import gzip, re, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, vstack, hstack, csr_matrix
from scipy.optimize import linprog


def local(s):
    return s.split("}")[-1]


def attrs(x):
    return {local(k): v for k, v in x.attrib.items()}


def formula(f):
    if not f or f in ("NA", "None") or re.sub("([A-Z][a-z]?)([0-9]*)", "", f):
        return None
    d = {}
    for e, n in re.findall("([A-Z][a-z]?)([0-9]*)", f):
        d[e] = d.get(e, 0) + int(n or 1)
    return d


class Network:

    def __init__(self, path, config=None):
        self.config = config or {}
        self.path = str(path)
        with (
            gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")
        ) as f:
            root = ET.parse(f).getroot()
        sp = [attrs(x) for x in root.iter() if local(x.tag) == "species"]
        self.species = sp
        self.sid = [x["id"] for x in sp]
        ix = {s: i for i, s in enumerate(self.sid)}
        self.formula = [formula(s.get("chemicalFormula")) for s in sp]
        self.carbon = np.array(
            [f.get("C", 0) if f is not None else np.nan for f in self.formula]
        )
        params = {
            x.get("id"): float(x.get("value"))
            for x in root.iter()
            if local(x.tag) == "parameter" and x.get("value") is not None
        }
        rx = [x for x in root.iter() if local(x.tag) == "reaction"]
        self.rid = []
        self.lb = []
        self.ub = []
        rr = []
        cc = []
        vv = []
        self.exchanges = {}
        self.boundaries = {}
        self.bio = None
        for j, x in enumerate(rx):
            a = attrs(x)
            self.rid.append(a["id"])
            self.lb.append(params[a["lowerFluxBound"]])
            self.ub.append(params[a["upperFluxBound"]])
            if a.get("sboTerm") == "SBO:0000629" or a["id"] == "R_bio1":
                self.bio = j
            pairs = {}
            for side in x:
                sign = {"listOfReactants": -1, "listOfProducts": 1}.get(local(side.tag))
                if sign is None:
                    continue
                for s in side:
                    k = ix[s.get("species")]
                    pairs[k] = pairs.get(k, 0) + sign * float(
                        s.get("stoichiometry", "1")
                    )
            pairs = {k: v for k, v in pairs.items() if v}
            for k, v in pairs.items():
                rr.append(k)
                cc.append(j)
                vv.append(v)
            if len(pairs) == 1:
                k, v = next(iter(pairs.items()))
                self.boundaries[j] = (k, v)
                if (
                    sp[k]["compartment"].lower().startswith("e")
                    or sp[k]["compartment"].lower().endswith("e")
                ) and a["id"].startswith(("R_EX_", "EX_")):
                    resource = re.sub(
                        "_[a-z][0-9]+$", "", re.sub("^M_", "", self.sid[k])
                    )
                    self.exchanges[resource] = (j, k, v)
        self.S = coo_matrix((vv, (rr, cc)), shape=(len(sp), len(rx))).tocsc()
        self.lb = np.array(self.lb)
        self.ub = np.array(self.ub)
        if self.bio is None:
            raise ValueError("No biomass reaction")
        self.pseudo = [i for i, s in enumerate(self.sid) if s.startswith("M_cpd11416_")]
        self.balance_rows = [i for i in range(len(sp)) if i not in self.pseudo]
        self.A = self.S[self.balance_rows, :]
        col = self.S[:, self.bio].toarray().ravel()
        self.biomass_unknown = [
            self.sid[i]
            for i in np.flatnonzero(col)
            if self.formula[i] is None and i not in self.pseudo
        ]
        self.biomass_C = -float(np.nansum(col * self.carbon))
        self.symbolic_residual = {}
        for sym in ["R", "X"]:
            z = sum(
                (
                    col[i] * f.get(sym, 0)
                    for i, f in enumerate(self.formula)
                    if f is not None
                )
            )
            if abs(z) > 1e-08:
                self.symbolic_residual[sym] = float(z)
        self.bad_reactions = []
        for j in range(len(rx)):
            if j == self.bio or j in self.boundaries:
                continue
            col = self.S[:, j]
            ids = col.indices
            vals = col.data
            unknown = any((self.formula[k] is None for k in ids))
            carbon = sum((v * self.carbon[k] for k, v in zip(ids, vals)))
            symbad = any(
                (
                    abs(
                        sum(
                            (
                                v * self.formula[k].get(s, 0)
                                for k, v in zip(ids, vals)
                                if self.formula[k] is not None
                            )
                        )
                    )
                    > 1e-07
                    for s in ["R", "X"]
                )
            )
            if unknown or symbad or abs(carbon) > 1e-07:
                self.bad_reactions.append(j)

    def bounds(self, rich=False, strict=True):
        lo = self.lb.copy()
        hi = self.ub.copy()
        if strict:
            lo[self.bad_reactions] = 0
            hi[self.bad_reactions] = 0
        for j, (k, s) in self.boundaries.items():
            lo[j] = hi[j] = 0
        for res, (j, k, s) in self.exchanges.items():
            c = self.carbon[k]
            if (
                not np.isfinite(c)
                or self.formula[k].get("R", 0)
                or self.formula[k].get("X", 0)
            ):
                continue
            inward = (
                self.config.get("network_inorganic_import_bound", 10)
                if rich or c == 0
                else 0
            )
            outward = self.config.get("network_export_bound", 1000)
            if s < 0:
                lo[j] = -inward / abs(s)
                hi[j] = outward / abs(s)
            else:
                lo[j] = -outward / abs(s)
                hi[j] = inward / abs(s)
        return (lo, hi)

    def transport_supported(self, resource):
        j, k, s = self.exchanges[resource]
        row = self.S.getrow(k)
        return any(
            (
                jj != j and (v < 0 and self.ub[jj] > 0 or (v > 0 and self.lb[jj] < 0))
                for jj, v in zip(row.indices, row.data)
            )
        )

    def precursor_indices(self):
        """Structural biomass precursors only.

        Energy carriers/cofactors are excluded: ATP alone carries ~92% of the
        nominal 'precursor carbon' in gapseq biomass reactions and is
        regenerated rather than retained, so leaving it in makes every yield
        and leakage calculation describe ATP turnover instead of biomass.
        """
        allowed = set(self.config.get("retention_precursors", []))
        cofactors = {
            re.sub("_[a-z][0-9]+$", "", x)
            for x in self.config.get(
                "retention_cofactor_exclusions",
                [
                    "cpd00001",
                    "cpd00002",
                    "cpd00003",
                    "cpd00004",
                    "cpd00005",
                    "cpd00006",
                    "cpd00008",
                    "cpd00009",
                    "cpd00010",
                    "cpd00067",
                    "cpd00011",
                    "cpd00012",
                ],
            )
        }
        max_coef = float(self.config.get("max_precursor_coefficient", 1.0))
        biomass = self.S[:, self.bio].toarray().ravel()
        out = []
        for k, s in enumerate(self.species):
            if s["compartment"] != "c0" or biomass[k] >= 0:
                continue
            coef = -float(biomass[k])
            if coef <= 0 or coef > max_coef:
                continue
            rid = re.sub("_[a-z][0-9]+$", "", re.sub("^M_", "", self.sid[k]))
            if rid in cofactors:
                continue
            if allowed and rid not in allowed:
                continue
            if self.formula[k] is None or self.carbon[k] <= 0:
                continue
            if self.formula[k].get("R", 0) or self.formula[k].get("X", 0):
                continue
            out.append(k)
        return out
