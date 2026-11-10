"""Step 5: prepare MiCRM inputs for each community (logic from hpc_batch/run_batch.py).

* binary uptake capability (0/1): a GEM can take up a carbon resource when the
  resource has an extracellular exchange, positive carbon, a supporting
  transport reaction, and is not a carbon sink (CO2/HCO3- by default)
* resource panel: union of resources taken up by the community's modeled
  genera; every panel resource is supplied (fed)
* K: substrate-level geometric mean of audited whole-cell uptake literature
  records (mM-C); other substrates use the geometric mean of those priors.
  The same K applies to every genus for a given resource.
* m: Madin cell size -> volume -> m = 0.01449 V^-0.09
* observed: all genus counts including unclassified ASVs (coverage denominator)
"""

import collections
import math
from pathlib import Path

from common import digest, read_json, write_csv, write_json
from kinetic_helpers import SUBSTRATES, gmean, resource_prior
from network import Network

EVIDENCE = Path(__file__).resolve().parent / "resources" / "kinetics_evidence.csv"


def inventory(catalog_rows, sinks):
    """accession -> sorted list of structurally importable carbon resources."""
    out, errors = {}, {}
    for row in catalog_rows:
        acc = row["accession"]
        if acc in out or acc in errors:
            continue
        xml = Path(row["xml_path"])
        cache = xml.parent / "uptake_inventory.json"
        key = dict(xml_sha256=digest(xml), sinks=sorted(sinks))
        cached = read_json(cache)
        if cached and cached.get("key") == key:
            out[acc] = cached["resources"]
            continue
        try:
            net = Network(xml)
            res = sorted(r for r, (_, k, _) in net.exchanges.items()
                         if r not in sinks and math.isfinite(net.carbon[k]) and net.carbon[k] > 0
                         and net.transport_supported(r))
            write_json(cache, dict(key=key, resources=res))
            out[acc] = res
        except Exception as e:
            errors[acc] = str(e)
    return out, errors


def kinetic_priors(evidence_path=EVIDENCE):
    from common import read_csv
    evidence = read_csv(evidence_path)
    priors = {r: resource_prior(r, evidence) for r in SUBSTRATES.values()}
    fallback = gmean([v["K"] for v in priors.values() if v])
    if not fallback:
        raise ValueError("No usable K evidence")
    return priors, fallback


def community_inputs(name, asvs, counts, models, maintenance, priors, fallback, out, min_genera):
    """models: genus -> (accession, [resources]). Returns summary dict (raises when unusable)."""
    by_asv = {r["ASV"]: r for r in asvs}
    observed = collections.Counter()
    for asv, c in counts.items():
        observed[(by_asv.get(asv, {}).get("Genus") or "__unclassified__")] += c
    eligible = {r["Genus"] for r in asvs if r.get("Eligible", "").upper() == "TRUE"}
    selected = sorted(g for g in eligible if g in models and observed[g] > 0 and models[g][1])
    if len(selected) < min_genera:
        raise ValueError(f"Fewer than {min_genera} positive-abundance genera with usable GEMs")
    resources = sorted(set().union(*(models[g][1] for g in selected)))
    write_csv(out / "genera.csv", [dict(genus=g, accession=models[g][0]) for g in selected], ["genus", "accession"])
    write_csv(out / "resources.csv", [dict(resource=r, fed=True) for r in resources], ["resource", "fed"])
    write_csv(out / "observed.csv", [dict(genus=g, count=c) for g, c in sorted(observed.items())], ["genus", "count"])
    write_csv(out / "maintenance.csv", [maintenance[g] for g in selected], ["genus", "m", "source", "volume_um3"])
    krows = [dict(genus=g, resource=r, K=priors[r]["K"] if priors.get(r) else fallback,
                  source="substrate_evidence_gmean" if priors.get(r) else "global_evidence_gmean")
             for g in selected for r in sorted(models[g][1])]
    write_csv(out / "kinetics.csv", krows, ["genus", "resource", "K", "source"])
    # binary (0/1) structural uptake matrix, for inspection
    write_csv(out / "binary_uptake.csv",
              [dict(genus=g, **{r: int(r in set(models[g][1])) for r in resources}) for g in selected],
              ["genus"] + resources)
    total = sum(observed.values())
    return dict(community=name, genera=len(selected), resources=len(resources),
                coverage=sum(observed[g] for g in selected) / total if total else 0)
