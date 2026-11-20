"""Step 5: prepare MiCRM inputs for each community (logic from hpc_batch/run_batch.py).

* binary uptake capability (0/1): a genus can take up a carbon resource when the
  final GEM has an exchange + transport for it, the transport is already in the
  draft GEM (not gap-filled), gapseq found a transporter sequence for it, it
  contains carbon and is not a carbon sink (CO2/HCO3- by default); see inventory()
* resource panel: union of resources taken up by the community's modeled
  genera; every panel resource may be supplied (amounts inferred in step 6)
* K: substrate-level geometric mean of audited whole-cell uptake literature
  records (mM-C); other substrates use the geometric mean of those priors.
  The same K applies to every genus for a given resource.
* m: Madin cell size -> volume -> m = 0.01449 V^-0.09
* observed: all genus counts including unclassified ASVs (coverage denominator)
* prior: each genus's mean relative abundance in the configured environment class
  of the MicrobeAtlas prior (matched by SILVA lineage, then by unique genus name;
  genera not found get the median of the found ones)
"""

import collections
import csv
import math
import re
from pathlib import Path

from common import digest, read_json, write_csv, write_json
from kinetic_helpers import SUBSTRATES, gmean, resource_prior
from network import Network

EVIDENCE = Path(__file__).resolve().parent / "resources" / "kinetics_evidence.csv"


def _importable(net, sinks):
    return {r for r, (_, k, _) in net.exchanges.items()
            if r not in sinks and math.isfinite(net.carbon[k]) and net.carbon[k] > 0
            and net.transport_supported(r)}


def _transporter_hits(folder):
    """Compounds with at least one transporter sequence hit in gapseq's
    find-transport table (genome-Transporter.tbl)."""
    path = Path(folder) / "genome-Transporter.tbl"
    if not path.exists():
        return None
    with path.open() as f:
        rows = [line for line in f if not line.startswith("#")]
    hits = set()
    for r in csv.DictReader(rows, delimiter="\t"):
        m = re.match(r"EX_(cpd\d+)_e0$", (r.get("exid") or "").strip())
        if m:
            hits.add(m.group(1))
    return hits


UPTAKE_RULE = "sequence_evidence_v1"


def inventory(catalog_rows, sinks):
    """accession -> sorted list of carbon resources the genus can import.

    A compound counts only when all three hold:
      1. the final GEM has an extracellular exchange and a transport reaction for it;
      2. that transport is already in gapseq's draft model (not added by gap-filling);
      3. gapseq's find-transport found a transporter sequence for it
         (genome-Transporter.tbl).
    Without 2-3, gap-filled models credit genera with uptake they have no
    genomic evidence for; in the MiCRM those extra niches made the genera with
    the most complete models win regardless of their real abundance
    (cf. misosoup: fundamental niche size tracks model size).
    """
    out, errors = {}, {}
    for row in catalog_rows:
        acc = row["accession"]
        if acc in out or acc in errors:
            continue
        xml = Path(row["xml_path"])
        draft = xml.parent / "genome-draft.xml"
        cache = xml.parent / "uptake_inventory.json"
        key = dict(xml_sha256=digest(xml), sinks=sorted(sinks), rule=UPTAKE_RULE)
        cached = read_json(cache)
        if cached and cached.get("key") == key:
            out[acc] = cached["resources"]
            continue
        try:
            hits = _transporter_hits(xml.parent)
            if hits is None or not draft.exists():
                raise FileNotFoundError("genome-draft.xml or genome-Transporter.tbl missing; rebuild this GEM")
            res = _importable(Network(xml), sinks) & _importable(Network(draft), sinks) & hits
            res = sorted(res)
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


def community_inputs(name, asvs, counts, models, maintenance, priors, fallback, out, min_genera, prior):
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
    # MicrobeAtlas prior (s, the genus's uptake capacity, is derived from it in step 6)
    lineages = collections.defaultdict(collections.Counter)
    for r in asvs:
        if r.get("Genus"):
            lineages[r["Genus"]][";".join(r.get(k) or "" for k in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))] += 1
    found = {g: prior.lookup(lineages[g].most_common(1)[0][0] if lineages[g] else "", g) for g in selected}
    values = [v for v, _ in found.values() if v is not None]
    fill = sorted(values)[len(values) // 2] if values else 1.0
    write_csv(out / "prior.csv",
              [dict(genus=g, prior=found[g][0] if found[g][0] is not None else fill,
                    match=found[g][1] if found[g][0] is not None else "missing_median")
               for g in selected], ["genus", "prior", "match"])
    # binary (0/1) structural uptake matrix, for inspection
    write_csv(out / "binary_uptake.csv",
              [dict(genus=g, **{r: int(r in set(models[g][1])) for r in resources}) for g in selected],
              ["genus"] + resources)
    total = sum(observed.values())
    return dict(community=name, genera=len(selected), resources=len(resources),
                coverage=sum(observed[g] for g in selected) / total if total else 0,
                prior_matched=f"{sum(v is not None for v, _ in found.values())}/{len(selected)}")
