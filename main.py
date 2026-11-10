#!/usr/bin/env python3
"""One-click workflow: raw 16S FASTQ -> MiCRM simulation vs observed abundance.

    python main.py                      # run everything
    python main.py --check              # only check tools and download databases
    python main.py --steps params,simulate --rerun   # redo parameters + simulation
    python main.py --communities S1,S2  # restrict to some communities

Steps
  1 denoise   DADA2 denoising + SILVA genus assignment         (workflow/denoise.R)
  2 taxa      collect genera and their ASV sequences across communities
  3 genomes   NCBI reference genome per genus (+related genus / family fallback)
  4 gems      gapseq GEM per selected genome (cached in gems/)
  5 params    binary uptake, resource panel, K (literature), m (cell size)
  6 simulate  U (min carbon-throughput LP), L (max output entropy), abundance
              calibration, BDF MiCRM, comparison bar chart
"""

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import collections
import concurrent.futures
import contextlib
import shutil
import sys
import tomllib
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "workflow"))

from common import digest, log, read_csv, read_json, write_csv, write_json  # noqa: E402
import databases  # noqa: E402
import denoise  # noqa: E402

STEPS = ["denoise", "taxa", "genomes", "gems", "params", "simulate"]


def load_config(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


# --------------------------------------------------------------------------- 2
def collect_taxa(dada2_dirs):
    """genus -> dict(asvs={ASV: seq}, family=most common SILVA family)."""
    info = {}
    families = collections.defaultdict(collections.Counter)
    for name, folder in sorted(dada2_dirs.items()):
        asvs, _ = denoise.load_community(folder)
        for r in asvs:
            if r.get("Eligible", "").upper() != "TRUE" or not r.get("Genus"):
                continue
            g = info.setdefault(r["Genus"], dict(asvs={}, family=""))
            g["asvs"][r["ASV"]] = r["Sequence"]
            if r.get("Family"):
                families[r["Genus"]][r["Family"]] += 1
    for g, c in families.items():
        info[g]["family"] = c.most_common(1)[0][0]
    return info


# --------------------------------------------------------------------------- 6
def simulate_one(job):
    """Run in a worker process: estimate U/L, calibrate, simulate, plot."""
    import argparse as _a
    import micrm
    name, inputs, catalog, cache, out, mcfg = job
    log_path = out.parent / "micrm.log"
    try:
        if out.exists():
            shutil.rmtree(out)
        with log_path.open("w") as f, contextlib.redirect_stdout(f):
            micrm.run(_a.Namespace(input=str(inputs), catalog=str(catalog), cache=str(cache), out=str(out),
                                   initial_total=float(mcfg["initial_total"]), feed_total=float(mcfg["feed_total"]),
                                   biomass=float(mcfg["biomass"]), turnover=float(mcfg["turnover"]),
                                   t_end=float(mcfg["t_end"]), title=name))
        return name, "completed", read_json(out / "metadata.json"), ""
    except Exception as e:
        with log_path.open("a") as f:
            f.write(traceback.format_exc())
        return name, "failed", None, str(e)[:300]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=str(ROOT / "config.toml"))
    p.add_argument("--steps", default=",".join(STEPS), help="comma-separated subset of: " + ",".join(STEPS))
    p.add_argument("--communities", default="", help="comma-separated community names (default: all)")
    p.add_argument("--rerun", action="store_true", help="redo simulations that already finished")
    p.add_argument("--check", action="store_true", help="check tools and download databases, then stop")
    a = p.parse_args()

    cfg = load_config(a.config)
    steps = [s.strip() for s in a.steps.split(",") if s.strip()]
    if bad := set(steps) - set(STEPS):
        p.error("unknown steps: " + ", ".join(sorted(bad)))
    results = ROOT / cfg["paths"]["results"]
    results.mkdir(parents=True, exist_ok=True)
    for key in ("raw_data", "databases", "genomes", "gems"):
        (ROOT / cfg["paths"][key]).mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- setup
    databases.check_tools(STEPS if a.check else steps)
    dbs = databases.ensure_databases(cfg, ROOT, STEPS if a.check else steps)
    if a.check:
        log("Tools and databases OK:", {k: str(v) for k, v in dbs.items()})
        return 0
    if dbs.get("gapseq_db"):
        cfg["gems"]["gapseq_db"] = dbs["gapseq_db"]
        os.environ["UL_GAPSEQ_DB"] = dbs["gapseq_db"]

    # ------------------------------------------------------- 1 read + denoise
    if "denoise" in steps:
        log("== Step 1: DADA2 denoising + genus assignment")
        dada2 = denoise.run(cfg, ROOT, dbs["silva"])
    else:
        dada2 = {d.parent.name: d for d in sorted(results.glob("*/01_dada2"))
                 if (d / "status.txt").exists() and (d / "status.txt").read_text().startswith("complete")}
    if a.communities:
        keep = {c.strip() for c in a.communities.split(",")}
        dada2 = {k: v for k, v in dada2.items() if k in keep}
    if not dada2:
        raise SystemExit("No successfully denoised communities (see results/<community>/01_dada2/).")
    log(f"{len(dada2)} denoised communities")

    # ----------------------------------------------------------- 2 genera
    genus_info = collect_taxa(dada2)
    log(f"== Step 2: {len(genus_info)} bacterial genera across communities")
    write_csv(results / "genera.csv",
              [dict(genus=g, family=i["family"], asvs=len(i["asvs"])) for g, i in sorted(genus_info.items())],
              ["genus", "family", "asvs"])

    # ---------------------------------------------------------- 3 genomes
    if "genomes" in steps:
        import genomes
        log("== Step 3: reference genomes from NCBI")
        refs = genomes.run(cfg, ROOT, genus_info)
        log(f"{len(refs)} genera have a reference genome")
    elif "gems" in steps:
        refs = {r["genus"]: r for r in read_csv(results / "genus_references.csv") if r["status"] == "ok"}

    # ------------------------------------------------------------- 4 GEMs
    if "gems" in steps:
        import gems
        log("== Step 4: gapseq GEM reconstruction")
        catalog, built = gems.run(cfg, ROOT, refs)
    else:
        catalog = results / "gem_catalog.csv"
        built = [r for r in read_csv(catalog) if r["status"] == "ok"]
    log(f"{len(built)} genera have a usable GEM")

    # ----------------------------------------------------------- 5 params
    status = {n: dict(community=n, status="denoised") for n in dada2}
    jobs = []
    if "params" in steps or "simulate" in steps:
        import maintenance
        import parameterize
        log("== Step 5: binary uptake, resource panel, K and m")
        pc = cfg["parameters"]
        inv, errors = parameterize.inventory(built, set(pc["carbon_sinks"]))
        if errors:
            write_json(results / "model_errors.json", errors)
        models = {r["genus"]: (r["accession"], inv[r["accession"]]) for r in built if r["accession"] in inv}
        write_json(results / "model_inventory.json",
                   {g: dict(accession=a_, resources=res) for g, (a_, res) in models.items()})
        priors, fallback = parameterize.kinetic_priors()
        m_rows = maintenance.estimate(dbs.get("madin") or databases.ensure_databases(cfg, ROOT, ["params"])["madin"],
                                      list(models), pc["missing_m_policy"])
        write_csv(results / "maintenance.csv", m_rows, ["genus", "volume_um3", "m", "source", "resolved"])
        mmap = {r["genus"]: r for r in m_rows}
        for name, folder in sorted(dada2.items()):
            inputs = results / name / "02_inputs"
            try:
                asvs, counts = denoise.load_community(folder)
                status[name].update(parameterize.community_inputs(
                    name, asvs, counts, models, mmap, priors, fallback, inputs, int(pc["min_genera"])))
                status[name]["status"] = "inputs_ready"
                jobs.append((name, inputs, catalog, ROOT / cfg["paths"]["gems"] / "parameter_cache",
                             results / name / "03_micrm", cfg["micrm"]))
            except Exception as e:
                status[name].update(status="skipped", reason=str(e))
                log("skip", name, e)

    # --------------------------------------------------------- 6 simulate
    if "simulate" in steps:
        log(f"== Step 6: U/L estimation + MiCRM for {len(jobs)} communities")
        figures = results / "figures"
        figures.mkdir(exist_ok=True)
        todo = []
        for job in jobs:
            name, out = job[0], job[4]
            meta = read_json(out / "metadata.json")
            same_inputs = meta and all(
                (job[1] / f).exists() and digest(job[1] / f) == h for f, h in meta.get("input_hashes", {}).items())
            if same_inputs and not a.rerun:
                status[name].update(status="completed", jsd=meta["jsd"])
            else:
                todo.append(job)
        workers = max(1, int(cfg["compute"]["parallel_communities"]))
        with concurrent.futures.ProcessPoolExecutor(workers) as pool:
            for name, state, meta, err in pool.map(simulate_one, todo):
                status[name]["status"] = state
                if meta:
                    status[name]["jsd"] = meta["jsd"]
                    shutil.copy2(results / name / "03_micrm" / "comparison.png", figures / f"{name}.png")
                if err:
                    status[name]["reason"] = err
                log("MiCRM", name, state, f"JSD={meta['jsd']:.3f}" if meta else err)

    rows = list(status.values())
    write_csv(results / "summary.csv", rows, ["community", "status", "genera", "resources", "coverage", "jsd", "reason"])
    counts = collections.Counter(r["status"] for r in rows)
    log("Done:", dict(counts), "->", results / "summary.csv", "| figures in", results / "figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
