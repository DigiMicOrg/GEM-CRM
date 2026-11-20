"""Genus x environment abundance prior from MicrobeAtlas (https://microbeatlas.org).

Built once, on the first run, into databases/microbeatlas/:
  taxon_env_prior.tsv   per SILVA genus and family and per environment class
                        (soil, soil;agricultural, animal;human, aquatic;marine, ...):
                        n_samples, mean relative abundance (zeros included), prevalence
  run_projects.tsv.gz   sequencing run -> study (project) for every MicrobeAtlas sample
  prior_meta.json       what was excluded when the table was built
The ~3.6 GB of downloads (OTU table, sample metadata, OTU->SILVA map) are deleted
after the build. Relative abundances use reads mapped to bacterial OTUs as the
denominator (unmapped reads, archaea, eukaryotes, chloroplasts and mitochondria
excluded) and only samples with >= 1000 such reads and a single environment.
Runs in the raw-data folder at build time, and every sample from their studies,
are excluded so the prior never contains the data it is applied to.
"""

import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from common import log

TABLE = "taxon_env_prior.tsv"
URLS = {
    "silva": "https://microbeatlas.org/downloads/mapref/otus.silva.csv",
    "env": "https://microbeatlas.org/downloads/samples/samples.env.info.tsv",
    "biom": "https://microbeatlas.org/downloads/otutables/samples-otus.97.mapped.metag.minfilter.refilt.biom.gz",
}
UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 "
      "(KHTML, like Gecko) Version/17.0 Safari/605.1.15")


def run_accessions(raw_dir):
    """Run accessions of the FASTQ files (ERR123_1.fastq.gz -> ERR123)."""
    out = set()
    for p in Path(raw_dir).rglob("*"):
        if p.is_file() and not p.name.startswith("."):
            out.add(p.name.split("_")[0].split(".")[0])
    return out


def _curl(url, dest):
    """curl with a browser user agent (the server rejects Python's) and resume support."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    part = dest.with_name(dest.name + ".part")
    log("Downloading", url)
    for attempt in range(5):
        r = subprocess.run(["curl", "-L", "--fail", "-sS", "-C", "-", "-A", UA, "--retry", "3",
                            "-o", str(part), url])
        if r.returncode == 0:
            part.replace(dest)
            return dest
        log(f"  curl exit {r.returncode}, retry {attempt + 1}/5")
        time.sleep(10)
    raise RuntimeError("MicrobeAtlas download failed: " + url)


def _free_gb(path):
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def build(folder, raw_dir, min_samples=30):
    """Download MicrobeAtlas data and write the prior table (one-off, ~10-30 min)."""
    try:
        import h5py
        import scipy.sparse as sp
    except ImportError as e:
        raise RuntimeError("h5py is required to build the MicrobeAtlas prior: "
                           "conda env update -f environment.yml (or pip install h5py)") from e
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    h5 = folder / "samples-otus.97.minfilter.biom.h5"
    if not h5.exists() and _free_gb(folder) < 16:
        raise RuntimeError(f"Building the MicrobeAtlas prior needs ~16 GB free disk "
                           f"(temporary); {_free_gb(folder):.1f} GB available")
    log("Building MicrobeAtlas genus x environment prior (one-off)")
    silva = _curl(URLS["silva"], folder / "otus.silva.csv")
    env = _curl(URLS["env"], folder / "samples.env.info.tsv")
    if not h5.exists():
        gz = _curl(URLS["biom"], folder / "samples-otus.97.minfilter.biom.gz")
        log("  decompressing OTU table")
        part = h5.with_name(h5.name + ".part")
        with gzip.open(gz, "rb") as i, open(part, "wb") as o:
            shutil.copyfileobj(i, o, 1 << 24)
        part.replace(h5)
        gz.unlink()
    with open(h5, "rb") as f:
        if f.read(8) != b"\x89HDF\r\n\x1a\n":
            raise RuntimeError("MicrobeAtlas OTU table is not HDF5 BIOM 2.1")

    otu_tax = {}
    with open(silva) as f:
        for row in csv.DictReader(f):
            ids = row["OTU"].split(";")
            if len(ids) == 3 and ids[-1].startswith("97_"):
                otu_tax[ids[-1]] = row["TaxSILVA"].split(";")

    runs = run_accessions(raw_dir)
    projects, sample_env = set(), {}
    with open(env) as f, gzip.open(folder / "run_projects.tsv.gz", "wt") as idx:   # pass 1
        rd = csv.reader(f, delimiter="\t")
        next(rd)
        for row in rd:
            run, project = row[0].split(".")[0], row[6] if len(row) > 6 else ""
            idx.write(f"{run}\t{project}\n")
            if run in runs and project:
                projects.add(project)
    with open(env) as f:                                                             # pass 2
        rd = csv.reader(f, delimiter="\t")
        next(rd)
        for row in rd:
            sid, envs = row[0], row[1] if len(row) > 1 else ""
            project = row[6] if len(row) > 6 else ""
            if not envs or sid.split(".")[0] in runs or project in projects:
                continue
            entries = [e.split(";") for e in envs.split("|") if e]
            tops = {e[0] for e in entries}
            if len(tops) != 1:
                continue  # several environments -> 'unknown' in MicrobeAtlas
            top = tops.pop()
            classes = [top]
            subs = {e[1] for e in entries if len(e) > 1}
            if len(subs) == 1:
                classes.append(f"{top};{subs.pop()}")
            sample_env[sid] = classes
    log(f"  {len(sample_env)} single-environment samples; excluded runs {sorted(runs)} and studies {sorted(projects)}")

    with h5py.File(h5, "r") as f:
        obs = [x.decode() if isinstance(x, bytes) else x for x in f["observation/ids"][:]]
        sids = [x.decode() if isinstance(x, bytes) else x for x in f["sample/ids"][:]]
        indptr = f["sample/matrix/indptr"][:]
        data, indices = f["sample/matrix/data"], f["sample/matrix/indices"]
        tax = [otu_tax.get(next((p for p in o.split(";") if p.startswith("97_")), o)) for o in obs]
        bact = np.array([bool(t) and t[0] == "Bacteria" and "Chloroplast" not in t and "Mitochondria" not in t
                         for t in tax], dtype=float)
        names = {}

        def tax_index(level):
            m = np.full(len(obs), -1, dtype=np.int64)
            for j, t in enumerate(tax):
                if t and len(t) > level and t[level]:
                    m[j] = names.setdefault((level, ";".join(t[:level + 1])), len(names))
            return m

        g_idx, f_idx = tax_index(5), tax_index(4)   # genus, family (0 = domain)
        ntax = len(names)
        env_names = sorted({c for v in sample_env.values() for c in v})
        ei = {e: k for k, e in enumerate(env_names)}
        S = np.zeros((len(env_names), ntax))
        P = np.zeros((len(env_names), ntax))
        N = np.zeros(len(env_names))
        cols = np.arange(len(obs))
        maps = []
        for idx in (g_idx, f_idx):
            ok = idx >= 0
            maps.append(sp.csr_matrix((np.ones(ok.sum()), (cols[ok], idx[ok])), shape=(len(obs), ntax)))
        step = 20000
        for s0 in range(0, len(sids), step):
            s1 = min(len(sids), s0 + step)
            a, b = indptr[s0], indptr[s1]
            X = sp.csr_matrix((data[a:b], indices[a:b], indptr[s0:s1 + 1] - a), shape=(s1 - s0, len(obs)))
            tot = X @ bact
            rows = [(ei[c], k - s0) for k in range(s0, s1) if tot[k - s0] >= 1000
                    for c in sample_env.get(sids[k], [])]
            if not rows:
                continue
            tot[tot == 0] = 1
            X = sp.diags(1 / tot) @ X
            er, sr = np.array(rows).T
            E = sp.csr_matrix((np.ones(len(er)), (er, sr)), shape=(len(env_names), s1 - s0))
            for M in maps:
                T = (X @ M).tocsr()
                S += (E @ T).toarray()
                P += (E @ (T > 0).astype(float)).toarray()
            N += np.bincount(er, minlength=len(env_names))
    inv = {v: k for k, v in names.items()}
    out = folder / TABLE
    part = out.with_name(out.name + ".part")
    with open(part, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["level", "taxon", "silva_lineage", "environment", "n_samples", "mean_rel_abundance", "prevalence"])
        for e, k in ei.items():
            if N[k] < min_samples:
                continue
            for t in np.flatnonzero(P[k] > 0):
                level, lineage = inv[t]
                w.writerow(["genus" if level == 5 else "family", lineage.split(";")[-1], lineage, e,
                            int(N[k]), f"{S[k, t] / N[k]:.6g}", f"{P[k, t] / N[k]:.6g}"])
    part.replace(out)
    (folder / "prior_meta.json").write_text(json.dumps(dict(
        source="MicrobeAtlas MAPv3, 97% OTUs, metag minfilter", urls=URLS, built=time.strftime("%Y-%m-%d %H:%M"),
        excluded_runs=sorted(runs), excluded_projects=sorted(projects),
        environments={e: int(N[k]) for e, k in ei.items() if N[k] >= min_samples}), indent=2))
    for name in ("otus.silva.csv", "samples.env.info.tsv", h5.name):
        (folder / name).unlink(missing_ok=True)
    log("  MicrobeAtlas prior ready:", out, f"({out.stat().st_size / 1e6:.0f} MB); raw downloads removed")
    return out


def ensure(folder, raw_dir, min_samples=30):
    table = Path(folder) / TABLE
    if table.exists() and (Path(folder) / "prior_meta.json").exists():
        return table
    return build(folder, raw_dir, min_samples)


class Prior:
    """Look up the MicrobeAtlas prior of SILVA genera in one environment class."""

    def __init__(self, table, environment):
        self.environment = environment
        self.by_lineage, self.by_name, envs = {}, {}, set()
        with open(table) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                envs.add(r["environment"])
                if r["environment"] != environment or r["level"] != "genus":
                    continue
                v = float(r["mean_rel_abundance"])
                self.by_lineage[r["silva_lineage"]] = v
                self.by_name.setdefault(r["taxon"], []).append(v)
        if not self.by_lineage:
            raise ValueError(f"Environment '{environment}' not in the MicrobeAtlas prior. Available: "
                             + ", ".join(sorted(envs)))
        meta = Path(table).parent / "prior_meta.json"
        self.meta = json.loads(meta.read_text()) if meta.exists() else {}

    def lookup(self, lineage, genus):
        """-> (mean relative abundance or None, how matched)."""
        if lineage in self.by_lineage:
            return self.by_lineage[lineage], "lineage"
        if len(self.by_name.get(genus, [])) == 1:
            return self.by_name[genus][0], "name"
        return None, "missing"

    def overlap_warning(self, folder, runs):
        """Studies of these runs that are inside the prior (built before they were added)."""
        index = Path(folder) / "run_projects.tsv.gz"
        if not index.exists():
            return []
        excluded = set(self.meta.get("excluded_projects", []))
        hits = set()
        with gzip.open(index, "rt") as f:
            for line in f:
                run, _, project = line.rstrip("\n").partition("\t")
                if run in runs and project and project not in excluded:
                    hits.add(project)
        return sorted(hits)
