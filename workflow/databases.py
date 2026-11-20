"""Tool check and automatic download of public reference databases.

Downloaded once into the databases/ folder and reused afterwards:
  * SILVA 138.2 DADA2 genus training set (taxonomy assignment)
  * Madin et al. 2020 condensed species trait table (cell size -> m)
  * gapseq sequence database (via `gapseq update-sequences`, optional)
  * MicrobeAtlas genus x environment abundance prior (see prior.py)
"""

import shutil
import subprocess
import urllib.request
from pathlib import Path

from common import command, digest, log, tool, write_json, read_json

REQUIRED_TOOLS = ["Rscript", "datasets", "barrnap", "mash", "gapseq"]


def check_tools(steps):
    need = set()
    if "denoise" in steps:
        need.add("Rscript")
    if "genomes" in steps:
        need |= {"datasets", "barrnap", "mash"}
    if "gems" in steps:
        need.add("gapseq")
    missing = [t for t in sorted(need) if not shutil.which(t)]
    if missing:
        raise RuntimeError(
            "Missing executables: " + ", ".join(missing)
            + ". Create and activate the environment: conda env create -f environment.yml; conda activate micrm16s"
        )
    if "denoise" in steps:
        ok = subprocess.run(
            [tool("Rscript"), "--vanilla", "-e",
             "quit(status=if(requireNamespace('dada2',quietly=TRUE))0 else 1)"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        ).returncode == 0
        if not ok:
            raise RuntimeError("R package dada2 is not installed in this environment")


def download(url, dest, label):
    dest = Path(dest)
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    log("Downloading", label, "->", dest)
    part = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "micrm-workflow"})
    with urllib.request.urlopen(req, timeout=300) as src, part.open("wb") as out:
        while chunk := src.read(1 << 20):
            out.write(chunk)
    if part.stat().st_size == 0:
        raise RuntimeError("Empty download: " + url)
    part.replace(dest)
    return dest


def ensure_databases(cfg, root, steps):
    db = root / cfg["paths"]["databases"]
    db.mkdir(parents=True, exist_ok=True)
    out = {}
    manifest = read_json(db / "manifest.json", {})
    if "denoise" in steps:
        d = cfg["denoise"]
        path = Path(d["silva_trainset"]) if d.get("silva_trainset") else db / Path(d["silva_url"].split("?")[0]).name
        if not path.is_absolute():
            path = root / path
        out["silva"] = download(d["silva_url"], path, "SILVA trainset")
        manifest.setdefault("silva", dict(url=d["silva_url"], file=str(path.name), sha256=digest(path)))
    if "params" in steps:
        p = cfg["parameters"]
        path = Path(p["madin_traits"]) if p.get("madin_traits") else db / "madin_condensed_species_NCBI.csv"
        if not path.is_absolute():
            path = root / path
        out["madin"] = download(p["madin_url"], path, "Madin trait table")
        manifest.setdefault("madin", dict(url=p["madin_url"], file=str(path.name), sha256=digest(path)))
        import prior
        pr = cfg.get("prior", {})
        out["prior"] = prior.ensure(db / "microbeatlas", root / cfg["paths"]["raw_data"], int(pr.get("min_samples", 30)))
    if "gems" in steps:
        g = cfg["gems"]
        if g.get("gapseq_db"):
            out["gapseq_db"] = str(Path(g["gapseq_db"]).expanduser().resolve())
        elif g.get("update_gapseq_sequences", True):
            marker = db / "gapseq_sequences.updated"
            if not marker.exists():
                log("Updating gapseq reference sequences (one-off)")
                try:
                    command([tool("gapseq"), "update-sequences"], db / "gapseq_update.log")
                    marker.write_text("ok\n")
                except Exception as e:  # gapseq may already ship its sequences
                    log("WARNING gapseq update-sequences failed; continuing with bundled data:", str(e)[:200])
    write_json(db / "manifest.json", manifest)
    return out
