"""Step 1: discover raw 16S FASTQ files and denoise each community with DADA2.

One FASTQ = one community.  Paired files (xxx_R1/xxx_R2, xxx_1/xxx_2) are
recognised; only the forward read is used, as in the original pipeline.
"""

import collections
import concurrent.futures
import gzip
import re
import subprocess
from pathlib import Path

from common import log, read_csv, tool, write_json

FASTQ = re.compile(r"\.(fastq|fq)(\.gz)?$", re.I)
PAIR = re.compile(r"^(?P<name>.+?)(?:[._-](?:R)?(?P<read>[12]))(?:_001)?$", re.I)


def discover(raw_dir):
    raw_dir = Path(raw_dir)
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"Raw data folder not found: {raw_dir}")
    found = {}
    for path in sorted(p for p in raw_dir.rglob("*") if p.is_file() and FASTQ.search(p.name)):
        stem = FASTQ.sub("", path.name)
        m = PAIR.match(stem)
        name, read = (m["name"], m["read"]) if m else (stem, "1")
        if read == "2":
            continue  # reverse read: not used (forward-only DADA2)
        if name in found:
            raise ValueError(f"Two forward FASTQ files map to community '{name}': {found[name]} and {path}")
        found[name] = path
    if not found:
        raise ValueError(f"No FASTQ files (*.fastq.gz, *.fq.gz, *.fastq, *.fq) in {raw_dir}")
    return found


def read_stats(path, limit=200000):
    """Modal read length and number of distinct quality symbols."""
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths = collections.Counter()
    quals = set()
    n = 0
    with opener(path, "rt") as f:
        while n < limit:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip()
            plus = f.readline()
            q = f.readline().strip()
            if not h.startswith("@") or not plus.startswith("+") or len(s) != len(q):
                raise ValueError(f"Malformed FASTQ record {n + 1} in {path}")
            lengths[len(s)] += 1
            if len(quals) < 64:
                quals.update(q)
            n += 1
    if n == 0:
        raise ValueError(f"Empty FASTQ: {path}")
    modal = max(lengths, key=lambda x: (lengths[x], x))
    return dict(reads_scanned=n, modal_length=modal, quality_symbols=len(quals))


def denoise_one(name, fastq, out, silva, d, threads):
    status = out / "status.txt"
    if status.exists() and status.read_text().startswith("complete") and (out / "asvs.csv").exists():
        return name, "cached"
    out.mkdir(parents=True, exist_ok=True)
    stats = read_stats(fastq)
    trunc = min(stats["modal_length"], int(d["max_trunc_len"]))
    # Near-constant quality scores (binned/stripped) break loess error fitting.
    model = "noqual" if stats["quality_symbols"] <= 2 else "loess"
    write_json(out / "read_stats.json", dict(stats, trunc_len=trunc, error_model=model, fastq=str(fastq)))
    script = Path(__file__).with_name("denoise.R")
    with (out / "dada2.log").open("w") as logf:
        code = subprocess.run(
            [tool("Rscript"), "--vanilla", str(script), str(fastq), str(out), name, str(silva),
             str(trunc), str(d.get("trim_left", 0)), model, str(threads), str(d["min_reads_after_filter"])],
            stdout=logf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
        ).returncode
    return name, ("ok" if code == 0 else "failed")


def run(cfg, root, silva):
    raw = root / cfg["paths"]["raw_data"]
    samples = discover(raw)
    log(f"Found {len(samples)} communities in {raw}")
    results = root / cfg["paths"]["results"]
    d = cfg["denoise"]
    workers = max(1, int(cfg["compute"]["parallel_denoise"]))
    threads = max(1, int(cfg["compute"]["threads"]) // workers)
    done = {}
    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(denoise_one, n, p, results / n / "01_dada2", silva, d, threads)
                   for n, p in samples.items()]
        for fut in concurrent.futures.as_completed(futures):
            name, state = fut.result()
            done[name] = state
            log("DADA2", name, state)
    return {n: results / n / "01_dada2" for n, s in done.items() if s in ("ok", "cached")}


def load_community(folder):
    """Return ASV table rows and per-ASV counts for a denoised community."""
    asvs = read_csv(folder / "asvs.csv")
    counts = {r["ASV"]: float(r["count"]) for r in read_csv(folder / "counts.csv")}
    return asvs, counts
