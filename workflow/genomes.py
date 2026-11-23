"""Step 3: find one NCBI reference genome for every 16S genus.

Order of attempts for a SILVA genus label (first that yields genomes wins):
  1. the exact genus name (NCBI taxonomy, restricted to Bacteria)
  2. explicit aliases for composite SILVA labels (e.g. Escherichia-Shigella)
  3. related genus names derived from the label ("[Eubacterium] xx group" ->
     Eubacterium, "Clostridium sensu stricto 1" -> Clostridium, hyphen parts)
  4. family-level surrogate: reference genomes of other genera of the same
     SILVA family (one per genus)
  5. otherwise the genus is skipped (reported as no_reference)

Within a tier genomes are searched as NCBI reference -> type material ->
any annotated assembly.  Among the candidates, barrnap extracts 16S copies and
Mash picks the genome closest to the genus' ASV sequences.

Everything is cached under genomes/ so reruns and new datasets reuse downloads.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from common import (command, digest, fasta, log, read_json, safe_name, text_hash,
                    tool, write_csv, write_fasta, write_json)

# Explicit composite-label mapping.
TAXON_ALIASES = {
    "Escherichia-Shigella": ["Escherichia", "Shigella"],
    "Burkholderia-Caballeronia-Paraburkholderia": ["Burkholderia", "Caballeronia", "Paraburkholderia"],
}

TIERS = [
    ("reference", ["--reference", "--assembly-level", "complete,chromosome"]),
    ("type", ["--from-type", "--annotated", "--assembly-level", "complete,chromosome", "--exclude-atypical"]),
    ("annotated", ["--annotated", "--assembly-level", "complete,chromosome,scaffold", "--exclude-atypical"]),
]

GENOME_FIELDS = ["genus", "accession", "source", "organism", "tier", "mash_distance",
                 "matched_asv", "candidates", "genome_path", "status", "note"]


class NCBI:
    def __init__(self, root, cfg):
        self.root = Path(root)
        self.api_key = cfg.get("ncbi_api_key") or os.environ.get("NCBI_API_KEY", "")
        self.max_candidates = int(cfg["max_candidates"])
        self.delay = 0.12 if self.api_key else 0.4
        for d in ("taxonomy", "summaries", "ncbi", "selection", "logs"):
            (self.root / d).mkdir(parents=True, exist_ok=True)

    # ---- taxonomy -------------------------------------------------------
    def taxid(self, name, rank):
        path = self.root / "taxonomy" / f"{safe_name(name)}__{rank}.json"
        cached = read_json(path)
        if cached is not None:
            return cached.get("taxid")
        term = f'"{name}"[Scientific Name] AND {rank}[Rank] AND txid2[Subtree]'
        params = dict(db="taxonomy", term=term, retmode="json", retmax=20)
        if self.api_key:
            params["api_key"] = self.api_key
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
        last = None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=45) as f:
                    data = json.load(f)["esearchresult"]
                last = None
                break
            except Exception as e:  # transient network/rate errors
                last = e
                time.sleep(1.5 * (attempt + 1))
        if last is not None:
            raise RuntimeError(f"NCBI taxonomy query failed for {name}: {last}")
        time.sleep(self.delay)
        ids = data.get("idlist", [])
        taxid = ids[0] if len(ids) == 1 else None
        write_json(path, dict(name=name, rank=rank, query=term, idlist=ids, taxid=taxid))
        return taxid

    # ---- genome listing -------------------------------------------------
    def _datasets(self, args, log_name):
        args = [tool("datasets")] + args + (["--api-key", self.api_key] if self.api_key else [])
        res = subprocess.run(args, capture_output=True, text=True, stdin=subprocess.DEVNULL)
        with (self.root / "logs" / log_name).open("a") as f:
            f.write(json.dumps(args[:-1] if self.api_key else args) + "\n" + res.stderr[-4000:] + "\n")
        return res

    def summary(self, taxid, tier, flags, limit):
        path = self.root / "summaries" / f"{taxid}__{tier}__{limit}.json"
        cached = read_json(path)
        if cached is not None:
            return cached
        res = self._datasets(["summary", "genome", "taxon", str(taxid)] + flags +
                             ["--limit", str(limit), "--as-json-lines"], "datasets_summary.log")
        rows = []
        for line in res.stdout.splitlines():
            line = line.strip()
            if not line.startswith("{"):
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            acc = rec.get("accession") or rec.get("current_accession")
            if acc and re.fullmatch(r"GC[AF]_\d+\.\d+", acc):
                rows.append(dict(accession=acc, organism=(rec.get("organism") or {}).get("organism_name", ""), tier=tier))
        if res.returncode == 0 or rows:
            write_json(path, rows)  # cache empty answers only when the call succeeded
        return rows

    def candidates(self, taxid, diversify=False):
        """First non-empty tier; with diversify keep one genome per genus (family level)."""
        limit = self.max_candidates * (10 if diversify else 1)
        for tier, flags in TIERS:
            rows = self.summary(taxid, tier, flags, limit)
            if diversify:
                seen, kept = set(), []
                for r in rows:
                    g = (r["organism"].split() or [""])[0]
                    if g and g not in seen:
                        seen.add(g)
                        kept.append(r)
                rows = kept
            if rows:
                return rows[: self.max_candidates]
        return []

    # ---- download + 16S --------------------------------------------------
    def genome(self, accession):
        folder = self.root / "ncbi" / accession
        target = folder / "genome.fna"
        if target.is_file() and target.stat().st_size:
            return target
        folder.mkdir(parents=True, exist_ok=True)
        archive = folder / "download.zip"
        res = self._datasets(["download", "genome", "accession", accession, "--include", "genome",
                              "--filename", str(archive), "--no-progressbar"], "datasets_download.log")
        if res.returncode or not zipfile.is_zipfile(archive):
            raise RuntimeError(f"Genome download failed for {accession}: {res.stderr[-300:]}")
        with zipfile.ZipFile(archive) as z:
            members = [m for m in z.namelist() if m.endswith("_genomic.fna") and f"/{accession}/" in m]
            if not members:
                raise RuntimeError(f"No genomic FASTA in download for {accession}")
            with z.open(members[0]) as src, target.with_suffix(".part").open("wb") as dst:
                shutil.copyfileobj(src, dst)
        target.with_suffix(".part").replace(target)
        archive.unlink()
        return target

    def sixteen(self, accession, genome):
        folder = self.root / "ncbi" / accession
        output = folder / "16s.fasta"
        if output.exists():
            return fasta(output)
        env = dict(os.environ)
        if sys.platform == "darwin":
            # barrnap 0.9 pipes nhmmer --tblout /dev/stdout, which can fail on
            # macOS; redirect that table to a regular file and replay it.
            real = tool("nhmmer")
            compat = folder / "compat"
            compat.mkdir(exist_ok=True)
            wrapper = compat / "nhmmer"
            wrapper.write_text(
                "#!" + sys.executable + "\n"
                "import os, subprocess, sys, tempfile\n"
                "args=sys.argv[1:]\n"
                "with tempfile.TemporaryDirectory() as d:\n"
                "    target=os.path.join(d,'table.txt')\n"
                "    redirect='--tblout' in args and args[args.index('--tblout')+1]=='/dev/stdout'\n"
                "    if redirect: args[args.index('--tblout')+1]=target\n"
                f"    result=subprocess.run([{real!r}]+args)\n"
                "    if redirect and os.path.exists(target):\n"
                "        with open(target) as f: sys.stdout.write(f.read())\n"
                "    sys.exit(result.returncode)\n")
            wrapper.chmod(0o755)
            env["PATH"] = str(compat) + os.pathsep + env["PATH"]
        gff = folder / "rrna.gff"
        with gff.open("w") as out, (folder / "barrnap.log").open("w") as err:
            subprocess.run([tool("barrnap"), "--kingdom", "bac", "--threads", "1", str(genome)],
                           stdout=out, stderr=err, check=True, env=env, stdin=subprocess.DEVNULL)
        seqs = fasta(genome)
        comp = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")
        records = {}
        for line in gff.read_text().splitlines():
            if line.startswith("#") or "16S_rRNA" not in line:
                continue
            f = line.split("\t")
            start, end = int(f[3]), int(f[4])
            s = seqs[f[0]][start - 1:end]
            if f[6] == "-":
                s = s.translate(comp)[::-1]
            records[f"{accession}__copy_{len(records) + 1}"] = s
        write_fasta(output, records)
        return records


def related_names(label):
    """Ordered groups of NCBI genus names to try for a SILVA genus label."""
    groups = [("genus", [label])]
    if label in TAXON_ALIASES:
        groups.append(("alias", TAXON_ALIASES[label]))
    cleaned = re.sub(r"[\[\]]", "", label)
    if "-" in cleaned:  # composite SILVA labels: pool every component genus
        parts = [p for p in cleaned.split("-") if re.fullmatch(r"[A-Z][a-z]+", p)]
        if len(parts) > 1:
            groups.append(("related_genus", parts))
    tokens = [t for t in re.findall(r"[A-Z][a-z]+", cleaned) if t != "Candidatus"]
    if tokens and tokens[0] != label:
        groups.append(("related_genus", [tokens[0]]))
    out, seen = [], set()
    for kind, names in groups:
        key = tuple(names)
        if key not in seen:
            seen.add(key)
            out.append((kind, names))
    return out


def mash_select(work, query, candidates, ncbi):
    """Return (best candidate, distance, ASV) by Mash distance of 16S copies."""
    ref_files, by_file = [], {}
    for c in candidates:
        try:
            g = ncbi.genome(c["accession"])
            rr = ncbi.sixteen(c["accession"], g)
        except Exception as e:
            log("  genome skipped", c["accession"], str(e)[:160])
            continue
        c["genome_path"] = str(g)
        if rr:
            p = work / f"{c['accession']}.fasta"
            write_fasta(p, rr)
            ref_files.append(p)
            by_file[str(p)] = c
    usable = [c for c in candidates if c.get("genome_path")]
    if not usable:
        return None, "", ""
    if not ref_files:
        return usable[0], "", ""  # no 16S annotated: keep the highest-tier genome
    command([tool("mash"), "sketch", "-k", "21", "-s", "1000", "-o", work / "reference"] + ref_files, work / "mash.log")
    command([tool("mash"), "sketch", "-i", "-k", "21", "-s", "1000", "-o", work / "query", query], work / "mash.log")
    res = subprocess.run([tool("mash"), "dist", str(work / "reference.msh"), str(work / "query.msh")],
                         capture_output=True, text=True, check=True)
    hits = [x.split("\t") for x in res.stdout.splitlines() if x.strip()]
    if not hits:
        return usable[0], "", ""
    hits.sort(key=lambda x: (float(x[2]), float(x[3]), x[0], x[1]))
    return by_file[hits[0][0]], hits[0][2], hits[0][1]


def resolve_genus(genus, info, ncbi, family_fallback):
    """info: dict(asvs={id: seq}, family=str)."""
    key = text_hash(genus + "|" + "|".join(sorted(info["asvs"].values())))
    work = ncbi.root / "selection" / f"{safe_name(genus)}__{key}"
    cached = read_json(work / "selection.json")
    if cached and cached.get("status") == "ok" and Path(cached["genome_path"]).exists():
        return cached
    work.mkdir(parents=True, exist_ok=True)
    query = work / "asvs.fasta"
    write_fasta(query, info["asvs"])
    cands, source = [], ""
    for kind, names in related_names(genus):
        pooled = []
        for n in names:
            t = ncbi.taxid(n, "genus")
            if t:
                pooled += ncbi.candidates(t)
        if pooled:
            cands = pooled[: ncbi.max_candidates * max(1, len(names))]
            source = kind if kind == "genus" else f"{kind}:{'/'.join(names)}"
            break
    if not cands and family_fallback and info.get("family"):
        t = ncbi.taxid(info["family"], "family")
        if t:
            cands = ncbi.candidates(t, diversify=True)
            source = f"family_surrogate:{info['family']}"
    row = dict(genus=genus, source=source, candidates=len(cands), status="no_reference",
               accession="", organism="", tier="", mash_distance="", matched_asv="", genome_path="", note="")
    if cands:
        best, dist, asv = mash_select(work, query, cands, ncbi)
        if best:
            row.update(accession=best["accession"], organism=best["organism"], tier=best["tier"],
                       mash_distance=dist, matched_asv=asv, genome_path=best["genome_path"], status="ok",
                       note="" if dist != "" else "no 16S copies annotated; first candidate used")
    write_json(work / "selection.json", row)
    return row


def run(cfg, root, genus_info):
    gcfg = cfg["genomes"]
    ncbi = NCBI(root / cfg["paths"]["genomes"], gcfg)
    rows = []
    for i, (genus, info) in enumerate(sorted(genus_info.items()), 1):
        try:
            row = resolve_genus(genus, info, ncbi, bool(gcfg.get("family_fallback", True)))
        except Exception as e:
            row = dict(genus=genus, status="failed", note=str(e)[:300])
        rows.append(row)
        log(f"Genome {i}/{len(genus_info)}", genus, row.get("status"), row.get("accession", ""), row.get("source", ""))
    write_csv(root / cfg["paths"]["results"] / "genus_references.csv", rows, GENOME_FIELDS)
    return {r["genus"]: r for r in rows if r.get("status") == "ok"}
