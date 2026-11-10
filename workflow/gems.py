"""Step 4: build one gapseq GEM per selected reference genome.

Stages (from the original vendor/run_pipeline.py):
  gapseq find -> find-transport -> draft -> medium -> fill
Each stage is checkpointed in gems/<accession>/, so an interrupted run resumes
and an existing GEM built by this workflow is reused.  Nothing is imported from
outside this folder.
"""

import concurrent.futures
import csv
import gzip
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

from common import command, digest, log, read_json, tool, write_csv, write_json
from network import Network
from ul_core import parse_sbml

CATALOG_FIELDS = ["genus", "accession", "xml_path", "status", "reference_source", "error"]


def tbl(path):
    name = str(path)
    opener = gzip.open(path, "rt", encoding="latin-1") if name.endswith(".gz") else Path(path).open(encoding="latin-1")
    delimiter = "," if name.removesuffix(".gz").endswith(".csv") else "\t"  # gapseq medium.csv is comma-separated
    with opener as f:
        return list(csv.DictReader((line for line in f if not line.startswith("#")), delimiter=delimiter))


def checked_table(path, required):
    if not Path(path).is_file() or not Path(path).stat().st_size:
        raise ValueError(f"Missing stage output: {path}")
    rows = tbl(path)
    if required and not set(required) <= set(rows[0].keys() if rows else []):
        raise ValueError(f"Invalid table header in {path}; expected {required}")
    return rows


def run_stage(folder, name, args, products, timeout):
    checkpoint = folder / (name + ".complete.json")
    signature = hashlib.sha256(json.dumps([str(a) for a in args]).encode()).hexdigest()
    old = read_json(checkpoint)
    if old and old.get("signature") == signature and all(
        (folder / f).is_file() and digest(folder / f) == sha for f, sha in old["outputs"].items()
    ):
        return
    with tempfile.TemporaryDirectory(prefix=name + "_", dir=folder) as temporary:
        work = Path(temporary)
        command(args, folder / (name + ".log"), work, timeout=timeout)
        for filename, columns in products.items():
            if filename.endswith(".xml"):
                species, reactions = parse_sbml(clean_copy(work / filename))
                if not species or not reactions:
                    raise ValueError("Empty metabolic model")
            elif filename.endswith((".tbl", ".csv")):
                checked_table(work / filename, columns)
            elif not (work / filename).is_file():
                raise ValueError(f"Missing stage output: {filename}")
        for filename in products:
            shutil.copy2(work / filename, folder / filename)
        for path in work.glob("genome*"):
            if path.is_file():
                shutil.copy2(path, folder / path.name)
    write_json(checkpoint, dict(signature=signature, outputs={f: digest(folder / f) for f in products}))


def clean_copy(xml):
    """gapseq occasionally writes non-UTF-8 bytes; parse a sanitised copy."""
    xml = Path(xml)
    data = xml.read_bytes()
    try:
        data.decode("utf-8")
        return xml
    except UnicodeDecodeError:
        out = xml.with_suffix(".clean.xml")
        out.write_bytes(data.decode("utf-8", "replace").encode())
        return out


def build(accession, genome, gem_root, threads, database, timeout):
    folder = Path(gem_root) / accession
    folder.mkdir(parents=True, exist_ok=True)
    final = folder / "model.xml"
    status = read_json(folder / "status.json", {})
    if status.get("status") == "ok" and final.exists() and status.get("genome_sha256") == digest(genome):
        return dict(status="ok", xml=str(final))
    local = folder / "genome.fna"
    if not local.exists() or digest(local) != digest(genome):
        shutil.copyfile(genome, local)
    gap = tool("gapseq")
    db = ["-D", database] if database else []
    p = "genome"
    run_stage(folder, "find",
              [gap, "find", "-p", "all", "-u", "all", "-b", "12", "-K", threads, "-A", "diamond"] + db + [local],
              {p + "-all-Reactions.tbl": ["rxn", "dbhit", "status"], p + "-all-Pathways.tbl": ["ID", "Prediction"]}, timeout)
    run_stage(folder, "find_transport",
              [gap, "find-transport", "-b", "12", "-K", threads, "-A", "diamond", local],
              {p + "-Transporter.tbl": ["exid", "rea"]}, timeout)
    run_stage(folder, "draft",
              [gap, "draft", "-r", folder / (p + "-all-Reactions.tbl"), "-t", folder / (p + "-Transporter.tbl"),
               "-p", folder / (p + "-all-Pathways.tbl")],
              {p + "-draft.xml": [], p + "-draft.RDS": []}, timeout)
    # A draft is not biomass-competent: predict a medium and gap-fill.
    run_stage(folder, "medium",
              [gap, "medium", "-m", folder / (p + "-draft.RDS"), "-p", folder / (p + "-all-Pathways.tbl"), "-o", "medium.csv"],
              {"medium.csv": ["compounds", "name", "maxFlux"]}, timeout)
    run_stage(folder, "fill",
              [gap, "fill", "-m", folder / (p + "-draft.RDS"), "-n", folder / "medium.csv", "-f", ".", "-z", "glpk"],
              {p + ".xml": []}, timeout)
    xml = clean_copy(folder / (p + ".xml"))
    species, _ = parse_sbml(xml)
    if not {"c0", "e0"} <= {s["compartment"] for s in species.values()}:
        raise ValueError("Unsupported model compartments: c0/e0 required")
    Network(xml)  # must expose a biomass reaction for the parameterisation
    shutil.copy2(xml, final)
    write_json(folder / "status.json", dict(status="ok", accession=accession, genome_sha256=digest(genome),
                                            xml_sha256=digest(final), gapseq=str(Path(gap).resolve())))
    return dict(status="ok", xml=str(final))


def run(cfg, root, references):
    """references: genus -> selection row (accession, genome_path, source)."""
    gem_root = root / cfg["paths"]["gems"]
    gem_root.mkdir(parents=True, exist_ok=True)
    g = cfg["gems"]
    database = g.get("gapseq_db", "")
    timeout = float(g["timeout_hours"]) * 3600 if g.get("timeout_hours") else None
    workers = max(1, int(cfg["compute"]["parallel_gems"]))
    threads = max(1, int(cfg["compute"]["threads"]) // workers)
    accessions = {r["accession"]: r["genome_path"] for r in references.values()}
    built = {}

    def job(acc):
        try:
            return acc, build(acc, accessions[acc], gem_root, threads, database, timeout)
        except Exception as e:
            write_json(gem_root / acc / "status.json", dict(status="failed", error=str(e)[:2000]))
            return acc, dict(status="failed", error=str(e)[:300])

    with concurrent.futures.ThreadPoolExecutor(workers) as pool:
        for n, (acc, res) in enumerate(pool.map(job, sorted(accessions)), 1):
            built[acc] = res
            log(f"GEM {n}/{len(accessions)}", acc, res["status"], res.get("error", "")[:160])
    rows = []
    for genus, ref in sorted(references.items()):
        res = built[ref["accession"]]
        rows.append(dict(genus=genus, accession=ref["accession"], xml_path=res.get("xml", ""),
                         status=res["status"], reference_source=ref.get("source", ""), error=res.get("error", "")))
    catalog = root / cfg["paths"]["results"] / "gem_catalog.csv"
    write_csv(catalog, rows, CATALOG_FIELDS)
    return catalog, [r for r in rows if r["status"] == "ok"]
