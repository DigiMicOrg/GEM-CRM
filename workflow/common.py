"""Shared helpers: CSV/JSON IO, hashing, logging and external command execution."""

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

csv.field_size_limit(10**9)


def log(*parts):
    stamp = time.strftime("%H:%M:%S")
    print(f"[{stamp}]", *parts, flush=True)


def read_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def read_json(path, default=None):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else default


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str))
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def text_hash(text, n=16):
    return hashlib.sha256(text.encode()).hexdigest()[:n]


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "x"


def tool(name):
    value = shutil.which(name)
    if not value:
        raise RuntimeError(
            f"Missing executable: {name}. Activate the conda environment from environment.yml."
        )
    return value


def child_environment():
    # gapseq R scripts interpret even RSTUDIO=0 as an interactive IDE.
    # Clean only the child environment; do not change the user's R session.
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("RSTUDIO") or key in {
            "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_SESSION_TMPDIR"
        }:
            env.pop(key, None)
    for key in ("R_ENVIRON", "R_ENVIRON_USER", "R_PROFILE", "R_PROFILE_USER"):
        env[key] = os.devnull
    return env


def command(args, log_path, cwd=None, timeout=None, check_output_text=True):
    """Run a CLI, append its output to log_path, and fail on hidden R errors."""
    args = [str(x) for x in args]
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write("\n" + json.dumps(args) + "\n")
        f.flush()
        start = f.tell()
        try:
            result = subprocess.run(
                args, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, env=child_environment(), timeout=timeout,
            )
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = "timeout"
    with log_path.open() as f:
        f.seek(start)
        output = f.read()
    inspected = output
    # gapseq 2.1 probes Zenodo even with a complete custom database and then
    # continues; only inspect what follows that message.
    custom_db = "Using custom database directory:"
    if custom_db in output:
        inspected = output.split(custom_db, 1)[1]
    bad = check_output_text and re.search(
        r"Execution halted|^Error in |illegal option|^Error:", inspected, re.M
    )
    if code != 0 or bad:
        tail = "\n".join(output.splitlines()[-18:])
        raise RuntimeError(
            f"{Path(args[0]).name} {args[1] if len(args) > 1 else ''} failed "
            f"(exit={code}). Log: {log_path}\n{tail}"
        )
    return output


def fasta(path):
    records = {}
    name, seq = None, []
    with Path(path).open() as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        records[name] = "".join(seq).upper()
    return records


def write_fasta(path, records):
    Path(path).write_text("".join(f">{k}\n{v}\n" for k, v in records.items()))


def python_executable():
    return sys.executable
