"""Step 6: chemical-potential U/L estimation, abundance calibration, BDF MiCRM and\ncomparison plot.  Copied from micrm_shared/run.py; only parameters made configurable."""

import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from estimate import Network, estimate
from bdf_solver import integrate


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def run(args):
    """args: namespace with input, catalog, cache, out, initial_total, feed_total
    and optional turnover, biomass, t_end, title."""
    source = Path(args.input).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    genera = read(source / "genera.csv")
    resources = read(source / "resources.csv")
    names = [x["genus"] for x in genera]
    compounds = [x["resource"] for x in resources]
    if len(names) != len(set(names)) or len(compounds) != len(set(compounds)):
        raise ValueError("Duplicate genus or resource")
    n, r = len(names), len(compounds)
    if n < 2 or r == 0:
        raise ValueError("Require two genera and at least one resource")
    gi = {g: i for i, g in enumerate(names)}
    ri = {c: a for a, c in enumerate(compounds)}
    observations = read(source / "observed.csv")
    counts = {x["genus"]: float(x["count"]) for x in observations}
    if (
        any(not np.isfinite(v) or v < 0 for v in counts.values())
        or sum(counts.values()) <= 0
    ):
        raise ValueError("Invalid observations")
    abundance = np.array([counts.get(g, 0) for g in names]) / sum(counts.values())
    if np.any(abundance <= 0):
        raise ValueError("Selected genera require positive counts")
    maintenance = {x["genus"]: float(x["m"]) for x in read(source / "maintenance.csv")}
    m = np.array([maintenance[g] for g in names])
    k = np.full((n, r), np.inf)
    if not np.all(np.isfinite(m)) or np.any(m <= 0):
        raise ValueError("Invalid m")
    for row in read(source / "kinetics.csv"):
        if row["genus"] in gi and row["resource"] in ri:
            value = float(row["K"])
            if not np.isfinite(value) or value <= 0:
                raise ValueError("Invalid K")
            k[gi[row["genus"]], ri[row["resource"]]] = value
    catalog = {x["accession"]: x for x in read(args.catalog)}
    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    code_hash = hashlib.sha256(
        b"".join(
            (Path(__file__).parent / x).read_bytes()
            for x in ["estimate.py", "network.py", "entropy_solver.py"]
        )
    ).hexdigest()
    u = np.zeros((n, r))
    leakage = np.zeros((n, r, r))
    loss = np.zeros((n, r))
    audit = []
    for i, genus in enumerate(genera):
        xml = Path(catalog[genus["accession"]]["xml_path"])
        if not xml.is_absolute():
            xml = Path(args.catalog).resolve().parent / xml
        digest = hashlib.sha256(xml.read_bytes()).hexdigest()
        net = None
        for a, compound in enumerate(compounds):
            if not np.isfinite(k[i, a]):
                continue
            key = hashlib.sha256(
                f"{digest}:{code_hash}:{compound}".encode()
            ).hexdigest()
            file = cache / (key + ".json")
            if file.exists():
                result = json.loads(file.read_text())
            else:
                if net is None:
                    net = Network(xml)
                result = estimate(net, compound)
                import tempfile, os

                with tempfile.NamedTemporaryFile(
                    mode="w", dir=cache, delete=False
                ) as f:
                    json.dump(result, f, allow_nan=False)
                    temporary = f.name
                os.replace(temporary, file)
            audit.append(dict(genus=names[i], resource=compound, **result))
            if result["status"] != "ok":
                continue
            u[i, a] = result["u"]
            for product, fraction in result["outputs"].items():
                if product in ri:
                    leakage[i, a, ri[product]] = fraction
                else:
                    loss[i, a] += fraction
    (out / "parameter_audit.json").write_text(json.dumps(audit, indent=2))
    # binary leakage capability (substrate -> released product), for inspection
    with (out / "binary_leakage.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["genus", "resource", "product", "in_panel", "base_fraction"])
        for a in audit:
            for product, fraction in (a.get("outputs") or {}).items():
                w.writerow([a["genus"], a["resource"], product, int(product in ri), fraction])
    if not np.any(u > 0):
        raise ValueError("No feasible uptake pairs")
    base_u = u.copy()
    preference = base_u / np.maximum(base_u.sum(axis=1, keepdims=True), 1e-300)
    demand = (preference * abundance[:, None]).sum(axis=0)
    weights = 0.1 / r + 0.9 * demand / demand.sum()
    organic = leakage.sum(axis=2)
    weighted = leakage * weights[None, None, :]
    scale = np.divide(
        organic,
        weighted.sum(axis=2),
        out=np.zeros_like(organic),
        where=weighted.sum(axis=2) > 0,
    )
    leakage = weighted * scale[:, :, None]
    u *= (abundance / abundance.mean())[:, None]
    if not np.allclose(leakage.sum(axis=2), organic, atol=1e-10):
        raise ValueError("Organic carbon mismatch")
    if np.any(leakage.sum(axis=2) + loss > 1 + 1e-8):
        raise ValueError("Carbon exceeds input")
    fed = np.array([x["fed"].lower() == "true" for x in resources])
    if not fed.any():
        raise ValueError("No supplied resources")
    initial = fed.astype(float) * args.initial_total / fed.sum()
    feed = fed.astype(float) * args.feed_total / fed.sum()
    params = dict(
        n=n,
        r=r,
        U=u.ravel(order="F").tolist(),
        L=leakage.ravel(order="F").tolist(),
        loss=loss.ravel(order="F").tolist(),
        m=m.tolist(),
        feed=feed.tolist(),
        initial=initial.tolist(),
        K=k.ravel(order="F").tolist(),
        turnover=float(getattr(args, "turnover", 0.1)),
        biomass=float(getattr(args, "biomass", 1.0)),
        t_end=float(getattr(args, "t_end", 300000.0)),
        atol=1e-8,
        rtol=1e-5,
        tolerance=1e-5,
    )
    result = integrate(params)
    biomass = np.asarray(result["state"][:n])
    if biomass.sum() <= 1e-12:
        raise ValueError("All modeled biomass extinct")
    predicted = biomass / biomass.sum()
    observed = abundance / abundance.sum()
    midpoint = (predicted + observed) / 2

    def divergence(x):
        mask = x > 0
        return float(np.sum(x[mask] * np.log(x[mask] / midpoint[mask])))

    jsd = 0.5 * (divergence(predicted) + divergence(observed))
    with (out / "abundance_comparison.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["genus", "observed", "predicted"])
        w.writerows(zip(names, observed, predicted))
    np.savez_compressed(
        out / "parameters.npz",
        U=u,
        L=leakage,
        loss=loss,
        K=k,
        m=m,
        base_U=base_u,
        R0=initial,
        feed=feed,
    )
    (out / "simulation.json").write_text(json.dumps(result))
    metadata = dict(
        coverage=float(abundance.sum()),
        jsd=jsd,
        time=result["time"],
        early_stopped=result["early_stopped"],
        leakage="abundance-reweighted; carbon totals preserved; full stoichiometric feasibility not re-certified",
        interpretation="within-sample calibration",
        code_hash=code_hash,
        arguments={k: str(v) for k, v in vars(args).items()},
    )
    metadata["input_hashes"] = {
        name: hashlib.sha256((source / name).read_bytes()).hexdigest()
        for name in [
            "genera.csv",
            "resources.csv",
            "kinetics.csv",
            "observed.csv",
            "maintenance.csv",
        ]
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = np.argsort(-observed)
    x = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(11, n * 0.14), 6))
    ax.bar(x - 0.2, observed[order], 0.4, label="observed 16S", color="#00BFC4")
    ax.bar(x + 0.2, predicted[order], 0.4, label="MiCRM", color="#F8766D")
    ax.set_xticks(x, [names[i] for i in order], rotation=90, fontsize=7)
    ax.set_ylabel("Relative abundance within modeled genera")
    ax.set_title(f"{getattr(args, 'title', None) or source.name} | coverage {abundance.sum():.1%} | JSD {jsd:.3f}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)
    print(json.dumps(metadata), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--initial-total", type=float, default=0.0001)
    parser.add_argument("--feed-total", type=float, required=True)
    args = parser.parse_args()
    if (
        not np.isfinite(args.initial_total)
        or not np.isfinite(args.feed_total)
        or args.initial_total < 0
        or args.feed_total <= 0
    ):
        parser.error(
            "Resource totals must be finite and nonnegative, with positive feed"
        )
    run(args)
