"""Step 6: genome-derived U/L, calibration, MiCRM simulation and comparison plot.

* U (min carbon-throughput LP) and L (max output entropy) per genus x resource,
  from the GEMs (estimate.py, network.py, entropy_solver.py).
* Calibration (calibrate.py):
  1. equal uptake budget: every genus has the same total uptake capacity;
  2. uptake capacity s_i of each genus from its typical abundance in this kind of
     environment (MicrobeAtlas prior, prior.py), never from this sample;
  3. nutrient supply rho inferred from this sample's abundances with U*s fixed.
  Genera represented by the same genome share one maintenance rate m.
* Simulation (crm_solver.py): MiCRM in log-biomass with an exact Jacobian, from
  uniform initial biomass to steady state; every genus's final per-capita growth
  rate is reported.
* Evaluation: leave-one-genus-out (each genus predicted with the supply inferred
  without its abundance) and a genome-only baseline (uniform medium, no abundance
  information). Scores go to metadata.json.
* Figure: grouped bars, observed 16S vs calibrated MiCRM.
"""

import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from estimate import Network, estimate
from calibrate import SUPPLY_ONLY, equal_budget, infer_environment, metrics
from crm_solver import CRM, simulate, leading_eigenvalue


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def load_and_estimate(args):
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
    prior_rows = {x["genus"]: x for x in read(source / "prior.csv")}
    missing = [g for g in names if g not in prior_rows]
    if missing:
        raise ValueError("No MicrobeAtlas prior for: " + ", ".join(missing))
    prior = np.array([float(prior_rows[g]["prior"]) for g in names])
    if not np.all(np.isfinite(prior)) or np.any(prior <= 0):
        raise ValueError("Invalid prior abundances")
    return dict(names=names, compounds=compounds, n=n, r=r, abundance=abundance, m=m, k=k,
                u=u, leakage=leakage, loss=loss, resources=resources, source=source, out=out,
                code_hash=code_hash, accessions=[x["accession"] for x in genera],
                prior=prior, prior_match=[prior_rows[g]["match"] for g in names])


def same_genome_same_m(accessions, m):
    """Genera represented by the same genome get the same maintenance rate
    (different m on identical uptake rows would exclude one of them by construction)."""
    m = np.array(m, float)
    groups = {}
    for i, a in enumerate(accessions):
        groups.setdefault(a, []).append(i)
    for idx in groups.values():
        m[idx] = m[idx].mean()
    return m


def run(args):
    defaults = dict(turnover=0.1, biomass=1.0, t_end=300000.0, detection_limit=1e-6,
                    resource_cap=0.2, leave_one_out=True, title=None, environment="")
    for key, value in defaults.items():
        if getattr(args, key, None) is None:
            setattr(args, key, value)
    e = load_and_estimate(args)
    out, names, compounds = e["out"], e["names"], e["compounds"]
    n, r = e["n"], e["r"]
    L, loss, k = e["leakage"], e["loss"], e["k"]
    m = same_genome_same_m(e["accessions"], e["m"])
    U0 = equal_budget(e["u"])                          # 1. equal uptake budget
    abundance = e["abundance"]
    observed = abundance / abundance.sum()
    prior = e["prior"] / e["prior"].sum()
    fed = np.array([x["fed"].lower() == "true" for x in e["resources"]])
    if not fed.any():
        raise ValueError("No supplied resources")
    F, tau, det, cap = float(args.feed_total), float(args.turnover), float(args.detection_limit), float(args.resource_cap)
    C0 = np.full(n, float(args.biomass))
    R0 = fed * args.initial_total / fed.sum()
    feed_uniform = fed * F / fed.sum()
    sim_kw = dict(t_end=float(args.t_end))
    predictions, runs, notes = {}, {}, {}

    def rel(C):
        return C / C.sum()

    # 2. uptake capacity s from the MicrobeAtlas prior (never from this sample)
    from_prior = infer_environment(U0, L, loss, m, k, prior, tau, F, cap=cap)
    s = from_prior["s"]
    U = U0 * s[:, None]
    # 3. nutrient supply rho from this sample, U held fixed
    env = infer_environment(U, L, loss, m, k, observed, tau, F, cap=cap, growth_weight=SUPPLY_ONLY)
    model = CRM(U, L, loss, m, k, env["rho"], tau)
    eigen = leading_eigenvalue(model, env["C"], env["R"])
    sim = simulate(model, C0, R0, **sim_kw)
    predictions["calibrated"] = rel(sim["C"])
    runs["calibrated"] = dict(time=sim["time"], converged=sim["converged"], steps=sim["steps"])
    print("calibrated done", flush=True)

    # baseline: equal-budget genome traits, uniform medium, no abundance information
    try:
        g = simulate(CRM(U0, L, loss, m, k, feed_uniform, tau), C0, R0, **sim_kw)
        predictions["gem_only"] = rel(g["C"])
        runs["gem_only"] = dict(time=g["time"], converged=g["converged"], steps=g["steps"])
    except Exception as ex:
        notes["gem_only"] = str(ex)[:300]

    # leave-one-genus-out: supply inferred without genus j's abundance, j simulated with
    # its prior-derived capacity; out-of-sample for j
    loo = np.full(n, np.nan)
    loo_converged = []
    if args.leave_one_out:
        for j in range(n):
            keep = np.arange(n) != j
            try:
                ej = infer_environment(U[keep], L[keep], loss[keep], m[keep], k[keep], observed[keep],
                                       tau, F, cap=cap, growth_weight=SUPPLY_ONLY)
                sj = simulate(CRM(U, L, loss, m, k, ej["rho"], tau), C0, R0, **sim_kw)
                loo[j] = sj["C"][j] / sj["C"].sum()
                loo_converged.append(bool(sj["converged"]))
            except Exception as ex:
                notes.setdefault("leave_one_out", []).append(f"{names[j]}: {str(ex)[:120]}")
            print(f"loo {j + 1}/{n}", flush=True)
        if np.isfinite(loo).any():
            predictions["leave_one_out"] = np.nan_to_num(loo, nan=0.0)

    # ---------------------------------------------------------------- outputs
    scores = {}
    for key, q in predictions.items():
        if key == "leave_one_out":
            ok = np.isfinite(loo)
            scores[key] = metrics(observed[ok], np.maximum(loo[ok], 0), det, normalize=False)
            for drop in ("jsd", "bray_curtis"):   # per-genus predictions from different runs
                scores[key][drop] = None
        else:
            scores[key] = metrics(observed, q, det)
    final_growth = model.growth(sim["R"])
    break_even = np.where(env["status"] == "corrected", "no", "yes")
    with (out / "abundance_comparison.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["genus", "observed", "predicted", "leave_one_out", "gem_only", "prior", "prior_match",
                    "uptake_capacity_s", "breaks_even_at_inferred_supply", "final_growth_rate", "final_state"])
        for i, gname in enumerate(names):
            w.writerow([gname, observed[i], predictions["calibrated"][i],
                        predictions["leave_one_out"][i] if "leave_one_out" in predictions else "",
                        predictions["gem_only"][i] if "gem_only" in predictions else "",
                        prior[i], e["prior_match"][i], s[i], break_even[i], final_growth[i],
                        "present" if predictions["calibrated"][i] >= det else "excluded"])
    with (out / "environment.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["resource", "fed", "rho_uniform", "rho_inferred", "R_steady"])
        for a, c in enumerate(compounds):
            w.writerow([c, bool(fed[a]), feed_uniform[a], env["rho"][a], env["R"][a]])
    np.savez_compressed(out / "parameters.npz", U=U, U_equal_budget=U0, U_genome=e["u"], s=s, L=L, loss=loss,
                        K=k, m=m, R0=R0, feed=env["rho"], feed_uniform=feed_uniform, prior=prior)
    (out / "simulation.json").write_text(json.dumps(dict(
        runs=runs, final_state_calibrated=sim["C"].tolist() + sim["R"].tolist(),
        leave_one_out_converged=loo_converged), indent=2))
    coverage = float(abundance.sum())
    metadata = dict(
        coverage=coverage,
        jsd=scores["calibrated"]["jsd"],
        metrics=scores,
        genera=n,
        resources=r,
        prior_environment=args.environment,
        prior_matched=int(sum(x != "missing_median" for x in e["prior_match"])),
        uptake_corrections=int((np.abs(s - 1) > 1e-6).sum()),
        uptake_capacity_range=[float(s.min()), float(s.max())],
        genera_not_breaking_even=[names[i] for i in np.flatnonzero(break_even == "no")],
        total_biomass=env["B"],
        supply_violation=env["rho_violation"],
        equilibrium_leading_eigenvalue=eigen,
        equilibrium_stable=bool(eigen < 1e-4 * float(np.max(m))),
        runs=runs,
        notes=notes,
        interpretation=("calibrated: equal uptake budget, uptake capacity s from the MicrobeAtlas prior, "
                        "supply rho inferred from this sample; leave_one_out: each genus predicted with the "
                        "supply inferred without its abundance; gem_only: equal-budget genome traits, uniform "
                        "medium, no abundances"),
        code_hash=e["code_hash"],
        arguments={kk: str(v) for kk, v in vars(args).items()},
    )
    metadata["input_hashes"] = {
        name: hashlib.sha256((e["source"] / name).read_bytes()).hexdigest()
        for name in ["genera.csv", "resources.csv", "kinetics.csv", "observed.csv", "maintenance.csv", "prior.csv"]
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str))
    plot(out / "comparison.png", names, observed, predictions, scores, det,
         f"{getattr(args, 'title', None) or e['source'].name}  |  {n} modelled genera  |  "
         f"JSD {scores['calibrated']['jsd']:.3f}", getattr(args, "plot_scale", None) or "linear")
    print(json.dumps(dict(jsd=metadata["jsd"], metrics=scores), default=str), flush=True)


STYLE = {"observed": "#2a78d6", "calibrated": "#e34948"}
LABEL = {"observed": "Observed 16S", "calibrated": "Calibrated MiCRM"}


def plot(path, names, observed, predictions, scores, det, title, scale="linear"):
    """Grouped bar chart, one pair of bars per modelled genus (sorted by
    observed): blue = observed 16S, red = calibrated MiCRM.  Linear y-axis by
    default (evenly spaced ticks); scale="log" is still available."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, muted, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": muted, "axes.labelcolor": ink,
                         "xtick.color": muted, "ytick.color": muted,
                         "axes.spines.top": False, "axes.spines.right": False})
    names = list(names)
    observed = np.asarray(observed, float)
    calibrated = np.asarray(predictions["calibrated"], float)
    n = len(names)
    order = np.argsort(-observed)
    x = np.arange(n)
    width = 0.4
    fig, ax = plt.subplots(figsize=(max(11, n * 0.36), 6), facecolor=surface)
    ax.set_facecolor(surface)
    ax.bar(x - width / 2, observed[order], width, color=STYLE["observed"], edgecolor=surface,
           linewidth=0.6, label=LABEL["observed"], zorder=2)
    ax.bar(x + width / 2, calibrated[order], width, color=STYLE["calibrated"], edgecolor=surface,
           linewidth=0.6, label=LABEL["calibrated"], zorder=2)
    if scale == "log":
        ax.set_yscale("log")
        ax.set_ylim(det, min(1.0, float(max(observed.max(), calibrated.max()))) * 2.5)
        ax.set_ylabel("Relative abundance (log scale)")
    else:
        ax.set_ylim(0, float(max(observed.max(), calibrated.max())) * 1.08)
        ax.set_ylabel("Relative abundance within modelled genera")
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_xticks(x, [names[i] for i in order], rotation=90, fontsize=8)
    ax.grid(axis="y", color=grid, lw=0.6, zorder=0)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title(title, loc="left", fontsize=11, color=ink)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor=surface)
    plt.close(fig)


def replot(out, title=None, scale="linear"):
    """Redraw comparison.png from a finished 03_micrm folder (no re-simulation)."""
    out = Path(out)
    rows = read(out / "abundance_comparison.csv")
    meta = json.loads((out / "metadata.json").read_text())
    names = [r["genus"] for r in rows]
    observed = np.array([float(r["observed"]) for r in rows])
    predictions = {"calibrated": np.array([float(r["predicted"]) for r in rows])}
    det = float(meta.get("arguments", {}).get("detection_limit", 1e-6))
    title = title or out.parent.name
    plot(out / "comparison.png", names, observed, predictions, meta.get("metrics", {}), det,
         f"{title}  |  {len(names)} modelled genera  |  JSD {meta['jsd']:.3f}", scale)
    return out / "comparison.png"



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--initial-total", type=float, default=0.0001)
    parser.add_argument("--feed-total", type=float, required=True)
    parser.add_argument("--no-loo", dest="leave_one_out", action="store_false")
    args = parser.parse_args()
    if (not np.isfinite(args.initial_total) or not np.isfinite(args.feed_total)
            or args.initial_total < 0 or args.feed_total <= 0):
        parser.error("Resource totals must be finite and nonnegative, with positive feed")
    run(args)
