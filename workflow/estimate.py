"""Free precursor pool and feasibility-selected retention, without growth maximization."""

import numpy as np
from scipy.sparse import coo_matrix, hstack, vstack, eye
from scipy.optimize import linprog
from network import Network
from entropy_solver import solve_entropy


def estimate(net, resource):
    if resource not in net.exchanges or not net.transport_supported(resource):
        return {"status": "no_import"}
    prec = net.precursor_indices()
    if not prec:
        return {"status": "no_precursors"}
    j, k, sign = net.exchanges[resource]
    carbon = net.carbon[k]
    if not np.isfinite(carbon) or carbon <= 0:
        return {"status": "unknown_carbon"}
    lo, hi = net.bounds(strict=True)
    lo[j] = hi[j] = 1 / (sign * carbon)
    lo[net.bio] = hi[net.bio] = 0
    nr = len(lo)
    np_ = len(prec)
    drains = coo_matrix(
        (-np.ones(np_), (prec, np.arange(np_))), shape=(len(net.sid), np_)
    ).tocsr()[net.balance_rows]
    A = hstack([net.A, drains], format="csr")
    nv = nr + np_
    outputs = [
        (r, jj, kk, ss)
        for r, (jj, kk, ss) in net.exchanges.items()
        if r != resource
        and np.isfinite(net.carbon[kk])
        and net.carbon[kk] > 0
        and net.formula[kk] is not None
        and not net.formula[kk].get("R", 0)
        and not net.formula[kk].get("X", 0)
    ]
    if not outputs:
        return {"status": "no_outputs"}
    oi = [x[1] for x in outputs]
    w = np.array([-x[3] * net.carbon[x[2]] for x in outputs])
    release = coo_matrix((w, (np.zeros(len(oi)), oi)), shape=(1, nv)).tocsr()
    retention = coo_matrix(
        (net.carbon[prec], (np.zeros(np_), nr + np.arange(np_))), shape=(1, nv)
    ).tocsr()
    eq = vstack([A, release, retention], format="csr")
    lower = np.r_[lo, np.zeros(np_)]
    upper = np.r_[hi, 1 / net.carbon[prec]]
    attempts = []
    lp = None
    # Prior-centred feasibility search; no yield or biomass maximization.
    for retained in [0.75, 0.5, 0.25, 0.1, 0.01, 0.001, 0.9, 0.99]:
        lam = 1 - retained
        rhs = np.r_[np.zeros(A.shape[0]), lam, retained]
        lp = linprog(
            np.zeros(nv),
            A_eq=eq,
            b_eq=rhs,
            bounds=list(zip(lower, upper)),
            method="highs",
        )
        attempts.append({"retained": retained, "status": int(lp.status)})
        if lp.success:
            break
    if not lp.success:
        return {"status": "no_feasible_retention", "attempts": attempts}
    cost = np.zeros(nv)
    for z in range(nr):
        if z == net.bio or z in net.boundaries:
            continue
        col = net.S[:, z]
        c = net.carbon[col.indices]
        if np.all(np.isfinite(c)):
            cost[z] = 0.5 * np.sum(np.abs(col.data) * c)
    # Minimize absolute carbon throughput for the selected retention constraint.
    zeros = coo_matrix(eq.shape)
    E = hstack([eq, zeros], format="csr")
    I = eye(nv, format="csr")
    inequalities = vstack([hstack([I, -I]), hstack([-I, -I])], format="csr")
    minimum = linprog(
        np.r_[np.zeros(nv), cost],
        A_ub=inequalities,
        b_ub=np.zeros(2 * nv),
        A_eq=E,
        b_eq=rhs,
        bounds=list(zip(lower, upper)) + [(0, None)] * nv,
        method="highs",
    )
    if not minimum.success:
        return {"status": "cost_failed", "attempts": attempts}
    entropy = solve_entropy(
        eq,
        rhs,
        lower,
        upper,
        oi,
        w,
        lam,
        [minimum.x[:nv]],
        tolerance=1e-5,
        max_iterations=100,
    )
    if entropy["status"] == "optimal":
        flux = entropy["solution"]
        mode = "maximum_entropy"
        gap = entropy["gap"]
    else:
        # Explicit feasible cost witness when entropy refinement is unresolved.
        flux = minimum.x[:nv]
        mode = "minimum_cost_feasible_fallback"
        gap = entropy.get("gap")
    release_values = w * flux[oi]
    residual = max(
        float(np.max(np.abs(eq @ flux - rhs))),
        float(np.max(lower - flux)),
        float(np.max(flux - upper)),
        max(0.0, -float(release_values.min())),
    )
    if residual > 1e-6:
        return {"status": "witness_failed", "residual": residual}
    return {
        "status": "ok",
        "u": 10 / (1 + float(minimum.fun)),
        "cost": float(minimum.fun),
        "retained": retained,
        "lambda": lam,
        "residual": residual,
        "allocation": mode,
        "entropy_status": entropy["status"],
        "entropy_gap": gap,
        "attempts": attempts,
        "outputs": {
            r: float(max(0.0, q))
            for (r, *_), q in zip(outputs, release_values)
            if q > 1e-10
        },
    }
