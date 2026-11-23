"""Maximum output entropy on a linear flux polytope, with an LP gap certificate.
Fully corrective conditional gradients add feasible flux vertices on demand;
this neither enumerates metabolic pathways nor maximizes growth.
"""

import numpy as np
from scipy.optimize import linprog, minimize


def solve_entropy(
    eq,
    rhs,
    lower,
    upper,
    indices,
    weights,
    lam,
    witnesses,
    tolerance=1e-6,
    max_iterations=300,
):
    indices = np.asarray(indices)
    weights = np.asarray(weights)
    vertices = [np.asarray(v).copy() for v in witnesses]
    theta = np.ones(len(vertices)) / len(vertices)
    bounds = list(zip(lower, upper))
    last_gap = np.inf
    for iteration in range(max_iterations):
        V = np.asarray(vertices)
        P = np.maximum(V[:, indices] * weights / lam, 0)

        def objective(t):
            p = t @ P
            return float(np.sum(p * np.log(np.maximum(p, 1e-16))))

        def jac(t):
            return P @ (np.log(np.maximum(t @ P, 1e-16)) + 1)

        correction = minimize(
            objective,
            theta,
            jac=jac,
            method="SLSQP",
            bounds=[(0, 1)] * len(theta),
            constraints=[
                {
                    "type": "eq",
                    "fun": lambda t: sum(t) - 1,
                    "jac": lambda t: np.ones(len(t)),
                }
            ],
            options={"ftol": 1e-12, "maxiter": 500},
        )
        if (
            np.all(np.isfinite(correction.x))
            and abs(sum(correction.x) - 1) < 1e-7
            and min(correction.x) >= -1e-9
        ):
            theta = np.maximum(correction.x, 0)
            theta /= sum(theta)
        current = theta @ V
        p = np.maximum(current[indices] * weights / lam, 1e-16)
        grad = -(np.log(p) + 1) / lam
        cost = np.zeros(len(lower))
        cost[indices] = -grad * weights
        oracle = linprog(cost, A_eq=eq, b_eq=rhs, bounds=bounds, method="highs")
        if not oracle.success:
            return {"status": "entropy_oracle_failed", "reason": oracle.message}
        gap = float(grad @ ((oracle.x[indices] - current[indices]) * weights))
        last_gap = gap
        if gap <= tolerance:
            return {
                "status": "optimal",
                "solution": current,
                "gap": max(0, gap),
                "iterations": iteration + 1,
            }
        # A new vertex may allow a better simultaneous allocation. Never equate
        # per-product individual maxima with one achievable allocation.
        vertices.append(oracle.x)
        theta = np.r_[theta, 0.0]
    return {
        "status": "entropy_gap_not_converged",
        "gap": last_gap,
        "iterations": max_iterations,
    }
