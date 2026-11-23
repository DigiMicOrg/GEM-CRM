"""Calibration of the MiCRM: uptake capacities from an external prior, nutrient
supply from the sample.

Assumptions and steps (see micrm.run):

1. Equal uptake budget (equal_budget): every genus gets the same total uptake
   capacity, spread over the carbon sources it can import. Being able to use
   many substrates makes a genus broad, not a bigger eater. Without this, GEM
   completeness (number of importable compounds) decides who wins.

2. Uptake capacity s_i from the MicrobeAtlas prior (infer_environment with the
   prior abundances): the smallest per-genus scaling of U such that the genus's
   *typical* abundance in this environment class (not this sample's) is a
   steady state. This is the genus's intrinsic competitiveness.

3. Nutrient supply rho from the sample (infer_environment with the observed
   abundances and growth_weight=SUPPLY_ONLY): U*s is held fixed; an LP finds
   resource levels at which every genus breaks even and the supply that makes
   the observed abundances a steady state, choosing the supply closest to a
   uniform medium. Genera that cannot break even jointly are reported, not
   corrected.

Steady state with every genus present:

  growth:  sum_a U_ia eps_ia x_a = m_i            x_a = R_a / (1 + R_a/K_a)
  supply:  rho_a = R_a + (B/tau) (P x)_a >= 0     P = net uptake per biomass

Why abundances are never put into U directly: in a consumer-resource model the
winner of a resource is the genus with the lowest R* (~ m / U), so scaling U by
abundance turns abundance into competitive rank and drives rare genera extinct.
"""

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack, identity
from scipy.stats import spearmanr


def _column_K(K):
    K = np.asarray(K, float)
    fin = np.where(np.isfinite(K) & (K > 0), K, np.inf)
    return fin.min(axis=0)


SUPPLY_ONLY = 1e8   # growth-slack weight that makes break-even (practically) a hard constraint


def equal_budget(U):
    """Rescale each genus's uptake row to the same total (the mean row total)."""
    U = np.asarray(U, float)
    rows = U.sum(axis=1)
    return U * np.where(rows > 0, rows.mean() / np.where(rows > 0, rows, 1), 0)[:, None]


def _lp(A, m, Q, xmax, f, total, growth_weight=1e4):
    """min  W_g sum|e|/m + W_v sum v/f + sum|Qx - f|/total
    s.t. A x + e+ - e- = m;  Q x + v >= 0;  0 <= x <= xmax."""
    n, r = A.shape
    nv = r + 2 * n + 2 * r  # x, e+, e-, v, t
    sx, sp, sm, sv, st = 0, r, r + n, r + 2 * n, 2 * r + 2 * n
    c = np.zeros(nv)
    c[sp:sm] = growth_weight / m
    c[sm:sv] = growth_weight / m
    c[sv:st] = 1e4 / f
    c[st:] = 1.0 / total
    In, Ir = identity(n, format="csr"), identity(r, format="csr")
    Zr = csr_matrix((r, r))
    Zrn = csr_matrix((r, n))
    Qs = csr_matrix(Q)
    A_eq = hstack([csr_matrix(A), In, -In, csr_matrix((n, r)), csr_matrix((n, r))], format="csr")
    A_ub = vstack([
        hstack([-Qs, Zrn, Zrn, -Ir, Zr]),   # -(Qx) - v <= 0
        hstack([Qs, Zrn, Zrn, Zr, -Ir]),    #  Qx - t <= f
        hstack([-Qs, Zrn, Zrn, Zr, -Ir]),   # -Qx - t <= -f
    ], format="csr")
    b_ub = np.r_[np.zeros(r), np.full(r, f), np.full(r, -f)]
    bounds = [(0, xm) for xm in xmax] + [(0, None)] * (nv - r)
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=m, bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError("Environment LP failed: " + res.message)
    return res.x[:r]


def infer_environment(U, L, loss, m, K, p, tau, feed_total, cap=0.2, iterations=15, growth_weight=1e4):
    """Return supply rho, resource levels R*, total biomass B and the per-genus
    uptake scaling s needed for p to be a steady state (1 = U already feasible).
    With growth_weight=SUPPLY_ONLY the LP keeps U fixed wherever possible and s
    only flags genera that cannot break even (callers must not apply it)."""
    U, L, loss, m, p = (np.asarray(a, float) for a in (U, L, loss, m, p))
    p = p / p.sum()
    n, r = U.shape
    eps = 1.0 - L.sum(axis=2) - loss
    A = U * eps
    Kc = _column_K(K)
    invK = np.where(np.isfinite(Kc), 1.0 / np.where(np.isfinite(Kc), Kc, 1.0), 0.0)
    f = feed_total / r
    # R <= cap*K keeps uptake in the unsaturated regime (x <= cap K / (1+cap))
    xmax = np.where(np.isfinite(Kc), cap * Kc / (1 + cap), feed_total)
    pu = p[:, None] * U
    P = np.diag(pu.sum(axis=0)) - np.einsum("ib,iba->ab", pu, L)
    B = tau * feed_total / max(float(p @ m), 1e-12)
    history = []
    for _ in range(iterations):
        Q = np.eye(r) + (B / tau) * P   # Q x ~ rho (uses R >= x)
        x = _lp(A, m, Q, xmax, f, feed_total, growth_weight)
        R = x / np.maximum(1.0 - x * invK, 1e-12)
        net = float((P @ x).sum())
        history.append(B)
        if net <= 0 or R.sum() >= feed_total:
            break
        B_new = tau * (feed_total - R.sum()) / net
        converged = abs(B_new - B) <= 1e-4 * B
        B = B_new
        if converged:
            break
    Q = np.eye(r) + (B / tau) * P
    x = _lp(A, m, Q, xmax, f, feed_total, growth_weight)
    R = x / np.maximum(1.0 - x * invK, 1e-12)
    rho_raw = R + (B / tau) * (P @ x)
    rho = np.maximum(rho_raw, 0.0)
    # Resources nobody in this set can take up are unconstrained by the data;
    # give them a typical supply rather than the uniform default, so a genus
    # left out (leave-one-out) is not handed a private niche by construction.
    idle = U.sum(axis=0) <= 0
    if idle.any() and (~idle).any():
        rho[idle] = float(np.median(rho[~idle]))
    supplied = A @ x
    s = np.where(supplied > 1e-15, m / np.maximum(supplied, 1e-300), np.nan)
    status = np.where(np.isnan(s), "no_usable_uptake",
                      np.where(np.abs(s - 1) <= 1e-6, "supported", "corrected"))
    s = np.where(np.isnan(s), 1.0, s)
    s[np.abs(s - 1) <= 1e-6] = 1.0
    return dict(x=x, R=R, rho=rho, rho_violation=float(np.maximum(-rho_raw, 0).sum()),
                B=float(B), C=B * p, s=s, status=status, B_history=history)


def metrics(observed, predicted, detection=1e-6, rare=0.01, normalize=True):
    o = np.asarray(observed, float)
    o = o / o.sum()
    q = np.asarray(predicted, float)
    if normalize and q.sum() > 0:
        q = q / q.sum()
    mid = (o + q) / 2

    def kl(a):
        k = a > 0
        return float(np.sum(a[k] * np.log(a[k] / mid[k])))

    kept = q >= detection
    rare_mask = o < rare
    lo, lq = np.log10(o), np.log10(np.maximum(q, detection))
    rho_s = spearmanr(lo, lq).statistic if len(o) > 2 and np.ptp(lq) > 0 else float("nan")
    return dict(
        jsd=0.5 * (kl(o) + kl(q)),
        bray_curtis=float(0.5 * np.abs(o - q).sum()),
        spearman_log=float(rho_s),
        median_abs_log10_ratio=float(np.median(np.abs(lo - lq))),
        retained=int(kept.sum()),
        retained_fraction=float(kept.mean()),
        rare_genera=int(rare_mask.sum()),
        rare_retained=int((kept & rare_mask).sum()),
        rare_retained_fraction=float((kept & rare_mask).sum() / rare_mask.sum()) if rare_mask.any() else float("nan"),
    )
