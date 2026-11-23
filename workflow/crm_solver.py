"""Carbon MiCRM integrated in log-biomass coordinates with an analytic Jacobian.

    uptake_ia = U_ia * C_i * f_a(R_a),  f_a(R) = R / (1 + R/K_ia)
    dC_i/dt   = sum_a uptake_ia * eps_ia - m_i * C_i
    dR_a/dt   = tau (rho_a - R_a) - sum_i uptake_ia + sum_ib uptake_ib L_iba

Biomass is integrated as ln C, so every genus keeps relative precision (on a
linear scale with an absolute tolerance, a genus at 1e-6 relative abundance
sits in the solver's noise floor). The Jacobian is exact and convergence is
judged per capita, so a genus is only called extinct when its own growth rate
is negative at the final resource levels.
"""

import numpy as np
from scipy.integrate import solve_ivp


class CRM:
    def __init__(self, U, L, loss, m, K, rho, tau):
        self.U = np.asarray(U, float)
        self.L = np.asarray(L, float)
        self.loss = np.asarray(loss, float)
        self.m = np.asarray(m, float)
        K = np.asarray(K, float)
        self.invK = np.where(np.isfinite(K) & (K > 0), 1.0 / np.where(K > 0, K, 1), 0.0)
        self.rho = np.asarray(rho, float)
        self.tau = float(tau)
        self.n, self.r = self.U.shape
        self.eps = 1.0 - self.L.sum(axis=2) - self.loss
        self.Lflat = self.L.reshape(self.n * self.r, self.r)

    # per-capita uptake G_ia and its derivative w.r.t. R_a
    def uptake(self, R):
        R = np.maximum(R, 0.0)
        d = 1.0 + R[None, :] * self.invK
        return self.U * R[None, :] / d, self.U / d**2

    def growth(self, R):
        """Per-capita net growth rate of every genus at resource levels R."""
        G, _ = self.uptake(R)
        return (G * self.eps).sum(axis=1) - self.m

    def rhs(self, t, y):
        n = self.n
        C = np.exp(np.minimum(y[:n], 700.0))
        R = y[n:]
        G, _ = self.uptake(R)
        dlnC = (G * self.eps).sum(axis=1) - self.m
        flux = G * C[:, None]
        dR = self.tau * (self.rho - np.maximum(R, 0.0)) - flux.sum(axis=0) + flux.ravel() @ self.Lflat
        return np.concatenate((dlnC, dR))

    def jac(self, t, y):
        n, r = self.n, self.r
        C = np.exp(np.minimum(y[:n], 700.0))
        R = y[n:]
        G, Gp = self.uptake(R)
        J = np.zeros((n + r, n + r))
        J[:n, n:] = Gp * self.eps
        # dR_a / d lnC_i = C_i * (-G_ia + sum_b G_ib L_iba)
        prod = np.einsum("ib,iba->ia", G, self.L)
        J[n:, :n] = ((-G + prod) * C[:, None]).T
        CGp = Gp * C[:, None]
        JRR = np.einsum("ib,iba->ab", CGp, self.L)
        JRR[np.diag_indices(r)] -= self.tau * (R >= 0) + CGp.sum(axis=0)
        J[n:, n:] = JRR
        return J


def simulate(model, C0, R0, t_end=3e5, chunk=2000.0, growth_tol=1e-7,
             resource_tol=1e-9, extinct=1e-12, rtol=1e-7, atol=1e-10):
    """Integrate until every surviving genus has |growth| < growth_tol and the
    resources have stopped moving, or until t_end."""
    n = model.n
    y = np.concatenate((np.log(np.maximum(np.asarray(C0, float), 1e-300)), np.asarray(R0, float)))
    t, steps, converged = 0.0, 0, False
    while t < t_end:
        t1 = min(t_end, t + chunk)
        sol = solve_ivp(model.rhs, (t, t1), y, method="BDF", jac=model.jac,
                        rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        y, t = sol.y[:, -1], t1
        steps += sol.t.size - 1
        C = np.exp(np.minimum(y[:n], 700.0))
        rel = C / C.sum()
        d = model.rhs(t, y)
        alive = rel > extinct
        ok_c = np.all(np.abs(d[:n][alive]) < growth_tol)
        ok_r = np.all(np.abs(d[n:]) < resource_tol + 1e-6 * np.abs(y[n:]))
        if ok_c and ok_r:
            converged = True
            break
        chunk = min(chunk * 2, 1e5)
    C = np.exp(np.minimum(y[:n], 700.0))
    R = np.maximum(y[n:], 0.0)
    if not np.all(np.isfinite(C)) or C.sum() <= 0:
        raise ValueError("Nonfinite or extinct community")
    return dict(C=C, R=R, time=t, steps=steps, converged=converged,
                growth=model.growth(R))


def leading_eigenvalue(model, C, R):
    """Largest real part of the Jacobian at a state, in (C, R) coordinates
    (absent genera contribute their own per-capita growth rate)."""
    n = model.n
    J = model.jac(0.0, np.concatenate((np.log(np.maximum(C, 1e-300)), R)))
    G, Gp = model.uptake(R)
    J[:n, :n] = np.diag(model.growth(R))
    J[:n, n:] = Gp * model.eps * C[:, None]
    J[n:, :n] = (-G + np.einsum("ib,iba->ia", G, model.L)).T
    return float(np.max(np.linalg.eigvals(J).real))
