"""Integrate the existing carbon MiCRM with BDF and an absolute-rate event."""

import json
import sys
import numpy as np
from scipy.integrate import solve_ivp


def integrate(p):
    n, r = int(p["n"]), int(p["r"])
    u = np.asarray(p["U"]).reshape((n, r), order="F")
    leakage = np.asarray(p["L"]).reshape((n, r, r), order="F")
    loss = np.asarray(p["loss"]).reshape((n, r), order="F")
    maintenance = np.asarray(p["m"])
    feed = np.asarray(p["feed"])
    initial = np.atleast_1d(p["initial"])
    k = None if p["K"] is None else np.asarray(p["K"]).reshape((n, r), order="F")
    retained = 1 - leakage.sum(axis=2) - loss
    leakage_flat = leakage.reshape((n * r, r))
    tolerance = float(p["tolerance"])
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Equilibrium tolerance must be finite and positive")

    def rhs(t, state):
        state = np.maximum(state, 0)
        biomass, resources = state[:n], state[n:]
        uptake = u * resources * biomass[:, None]
        if k is not None:
            uptake = uptake / (1 + resources / k)
        dc = (uptake * retained).sum(axis=1) - maintenance * biomass
        dr = (
            p["turnover"] * (feed - resources)
            - uptake.sum(axis=0)
            + uptake.ravel() @ leakage_flat
        )
        return np.concatenate((dc, dr))

    def equilibrium(t, state):
        return np.max(np.abs(rhs(t, state))) - tolerance

    equilibrium.terminal = True
    equilibrium.direction = -1
    y0 = np.concatenate((np.full(n, p["biomass"]), initial))
    if equilibrium(0, y0) <= 0:
        times, states, early = np.array([0.0]), y0[:, None], True
    else:
        solution = solve_ivp(
            rhs,
            (0, p["t_end"]),
            y0,
            method="BDF",
            atol=p["atol"],
            rtol=p["rtol"],
            events=equilibrium,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        times, states, early = solution.t, solution.y, solution.status == 1
    if not np.all(np.isfinite(states)):
        raise ValueError("Nonfinite integration output")
    states = np.maximum(states, 0)
    return dict(
        time=float(times[-1]),
        state=states[:, -1].tolist(),
        max_resource=states[n:].max(axis=1).tolist(),
        early_stopped=early,
        steps=len(times) - 1,
        trajectory=np.column_stack((times, states.T)).tolist(),
    )


if __name__ == "__main__":
    with open(sys.argv[1]) as stream:
        parameters = json.load(stream)
    result = integrate(parameters)
    with open(sys.argv[2], "w") as stream:
        json.dump(result, stream, allow_nan=False)
