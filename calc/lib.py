"""Shared helpers for the CW / pulsed ODMR / Ramsey comparison, built on nvkit's protocols module.
Global search (differential evolution) + Powell polish on nvkit's own knobs,
then every optimum is re-evaluated with finer T2* quadratures (41 and 61 nodes)."""
import json, os, sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import nvkit
from nvkit import core, protocols as P

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
GAMMA = 2 * np.pi * core.gamma_e


def make(T2s, nv0=False, n_nodes=15, **kw):
    """Paper setup: B = 10 mT, T1 = 5 ms, T2 = max(5 us, 2 T2*), eta_det = 0.05, dead 1 us."""
    base = dict(B=10e-3, T1=5e-3, T2=max(5e-6, 2 * T2s), T2_star=T2s, eta_det=0.05, t_dead=1e-6,
                max_rabi=10e6, charge=nvkit.ChargeModel() if nv0 else None, n_nodes=n_nodes)
    base.update(kw)
    return P._make_setup(**base)


def with_nodes(st, n):
    st.nodes, st.wts = core.gauss_hermite_detunings(st.sigma, n)


def optimise(protocol, T2s, nv0=False, n_search=9, seed=1, popsize=6, maxiter=30, **kw):
    t0 = time.time()
    kw_pair = kw.pop('pair', None)
    x0_from = kw.pop('x0_from', None)
    opt = kw.pop('opt', {})
    seed = opt.get('seed', seed); popsize = opt.get('popsize', popsize); maxiter = opt.get('maxiter', maxiter)
    n_search = opt.get('n_search', n_search); free = opt.get('free_ramsey', False)
    st = make(T2s, nv0, max(kw_pair) if kw_pair else n_search, **kw)
    knobs = P._knobs(protocol, st)
    fixed = {"t_dead": st.t_dead} if protocol != "cw" else {}
    if protocol == "ramsey":
        fixed.update(rabi=st.max_rabi, detuning=0.0)
        if free:      # also search the Rabi frequency and the drive detuning
            knobs = knobs + [("rabi", 0.2e6, st.max_rabi, True), ("detuning", -3e6, 3e6, False)]
            fixed = {"t_dead": st.t_dead}
    nev = [0]

    pair = kw_pair
    if pair:
        quads = [core.gauss_hermite_detunings(st.sigma, n) for n in pair]

    def cost(x):
        nev[0] += 1
        try:
            s_ = P._from_x(protocol, st, x, knobs, fixed)
            if pair:      # worse of two quadratures: an artefact of one node set cannot win
                e = 0.0
                for q in quads:
                    st.nodes, st.wts = q
                    e = max(e, P._evaluate(protocol, st, s_).eta)
            else:
                e = P._evaluate(protocol, st, s_).eta
        except Exception:
            return 1e3
        return float(np.log(e)) if np.isfinite(e) and e > 0 else 1e3

    start = P._default_settings(protocol, st)
    start["t_laser"] = max(start.get("t_laser", 0.0), start.get("t_read", 0.0) + 50e-9)
    x0 = P._to_x(knobs, start)
    nk = len(knobs)
    if opt.get('no_x0'):
        x0 = None
    if protocol == "cw":
        popsize, maxiter = 20, 80
    if x0_from:
        ws = json.load(open(os.path.join(OUT, x0_from + ".json")))["settings"]
        x0 = P._to_x(knobs, ws)
        maxiter, popsize = 1, 2
    r = differential_evolution(cost, [(0, 1)] * nk, seed=seed, popsize=popsize, maxiter=maxiter, tol=1e-4,
                               x0=x0, polish=False, init="sobol", mutation=(0.4, 1.0), recombination=0.8)
    xb, cb = r.x, r.fun
    # insurance against a missed basin: a local search from nvkit's default start
    rp = minimize(cost, x0 if x0 is not None else xb, method="Powell", bounds=[(0, 1)] * nk, options={"xtol": 1e-4, "ftol": 1e-6, "maxiter": 400} if x0_from else {"xtol": 1e-3, "ftol": 1e-4, "maxiter": 200})
    if rp.fun < cb:
        xb, cb = np.clip(rp.x, 0, 1), rp.fun
    for _ in range(3):
        rp = minimize(cost, xb, method="Powell", bounds=[(0, 1)] * nk, options={"xtol": 1e-4, "ftol": 1e-6, "maxiter": 400})
        if rp.fun < cb - 1e-6:
            xb, cb = np.clip(rp.x, 0, 1), rp.fun
        else:
            break
    s = P._from_x(protocol, st, xb, knobs, fixed)
    res = {"protocol": protocol, "T2s": T2s, "nv0": nv0, "kw": {k: v for k, v in kw.items()},
           "settings": s, "eta_search": float(np.exp(cb)) * 1e9, "nev": nev[0]}
    if protocol == "cw":
        ev = P._evaluate(protocol, st, s)
        res["eta"] = ev.eta * 1e9
        res["t_shot"] = 1.0
    else:
        etas = {}
        for n in (41, 61):
            with_nodes(st, n)
            ev = P._evaluate(protocol, st, s)
            etas[n] = ev.eta * 1e9
        res["eta_41"], res["eta_61"] = etas[41], etas[61]
        res["eta"] = etas[61]
        res["t_shot"] = ev.t_shot
        res["rabi"] = 0.5 / s["t_p"] if protocol == "pulsed_odmr" else s["rabi"]
    res["counts"] = {k: float(v) for k, v in ev.counts.items()}
    res["time_s"] = time.time() - t0
    return res


def run_task(task):
    name, protocol, T2s, nv0, kw = task
    path = os.path.join(OUT, name + ".json")
    if os.path.exists(path):
        return name
    best = None
    for seed in ((11, 12, 13, 14, 15) if protocol == "cw" else (1,)):
        res = optimise(protocol, T2s, nv0, seed=seed, **{k: (dict(v) if isinstance(v, dict) else v) for k, v in kw.items()})
        if best is None or res["eta"] < best["eta"]:
            best = res
    best["name"] = name
    with open(path, "w") as f:
        json.dump(best, f, indent=1, default=float)
    print(f"{name}: eta={best['eta']:.1f} search={best['eta_search']:.1f} nev={best['nev']} t={best['time_s']:.0f}s", flush=True)
    return name


# ---------------- CW helpers (exact T2* average, vectorised over detuning) ----------------
def cw_rate(st, laser, rabi, det_hz, h=0.0):
    """Emitted photon rate (s^-1, before eta_det) at detunings det_hz (array, Hz from the line)."""
    eng = st.engine
    d = 2 * np.pi * (st.line + np.asarray(det_hz, float))
    R = 0.0
    for nl in eng.branches:
        R = R + eng.weights[nl] * np.asarray(eng.cw_average(nl, laser, 2 * np.pi * rabi, d, st.sigma, h), float)
    return R


def cw_eta(st, laser, rabi, det_hz):
    """CRB sensitivity (T/sqrt(Hz)) at an array of detunings."""
    h = P._DB_STEP
    R0 = cw_rate(st, laser, rabi, det_hz)
    dR = (cw_rate(st, laser, rabi, det_hz, h) - cw_rate(st, laser, rabi, det_hz, -h)) / (2 * h)
    a = st.n_nv * st.eta_det
    mu = a * R0 + st.background
    F = (a * dR) ** 2 / mu
    return 1.0 / np.sqrt(np.maximum(F, 1e-300))


def cw_best_detuning(st, laser, rabi):
    det = np.concatenate([np.geomspace(2e3, 30e6, 500)])
    e = cw_eta(st, laser, rabi, det)
    i = int(np.argmin(e))
    lo, hi = det[max(i - 1, 0)], det[min(i + 1, len(det) - 1)]
    fine = np.linspace(lo, hi, 200)
    e2 = cw_eta(st, laser, rabi, fine)
    j = int(np.argmin(e2))
    return float(e2[j]), float(fine[j])


def slope_formula(st, laser, rabi):
    """Textbook eta = 0.77 dnu / (gamma C sqrt(R0)) with width, contrast and rate read off
    the computed spectrum: R0 = detected rate without MW, C = depth of the deepest dip,
    dnu = full width of that dip at half depth."""
    det = np.linspace(-15e6, 15e6, 12001)
    R = cw_rate(st, laser, rabi, det)
    R0 = float(cw_rate(st, laser, 0.0, np.array([0.0]))[0])
    i = int(np.argmin(R))
    C = (R0 - R[i]) / R0
    half = R0 - 0.5 * (R0 - R[i])
    l = i
    while l > 0 and R[l] < half:
        l -= 1
    r = i
    while r < len(R) - 1 and R[r] < half:
        r += 1
    dnu = det[r] - det[l]
    a = st.n_nv * st.eta_det
    eta = 0.77 * dnu / (core.gamma_e * C * np.sqrt(a * R0))
    return float(eta), float(dnu), float(C), float(a * R0)
