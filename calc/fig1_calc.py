import os, json
os.environ["OMP_NUM_THREADS"] = "1"
import multiprocessing as mp
import numpy as np
import lib

T2S = 2e-6
RABI = np.geomspace(0.05e6, 10e6, 16)
LASER = np.geomspace(1e-3, 1.0, 12)
CUT = np.geomspace(0.3e6, 10e6, 9)

def row(args):
    kind, nv0, laser, rabis = args
    st = lib.make(T2S, nv0)
    out = []
    for r in rabis:
        e, d = lib.cw_best_detuning(st, laser, r)
        rec = {"laser": float(laser), "rabi": float(r), "eta": e, "det": d}
        if kind != "map0":
            s, dnu, C, R0 = lib.slope_formula(st, laser, r)
            rec.update(slope=s, dnu=dnu, C=C, R0=R0)
        out.append(rec)
    return kind, nv0, out

if __name__ == "__main__":
    jobs = [("map", False, l, RABI) for l in LASER] + [("map0", True, l, RABI) for l in LASER]
    with mp.Pool(4) as pool:
        res = pool.map(row, jobs)
    data = {"T2s": T2S, "rabi": RABI.tolist(), "laser": LASER.tolist(),
            "map": [r[2] for r in res if r[0] == "map"], "map_nv0": [r[2] for r in res if r[0] == "map0"]}
    # optima (global search) and 1D cuts
    for nv0 in (False, True):
        r = lib.optimise("cw", T2S, nv0)
        data["opt_nv0" if nv0 else "opt"] = {"eta": r["eta"], **r["settings"]}
    lc = data["opt"]["laser"]
    data["cut"] = row(("cut", False, lc, CUT))[2]
    data["cut_nv0"] = row(("cut", True, 1.0, CUT))[2]
    data["cut_fine"] = row(("cut", False, lc, np.geomspace(0.3e6, 10e6, 41)))[2]
    json.dump(data, open(os.path.join(lib.OUT, "fig1.json"), "w"), indent=1)
    print("fig1 done", data["opt"], data["opt_nv0"])
