import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import multiprocessing as mp
import numpy as np
import lib

T2S = [0.25e-6, 0.30e-6, 0.35e-6, 0.40e-6] + list(np.geomspace(0.465e-6, 10e-6, 8))
TAB = 1.874e-6
tasks = []
# table + robustness first (most quoted numbers), then the sweep (long T2* first: slowest)
for p in ("ramsey", "pulsed_odmr", "cw"):
    for nv0 in (False, True):
        tasks.append((f"tab_{p}_{int(nv0)}", p, TAB, nv0, {}))
for p in ("ramsey", "pulsed_odmr", "cw"):
    for nv0 in (False, True):
        for i, T in sorted(enumerate(T2S), key=lambda t: -t[1]):
            tasks.append((f"f2_{p}_{int(nv0)}_{i:02d}", p, float(T), nv0, {}))
ROB = {"dead03": dict(t_dead=0.3e-6), "dead3": dict(t_dead=3e-6), "eta001": dict(eta_det=0.01),
       "eta02": dict(eta_det=0.2), "bg5k": dict(background=5000.0), "bins20": dict(bins=20),
       "rabi2": dict(max_rabi=2e6)}
for key, kw in ROB.items():
    for p in ("ramsey", "pulsed_odmr", "cw"):
        tasks.append((f"rob_{key}_{p}", p, TAB, False, kw))

for p in ("ramsey", "pulsed_odmr"):
    for nv0 in (False, True):
        for i in (4, 3, 2, 1, 0):
            tasks.append((f"s2_{p}_{int(nv0)}_{i:02d}", p, float(T2S[i]), nv0, {"pair": (21, 31)}))

for i in (4, 3, 2, 1, 0):
    tasks.append((f"s3_pulsed_odmr_1_{i:02d}", "pulsed_odmr", float(T2S[i]), True, {"pair": (21, 31), "x0_from": f"s2_pulsed_odmr_0_{i:02d}"}))

for i in range(12):
    for nv0 in (False, True):
        tasks.append((f"s4_cw_{int(nv0)}_{i:02d}", "cw", float(T2S[i]), nv0, {}))
    for p in ("ramsey", "pulsed_odmr"):
        if i >= 5:
            tasks.append((f"s4_{p}_1_{i:02d}", p, float(T2S[i]), True, {"x0_from": f"f2_{p}_0_{i:02d}"}))
        elif p == "ramsey":
            tasks.append((f"s4_{p}_1_{i:02d}", p, float(T2S[i]), True, {"pair": (21, 31), "x0_from": f"s2_{p}_0_{i:02d}"}))

for i in range(5, 12):
    for p in ("ramsey", "pulsed_odmr"):
        src = min((f"{t}_{p}_1_{i:02d}" for t in ("f2", "s4") if os.path.exists(os.path.join(lib.OUT, f"{t}_{p}_1_{i:02d}.json"))),
                  key=lambda n: __import__("json").load(open(os.path.join(lib.OUT, n + ".json")))["eta"])
        tasks.append((f"s5_{p}_0_{i:02d}", p, float(T2S[i]), False, {"x0_from": src}))

for p in ("ramsey", "pulsed_odmr"):
    for nv0 in (False, True):
        for sd in (7, 23):
            tasks.append((f"v1_{p}_{int(nv0)}_s{sd}", p, TAB, nv0, {"opt": dict(seed=sd, popsize=10, maxiter=40, n_search=15, no_x0=True)}))
for nv0 in (False, True):
    tasks.append((f"v2_ramseyfree_{int(nv0)}", "ramsey", TAB, nv0, {"opt": dict(seed=5, popsize=8, maxiter=40, n_search=15, free_ramsey=True)}))
for i in (5, 9, 11):
    tasks.append((f"v2_ramseyfree_0_T{i:02d}", "ramsey", float(T2S[i]), False, {"opt": dict(seed=5, popsize=8, maxiter=40, n_search=15, free_ramsey=True)}))

if __name__ == "__main__":
    sel = sys.argv[1] if len(sys.argv) > 1 else ""
    todo = sorted([t for t in tasks if t[0].startswith(sel)], key=lambda t: -int(t[3]))
    print(len(todo), "tasks", flush=True)
    with mp.Pool(4, maxtasksperchild=1) as pool:
        for n in pool.imap_unordered(lib.run_task, todo, chunksize=1):
            pass
    print("ALL DONE", flush=True)
