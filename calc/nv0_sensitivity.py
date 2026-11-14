"""NV0 sensitivity check at T2*=1.874 us: ionisation ratio beta_ion x0.5/x2 and steady NV- fraction 0.65/0.85,
warm-started from the table optima (tab_*_1)."""
import os, sys, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import multiprocessing as mp
import nvkit, lib
TAB = 1.874e-6
def job(a):
    p, key, tag, b = a
    out = os.path.join(lib.OUT, f'nv0{tag}_{p}_{b:g}.json')
    if os.path.exists(out): return out
    best = None
    for seed in ((11, 12) if p == 'cw' else (1,)):
        r = lib.optimise(p, TAB, True, seed=seed, charge=nvkit.ChargeModel(**{key: b}), x0_from=f'tab_{p}_1')
        if best is None or r['eta'] < best['eta']: best = r
    best['kw'] = {key: b}; best['name'] = os.path.basename(out)[:-5]
    json.dump(best, open(out, 'w'), indent=1, default=float)
    print(p, b, best['eta'], flush=True)
    return out
if __name__ == '__main__':
    with mp.Pool(4) as pool:
        list(pool.imap_unordered(job, [(p, k, t, b) for k, t, bs in (('beta_ion', 'b', (0.025, 0.1)), ('f_minus', 'f', (0.65, 0.85))) for b in bs for p in ('ramsey', 'pulsed_odmr', 'cw')]))
