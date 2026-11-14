# NV magnetometry: CW ODMR vs pulsed ODMR vs Ramsey

Code, raw results and supplementary material for the paper

> *A photophysical model for numerical simulation of magnetometry with nitrogen-vacancy centers in diamond*,
> 33rd Iranian Conference on Optics and Photonics (ICOP 2026).

The DC magnetic sensitivity of CW ODMR, single-tone pulsed ODMR and Ramsey interferometry is compared with one
Lindblad model of the NV⁻ centre (¹⁴N hyperfine structure, optional NV⁰ charge dynamics) and one figure of merit:
the Cramér–Rao bound of the detected photon counts per unit of total cycle time. Every setting of every protocol is
globally optimised under the same caps on the peak laser excitation rate and the Rabi frequency.

* **[SUPPLEMENT.md](SUPPLEMENT.md)**: Hamiltonian, all photophysical and charge-state rates, T₂* convention,
  Fisher-information treatment of the nuclear mixture, optimiser settings and bounds, optimal settings, convergence
  and sensitivity checks.
* **`calc/`**: scripts; **`calc/out/`**: every optimisation result (JSON); **`figures/`**: Fig. 1 and Fig. 2.

## Requirements

* [nvkit](https://github.com/ghmhrdd/nvkit) at commit `23f6d18` (version 0.2.0)
* Python ≥ 3.11, NumPy, SciPy, QuTiP 5, Matplotlib

```bash
pip install -e <nvkit checkout at commit 23f6d18>
cd calc
python run_all.py                  # all optimisations → out/*.json (≈ 2.5 h on 4 cores; existing files are skipped)
python fig1_calc.py && python plot_fig1.py
python plot_fig2.py                # best run per point → out/fig2_summary.json and figures/fig2.*
python nv0_sensitivity.py          # NV⁰ assumption checks
```

## Files

| Script | Purpose |
|---|---|
| `calc/lib.py` | setup (`make`), optimiser (`optimise`), CW helpers and the slope formula |
| `calc/run_all.py [prefix]` | all optimisation tasks (optionally only names starting with `prefix`) |
| `calc/fig1_calc.py` | CW sensitivity map and cuts → `out/fig1.json` |
| `calc/plot_fig1.py`, `calc/plot_fig2.py` | figures |
| `calc/nv0_sensitivity.py` | ionisation ratio ×0.5 / ×2 and steady NV⁻ fraction 0.65 / 0.85 |

Result files in `calc/out/` are named `<set>_<protocol>_<NV⁰: 0|1>_<T₂* index>`:

| Set | Content |
|---|---|
| `tab_*` | Table 1 (T₂* = 1.874 µs) |
| `f2_*`, `s2_*` … `s5_*` | Fig. 2 sweep and its refinements (finer quadrature for T₂* < 0.5 µs, warm starts); `plot_fig2.py` keeps the best per point |
| `rob_*` | robustness: dead time 0.3 / 3 µs, η_det 0.01 / 0.2, 5000 s⁻¹ background, 20 read bins, Rabi cap 2 MHz |
| `v1_*`, `v2_*` | independent-seed reproductions; Ramsey with free Rabi frequency and detuning |
| `nv0b_*`, `nv0f_*` | NV⁰ assumption checks |
| `fig1*.json`, `fig2_summary.json` | data plotted in Fig. 1 and Fig. 2 |

T₂* grid (µs): 0.25, 0.30, 0.35, 0.40, then 8 log-spaced points from 0.465 to 10. Each JSON holds the optimal
settings, the sensitivity in nT/√Hz (61-node quadrature; 41-node value alongside), cycle time and counts.
