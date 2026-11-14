# Supplementary material: model, parameters and reproducibility

Supplement to the paper *A photophysical model for numerical simulation of magnetometry with nitrogen-vacancy centers in diamond* (33rd Iranian Conference on Optics and Photonics, ICOP 2026).
Every number below is the value used in the calculations (nvkit `0.2.0`, commit `23f6d18`; Python 3.13, NumPy 2.5, SciPy 1.18, QuTiP 5.3).

## S1. Hamiltonian

Basis per ¹⁴N branch m_I ∈ {+1, 0, −1}: NV⁻ ground triplet |g, m_s⟩, excited triplet |e, m_s⟩, one effective intermediate (singlet) level |s⟩, plus, with charge dynamics, NV⁰ ground |n0g⟩ and excited |n0e⟩ (9 levels per branch).

Ground-state Hamiltonian (Hz), field B along the NV axis:

$$
H_g/h = D\,S_z^2 + \gamma_e B\,S_z + A_\parallel S_z I_z + Q\,I_z^2 - \gamma_n B\,I_z
$$

| Symbol | Value | Note |
|---|---|---|
| D | 2.87 GHz | ground-state zero-field splitting |
| γ_e | 28.03 GHz/T | g = 2.0028 |
| A∥ | −2.14 MHz | secular (longitudinal) ¹⁴N hyperfine |
| A⊥ | −2.70 MHz | **dropped** (secular approximation): couples levels ≈ D apart, shift ≈ A⊥²/D ≈ 2.5 kHz |
| Q | −5.01 MHz | ¹⁴N quadrupole |
| γ_n | 3.077 MHz/T | ¹⁴N gyromagnetic ratio |
| E (strain) | 0 | |
| B | 10 mT | bias; ms = 0 ↔ −1 line at D − γ_e B ≈ 2589.7 MHz, hyperfine lines at −A∥ m_I (2.14 MHz apart) |

Microwave drive, rotating-wave approximation, frame rotating at ν_mw = D + δ/2π:

$$
H_{\rm rot}/\hbar = 2\pi\left(H_g/h - \nu_{\rm mw} S_z^2\right) + \frac{\Omega_R}{\sqrt2}\,S_x
$$

Ω_R is the on-resonance Rabi (angular) frequency of each ms = 0 ↔ ±1 transition; the "Rabi frequency" quoted in the paper is Ω_R/2π. The ms = +1 transition is in H but lies ≈ 560 MHz away at 10 mT, so it is not excited. Each m_I branch is solved separately and carries its own detuning δ + 2πA∥m_I (this includes off-resonant driving of the outer hyperfine lines *during* every pulse; nothing is idealised as an instantaneous rotation). The excited state carries no coherent dynamics (only rates, below); ¹⁴N hyperfine is neglected there (far from the 51 mT excited-state anticrossing).

## S2. Incoherent processes (Lindblad jump operators)

Optical / ISC model of Robledo et al., NJP 13, 025013 (2011). Laser enters only through the pump rate k_p (spin-conserving).

| Process | Jump | Rate |
|---|---|---|
| optical excitation | g, m_s → e, m_s | k_p = laser × 16.9 µs⁻¹ (cap ≈ ¼ saturation; peak value during the laser pulse) |
| radiative decay | e, m_s → g, m_s | k_r = 65 µs⁻¹ |
| ISC from e, ±1 | → s | 79 µs⁻¹ (each) |
| ISC from e, 0 | → s | 11 µs⁻¹ |
| singlet decay | s → g, 0 | 3 µs⁻¹ |
| singlet decay | s → g, ±1 | 1 µs⁻¹ (each) |
| spin relaxation | g,0 ↔ g,±1 | 1/(3T₁) each direction (T₁ = relaxometry time, 5 ms) |
| pure dephasing | two operators √γ_φ (\|0⟩⟨0\| − \|∓1⟩⟨∓1\|) | γ_φ = 1/(2.5 T₂), so the coherence decays at 1/T₂ (T₂ = Hahn-echo time = max(5 µs, 2T₂*)) |

Laser-off (dark) intervals keep relaxation, dephasing and the decay of e and s.

## S3. NV⁰ charge model (`nvkit.ChargeModel()` defaults)

| Process | Jump | Rate at the laser cap (k_p = 16.9 µs⁻¹) |
|---|---|---|
| ionisation | e, m_s (NV⁻) → n0g | β_ion k_p, β_ion = 0.05 → 0.845 µs⁻¹ |
| NV⁰ excitation | n0g → n0e | k_p (abs0 = 1) |
| NV⁰ radiative decay | n0e → n0g | k_r0 = 1/19 ns = 52.6 µs⁻¹ |
| recombination | n0e → g, m_s (NV⁻), ⅓ each | β_rec k_p, β_rec = 0.0667 → 1.13 µs⁻¹ |

β_rec is calibrated so that the laser-on steady NV⁻ fraction at the cap is f₋ = 0.75 (Aslam et al., NJP 15, 013064 (2013), order of magnitude). NV⁰ photons are detected with relative efficiency η₀ = 0.4 (red tail through a ≈650 nm long-pass). The ¹⁴N spin is a spectator of charge transitions; NV⁰ is inert in the dark.

**Sensitivity to these assumptions** (T₂* = 1.874 µs, re-optimised; `calc/nv0_sensitivity.py`, `calc/out/nv0b_*`, `calc/out/nv0f_*`). NV⁰ penalty = η(with NV⁰)/η(NV⁻ only) − 1:

| Assumption | Ramsey | pulsed | CW | Ramsey/pulsed |
|---|---|---|---|---|
| β_ion = 0.025 | 39 % | 36 % | 33 % | 3.33 |
| β_ion = 0.05 (paper) | 42 % | 39 % | 35 % | 3.33 |
| β_ion = 0.10 | 48 % | 45 % | 40 % | 3.32 |
| f₋ = 0.85 | 24 % | 23 % | 20 % | 3.35 |
| f₋ = 0.65 | 64 % | 60 % | 54 % | 3.31 |

The penalty is set mainly by the steady NV⁻ fraction, not by the ionisation rate; the ranking and ratios do not change.

## S4. Inhomogeneous dephasing (T₂*)

Quasi-static detuning δ_inh ~ Gaussian with standard deviation σ = √2/T₂* (rad/s), i.e. a free-induction envelope exp[−(t/T₂*)²] and FWHM linewidth 0.53/T₂* (Hz × 2π convention: FWHM = 2√(2 ln2)σ/2π). T₂* is therefore the Gaussian (inhomogeneous) part only; with the Markovian T₂ above the total envelope reaches 1/e at 0.93 T₂* for T₂* = 1.9 µs. The average uses Gauss–Hermite quadrature (search with 9 nodes, or the worse of 21 and 31 nodes for T₂* < 0.5 µs; every optimum re-evaluated with 41 and 61 nodes, change < 1 %). CW uses the exact (Faddeeva) average.

## S5. Measurement model and figure of merit

* Sequences (periodic steady state, the cycle repeated indefinitely):
  * pulsed ODMR: MW π pulse (length t_p, detuning δ) → laser t_laser (read the first t_read) → dark t_dead;
  * Ramsey: π/2 → free evolution τ → π/2 with phase φ → laser t_laser (read the first t_read) → dark t_dead; Rabi at the cap (π/2 = 25 ns), drive on the central hyperfine line (δ = 0);
  * CW: laser and MW on, steady state.
* Detected photon rate = η_det [k_r P(e) + η₀ k_r0 P(n0e)] + background, η_det = 0.05, background 0 (5000 s⁻¹ in the robustness test), dead time 1 µs.
* **Nuclear mixture:** expected counts N_k and their field derivatives ∂_B N_k are averaged over the three m_I branches (weight ⅓ each, unpolarised ¹⁴N) and over the T₂* nodes **before** the Fisher information is formed, so F is that of the photon distribution actually detected (no overestimate from summing per-branch Fisher information).
* Poisson statistics (< 1 detected photon per cycle): F = Σ_k (∂_B N_k)² / (N_k + N_bg,k) per cycle, η = √(T_cyc / F). ∂_B is analytic (finite-difference check: relative difference ≈ 2×10⁻⁶).
* η is a local Cramér–Rao bound: a small δB around the 10 mT bias with all other parameters known; the working point (detuning, phase, τ) is optimised, which places it at maximum slope.

## S6. Optimisation

Knobs and bounds (nvkit `protocols._knobs`):

| Protocol | Knobs (bounds; log-scaled unless noted) |
|---|---|
| CW | laser 10⁻³–1 (fraction of cap), Rabi 10 kHz–cap, detuning 1 kHz–30 MHz |
| pulsed | t_p from 0.5/Rabi_cap to min(4T₂*, T₂/2, 50 µs, quadrature limit); detuning 1 kHz–30 MHz; laser 0.01–1; t_read 20 ns–3 µs; t_laser − t_read 5 ns–20 µs |
| Ramsey | τ 10 ns–min(3T₂*, T₂, quadrature limit); phase 0–360° (linear); laser, t_read, t_laser as pulsed |

Global search: SciPy `differential_evolution` (Sobol init, popsize 6, maxiter 30, mutation 0.4–1.0, recombination 0.8, tol 10⁻⁴, seed 1) on log η, then Powell from nvkit's default settings, then Powell polishing (xtol 10⁻⁴, ftol 10⁻⁶, up to 3 rounds). CW: 5 seeds, popsize 20, maxiter 80. Fig. 2 points were refined by warm starts from neighbouring optima and the best run per point is kept (`plot_fig2.py`).

Convergence evidence: the Table-1 pulsed and Ramsey optima (with and without NV⁰) were reproduced to 0.1 nT/√Hz by two independent searches (seeds 7 and 23, popsize 10, maxiter 40, 15 nodes, no start point; `calc/out/v1_*`); freeing the Ramsey Rabi frequency and drive detuning never beat Rabi at the cap and zero detuning (`calc/out/v2_*`).

## S7. Optimal settings at T₂* = 1.874 µs (no NV⁰)

| | CW | pulsed | Ramsey |
|---|---|---|---|
| η (nT/√Hz) | 1845 | 639 | 189 |
| laser (fraction of cap) | 0.70 | 1.00 | 1.00 |
| Rabi Ω_R/2π (MHz) | 2.26 | 0.30 | 10 (cap) |
| detuning from central line (MHz) | +3.59 | +0.25 | 0 |
| t_p or τ (µs) | — | 1.65 | 0.904 |
| Ramsey 2nd-pulse phase | — | — | 268.9° (quadrature) |
| t_laser / t_read (µs) | — | 0.49 / 0.46 | 0.46 / 0.43 |
| cycle (µs) | — | 3.14 | 2.41 |

**Ramsey working point stability:** η rises by 1.5 % for φ ± 10°, 6 % for ± 20°; by 14 % for τ ± 5 % and 78 % for τ ± 10 % (τ sits on the second ¹⁴N hyperfine revival, 2 × 463 ns).

## S8. Reproducing

```bash
pip install -e <nvkit checkout at commit 23f6d18>
cd calc
python run_all.py            # all optimisations → out/*.json (≈ 2.5 h on 4 cores)
python fig1_calc.py && python plot_fig1.py
python plot_fig2.py          # also writes out/fig2_summary.json
python nv0_sensitivity.py    # NV0 sensitivity table (Sec. S3)
```

See `README.md` for the file map.
