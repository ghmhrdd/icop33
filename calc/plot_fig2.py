import json, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "DejaVu Serif", "mathtext.fontset": "dejavuserif", "font.size": 8,
                     "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": False, "ytick.right": True, "pdf.fonttype": 42})

def best(p, nv0, i):
    c = [json.load(open(f)) for f in [f"out/{t}_{p}_{nv0}_{i:02d}.json" for t in ("f2", "s2", "s3", "s4", "s5")] if os.path.exists(f)]
    return min(c, key=lambda r: r["eta"])

N = 12
D = {(p, nv0): [best(p, nv0, i) for i in range(N)] for p in ("cw", "pulsed_odmr", "ramsey") for nv0 in (0, 1)}
T = np.array([r["T2s"] for r in D[("cw", 0)]]) * 1e6
E = {k: np.array([r["eta"] for r in v]) for k, v in D.items()}
json.dump({"T2s_us": T.tolist(), **{f"{p}_{nv0}": E[(p, nv0)].tolist() for (p, nv0) in E}}, open("out/fig2_summary.json", "w"), indent=1)

COL = {"cw": "#2a78d6", "pulsed_odmr": "#eb6834", "ramsey": "#1baf7a"}
MK = {"cw": "o", "pulsed_odmr": "s", "ramsey": "^"}
LAB = {"cw": "CW ODMR", "pulsed_odmr": "pulsed ODMR", "ramsey": "Ramsey"}
fig = plt.figure(figsize=(3.17, 4.33))
ax = fig.add_axes([0.20, 0.42, 0.76, 0.47]); bx = fig.add_axes([0.20, 0.10, 0.76, 0.28], sharex=ax)
for p in ("cw", "pulsed_odmr", "ramsey"):
    ax.plot(T, E[(p, 1)], "-" + MK[p], color=COL[p], ms=3.5, lw=1.4, label=LAB[p])
    ax.plot(T, E[(p, 0)], "--" + MK[p], color=COL[p], ms=3.5, lw=1.2, mfc="w")
ax.plot([], [], "-o", color="0.6", ms=3.5, label=r"solid: with NV$^0$")
ax.plot([], [], "--o", color="0.6", ms=3.5, mfc="w", label=r"dashed: NV$^-$ only")
ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(True, lw=0.3, alpha=0.5)
ax.set_ylabel(r"$\eta_{\rm opt}$ (nT/$\sqrt{\rm Hz}$)", labelpad=1)
ax.legend(fontsize=6.3, frameon=False, loc="lower left"); plt.setp(ax.get_xticklabels(), visible=False)
tP = 0.53 / 0.3
ax.axvline(tP, color="0.4", lw=0.7, ls=":"); bx.axvline(tP, color="0.4", lw=0.7, ls=":")
ax.set_ylim(50, 6500)
ax.text(tP * 1.07, 0.975, "Poulsen et al.\n0.3 MHz", transform=ax.get_xaxis_transform(), fontsize=5.5, color="0.35", va="top")
top = ax.secondary_xaxis("top", functions=(lambda t: 0.53 / t, lambda L: 0.53 / L))
top.set_xlabel(r"inhomogeneous FWHM $L_{\rm IHB}$ (MHz)", labelpad=3)
ax.text(0.03, 0.97, "(a)", transform=ax.transAxes, va="top")
for nv0, ls, mfc in ((1, "-", None), (0, "--", "w")):
    bx.plot(T, E[("cw", nv0)] / E[("pulsed_odmr", nv0)], ls + "s", color=COL["pulsed_odmr"], ms=3.5, lw=1.3, mfc=mfc)
    bx.plot(T, E[("pulsed_odmr", nv0)] / E[("ramsey", nv0)], ls + "^", color=COL["ramsey"], ms=3.5, lw=1.3, mfc=mfc)
bx.axhline(1, color="k", lw=0.5); bx.grid(True, lw=0.3, alpha=0.5)
bx.set_xlabel(r"$T_2^*$ ($\mu$s)", labelpad=1); bx.set_ylabel("sensitivity gain", labelpad=1)
bx.text(0.03, 0.93, "(b)", transform=bx.transAxes, va="top")
ym = max((E[("cw", 1)] / E[("pulsed_odmr", 1)]).max(), (E[("pulsed_odmr", 0)] / E[("ramsey", 0)]).max())
bx.set_ylim(0, ym * 1.08)
bx.text(0.36, 0.86, r"$\eta_{\rm CW}/\eta_{\rm pulsed}$", transform=bx.transAxes, color=COL["pulsed_odmr"], fontsize=6.5)
bx.text(0.62, 0.30, r"$\eta_{\rm pulsed}/\eta_{\rm Ramsey}$", transform=bx.transAxes, color=COL["ramsey"], fontsize=6.5)
OUT_FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(OUT_FIG, exist_ok=True)
for ext, kw in (("pdf", {}), ("jpg", {"dpi": 600, "pil_kwargs": {"quality": 95}}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(OUT_FIG, "fig2." + ext), **kw)
np.set_printoptions(linewidth=200, precision=3, suppress=True)
print("T", T)
for k in E: print(k, np.round(E[k]))
for nv0 in (0, 1):
    print("CW/pulsed", nv0, E[("cw", nv0)] / E[("pulsed_odmr", nv0)]); print("pulsed/Ramsey", nv0, E[("pulsed_odmr", nv0)] / E[("ramsey", nv0)])
for p in E: print("NV0 penalty", p[0], E[(p[0], 1)] / E[(p[0], 0)]) if p[1] == 0 else None
