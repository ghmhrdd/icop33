import json, os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import lib, fig1_calc

plt.rcParams.update({"font.family": "DejaVu Serif", "mathtext.fontset": "dejavuserif", "font.size": 8,
                     "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42})
d = json.load(open("out/fig1.json"))
p2 = "out/fig1_cut_nv0.json"
if not os.path.exists(p2):
    json.dump(fig1_calc.row(("cut", True, d["opt_nv0"]["laser"], fig1_calc.CUT))[2], open(p2, "w"))
cut_nv0 = json.load(open(p2))
rabi = np.array(d["rabi"]) / 1e6; laser = np.array(d["laser"])
m = np.array([[r["eta"] for r in row] for row in d["map"]]) * 1e6   # uT
o, o0 = d["opt"], d["opt_nv0"]
emin = o["eta"] / 1e3

BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#555555"
fig = plt.figure(figsize=(3.32, 4.1))
ax = fig.add_axes([0.17, 0.57, 0.58, 0.40]); cax = fig.add_axes([0.77, 0.57, 0.035, 0.40])
pm = ax.pcolormesh(rabi, laser, m, norm=LogNorm(vmin=emin * 0.98, vmax=60), cmap="Blues_r", shading="nearest", rasterized=True)
cs = ax.contour(rabi, laser, m / emin, levels=[1.1, 1.5, 2, 4], colors="w", linewidths=0.7)
ax.clabel(cs, fmt=lambda v: f"×{v:g}", fontsize=5.5, inline=True)
ax.axhline(o["laser"], color=ORANGE, lw=0.8, ls="--")
ax.axhline(o0["laser"], color=ORANGE, lw=0.8, ls="-")
ax.plot(o["rabi"] / 1e6, o["laser"], marker="*", ms=9, mfc=ORANGE, mec="w", mew=0.6, ls="none")
ax.plot(o0["rabi"] / 1e6, o0["laser"], marker="*", ms=9, mfc="none", mec=ORANGE, mew=0.9, ls="none")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"Rabi frequency $\Omega_R/2\pi$ (MHz)", labelpad=1); ax.set_ylabel("laser pump / cap", labelpad=1)
cb = fig.colorbar(pm, cax=cax); cb.ax.minorticks_off(); cb.set_ticks([2, 5, 10, 20, 50]); cb.set_ticklabels(["2", "5", "10", "20", "50"])
cb.set_label(r"$\eta_{\rm CW}$ ($\mu$T/$\sqrt{\rm Hz}$), NV$^-$ only", labelpad=2)
ax.text(0.04, 0.06, "(a)", transform=ax.transAxes, bbox=dict(fc="w", ec="none", pad=1.5, alpha=0.85))

bx = fig.add_axes([0.17, 0.10, 0.70, 0.36])
c = d["cut"]
x = np.array([r["rabi"] for r in c]) / 1e6
bx.plot(x, [r["eta"] * 1e6 for r in cut_nv0], "-o", color=BLUE, ms=3.5, lw=1.4, label=r"CRB, with NV$^0$")
bx.plot(x, [r["eta"] * 1e6 for r in c], "--o", color=BLUE, ms=3.5, lw=1.4, mfc="w", label=r"CRB, NV$^-$ only")
bx.plot(x, [r["slope"] * 1e6 for r in c], ":D", color=GREY, ms=3, lw=1.2, label=r"slope formula, NV$^-$ only")
bx.set_xscale("log"); bx.set_yscale("log")
bx.set_yticks([2, 5, 10, 20, 50]); bx.set_yticklabels(["2", "5", "10", "20", "50"]); bx.minorticks_off()
bx.set_xticks([1, 10]); bx.set_xticklabels(["1", "10"])
bx.set_xlabel(r"Rabi frequency $\Omega_R/2\pi$ (MHz)", labelpad=1)
bx.set_ylabel(r"$\eta_{\rm CW}$ ($\mu$T/$\sqrt{\rm Hz}$)", labelpad=1)
bx.grid(True, lw=0.3, alpha=0.5); bx.legend(fontsize=5.8, frameon=False, loc="upper center")
bx.text(0.04, 0.06, "(b)", transform=bx.transAxes)
OUT_FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
os.makedirs(OUT_FIG, exist_ok=True)
for ext, kw in (("pdf", {}), ("jpg", {"dpi": 600, "pil_kwargs": {"quality": 95}}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(OUT_FIG, "fig1." + ext), **kw)
print([(round(r["rabi"]/1e6,2), round(r["eta"]*1e9)) for r in cut_nv0])
