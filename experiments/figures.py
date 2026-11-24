"""
figures.py -- the report's figures, drawn only from results/*.csv (needs matplotlib).

Run from the project root after the experiments:   python experiments/figures.py
Writes figures/*.png. A figure whose CSV file is missing is skipped.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import FIGURES, RESULTS, read_csv, fit_slope  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#6250d6"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
TK = {(1, 1): C1, (1, 2): C3, (2, 1): C2}
plt.rcParams.update({"font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.titlesize": 11,
                     "axes.titlelocation": "left", "legend.frameon": False, "legend.fontsize": 8.5,
                     "axes.axisbelow": True, "savefig.dpi": 160})


def have(*names):
    missing = [n for n in names if not os.path.exists(os.path.join(RESULTS, n))]
    if missing:
        print("  skipped (missing " + ", ".join(missing) + ")")
    return not missing


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, name))
    plt.close(fig)
    print("  wrote figures/" + name)


def by_tk(rows):
    out = {}
    for r in rows:
        out.setdefault((int(r["t"]), int(r["k"])), []).append(r)
    return out


# ================================================================ Theorem 1
def theorem1():
    if not have("worst_case_mc.csv", "worst_case_exhaustive.csv"):
        return
    mc, ex = by_tk(read_csv("worst_case_mc.csv")), by_tk(read_csv("worst_case_exhaustive.csv"))

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, tk in zip(axes, [(1, 1), (1, 2), (2, 1)]):
        t, k = tk
        rs, es = mc.get(tk, []), ex.get(tk, [])
        ns = sorted({r["n"] for r in rs} | {r["n"] for r in es})
        ax.plot(ns, [4 * t * math.log2(n) for n in ns], "--", color=INK2, lw=1.6, label="4t·log n")
        if es:
            ax.plot([r["n"] for r in es], [r["bits_max"] for r in es], "s", color=C2, ms=7, label="max over all x (exhaustive)")
        if rs:
            ax.errorbar([r["n"] for r in rs], [r["bits_mean"] for r in rs], yerr=[r["bits_ci"] for r in rs], fmt="o-",
                        color=C1, lw=2, ms=6, capsize=3, label="mean ± 95% CI (Monte Carlo)")
            ax.plot([r["n"] for r in rs], [r["bits_max"] for r in rs], "^", color=C3, ms=7, label="max over samples")
        ax.set_xscale("log", base=2)
        ax.set_title(f"t = {t}, k = {k}")
        ax.set_xlabel("n (bits)")
    axes[0].set_ylabel("tag length |Enc(x)| (bits)")
    axes[0].legend(loc="upper left")
    fig.suptitle("Theorem 1: measured tag length vs the leading term 4t·log n", x=0.01, ha="left", fontweight="bold")
    save(fig, "fig01_thm1_length.png")

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for tk, rs in mc.items():
        t, k = tk
        ns = [r["n"] for r in rs]
        b, a0 = fit_slope([math.log2(n) for n in ns], [math.log2(r["ball_mean"]) for r in rs])
        ax.loglog(ns, [r["ball_mean"] for r in rs], "o-", color=TK[tk], lw=2, label=f"t={t}, k={k}: slope {b:.2f} (2t = {2 * t})")
        ax.loglog(ns, [rs[0]["ball_mean"] * (n / ns[0]) ** (2 * t) for n in ns], ":", color=TK[tk], lw=1.3)
    ax.set_xlabel("n"); ax.set_ylabel("|C(x)| (mean over random x)")
    ax.set_title("Lemma 3: |C(x)| grows like n^(2t) (dotted: slope 2t)")
    ax.legend()
    save(fig, "fig02_thm1_ball.png")

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.4))
    for tk, rs in mc.items():
        t, k = tk
        ns = [r["n"] for r in rs]
        for ax, key, theory in ((axes[0], "enc_s", 2 * t + 1), (axes[1], "dec_s", t + 1)):
            b, _ = fit_slope([math.log2(n) for n in ns], [math.log2(r[key]) for r in rs])
            ax.loglog(ns, [r[key] * 1e3 for r in rs], "o-", color=TK[tk], lw=2, label=f"t={t}, k={k}: slope {b:.2f} (≤ {theory})")
    axes[0].set_title("Encoding time per string, O(n^(2t+1))")
    axes[1].set_title("Decoding time per string, O(n^(t+1))")
    for ax in axes:
        ax.set_xlabel("n"); ax.set_ylabel("milliseconds"); ax.legend()
    save(fig, "fig03_thm1_time.png")


FIGS = [theorem1]

if __name__ == "__main__":
    for f in FIGS:
        print(f.__name__)
        f()
