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


# ================================================================ Theorem 2
def theorem2():
    from docex.schemes.average_case import density, hash_h, psi
    if have("lemma4_density.csv"):
        fig, ax = plt.subplots(figsize=(6.6, 4.4))
        groups = {}
        for r in read_csv("lemma4_density.csv"):
            groups.setdefault((int(r["k"]), int(r["n"])), []).append(r)
        for (k, n), rs, c in zip(groups, groups.values(), [C1, C2, C3]):
            d = [r["delta"] for r in rs]
            ax.errorbar(d, [r["nondense"] for r in rs], yerr=[[max(0, r["nondense"] - r["ci_lo"]) for r in rs],
                        [max(0, r["ci_hi"] - r["nondense"]) for r in rs]], fmt="o-", color=c, lw=2, capsize=3,
                        label=f"measured, k={k}, n={n}")
            ax.plot(d, [r["bound"] for r in rs], "--", color=c, lw=1.4, label=f"Lemma 4 bound, k={k}, n={n}")
        ax.set_xlabel("δ (window length)"); ax.set_ylabel("Pr[x not (p, δ)-dense]")
        ax.set_title("Lemma 4: non-dense strings are rare (bound valid but loose)")
        ax.legend(fontsize=7.8)
        save(fig, "fig04_lemma4_density.png")

    if have("lemma5_one_edit.csv"):
        fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.4))
        groups = {}
        for r in read_csv("lemma5_one_edit.csv"):
            groups.setdefault(int(r["k"]), []).append(r)
        for k, rs in groups.items():
            c = C1 if k == 1 else C2
            x = [r["log2n"] for r in rs]
            axes[0].errorbar(x, [100 * r["success"] for r in rs], yerr=[[max(0, 100 * (r["success"] - r["ci_lo"])) for r in rs],
                             [max(0, 100 * (r["ci_hi"] - r["success"])) for r in rs]], fmt="o-", color=c, lw=2, capsize=3,
                             label=f"k = {k}")
            axes[1].plot(x, [100 * min(1, r["window_frac"]) for r in rs], "o-", color=c, lw=2, label=f"k = {k}")
        axes[0].set_title("Lemma 5: x recovered from (y, h, ψ) after one edit"); axes[0].set_ylabel("recovered (%)")
        axes[1].set_title("Window located by h, as a share of x"); axes[1].set_ylabel("median window / n (%)")
        for ax in axes:
            ax.set_xlabel("log₂ n"); ax.set_ylim(0, 105); ax.legend()
        save(fig, "fig05_lemma5_one_edit.png")

    if not have("average_case_mc.csv"):
        return
    rows = by_tk(read_csv("average_case_mc.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    for tk, rs in rows.items():
        t, k = tk
        x = [r["log2n"] for r in rs]
        b, _ = fit_slope(x, [r["s_bar"] for r in rs])
        axes[0].errorbar(x, [r["s_bar"] for r in rs], yerr=[r["s_ci"] for r in rs], fmt="o-", color=TK[tk], lw=2, capsize=3,
                         label=f"t={t}, k={k}: slope {b:.2f} (≤ {4 * t - 1})")
        axes[0].plot([x[0], x[-1]], [rs[0]["s_bar"], rs[0]["s_bar"] + (4 * t - 1) * (x[-1] - x[0])], ":", color=TK[tk])
        axes[1].plot(x, [r["excess_per_logn"] for r in rs], "o-", color=TK[tk], lw=2, label=f"t={t}, k={k}")
    axes[0].set_title("Average tag s̄ vs log n (dotted: slope 4t−1)"); axes[0].set_ylabel("s̄ = E|Enc(x)| (bits)")
    axes[1].set_title("Remainder is o(log n): [s̄ − (4t−1)·log n] / log n falls")
    for ax in axes:
        ax.set_xlabel("log₂ n"); ax.legend()
    save(fig, "fig06_thm2_length.png")

    rs = rows.get((1, 1), [])
    if rs:
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        ns = [r["n"] for r in rs]
        ax.loglog(ns, [r["ball_full"] for r in rs], "o-", color=C1, lw=2, label="|C(x)|: full confusion ball")
        ax.loglog(ns, [r["ball_h"] for r in rs], "o-", color=C2, lw=2, label="|C̃_h(x)|: dense, same h")
        ax.loglog(ns, [r["ball_restricted"] for r in rs], "o-", color=C3, lw=2, label="|C̃(x)|: dense, same h and ψ")
        ax.loglog(ns, ns, "--", color=INK2, lw=1.2, label="n^(2t−1) = n (Lemma 6 bound)")
        ax.set_xticks(ns, [str(int(n)) for n in ns]); ax.minorticks_off()
        ax.set_xlabel("n"); ax.set_ylabel("strings Bob must rule out")
        ax.set_title("Lemma 6: the filters shrink the confusion ball (t = k = 1)")
        ax.legend()
        save(fig, "fig07_thm2_ball.png")

        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        xs = list(range(len(rs)))
        ax.bar([i - .2 for i in xs], [r["thm1_mean"] for r in rs], .38, color=C1, label="Theorem 1 tag (r, a)")
        bottom = [0.0] * len(rs)
        for key, col, lab in (("bits_ra", "#9ec5f4", "Theorem 2: (r, a)"), ("bits_flag", C4, "flag"),
                              ("bits_h", C2, "h(x)"), ("bits_psi", C3, "ψ(x)")):
            v = [r[key] for r in rs]
            ax.bar([i + .2 for i in xs], v, .38, bottom=bottom, color=col, label=lab, edgecolor="white", linewidth=1)
            bottom = [b + q for b, q in zip(bottom, v)]
        ax.set_xticks(xs, [f"n={int(r['n'])}" for r in rs]); ax.grid(axis="x", visible=False)
        ax.set_ylabel("bits"); ax.set_title("Where the bits go (t = k = 1)")
        ax.legend(ncol=3, fontsize=8, loc="upper left"); ax.set_ylim(0, max(bottom) * 1.35)
        save(fig, "fig08_thm2_breakdown.png")

    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    e = list(range(4, 161, 2))
    model = []
    for q in e:
        n = 2.0 ** q
        p = density.lemma4_bound(n, 1, density.delta_paper(n, 1))
        dense = 1 + hash_h.h_bits(n, 1) + psi.psi_bits(2 ** q, 1) + 2 * q
        model.append(p * (1 + 4 * q) + (1 - p) * dense)
    worst = [4 * q for q in e]
    ax.plot(e, worst, color=C1, lw=2, label="Theorem 1: 4t·log n")
    ax.plot(e, model, color=C2, lw=2, label="Theorem 2 with our |h|, |ψ|")
    ax.plot(e, [3 * q for q in e], "--", color=C3, lw=2, label="Theorem 2 as stated: (4t−1)·log n")
    cross = next(q for q, a, b in zip(e, model, worst) if a < b)
    ax.axvline(cross, color=INK2, lw=1, ls=":")
    ax.annotate(f"break-even n ≈ 2^{cross}", (cross, 4 * cross), xytext=(cross + 6, 4 * cross - 150), fontsize=9.5,
                arrowprops=dict(arrowstyle="-", color=INK2))
    ax.set_xlabel("log₂ n"); ax.set_ylabel("expected tag length (bits)")
    ax.set_title("Theorem 2 overtakes Theorem 1 only after its constant is paid (t = k = 1)")
    ax.legend(loc="upper left")
    save(fig, "fig09_thm2_break_even.png")


FIGS = [theorem1, theorem2]

if __name__ == "__main__":
    for f in FIGS:
        print(f.__name__)
        f()
