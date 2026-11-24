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
            ax.plot([r["n"] for r in rs], [r["bits_mean"] for r in rs], "o-",
                    color=C1, lw=2, ms=6, label="mean (Monte Carlo)")
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
    if not have("average_case_mc.csv"):
        return
    rows = by_tk(read_csv("average_case_mc.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))
    for tk, rs in rows.items():
        t, k = tk
        x = [r["log2n"] for r in rs]
        b, _ = fit_slope(x, [r["s_bar"] for r in rs])
        axes[0].plot(x, [r["s_bar"] for r in rs], "o-", color=TK[tk], lw=2,
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


# ================================================================ comparison (synthetic data)
SCHEME_COLOUR = {"worst-case (Thm 1)": C1, "average-case (Thm 2)": C2, "distributed colouring [7]": C5}


def comparison():
    if not have("compare.csv"):
        return
    groups = {}
    for r in read_csv("compare.csv"):
        groups.setdefault(r["scheme"], []).append(r)

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ns = sorted({r["n"] for rs in groups.values() for r in rs})
    ax.plot(ns, [4 * math.log2(n) for n in ns], "--", color=INK2, lw=1.4, label="4t·log n")
    for name, rs in groups.items():
        ax.plot([r["n"] for r in rs], [r["bits_mean"] for r in rs], "o-", color=SCHEME_COLOUR.get(name, INK), lw=2, ms=6,
                label=name)
    ax.set_xlabel("n (bits per block)"); ax.set_ylabel("tag length (bits)")
    ax.set_title("Redundancy of the three schemes (t = k = 1)")
    ax.legend(loc="center right")
    save(fig, "fig10_compare_redundancy.png")

    fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.6))
    for ax, key, title in ((axes[0], "enc_s", "Encoding time per block"), (axes[1], "dec_s", "Decoding time per block")):
        for name, rs in groups.items():
            xs = [r["n"] for r in rs]
            b, _ = fit_slope([math.log2(x) for x in xs], [math.log2(max(r[key], 1e-7)) for r in rs])
            ax.loglog(xs, [r[key] * 1e3 for r in rs], "o-", color=SCHEME_COLOUR.get(name, INK), lw=2,
                      label=f"{name}: slope {b:.1f}")
        ax.set_title(title); ax.set_xlabel("n (bits per block)"); ax.set_ylabel("milliseconds")
        ax.set_xticks(ns, [str(int(n)) for n in ns]); ax.minorticks_off()
        ax.legend(fontsize=8)
    save(fig, "fig11_compare_cost.png")


# ================================================================ the two hand-edited documents
def documents():
    if not have("documents_blocks.csv"):
        return
    with open(os.path.join(RESULTS, "documents_blocks.csv"), encoding="utf-8") as fh:
        import csv
        rows = list(csv.DictReader(fh))
    docs = sorted({r["document"] for r in rows})
    names = list(SCHEME_COLOUR)
    fig, axes = plt.subplots(len(docs), 1, figsize=(7.2, 3.3 * len(docs)), squeeze=False)
    for ax, doc in zip(axes[:, 0], docs):
        rs = [r for r in rows if r["document"] == doc]
        blocks = sorted({int(r["block"]) for r in rs})
        labels = []
        for j in blocks:
            r = next(q for q in rs if int(q["block"]) == j)
            lx, ly = len(r["original"]), len(r["edited"])
            kind = "–" if r["original"] == r["edited"] else ("ins" if ly > lx else "del" if ly < lx else "flip")
            labels.append(f"{j + 1}\n{kind}")
        w = 0.26
        for i, name in enumerate(names):
            vals = [float(next(q["bits"] for q in rs if int(q["block"]) == j and q["scheme"] == name)) for j in blocks]
            ok = [next(q["recovered"] for q in rs if int(q["block"]) == j and q["scheme"] == name) == "1" for j in blocks]
            ax.bar([j + (i - 1) * w for j in blocks], vals, w, color=SCHEME_COLOUR[name], label=name,
                   edgecolor="white", linewidth=0.8)
            for j, good in zip(blocks, ok):
                if not good:
                    ax.text(j + (i - 1) * w, 2, "✗", ha="center", color=INK)
        ax.axhline(16, color=INK2, ls="--", lw=1.2)
        ax.text(blocks[-1] + 0.45, 17, "block = 16 bits", fontsize=8, color=INK2, ha="right")
        ax.set_xticks(blocks, labels)
        ax.set_ylabel("tag bits"); ax.set_title(f"{doc}: tag per block (line, edit) — every block recovered")
        ax.grid(axis="x", visible=False)
    axes[0, 0].legend(ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.32))
    save(fig, "fig15_documents.png")


# ================================================================ extension: real version-control edits
def extension():
    if not have("ext_fit.csv", "ext_pairs.csv"):
        return
    fit = read_csv("ext_fit.csv")
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    for T, c in zip([1, 2, 3, 5], [C1, C2, C3, C5]):
        ax.semilogx([r["k"] for r in fit], [100 * r[f"t_le_{T}"] for r in fit], "o-", color=c, lw=2, label=f"t ≤ {T}")
    ax.axvline(8, color=INK2, ls=":", lw=1)
    ax.text(8.6, 92, "1 byte", fontsize=8, color=INK2)
    ax.set_xlabel("k (bits per substring edit)"); ax.set_ylabel("changes that fit the model (%)")
    ax.set_title("Do real commits fit 't edits of ≤ k bits'?"); ax.set_ylim(0, 100); ax.legend()
    save(fig, "fig12_real_fit.png")

    pairs = read_csv("ext_pairs.csv")
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 7.0))
    axes[0].hist([min(r["hunks"], 60) for r in pairs], bins=30, color=C1)
    axes[0].set_xlabel("changed regions per file change (capped at 60)"); axes[0].set_ylabel("file changes")
    axes[0].set_title("How many edits does a change make?")
    axes[1].hist([r["max_hunk_bits"] for r in pairs if r["max_hunk_bits"] > 0], bins=[2 ** i for i in range(0, 19)], color=C2)
    axes[1].set_xscale("log"); axes[1].set_xlabel("largest changed region (bits)"); axes[1].set_ylabel("file changes")
    axes[1].set_title("How long is the largest edit?")
    save(fig, "fig13_real_hunks.png")

    if have("ext_tag.csv"):
        tag = [r for r in read_csv("ext_tag.csv") if r["k"] > 0]
        if tag:
            fig, ax = plt.subplots(figsize=(6.6, 4.4))
            allx = []
            for k, c in zip((8, 64, 512), (C3, C4, C2)):
                rs = [r for r in tag if int(r["k"]) == k]
                if rs:
                    ax.loglog([r["n_bits"] for r in rs], [r["est_tag"] for r in rs], "o", color=c, ms=6, label=f"estimated tag, k = {k}")
                    allx += [r["n_bits"] for r in rs]
            ax.loglog([min(allx), max(allx)], [min(allx), max(allx)], "--", color=INK2, lw=1.4, label="n (send the whole file)")
            ax.set_xlabel("file size n (bits)"); ax.set_ylabel("bits")
            ax.set_title("Estimated tag for changes with t ≤ 3 (Lemma 3 bound)"); ax.legend(fontsize=8)
            save(fig, "fig14_real_tag.png")


FIGS = [theorem1, theorem2, comparison, documents, extension]

if __name__ == "__main__":
    for f in FIGS:
        print(f.__name__)
        f()
