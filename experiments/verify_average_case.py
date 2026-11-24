"""
verify_average_case.py -- verification of Theorem 2 and of Lemmas 4-6 (report, Sec. III-A2).

Run from the project root:   python3 experiments/verify_average_case.py [--quick]
Writes results/lemma4_density.csv, results/lemma5_one_edit.csv, results/average_case_mc.csv.

  Lemma 4   empirical Pr[x not dense] (window form) with a Wilson interval, against the union bound
  Lemma 5   one edit undone from (y, h(x), psi(x)) alone: success rate and located window
  Theorem 2 average tag length s̄(n) with its breakdown (flag, h, psi, (r,a)); |C|, |C~_h|, |C~| (Lemma 6);
            the Theorem 1 tag on the same strings; decoding success
A bound "(4t-1) log n + o(log n)" cannot be checked at one n, so the checks are:
  C1 every tag decodes   C2 (r,a) <= (4t-2) log n   C3 |h| - log n is constant
  C4 fitted slope of s̄ against log n <= 4t-1         C5 [s̄ - (4t-1) log n] / log n decreases
"""
import math
import random

from common import QUICK, Timer, mean_ci, wilson, fit_slope, write_csv
from docex.core import edits
from docex.schemes.worst_case import WorstCaseScheme
from docex.schemes.average_case import density, hash_h, psi
from docex.schemes.average_case.scheme import AverageCaseScheme

rng = random.Random(2)


def lemma4():
    rows, trials = [], (300 if QUICK else 2000)
    cases = [(1, 256, [8, 10, 12, 14, 16, 18, 20, 24, 28]), (1, 1024, [10, 12, 14, 16, 20, 24, 28, 32]),
             (2, 1024, [32, 40, 48, 56, 64, 80, 96])]
    for k, n, deltas in cases:
        if QUICK and n > 256:
            continue
        xs = [edits.random_string(n, rng) for _ in range(trials)]
        for d in deltas:
            bad = sum(not density.is_window_dense(x, k, d) for x in xs)
            p, lo, hi = wilson(bad, trials)
            rows.append(dict(k=k, n=n, delta=d, trials=trials, nondense=p, ci_lo=lo, ci_hi=hi,
                             bound=density.lemma4_bound(n, k, d)))
            print(f"  Lemma 4 k={k} n={n} δ={d:3d}: Pr[not dense] {p:.4f} [{lo:.4f}, {hi:.4f}], bound {rows[-1]['bound']:.4f}",
                  flush=True)
    d = density.delta_paper(1024, 1)
    print(f"  at the paper's δ = {d:.0f} (n = 1024, k = 1) the bound is {density.lemma4_bound(1024, 1, d):.1e}")
    return rows


def lemma5():
    rows = []
    grid = [(1, [2 ** 12, 2 ** 14, 2 ** 16, 2 ** 18]), (2, [2 ** 14, 2 ** 16, 2 ** 18])]
    if QUICK:
        grid = [(1, [2 ** 12, 2 ** 14]), (2, [2 ** 14])]
    for k, ns in grid:
        for n in ns:
            trials = 10 if QUICK else (60 if n <= 2 ** 16 else 25)
            ok, widths, fails = 0, [], 0
            for _ in range(trials):
                x = edits.random_string(n, rng)
                y = edits.random_edit(x, k, rng)
                while y == x:
                    y = edits.random_edit(x, k, rng)
                try:
                    xh, info = psi.decode_one_edit(y, hash_h.h_value(x, k, n), psi.psi_value(x, k), n, k, return_info=True)
                    ok += (xh == x)
                    widths.append(info["window"])
                except ValueError:
                    fails += 1
            p, lo, hi = wilson(ok, trials)
            med = sorted(widths)[len(widths) // 2] if widths else float("nan")
            rows.append(dict(k=k, n=n, log2n=math.log2(n), trials=trials, recovered=ok, success=p, ci_lo=lo, ci_hi=hi,
                             window_median=med, window_frac=med / n, failures=fails,
                             h_bits=hash_h.h_bits(n, k), psi_bits=psi.psi_bits(n, k)))
            print(f"  Lemma 5 k={k} n=2^{int(math.log2(n))}: recovered {ok}/{trials} [{lo:.2f}, {hi:.2f}], window ≈ {med} bits"
                  f" ({med / n:.2f} n), failures {fails}", flush=True)
    return rows


GRID = {(1, 1): [(16, 40), (32, 40), (64, 30), (128, 20), (256, 10)],
        (2, 1): [(8, 30), (12, 30), (16, 20), (20, 10), (24, 6)],
        (1, 2): [(16, 30), (32, 30), (64, 20), (128, 10)]}
if QUICK:
    GRID = {(1, 1): [(16, 10), (32, 10), (64, 6), (128, 3)], (2, 1): [(8, 8), (12, 6), (16, 4)]}


def theorem2():
    rows = []
    for (t, k), grid in GRID.items():
        avg, worst = AverageCaseScheme(k, t), WorstCaseScheme(k, t)
        for n, nx in grid:
            per, ok, tot = [], 0, 0
            with Timer() as tm:
                for _ in range(nx):
                    x = edits.random_string(n, rng)
                    tag, info = avg.encode(x, return_info=True)
                    wtag, winfo = worst.encode(x, return_info=True)
                    for _ in range(3 if QUICK else 5):
                        ok += avg.decode(edits.random_edits(x, k, t, rng), tag, n) == x
                        tot += 1
                    per.append(dict(bits=tag.bits, **{f"bits_{q}": v for q, v in tag.parts.items()}, dense=info["dense"],
                                    ball_h=info["ball_h"], ball_restricted=info["ball_restricted"],
                                    ball_full=winfo["ball"], thm1_bits=wtag.bits))
            m = lambda key: mean_ci([q[key] for q in per])[0]
            ln = math.log2(n)
            s_bar, s_ci = mean_ci([q["bits"] for q in per])
            row = dict(t=t, k=k, n=n, log2n=ln, strings=nx, frac_dense=m("dense"), s_bar=s_bar, s_ci=s_ci,
                       bits_flag=m("bits_flag"), bits_h=m("bits_h"), bits_psi=m("bits_psi"), bits_ra=m("bits_ra"),
                       thm1_mean=m("thm1_bits"), ball_full=m("ball_full"), ball_h=m("ball_h"),
                       ball_restricted=m("ball_restricted"), bound_thm2=(4 * t - 1) * ln, bound_ra=(4 * t - 2) * ln,
                       decoded=ok, trials=tot, seconds=tm.s)
            row["excess_per_logn"] = (row["s_bar"] - row["bound_thm2"]) / ln
            row["h_minus_logn"] = row["bits_h"] - ln
            rows.append(row)
            print(f"  Thm 2 t={t} k={k} n={n:4d}: s̄ {s_bar:5.1f} = 1 + h {row['bits_h']:.0f} + ψ {row['bits_psi']:.0f} + (r,a) "
                  f"{row['bits_ra']:.1f} | Thm 1 {row['thm1_mean']:.1f} | |C| {row['ball_full']:.0f} → |C~_h| "
                  f"{row['ball_h']:.2f} → |C~| {row['ball_restricted']:.2f} | decoded {ok}/{tot} | {tm.s:.0f}s", flush=True)
    return rows


def checks(rows):
    print("\n  checks")
    for (t, k) in GRID:
        rs = [r for r in rows if (r["t"], r["k"]) == (t, k)]
        b, _ = fit_slope([r["log2n"] for r in rs], [r["s_bar"] for r in rs])
        ep = [r["excess_per_logn"] for r in rs]
        hm = [r["h_minus_logn"] for r in rs]
        for name, good in [("C1 every tag decodes", all(r["decoded"] == r["trials"] for r in rs)),
                           ("C2 (r,a) <= (4t-2) log n", all(r["bits_ra"] <= r["bound_ra"] for r in rs)),
                           ("C3 |h| - log n constant", max(hm) - min(hm) <= 2),
                           (f"C4 slope {b:.2f} <= 4t-1 = {4 * t - 1}", b <= 4 * t - 1),
                           ("C5 excess/log n decreasing " + " -> ".join(f"{e:.1f}" for e in ep),
                            all(e2 <= e1 for e1, e2 in zip(ep, ep[1:])))]:
            print(f"    t={t} k={k} [{'PASS' if good else 'FAIL'}] {name}")


if __name__ == "__main__":
    print("Lemma 4")
    write_csv("lemma4_density.csv", lemma4())
    print("Lemma 5")
    write_csv("lemma5_one_edit.csv", lemma5())
    print("Theorem 2")
    rows = theorem2()
    write_csv("average_case_mc.csv", rows)
    checks(rows)
