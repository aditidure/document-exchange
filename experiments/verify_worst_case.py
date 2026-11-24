"""
verify_worst_case.py -- verification of Theorem 1 (report, Sec. III-A1).

Run from the project root:   python3 experiments/verify_worst_case.py [--quick]
Writes results/worst_case_mc.csv and results/worst_case_exhaustive.csv.

  Monte Carlo   random x, t random edits per trial: tag bits (mean ± CI, max), |C(x)|, a_x, encode/decode time,
                decoding success with a Wilson interval
  Exhaustive    ALL 2^n strings at small n: the true worst case max_x |Enc(x)|, and min/mean/max |C(x)|
Checks printed at the end:
  every tag decodes;  max |Enc(x)| against 4t log n;  slope of log|C| vs log n (Lemma 3: 2t);
  bits added per doubling of n (<= 4t);  encode/decode time slopes (<= 2t+1, <= t+1).
"""
import math
import random

from common import QUICK, Timer, mean_ci, wilson, fit_slope, write_csv
from docex.core import edits
from docex.schemes.worst_case import WorstCaseScheme

rng = random.Random(1)

MC = {  # (t, k): [(n, strings, edit patterns per string)]
    (1, 1): [(16, 40, 5), (32, 40, 5), (64, 30, 5), (128, 20, 5), (256, 10, 5)],
    (1, 2): [(16, 40, 5), (24, 40, 5), (32, 30, 5), (48, 20, 5), (64, 10, 5)],
    (2, 1): [(8, 40, 5), (12, 40, 5), (16, 30, 5), (20, 15, 5), (24, 8, 5)],
}
EXHAUSTIVE = {(1, 1): [8, 10, 12, 14, 16], (1, 2): [8, 10, 12], (2, 1): [6, 8, 10]}
if QUICK:
    MC = {(1, 1): [(16, 10, 3), (32, 10, 3), (64, 6, 3)], (2, 1): [(8, 10, 3), (12, 6, 3)]}
    EXHAUSTIVE = {(1, 1): [8, 10, 12], (2, 1): [6, 8]}


def bound(n, t):
    return 4 * t * math.log2(n)


def monte_carlo():
    rows = []
    for (t, k), grid in MC.items():
        scheme = WorstCaseScheme(k, t)
        for n, nx, nd in grid:
            bits, balls, mods, tenc, tdec, ok, tot = [], [], [], [], [], 0, 0
            for _ in range(nx):
                x = edits.random_string(n, rng)
                with Timer() as tm:
                    tag, info = scheme.encode(x, return_info=True)
                tenc.append(tm.s)
                bits.append(tag.bits); balls.append(info["ball"]); mods.append(tag.a)
                for _ in range(nd):
                    y = edits.random_edits(x, k, t, rng)
                    with Timer() as tm:
                        ok += scheme.decode(y, tag, n) == x
                    tdec.append(tm.s)
                    tot += 1
            mb, cb = mean_ci(bits)
            p, lo, hi = wilson(ok, tot)
            row = dict(t=t, k=k, n=n, strings=nx, bits_mean=mb, bits_ci=cb, bits_max=max(bits), bound=bound(n, t),
                       ball_mean=mean_ci(balls)[0], a_mean=mean_ci(mods)[0], enc_s=mean_ci(tenc)[0],
                       dec_s=mean_ci(tdec)[0], decoded=ok, trials=tot, success_lo=lo)
            rows.append(row)
            print(f"  MC  t={t} k={k} n={n:4d}: bits {mb:5.1f} ± {cb:.1f} (max {row['bits_max']}) vs 4t·log n {row['bound']:5.1f}"
                  f" | |C| {row['ball_mean']:8.0f} | a {row['a_mean']:7.1f} | enc {row['enc_s'] * 1e3:8.1f} ms"
                  f" dec {row['dec_s'] * 1e3:6.2f} ms | decoded {ok}/{tot}", flush=True)
    return rows


def exhaustive():
    rows = []
    for (t, k), ns in EXHAUSTIVE.items():
        scheme = WorstCaseScheme(k, t)
        for n in ns:
            bits, balls = [], []
            with Timer() as tm:
                for x in edits.all_bits(n):
                    tag, info = scheme.encode(x, return_info=True)
                    bits.append(tag.bits); balls.append(info["ball"])
            worst = edits.all_bits(n)[max(range(len(balls)), key=balls.__getitem__)]
            row = dict(t=t, k=k, n=n, strings=2 ** n, bits_mean=sum(bits) / len(bits), bits_max=max(bits),
                       bound=bound(n, t), ball_min=min(balls), ball_mean=sum(balls) / len(balls), ball_max=max(balls),
                       largest_ball_x=worst, seconds=tm.s)
            rows.append(row)
            print(f"  ALL t={t} k={k} n={n:3d}: max bits {row['bits_max']} vs {row['bound']:.1f} | |C| min/mean/max "
                  f"{row['ball_min']}/{row['ball_mean']:.1f}/{row['ball_max']} (largest at {worst}) | {tm.s:.1f}s", flush=True)
    return rows


def checks(mc, ex):
    print("\n  checks")
    print(f"    every tag decoded: {sum(r['decoded'] for r in mc)}/{sum(r['trials'] for r in mc)}")
    for r in ex:
        if r["bits_max"] > r["bound"]:
            print(f"    t={r['t']} k={r['k']} n={r['n']}: worst tag {r['bits_max']} exceeds the leading term "
                  f"{r['bound']:.1f} by {r['bits_max'] - r['bound']:.1f} bits (allowed by the o(log n) term)")
    for (t, k) in MC:
        rs = [r for r in mc if (r["t"], r["k"]) == (t, k)]
        lx = [math.log2(r["n"]) for r in rs]
        print(f"    t={t} k={k}: slope log|C| {fit_slope(lx, [math.log2(r['ball_mean']) for r in rs])[0]:.2f} (2t={2 * t})"
              f" | bits/doubling {fit_slope(lx, [r['bits_mean'] for r in rs])[0]:.2f} (<= {4 * t})"
              f" | enc slope {fit_slope(lx, [math.log2(r['enc_s']) for r in rs])[0]:.2f} (<= {2 * t + 1})"
              f" | dec slope {fit_slope(lx, [math.log2(r['dec_s']) for r in rs])[0]:.2f} (<= {t + 1})")


if __name__ == "__main__":
    print("Theorem 1, Monte Carlo")
    mc = monte_carlo()
    write_csv("worst_case_mc.csv", mc)
    print("Theorem 1, exhaustive")
    ex = exhaustive()
    write_csv("worst_case_exhaustive.csv", ex)
    checks(mc, ex)
