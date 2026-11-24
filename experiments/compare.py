"""
compare.py -- comparison of redundancy and complexity on synthetic data (report, Sec. III-D).

Run from the project root:   python3 experiments/compare.py [--quick]
Writes results/compare.csv.

For t = k = 1 and block lengths n = 8 ... 16, the three schemes run through the common interface on the SAME random
strings and the SAME random edits:
    worst-case (Theorem 1), average-case (Theorem 2), distributed colouring [7]
Recorded per scheme and n: tag bits (mean, max), encoding and decoding time per string, decoding success.
"""
import math
import random

from common import QUICK, Timer, mean_ci, write_csv
from docex.core import edits
from docex.schemes.worst_case import WorstCaseScheme
from docex.schemes.average_case.scheme import AverageCaseScheme
from docex.schemes.li21 import Li21Scheme

rng = random.Random(11)
T, K = 1, 1
NS = [8, 10, 12, 14, 16] if not QUICK else [8, 10]
STRINGS, PATTERNS = (15, 3) if not QUICK else (4, 2)


def run():
    rows = []
    for n in NS:
        with Timer() as setup:
            schemes = [WorstCaseScheme(K, T), AverageCaseScheme(K, T), Li21Scheme(K, T, n)]
        xs = [edits.random_string(n, rng) for _ in range(STRINGS)]
        ys = [[edits.random_edits(x, K, T, rng) for _ in range(PATTERNS)] for x in xs]
        for s in schemes:
            bits, tenc, tdec, ok, tot = [], [], [], 0, 0
            for x, yl in zip(xs, ys):
                with Timer() as tm:
                    tag = s.encode(x)
                tenc.append(tm.s)
                bits.append(tag.bits)
                for y in yl:
                    with Timer() as tm:
                        ok += s.decode(y, tag, n) == x
                    tdec.append(tm.s)
                    tot += 1
            row = dict(scheme=s.name, n=n, t=T, k=K, bits_mean=mean_ci(bits)[0], bits_max=max(bits),
                       enc_s=mean_ci(tenc)[0], dec_s=mean_ci(tdec)[0], decoded=ok, trials=tot,
                       bound_4tlogn=4 * T * math.log2(n),
                       setup_s=setup.s if isinstance(s, Li21Scheme) else 0.0,
                       delta=s.delta if isinstance(s, Li21Scheme) else "")
            rows.append(row)
            print(f"  n={n:2d} {s.name:26s}: bits {row['bits_mean']:5.1f} (max {row['bits_max']}) | enc {row['enc_s'] * 1e3:8.2f} ms"
                  f" | dec {row['dec_s'] * 1e3:8.2f} ms | decoded {ok}/{tot}", flush=True)
    return rows


if __name__ == "__main__":
    print(f"Comparison on synthetic data (t = {T}, k = {K})")
    write_csv("compare.csv", run())
