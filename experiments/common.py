"""
common.py -- shared helpers for the experiments: paths, Monte Carlo statistics, CSV files, timing.

Monte Carlo: draw x uniformly from {0,1}^n, draw t random k-substring edits, run the scheme; repeat and summarise
  mean_ci(values)       sample mean and 95% half-width 1.96 s / sqrt(N)
  wilson(successes, N)  95% Wilson interval for a success rate (stays inside [0, 1], also at 100%)
  fit_slope(xs, ys)     least-squares slope (used on log-log and log-linear data)
"""
import csv
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)

QUICK = "--quick" in sys.argv


def mean_ci(values):
    v = list(values)
    m = sum(v) / len(v)
    if len(v) < 2:
        return m, 0.0
    s = math.sqrt(sum((a - m) ** 2 for a in v) / (len(v) - 1))
    return m, 1.96 * s / math.sqrt(len(v))


def wilson(successes, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = successes / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, max(0.0, centre - half), min(1.0, centre + half)


def fit_slope(xs, ys):
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    b = sum((a - xm) * (c - ym) for a, c in zip(xs, ys)) / sum((a - xm) ** 2 for a in xs)
    return b, ym - b * xm


def write_csv(name, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    path = os.path.join(RESULTS, name)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    return path


def read_csv(name):
    with open(os.path.join(RESULTS, name), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


class Timer:
    def __enter__(self):
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.s = time.perf_counter() - self.t0
