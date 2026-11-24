"""
density.py -- pattern-dense strings (Definition 4) and the bound of Lemma 4 (report, Sec. II-B2).

The pattern p = 0^k 1^k acts as a fence post: cutting x at the left edge of every occurrence gives "parts".
Two forms of density are used:
  * window form (Definition 4 of the report): every length-delta substring of x contains p;
  * parts form (used by [8], which Lemma 5 relies on): every part has length <= delta.
A window-dense string is parts-dense, so the parts form is the more permissive one; the scheme uses it.

Lemma 4: for uniform x and delta = alpha log n,
    Pr[x not dense] <= (n - delta + 1)(1 - 2^-l)^floor(delta / l),   l = |p| = 2k,
which is O(1/n^3) for alpha = k 2^(2k+3).
"""
import math


def pattern(k):
    return '0' * k + '1' * k


def delta_paper(n, k):
    """delta = k 2^(2k+3) log2 n. Larger than n for n < 256 (k = 1) and n < 2951 (k = 2)."""
    return k * 2 ** (2 * k + 3) * math.log2(n)


def partition(x, k):
    """Cut x at the left edge of every occurrence of 0^k 1^k (the pattern cannot overlap itself).
    Returns (parts, start positions). The first part is the prefix before the first occurrence (possibly empty)."""
    p = pattern(k)
    starts, i = [0], x.find(p)
    while i != -1:
        starts.append(i)
        i = x.find(p, i + len(p))
    ends = starts[1:] + [len(x)]
    return [x[s:e] for s, e in zip(starts, ends)], starts


def is_dense(x, k, delta):
    """Parts form: every part has length <= delta."""
    parts, _ = partition(x, k)
    return max(len(u) for u in parts) <= delta


def is_window_dense(x, k, delta):
    """Window form (Definition 4): every length-delta window contains 0^k 1^k (always true when delta >= n)."""
    p, d = pattern(k), int(delta)
    if d >= len(x):
        return True
    occ = [i for i in range(len(x) - len(p) + 1) if x.startswith(p, i)]
    j = 0
    for s in range(len(x) - d + 1):          # window [s, s+d) needs an occurrence starting in [s, s+d-|p|]
        while j < len(occ) and occ[j] < s:
            j += 1
        if j == len(occ) or occ[j] > s + d - len(p):
            return False
    return True


def lemma4_bound(n, k, delta):
    """Union bound of Lemma 4, capped at 1 (0 when delta >= n: then every string is dense)."""
    if delta >= n:
        return 0.0
    l = 2 * k
    return min(1.0, (n - delta + 1) * (1 - 2.0 ** -l) ** math.floor(delta / l))
