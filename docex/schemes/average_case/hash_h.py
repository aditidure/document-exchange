"""
hash_h.py -- the hash h of Li et al. [8] and Bob's localisation (the "where" half of Lemma 5).

1. Cut x into parts at every 0^k 1^k (density.partition).
2. Map each part u to one integer  f(u) = d(u) + B n_1k(u) + C |u|,   B = 3 2^k,  C = 40 k 2^k:
       d(u)     parity of every k-th bit inside u, read as a k-bit number   (sees substitutions)
       n_1k(u)  number of overlapping runs 1^k in u                         (sees a shifted pattern)
       |u|      length of u                                                 (sees insertions and deletions)
   The weights B and C keep the three effects from cancelling.
3. h(x) = ( VT(z) mod 2000 C n,  sum(z) mod 10 C k,  #parts mod 5 ),  z = (f(u_1), f(u_2), ...),
   VT(z) = sum_i i z_i  (Varshamov-Tenengolts weighted sum).
One edit changes only a few consecutive entries of z, so comparing VT(z), sum(z) and #parts with Bob's own string
reveals WHERE the edit is (locate). |h| = log n + O(1) bits.
"""
import math
import re

from docex.schemes.average_case.density import partition


def weights(k):
    """(B, C) of the part map."""
    return 3 * 2 ** k, 40 * k * 2 ** k


def _parity_digest(u, k):
    return int(''.join(str(u[j::k].count('1') & 1) for j in range(k)), 2)


def _runs_1k(u, k):
    return sum(max(0, len(r) - k + 1) for r in re.findall('1+', u))


def part_integers(x, k):
    """z = (f(u_1), f(u_2), ...), with the part start positions and the parts themselves."""
    B, C = weights(k)
    parts, starts = partition(x, k)
    return [_parity_digest(u, k) + B * _runs_1k(u, k) + C * len(u) for u in parts], starts, parts


def vt(z):
    return sum((i + 1) * v for i, v in enumerate(z))


def h_from_integers(z, n, k):
    _, C = weights(k)
    return (vt(z) % (2000 * C * n), sum(z) % (10 * C * k), len(z) % 5)


def h_value(x, k, n=None):
    n = len(x) if n is None else n
    z, _, _ = part_integers(x, k)
    return h_from_integers(z, n, k)


def h_bits(n, k):
    _, C = weights(k)
    return math.ceil(math.log2(2000 * C * n)) + math.ceil(math.log2(10 * C * k)) + math.ceil(math.log2(5))


def locate(y, hx, n, k, A, K=3):
    """Bob's localisation from y and h(x): the interval of y (in parts and in bits) that must contain the edit.
    A bounds the part integers of x and y. Returns (part_lo, part_hi, bit_lo, bit_hi, #parts of y) or None."""
    _, C = weights(k)
    w, starts, _ = part_integers(y, k)
    vt_mod, sum_mod, np_mod = hx
    m = 2000 * C * n
    d_vt = (vt_mod - vt(w)) % m                    # VT(z) - VT(w), recovered from its residue
    if d_vt > m // 2:
        d_vt -= m
    d_len = (np_mod - len(w)) % 5                  # change in the number of parts (|.| <= 2)
    if d_len > 2:
        d_len -= 5
    mod_sum = 10 * C * k
    d_sum = (sum_mod - sum(w)) % mod_sum           # change in the sum of part integers
    if d_sum > mod_sum // 2:
        d_sum -= mod_sum
    L = len(w)
    tail = [0] * (L + 2)                           # tail[v] = sum of w_u for u >= v
    for v in range(L, 0, -1):
        tail[v] = tail[v + 1] + w[v - 1]
    tol = 3 * K * K * A                            # an edit at part j gives |d_vt - (d_len tail[j] + j d_sum)| <= tol
    cand = [v for v in range(1, L + 2) if abs(d_len * tail[v] + v * d_sum - d_vt) <= tol]
    if not cand:
        return None
    lo, hi = min(cand), max(cand)
    bit_lo = starts[lo - 1] if lo - 1 < len(starts) else len(y)
    bit_hi = starts[hi] if hi < len(starts) else len(y)
    return lo, hi, bit_lo, bit_hi, L
