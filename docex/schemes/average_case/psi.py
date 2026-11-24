"""
psi.py -- psi(x) and the one-edit decoder (the "what" half of Lemma 5, report Sec. II-B2 and III-A2).

Lemma 5: for a dense x, (h(x), psi(x)) lets Bob undo ANY single k-substring edit, with |h| <= log n and
|psi| = o(log n). h locates the edit (hash_h.locate); psi recovers the bits.

Our psi(x) = (A_x, S_0, S_1):
  A_x  largest part integer max_i f(u_i): gives Bob a valid bound for locate(); A_x <= 8 delta C (Lemma 7 of [8]),
       so it needs log2(8 C delta) = O(log log n) bits.
  S_c  sum of int(block_j) mod P over the blocks j = c (mod 2), blocks of L bits, P a fixed prime with
       ceil(log2 L) + 12 bits.
Recovery: the located window touches at most two consecutive blocks (one even, one odd). Every other block is
known from y, so subtracting its residue from S_0 or S_1 leaves the residues of the damaged blocks of x; Bob tries
every single edit of the window and keeps the one whose blocks match.

Choices that differ from [8] (stated in the report): the block checksum is a residue modulo a fixed prime instead of
the hash of [4], so recovery holds with high probability rather than always; the block length L is a parameter.
The Theorem 2 decoder uses psi only as a filter, so its correctness never depends on this recovery.
"""
import math

from docex.core.edits import all_bits
from docex.schemes.average_case.hash_h import part_integers, h_value, locate, weights

MARGIN = 12                       # extra prime bits beyond log2 L
DEFAULT_L = {1: 8192, 2: 16384}   # block length per k (from the measured localisation windows)


def _is_prime(m):
    if m < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if m % p == 0:
            return m == p
    d, s = m - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:                                 # deterministic Miller-Rabin for m < 3.3e24
        y = pow(a, d, m)
        if y in (1, m - 1):
            continue
        for _ in range(s - 1):
            y = y * y % m
            if y == m - 1:
                break
        else:
            return False
    return True


def psi_params(n, k, L0=None):
    """Block length L = min(L0, n) and the prime P >= 2^(ceil(log2 L) + MARGIN)."""
    L0 = DEFAULT_L.get(k, 16384) if L0 is None else L0
    L = max(1, min(L0, n))
    P = 1 << (max(1, math.ceil(math.log2(L))) + MARGIN)
    while not _is_prime(P):
        P += 1
    return L, P


def amax_bits(n, k):
    """Field size of A_x: log2(8 C delta) with the paper's delta = k 2^(2k+3) log2 n."""
    _, C = weights(k)
    delta = k * 2 ** (2 * k + 3) * max(1.0, math.log2(n))
    return math.ceil(math.log2(8 * C * delta))


def psi_bits(n, k, L0=None):
    L, P = psi_params(n, k, L0)
    sums = 1 if math.ceil(n / L) == 1 else 2
    return amax_bits(n, k) + sums * math.ceil(math.log2(P))


def psi_value(x, k, L0=None):
    L, P = psi_params(len(x), k, L0)
    z, _, _ = part_integers(x, k)
    res = [int(x[j:j + L], 2) % P for j in range(0, len(x), L)]
    return (max(z), sum(res[0::2]) % P, sum(res[1::2]) % P)


def _one_edit_variants(s, k, m):
    """All length-m strings within ONE k-substring edit of s."""
    out, ls = set(), len(s)
    for i in range(ls + 1):
        for len_u in range(min(k, ls - i) + 1):
            len_v = m - (ls - len_u)
            if 0 <= len_v <= k:
                for v in all_bits(len_v):
                    out.add(s[:i] + v + s[i + len_u:])
    return out


def decode_one_edit(y, hx, psix, n, k, L0=None, return_info=False):
    """Lemma 5: recover a dense x (length n) from y in B_1(x), h(x) and psi(x). Raises ValueError if it cannot."""
    L, P = psi_params(n, k, L0)
    amax, S0, S1 = psix
    info = {"window": 0, "candidates": 1}
    if len(y) == n and h_value(y, k, n) == hx and psi_value(y, k, L0) == psix:
        return (y, info) if return_info else y               # the edit changed nothing visible
    w, _, _ = part_integers(y, k)
    loc = locate(y, hx, n, k, max(amax, max(w)))
    if loc is None:
        raise ValueError("localisation failed")
    lo, hi_y = loc[2], loc[3]
    hi_x = n - (len(y) - hi_y)                                # window end in x (the suffix after it is shared)
    if hi_x < lo:
        raise ValueError("inconsistent window")
    info["window"] = hi_y - lo
    jlo = lo // L
    jhi = max(jlo, (hi_x - 1) // L) if hi_x > lo else jlo
    if jhi - jlo > 1:
        raise ValueError(f"window spans {jhi - jlo + 1} blocks (more than 2)")
    shift = len(y) - n
    target = [S0, S1]
    for j in range(math.ceil(n / L)):                         # subtract the blocks that the edit did not touch
        if jlo <= j <= jhi:
            continue
        a, b = j * L, min((j + 1) * L, n)
        block = y[a:b] if b <= lo else y[a + shift:b + shift]
        target[j % 2] = (target[j % 2] - int(block, 2) % P) % P
    ra, rb = jlo * L, min((jhi + 1) * L, n)
    left, right, seg = y[ra:lo], y[hi_x + shift:rb + shift], y[lo:hi_y]
    cands = _one_edit_variants(seg, k, hi_x - lo)
    info["candidates"] = len(cands)
    hits = []
    for c in cands:
        region = left + c + right
        if all(int(region[j * L - ra:min((j + 1) * L, n) - ra], 2) % P == target[j % 2] for j in range(jlo, jhi + 1)):
            hits.append(c)
    if len(hits) != 1:
        full = [z for z in (y[:lo] + c + y[hi_y:] for c in hits) if h_value(z, k, n) == hx]
        if len(full) != 1:
            raise ValueError(f"{len(hits)} candidates match psi")
        return (full[0], info) if return_info else full[0]
    xhat = y[:lo] + hits[0] + y[hi_y:]
    return (xhat, info) if return_info else xhat
