"""
syndrome.py -- syndrome compression (Lemma 1, report Sec. II-B).

Alice does not send the long label f(x); she sends f(x) mod a for a small modulus a that still separates x from
every string in its confusion ball:
        f(x) mod a != f(y) mod a      for all y in C(x), y != x.
a depends on x, so it is sent too: the tag is (f(x) mod a, a).

From the proof of Lemma 1: a fails exactly when it divides some difference d_y = f(x) - f(y); each d_y has at most
R - 1 distinct prime factors, so fewer than M R primes are bad (M = |C(x)| - 1) and a_x <= 2 M R ln(M R).
"""
import math


def is_valid_modulus(a, f_x, neighbour_labels):
    """True if f(x) mod a differs from f(y) mod a for every neighbour label f(y)."""
    return all((f_x - f_y) % a for f_y in neighbour_labels)


def find_modulus(f_x, neighbour_labels):
    """Smallest a >= 2 separating f(x) from every neighbour label (Lemma 1).

    Each candidate a is checked against ALL neighbours. The reference paper raises a "at no extra cost" whenever a
    collision appears while labelling; that shortcut would need the earlier neighbours to be re-checked against the
    new a, so we do the full check (correct, possibly slower)."""
    diffs = [f_x - f_y for f_y in neighbour_labels]
    a = 2
    while not all(d % a for d in diffs):
        a += 1
    return a


def tag_bits(r, a):
    """Bits needed to send (r, a) in fixed-width fields: r < a takes ceil(log2 a) bits, a takes ceil(log2(a+1))."""
    return math.ceil(math.log2(a)) + math.ceil(math.log2(a + 1))


def lemma1_bound(M, R):
    """The proof's guarantee a_x <= 2 M R ln(M R) (valid for M R >= 6)."""
    K = max(6, M * R)
    return 2 * K * math.log(K)


def distinct_prime_factors(d):
    """Distinct prime factors of |d| (used by the tests of the bad-prime counting step)."""
    d, p, out = abs(d), 2, []
    while p * p <= d:
        if d % p == 0:
            out.append(p)
            while d % p == 0:
                d //= p
        p += 1
    if d > 1:
        out.append(d)
    return out
