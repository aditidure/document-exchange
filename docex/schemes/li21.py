"""
li21.py -- the higher-complexity scheme: distributed graph colouring (Li, Gabrys, Farnoud, arXiv:2512.04197;
reference [7] of the report; Corollary 34 of the paper for short substring edits).

Proposition 1 of the report: a valid scheme is a proper colouring of the confusion graph G (vertices {0,1}^n,
x ~ x' iff x' in C(x)). The paper computes the colour of ONE string with a two-round version of Linial's algorithm:

  Cover-free family (Lemma 13): for a prime Q and degree b, string/colour i gets the polynomial g_i over F_Q whose
      coefficients are the base-Q digits of i, and the set F_i = {(alpha, g_i(alpha)) : alpha in F_Q}.
      If b * r < Q, no F_i is covered by r others.
  Witness (Lemma 13):  W(i0, {i1..ir}) = (alpha, g_i0(alpha)) for the smallest alpha where g_i0 differs from every
      g_ij (it exists because prod_j (g_i0 - g_ij) is a nonzero polynomial of degree <= b r < Q).
  Recolouring (Lemma 15): new colour of v = W(old colour of v, old colours of the neighbours of v);
      adjacent vertices get different new colours.
  Two rounds (Lemma 16.2):
      round 1: old colour = the string itself, N1 = 2^n, b1 = ceil(log2 N1) = n,     Q1 = prime > Delta b1
      round 2: old colour = round-1 colour,  N2 = Q1^2, b2 = ceil(log2 N2),           Q2 = prime > Delta b2
      tag = round-2 colour, an integer below Q2^2:  ceil(log2 Q2^2) = 2 log Delta + O(log log n) bits.
  Decoding: list Bob's suspects; they all share y, so they are pairwise adjacent and have distinct colours;
      the one with the received colour is x.

Choices (stated in the report):
  * Delta is the EXACT maximum degree of G, computed once per (n, k, t) by enumerating all 2^n strings;
  * Q is the smallest prime above Delta b (the paper allows the smallest prime power);
  * round-1 colours are cached, because the colour of a string needs the round-1 colours of all its neighbours.
Complexity: one colour needs the round-1 colours of Delta neighbours, each needing that neighbour's confusion ball,
i.e. about Delta^2 ball enumerations, polynomial in n for constant t, k (the paper states O(n^8t) encoding).
"""
import math
from dataclasses import dataclass
from functools import lru_cache

from docex.core.edits import confusion_ball, suspects, all_bits
from docex.schemes.base import Scheme, Tag, DecodingError


def _is_prime(m):
    if m < 2:
        return False
    if m % 2 == 0:
        return m == 2
    i = 3
    while i * i <= m:
        if m % i == 0:
            return False
        i += 2
    return True


def smallest_prime_above(m):
    p = m + 1
    while not _is_prime(p):
        p += 1
    return p


@lru_cache(maxsize=None)
def max_degree(n, k, t):
    """Exact maximum degree Delta of the confusion graph on {0,1}^n (enumerates all 2^n strings once)."""
    return max(len(confusion_ball(x, k, t)) - 1 for x in all_bits(n))


def poly_value(i, Q, alpha):
    """g_i(alpha) over F_Q, the coefficients of g_i being the base-Q digits of i (lowest degree first)."""
    value, power = 0, 1
    while i:
        i, digit = divmod(i, Q)
        value = (value + digit * power) % Q
        power = power * alpha % Q
    return value


def witness(colour, neighbour_colours, Q):
    """Lemma 13: (alpha, g_colour(alpha)) for the smallest alpha with g_colour(alpha) != g_c(alpha) for every
    neighbour colour c, packed as the integer alpha * Q + g_colour(alpha) in [0, Q^2)."""
    others = list(neighbour_colours)
    for alpha in range(Q):
        g = poly_value(colour, Q, alpha)
        if all(poly_value(c, Q, alpha) != g for c in others):
            return alpha * Q + g
    raise RuntimeError("no witness: b r < Q is violated")


@dataclass
class ColouringTag(Tag):
    colour: int
    bits_len: int

    @property
    def parts(self):
        return {"colour": self.bits_len}


class Li21Scheme(Scheme):
    name = "distributed colouring [7]"

    def __init__(self, k, t, n):
        super().__init__(k, t)
        self.n = n
        self.delta = max_degree(n, k, t)
        self.b1 = n
        self.Q1 = smallest_prime_above(self.delta * self.b1)
        self.b2 = math.ceil(math.log2(self.Q1 * self.Q1))
        self.Q2 = smallest_prime_above(self.delta * self.b2)
        self.tag_len = math.ceil(math.log2(self.Q2 * self.Q2))
        self._neighbours = lru_cache(maxsize=None)(self._neighbours_uncached)
        self.round1 = lru_cache(maxsize=None)(self._round1_uncached)

    def _neighbours_uncached(self, x):
        return tuple(z for z in confusion_ball(x, self.k, self.t) if z != x)

    def _round1_uncached(self, x):
        """Round 1: old colour = the string itself (an integer below 2^n)."""
        return witness(int(x, 2), (int(z, 2) for z in self._neighbours(x)), self.Q1)

    def colour(self, x):
        """Round 2: old colour = the round-1 colour."""
        return witness(self.round1(x), (self.round1(z) for z in self._neighbours(x)), self.Q2)

    def encode(self, x, return_info=False):
        if len(x) != self.n:
            raise ValueError(f"this scheme was built for n = {self.n}")
        tag = ColouringTag(self.colour(x), self.tag_len)
        if return_info:
            return tag, {"delta": self.delta, "Q1": self.Q1, "Q2": self.Q2}
        return tag

    def decode(self, y, tag, n):
        hits = [z for z in suspects(y, n, self.k, self.t) if self.colour(z) == tag.colour]
        if len(hits) != 1:
            raise DecodingError(f"{len(hits)} candidates have the received colour")
        return hits[0]
