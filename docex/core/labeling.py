"""
labeling.py -- the labeling function f of Lemma 2 (report, Sec. II-B).

What Lemma 1 needs from f:
    (1) f(c) != f(c') for every c' confusable with c     -> labels separate confusable pairs
    (2) labels are integers below 2^R, R = O(log n log log n)

Lemma 2 ([7] of the reference paper) gives such an f with
    R <= [(tau^2 + 1)(2 tau^2 + 1) + 2 tau^2 (tau - 1)] log n + o(log n),   tau = 2tk.
The construction is only cited in the reference paper, and it is long: 657 log n bits for t = 1, k = 2, which is
more than n at every n we can simulate.

Choice made here: f(x) = x read as a binary number (R = n).
  * injective, so (1) holds for ALL pairs, in particular for confusable ones -> decoding is always correct;
  * violates (2), so the GUARANTEE of Lemma 1 on a_x is not inherited; the tag lengths are measured instead.
"""
import math


def label_int(x):
    """f(x) = x as an integer (0 for the empty string)."""
    return int(x, 2) if x else 0


def lemma2_label_bits(n, t, k):
    """Leading term of R in Lemma 2 with tau = 2tk (the o(log n) term is not specified)."""
    tau = 2 * t * k
    return ((tau ** 2 + 1) * (2 * tau ** 2 + 1) + 2 * tau ** 2 * (tau - 1)) * math.log2(n)
