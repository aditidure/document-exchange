r"""
worst_case.py -- Theorem 1: the low-complexity worst-case scheme (report, Sec. II-B1 and III-A1).

    Alice:  C(x) = B_2t(x) ∩ {0,1}^n                                   (edits.confusion_ball, Lemma 3)
            a_x  = smallest a with f(x) mod a != f(y) mod a, all y in C(x)\{x}   (syndrome.find_modulus, Lemma 1)
            Enc(x) = (f(x) mod a_x, a_x)
    Bob:    S = B_t(y) ∩ {0,1}^n                                       (edits.suspects)
            output the unique z in S with f(z) mod a_x = r
Correctness: every z in S shares the descendant y with x, so S ⊆ C(x); a_x separates x from all of C(x).
Encoding enumerates |C(x)| = O(n^2t) strings and labels each in O(n): O(n^(2t+1)).
Decoding enumerates |B_t(y)| = O(n^t) strings and labels each in O(n): O(n^(t+1)).
"""
from dataclasses import dataclass

from docex.core.edits import confusion_ball, suspects
from docex.core.labeling import label_int
from docex.core.syndrome import find_modulus, tag_bits
from docex.schemes.base import Scheme, Tag, DecodingError


@dataclass
class WorstCaseTag(Tag):
    r: int          # f(x) mod a
    a: int          # the modulus a_x

    @property
    def parts(self):
        return {"ra": tag_bits(self.r, self.a)}


class WorstCaseScheme(Scheme):
    name = "worst-case (Thm 1)"

    def __init__(self, k, t, label=label_int):
        super().__init__(k, t)
        self.f = label

    def encode(self, x, return_info=False):
        ball = confusion_ball(x, self.k, self.t)
        fx = self.f(x)
        a = find_modulus(fx, [self.f(z) for z in ball if z != x])
        tag = WorstCaseTag(fx % a, a)
        if return_info:
            return tag, {"ball": len(ball)}
        return tag

    def decode(self, y, tag, n):
        hits = [z for z in suspects(y, n, self.k, self.t) if self.f(z) % tag.a == tag.r]
        if len(hits) != 1:
            raise DecodingError(f"{len(hits)} candidates match the tag")
        return hits[0]
