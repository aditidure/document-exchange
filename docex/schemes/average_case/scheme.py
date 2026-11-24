"""
scheme.py -- Theorem 2: the average-case scheme (report, Sec. II-B2 and III-A2).

    x not dense:  Enc(x) = (0, f(x) mod a, a)               a separates x from the full ball C(x)   (Theorem 1)
    x dense:      Enc(x) = (1, h(x), psi(x), f(x) mod a, a)  a separates x from the restricted ball
                  C~(x) = { z in C(x) : z dense, h(z) = h(x), psi(z) = psi(x) }      (Definition 5, Lemma 6)
Bob reads the flag. In the dense branch he keeps the suspects that are dense with the same h and psi; they all lie
in C~(x), so the residue picks out x. The first bit is the partition index of Definition 3.

delta defaults to the paper's value k 2^(2k+3) log n capped at n: below n = 256 (k = 1) every string is then dense,
which is the regime of Lemma 4, so the dense branch is the one measured.
"""
from dataclasses import dataclass
from typing import Optional, Tuple

from docex.core.edits import confusion_ball, suspects
from docex.core.labeling import label_int
from docex.core.syndrome import find_modulus, tag_bits
from docex.schemes.base import Scheme, Tag, DecodingError
from docex.schemes.average_case.density import is_dense, delta_paper
from docex.schemes.average_case.hash_h import part_integers, h_value, h_from_integers, h_bits
from docex.schemes.average_case.psi import psi_value, psi_bits


@dataclass
class AverageCaseTag(Tag):
    flag: int                     # 1 = dense branch, 0 = Theorem 1 branch
    h: Optional[Tuple]
    psi: Optional[Tuple]
    r: int
    a: int
    n: int
    k: int

    @property
    def parts(self):
        if self.flag:
            return {"flag": 1, "h": h_bits(self.n, self.k), "psi": psi_bits(self.n, self.k), "ra": tag_bits(self.r, self.a)}
        return {"flag": 1, "h": 0, "psi": 0, "ra": tag_bits(self.r, self.a)}


def scheme_delta(n, k):
    return min(n, delta_paper(n, k))


class AverageCaseScheme(Scheme):
    name = "average-case (Thm 2)"

    def __init__(self, k, t, delta=None, label=label_int):
        super().__init__(k, t)
        self.delta, self.f = delta, label

    def _delta(self, n):
        return scheme_delta(n, self.k) if self.delta is None else self.delta

    def restricted_ball(self, x):
        """(C~(x), C~_h(x)): confusable strings that are dense with the same h, and among those the ones with the same
        psi. One partition per candidate: density and h are read from the same part integers."""
        n, k, delta = len(x), self.k, self._delta(len(x))
        hx, px = h_value(x, k, n), psi_value(x, k)
        same_h, same_h_psi = [], []
        for z in confusion_ball(x, k, self.t):
            zi, _, parts = part_integers(z, k)
            if max(len(u) for u in parts) > delta or h_from_integers(zi, n, k) != hx:
                continue
            same_h.append(z)
            if psi_value(z, k) == px:
                same_h_psi.append(z)
        return same_h_psi, same_h

    def encode(self, x, return_info=False):
        n, fx = len(x), self.f(x)
        if is_dense(x, self.k, self._delta(n)):
            ball, ball_h = self.restricted_ball(x)
            a = find_modulus(fx, [self.f(z) for z in ball if z != x])
            tag = AverageCaseTag(1, h_value(x, self.k, n), psi_value(x, self.k), fx % a, a, n, self.k)
            info = {"dense": True, "ball_h": len(ball_h), "ball_restricted": len(ball)}
        else:
            ball = confusion_ball(x, self.k, self.t)
            a = find_modulus(fx, [self.f(z) for z in ball if z != x])
            tag = AverageCaseTag(0, None, None, fx % a, a, n, self.k)
            info = {"dense": False, "ball_h": len(ball), "ball_restricted": len(ball)}
        return (tag, info) if return_info else tag

    def decode(self, y, tag, n):
        cands = suspects(y, n, self.k, self.t)
        if tag.flag:
            delta = self._delta(n)
            cands = [z for z in cands if is_dense(z, self.k, delta) and h_value(z, self.k, n) == tag.h
                     and psi_value(z, self.k) == tag.psi]
        hits = [z for z in cands if self.f(z) % tag.a == tag.r]
        if len(hits) != 1:
            raise DecodingError(f"{len(hits)} candidates match the tag")
        return hits[0]
