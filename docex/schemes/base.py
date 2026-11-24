"""
base.py -- the interface shared by every document exchange scheme.

A scheme is built for fixed parameters (t, k) and offers
    encode(x)        -> Tag     Alice: the message sent for her version x
    decode(y, tag, n) -> x      Bob:   recover x (length n) from his version y and the tag
and every Tag reports its length in bits (the redundancy). Experiments loop over schemes through this interface only,
so the worst-case, average-case and higher-complexity schemes are measured in exactly the same way.
"""


class DecodingError(ValueError):
    """Raised when Bob's candidate set does not contain exactly one string matching the tag."""


class Tag:
    """Base class of a transmitted message; subclasses define `parts` (field name -> bits)."""

    @property
    def parts(self):
        raise NotImplementedError

    @property
    def bits(self):
        return sum(self.parts.values())


class Scheme:
    name = "scheme"

    def __init__(self, k, t):
        self.k, self.t = k, t

    def encode(self, x):
        raise NotImplementedError

    def decode(self, y, tag, n):
        raise NotImplementedError

    def __repr__(self):
        return f"{self.name}(k={self.k}, t={self.t})"
