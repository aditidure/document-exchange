"""
edits.py -- the edit model (report, Sec. II-A: Definitions 1-2 and Eq. (1)).

Strings are Python str over '0'/'1'.

Definition 1  A k-substring edit deletes a block u of x (|u| <= k) and inserts a block v (|v| <= k) at the same
              position. The edit ball B_t(x) is every string reachable from x by AT MOST t such edits.
Definition 2  x and x' (same length n) are confusable if B_t(x) and B_t(x') intersect; the confusion ball C_t(x)
              is the set of all such x'.
Eq. (1)       C_t(x) = B_2t(x) restricted to length n.
              An edit u -> v is undone by the edit v -> u (also a k-substring edit). So a common descendant z of x
              and x' gives a path of <= 2t edits from x to x', and conversely the string after the first t edits of
              such a path is a common descendant.

Functions
  all_bits(L)                 every binary string of length L (the possible inserted blocks v)
  apply_edit(x, i, len_u, v)  one k-substring edit
  one_edit_ball(s, k)         B_1(s)
  edit_ball(s, k, t)          B_t(s), built one edit at a time
  confusion_ball(x, k, t)     C_t(x), via Eq. (1)
  suspects(y, n, k, t)        Bob's candidate set B_t(y) ∩ {0,1}^n
  random_edit / random_edits  the synthetic channel used in the experiments
"""
from functools import lru_cache
from itertools import product


@lru_cache(maxsize=None)
def all_bits(length):
    """All 2^length binary strings of the given length (cached: the same lengths are needed over and over)."""
    return tuple(''.join(bits) for bits in product('01', repeat=length))


def apply_edit(x, i, len_u, v):
    """One k-substring edit: delete x[i : i+len_u] and insert v at position i (0-indexed)."""
    if not (0 <= i <= len(x) and 0 <= len_u and i + len_u <= len(x)):
        raise ValueError(f"edit out of range: i={i}, |u|={len_u}, |x|={len(x)}")
    return x[:i] + v + x[i + len_u:]


def one_edit_ball(s, k, lengths=None):
    """B_1(s): every string reachable from s by ONE k-substring edit (s itself included, via u = v).

    An edit is fixed by its position i, the deleted length |u| in 0..k and the inserted block v (|v| in 0..k).
    `lengths`, if given, keeps only results whose length is in that set (saves memory on the last layer)."""
    n, out = len(s), set()
    for i in range(n + 1):
        for len_u in range(min(k, n - i) + 1):
            before, after = s[:i], s[i + len_u:]
            for len_v in range(k + 1):
                if lengths is not None and (n - len_u + len_v) not in lengths:
                    continue
                for v in all_bits(len_v):
                    out.add(before + v + after)
    return out


def edit_ball(s, k, t, final_len=None):
    """B_t(s): at most t k-substring edits, applied one layer at a time.

    With final_len, only strings of that length are returned. Intermediate strings may have any length
    (a deletion followed by an insertion can return to the original length)."""
    layer = {s}
    for r in range(t):
        last = (r == t - 1)
        keep = {final_len} if (last and final_len is not None) else None
        new = set()
        for z in layer:
            new |= one_edit_ball(z, k, lengths=keep)
        layer |= new                      # "at most t": strings reached with fewer edits stay in
    if final_len is not None:
        layer = {z for z in layer if len(z) == final_len}
    return layer


def confusion_ball(x, k, t):
    """C_t(x) = B_2t(x) ∩ {0,1}^n, by Eq. (1). Includes x itself."""
    return edit_ball(x, k, 2 * t, final_len=len(x))


def suspects(y, n, k, t):
    """Bob's candidates: every length-n string within t edits of y, i.e. B_t(y) ∩ {0,1}^n."""
    return edit_ball(y, k, t, final_len=n)


# ---------------------------------------------------------------- the synthetic channel
def random_edit(x, k, rng):
    """One random k-substring edit: random |u| in 0..k, random position, random |v| in 0..k, random bits of v."""
    len_u = rng.randint(0, min(k, len(x)))
    i = rng.randint(0, len(x) - len_u)
    len_v = rng.randint(0, k)
    v = ''.join(rng.choice('01') for _ in range(len_v))
    return apply_edit(x, i, len_u, v)


def random_edits(x, k, t, rng):
    """t successive random k-substring edits: the version y held by Bob."""
    for _ in range(t):
        x = random_edit(x, k, rng)
    return x


def random_string(n, rng):
    """A uniformly random n-bit string."""
    return format(rng.getrandbits(n), f'0{n}b') if n > 0 else ''
