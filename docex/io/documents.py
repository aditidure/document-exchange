"""
documents.py -- binary documents split into fixed-size blocks (report, Sec. III-C).

A document is a text file of '0'/'1' characters (line breaks are ignored when reading the bits). It is cut into
consecutive blocks of b bits; every block is encoded separately and Bob decodes block j of his version with tag j.
Because an insertion or a deletion would shift every later block, the edited version must keep the length, and
each block must be reachable from the original block by at most t length-preserving edits of at most k bits.

  read_bits(path)                         the document as one bit string
  split_blocks(bits, b)                   consecutive b-bit blocks (the length must be a multiple of b)
  min_substitution_edits(x, y, k)         fewest length-preserving k-bit edits turning block x into block y
  classify_blocks(xs, ys, t, k)           per block: 'unchanged', 'fits' or 'violation'
"""


def read_bits(path):
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    bits = "".join(text.split())
    bad = set(bits) - {"0", "1"}
    if bad:
        raise ValueError(f"{path}: only 0 and 1 are allowed, found {sorted(bad)}")
    return bits


def split_blocks(bits, b):
    if len(bits) % b:
        raise ValueError(f"document length {len(bits)} is not a multiple of the block size {b}")
    return [bits[i:i + b] for i in range(0, len(bits), b)]


def min_substitution_edits(x, y, k):
    """Fewest edits, each replacing at most k consecutive bits by as many bits, that turn x into y (|x| = |y|).
    Greedy: start a window of k bits at the first differing position not yet covered."""
    if len(x) != len(y):
        raise ValueError("blocks of different lengths")
    diff = [i for i in range(len(x)) if x[i] != y[i]]
    edits, covered_to = 0, -1
    for i in diff:
        if i > covered_to:
            edits += 1
            covered_to = i + k - 1
    return edits


def classify_blocks(xs, ys, t, k):
    """For each pair of blocks: 'unchanged', 'fits' (at most t edits of at most k bits) or 'violation'."""
    if len(xs) != len(ys):
        raise ValueError(f"the versions have {len(xs)} and {len(ys)} blocks")
    out = []
    for x, y in zip(xs, ys):
        if x == y:
            out.append("unchanged")
        elif len(x) == len(y) and min_substitution_edits(x, y, k) <= t:
            out.append("fits")
        else:
            out.append("violation")
    return out
