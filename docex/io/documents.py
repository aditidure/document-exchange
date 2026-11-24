"""
documents.py -- binary documents split into blocks at line breaks (report, Sec. III-C).

A document is a text file of '0'/'1' characters, one block per line. Every block is encoded separately and Bob
decodes line j of his version with tag j. Because the blocks are delimited by line breaks rather than by bit
positions, an insertion or a deletion only changes the length of its own line, and every other block still lines up.
Each edited line must be reachable from the original line by at most t k-substring edits (insertions, deletions and
substitutions allowed); the original block length n is known to Bob.

  read_blocks(path)                  the document's blocks (one per non-empty line)
  classify_blocks(xs, ys, k, t)      per block: 'unchanged', 'fits' (y in B_t(x)) or 'violation'
"""
from docex.core.edits import edit_ball


def read_blocks(path):
    with open(path, encoding="utf-8") as fh:
        lines = [line.strip() for line in fh]
    blocks = [line for line in lines if line]
    for j, block in enumerate(blocks):
        bad = set(block) - {"0", "1"}
        if bad:
            raise ValueError(f"{path}, line {j + 1}: only 0 and 1 are allowed, found {sorted(bad)}")
    return blocks


def classify_blocks(xs, ys, k, t):
    """For each pair of blocks: 'unchanged', 'fits' (the edited block is within t k-substring edits of the
    original, i.e. in B_t(x)) or 'violation'."""
    if len(xs) != len(ys):
        raise ValueError(f"the versions have {len(xs)} and {len(ys)} lines: line breaks must not be added or removed")
    out = []
    for x, y in zip(xs, ys):
        if x == y:
            out.append("unchanged")
        elif y in edit_ball(x, k, t):
            out.append("fits")
        else:
            out.append("violation")
    return out
