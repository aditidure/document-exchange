"""
diffs.py -- real edits seen through the k-substring edit model (report, Sec. IV).

Given the old and new version of a file (bytes):
  byte_hunks(old, new)      changed regions (old_start, old_end, new_start, new_end) in BYTES: a line-level diff,
                            then a byte-level diff inside each changed block (blocks above BYTE_CAP stay coarse)
  bit_hunks(old, new, hb)   the same regions in BITS, after trimming the bits that agree at both ends
  t_min(hunks, k)           fewest k-substring edits that turn old into new: a group of consecutive hunks spanning
                            s_old bits of old and s_new bits of new needs ceil(max(s_old, s_new) / k) edits; the best
                            grouping is found by dynamic programming (groups of at most MAX_GROUP hunks)
  line_changes(old, new)    the edit seen line by line: changed line pairs, and lines added or removed
"""
from difflib import SequenceMatcher

BYTE_CAP = 1500
MAX_GROUP = 40


def bits_of(data):
    return ''.join(format(c, '08b') for c in data)


def byte_hunks(old, new):
    xl, yl = old.splitlines(keepends=True), new.splitlines(keepends=True)
    xo, yo = [0], [0]
    for line in xl:
        xo.append(xo[-1] + len(line))
    for line in yl:
        yo.append(yo[-1] + len(line))
    out = []
    for tag, i1, i2, j1, j2 in SequenceMatcher(None, xl, yl, autojunk=False).get_opcodes():
        if tag == 'equal':
            continue
        xa, xb, ya, yb = xo[i1], xo[i2], yo[j1], yo[j2]
        if xb - xa <= BYTE_CAP and yb - ya <= BYTE_CAP:
            for t2, a1, a2, b1, b2 in SequenceMatcher(None, old[xa:xb], new[ya:yb], autojunk=False).get_opcodes():
                if t2 != 'equal':
                    out.append((xa + a1, xa + a2, ya + b1, ya + b2))
        else:
            out.append((xa, xb, ya, yb))
    return out


def bit_hunks(old, new, hb):
    res = []
    for xa, xb, ya, yb in hb:
        bx, by = bits_of(old[xa:xb]), bits_of(new[ya:yb])
        m, p = min(len(bx), len(by)), 0
        while p < m and bx[p] == by[p]:
            p += 1
        s = 0
        while s < m - p and bx[-1 - s] == by[-1 - s]:
            s += 1
        res.append((xa * 8 + p, xb * 8 - s, ya * 8 + p, yb * 8 - s))
    return [h for h in res if h[1] > h[0] or h[3] > h[2]]


def hunk_size(h):
    return max(h[1] - h[0], h[3] - h[2])


def t_min(hunks, k):
    m = len(hunks)
    if m == 0:
        return 0
    best = [0] + [10 ** 9] * m
    for j in range(1, m + 1):
        for i in range(j - 1, max(-1, j - 1 - MAX_GROUP), -1):
            span = max(hunks[j - 1][1] - hunks[i][0], hunks[j - 1][3] - hunks[i][2], 1)
            best[j] = min(best[j], best[i] + -(-span // k))
    return best[m]


def line_changes(old, new):
    """(changed line pairs, lines added, lines removed). A 'replace' block with as many old as new lines is read as
    pairs of changed lines; any surplus counts as added or removed lines."""
    xl, yl = old.splitlines(), new.splitlines()
    pairs, added, removed = [], 0, 0
    for tag, i1, i2, j1, j2 in SequenceMatcher(None, xl, yl, autojunk=False).get_opcodes():
        if tag == 'insert':
            added += j2 - j1
        elif tag == 'delete':
            removed += i2 - i1
        elif tag == 'replace':
            n = min(i2 - i1, j2 - j1)
            pairs += list(zip(xl[i1:i1 + n], yl[j1:j1 + n]))
            added += (j2 - j1) - n
            removed += (i2 - i1) - n
    return pairs, added, removed
