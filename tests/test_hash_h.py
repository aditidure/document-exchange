"""
test_hash_h.py -- the hash h and Bob's localisation (docex/schemes/average_case/hash_h.py).

Run from the project root:  python3 -m unittest tests.test_hash_h -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes.average_case import hash_h  # noqa: E402

rng = random.Random(4)


def edit_window(x, y):
    """Bits of y that differ from x after removing the common prefix and suffix."""
    p = 0
    while p < min(len(x), len(y)) and x[p] == y[p]:
        p += 1
    s = 0
    while s < min(len(x), len(y)) - p and x[-1 - s] == y[-1 - s]:
        s += 1
    return p, len(y) - s


class HashH(unittest.TestCase):
    def test_runs_example_of_li22(self):
        # [8]: "111100111" with k = 2 has 5 overlapping runs 11
        self.assertEqual(hash_h._runs_1k("111100111", 2), 5)

    def test_length_is_log_n_plus_constant(self):
        # |h| - log2 n is the same constant for every n (here k = 1)
        diffs = {hash_h.h_bits(2 ** e, 1) - e for e in (8, 12, 16, 20)}
        self.assertEqual(len(diffs), 1)

    def test_h_detects_a_change(self):
        # an edit that changes x changes h (on these random cases)
        for _ in range(200):
            x = edits.random_string(256, rng)
            y = edits.random_edit(x, 1, rng)
            if y != x and len(y) == len(x):
                self.assertNotEqual(hash_h.h_value(x, 1), hash_h.h_value(y, 1, 256))

    def test_located_window_contains_the_edit(self):
        for k, n in [(1, 1024), (2, 1024)]:
            for _ in range(60):
                x = edits.random_string(n, rng)
                y = edits.random_edit(x, k, rng)
                if y == x:
                    continue
                z, _, _ = hash_h.part_integers(x, k)
                w, _, _ = hash_h.part_integers(y, k)
                loc = hash_h.locate(y, hash_h.h_value(x, k, n), n, k, A=max(max(z), max(w)))
                self.assertIsNotNone(loc)
                lo, hi = edit_window(x, y)
                self.assertTrue(loc[2] <= lo and hi <= loc[3], (loc, lo, hi))


if __name__ == "__main__":
    unittest.main()
