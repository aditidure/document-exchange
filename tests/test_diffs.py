"""
test_diffs.py -- real edits through the k-substring edit model (docex/io/diffs.py) and the bounds (core/bounds.py).

Run from the project root:  python3 -m unittest tests.test_diffs -v
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import bounds  # noqa: E402
from docex.io import diffs  # noqa: E402


class Diffs(unittest.TestCase):
    def test_one_character_substitution(self):
        old, new = b"value = 10\n", b"value = 11\n"          # '0' (0x30) -> '1' (0x31): one bit
        hs = diffs.bit_hunks(old, new, diffs.byte_hunks(old, new))
        self.assertEqual(len(hs), 1)
        self.assertEqual(diffs.hunk_size(hs[0]), 1)
        self.assertEqual(diffs.t_min(hs, 1), 1)

    def test_inserted_character(self):
        old, new = b"abc\n", b"abXc\n"
        hs = diffs.bit_hunks(old, new, diffs.byte_hunks(old, new))
        self.assertEqual(diffs.t_min(hs, 8), 1)               # one inserted byte = one 8-bit edit

    def test_t_min_groups_close_hunks(self):
        hs = [(0, 1, 0, 1), (3, 4, 3, 4), (100, 101, 100, 101)]
        self.assertEqual(diffs.t_min(hs, 1), 3)               # three 1-bit edits
        self.assertEqual(diffs.t_min(hs, 4), 2)               # the first two fit one 4-bit window
        self.assertEqual(diffs.t_min([], 4), 0)

    def test_line_changes(self):
        old = b"a\nb\nc\n"
        new = b"a\nB\nc\nd\n"
        pairs, added, removed = diffs.line_changes(old, new)
        self.assertEqual((pairs, added, removed), ([(b"b", b"B")], 1, 0))


class Bounds(unittest.TestCase):
    def test_leading_terms(self):
        self.assertEqual(bounds.theorem1(16, 1), 16)
        self.assertEqual(bounds.theorem2(16, 2), 28)

    def test_lemma3_bound_value(self):
        # prod (n + ik) for i = 0..2 with n = 8, k = 1: 8 * 9 * 10 = 720; times 1^2 * 3^2 = 9
        self.assertAlmostEqual(bounds.lemma3_log_ball(8, 1, 1), math.log2(720 * 9))


if __name__ == "__main__":
    unittest.main()
