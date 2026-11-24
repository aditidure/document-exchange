"""
test_density.py -- pattern density (docex/schemes/average_case/density.py).

Run from the project root:  python3 -m unittest tests.test_density -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes.average_case import density  # noqa: E402

rng = random.Random(3)


class Density(unittest.TestCase):
    def test_partition_example_of_li22(self):
        # example of [8] with k = 3, p = 000111
        x = "01001" + "000111" + "100" + "000111" + "101001" + "000111" + "1"
        parts, starts = density.partition(x, 3)
        self.assertEqual(parts, ["01001", "000111100", "000111101001", "0001111"])
        self.assertEqual("".join(parts), x)

    def test_window_form_by_hand(self):
        # k = 1, p = 01, delta = 4 (the presentation's examples)
        self.assertTrue(density.is_window_dense("011010010110", 1, 4))
        self.assertFalse(density.is_window_dense("111100001011", 1, 4))

    def test_window_dense_implies_parts_dense(self):
        for _ in range(300):
            x = edits.random_string(40, rng)
            if density.is_window_dense(x, 1, 6):
                self.assertTrue(density.is_dense(x, 1, 6))

    def test_window_form_matches_brute_force(self):
        p = density.pattern(1)
        for _ in range(200):
            x, d = edits.random_string(30, rng), rng.randint(3, 12)
            brute = all(p in x[s:s + d] for s in range(len(x) - d + 1))
            self.assertEqual(density.is_window_dense(x, 1, d), brute)

    def test_lemma4_bound_at_paper_delta(self):
        # paper's delta: the bound is far below 1/n^3
        for n, k in [(2 ** 10, 1), (2 ** 16, 1)]:
            self.assertLess(density.lemma4_bound(n, k, density.delta_paper(n, k)), n ** -3)


if __name__ == "__main__":
    unittest.main()
