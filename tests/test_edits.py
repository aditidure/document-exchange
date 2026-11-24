"""
test_edits.py -- checks of the edit model (docex/core/edits.py) against the definitions.

Run from the project root:  python3 -m unittest tests.test_edits -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402

rng = random.Random(0)


class OneEdit(unittest.TestCase):
    def test_paper_example_1(self):
        # [23], Example 1: x -> y by a 5-substring edit, u = 0010 replaced by v = 10011 at position 5 (1-indexed)
        self.assertEqual(edits.apply_edit("11100010100110", 4, 4, "10011"), "111010011100110")

    def test_out_of_range_edit_is_rejected(self):
        with self.assertRaises(ValueError):
            edits.apply_edit("0101", 3, 2, "1")

    def test_one_edit_ball_by_hand(self):
        # from "0" with k = 1: keep it, delete it, substitute it, or insert one bit before/after it
        self.assertEqual(edits.one_edit_ball("0", 1), {"0", "", "1", "00", "10", "01"})

    def test_one_edit_count_bound(self):
        # Lemma 3 proof: |B_1(s)| <= (m + 1)(k + 1)(2^(k+1) - 1) for a string of length m
        for m, k in [(6, 1), (6, 2), (5, 3)]:
            for _ in range(5):
                s = edits.random_string(m, rng)
                self.assertLessEqual(len(edits.one_edit_ball(s, k)), (m + 1) * (k + 1) * (2 ** (k + 1) - 1))


class Balls(unittest.TestCase):
    def test_ball_contains_centre_and_grows_with_t(self):
        s = "011010"
        b1, b2 = edits.edit_ball(s, 1, 1), edits.edit_ball(s, 1, 2)
        self.assertIn(s, b1)
        self.assertTrue(b1 <= b2)

    def test_random_edits_stay_in_the_ball(self):
        for _ in range(30):
            x = edits.random_string(8, rng)
            y = edits.random_edits(x, 2, 2, rng)
            self.assertIn(y, edits.edit_ball(x, 2, 2))

    def test_confusion_ball_equals_definition(self):
        # Eq. (1): C_t(x) = B_2t(x) ∩ {0,1}^n must equal {x' : B_t(x) ∩ B_t(x') ≠ ∅}, for every x
        for n, k, t in [(6, 1, 1), (5, 2, 1), (5, 1, 2)]:
            balls = {x: edits.edit_ball(x, k, t) for x in edits.all_bits(n)}
            for x in edits.all_bits(n):
                by_definition = {x2 for x2 in edits.all_bits(n) if balls[x] & balls[x2]}
                self.assertEqual(edits.confusion_ball(x, k, t), by_definition, (n, k, t, x))

    def test_confusion_is_symmetric(self):
        for x in edits.all_bits(7):
            for z in edits.confusion_ball(x, 1, 1):
                self.assertIn(x, edits.confusion_ball(z, 1, 1))


class BobsSuspects(unittest.TestCase):
    def test_x_is_a_suspect_and_suspects_are_confusable(self):
        # Theorem 1 correctness: x ∈ B_t(y) ∩ {0,1}^n (reversibility) and that set ⊆ C_t(x)
        for n, k, t in [(10, 1, 1), (8, 2, 1), (8, 1, 2)]:
            for _ in range(15):
                x = edits.random_string(n, rng)
                y = edits.random_edits(x, k, t, rng)
                s = edits.suspects(y, n, k, t)
                self.assertIn(x, s)
                self.assertTrue(s <= edits.confusion_ball(x, k, t))


if __name__ == "__main__":
    unittest.main()
