"""
test_average_case.py -- Theorem 2 scheme (docex/schemes/average_case/scheme.py).

Run from the project root:  python3 -m unittest tests.test_average_case -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes.average_case.scheme import AverageCaseScheme  # noqa: E402

rng = random.Random(6)


class AverageCase(unittest.TestCase):
    def test_round_trip_both_branches(self):
        # a small delta makes some strings non-dense, so both branches are exercised
        for n, k, t, delta in [(16, 1, 1, 6), (12, 1, 2, 5)]:
            s = AverageCaseScheme(k, t, delta=delta)
            seen = set()
            for _ in range(14):
                x = edits.random_string(n, rng)
                tag, info = s.encode(x, return_info=True)
                seen.add(info["dense"])
                for _ in range(4):
                    self.assertEqual(s.decode(edits.random_edits(x, k, t, rng), tag, n), x)
            self.assertEqual(seen, {True, False})

    def test_paper_delta_makes_every_short_string_dense(self):
        s = AverageCaseScheme(1, 1)
        for _ in range(10):
            self.assertTrue(s.encode(edits.random_string(32, rng), return_info=True)[1]["dense"])

    def test_restricted_ball_is_inside_confusion_ball_and_contains_x(self):
        s = AverageCaseScheme(1, 1)
        for _ in range(5):
            x = edits.random_string(32, rng)
            ball, ball_h = s.restricted_ball(x)
            self.assertIn(x, ball)
            self.assertTrue(set(ball) <= set(ball_h) <= edits.confusion_ball(x, 1, 1))

    def test_survivors_lie_in_restricted_ball(self):
        # Bob's suspects that pass the filters are a subset of C~(x) (correctness step of Theorem 2)
        s = AverageCaseScheme(1, 1)
        for _ in range(5):
            x = edits.random_string(32, rng)
            tag = s.encode(x)
            y = edits.random_edits(x, 1, 1, rng)
            ball = set(s.restricted_ball(x)[0])
            from docex.schemes.average_case.density import is_dense
            from docex.schemes.average_case.hash_h import h_value
            from docex.schemes.average_case.psi import psi_value
            survivors = {z for z in edits.suspects(y, 32, 1, 1) if is_dense(z, 1, 32) and h_value(z, 1, 32) == tag.h
                         and psi_value(z, 1) == tag.psi}
            self.assertIn(x, survivors)
            self.assertTrue(survivors <= ball)


if __name__ == "__main__":
    unittest.main()
