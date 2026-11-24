"""
test_worst_case.py -- Theorem 1 scheme (docex/schemes/worst_case.py).

Run from the project root:  python3 -m unittest tests.test_worst_case -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes.worst_case import WorstCaseScheme  # noqa: E402

rng = random.Random(2)


class WorstCase(unittest.TestCase):
    def test_slides_example(self):
        # running example of the slides: x = 10110, t = k = 1 -> tag (9, 13); Bob holds y = 1010
        s = WorstCaseScheme(k=1, t=1)
        tag = s.encode("10110")
        self.assertEqual((tag.r, tag.a), (9, 13))
        self.assertEqual(s.decode("1010", tag, 5), "10110")

    def test_round_trip_random(self):
        for n, k, t in [(12, 1, 1), (12, 2, 1), (9, 1, 2)]:
            s = WorstCaseScheme(k, t)
            for _ in range(10):
                x = edits.random_string(n, rng)
                tag = s.encode(x)
                for _ in range(5):
                    self.assertEqual(s.decode(edits.random_edits(x, k, t, rng), tag, n), x)

    def test_round_trip_every_string(self):
        # exhaustive: every x of length 8 and every y one edit away decodes (t = k = 1)
        s = WorstCaseScheme(1, 1)
        for x in edits.all_bits(8):
            tag = s.encode(x)
            for y in edits.edit_ball(x, 1, 1):
                self.assertEqual(s.decode(y, tag, 8), x)

    def test_bits_counted(self):
        tag = WorstCaseScheme(1, 1).encode("10110")
        self.assertEqual(tag.bits, 8)


if __name__ == "__main__":
    unittest.main()
