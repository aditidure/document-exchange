"""
test_li21.py -- the distributed-colouring scheme of [7] (docex/schemes/li21.py).

Run from the project root:  python3 -m unittest tests.test_li21 -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes import li21  # noqa: E402

rng = random.Random(7)


class CoverFreeWitness(unittest.TestCase):
    def test_witness_point_separates(self):
        # Lemma 13: the witness alpha gives g_c(alpha) different from every other polynomial's value
        Q = li21.smallest_prime_above(50 * 6)
        for _ in range(50):
            c = rng.randrange(2 ** 12)
            others = [rng.randrange(2 ** 12) for _ in range(50)]
            others = [o for o in others if o != c]
            w = li21.witness(c, others, Q)
            alpha, g = divmod(w, Q)
            self.assertEqual(g, li21.poly_value(c, Q, alpha))
            self.assertTrue(all(li21.poly_value(o, Q, alpha) != g for o in others))

    def test_max_degree_matches_enumeration(self):
        n = 8
        self.assertEqual(li21.max_degree(n, 1, 1),
                         max(len(edits.confusion_ball(x, 1, 1)) - 1 for x in edits.all_bits(n)))


class Li21Scheme(unittest.TestCase):
    def test_proper_colouring_on_all_strings(self):
        # every pair of confusable strings gets different colours (n = 8, t = k = 1, all 256 strings)
        s = li21.Li21Scheme(k=1, t=1, n=8)
        colours = {x: s.colour(x) for x in edits.all_bits(8)}
        for x in edits.all_bits(8):
            for z in edits.confusion_ball(x, 1, 1):
                if z != x:
                    self.assertNotEqual(colours[x], colours[z])

    def test_round_trip(self):
        for n in (8, 12):
            s = li21.Li21Scheme(k=1, t=1, n=n)
            for _ in range(6):
                x = edits.random_string(n, rng)
                tag = s.encode(x)
                for _ in range(3):
                    self.assertEqual(s.decode(edits.random_edits(x, 1, 1, rng), tag, n), x)

    def test_tag_length(self):
        s = li21.Li21Scheme(k=1, t=1, n=8)
        self.assertEqual((s.delta, s.Q1, s.Q2, s.tag_len), (64, 521, 1217, 21))


if __name__ == "__main__":
    unittest.main()
