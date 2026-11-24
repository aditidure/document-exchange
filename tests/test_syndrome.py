"""
test_syndrome.py -- labeling function and syndrome compression (docex/core/labeling.py, syndrome.py).

Run from the project root:  python3 -m unittest tests.test_syndrome -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits, labeling, syndrome  # noqa: E402

rng = random.Random(1)


class Labeling(unittest.TestCase):
    def test_separates_every_pair(self):
        # property (1) of Lemma 1 holds for all pairs, hence for confusable ones
        labels = [labeling.label_int(x) for x in edits.all_bits(10)]
        self.assertEqual(len(set(labels)), len(labels))

    def test_lemma2_length_example(self):
        # t = 1, k = 2 -> tau = 4 -> 657 log n
        self.assertAlmostEqual(labeling.lemma2_label_bits(1024, 1, 2), 657 * 10)


class SyndromeCompression(unittest.TestCase):
    def test_stage1_toy_example(self):
        # labels 100 (x) vs 103, 97, 110, 112: moduli 2..6 fail, a = 7 works, residue 2
        a = syndrome.find_modulus(100, [103, 97, 110, 112])
        self.assertEqual(a, 7)
        self.assertEqual(100 % a, 2)

    def test_modulus_is_valid_and_smallest(self):
        for n in (8, 10, 12):
            for _ in range(10):
                x = edits.random_string(n, rng)
                nb = [labeling.label_int(z) for z in edits.confusion_ball(x, 1, 1) if z != x]
                fx = labeling.label_int(x)
                a = syndrome.find_modulus(fx, nb)
                self.assertTrue(syndrome.is_valid_modulus(a, fx, nb))
                self.assertFalse(any(syndrome.is_valid_modulus(b, fx, nb) for b in range(2, a)))

    def test_bad_prime_counting(self):
        # proof of Lemma 1: each difference has at most R - 1 distinct prime factors (labels below 2^R)
        R = 12
        for _ in range(200):
            a, b = rng.randrange(2 ** R), rng.randrange(2 ** R)
            if a != b:
                self.assertLessEqual(len(syndrome.distinct_prime_factors(a - b)), R - 1)

    def test_found_modulus_respects_lemma1_bound(self):
        n = 12
        for _ in range(10):
            x = edits.random_string(n, rng)
            nb = [labeling.label_int(z) for z in edits.confusion_ball(x, 1, 1) if z != x]
            a = syndrome.find_modulus(labeling.label_int(x), nb)
            self.assertLessEqual(a, syndrome.lemma1_bound(len(nb), n))

    def test_tag_bits(self):
        self.assertEqual(syndrome.tag_bits(9, 13), 4 + 4)


if __name__ == "__main__":
    unittest.main()
