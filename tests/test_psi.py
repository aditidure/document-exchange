"""
test_psi.py -- psi and the one-edit decoder of Lemma 5 (docex/schemes/average_case/psi.py).

Run from the project root:  python3 -m unittest tests.test_psi -v
"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.core import edits  # noqa: E402
from docex.schemes.average_case import hash_h, psi  # noqa: E402

rng = random.Random(5)


class Psi(unittest.TestCase):
    def test_prime_is_prime(self):
        for n in (64, 1024, 2 ** 16):
            L, P = psi.psi_params(n, 1)
            self.assertTrue(psi._is_prime(P))
            self.assertTrue(all(P % d for d in range(2, 1000) if d < P))

    def test_psi_is_short(self):
        # |psi| grows far slower than log n: the same size for n = 2^16 and n = 2^20 (L capped)
        self.assertLessEqual(psi.psi_bits(2 ** 20, 1) - psi.psi_bits(2 ** 16, 1), 2)

    def test_one_edit_recovered(self):
        # Lemma 5 end to end: x recovered from (y, h(x), psi(x)) after one edit (n = 2^13, k = 1)
        n, k, ok, trials = 2 ** 13, 1, 0, 8
        for _ in range(trials):
            x = edits.random_string(n, rng)
            y = edits.random_edit(x, k, rng)
            try:
                ok += psi.decode_one_edit(y, hash_h.h_value(x, k, n), psi.psi_value(x, k), n, k) == x
            except ValueError:
                pass
        self.assertGreaterEqual(ok, trials - 1)

    def test_never_returns_a_wrong_string(self):
        # when the decoder answers, the answer is x (failures are reported, not hidden)
        n, k = 2 ** 12, 1
        for _ in range(10):
            x = edits.random_string(n, rng)
            y = edits.random_edit(x, k, rng)
            try:
                self.assertEqual(psi.decode_one_edit(y, hash_h.h_value(x, k, n), psi.psi_value(x, k), n, k), x)
            except ValueError:
                pass


if __name__ == "__main__":
    unittest.main()
