"""
test_documents.py -- binary documents and blocks (docex/io/documents.py).

Run from the project root:  python3 -m unittest tests.test_documents -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.io import documents  # noqa: E402


class Documents(unittest.TestCase):
    def test_read_ignores_line_breaks_and_rejects_other_characters(self):
        with tempfile.TemporaryDirectory() as d:
            good, bad = os.path.join(d, "a.txt"), os.path.join(d, "b.txt")
            open(good, "w").write("0101\n1100\n")
            open(bad, "w").write("0102\n")
            self.assertEqual(documents.read_bits(good), "01011100")
            with self.assertRaises(ValueError):
                documents.read_bits(bad)

    def test_split_blocks(self):
        self.assertEqual(documents.split_blocks("00011011", 2), ["00", "01", "10", "11"])
        with self.assertRaises(ValueError):
            documents.split_blocks("000", 2)

    def test_min_substitution_edits(self):
        self.assertEqual(documents.min_substitution_edits("0000", "0000", 1), 0)
        self.assertEqual(documents.min_substitution_edits("0000", "0100", 1), 1)
        self.assertEqual(documents.min_substitution_edits("0000", "0110", 1), 2)
        self.assertEqual(documents.min_substitution_edits("0000", "0110", 2), 1)
        self.assertEqual(documents.min_substitution_edits("000000", "100001", 2), 2)

    def test_classify(self):
        xs = ["0000", "0000", "0000"]
        ys = ["0000", "0010", "0110"]
        self.assertEqual(documents.classify_blocks(xs, ys, t=1, k=1), ["unchanged", "fits", "violation"])

    def test_shipped_documents_are_valid(self):
        root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "documents")
        for name in ("doc1", "doc2"):
            bits = documents.read_bits(os.path.join(root, name, "original.txt"))
            self.assertEqual(len(documents.split_blocks(bits, 16)), 8)


if __name__ == "__main__":
    unittest.main()
