"""
test_documents.py -- binary documents with one block per line (docex/io/documents.py).

Run from the project root:  python3 -m unittest tests.test_documents -v
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docex.io import documents  # noqa: E402


class Documents(unittest.TestCase):
    def test_read_blocks_one_per_line(self):
        with tempfile.TemporaryDirectory() as d:
            good, bad = os.path.join(d, "a.txt"), os.path.join(d, "b.txt")
            open(good, "w").write("0101\n110\n\n10111\n")
            open(bad, "w").write("0102\n")
            self.assertEqual(documents.read_blocks(good), ["0101", "110", "10111"])
            with self.assertRaises(ValueError):
                documents.read_blocks(bad)

    def test_classify_with_insertions_and_deletions(self):
        x = "1001010011000110"
        ys = [x,                          # unchanged
              "1001010111000110",         # one bit flipped
              "10010100110001100",        # one bit inserted (17 bits)
              "100101001100011",          # one bit deleted (15 bits)
              "1001010111000111"]         # two bits flipped far apart: two edits
        self.assertEqual(documents.classify_blocks([x] * 5, ys, k=1, t=1),
                         ["unchanged", "fits", "fits", "fits", "violation"])

    def test_line_count_must_match(self):
        with self.assertRaises(ValueError):
            documents.classify_blocks(["01", "10"], ["01"], k=1, t=1)

    def test_shipped_documents_are_valid(self):
        root = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "documents")
        for name in ("doc1", "doc2"):
            blocks = documents.read_blocks(os.path.join(root, name, "original.txt"))
            self.assertEqual([len(b) for b in blocks], [16] * 8)


if __name__ == "__main__":
    unittest.main()
