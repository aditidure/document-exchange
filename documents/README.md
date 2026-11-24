# Hand-edited documents (Implementation, part C)

Each document has two versions:

    documents/docN/original.txt    the version Alice holds (x): 8 lines of 16 bits = 128 bits
    documents/docN/edited.txt      the version Bob holds (y): edit THIS file by hand

Every line is one block of 16 bits, and every scheme is run per block with t = 1, k = 1.

## Rules for editing `edited.txt`

* Change **at most one bit per line**: turn one `0` into `1` or one `1` into `0`.
* Keep every line **exactly 16 characters** of `0`/`1` (no inserting or deleting bits), and keep **8 lines**.
* You may leave lines unchanged. Editing several lines is fine (one bit each).

Why: with fixed-size blocks an insertion or deletion would shift every later block, so the edits must keep the
length; with t = k = 1 each block can absorb one edit of one bit. Edits that break these rules are reported as
violations of the edit model by `experiments/run_documents.py`, not as decoder failures.
