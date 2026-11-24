# Hand-edited documents (Implementation, part C)

Each document has two versions:

    documents/docN/original.txt    the version Alice holds (x): 8 lines of 16 bits
    documents/docN/edited.txt      the version Bob holds (y): edit THIS file by hand

Every line is one block, and every scheme is run per block with t = 1, k = 1.

## Rules for editing `edited.txt`

In each line make **at most one edit of one bit**:

* **flip** a bit (`0` -> `1` or `1` -> `0`), or
* **insert** one bit anywhere in the line (the line becomes 17 bits), or
* **delete** one bit (the line becomes 15 bits).

Keep the **8 lines**: do not add or remove line breaks. You may leave lines unchanged, and edit as many lines as you
like (one edit each).

Why: blocks are delimited by line breaks, so an insertion or deletion only changes its own line and the other
blocks still line up; with t = k = 1 each block can absorb one edit of one bit. A line that breaks the rule (e.g.
two edits) is reported as a violation of the edit model by `experiments/run_documents.py`, not as a decoder failure.
