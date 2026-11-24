# Hand-edited documents (Implementation, part C)

Put each document in its own folder with two versions:

    documents/doc1/original.txt    the version Alice holds (x)
    documents/doc1/edited.txt      the version Bob holds (y)

Edit `edited.txt` by hand. The documents are split into fixed-size blocks and every scheme is run per block
(see `docex/io/documents.py`, added in a later commit).
