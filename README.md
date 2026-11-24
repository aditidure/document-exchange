# Codes for Document Exchange under Multiple Substring Edits

EE5350 Project 12. Reference paper: Narayanan, Ramkumar, Bitar, Wachter-Zeh, *Coding schemes for document exchange
under multiple substring edits*, arXiv:2601.18441.

The implementation is built one commit at a time; every commit adds one component with its tests.

## Layout

```
docex/                      library
  core/                     edit model, labeling function, syndrome compression, bounds
  schemes/                  worst_case (Theorem 1), average_case/ (Theorem 2), li21 (higher-complexity scheme)
  io/                       documents: text <-> bits, blocks
experiments/                one script per question; writes CSVs to results/
documents/                  the two hand-edited documents
tests/                      unit tests, one file per component
results/  figures/          generated outputs
report/                     IEEE-style report (LaTeX)
```

## Running

Python 3.9 or newer, standard library only (matplotlib only for figures).

```bash
python3 -m unittest discover -s tests -v      # all tests
```
