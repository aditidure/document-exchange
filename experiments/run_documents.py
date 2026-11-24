"""
run_documents.py -- the three schemes on the two hand-edited documents (report, Sec. III-C).

Run from the project root:   python3 experiments/run_documents.py
Writes results/documents_blocks.csv (one row per document, block and scheme) and results/documents_summary.csv.

For every document in documents/ (original.txt = Alice's x, edited.txt = Bob's y):
  1. read both versions, one block per line (B bits per original block);
  2. classify each block against the edit model: unchanged / fits (edited line in B_t(original line)) / violation;
  3. for every scheme and every block: encode the original block, decode from the edited block, and record the tag
     length, the encoding and decoding time and whether the original block was recovered.
Blocks that violate the edit model are still run, but are reported separately: a failure there is a departure from
the model, not a decoder failure.
"""
import os

from common import ROOT, Timer, write_csv
from docex.io.documents import read_blocks, classify_blocks
from docex.schemes.base import DecodingError
from docex.schemes.worst_case import WorstCaseScheme
from docex.schemes.average_case.scheme import AverageCaseScheme
from docex.schemes.li21 import Li21Scheme

B, T, K = 16, 1, 1
DOCS = os.path.join(ROOT, "documents")


def main():
    schemes = [WorstCaseScheme(K, T), AverageCaseScheme(K, T), Li21Scheme(K, T, B)]
    rows, summary = [], []
    for name in sorted(d for d in os.listdir(DOCS) if os.path.isdir(os.path.join(DOCS, d))):
        xs = read_blocks(os.path.join(DOCS, name, "original.txt"))
        ys = read_blocks(os.path.join(DOCS, name, "edited.txt"))
        status = classify_blocks(xs, ys, K, T)
        print(f"{name}: {len(xs)} blocks of {B} bits | unchanged {status.count('unchanged')}, fits {status.count('fits')},"
              f" violations {status.count('violation')}")
        for s in schemes:
            total_bits = enc_total = dec_total = 0.0
            recovered = 0
            for j, (x, y, st) in enumerate(zip(xs, ys, status)):
                with Timer() as te:
                    tag = s.encode(x)
                try:
                    with Timer() as td:
                        ok = s.decode(y, tag, B) == x
                    dec_s = td.s
                except DecodingError:
                    ok, dec_s = False, float("nan")
                rows.append(dict(document=name, block=j, status=st, scheme=s.name, bits=tag.bits, enc_s=te.s, dec_s=dec_s,
                                 recovered=int(ok), original=x, edited=y))
                total_bits += tag.bits
                enc_total += te.s
                dec_total += dec_s if dec_s == dec_s else 0.0
                recovered += ok
            in_model = [st != "violation" for st in status]
            summary.append(dict(document=name, scheme=s.name, blocks=len(xs), document_bits=len(xs) * B,
                                total_tag_bits=total_bits, recovered=recovered,
                                recovered_in_model=sum(r for r, m in zip([row["recovered"] for row in rows[-len(xs):]], in_model) if m),
                                blocks_in_model=sum(in_model), enc_total_s=enc_total, dec_total_s=dec_total))
            r = summary[-1]
            print(f"  {s.name:26s}: total tag {r['total_tag_bits']:5.0f} bits for {r['document_bits']} document bits | "
                  f"recovered {r['recovered']}/{r['blocks']} (in-model {r['recovered_in_model']}/{r['blocks_in_model']}) | "
                  f"enc {enc_total:.2f}s dec {dec_total:.2f}s", flush=True)
    write_csv("documents_blocks.csv", rows)
    write_csv("documents_summary.csv", summary)


if __name__ == "__main__":
    main()
