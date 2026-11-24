"""
extension_commits.py -- Extension: do real version-control edits fit the multiple-substring-edit model?
(report, Sec. IV)

Run from the project root:
    python3 experiments/extension_commits.py --repo PATH [--repo PATH2 ...] [--revs REVS]
    python3 experiments/extension_commits.py            (asks for the repository and the revisions)
REVS is any git revision selection (default: --all, every commit on every branch).

Only statistics are written (results/ext_*.csv): file type, sizes, edit counts. No file contents and no file names.

For every non-merge commit and every MODIFIED text file (old version = parent, new version = commit):
  1. changed regions in bits (line diff, then byte diff inside changed blocks)            -> io/diffs.py
  2. t_min(k): fewest k-substring edits explaining the change, for k = 1 bit ... 512 bytes;
     the commit FITS (t, k) if t_min(k) <= t. Pairs that do not fit are model violations, not decoder failures.
  3. the line view used for the documents (Sec. III-C): lines added/removed, and whether each changed line is a
     single edit of at most k bits
  4. density of the old file (Definition 4, paper's delta) for k = 1, 2: the precondition of Theorem 2
  5. estimated tag for the whole file, 2 log2 of the Lemma 3 bound (files are far too long to encode exactly)
  6. real run: 16-bit windows around real single-bit edits, decoded by all three schemes (t = k = 1)
"""
import argparse
import os
import random
import statistics as st
import subprocess

from common import QUICK, write_csv
from docex.core import bounds
from docex.core.edits import apply_edit
from docex.io import diffs
from docex.schemes.average_case.density import partition, delta_paper
from docex.schemes.base import DecodingError
from docex.schemes.worst_case import WorstCaseScheme
from docex.schemes.average_case.scheme import AverageCaseScheme
from docex.schemes.li21 import Li21Scheme

KS = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096]
TS = [1, 2, 3, 5]
MAX_BYTES = 256 * 1024
WINDOW = 16


def git(repo, *args):
    r = subprocess.run(["git", "-C", repo, *args], capture_output=True)
    return r.stdout if r.returncode == 0 else None


def collect(repo, revs):
    """(commit, extension, old bytes, new bytes) for every modified text file in the selected history."""
    commits = (git(repo, "rev-list", "--no-merges", *revs.split()) or b"").decode().split()
    pairs = []
    for c in commits:
        listing = (git(repo, "diff-tree", "--no-commit-id", "--name-status", "-r", "-M0", c) or b"").decode("utf-8", "replace")
        for line in listing.splitlines():
            parts = line.split("\t")
            if len(parts) != 2 or parts[0] != "M":
                continue
            path = parts[1]
            old, new = git(repo, "show", f"{c}^:{path}"), git(repo, "show", f"{c}:{path}")
            if old is None or new is None or old == new or b"\x00" in old + new:
                continue                                   # missing or binary
            if not (1 <= len(old) <= MAX_BYTES and 1 <= len(new) <= MAX_BYTES):
                continue
            ext = os.path.splitext(path)[1].lower() or os.path.basename(path)
            pairs.append((c[:8], ext, old, new))
    return pairs


def analyse(label, pairs):
    rows = []
    for commit, ext, old, new in pairs:
        hs = diffs.bit_hunks(old, new, diffs.byte_hunks(old, new))
        line_pairs, added, removed = diffs.line_changes(old, new)
        n = 8 * len(old)
        row = dict(repo=label, commit=commit, ext=ext, n_bits=n, hunks=len(hs),
                   max_hunk_bits=max((diffs.hunk_size(h) for h in hs), default=0),
                   changed_bits=sum(diffs.hunk_size(h) for h in hs),
                   lines_changed=len(line_pairs), lines_added=added, lines_removed=removed)
        for k in KS:
            row[f"t_k{k}"] = diffs.t_min(hs, k)
        for k in (8, 64):                                    # changed lines that are ONE edit of <= k bits
            row[f"lines_single_edit_k{k}"] = sum(
                diffs.t_min(diffs.bit_hunks(a, b, diffs.byte_hunks(a, b)), k) <= 1 for a, b in line_pairs)
        bits = diffs.bits_of(old)
        for k in (1, 2):
            longest = max(len(u) for u in partition(bits, k)[0])
            row[f"maxpart_k{k}"] = longest
            row[f"dense_k{k}"] = int(longest <= delta_paper(n, k))
        rows.append(row)
    return rows


def real_windows(pairs, rng, want):
    """16-bit windows of the old file around isolated real single-bit edits, with the matching window of the new
    file; returns (window_x, window_y, kind) where kind is flip / ins / del."""
    found = []
    for commit, ext, old, new in pairs:
        hs = diffs.bit_hunks(old, new, diffs.byte_hunks(old, new))
        xb, yb = diffs.bits_of(old), diffs.bits_of(new)
        for (xa, xe, ya, ye) in hs:
            if max(xe - xa, ye - ya) != 1:
                continue
            s = min(max(0, xa - rng.randint(2, WINDOW - 3)), len(xb) - WINDOW)
            if s < 0 or not (s <= xa and xe <= s + WINDOW):
                continue
            xw = xb[s:s + WINDOW]
            yw = apply_edit(xw, xa - s, xe - xa, yb[ya:ye])
            s2 = ya - (xa - s)
            if yb[s2:s2 + len(yw)] != yw:                      # another change inside the window: skip
                continue
            kind = "flip" if xe - xa == ye - ya else ("ins" if ye > ya + (xe - xa) else "del")
            found.append((xw, yw, kind))
    rng.shuffle(found)
    return found[:want]


def run_windows(windows):
    schemes = [WorstCaseScheme(1, 1), AverageCaseScheme(1, 1), Li21Scheme(1, 1, WINDOW)]
    rows = []
    for x, y, kind in windows:
        for s in schemes:
            tag = s.encode(x)
            try:
                ok = s.decode(y, tag, WINDOW) == x
            except DecodingError:
                ok = False
            rows.append(dict(kind=kind, scheme=s.name, bits=tag.bits, recovered=int(ok)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", action="append", help="path of a local git repository (repeatable)")
    ap.add_argument("--revs", default="--all", help="git revision selection (default --all)")
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    repos = args.repo or [input("Path to a local git repository: ").strip()]
    revs = args.revs if args.repo else (input("Revisions to analyse [--all]: ").strip() or "--all")

    all_pairs, rows = [], []
    for repo in repos:
        label = os.path.basename(os.path.abspath(repo))
        pairs = collect(repo, revs)
        print(f"{label}: {len(pairs)} modified text files", flush=True)
        all_pairs += pairs
        rows += analyse(label, pairs)
    write_csv("ext_pairs.csv", rows)
    N = len(rows)

    fit = [dict(k=k, **{f"t_le_{T}": sum(r[f"t_k{k}"] <= T for r in rows) / N for T in TS}) for k in KS]
    write_csv("ext_fit.csv", fit)
    print(f"\n{N} file changes; median file {st.median(r['n_bits'] for r in rows) / 8:.0f} bytes, "
          f"median {st.median(r['hunks'] for r in rows)} changed regions and {st.median(r['changed_bits'] for r in rows)} changed bits")
    print("share of changes that fit t edits of <= k bits")
    print("     k  " + "  ".join(f"t<={T}" for T in TS))
    for f in fit:
        print(f"  {f['k']:5d}  " + "  ".join(f"{f[f't_le_{T}']:5.2f}" for T in TS))

    keep_breaks = sum(r["lines_added"] == 0 and r["lines_removed"] == 0 for r in rows)
    changed = sum(r["lines_changed"] for r in rows)
    print(f"\nline view: {keep_breaks}/{N} changes keep every line break; of {changed} changed lines, "
          f"{sum(r['lines_single_edit_k8'] for r in rows)} are one edit of <= 8 bits and "
          f"{sum(r['lines_single_edit_k64'] for r in rows)} one edit of <= 64 bits")
    for k in (1, 2):
        print(f"density (k={k}, paper's delta): {sum(r[f'dense_k{k}'] for r in rows)}/{N} dense")

    tags = []
    for k in (8, 64, 512):
        for r in rows:
            t = r[f"t_k{k}"]
            if 1 <= t <= 3:
                tags.append(dict(k=k, t=t, n_bits=r["n_bits"], est_tag=bounds.estimated_tag(r["n_bits"], t, k),
                                 four_t_log_n=bounds.theorem1(r["n_bits"], t)))
    write_csv("ext_tag.csv", tags or [dict(k=0, t=0, n_bits=0, est_tag=0, four_t_log_n=0)])
    for k in (8, 64, 512):
        sel = [q for q in tags if q["k"] == k]
        if sel:
            print(f"k={k:4d}: {len(sel)} changes with t<=3; estimated tag median {st.median(q['est_tag'] for q in sel):.0f} bits "
                  f"(4t log n {st.median(q['four_t_log_n'] for q in sel):.0f}, file {st.median(q['n_bits'] for q in sel):.0f} bits)")

    windows = real_windows(all_pairs, random.Random(5), 8 if QUICK else 40)
    print(f"\nreal single-bit edits usable as 16-bit windows: {len(windows)} "
          f"({', '.join(f'{k} {sum(w[2] == k for w in windows)}' for k in ('flip', 'ins', 'del'))})")
    wrows = run_windows(windows)
    if wrows:
        write_csv("ext_windows.csv", wrows)
        for name in dict.fromkeys(r["scheme"] for r in wrows):
            rs = [r for r in wrows if r["scheme"] == name]
            print(f"  {name:26s}: recovered {sum(r['recovered'] for r in rs)}/{len(rs)}, tag {st.mean(r['bits'] for r in rs):.1f} bits")


if __name__ == "__main__":
    main()
