#!/usr/bin/env python3
"""Prints the counts that panels C and D of chry_composite.R quote.

Reads the pileups make_inputs.sh leaves in the work directory:

  <sample>.msy.pile   chrY:13,000,000-20,000,000, position/ref/depth/bases
  <sample>.xtr.pile   chrY:3,050,000-6,234,000, position/depth/non-reference count

  usage:  python3 analysis.py [workdir]        (default: ./work)
"""

import collections
import sys

WORK = sys.argv[1] if len(sys.argv) > 1 else "work"


def consensus(sample, mindepth=1):
    """Majority base per position, for positions covered at least mindepth."""
    out = {}
    for line in open(f"{WORK}/{sample}.msy.pile"):
        pos, ref, depth, bases = line.rstrip("\n").split("\t")
        if int(depth) < mindepth:
            continue
        counts, i = collections.Counter(), 0
        while i < len(bases):
            ch = bases[i]
            if ch in ".,":
                counts[ref.upper()] += 1
                i += 1
            elif ch.upper() in "ACGTN":
                counts[ch.upper()] += 1
                i += 1
            elif ch == "^":          # read start, next char is mapping quality
                i += 2
            elif ch in "+-":         # indel, e.g. +4ACGT
                j, num = i + 1, ""
                while j < len(bases) and bases[j].isdigit():
                    num += bases[j]
                    j += 1
                i = j + int(num or 0)
            else:
                i += 1
        if counts:
            out[int(pos)] = counts.most_common(1)[0][0]
    return out


def sites(sample, mindepth=3, minalt=2, frac=0.8):
    """Positions where the sample's reads disagree with the chrY reference."""
    called, covered = set(), {}
    for line in open(f"{WORK}/{sample}.xtr.pile"):
        pos, depth, alt = (int(x) for x in line.split())
        if depth < mindepth:
            continue
        covered[pos] = (depth, alt)
        if alt >= minalt and alt / depth >= frac:
            called.add(pos)
    return called, covered


WINDOW = 7_000_000
xx = consensus("NA12878_mother")
men = {
    "NA12882 her son": "NA12882_son",
    "NA12889 his grandfather": "NA12889_unrel_male",
    "NA12891 an unrelated Y": "NA12891_her_father",
}
con = {name: consensus(s, mindepth=5) for name, s in men.items()}

print("panel D — coverage of chrY:13-20 Mb")
print(f"  {'NA12878 (XX)':24s}{len(xx):>9,} positions  {100 * len(xx) / WINDOW:5.1f}%")
for name, c in con.items():
    print(f"  {name:24s}{len(c):>9,} positions  {100 * len(c) / WINDOW:5.1f}%")

print(f"\npanel C — scored over the {len(xx):,} positions the XX sample covers")
pairs = [
    ("NA12882 & NA12889 — one Y, two generations apart", con["NA12882 her son"], con["NA12889 his grandfather"]),
    ("NA12882 & NA12891 — two men, different Ys", con["NA12882 her son"], con["NA12891 an unrelated Y"]),
    ("NA12878 & her son NA12882 — no Y in common", xx, con["NA12882 her son"]),
    ("NA12878 & NA12891 — no Y in common", xx, con["NA12891 an unrelated Y"]),
]
for label, a, b in pairs:
    both = [p for p in xx if p in a and p in b]
    diff = sum(1 for p in both if a[p] != b[p])
    print(f"  {label:52s} {100 * (1 - diff / len(both)):8.4f}%  {diff:>5,} differences")

print("\ndirection test — non-reference alleles in the X-transposed region")
her, _ = sites("NA12878_mother")
print(f"  sites called from HER reads ({len(her):,}), carried by:")
for label, s in [("her son NA12882", "NA12882_son"), ("her son NA12886", "NA12886_son2"),
                 ("NA12889", "NA12889_unrel_male"), ("NA12891 her father", "NA12891_her_father"),
                 ("NA12879 her daughter", "NA12879_daughter")]:
    _, cov = sites(s)
    shared = [p for p in her if p in cov]
    carry = [p for p in shared if cov[p][1] / cov[p][0] >= 0.5]
    print(f"    {label:24s} {100 * len(carry) / max(1, len(shared)):5.1f}%  ({len(carry):,}/{len(shared):,})")

his, _ = sites("NA12882_son", mindepth=5, minalt=3)
_, hercov = sites("NA12878_mother")
shared = [p for p in his if p in hercov]
carry = [p for p in shared if hercov[p][1] / hercov[p][0] >= 0.5]
print(f"  sites called from HIS reads ({len(his):,}), carried by:")
print(f"    {'NA12878 his mother':24s} {100 * len(carry) / max(1, len(shared)):5.1f}%  ({len(carry):,}/{len(shared):,})")
