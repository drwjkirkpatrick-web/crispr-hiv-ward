#!/usr/bin/env python3
"""Cross-subtype conservation: star-align every strain to HXB2 and build
per-position base maps, then a per-guide exact window+PAM match fraction.

Outputs:
  results/conservation.npz      (per-strain base-at-HXB2-position arrays, codes
                                 0..3, 4 = gap/ambiguous)
  results/conservation_profile.tsv (pos, hxb2_base, match_frac, major_frac)
Run: /usr/bin/python3 src/conservation.py
"""
import os
import numpy as np
from Bio import Align
from Bio import SeqIO

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HIV = os.path.join(HERE, "data", "hiv")
RESULTS = os.path.join(HERE, "results")
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
COMP = {0: 2, 1: 3, 2: 0, 3: 1}  # complement codes


def load_seq(path):
    return str(next(SeqIO.parse(path, "fasta")).seq).upper()


def code_seq(seq):
    arr = np.full(len(seq), 4, dtype=np.uint8)
    for i, ch in enumerate(seq):
        if ch in CODE:
            arr[i] = CODE[ch]
    return arr


def strain_maps(hxb2, strains):
    """Align each strain to HXB2; return dict name -> base-at-hxb2-position."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    maps = {}
    for name, sseq in strains:
        sarr = code_seq(sseq)
        aln = aligner.align(hxb2, sseq)[0]
        m = np.full(len(hxb2), 4, dtype=np.uint8)
        for (s1a, s1b), (s2a, s2b) in aln.aligned:
            n = s1b - s1a
            m[s1a:s1b] = sarr[s2a:s2b]
            if s2b - s2a != n:
                # paired blocks are equal-length by construction; guard anyway
                n = min(s1b - s1a, s2b - s2a)
                m[s1a:s1a + n] = sarr[s2a:s2a + n]
        maps[name] = m
        pid = float((m == code_seq(hxb2)).mean())
        print(f"  {name}: aligned, {pid:.3%} per-position identity to HXB2")
    return maps


def guide_conservation(maps, hxb2_arr, pos1, length, strand, pam_constrained):
    """Fraction of strains with an exact window+PAM match to the HXB2-derived
    guide. pam_constrained = 1-based HXB2 positions whose base must be
    conserved (the fixed PAM bases; window positions are all constrained)."""
    idx = np.arange(pos1 - 1, pos1 - 1 + length)
    win = hxb2_arr[idx]
    n = 0
    for m in maps.values():
        ok = bool((m[idx] == win).all())
        if ok:
            for p in pam_constrained:
                if m[p - 1] != hxb2_arr[p - 1]:
                    ok = False
                    break
        n += int(ok)
    return n / len(maps)


def main():
    os.makedirs(RESULTS, exist_ok=True)
    hxb2 = load_seq(os.path.join(HIV, "K03455_HXB2.fasta"))
    hxb2_arr = code_seq(hxb2)
    strains = []
    tsv = os.path.join(HIV, "strains.tsv")
    with open(tsv) as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            acc = f[header.index("accession")]
            path = os.path.join(HIV, "strains", acc + ".fasta")
            if os.path.exists(path):
                strains.append((acc, load_seq(path)))
    # always include NL4-3 (subtype B, the common lab strain guides are
    # validated against) — its sequence was extracted from M19921.
    nl = os.path.join(HIV, "NL4-3_proval.fasta")
    if os.path.exists(nl):
        strains.append(("NL4-3", load_seq(nl)))
    print(f"conservation over {len(strains)} strains (HXB2 = reference)")
    maps = strain_maps(hxb2, strains)
    np.savez_compressed(os.path.join(RESULTS, "conservation.npz"),
                        hxb2=hxb2_arr, **{f"strain__{k}": v for k, v in maps.items()})
    match = np.mean([m == hxb2_arr for m in maps.values()], axis=0)
    major = np.zeros(len(hxb2))
    for i in range(len(hxb2)):
        col = [m[i] for m in maps.values()]
        vals = [c for c in col if c != 4]
        if vals:
            major[i] = max(vals.count(v) for v in set(vals)) / len(vals)
    with open(os.path.join(RESULTS, "conservation_profile.tsv"), "w") as fh:
        fh.write("pos\thxb2_base\tmatch_frac\tmajor_frac\n")
        bases = "ACGTN"
        for i in range(len(hxb2)):
            fh.write(f"{i+1}\t{bases[hxb2_arr[i]]}\t{match[i]:.4f}\t{major[i]:.4f}\n")
    cov = float((match >= 0.9).mean())
    print(f"per-position identity >=0.9 over {cov:.1%} of the genome")


if __name__ == "__main__":
    main()