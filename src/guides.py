#!/usr/bin/env python3
"""Enumerate and filter CRISPR guides across the HXB2 genome.

Conventions (identical to src/crispr_scan.c spec):
  guide sequence = protospacer as read 5'->3' on the PAM-bearing strand,
  immediately 5' of the PAM.
  SpCas9: 20nt guide, PAM NGG. Plus site: plus[p..p+19]=guide, plus[p+21..p+22]=GG.
          Minus site: plus[p..p+19]=revcomp(guide), plus[p-3..p-2]=CC.
  SaCas9: 21nt guide, PAM NNGRRT (canonical). Plus: base(p+23)=G, base(p+24) in AG,
          base(p+25) in AG, base(p+26)=T. Minus: base(p-3)=C, base(p-4) in CT,
          base(p-5) in CT, base(p-6)=A.
  (positions above are 0-based; pos in outputs is 1-based window start)

Filters: GC in [0.40, 0.65]; no run of 4 identical bases (covers U6-terminating
TTTT); guides must not start with TTT.

Outputs:
  results/guides_sp_all.tsv, results/guides_sa_all.tsv  (all enumerated+filtered)
  results/scan_input.tsv  (top Sp + all Sa + optional literature guides)
Run: /usr/bin/python3 src/guides.py [--top-sp 150]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HIV = os.path.join(HERE, "data", "hiv")
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, os.path.join(HERE, "src"))
from regions import regions_for  # noqa: E402

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def load_hxb2():
    with open(os.path.join(HIV, "K03455_HXB2.fasta")) as fh:
        seq = "".join(l.strip() for l in fh if not l.startswith(">")).upper()
    return seq, np.array([CODE.get(c, 4) for c in seq], dtype=np.uint8)


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def pam_constrained_sp(p, L):
    """0-based positions whose base is a fixed PAM base (must be conserved)."""
    return ([p + L + 1, p + L + 2],          # plus: G,G
            [p - 2, p - 3])                   # minus: C,C


def enumerate_sp(seq):
    """Yield (pos1, strand, guide, pam_constrained_0based)."""
    L = 20
    for p in range(0, len(seq) - L - 2):
        if seq[p + L + 1] == "G" and seq[p + L + 2] == "G":  # plus NGG
            yield (p + 1, "+", seq[p:p + L], pam_constrained_sp(p, L)[0])
        if p - 3 >= 0 and seq[p - 3] == "C" and seq[p - 2] == "C":  # minus NGG
            yield (p + 1, "-", revcomp(seq[p:p + L]), pam_constrained_sp(p, L)[1])


def enumerate_sa(seq):
    L = 21
    for p in range(0, len(seq) - L - 5):
        if (seq[p + L + 2] == "G" and seq[p + L + 3] in "AG"
                and seq[p + L + 4] in "AG" and seq[p + L + 5] == "T"):
            yield (p + 1, "+", seq[p:p + L],
                   [p + L + 2, p + L + 3, p + L + 4, p + L + 5])
        if (p - 6 >= 0 and seq[p - 3] == "C" and seq[p - 4] in "CT"
                and seq[p - 5] in "CT" and seq[p - 6] == "A"):
            yield (p + 1, "-", revcomp(seq[p:p + L]),
                   [p - 3, p - 4, p - 5, p - 6])


def pass_filters(guide):
    gc = sum(guide.count(b) for b in "GC") / len(guide)
    if not 0.40 <= gc <= 0.65:
        return False, gc
    for b in "ACGT":
        if b * 4 in guide:
            return False, gc  # homopolymer run >=4 (incl. TTTT U6 terminator)
    if guide.startswith("TTT"):
        return False, gc
    return True, gc


def conservation_frac(maps, hxb2_arr, pos1, length, pam_constrained):
    """Fraction of strains with an exact window+PAM match to the HXB2-derived
    guide. pam_constrained = 0-based positions of fixed PAM bases."""
    idx = np.arange(pos1 - 1, pos1 - 1 + length)
    win = hxb2_arr[idx]
    n = 0
    for m in maps.values():
        ok = bool((m[idx] == win).all())
        if ok and pam_constrained:
            for c in pam_constrained:
                if m[c] != hxb2_arr[c]:
                    ok = False
                    break
        n += int(ok)
    return n / len(maps)


def ltr_presence(seq, guide, system="SPCAS9"):
    """Exact window+PAM occurrences of guide in the full genome, flagged by
    whether they fall in the 5' LTR (<=634) or 3' LTR (>=9086)."""
    hits = {"5": False, "3": False, "sites": 0}
    L = len(guide)
    rc = revcomp(guide)

    def plus_pam(p):
        if system == "SPCAS9":
            return seq[p + L + 1: p + L + 3] == "GG"
        return (seq[p + L + 2] == "G" and seq[p + L + 3] in "AG"
                and seq[p + L + 4] in "AG" and seq[p + L + 5] == "T")

    def minus_pam(p):
        if system == "SPCAS9":
            return p - 3 >= 0 and seq[p - 3: p - 1] == "CC"
        return (p - 6 >= 0 and seq[p - 3] == "C" and seq[p - 4] in "CT"
                and seq[p - 5] in "CT" and seq[p - 6] == "A")

    for p in range(0, len(seq) - L):
        hit = False
        if seq[p:p + L] == guide and plus_pam(p):
            hit = True
        elif seq[p:p + L] == rc and minus_pam(p):
            hit = True
        if hit:
            hits["sites"] += 1
            if p + 1 <= 634:
                hits["5"] = True
            if p + 1 >= 9086:
                hits["3"] = True
    return hits


def build(system, seq, hxb2_arr, maps):
    enum = enumerate_sp if system == "SPCAS9" else enumerate_sa
    rows = []
    for pos1, strand, guide, pam_c in enum(seq):
        ok, gc = pass_filters(guide)
        cons = conservation_frac(maps, hxb2_arr, pos1, len(guide), pam_c)
        primary, labels = regions_for(pos1, len(guide))
        row = {"system": system, "pos": pos1, "strand": strand, "guide": guide,
               "gc": round(gc, 3), "pass": ok, "cons": round(cons, 4),
               "primary": primary or "", "labels": ";".join(labels)}
        if pos1 <= 634 or pos1 >= 9086 - len(guide):
            h = ltr_presence(seq, guide, system)
            row["in_5ltr"], row["in_3ltr"] = h["5"], h["3"]
        else:
            row["in_5ltr"], row["in_3ltr"] = "", ""
        rows.append(row)
    return rows


def write_tsv(path, rows):
    cols = ["guide_id", "system", "sequence", "pos", "strand", "gc", "pass",
            "cons", "primary", "labels", "in_5ltr", "in_3ltr"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r = dict(r, guide_id=f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}"
                               f"_{r['pos']}{r['strand']}")
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main():
    top_sp = 150
    if "--top-sp" in sys.argv:
        top_sp = int(sys.argv[sys.argv.index("--top-sp") + 1])
    os.makedirs(RESULTS, exist_ok=True)
    seq, hxb2_arr = load_hxb2()
    npz = np.load(os.path.join(RESULTS, "conservation.npz"))
    maps = {k[len("strain__"):]: npz[k] for k in npz.files if k.startswith("strain__")}
    print(f"loaded {len(maps)} strain maps")

    sp = build("SPCAS9", seq, hxb2_arr, maps)
    sa = build("SACAS9", seq, hxb2_arr, maps)
    write_tsv(os.path.join(RESULTS, "guides_sp_all.tsv"), sp)
    write_tsv(os.path.join(RESULTS, "guides_sa_all.tsv"), sa)
    print(f"SPCAS9: {len(sp)} sites ({sum(r['pass'] for r in sp)} pass filters)")
    print(f"SACAS9: {len(sa)} sites ({sum(r['pass'] for r in sa)} pass filters)")

    good_sp = sorted([r for r in sp if r["pass"]],
                     key=lambda r: (-r["cons"], abs(r["gc"] - 0.5), r["pos"]))
    good_sa = sorted([r for r in sa if r["pass"]],
                     key=lambda r: (-r["cons"], abs(r["gc"] - 0.5), r["pos"]))
    sel = good_sp[:top_sp] + good_sa[:40]

    # Optional literature benchmark guides for comparison
    lit_path = os.path.join(HERE, "docs", "literature_guides.tsv")
    lit = []
    if os.path.exists(lit_path):
        with open(lit_path) as fh:
            fh.readline()
            for line in fh:
                f = line.strip().split("\t")
                if len(f) >= 3:
                    lit.append({"guide_id": f[0], "system": f[1].upper(),
                                "sequence": f[2].upper()})
        print(f"loaded {len(lit)} literature benchmark guides")

    out = os.path.join(RESULTS, "scan_input.tsv")
    with open(out, "w") as fh:
        fh.write("#guide_id\tsystem\tsequence\n")
        for r in sel:
            fh.write(f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}"
                     f"_{r['pos']}{r['strand']}\t{r['system']}\t{r['guide']}\n")
        for r in lit:
            fh.write(f"{r['guide_id']}\t{r['system']}\t{r['sequence']}\n")
    print(f"wrote {out} ({len(sel) + len(lit)} guides to scan)")


if __name__ == "__main__":
    main()