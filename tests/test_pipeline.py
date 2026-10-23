#!/usr/bin/env python3
"""Unit tests for the crispr-hiv-ward Python pipeline.
Run: /usr/bin/python3 tests/test_pipeline.py
Plain asserts, no pytest dependency. Exits non-zero on failure.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))

from regions import regions_for  # noqa: E402
import guides as G  # noqa: E402

PASS, FAIL = 0, 0


def check(name, cond):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"PASS {name}")
    else:
        FAIL += 1
        print(f"FAIL {name}")


def test_regions():
    primary, labels = regions_for(454, 20)
    check("regions: TAR precedence at 454", primary == "TAR" and "R_5" in labels)
    primary, _ = regions_for(700, 20)
    check("regions: intergenic at 700", primary is None)
    primary, labels = regions_for(1000, 20)
    check("regions: p17 (LANL subregion) at 1000", primary == "p17" and "gag" in labels)
    primary, _ = regions_for(540, 20)
    check("regions: R_5 beats U5_5 at 540", primary == "R_5")
    primary, _ = regions_for(556, 20)
    check("regions: U5_5 at 556", primary == "U5_5")
    primary, _ = regions_for(9086, 20)
    check("regions: 3' LTR U3 at 9086", primary == "U3_3")
    primary, labels = regions_for(9086 + 400, 20)
    check("regions: U3_3 mid", primary == "U3_3")


def test_enumeration_sp():
    seq = "ACGTACGTACGTACGTACGTAGG"  # plus site: guide + AGG
    sites = list(G.enumerate_sp(seq))
    plus = [s for s in sites if s[2] == "ACGTACGTACGTACGTACGT"]
    check("enum sp: plus site found at pos 1, strand +",
          len(plus) == 1 and plus[0][0] == 1 and plus[0][1] == "+")
    check("enum sp: plus PAM constrained = 0-based [21, 22]",
          sorted(plus[0][3]) == [21, 22])
    g = "AAACCCGGGTTTACGTACGT"
    rc = G.revcomp(g)
    check("enum sp: revcomp worked example", rc == "ACGTACGTAAACCCGGGTTT")
    seq2 = "TTAC" + "CC" + "A" + rc + "GGGGGGGGGG"
    sites2 = list(G.enumerate_sp(seq2))
    minus = [s for s in sites2 if s[2] == g]
    # window starts right after 'TTACCCA' -> 0-based p=7 -> 1-based pos 8
    check("enum sp: minus site found (worked example)",
          len(minus) == 1 and minus[0][1] == "-" and minus[0][0] == 8)
    check("enum sp: minus PAM constrained = 0-based [4, 5]",
          sorted(minus[0][3]) == [4, 5])
    seq3 = "TT" + rc + "CC" + "AAAAAAAAA"
    sites3 = [s for s in G.enumerate_sp(seq3) if s[2] == g]
    check("enum sp: downstream CC is not a minus site", len(sites3) == 0)


def test_enumeration_sa():
    g = "ACGTACGTACGTACGTACGTA"
    seq = g + "AAGAAT" + "TTTTTTT"
    sites = list(G.enumerate_sa(seq))
    hits = [s for s in sites if s[2] == g]
    check("enum sa: plus NNGAAT site (canonical NNGRRT)",
          len(hits) == 1 and hits[0][1] == "+")
    seq2 = g + "AAGAAC"  # last PAM base C: design enumeration accepts only NNGRRT
    sites2 = [s for s in G.enumerate_sa(seq2) if s[2] == g]
    check("enum sa: NNGAAC not a design site", len(sites2) == 0)


def test_filters():
    ok, gc = G.pass_filters("ACGTACGTACGTACGTACGT")
    check("filter: balanced guide passes", ok and abs(gc - 0.5) < 0.01)
    ok, _ = G.pass_filters("AAATTTTACGCGCGCGCGC")
    check("filter: TTTT rejected", not ok)
    ok, _ = G.pass_filters("AAAGGGGACACACACACA")
    check("filter: GGGG rejected", not ok)
    ok2, gc2 = G.pass_filters("ATATATATATATATATAT")
    check("filter: low-GC rejected", not ok2 and gc2 < 0.4)
    ok3, _ = G.pass_filters("GGGGTTAAGGTTAACCGGTT")
    check("filter: GGGG rejected (high GC case)", not ok3)


def test_conservation():
    hxb2_arr = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                         0, 1, 2, 3, 0, 0, 0, 0], dtype=np.uint8)
    maps = {"s1": hxb2_arr.copy(), "s2": hxb2_arr.copy(),
            "s3": hxb2_arr.copy()}
    maps["s3"][0] = 2  # one mismatch at window pos 1
    frac = G.conservation_frac(maps, hxb2_arr, 1, 20, [21, 22])
    check("conservation: 2/3 strains exact", abs(frac - 2 / 3) < 1e-9)
    maps["s2"][21] = 1  # PAM G->C in strain 2 (0-based 21 = pos 22)
    frac = G.conservation_frac(maps, hxb2_arr, 1, 20, [21, 22])
    check("conservation: PAM mutation disqualifies", frac == 1 / 3)


def test_ltr_presence():
    g = "ACGTACGTACGTACGTACGT"
    seq = "N" * 620 + g + "AGG" + "N" * 8442 + g + "AGG" + "N" * 700
    h = G.ltr_presence(seq, g)
    check("ltr presence: both 5' and 3' copies detected",
          h["5"] and h["3"] and h["sites"] == 2)
    seq2 = "N" * 620 + g + "AGG" + "N" * 9140
    h2 = G.ltr_presence(seq2, g)
    check("ltr presence: only 5' copy", h2["5"] and not h2["3"] and h2["sites"] == 1)
    # SaCas9 guide: PAM must be NNGRRT, not NGG
    gsa = "ACGTACGTACGTACGTACGTA"
    seqsa = "N" * 620 + gsa + "AAGAAT" + "N" * 9140
    hsa = G.ltr_presence(seqsa, gsa, system="SACAS9")
    check("ltr presence: Sa PAM accepted for Sa guide",
          hsa["5"] and hsa["sites"] == 1)
    hsa_bad = G.ltr_presence("N" * 620 + gsa + "AAGAAG" + "N" * 9140, gsa,
                             system="SACAS9")
    check("ltr presence: Sa wrong PAM rejected", hsa_bad["sites"] == 0)


def main():
    test_regions()
    test_enumeration_sp()
    test_enumeration_sa()
    test_filters()
    test_conservation()
    test_ltr_presence()
    print(f"\n{PASS} passed, {FAIL} failed")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()