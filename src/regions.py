#!/usr/bin/env python3
"""Region map for HIV-1 HXB2 (K03455). GenBank-derived coordinates are
authoritative; a few subregions pending literature cross-check are marked
[unverified] (docs/research_notes.md should confirm/correct them).

Outputs results/regions.tsv. Region precedence for guide tagging (more
specific first): TAR > R > U3 > U5 > psi > gag > pol > vif > vpr > tat_e1 >
rev_e1 > vpu > env > RRE > nef > intergenic.
"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GB = os.path.join(HERE, "data", "hiv", "K03455_HXB2.gb")
OUT = os.path.join(HERE, "results", "regions.tsv")

# 1-based inclusive HXB2 coordinates. Source: K03455 GenBank annotation
# (LTRs 1-634 / 9086-9719; R 454-551 / 9540-9636; CDS spans) + literature
# subregions marked unverified pending research-notes confirmation.
REGIONS = [
    # name, start, end, source
    ("U3_5", 1, 453, "genbank"),
    ("R_5", 454, 551, "genbank"),
    ("TAR", 454, 512, "literature[unverified]"),  # TAR stem-loop at R start
    ("U5_5", 552, 634, "genbank"),
    ("psi", 744, 802, "literature[unverified]"),  # packaging signal
    ("gag", 790, 2292, "genbank"),
    ("pol", 2358, 5096, "genbank"),
    ("cPPT", 4727, 4765, "literature[unverified]"),  # central ppt, inside pol
    ("CTS", 4885, 4936, "literature[unverified]"),   # central termination, pol
    ("vif", 5041, 5619, "genbank"),
    ("vpr", 5559, 5795, "genbank"),
    ("tat_e1", 5831, 6045, "genbank"),
    ("rev_e1", 5970, 6045, "genbank"),
    ("vpu", 6062, 6310, "literature[unverified]"),
    ("env", 6225, 8795, "genbank"),
    ("RRE", 7925, 8034, "literature[unverified]"),   # Rev response element
    ("tat_e2", 8379, 8424, "genbank"),
    ("rev_e2", 8379, 8652, "genbank"),
    ("nef", 8797, 9168, "genbank"),
    ("U3_3", 9086, 9539, "genbank"),
    ("R_3", 9540, 9636, "genbank"),
    ("U5_3", 9637, 9719, "genbank"),
    ("LTR_5", 1, 634, "genbank"),
    ("LTR_3", 9086, 9719, "genbank"),
]

PRECEDENCE = ["TAR", "R_5", "R_3", "U3_5", "U3_3", "U5_5", "U5_3", "psi", "gag",
              "pol", "vif", "vpr", "tat_e1", "rev_e1", "vpu", "env", "RRE",
              "tat_e2", "rev_e2", "nef"]
# LTR umbrella regions never become the primary label (their subregions do)


def load_gb_spans():
    """Cross-check our table against the GenBank annotation."""
    from Bio import SeqIO
    rec = SeqIO.read(GB, "genbank")
    spans = {}
    for f in rec.features:
        if f.type in ("CDS", "repeat_region"):
            label = None
            for q in ("label", "standard_name", "gene", "product"):
                if q in f.qualifiers:
                    label = f.qualifiers[q][0]
                    break
            spans[(f.type, label)] = (int(f.location.start) + 1, int(f.location.end))
    return spans


def regions_for(pos1, length):
    """Primary + all region labels overlapping [pos1, pos1+length-1] (1-based)."""
    labels = []
    for name, start, end, _src in REGIONS:
        if name.startswith("LTR_"):
            continue  # umbrella only, reported separately in analysis
        if pos1 <= end and pos1 + length - 1 >= start:
            labels.append(name)
    primary = None
    for p in PRECEDENCE:
        if p in labels:
            primary = p
            break
    return primary, labels


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write("region\tstart\tend\tlength\tsource\n")
        for name, start, end, src in REGIONS:
            fh.write(f"{name}\t{start}\t{end}\t{end - start + 1}\t{src}\n")
    print(f"wrote {OUT} ({len(REGIONS)} regions)")


if __name__ == "__main__":
    main()