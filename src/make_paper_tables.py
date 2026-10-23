#!/usr/bin/env python3
"""Generate LaTeX tables + number macros for the paper from results TSVs.
Outputs: paper/numbers.tex, paper/tables/*.tex
Run: /usr/bin/python3 src/make_paper_tables.py
"""
import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
PAPER = os.path.join(HERE, "paper")
sys.path.insert(0, os.path.join(HERE, "src"))
from regions import REGIONS  # noqa: E402

# LANL HXB2-map per-region mean Shannon entropy (2022 DNA filtered web
# alignment), transcribed verbatim from the fetched map export [51].
LANL_ENTROPY = {
    "p17": 0.3664, "p24": 0.2073, "protease": 0.2002, "RT": 0.2022,
    "RNaseH": 0.2214, "IN": 0.1598, "cPPT": 0.0112, "vif": 0.2600,
    "vpr": 0.2616, "tat_e1": 0.3129, "rev_e1": 0.3299, "vpu": 0.4378,
    "gp120": 0.4433, "RRE": 0.2156, "gp41": 0.3234,
}


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip() and not line.startswith("#")]


def esc(s):
    return str(s).replace("_", r"\_")


def main():
    os.makedirs(os.path.join(PAPER, "tables"), exist_ok=True)
    seq = "".join(l.strip() for l in open(os.path.join(HERE, "data", "hiv",
                                                       "K03455_HXB2.fasta"))
                  if not l.startswith(">")).upper()
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def site_pam(gid):
        """(protospacer, PAM) as read 5'->3' on the PAM-bearing strand."""
        pos = int(gid[3:].rstrip("+-"))
        strand = gid[-1]
        sa = gid.startswith("SA")
        L = 21 if sa else 20
        p = pos - 1
        if strand == "+":
            pam_len = 6 if sa else 3
            return seq[p:p + L], seq[p + L:p + L + pam_len]
        pam_len = 6 if sa else 3
        return rc[p:p + L], rc[p - pam_len:p][::-1]

    guides_all = read_tsv(os.path.join(RESULTS, "guides_sp_all.tsv")) + \
        read_tsv(os.path.join(RESULTS, "guides_sa_all.tsv"))
    scores = {r["guide_id"]: r
              for r in read_tsv(os.path.join(RESULTS, "guide_scores.tsv"))}
    scan_in = {r["guide_id"]: r
               for r in read_tsv(os.path.join(RESULTS, "scan_input.tsv"))}

    # median conservation per primary region (filtered guides)
    cons = defaultdict(list)
    for r in guides_all:
        if r["pass"] == "True" and r["primary"]:
            cons[r["primary"]].append(float(r["cons"]))
    med = {k: sorted(v)[len(v) // 2] for k, v in cons.items() if v}

    spans = {name: (start, end, src) for name, start, end, src in REGIONS}

    # --- Table: region map + conservation + human similarity ---
    rr = {r["region"]: r
          for r in read_tsv(os.path.join(RESULTS, "region_ranking.tsv"))}
    order = ["TAR", "R_5", "U3_5", "U5_5", "psi", "p17", "p24", "protease",
             "RT", "RNaseH", "IN", "cPPT", "vif", "vpr", "tat_e1", "rev_e1",
             "vpu", "gp120", "RRE", "gp41", "nef", "U3_3", "R_3", "U5_3"]
    rows = []
    for name in order:
        if name not in spans:
            continue
        start, end, src = spans[name]
        rows.append((name, f"{start}--{end}", src.replace("literature[unverified]", "lit.$^{\\dagger}$"),
                     med.get(name, ""), LANL_ENTROPY.get(name, ""),
                     rr.get(name, {}).get("tierAplus", ""),
                     rr.get(name, {}).get("human_k20_hits", ""),
                     rr.get(name, {}).get("txome_k20_hits", "")))
    with open(os.path.join(PAPER, "tables", "regions.tex"), "w") as fh:
        fh.write("\\begin{tabular}{llllrrrr}\n\\toprule\n"
                 "Region & HXB2 span & Source & Med.\\ cons. & LANL ent. & "
                 "A+ guides & hu\\,20-mer & tx\\,20-mer \\\\\n\\midrule\n")
        for name, span, src, m, ent, ap, hu, tx in rows:
            m = f"{m:.2f}" if m != "" else "--"
            ent = f"{ent:.4f}" if ent != "" else "--"
            fh.write(f"{esc(name)} & {span} & {src} & {m} & {ent} & {ap or 0} & "
                     f"{hu} & {tx} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")

    # --- Table: top tier A+ guides ---
    aplus = [r for r in scores.values() if r["tier"] == "A+"]
    aplus.sort(key=lambda r: (-float(r["cons"]), int(r["off_mm4"])))
    with open(os.path.join(PAPER, "tables", "topguides.tex"), "w") as fh:
        fh.write("\\begin{tabular}{llllrrrr}\n\\toprule\n"
                 "ID & Region & Position & Protospacer ($5'\\to3'$) & PAM & "
                 "System & Cons. & mm$\\le$4 \\\\\n\\midrule\n")
        for r in aplus[:14]:
            gid = r["guide_id"]
            proto, pam = site_pam(gid)
            pos = f"{r['pos']}\\,{r['strand']}"
            fh.write(f"{esc(gid)} & {esc(r['region'])} & {pos} & "
                     f"\\texttt{{{proto}}} & \\texttt{{{pam}}} & "
                     f"{'SaCas9' if gid.startswith('SA') else 'SpCas9'} & "
                     f"{float(r['cons']):.2f} & {r['off_mm4']} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")

    # --- Table: recommended strategies ---
    strat = read_tsv(os.path.join(RESULTS, "recommended_strategy.tsv"))
    with open(os.path.join(PAPER, "tables", "strategy.tex"), "w") as fh:
        fh.write("\\begin{tabular}{p{2.6cm}p{5.6cm}p{2.2cm}p{1.8cm}}\n\\toprule\n"
                 "Strategy & Guides (protospacer/PAM) & Excised span & "
                 "Min.\\ cons. \\\\\n\\midrule\n")
        for s in strat:
            ids = [s.get("guide_a", ""), s.get("guide_b", ""),
                   s.get("guide_c", "")]
            ids = [i for i in ids if i]
            parts = []
            for gid in ids:
                proto, pam = site_pam(gid)
                parts.append(f"\\texttt{{{proto}}}--\\texttt{{{pam}}} "
                             f"({esc(gid)})")
            span = s.get("excised_fragment_bp", "")
            span = (f"{int(span):,}\\,bp" if span
                    else "multi-fragment (three cut sites)")
            cons = min(float(s[k]) for k in ("cons_a", "cons_b", "cons_c")
                       if s.get(k))
            fh.write(f"{s['strategy'].split(':')[0]} & "
                     + " \\newline ".join(parts) + f" & {span} & "
                     f"{cons:.2f} \\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n")

    # --- numbers.tex macros ---
    k20 = read_tsv(os.path.join(RESULTS, "kmap_genome_k20.tsv"))
    k30 = read_tsv(os.path.join(RESULTS, "kmap_genome_k30.tsv"))
    ktx = read_tsv(os.path.join(RESULTS, "kmap_txome_k20.tsv"))
    d20 = len({r["query_pos"] for r in k20})
    d30 = len({r["query_pos"] for r in k30})
    dtx = len({(r["query_name"], r["query_pos"]) for r in ktx})
    n_strains = 24
    tiers = defaultdict(int)
    for r in scores.values():
        tiers[r["tier"]] += 1
    strat0 = strat[0] if strat else {}
    excise = int(strat0.get("excised_fragment_bp", 0))
    macros = {
        "nstrains": n_strains, "nguides": len(scan_in),
        "nsp": sum(1 for r in scan_in.values() if r["system"] == "SPCAS9"),
        "nsa": sum(1 for r in scan_in.values() if r["system"] == "SACAS9"),
        "nexact": 0, "naplus": tiers["A+"], "na": tiers["A"], "nb": tiers["B"],
        "ktwenty": d20, "kthirty": d30, "ktx": dtx,
        "excisebp": f"{excise:,}", "guidesall": len(guides_all),
    }
    with open(os.path.join(PAPER, "numbers.tex"), "w") as fh:
        for k, v in macros.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
    print(f"wrote paper/numbers.tex + 3 tables; tiers={dict(tiers)}, "
          f"k20 distinct pos={d20}, k30={d30}, tx={dtx}")


if __name__ == "__main__":
    main()