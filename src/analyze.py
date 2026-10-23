#!/usr/bin/env python3
"""Aggregate scanner outputs into tiered guide scores, a region ranking, and a
recommended multiplex excision strategy.

Safety tiers (empirically calibrated: a random 20-mer has ~30 by-chance <=3mm
sites genome-wide, so 'zero <=3mm' is not a meaningful bar):
  A+ : no genomic off-target with <=2 mismatches AND no transcriptome hit <=1 mm
  A  : no genomic off-target with <=1 mismatch AND no exact transcriptome hit
  B  : everything else (has a <=1 mm genomic off-target)
Cleavage at >=3 mismatches is weak; mm<=2 sites are the risk band.

Inputs:
  results/scan_input.tsv, offtarget_human.tsv, offtarget_hiv.tsv,
  offtarget_lowmm_txome.tsv, guides_sp_all.tsv, guides_sa_all.tsv,
  kmap_genome_k20.tsv, kmap_genome_k30.tsv, kmap_txome_k20.tsv
Outputs:
  results/guide_scores.tsv, results/region_ranking.tsv,
  results/txome_hits_summary.tsv, results/recommended_strategy.tsv
Run: /usr/bin/python3 src/analyze.py
"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, os.path.join(HERE, "src"))
from regions import REGIONS  # noqa: E402


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip() and not line.startswith("#")]


def gene_name(chrom_field):
    """GENCODE transcript FASTA header: ENST|ENSG|OTTHUMT|OTTHUMG|GENE|..."""
    for f in chrom_field.split("|"):
        if not f or f == "-" or f.startswith(("ENST", "ENSG", "OTTH")):
            continue
        return f
    return chrom_field[:30]


def pam_class(system, pam):
    if system == "SPCAS9":
        return "canonical" if pam[1:] == "GG" else "noncanonical"
    return "canonical" if pam[5] == "T" else "noncanonical"


def main():
    scan_in = read_tsv(os.path.join(RESULTS, "scan_input.tsv"))
    guides = {g["guide_id"]: g for g in scan_in}
    allrows = read_tsv(os.path.join(RESULTS, "guides_sp_all.tsv")) + \
        read_tsv(os.path.join(RESULTS, "guides_sa_all.tsv"))
    by_id = {f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}_"
             f"{r['pos']}{r['strand']}": r for r in allrows}

    # on-target validation vs HXB2
    ontarget = {}
    for h in read_tsv(os.path.join(RESULTS, "offtarget_hiv.tsv")):
        if h["mismatches"] == "0":
            ontarget.setdefault(h["guide_id"], []).append(
                (h["chrom"], h["pos"], h["strand"]))
    bad = [gid for gid in guides if gid not in ontarget]
    print(f"on-target validation: {len(guides) - len(bad)}/{len(guides)} "
          f"hit HXB2 exactly at the expected coordinate")
    if bad:
        print("  no exact HXB2 hit (literature benchmarks, expected):", bad[:10])

    # human genomic off-targets
    hits = read_tsv(os.path.join(RESULTS, "offtarget_human.tsv"))
    per = {gid: {"mm0": 0, "mm1": 0, "mm2": 0, "mm3": 0, "mm4": 0,
                 "mm1_can": 0, "mm2_can": 0} for gid in guides}
    for h in hits:
        gid = h["guide_id"]
        if gid not in per:
            continue
        mm = int(h["mismatches"])
        per[gid][f"mm{mm}"] += 1
        if mm in (1, 2) and pam_class(h["system"], h["pam"]) == "canonical":
            per[gid][f"mm{mm}_can"] += 1
    n_exact = sum(p["mm0"] for p in per.values())
    print(f"exact (mm=0) human off-targets across all {len(guides)} guides: {n_exact}")

    # transcriptome hits for the mm<=2 genomic windows
    tx = read_tsv(os.path.join(RESULTS, "offtarget_lowmm_txome.tsv"))
    tx_per = {gid: {"tx0": 0, "tx1": 0, "genes": []} for gid in guides}
    tx_rows_out = []
    for h in tx:
        base = h["guide_id"].split("_", 1)[1] if "_" in h["guide_id"] else ""
        if base not in tx_per:
            continue
        mm = int(h["mismatches"])
        gene = gene_name(h["chrom"])
        tx_per[base][f"tx{mm}"] += 1
        if gene:
            tx_per[base]["genes"].append((gene, mm, h["chrom"], h["pos"]))
        if mm <= 1:
            tx_rows_out.append((base, gene, mm, h["pos"], h["strand"],
                                h["pam"], h["plus_window"]))
    with open(os.path.join(RESULTS, "txome_hits_summary.tsv"), "w") as fh:
        fh.write("guide_id\tgene\tmm\ttx_pos\ttx_strand\tpam\twindow\n")
        for r in sorted(tx_rows_out):
            fh.write("\t".join(str(x) for x in r) + "\n")

    # guide scores + tiers
    scored = []
    for gid, g in guides.items():
        p = per[gid]
        t = tx_per.get(gid, {"tx0": 0, "tx1": 0, "genes": []})
        meta = by_id.get(gid, {})
        if p["mm0"] == 0 and p["mm1"] == 0 and p["mm2"] == 0 and t["tx0"] == 0 \
                and t["tx1"] == 0:
            tier = "A+"
        elif p["mm0"] == 0 and p["mm1"] == 0 and t["tx0"] == 0:
            tier = "A"
        else:
            tier = "B"
        genes = ";".join(sorted({g0 for g0, mm, *_ in t["genes"] if mm <= 1}))
        scored.append({
            "guide_id": gid, "system": g["system"], "sequence": g["sequence"],
            "pos": meta.get("pos", ""), "strand": meta.get("strand", ""),
            "region": meta.get("primary", "LITERATURE"),
            "cons": meta.get("cons", ""), "gc": meta.get("gc", ""),
            "tier": tier,
            "off_mm0": p["mm0"], "off_mm1": p["mm1"], "off_mm2": p["mm2"],
            "off_mm3": p["mm3"], "off_mm4": p["mm4"],
            "off_mm1_can": p["mm1_can"], "off_mm2_can": p["mm2_can"],
            "tx_mm0": t["tx0"], "tx_mm1": t["tx1"], "tx_genes": genes,
            "hxb2_ontarget": ";".join(f"{p0}{s}" for _c, p0, s
                                      in ontarget.get(gid, [])),
        })
    trank = {"A+": 0, "A": 1, "B": 2}
    scored.sort(key=lambda r: (trank[r["tier"]],
                               -float(r["cons"] or 0), int(r["off_mm2"])))
    with open(os.path.join(RESULTS, "guide_scores.tsv"), "w") as fh:
        cols = list(scored[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in scored:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    tiers = {}
    for r in scored:
        tiers[r["tier"]] = tiers.get(r["tier"], 0) + 1
    print(f"tier counts: {tiers}")

    # region ranking
    region_stats = {}
    for r in allrows:
        if r["pass"] != "True":
            continue
        name = r["primary"]
        if not name:
            continue
        d = region_stats.setdefault(
            name, {"n": 0, "cons": [], "scanned": 0, "aplus": 0, "a": 0})
        d["n"] += 1
        d["cons"].append(float(r["cons"]))
        gid = f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}_{r['pos']}{r['strand']}"
        if gid in {s["guide_id"] for s in scored}:
            row = next(s for s in scored if s["guide_id"] == gid)
            d["scanned"] += 1
            d["aplus"] += int(row["tier"] == "A+")
            d["a"] += int(row["tier"] == "A")

    def region_kmer_share(path):
        hits_per, total = {}, {}
        for r in read_tsv(path):
            p = int(r["query_pos"])
            for name, start, end, _ in REGIONS:
                if start <= p <= end:
                    hits_per.setdefault(name, set()).add(p)
        for name, start, end, _ in REGIONS:
            total[name] = max(0, end - start + 1 - 19)
        return {k: len(hits_per.get(k, ())) for k in total}, total

    k20, k20t = region_kmer_share(os.path.join(RESULTS, "kmap_genome_k20.tsv"))
    k30, k30t = region_kmer_share(os.path.join(RESULTS, "kmap_genome_k30.tsv"))
    txx, _ = region_kmer_share(os.path.join(RESULTS, "kmap_txome_k20.tsv"))

    rows = []
    for name, d in region_stats.items():
        rows.append({
            "region": name, "n_guides_pass": d["n"],
            "median_cons": round(sorted(d["cons"])[len(d["cons"]) // 2], 3),
            "scanned": d["scanned"],
            "tierAplus": d["aplus"], "tierA_or_better": d["a"],
            "human_k20_hits": k20.get(name, 0),
            "human_k30_hits": k30.get(name, 0),
            "txome_k20_hits": txx.get(name, 0),
        })
    rows.sort(key=lambda r: (-r["tierAplus"], -r["tierA_or_better"],
                             -r["median_cons"]))
    with open(os.path.join(RESULTS, "region_ranking.tsv"), "w") as fh:
        cols = list(rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # recommended multiplex strategy (tier A+ guides, fallback tier A)
    def best_in(pred, tiers=("A+",)):
        cand = [r for r in scored if pred(r["region"]) and r["tier"] in tiers
                and float(r["cons"] or 0) > 0]
        return max(cand, key=lambda r: (float(r["cons"]), -int(r["off_mm4"]))) \
            if cand else None

    def best_any(pred):
        return best_in(pred, ("A+", "A"))

    ltr = best_in(lambda r: r in ("TAR", "R_5", "R_3", "U3_5", "U3_3",
                                  "U5_5", "U5_3")) or \
        best_any(lambda r: r in ("TAR", "R_5", "R_3", "U3_5", "U3_3",
                                 "U5_5", "U5_3"))

    def internal_far(regions):
        """Best internal guide >=2000 bp from the LTR cut (real excision span)."""
        lpos = int(ltr["pos"]) if ltr else 0
        for tiers in (("A+",), ("A+", "A")):
            cand = [r for r in scored if r["region"] in regions
                    and r["tier"] in tiers and float(r["cons"] or 0) > 0
                    and (not lpos or abs(int(r["pos"]) - lpos) >= 2000)]
            if cand:
                return max(cand, key=lambda r: (float(r["cons"]),
                                                -int(r["off_mm4"])))
        return None

    internal = internal_far(("gag", "p17", "p24", "pol", "RT", "RNaseH", "IN",
                             "protease", "psi"))
    internal2 = internal_far(("env", "gp120", "gp41", "vif", "nef", "RRE"))
    strat = []
    if ltr and internal:
        lo, hi = sorted((int(ltr["pos"]), int(internal["pos"])))
        strat.append({
            "strategy": "dual excision: LTR guide + internal guide",
            "guide_a": ltr["guide_id"], "seq_a": ltr["sequence"],
            "guide_b": internal["guide_id"], "seq_b": internal["sequence"],
            "tier_a": ltr["tier"], "tier_b": internal["tier"],
            "cons_a": ltr["cons"], "cons_b": internal["cons"],
            "excised_fragment_bp": hi - lo,
            "rationale": "LTR guide cuts both proviral LTRs (identical in vivo) "
                         "and the internal guide is an independent conserved "
                         "target that mitigates NHEJ-mediated escape."})
    if ltr and internal and internal2 and internal2["guide_id"] != internal["guide_id"]:
        strat.append({
            "strategy": "triple-guide escape-hardened excision",
            "guide_a": ltr["guide_id"], "seq_a": ltr["sequence"],
            "guide_b": internal["guide_id"], "seq_b": internal["sequence"],
            "guide_c": internal2["guide_id"], "seq_c": internal2["sequence"],
            "tier_a": ltr["tier"], "tier_b": internal["tier"],
            "tier_c": internal2["tier"],
            "cons_a": ltr["cons"], "cons_b": internal["cons"],
            "cons_c": internal2["cons"], "excised_fragment_bp": "",
            "rationale": "Multiplexing independent conserved sites reduces "
                         "viral escape (Wang 2016)."})
    with open(os.path.join(RESULTS, "recommended_strategy.tsv"), "w") as fh:
        if strat:
            cols = ["strategy", "guide_a", "seq_a", "guide_b", "seq_b",
                    "guide_c", "seq_c", "tier_a", "tier_b", "tier_c",
                    "cons_a", "cons_b", "cons_c", "excised_fragment_bp",
                    "rationale"]
            fh.write("\t".join(cols) + "\n")
            for s in strat:
                fh.write("\t".join(str(s.get(c, "")) for c in cols) + "\n")

    # console summary
    print("\n=== region ranking (top 12) ===")
    for r in rows[:12]:
        print(f"  {r['region']:<8} cons={r['median_cons']:<6} A+={r['tierAplus']:<3}"
              f"A={r['tierA_or_better']:<3} hu20={r['human_k20_hits']:<4}"
              f"hu30={r['human_k30_hits']} tx20={r['txome_k20_hits']}")
    print("\n=== top tier A+/A guides ===")
    for r in [r for r in scored if r["tier"] in ("A+", "A")][:15]:
        print(f"  {r['guide_id']:<10} {r['tier']:<3} {r['region']:<8} "
              f"cons={r['cons']:<6} mm2={r['off_mm2']:<4} {r['sequence']}")
    print(f"\nstrategy options: {len(strat)}")
    for s in strat:
        print(f"  {s['strategy']}: {s['guide_a']} + {s['guide_b']}"
              + (f" + {s['guide_c']}" if 'guide_c' in s else ""))


if __name__ == "__main__":
    main()