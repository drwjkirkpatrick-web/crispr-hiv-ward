#!/usr/bin/env python3
"""Aggregate scanner outputs into guide scores, a region ranking, and a
recommended multiplex excision strategy.

Inputs (produced by earlier stages):
  results/scan_input.tsv            guides scanned
  results/offtarget_human.tsv       scanner output vs GRCh38
  results/offtarget_hiv.tsv         scanner output vs HXB2 (on-target validation)
  results/guides_sp_all.tsv, guides_sa_all.tsv
  results/kmap_genome_k20.tsv, kmap_genome_k30.tsv, kmap_txome_k20.tsv

Outputs:
  results/guide_scores.tsv, results/region_ranking.tsv,
  results/recommended_strategy.tsv
Run: /usr/bin/python3 src/analyze.py
"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
sys.path.insert(0, os.path.join(HERE, "src"))
from regions import REGIONS  # noqa: E402

GENOME_LEN = 9719


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip() and not line.startswith("#")]


def pam_class(system, pam):
    if system == "SPCAS9":
        return "canonical" if pam[1:] == "GG" else "noncanonical"
    return "canonical" if pam[5] == "T" else "noncanonical"


def main():
    scan_in = read_tsv(os.path.join(RESULTS, "scan_input.tsv"))
    guides = {}
    for g in scan_in:
        guides[g["guide_id"]] = g
    allrows = read_tsv(os.path.join(RESULTS, "guides_sp_all.tsv")) + \
        read_tsv(os.path.join(RESULTS, "guides_sa_all.tsv"))
    by_id = {f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}_{r['pos']}{r['strand']}":
             r for r in allrows}

    # --- on-target validation vs HXB2 ---
    hiv_hits = read_tsv(os.path.join(RESULTS, "offtarget_hiv.tsv"))
    ontarget = {}
    for h in hiv_hits:
        if int(h["mismatches"]) == 0:
            ontarget.setdefault(h["guide_id"], []).append((h["chrom"], h["pos"],
                                                           h["strand"]))
    bad = [gid for gid in guides if gid not in ontarget]
    print(f"on-target validation: {len(guides) - len(bad)}/{len(guides)} guides hit "
          f"HXB2 exactly at the expected coordinate")
    if bad:
        print("  MISSING on-target (excluded):", bad[:10])

    # --- human off-target counts ---
    hits = read_tsv(os.path.join(RESULTS, "offtarget_human.tsv"))
    per = {gid: {"le3": 0, "eq4": 0, "le3_can": 0, "le3_noncan": 0,
                 "n4_can": 0, "n4_noncan": 0} for gid in guides}
    for h in hits:
        gid = h["guide_id"]
        if gid not in per:
            continue
        mm = int(h["mismatches"])
        cls = pam_class(h["system"], h["pam"])
        if mm <= 3:
            per[gid]["le3"] += 1
            per[gid]["le3_can" if cls == "canonical" else "le3_noncan"] += 1
        elif mm == 4:
            per[gid]["eq4"] += 1
            per[gid]["n4_can" if cls == "canonical" else "n4_noncan"] += 1

    # --- region similarity from k-maps ---
    def region_kmer_share(path):
        rows = read_tsv(path)
        hits_per_region = {}
        for r in rows:
            p = int(r["query_pos"])
            for name, start, end, _ in REGIONS:
                if start <= p <= end:
                    hits_per_region.setdefault(name, set()).add(p)
        hits_per_region = {k: len(v) for k, v in hits_per_region.items()}
        total = {}
        for name, start, end, _ in REGIONS:
            total[name] = max(0, end - start + 1 - 19)
        return {k: hits_per_region.get(k, 0) for k in total}, total

    k20_hits, k20_total = region_kmer_share(os.path.join(RESULTS, "kmap_genome_k20.tsv"))
    k30_hits, k30_total = region_kmer_share(os.path.join(RESULTS, "kmap_genome_k30.tsv"))
    tx_hits, _ = region_kmer_share(os.path.join(RESULTS, "kmap_txome_k20.tsv"))

    # --- region ranking ---
    region_stats = {}
    for r in allrows:
        if r["pass"] != "True":
            continue
        name = r["primary"]
        if not name:
            continue
        d = region_stats.setdefault(name, {"n_guides": 0, "cons": [],
                                           "zero_le3": 0, "scanned": 0})
        d["n_guides"] += 1
        d["cons"].append(float(r["cons"]))
        gid = f"{('SP' if r['system'] == 'SPCAS9' else 'SA')}_{r['pos']}{r['strand']}"
        if gid in per:
            d["scanned"] += 1
            if per[gid]["le3"] == 0:
                d["zero_le3"] += 1
    rows = []
    for name, d in region_stats.items():
        k20_share = k20_hits.get(name, 0) / max(1, k20_total.get(name, 1))
        k30_share = k30_hits.get(name, 0) / max(1, k30_total.get(name, 1))
        rows.append({
            "region": name,
            "n_guides_pass": d["n_guides"],
            "median_cons": sorted(d["cons"])[len(d["cons"]) // 2],
            "scanned": d["scanned"],
            "zero_le3_guides": d["zero_le3"],
            "human_k20_share": round(k20_share, 5),
            "human_k30_share": round(k30_share, 5),
            "txome_k20_hits": tx_hits.get(name, 0),
        })
    rows.sort(key=lambda r: (-r["zero_le3_guides"], -r["median_cons"]))
    with open(os.path.join(RESULTS, "region_ranking.tsv"), "w") as fh:
        cols = list(rows[0].keys()) if rows else ["region"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # --- guide scores ---
    scored = []
    for gid, g in guides.items():
        p = per[gid]
        meta = by_id.get(gid, {})
        safe = p["le3"] == 0
        scored.append({
            "guide_id": gid, "system": g["system"], "sequence": g["sequence"],
            "pos": meta.get("pos", ""), "strand": meta.get("strand", ""),
            "region": meta.get("primary", "LITERATURE"),
            "cons": meta.get("cons", ""),
            "gc": meta.get("gc", ""),
            "off_le3": p["le3"], "off_le3_can": p["le3_can"],
            "off_le3_noncan": p["le3_noncan"],
            "off_eq4": p["eq4"],
            "safe_le3": safe,
            "hxb2_ontarget_pos": ";".join(f"{c}:{p0}:{s}" for c, p0, s
                                          in ontarget.get(gid, [])),
        })
    scored.sort(key=lambda r: (r["safe_le3"] is not True, -float(r["cons"] or 0),
                               int(r["off_eq4"])))
    with open(os.path.join(RESULTS, "guide_scores.tsv"), "w") as fh:
        cols = list(scored[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in scored:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # --- recommended multiplex strategy ---
    def best_in(pred, need_safe=True):
        cand = [r for r in scored if pred(r) and (r["safe_le3"] or not need_safe)
                and float(r["cons"] or 0) > 0]
        if not cand:
            return None
        return max(cand, key=lambda r: (float(r["cons"]), -int(r["off_eq4"])))

    safe = [r for r in scored if r["safe_le3"] and float(r["cons"] or 0) > 0.5]
    ltr = best_in(lambda r: r["region"] in ("TAR", "R_5", "R_3", "U3_5", "U3_3",
                                            "U5_5", "U5_3"))
    internal = best_in(lambda r: r["region"] in ("gag", "pol", "psi"))
    internal2 = best_in(lambda r: r["region"] in ("env", "vif", "nef", "RRE"))
    strat = []
    if ltr and internal:
        span = (int(ltr["pos"]), int(internal["pos"]))
        lo, hi = min(span), max(span)
        strat.append({
            "strategy": "dual excision: LTR guide + internal guide",
            "guide_a": ltr["guide_id"], "seq_a": ltr["sequence"],
            "guide_b": internal["guide_id"], "seq_b": internal["sequence"],
            "cons_a": ltr["cons"], "cons_b": internal["cons"],
            "excised_fragment_bp": hi - lo,
            "rationale": "Two independent conserved targets; LTR guide cuts both "
                         "proviral LTRs (they are identical in vivo), excising "
                         "most of the provirus; internal guide mitigates escape."})
    if ltr and internal2 and internal and internal2["guide_id"] != internal["guide_id"]:
        strat.append({
            "strategy": "triple-guide escape-hardened excision",
            "guide_a": ltr["guide_id"], "seq_a": ltr["sequence"],
            "guide_b": internal["guide_id"], "seq_b": internal["sequence"],
            "guide_c": internal2["guide_id"], "seq_c": internal2["sequence"],
            "cons_a": ltr["cons"], "cons_b": internal["cons"],
            "cons_c": internal2["cons"], "excised_fragment_bp": "",
            "rationale": "Multiplexing independent conserved sites reduces "
                         "NHEJ-mediated escape (Wang 2016)."})
    with open(os.path.join(RESULTS, "recommended_strategy.tsv"), "w") as fh:
        if strat:
            cols = list(strat[0].keys())
            fh.write("\t".join(cols) + "\n")
            for s in strat:
                fh.write("\t".join(str(s.get(c, "")) for c in cols) + "\n")
    print("\n=== top safe guides (no <=3mm human off-target) ===")
    for r in [r for r in scored if r["safe_le3"]][:12]:
        print(f"  {r['guide_id']}\t{r['region']}\tcons={r['cons']}\t"
              f"off4={r['off_eq4']}\t{r['sequence']}")
    print(f"\nsafe guides: {sum(1 for r in scored if r['safe_le3'])}/{len(scored)}")
    if strat:
        print(f"strategy file written: {len(strat)} options")
    else:
        print("WARNING: no safe conserved pair found — see guide_scores.tsv")


if __name__ == "__main__":
    main()