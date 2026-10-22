# Plan: crispr-hiv-ward in-silico CRISPR guide design vs HIV-1 (2026-10-22)

Goal: identify the best HIV-1 target region (conserved, least similar to human
genes, evidence-backed) and produce off-target-screened SpCas9/SaCas9 guide sets
plus a recommended multiplex excision strategy. Computational only. No GitHub push.

## Stages
1. **Research** (subagent A → docs/research_notes.md, verified citations):
   published targets, EBT-101 design/status, escape, off-target data, conservation,
   HXB2 gene map, design constraints, 2023-26 landscape.
2. **Scanner** (subagent B → src/crispr_scan.c + tests): pack FASTA→.cgb (4-bit),
   exhaustive off-target scan (≤4 mm; Sp: NGG+NAG; Sa: NNGRRT+NNGRRV), exact k-map.
3. **Data** (parent): HXB2/NL4-3 refs; ~24 diverse complete genomes via NCBI
   eutils; GRCh38 primary (Ensembl, fallbacks NCBI/UCSC); GENCODE v44 pc transcripts.
4. **Conservation**: star-align strains to HXB2 (Biopython PairwiseAligner);
   per-position identity; per-guide exact window+PAM match fraction across strains.
5. **Guide enumeration**: all SpCas9 NGG + SaCas9 NNGRRT sites both strands;
   filters (GC 40–65%, no TTTT / no ≥4 homopolymer); region tagging (genbank-derived
   coordinates, precedence TAR>R>U3>U5>psi>gag>pol>vif>vpr>tat_e1>rev_e1>vpu>env>RRE>…);
   LTR guides flagged for presence in both 5' and 3' LTR.
6. **Scans**: pack GRCh38 (+ txome, + HXB2); scan top ~150 Sp + ~40 Sa guides +
   literature benchmark guides (child A may supply sequences); k-map HIV 20/30-mers
   vs human genome + transcriptome.
7. **Analysis**: per-guide off-target counts (≤3 mm = disqualifying), region
   similarity shares, region ranking, top guides, recommended dual/triple-guide
   excision pairs (LTR guide + internal guide, escape mitigation).
8. **Deliverables**: docs/design_report.md + README.md + results/*.tsv; full test
   suite; local git commits (no remote).

## Verification gates
- Scanner: unit tests (hand-derived strand/PAM vectors) + on-target validation:
  every designed guide must hit HXB2 .cgb at its expected coordinate, mm=0.
- Pipeline tests: enumeration/filter/conservation/analysis unit tests.
- Citations: sources.py verify --evidence; parent spot-checks key claims.
- Report: no guide recommended with any ≤3-mm human off-target.