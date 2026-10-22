# Project State: crispr-hiv-ward

> **Last updated:** 2026-10-22
> **Current phase:** build
> **Overall health:** green

---

## 1. Goal (1–2 sentences)
In-silico design of CRISPR guide RNAs that cleave and inactivate the integrated
HIV-1 proviral genome, identifying the target region that is (a) conserved across
HIV-1 strains, (b) least similar to human genes, and (c) grounded in published
evidence. Computational/literature work only; no wet-lab, no GitHub push yet.

## 2. Current Status
### Done
- [x] Web research (first pass): Khalili excision strategy, EBT-101 = AAV9-SaCas9-GagD-LTR1,
      escape (Wang 2016), off-target data, Cas13 work
- [x] Scaffold + git init; HXB2 (K03455) FASTA + GenBank fetched; HXB2 = 9719 bp
- [x] HXB2 feature table parsed (LTRs 1–634 / 9086–9719, R 454–551, gag 790–2292,
      pol 2358–5096, vif 5041–5619, vpr 5559–5795, tat 5831–8424, rev 5970–8653,
      env 6225–8795, nef 8797–9168)
- [x] GRCh38 + GENCODE download started (background, proc_e1277f65600b)

### In Progress
- [ ] Subagent A: docs/research_notes.md (literature, verified citations)
- [ ] Subagent B: src/crispr_scan.c + tests (off-target scanner, k-mer map)
- [ ] GRCh38 + transcriptome download (background)
- [ ] Python pipeline (fetch/conservation/guides/analyze) — being written

### Not Started
- [ ] Strain set fetch + conservation analysis
- [ ] Guide enumeration + off-target scans vs human genome
- [ ] Region ranking + design report + README
- [ ] Full test suite green + local commits (NO push)

## 3. Architecture & Key Decisions
| Decision | Rationale | Date |
|---|---|---|
| SpCas9 (NGG) + SaCas9 (NNGRRT) design; Cas13 covered in report only | SaCas9 fits AAV (EBT-101 precedent); SpCas9 is the literature workhorse | 2026-10-22 |
| Off-target scan: custom C tool, ≤4 mismatches, NGG+NAG (Sp) / NNGRRT+NNGRRV (Sa) | Cas-OFFiner-equivalent exhaustive scan, self-contained on Jetson | 2026-10-22 |
| Human similarity: exact 20-mer k-map vs GRCh38 + GENCODE pc transcripts | quantitative "least similar to human genes" evidence | 2026-10-22 |
| Conservation: star-align ~24 diverse subtype genomes to HXB2 | guide must work across group M; escape resistance | 2026-10-22 |
| GRCh38 Ensembl primary assembly (no alts), NCBI/UCSC fallback | standard for off-target analysis | 2026-10-22 |

## 4. Blockers & Risks
- **Risk:** GENCODE/Ensembl URLs stale → fallback chain in download script.
- **Risk:** scanner conventions (window/PAM strand rules) → unit tests + on-target
  validation vs HXB2 itself before trusting human-genome results.
- **Risk:** strain set may skew subtype B → label and report diversity honestly.

## 5. Next Step (only ONE)
> **Next:** Write Python pipeline (fetch_hiv, regions, conservation, guides, analyze) + tests, then run strain fetch and conservation as soon as it lands.

## 6. Environment & Tooling Notes
- Use `/usr/bin/python3` (has biopython 1.88 + numpy 2.2.6); the venv `python3`
  (hermes) has no pip/biopython.
- gcc 11.4 ARM64, 6 cores, ~4 GB RAM available; keep scans < 2 GB RAM.
- Subagents own: A → docs/research_notes.md + docs/citations/; B → src/crispr_scan.c,
  src/Makefile, tests/test_scanner.py. Parent owns all other files.
- No GitHub push (user instruction).

## 7. Recent Session Log
- 2026-10-22: Project created; initial research; scaffold; downloads + 2 subagents
  (literature, C scanner) dispatched in parallel.

## 8. References
- Plan: `.hermes/plans/2026-10-22_crispr-hiv-ward.md`
- Research notes: `docs/research_notes.md` (subagent A)
- Design report: `docs/design_report.md` (final deliverable)