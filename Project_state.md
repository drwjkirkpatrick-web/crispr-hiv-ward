# Project State: crispr-hiv-ward

> **Last updated:** 2026-10-22
> **Current phase:** polish (v1.0 research complete)
> **Overall health:** green

---

## 1. Goal (1–2 sentences)
In-silico design of CRISPR guide RNAs that cleave and inactivate the
integrated HIV-1 proviral genome, identifying the target region that is (a)
conserved across HIV-1 strains, (b) least similar to human genes, and (c)
grounded in published evidence. Computational/literature work only; no
wet-lab; NOT pushed to GitHub yet.

## 2. Current Status
### Done
- [x] Full pipeline: strain fetch (24 diverse genomes), conservation (25
      strains), guide enumeration (1,365 sites -> 668 filtered -> 190 scanned)
- [x] C scanner `crispr_scan` (pack/scan/kmap), 13/13 tests, ~5 min GRCh38
- [x] On-target validation 190/190 exact at expected HXB2 coordinates
- [x] Human off-target scan: 0 exact off-targets; tiers A+/A/B = 63/102/25
- [x] k-mer similarity maps: 0 shared 30-mers, 261 20-mer positions, 9
      transcriptome positions (NOTCH3 trap found in psi leader guides)
- [x] Region ranking (LANL-verified coordinates) + dual/triple strategy
- [x] Paper `paper/hiv-crispr-ward.pdf` (9 pp, 0 overfull boxes, 12/12
      rendered-PDF claim checks, 21 cited sources / 25 rendered quotes)
- [x] README.md; 39/39 tests green; all work committed locally (7 commits)

### In Progress
- [x] GitHub push: drwjkirkpatrick-web/crispr-hiv-ward (public, main, 9 commits)

### Not Started
- [ ] Benchmark published guides (GagD/LTR1) through the tier system
- [ ] Cas13 crRNA design reusing the same machinery

## 3. Architecture & Key Decisions
| Decision | Rationale | Date |
|---|---|---|
| SpCas9 (NGG) + SaCas9 (NNGRRT) design; Cas13 covered in paper only | SaCas9 fits AAV (EBT-101 precedent); SpCas9 is the literature workhorse | 2026-10-22 |
| Custom C scanner, exhaustive <=4 mm, both strands | Cas-OFFinder-equivalent, self-contained on Jetson | 2026-10-22 |
| Tiered safety (A+/A/B), not "zero <=3 mm" | ~30 by-chance <=3 mm sites per 20-mer genome-wide; <=2 mm is the risk band | 2026-10-22 |
| Transcriptome screen of mm<=2 windows | caught psi-leader guides exactly matching human NOTCH3 | 2026-10-22 |
| Regions from GenBank K03455 + LANL HXB2 map | both verified sources; psi literature-only | 2026-10-22 |

## 4. Blockers & Risks
- **Risk:** in-silico tiers are necessary-not-sufficient; bulge off-targets
  not scanned -> CHOPCHOP/CRISPOR re-score before any wet-lab use.
- **Risk:** strain set subtype labels sparse (2/24 records) -> geographic
  diversity is the honesty metric.

## 5. Next Step (only ONE)
> **Next:** Optional: benchmark published guides (GagD/LTR1) through the tier system.

## 6. Environment & Tooling Notes
- Use `/usr/bin/python3` (biopython 1.88 + numpy 2.2.6); the venv `python3`
  (hermes) has no pip/biopython.
- LaTeX: pdflatex + lmodern + microtype[expansion=false] (bitmap-font
  expansion crashes otherwise); require `grep Overfull *.log` = 0.
- Citation ledger: docs/citations/ledger.json (58 sources); render Sources
  blocks only via src/make_sources_block.py — never hand-typed.
- NCBI E-utilities: pace 0.5 s/request, back off on HTTP 429.

## 7. Recent Session Log
- 2026-10-22: Full research cycle completed step-by-step: pipeline, scanner,
  scans, analysis, LANL region verification, 9-page cited paper, README;
  9 commits pushed to github.com/drwjkirkpatrick-web/crispr-hiv-ward.
- 2026-10-22: Two subagents (literature, C scanner) died on provider 429
  limits; parent absorbed both workstreams (scanner was ~90% complete).

## 8. References
- Plan: `.hermes/plans/2026-10-22_crispr-hiv-ward.md`
- Paper: `paper/hiv-crispr-ward.pdf` (+ .tex, generated sources block)
- Citation ledger: `docs/citations/ledger.json`