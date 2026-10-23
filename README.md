# crispr-hiv-ward

[![Version](https://img.shields.io/github/v/tag/drwjkirkpatrick-web/crispr-hiv-ward?label=version&sort=semver)](https://github.com/drwjkirkpatrick-web/crispr-hiv-ward/releases)
[![Paper](https://img.shields.io/badge/paper-PDF-red)](paper/hiv-crispr-ward.pdf)
[![Tests](https://img.shields.io/badge/tests-39%2F39-brightgreen)](tests/)
[![Platform](https://img.shields.io/badge/platform-NVIDIA%20Jetson%20ARM64-green)]()
[![Python](https://img.shields.io/badge/python-3.10+-blue)]()
[![License](https://img.shields.io/badge/license-TBD-lightgrey)]()

In-silico design of CRISPR guide RNAs that target and inactivate the
integrated HIV-1 proviral genome — selecting the target regions that are
most conserved across viral strains and least similar to human genes.

**Status: research complete, v1.0.1.**

## Note for researchers: p24 (capsid) — the go-to internal target going forward

Across every metric in this study, the **p24 capsid region of gag (HXB2
1,186–1,878) is the strongest internal target for future guide development**:

- **Largest safe pool of any internal region:** 12 of its 34 guides are
  tier A+ (no human-genomic off-target with ≤2 mismatches *and* no
  transcriptome hit with ≤1 mismatch) — see `results/guide_scores.tsv`.
- **Least human-like coding region:** only 7 exact 20-mer positions shared
  with the entire human genome, and **zero** shared with human
  protein-coding transcripts — compare RT (53 shared) and integrase (40),
  which carry the endogenous-retrovirus footprint.
- **Independently corroborated conservation:** the LANL HXB2 entropy map
  scores p24 at 0.2073 mean Shannon entropy, among the lowest of all coding
  regions (integrase 0.1598, cPPT 0.0112) — i.e., strongly functionally
  constrained.

Top tier-A+ p24 guides (protospacer + PAM, HXB2 coordinates, SpCas9/SaCas9):

| Guide | HXB2 pos | Protospacer + PAM | Cons. (25 strains) |
|---|---|---|---|
| SA_1397+ | 1,397 | `AGACCATCAATGAGGAAGCT`–`GCA` | 0.36 |
| SA_1402+ | 1,402 | `ATCAATGAGGAAGCTGCAGA`–`ATG` | 0.32 |
| SP_1500+ | 1,500 | `AGGAACTACTAGTACCCTTC`–`AGG` | 0.28 |
| SP_1456− | 1,456 | `ATTGCACCAGGCCAGATGAG`–`AGA` | 0.28 |

**Recommended going forward:** pair one 5′LTR-U5 guide (dual-cut excision
geometry) with **two or more p24 guides** rather than a single internal
guide — this maximizes simultaneous-escape requirements at the most
functionally constrained, least human-similar internal region. Before any
wet-lab use: (1) re-score finalists with bulge-aware off-target tools
(CRISPOR/CHOPCHOP/Cas-OFFinder with DNA/RNA bulges), (2) confirm guide
activity across a patient-derived panel spanning subtype B plus global
subtypes (our conservation scoring already covers 25 genomes from 10+
countries), and (3) deep-sequence the capsid target region to quantify
pre-existing variation. Full per-guide data: `results/guide_scores.tsv`,
`results/region_ranking.tsv`, and the paper (§ Results).

## Headline findings

- **Zero exact-match off-targets.** Across all 190 screened guides, not one
  has a single exact (mm=0) off-target site anywhere in the 3.1-billion-base
  human genome (GRCh38, both strands, canonical + non-canonical PAMs,
  scanned to 4 mismatches).
- **HIV is almost perfectly alien to human DNA.** Not a single exact 30-mer
  is shared between HXB2 and the entire human genome; only 261 of 9,700
  possible 20-mer positions (2.7%) share a 20-mer with it. Sharing
  concentrates in *pol* and *env* (RT 53, gp120 42, integrase 40 shared
  20-mer positions) — the endogenous-retrovirus footprint.
- **The conservation trap.** The most conserved guides in the genome (ψ
  leader, 0.92–0.96 across 25 strains) exactly match human *NOTCH3*
  transcripts and were excluded by the transcriptome screen. Conservation
  alone is not a safe design criterion.
- **The sweet spot: 5′ LTR U5.** 0.72–0.88 conservation, one human-shared
  20-mer genome-wide, dual-cut excision geometry (in vivo proviral LTRs are
  identical), five tier-A+ guides. Best internal partner: p24 (12 tier-A+
  guides) or integrase (SA_5067+, 0.48).
- **Recommended design** (all tier A+ — no off-target ≤2 mm, no
  transcriptome hit ≤1 mm):
  - Dual: `TCTAGCAGTGGCGCCCGAAC`/AGG (5′LTR U5, SP_628+) +
    `TGTGTGGCAAGTAGACAGGAT`/GAGGAT (integrase, SA_5067+) — 4,439 bp excision.
  - Triple (escape-hardened): + `GGAGCAGCAGGAAGCACTAT`/GGG (RRE, SP_7794+).

## The paper

`paper/hiv-crispr-ward.pdf` — 9-page research paper with methods, all
tables, an explicit limitations section, and a Sources block with 44
verbatim quotes rendered mechanically from the citation ledger
(`docs/citations/ledger.json`). Rebuild: `cd paper && pdflatex
hiv-crispr-ward.tex` (run twice).

## Pipeline

```bash
bash src/run_pipeline.sh fetch         # NCBI E-utilities strain set
bash src/run_pipeline.sh regions       # HXB2 region table (GenBank + LANL map)
bash src/run_pipeline.sh conservation  # 25-strain star alignment + per-guide scores
bash src/run_pipeline.sh guides        # enumerate + filter SpCas9/SaCas9 guides
bash src/run_pipeline.sh pack          # GRCh38 + transcriptome + HXB2 -> packed binary
bash src/run_pipeline.sh scan_hiv      # on-target validation (190/190 exact)
bash src/run_pipeline.sh scan_human    # genome-wide off-target scan (~5 min)
bash src/run_pipeline.sh kmap          # human-similarity k-mer maps
bash src/run_pipeline.sh analyze       # tiers, region ranking, strategy
bash src/run_pipeline.sh tests         # 39 tests (Python pipeline + C scanner)
```

Requires: `/usr/bin/python3` (biopython, numpy), gcc, ~4 GB free RAM, ~6 GB
disk for the human genome. The C scanner (`src/crispr_scan.c`) packs FASTA
to a 4-bit binary format and exhaustively scans ~3.1 Gbp × 190 guides in
~5 minutes on 6 ARM cores.

## Repository layout

```
paper/            LaTeX paper + generated tables + PDF
src/              pipeline (Python), C scanner, table/source generators
tests/            39 unit tests (26 pipeline + 13 scanner)
results/          guides, off-target tables, region ranking, strategy
docs/             citation ledger + fetched page texts (evidence base)
data/hiv/         HXB2, NL4-3, 24 strain genomes (data/human/ is gitignored)
```

## Safety and scope

Computational/literature research only — no wet-lab work, no viral cultures,
no delivery vectors. Guide designs are *candidates*; any experimental use
requires institutional biosafety review, up-to-date off-target re-screening
(including DNA/RNA bulges), and clinical oversight. This repository is not
medical advice. See the paper's Limitations section.

## Provenance

Every factual claim about prior work traces to a registered source with a
verbatim quote: `docs/citations/ledger.json` (58 sources, fetched pages in
`docs/citations/pages/`). Key anchors: Ebina 2013, Hu 2014, Yin 2017
(saCas9 multiplex excision), Dash 2019, Wang 2016 (escape), J Virol 2019
(the 3.60% off-target figure), the EBT-101 trial record NCT05144386, and
the LANL HXB2 annotated map.
