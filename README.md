# crispr-hiv-ward

In-silico design of CRISPR guide RNAs that target and inactivate the
integrated HIV-1 proviral genome — selecting the target regions that are
most conserved across viral strains and least similar to human genes.

**Status: research complete, v1.0 — local repository only (not on GitHub).**

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
