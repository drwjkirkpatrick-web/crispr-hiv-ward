#!/bin/bash
# Canonical pipeline order for crispr-hiv-ward.
# Use /usr/bin/python3 (biopython + numpy live there); the venv python3 has no pip.
set -e
cd "$(dirname "$0")/.."
PY=/usr/bin/python3

stage="${1:-all}"

run_stage() { [ "$stage" = "all" ] || [ "$stage" = "$1" ]; }

if run_stage fetch; then
  $PY src/fetch_hiv.py
fi
if run_stage regions; then
  $PY src/regions.py
fi
if run_stage conservation; then
  $PY src/conservation.py
fi
if run_stage guides; then
  $PY src/guides.py --top-sp 150
fi

# --- stages below need the C scanner (make -C src) and the human genome ---
if run_stage pack; then
  gunzip -kc data/human/GRCh38.primary.fa.gz > data/human/GRCh38.primary.fa
  ./src/crispr_scan pack data/human/GRCh38.primary.fa data/human/GRCh38.cgb
  ./src/crispr_scan stat data/human/GRCh38.cgb
  gunzip -kc data/human/gencode.v44.pc_transcripts.fa.gz > data/human/pc_txome.fa
  ./src/crispr_scan pack data/human/pc_txome.fa data/human/pc_txome.cgb
  ./src/crispr_scan pack data/hiv/K03455_HXB2.fasta data/hiv/HXB2.cgb
fi
if run_stage scan_hiv; then
  # on-target validation: every designed guide must hit HXB2 exactly
  ./src/crispr_scan scan data/hiv/HXB2.cgb results/scan_input.tsv 4 results/offtarget_hiv.tsv 6
fi
if run_stage scan_human; then
  ./src/crispr_scan scan data/human/GRCh38.cgb results/scan_input.tsv 4 results/offtarget_human.tsv 6
fi
if run_stage kmap; then
  ./src/crispr_scan kmap data/human/GRCh38.cgb data/hiv/K03455_HXB2.fasta 20 results/kmap_genome_k20.tsv 6
  ./src/crispr_scan kmap data/human/GRCh38.cgb data/hiv/K03455_HXB2.fasta 30 results/kmap_genome_k30.tsv 6
  ./src/crispr_scan kmap data/human/pc_txome.cgb data/hiv/K03455_HXB2.fasta 20 results/kmap_txome_k20.tsv 6
fi
if run_stage analyze; then
  $PY src/analyze.py
fi
if run_stage tests; then
  $PY tests/test_pipeline.py
  $PY tests/test_scanner.py
fi
echo "pipeline stage '$stage' complete"