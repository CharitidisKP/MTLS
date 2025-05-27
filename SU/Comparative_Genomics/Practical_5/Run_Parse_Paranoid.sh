#!/usr/bin/env bash
set -euo pipefail

# Base project directory
MAIN_DIR="/home/charikonst_compgen2025/Practical_5"

# Invoke the parser
"$MAIN_DIR/Scripts/Parse_Paranoid.py" \
  --inp_dir  "$MAIN_DIR/Inparanoid_Results" \
  --ref_species "9.fa.pfa" \
  --out_pairs "$MAIN_DIR/Inparanoid_Results/Ortholog_pairs.tsv" \
  --out_clusters "$MAIN_DIR/Inparanoid_Results/Paranoid_clusters.tsv"