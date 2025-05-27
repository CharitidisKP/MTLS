#!/usr/bin/env bash
set -euo pipefail

SCRIPT="/home/charikonst_compgen2025/Practical_5/Scripts/Extract_BBH.py"
OUTDIR="/home/charikonst_compgen2025/Practical_5/BBH_Results"

# Reference is "9", queries are 15, 21 and 25
for Q in 15 21 25; do
  Q2R=../Pairwise_Blastp_Results/${Q}_as_Query/${Q}_vs_9_blast.tsv
  R2Q=../Pairwise_Blastp_Results/9_as_Query/9_vs_${Q}_blast.tsv
  OUT=$OUTDIR/Query_${Q}_vs_Ref.tsv

  echo ">> processing Q=$Q"
  "$SCRIPT" "$Q2R" "$R2Q" "$OUT"
done