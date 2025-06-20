#!/usr/bin/env bash
set -euo pipefail

## Set up the directory paths ##
PFA_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/PFAs"
REF_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Protein_seqs"
RESULT_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Blast_Results"
OUT_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Comparison_Results"
# ensure it exists
mkdir -p "$RESULT_DIR"
mkdir -p "$OUT_DIR" 

## Minimum fraction of both sequences covered by the alignment ##
MIN_COV=0.8     

## Loop over each blast result file ##
for raw in "$RESULT_DIR"/*_raw.tsv; do
  genome=$(basename "$raw" _raw.tsv)          
  hits="${OUT_DIR}/${genome}_hits.tsv"
  summary="${OUT_DIR}/${genome}_summary.txt"

  echo "Evaluating genome $genome"

  ## Filter by coverage on both query and reference ##
  echo "raw = [$raw]"
  echo "hits out= [$hits]"
  awk '{ C='"$MIN_COV"'; qcov=$3/$4; scov=$3/$5; if (qcov>=C && scov>=C) print }' "$raw" > "$hits"

  ## Get TP, FP and FN ##
  total_q=$(grep -c '^>' "$PFA_DIR/${genome}_orfs.pfa")
  tp=$(cut -f1 "$hits" | sort -u | wc -l | xargs)
  fp=$(( total_q - tp ))

  total_ref=$(grep -c '^>' "$REF_DIR/${genome}.fa.pfa")
  recov=$(cut -f2 "$hits" | sort -u | wc -l | xargs)
  fn=$(( total_ref - recov ))

  ## Present the statistics in a clean way
  prec=$(awk -v "t=$tp" -v "f=$fp" 'BEGIN{print (t+f>0? t/(t+f): 0)}')
  rec=$(awk -v "t=$tp" -v "n=$fn" 'BEGIN{print (t+n>0? t/(t+n): 0)}')
  f1=$(awk -v "p=$prec" -v "r=$rec" 'BEGIN{print (p+r>0? 2*p*r/(p+r): 0)}')

  ## Write summary ##
  {
    echo "Genome: $genome"
    echo "Predicted ORFs: $total_q"
    echo "Reference proteins: $total_ref"
    echo "TP=$tp, FP=$fp, FN=$fn"
    printf "Precision=%.3f, Recall=%.3f, F1=%.3f\n" "$prec" "$rec" "$f1"
  } > "$summary"

  ## Sanity check ##
  echo "Hits: $hits"
  echo "Summary: $summary"
  echo
done