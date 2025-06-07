#!/usr/bin/env bash
set -euo pipefail

## Set up the directory paths
PFA_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/PFAs"
DB_BASE="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Databases"
REF_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Protein_seqs"
RESULT_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Blast_Results"

## Cutoffs ##
EVAL=1e-5
MIN_COV=0.8    

cd "$PFA_DIR"

for qry in *_orfs.pfa; do
  genome="${qry%_orfs.pfa}"          # e.g. "1"
  db_dir="${DB_BASE}/${genome}_db"    # e.g. ".../Databases/1_db"
  db_name="${db_dir}/${genome}"    # BLAST DB prefix
  ref_faa="${REF_DIR}/${genome}.fa.pfa"

  echo "=== Processing genome $genome ==="

  # designate output paths
  raw_out="${RESULT_DIR}/${genome}_raw.tsv"
  hits_out="${RESULT_DIR}/${genome}_hits.tsv"
  summary_out="${RESULT_DIR}/${genome}_summary.txt"

  # 1) BLASTP: predictions vs. reference database
  blastp \
    -query "$qry" \
    -db   "$db_name" \
    -evalue $EVAL \
    -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" \
    > "$raw_out"

#   # 2) filter by coverage on both query & subject
# awk -v C="$MIN_COV" '{ qcov = $3/$4; scov = $3/$5; if (qcov>=C && scov>=C) print }' \
#     "$raw_out" > "$hits_out"

#   # 3) count TP, FP, FN
#   total_q=$(grep -c '^>' "$qry")
#   tp=$(cut -f1 "$hits_out" | sort -u | wc -l)
#   fp=$(( total_q - tp ))

#   total_ref=$(grep -c '^>' "$ref_faa")
#   recov=$(cut -f2 "$hits_out" | sort -u | wc -l)
#   fn=$(( total_ref - recov ))

#   # 4) precision, recall, F1
#   prec=$(awk -v t=$tp -v f=$fp 'BEGIN{print (t/(t+f))}')
#   rec=$(awk -v t=$tp -v m=$fn 'BEGIN{print (t/(t+m))}')
#   f1=$(awk -v p=$prec -v r=$rec 'BEGIN{print (2*p*r/(p+r))}')

#   # 5) write summary
#   {
#     echo "Genome: $genome"
#     echo "Predicted ORFs: $total_q"
#     echo "Reference proteins: $total_ref"
#     echo "TP=$tp, FP=$fp, FN=$fn"
#     printf "Precision=%.3f, Recall=%.3f, F1=%.3f\n" "$prec" "$rec" "$f1"
#   } > "$summary_out"

  echo
done