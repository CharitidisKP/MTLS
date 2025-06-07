#!/usr/bin/env bash
set -euo pipefail

## Set up the directory paths ##
PFA_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/PFAs"
DB_BASE="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Databases"
REF_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Protein_seqs"
RESULT_DIR="/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Project/Blast_Results"

## Cutoffs ##
EVAL=1e-5
MIN_COV=0.8    

cd "$PFA_DIR"

## Loop over the orf_pfas, the references and the databases ##
for qry in *_orfs.pfa; do
  genome="${qry%_orfs.pfa}"          
  db_dir="${DB_BASE}/${genome}_db"    
  db_name="${db_dir}/${genome}"    
  ref_faa="${REF_DIR}/${genome}.fa.pfa"

  echo "---------- Processing genome $genome ----------"

  ## Output paths ##
  raw_out="${RESULT_DIR}/${genome}_raw.tsv"
  hits_out="${RESULT_DIR}/${genome}_hits.tsv"
  summary_out="${RESULT_DIR}/${genome}_summary.txt"

  ## Actual blast command ##
  blastp \
    -query "$qry" \
    -db   "$db_name" \
    -evalue $EVAL \
    -outfmt "6 qseqid sseqid length qlen slen evalue bitscore" \
    > "$raw_out"

  echo
done