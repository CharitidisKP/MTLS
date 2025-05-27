#!/usr/bin/env bash
set -euo pipefail

## Paths for the pfa files, the databases and the output results ##
PROTEOME_DIR="/home/charikonst_compgen2025/Practical_5/Protein_seqs"
DB_ROOT="/home/charikonst_compgen2025/Practical_5/Databases"
OUT_DIR="/home/charikonst_compgen2025/Practical_5/Pairwise_Blastp_Results"

## Loop all proteome files as queries ##
for query_fa in "$PROTEOME_DIR"/*.pfa; do
  qbase=$(basename "$query_fa" .fa.pfa)

  ## Loop all proteome database folders ##
  for dbdir in "$DB_ROOT"/*_db; do
    dbbase=$(basename "$dbdir" _db)

    # Optional: skip self-vs-self
    if [[ "$qbase" == "$dbbase" ]]; then
      continue
    fi

    ## Prepare subdirectories for the pairwise comparisons ##
    OUT_sub_DIR="$OUT_DIR/${qbase}_as_Query"
    ## Make them ##
    mkdir -p "$OUT_sub_DIR"

    ## Prepare output names and paths ##
    out_tsv="$OUT_sub_DIR/${qbase}_vs_${dbbase}_blast.tsv"
    ## Sanity check ##
    echo "Files $qbase vs $dbbase have been placed in: $out_tsv"

    ## Run the BLASTp command ##
    blastp \
      -query  "$query_fa" \
      -db     "$dbdir/$dbbase" \
      -outfmt 6 \
      -out    "$out_tsv"
  done
done

## Sanity check number 2 ##
echo "Pairwise BLASTP results in $OUT_DIR/"
