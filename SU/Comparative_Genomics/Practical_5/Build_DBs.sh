#!/usr/bin/env bash
set -euo pipefail

## Protein files directory ##
SRC_DIR="/home/charikonst_compgen2025/Practical_5/Protein_seqs"
## Databases directory ##
DB_DIR="/home/charikonst_compgen2025/Practical_5/Databases"

for file in "$SRC_DIR"/*.pfa; do
    ## Clean up the names ##
    base=$(basename "$file" .fa.pfa)

    ## Make a subdirectory for each pfa database ##
    DB_sub_DIR="$DB_DIR/${base}_db"
    mkdir -p "$DB_sub_DIR"

    ## Create the database in each subfolder ##
    echo ">> makeblastdb -in $file -dbtype prot -out $DB_sub_DIR/$base"
    makeblastdb -in "$file" -dbtype prot -out "$DB_sub_DIR/$base"

done

echo "All databases built within $DB_DIR"