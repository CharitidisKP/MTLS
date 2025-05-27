#!/usr/bin/env bash
set -euo pipefail

## Directory locations ##
INP_DIR="/home/charikonst_compgen2025/Practical_5/inparanoid"
REF_FASTA="/home/charikonst_compgen2025/Practical_5/Protein_seqs/9.fa.pfa"  
QUERY_DIR="/home/charikonst_compgen2025/Practical_5/Protein_seqs"       
OUT_ROOT="/home/charikonst_compgen2025/Practical_5/Inparanoid_Results"          

## Move to the inparanoid directory ##
cd "$INP_DIR"    

## Iterate over the pfa files ##
for Q in "$QUERY_DIR"/*.pfa; do
  QBASE=$(basename "$Q" .pfa)
  ## Skip the reference ##
  if [[ "$Q" == "$REF_FASTA" ]]; then
    continue
  fi

  echo "Running InParanoid: Reference vs $QBASE"
  perl inparanoid.pl "$REF_FASTA" "$Q"

  ## InParanoid writes its tables into ./output/ ##
  ## Pick up the default table rename and move it: ##
  REFBASE=$(basename "$REF_FASTA" .pfa)

  ## Move each output file ##
  for outf in output/*; do
    fname=$(basename "$outf")
    mv "$outf" "$OUT_ROOT/$fname"
    echo "Moved $fname"
  done

  ## Clean out the output directory ##
  rm -f output/*   
done

echo "Done! All InParanoid outputs are in $OUT_ROOT/"