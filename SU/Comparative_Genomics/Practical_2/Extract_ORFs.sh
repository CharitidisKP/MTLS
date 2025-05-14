#!/usr/bin/env bash
set -euo pipefail

INPUT="$HOME/Practical_2/Fasta_Files"
GLIMMER_DIR="$HOME/Practical_2/Glimmer_Results"
OUTPUT="$HOME/Practical_2/ORF_Seqs"
SCRIPTS="$HOME/Practical_2/Scripts"

mkdir -p "$OUTPUT"

for fa in "$INPUT"/*.fa; do
  base=$(basename "$fa" .fa)

  glim="$GLIMMER_DIR/$base/${base}.glimmer.predict"
  if [[ ! -f "$glim" ]]; then
    echo "Prediction file not found for $base, skipping."
    continue
  fi

  # make an output folder for this sample
  outdir="$OUTPUT/$base"
  mkdir -p "$outdir"

  # run your extraction script
  python3 $SCRIPTS/parseGlimmer.py \
    "$fa" \
    "$glim" \
    --translate ## --translate for the protein output ##

  # move whatever .nfa or .pfa files it produced into the sample folder
  mv *.nfa *.pfa "$outdir"/

  echo "Extracted ORFs for $base → $outdir/"
done

echo "All done!  Check $OUTPUT."