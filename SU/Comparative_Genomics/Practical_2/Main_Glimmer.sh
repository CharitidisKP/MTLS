#!/usr/bin/env bash
set -euo pipefail

# define input / output
INPUT="$HOME/Practical_2/Fasta_Files"
RESULTS="$HOME/Practical_2/Glimmer_Results"

# make top-level results dir
mkdir -p "$RESULTS"

# go to where your .fa files live
cd "$INPUT"

# loop!
for fa in *.fa; do
  # strip off the .fa to get your sample prefix
  base="${fa%.fa}"

  # make a per-sample directory under RESULTS
  outdir="$RESULTS/$base"
  mkdir -p "$outdir"

  # 1) find long ORFs
  tigr-glimmer long-orfs -n -t 1.15 \
      "$fa" \
      "$outdir/${base}.long-orf-coords"

  # 2) extract the ORF sequences
  tigr-glimmer extract -t \
      "$fa" \
      "$outdir/${base}.long-orf-coords" \
    > "$outdir/${base}.longorf"

  # 3) build the ICM model
  tigr-glimmer build-icm -r \
      "$outdir/${base}.icm" \
    < "$outdir/${base}.longorf"

  # 4) run Glimmer3 gene prediction
  tigr-glimmer glimmer3 -o50 -g110 -t30 \
      "$fa" \
      "$outdir/${base}.icm" \
      "$outdir/${base}.glimmer"

  echo "✓ $base done ⇒ $outdir/"
done

echo "All samples processed. Check your results under $RESULTS/" 