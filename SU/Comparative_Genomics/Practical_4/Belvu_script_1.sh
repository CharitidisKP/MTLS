#!/usr/bin/env bash
cd ~/Practical_4/Cluster_Alignments

# make a place to collect your Newick trees
mkdir -p trees_NJ

for aln in cluster_*_aligned.fasta; do
  base="${aln%.fasta}"
  belvu \
    -quiet \
    -o tree \
    "$aln" \
    > "trees_NJ/${base}.nwk"
done