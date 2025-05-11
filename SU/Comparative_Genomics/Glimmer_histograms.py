#!/usr/bin/env python3
"""
Plot ORF-length histograms for all .glimmer.predict files under a directory tree.
"""
import argparse
import os
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot ORF-length histograms for all .glimmer.predict files from subfolders."
    )
    parser.add_argument(
        "-i", "--input-dir",
        default=".",
        help="Root directory containing subfolders with .glimmer.predict files."
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="Plots",
        help="Directory to save histogram PNGs."
    )
    parser.add_argument(
        "-b", "--bins",
        type=int,
        default=50,
        help="Number of bins for the histogram."
    )
    parser.add_argument(
        "-m", "--max-length",
        type=float,
        default=500000,
        help="Maximum ORF length (bp) to include in the histogram."
    )
    return parser.parse_args()


def plot_histogram(input_file, output_file, bins, max_length):
    """
    Read a .glimmer.predict file, compute ORF lengths, and save a histogram.
    """
    orf_lengths = []
    with open(input_file, "r") as fh:
        # skip header
        next(fh, None)
        for line in fh:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    length = abs(float(parts[2]) - float(parts[1]))
                    if length <= max_length:
                        orf_lengths.append(length)
                except ValueError:
                    continue

    if not orf_lengths:
        print(f"No ORFs found in {input_file}, skipping.")
        return

    # plot
    plt.figure(figsize=(10, 6))
    plt.hist(orf_lengths, bins=bins)
    plt.title(f"ORF Lengths from {os.path.basename(input_file)}")
    plt.xlabel("ORF Length (bp)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Histogram saved to: {output_file}")


def main():
    args = parse_args()

    # ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # walk input tree
    for root, _, files in os.walk(args.input_dir):
        for fname in files:
            if fname.endswith(".glimmer.predict"):
                infile = os.path.join(root, fname)
                # derive sample prefix (remove .glimmer.predict)
                prefix = fname[:-len(".glimmer.predict")]
                out_name = f"ORFlength_{prefix}_hist.png"
                out_path = os.path.join(args.output_dir, out_name)
                plot_histogram(
                    infile,
                    out_path,
                    bins=args.bins,
                    max_length=args.max_length
                )


if __name__ == "__main__":
    main()
