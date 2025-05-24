#!/usr/bin/env python3
import argparse
import glob
import os

def load_ref2queries(bbh_dir):
    """
    Scan all .tsv in bbh_dir matching Query_*_vs_Ref.tsv,
    skip the header line, and build { ref_orf: [qry_orf, ...], ... }.
    """
    ref2q = {}
    pattern = os.path.join(bbh_dir, "Query_*_vs_Ref.tsv")
    for path in glob.glob(pattern):
        with open(path) as f:
            next(f)  # skip header
            for line in f:
                ref, qry = line.rstrip().split('\t')
                ref2q.setdefault(ref, []).append(qry)
    return ref2q

def main():
    p = argparse.ArgumentParser(
        description="Group BBHs by reference ORF"
    )
    p.add_argument(
        "bbh_dir",
        help="Directory containing your Query_*_vs_Ref.tsv files"
    )
    args = p.parse_args()

    bbh_dir = args.bbh_dir
    if not os.path.isdir(bbh_dir):
        print(f"ERROR: '{bbh_dir}' is not a directory.")
        exit(1)

    ref2q = load_ref2queries(bbh_dir)

    # Print out in grouped form
    for ref, qlist in sorted(ref2q.items()):
        print(ref, *qlist, sep="\t")

if __name__ == "__main__":
    main()