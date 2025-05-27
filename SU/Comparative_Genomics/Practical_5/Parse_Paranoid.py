#!/usr/bin/env python3
import argparse
import glob
import os
import itertools

def parse_inparanoid_file(path):
    """
    Parse one InParanoid output file.
    Returns list of tuples: (group_id, species, protein_name)
    Assumes columns: Group-id, Score, Species, Confidence-score, Protein-name
    """
    entries = []
    with open(path) as f:
        for line in f:
            cols = line.rstrip().split()
            if len(cols) < 5 or cols[0].startswith("#"):
                continue
            group_id = cols[0]
            species = cols[2]
            protein = cols[4]
            entries.append((group_id, species, protein))
    return entries

def extract_pairs(entries):
    """
    Given a list of (group, species, protein),
    group by group_id and yield all inter-species pairs:
    (group, sp1, prot1, sp2, prot2)
    """
    pairs = []
    # group entries by group_id
    by_group = {}
    for gid, sp, pr in entries:
        by_group.setdefault(gid, []).append((sp, pr))
    # for each group, form all unordered pairs across species
    for gid, members in by_group.items():
        for (sp1, pr1), (sp2, pr2) in itertools.combinations(members, 2):
            if sp1 != sp2:
                pairs.append((gid, sp1, pr1, sp2, pr2))
    return pairs

def main():
    p = argparse.ArgumentParser(
        description="Parse InParanoid tables into ortholog pairs and clusters"
    )
    p.add_argument("--inp_dir", required=True,
                   help="Directory with InParanoid output files")
    p.add_argument("--ref_species", required=True,
                   help="Reference species identifier (must match the 'Species' field in files)")
    p.add_argument("--out_pairs", default="ortholog_pairs.tsv",
                   help="Output TSV of all inter-species ortholog pairs")
    p.add_argument("--out_clusters", default="clusters.tsv",
                   help="Output TSV of reference→comma-list-of-query clusters")
    args = p.parse_args()

    ## Parse entries from all files ##
    all_entries = []
    pattern = os.path.join(args.inp_dir, "*")
    for path in glob.glob(pattern):
        if os.path.isdir(path):
            continue
        all_entries.extend(parse_inparanoid_file(path))

    ## Extract inter-species pairs ##
    pairs = extract_pairs(all_entries)

    ## Write out all pairs ##
    with open(args.out_pairs, "w") as f:
        f.write("GroupID\tSpecies1\tProtein1\tSpecies2\tProtein2\n")
        for gid, sp1, pr1, sp2, pr2 in pairs:
            f.write(f"{gid}\t{sp1}\t{pr1}\t{sp2}\t{pr2}\n")
    print(f"[ok] Wrote {len(pairs)} pairs to {args.out_pairs}")

    ## Cluster around the reference species ##
    clusters = {}
    for gid, sp1, pr1, sp2, pr2 in pairs:
        if sp1 == args.ref_species:
            clusters.setdefault(pr1, set()).add(pr2)
        if sp2 == args.ref_species:
            clusters.setdefault(pr2, set()).add(pr1)

    ## Write clusters ##
    with open(args.out_clusters, "w") as f:
        f.write("ReferenceProtein\tQueryProteins\n")
        for ref, qset in sorted(clusters.items()):
            qlist = ",".join(sorted(qset))
            f.write(f"{ref}\t{qlist}\n")
    print(f"[ok] Wrote {len(clusters)} clusters to {args.out_clusters}")

if __name__ == "__main__":
    main()