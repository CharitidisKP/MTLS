#!/usr/bin/env python3
import sys

def parse_best_hits(blast_file):
    """
    From an outfmt6 BLASTP file, keep for each subject (col2) the one
    query (col1) with smallest E-value (col11) and, if tied, highest bitscore (col12).
    Returns dict: subject_id -> best_query_id
    """
    best_hits = {}
    with open(blast_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            qid, sid = cols[0], cols[1]
            evalue, bitscore = float(cols[10]), float(cols[11])
            # first hit or strictly better?
            if sid not in best_hits or \
               (evalue < best_hits[sid][1]) or \
               (evalue == best_hits[sid][1] and bitscore > best_hits[sid][2]):
                best_hits[sid] = (qid, evalue, bitscore)
    # collapse to subject -> query
    return {sid: data[0] for sid, data in best_hits.items()}

def main():
    if len(sys.argv) != 4:
        print("Usage: python bbh_parser.py <Q_vs_R.tsv> <R_vs_Q.tsv> <output.tsv>")
        sys.exit(1)

    q2r_file, r2q_file, out_file = sys.argv[1:]
    # parse best hits in each direction
    q2r = parse_best_hits(q2r_file)  # R_subject -> Q_query
    r2q = parse_best_hits(r2q_file)  # Q_subject -> R_query

    # find BBHs
    bbh = []
    for r_subj, q_query in q2r.items():
        if r2q.get(q_query) == r_subj:
            bbh.append((r_subj, q_query))

    # write
    with open(out_file, 'w') as out:
        out.write("Reference\tQuery\n")
        for r, q in sorted(bbh):
            out.write(f"{r}\t{q}\n")

    print(f"[bbh_parser] {len(bbh)} BBH pairs written to '{out_file}'")

if __name__ == "__main__":
    main()