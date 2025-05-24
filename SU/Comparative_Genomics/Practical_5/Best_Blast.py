import sys

def parse_blast_best_hits(blast_file, output_file):
    best_hits = {}

    with open(blast_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue  # Skip comments or empty lines

            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue  # Skip malformed lines

            query_id = parts[0]
            subject_id = parts[1]
            evalue = float(parts[10])
            bitscore = float(parts[11])

            # Store best hit for each subject (reference protein)
            if subject_id not in best_hits:
                best_hits[subject_id] = (query_id, evalue, bitscore)
            else:
                # Replace if e-value is better or equal with higher bitscore
                _, prev_evalue, prev_bitscore = best_hits[subject_id]
                if (evalue < prev_evalue) or (evalue == prev_evalue and bitscore > prev_bitscore):
                    best_hits[subject_id] = (query_id, evalue, bitscore)

    # Write results
    with open(output_file, 'w') as out:
        for subject_id in sorted(best_hits.keys()):
            query_id = best_hits[subject_id][0]
            out.write(f"{subject_id} {query_id}\n")

    print(f"Parsed best hits written to '{output_file}'.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python best_hits_parser.py <blast_output.txt> <output.txt>")
        sys.exit(1)

    blast_file = sys.argv[1]
    output_file = sys.argv[2]

    parse_blast_best_hits(blast_file, output_file)

