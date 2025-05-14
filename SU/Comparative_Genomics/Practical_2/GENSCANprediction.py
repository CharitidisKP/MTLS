import re
import argparse

def extract_sequences(genscan_output_file, aa_output_file, nt_output_file):
    with open(genscan_output_file, 'r') as infile:
        content = infile.read()

    # Regular expressions to match amino acid and nucleotide entries
    aa_entries = re.findall(
        r'>([^|]+)\|GENSCAN_predicted_peptide_\d+\|(\d+)_aa\n([A-Z\n]+)', content)
    nt_entries = re.findall(
        r'>([^|]+)\|GENSCAN_predicted_CDS_\d+\|(\d+)_bp\n([acgt\n]+)', content, re.IGNORECASE)

    # Write amino acid sequences
    with open(aa_output_file, 'w') as aa_out:
        for seq_id, length, seq in aa_entries:
            seq_clean = ''.join(seq.split())
            aa_out.write(f'>{seq_id}|{length}_aa\n')
            for i in range(0, len(seq_clean), 60):
                aa_out.write(seq_clean[i:i+60] + '\n')

    # Write nucleotide sequences
    with open(nt_output_file, 'w') as nt_out:
        for seq_id, length, seq in nt_entries:
            seq_clean = ''.join(seq.split())
            nt_out.write(f'>{seq_id}|{length}_bp\n')
            for i in range(0, len(seq_clean), 60):
                nt_out.write(seq_clean[i:i+60] + '\n')

    print(f"Amino acid sequences written to: {aa_output_file}")
    print(f"Nucleotide sequences written to: {nt_output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract AA and NT sequences from GENSCAN output")
    parser.add_argument("genscan_output_file", help="Input GENSCAN output file")
    parser.add_argument("aa_output_file", help="Output file for amino acid sequences")
    parser.add_argument("nt_output_file", help="Output file for nucleotide sequences")
    args = parser.parse_args()

    extract_sequences(args.genscan_output_file, args.aa_output_file, args.nt_output_file)
