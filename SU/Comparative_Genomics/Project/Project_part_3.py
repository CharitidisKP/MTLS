import os

def read_fasta(path):
    """
    Read a FASTA file and return a list of tuples: (header, sequence).
    Removes whitespace so that each sequence is continuous.
    """
    records = []
    with open(path) as fh:
        header = None
        seq_lines = []
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                # Save the previous record
                if header is not None:
                    sequence = "".join(seq_lines).replace(" ", "").upper()
                    records.append((header, sequence))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        # Don’t forget the last record
        if header is not None:
            sequence = "".join(seq_lines).replace(" ", "").upper()
            records.append((header, sequence))
    return records


Fasta_Path = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Fasta_Files/")
my_file = "1.fa"

full_path = os.path.join(Fasta_Path, my_file)
if not os.path.isfile(full_path):
    raise FileNotFoundError(f"File not found: {full_path}")

Yeet = read_fasta(full_path)


import os

def load_fasta_sequence(path):
    """
    Read a FASTA file and return the concatenated sequence (uppercase, no spaces).
    Assumes a single sequence per file.
    """
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip())
    return "".join(seq_parts).upper()

# Directory containing FASTA files
dir_path = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Fasta_Files/")

# Dictionary to store sequences: key = filename base (without extension), value = sequence string
seq_dict = {}

for filename in os.listdir(dir_path):
    if filename.endswith(".fa") or filename.endswith(".fasta"):
        base_name, _ = os.path.splitext(filename)
        full_path = os.path.join(dir_path, filename)
        if os.path.isfile(full_path):
            seq_dict[base_name] = load_fasta_sequence(full_path)

# Example usage: print base name and sequence length
for key, seq in seq_dict.items():
    print(f"{key}: {len(seq)} bases")

load_fasta_sequence(dir_path)