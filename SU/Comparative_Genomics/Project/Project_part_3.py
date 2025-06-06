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

# Dictionary to store sequences: key = filename base, value = sequence string
seq_dict = {}

for filename in os.listdir(dir_path):
    if filename.endswith(".fa") or filename.endswith(".fasta"):
        base_name, _ = os.path.splitext(filename)
        full_path = os.path.join(dir_path, filename)
        if os.path.isfile(full_path):
            seq_dict[base_name] = load_fasta_sequence(full_path)

# Now print each sequence length
for key, seq in seq_dict.items():
    print(f"{key}: {len(seq)} bases")




