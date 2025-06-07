#!/usr/bin/env python3
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# adjust these two paths to where your .fa files live and where you want the .pfa files
in_dir  = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Exercises/Project/ORFs")
out_dir = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Exercises/Project/PFAs")

for fn in os.listdir(in_dir):
    if not fn.endswith("_orfs.fa"):
        continue

    in_path  = os.path.join(in_dir,  fn)
    out_path = os.path.join(out_dir, fn.replace("_orfs.fa", "_orfs.pfa"))

    prot_records = []
    for nuc_rec in SeqIO.parse(in_path, "fasta"):
        # translate up to the first stop codon, drop the '*' at end
        aa_seq = nuc_rec.seq.translate(table=11, to_stop=True)
        prot_records.append(
            SeqRecord(
                aa_seq,
                id=nuc_rec.id,
                description=""    # you can copy nuc_rec.description if you like
            )
        )

    SeqIO.write(prot_records, out_path, "fasta")
    print(f"Wrote {len(prot_records)} proteins to {out_path}")

