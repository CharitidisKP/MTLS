#!/usr/bin/env python3
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

## Set up the directory paths ##
in_dir  = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Exercises/Project/ORFs")
out_dir = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Exercises/Project/PFAs")

## Loop over the ORF files ##
for fn in os.listdir(in_dir):
    if not fn.endswith("_orfs.fa"):
        continue

    in_path  = os.path.join(in_dir,  fn)
    out_path = os.path.join(out_dir, fn.replace("_orfs.fa", "_orfs.pfa"))

    prot_records = []
    ## Translate each instance in the fastas ##
    for nuc_rec in SeqIO.parse(in_path, "fasta"):
        ## Translate up to the first stop codon ##
        aa_seq = nuc_rec.seq.translate(table=11, to_stop=True)
        prot_records.append(
            SeqRecord(
                aa_seq,
                id=nuc_rec.id,
                description=""    
            )
        )

    ## Save the pfa files ##
    SeqIO.write(prot_records, out_path, "fasta")
    print(f"Wrote {len(prot_records)} proteins to {out_path}")

