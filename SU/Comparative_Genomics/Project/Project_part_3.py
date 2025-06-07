#! /usr/bin/env python3
import os 
from textwrap import wrap

## Parameters ##
Dir_Path = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Fasta_Files/")
Out_Path = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Exercises/Project/ORFs")

if not os.path.isdir(Out_Path):
    os.makedirs(Out_Path)

Min_AA_len = 35
Stop_Codons = {"TAA", "TAG", "TGA"}
Start_Codon = {"ATG"}

## Load the data ##
def Load_Fasta(Path):
    Sub_seqs = []
    with open(Path) as current_file:
        for line in current_file:
            if line.startswith(">"):
                continue
            Sub_seqs.append(line.strip())
    return "".join(Sub_seqs).upper()            

## Get the reverse complement ##
def Reverse_Comp(Seq):
    table = str.maketrans("ACGTacgtnN", "TGCAtgcanN")
    return Seq.translate(table)[::-1]

## Function to get the ORFs ##
def ORF_Sniffer(Seq, Frame, Min_len):
    ORFs = []
    i = Frame
    L = len(Seq)


    while i + 3 <= L:
        codon = Seq[i : i + 3]
        if codon == "ATG":

            j = i + 3
            while j + 3 <= L:
                stop_codon = Seq[j : j + 3]
                if stop_codon in Stop_Codons:
                    orf_nt_len = (j + 3) - i
                    if (orf_nt_len // 3) >= Min_len:
                        ORFs.append((i, j + 3, Seq[i : j + 3]))
                    break
                j += 3
            i += 3
        else:
            i += 3

    return ORFs
    
## Retrieve the ORF from the sequences ##
def ORF_Getter(Seq):
    Min_Nt_len = Min_AA_len * 3
    L = len(Seq)

    ## Main strand ##
    Forward_ORFs = []
    for Frame in (0, 1, 2): 
        hits = ORF_Sniffer(Seq, Frame, Min_Nt_len)

        for (s, e, nt_seq) in hits:
            Forward_ORFs.append({
                "Start": s, 
                "End": e, 
                "Strand": "Forward", 
                "Sequence": Seq[s:e]
            } )

    ## Reverse strand ##
    Reverse_ORFs = []
    Rev = Reverse_Comp(Seq)

    for Frame in (0, 1, 2): 
        hits = ORF_Sniffer(Rev, Frame, Min_Nt_len)

        for (s_r, e_r, nt_seq_r) in hits:
            orig_start = L - e_r
            orig_end = L - s_r
            Reverse_ORFs.append({
                "Start": orig_start,
                "End": orig_end,
                "Strand": "Reverse",
                "Sequence": nt_seq_r
            })

    return Forward_ORFs, Reverse_ORFs

## Filter the ORFs based on nesting ##
def Filter_ORFs(ORF_list):
    ORF_list_sorted = sorted(ORF_list, key = lambda x: (x["End"] - x["Start"]), reverse = True)
    kept = []

    for Current_ORF in ORF_list_sorted:
        s_c, e_c = Current_ORF["Start"], Current_ORF["End"]
        is_Nested = False
        for bigger in kept: 
            s_b, e_b = bigger["Start"], bigger["End"]
            if (s_b <= s_c) and (e_c <= e_b): 
                is_Nested = True
                break 
        if not is_Nested:
            kept.append(Current_ORF)
    return kept

## Initiate main ##
seq_dict = {}

for filename in os.listdir(Dir_Path):
    if filename.lower().endswith((".fa", ".fasta")):
        base_name, _ = os.path.splitext(filename)
        full_path = os.path.join(Dir_Path, filename)

        if os.path.isfile(full_path):
            seq = Load_Fasta(full_path)
            seq_dict[base_name] = seq

for genome_name, dna_seq in seq_dict.items():
    f_orfs, r_orfs = ORF_Getter(dna_seq)

    f_filtered = Filter_ORFs(f_orfs)
    r_filtered = Filter_ORFs(r_orfs)

    all_orfs = sorted(f_filtered + r_filtered, key=lambda x: (x["Start"], x["End"], x["Strand"]))

    out_fname = os.path.join(Out_Path, f"{genome_name}_orfs.fa")
    with open(out_fname, "w") as outfh:
        for idx, orf in enumerate(all_orfs, start = 1):
            s, e, strand = orf["Start"], orf["End"], orf["Strand"]

            ## format: genome_ORF0001_123-456_forward ##
            strand_label = "forward" if strand == "Forward" else "reverse"
            header = f">{genome_name}_ORF{idx:04d}_{s+1}-{e}_{strand_label}"
            outfh.write(header + "\n")

            ## Wrap the sequence to 60 nt/line ##
            for line in wrap(orf["Sequence"], 60):
                outfh.write(line + "\n")

    print(f"Wrote {len(all_orfs)} ORFs to {out_fname}")

