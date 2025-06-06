#! /usr/bin/env python3
import os 
from textwrap import wrap

## Parameters ##
Dir_Path = os.path.expanduser("~/MTLS/SU/Comparative_Genomics/Fasta_Files/")

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

def ORF_Sniffer(Seq, Frame, Min_len):
    ORFs = []
    i = Frame
    L = len(Seq)

    while i + 3 <= L:
        codon = Seq[i:i+3]
        if codon in Start_Codon:
            for j in range(i, L, 3):
                stop_codon = Seq[j:j+3]
            if stop_codon in Stop_Codons:
                orf_len = j + 3 - i
                if orf_len // 3 >= Min_len:
                    ORFs.append((i, j + 3, Seq[i:j+3]))
                break
        i += 3
        return ORFs
    

def ORF_Getter(Seq, Frame, Min_Nt_len):
    Min_Nt_len = Min_AA_len * 3
    L = len(Seq)

    ## Main strand ##
    Forward_ORFs = []
    for frame in (0, 1, 2): 
        hits = ORF_Sniffer(Seq, Frame, Min_Nt_len)

        for (s, e) in hits:
            Forward_ORFs.append({
                "Start": s, 
                "End": e, 
                "Strand": "Forward", 
                "Sequence": Seq[s:e]
            } )

    ## Reverse strand ##
    Reverse_ORFs = []
    Rev = Reverse_Comp(Seq)

    for frame in (0, 1, 2): 
        hits = ORF_Sniffer(Rev, Frame, Min_Nt_len)

        for (s_r, e_r) in hits:
            Reverse_ORFs.append({
                "Start": L - s_r, 
                "End": L - e_r, 
                "Strand": "Reverse", 
                "Sequence": Seq[L - s_r:L - e_r]
            } )


def Filter_ORFs(ORF_list):
    ORF_list_sorted = sorted(orf_list, key = lambda x: (x["end"] - x["start"]), reverse = True)
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
