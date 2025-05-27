from matplotlib import pyplot as plt
from matplotlib_venn import venn2

def load_refs_simple(path):
    refs = set()
    with open(path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            refs.add(line.split('\t')[0])
    return refs

def load_refs_paranoid(path):
    refs = set()
    with open(path) as f:
        next(f)
        for line in f:
            line=line.strip()
            if line:
                refs.add(line.split('\t')[0])
    return refs

bbh_refs = load_refs_simple("/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Practical_5/BBH_clusters.tsv")
inp_refs = load_refs_paranoid("/Users/charitidisk/MTLS/SU/Comparative_Genomics/Exercises/Practical_5/Inparanoid_Results/Paranoid_clusters.tsv")

plt.figure()
venn2([bbh_refs, inp_refs], set_labels=('BBH clusters', 'InParanoid clusters'))
plt.title('Overlap of Ortholog Clusters by Method')
plt.show()
