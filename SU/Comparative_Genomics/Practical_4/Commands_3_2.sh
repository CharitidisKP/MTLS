## Commands for exercise 3.2 ##
cd ~/Practical_4/Cluster_Alignments/ ## not necesary but yeah ##

~/Practical_4/Belvu_script_1.sh

nano rename.sed

## Insert the following into rename.sed without the coments: ##
## rename.sed
#s|15\.fa_orf[^:]*:([0-9]+\.[0-9]+)|C_trachomatis:\1|g
#s|9\.fa_orf[^:]*:([0-9]+\.[0-9]+)|S_epidermidis:\1|g
#s|21\.fa_orf[^:]*:([0-9]+\.[0-9]+)|S_pyogenes:\1|g
#s|25\.fa_orf[^:]*:([0-9]+\.[0-9]+)|L_gelidum:\1|g
#s|1\.fa_orf[^:]*:([0-9]+\.[0-9]+)|S_cerevisiae:\1|g

cd trees_NJ

for f in cluster_*_aligned.nwk; do   
	sed -E -f rename.sed "$f" 
> trees_named_fixed/"$f"; 
done

cd trees_named
cat cluster_* > intree
grep -v '^$' intree > intree.tmp && mv intree.tmp intree 
 
grep -c ';' intree
phylip consense

