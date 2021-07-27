#!/usr/bin/bash
#SBATCH -p intel -n 2 -N 1 --mem 16gb --out iqtree.%A.log

module load IQ-TREE
NAMES=../jgi_names.tab
NUM=$(wc -l ../prefix.tab | awk '{print $1}')
PREF=Suillus
ALN=../$PREF.${NUM}_taxa.JGI_1086.aa.fasaln
PART=../$PREF.${NUM}_taxa.JGI_1086.aa.partitions.txt
iqtree -nt 2 -alrt 1000 -bb 1000 -s $ALN -pre $PREF.${NUM}_taxa -m TESTMERGE
