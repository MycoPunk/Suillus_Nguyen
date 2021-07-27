#!/usr/bin/bash
#SBATCH --ntasks 10 --mem 18G --time 5:00:00 -p short -N 1

module load hmmer/3
module unload miniconda2/4.4.10
module load python/3
module unload perl
module load parallel

if [ ! -f config.txt ]; then
	echo "Need config.txt for PHYling"
	exit
fi

source config.txt
if [ ! -z $PREFIX ]; then
	rm -rf aln/$PREFIX
fi

./PHYling_unified/PHYling init
./PHYling_unified/PHYling search 
./PHYling_unified/PHYling aln -c 
pushd phylo
sbatch -p short --time 2:00:00 fast_run.sh
