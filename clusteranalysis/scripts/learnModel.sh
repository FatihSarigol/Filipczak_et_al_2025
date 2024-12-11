#!/bin/bash
#
#SBATCH --job-name=learn
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --output=learn.out
#SBATCH --mail-type=ALL
#SBATCH --qos=medium
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load java/1.8.0_212

for binSize in 5000 20000 50000 100000 500000;
do
	for design in `ls designs/*`;
	do
		for nstates in 5 6 7 8 9 10;
		do
			java -mx4000M -jar /users/daniel.malzl/ChromHMM/ChromHMM.jar LearnModel \
				-b ${binSize} \
				-p 8 \
				-r 1000 \
				-noenrich \
				-l resource/grcm38_chromsizes.tsv \
				-printstatebyline \
				binaries/`basename ${design%.*}`_binary_${binSize} \
				models/`basename ${design%.*}`_model_${binSize}_${nstates} \
				$nstates \
				mm10
		done
	done
done
