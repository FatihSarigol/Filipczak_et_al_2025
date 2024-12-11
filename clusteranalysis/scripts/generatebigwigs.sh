#!/bin/bash
#
#SBATCH --job-name=bigwiggen
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --output=bigwiggen.out
#SBATCH --mail-type=ALL
#SBATCH --time=00-08:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

for bam in `ls renamed_bams/*bam`;
do
	for bs in 5000 20000 50000 100000 500000;
	do
		echo processing ${bam}
		bamCoverage -b ${bam} -o bigwigs/`basename ${bam%.*}`_${bs}.bw -of bigwig -bs $bs -p 8 --ignoreDuplicates 
	done
done
