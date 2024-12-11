#!/bin/bash
#
#SBATCH --job-name=binarization
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --output=binarization.out
#SBATCH --mail-type=ALL
#SBATCH --qos=medium
#SBATCH --time=01-00:00:00
#SBATCH --mail-user=daniel.malzl@imp.ac.at

module load java/1.8.0_212

for binSize in 5000 20000 50000 100000 500000;
do
	for design in `ls designs/*`;
	do
		java -mx4000M -jar /users/daniel.malzl/ChromHMM/ChromHMM.jar BinarizeBam \
			-mixed \
			-b ${binSize} \
			resource/grcm38_chromsizes.tsv renamed_bams \
			$design \
			binaries/`basename ${design%.*}`_binary_${binSize}
	done
done
