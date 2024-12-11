for bam in `ls /groups/pavri/bioinfo/daniel/lamina/RNA/*rmdups.bam`;
do 
	condition=$( basename ${bam%.rmdups.bam} )
	echo $condition
	ln -s $bam renamed_bams/rnaseq_${condition}.bam
	ln -s ${bam}.bai renamed_bams/rnaseq_${condition}.bam.bai
done

n=$( wc -l chip_sample_id_map.txt | cut -d ' ' -f 1 )
for i in `seq 1 $n`;
do
	line=$( cat chip_sample_id_map.txt | awk -v var=$i 'NR==var {print}' ) 
	IFS=',' read sample_id condition protein cycles <<< "${line}"

	ln -s /groups/pavri/bioinfo/daniel/lamina/CHIP${cycles}cyc/${sample_id}_Filtered.rmdups.bam renamed_bams/chipseq_${protein}_${condition}_${cycles}.bam
	ln -s /groups/pavri/bioinfo/daniel/lamina/CHIP${cycles}cyc/${sample_id}_Filtered.rmdups.bam.bai renamed_bams/chipseq_${protein}_${condition}_${cycles}.bam.bai
done

ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210460_S41_R1_001.sort.rmdups.bam renamed_bams/atacseq_WT1.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210460_S41_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_WT1.bam.bai
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210461_S42_R1_001.sort.rmdups.bam renamed_bams/atacseq_WT2.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210461_S42_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_WT2.bam.bai
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210462_S43_R1_001.sort.rmdups.bam renamed_bams/atacseq_WT3.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210462_S43_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_WT3.bam.bai
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210463_S44_R1_001.sort.rmdups.bam renamed_bams/atacseq_KO1.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210463_S44_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_KO1.bam.bai
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210464_S45_R1_001.sort.rmdups.bam renamed_bams/atacseq_KO2.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210464_S45_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_KO2.bam.bai
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210465_S46_R1_001.sort.rmdups.bam renamed_bams/atacseq_KO3.bam
ln -s /groups/pavri/bioinfo/daniel/lamina/ATAC/210465_S46_R1_001.sort.rmdups.bam.bai renamed_bams/atacseq_KO3.bam.bai
