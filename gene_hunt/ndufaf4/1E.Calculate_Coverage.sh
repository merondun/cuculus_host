for RUN in $(ls *bam | sed 's/.bam//g'); do 

	cds=/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/cuculiformes_NDUFAF4_blasting/Alignment_Based_Approach/CDS_Regions_N3_BothParalogs.bed

	# Calculate coverage MQ50 
	mosdepth --threads 3 --mapq 50 --by ${cds} --fast-mode --no-per-base ../coverage/${RUN}_cds50 ${RUN}.bam

done 

# After, in /coverage/, merge into single file:
for i in $(ls *_cds50.regions.bed | sed 's/_cds50.regions.bed//g'); do awk -v i=${i} '{OFS="\t"}{print $1, $2, $3, $4, i}' ${i}_cds50.regions.bed > ${i}.cds50.cov ; done
cat *.cds50.cov > NDUFAF4_MQ50_Coverage_2024SEPT07.txt
# I then add e.g. species / tech in excel manually since there are not so many fields 
