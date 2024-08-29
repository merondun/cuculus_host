# Distance-based Redundancy Analysis of gens 

Determine if e.g. dbrda(Egg ~ Geography), using continuous distances. 

Inputs are autosomal distance, mtDNA fasta to calculate distance, and latitude/longitude from metadata.

`autos_canorus_LD.pdist`: Autosomal distance matrix, calcuated as `~/modules/VCF2Dis-1.50/bin/VCF2Dis -InPut autos_canorus_LD.vcf.gz -OutPut autos_canorus_LD.pdist `. 

`chr_MT_All.SNP.DP3-AC1-MQ40.min4.fasta`: mtDNA fasta for all samples.

`20240828_dbRDA_Results.txt`: Results. 

