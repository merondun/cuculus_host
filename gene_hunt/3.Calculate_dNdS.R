#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(GenomicFeatures)

chr = args[1]
txdb = makeTxDbFromGFF("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR-chr_MT.gff")
fasta_seq <- readDNAStringSet("/dss/dsslegfs01/pr53da/pr53da-dss-0021/assemblies/Cuculus.canorus/VGP.bCucCan1.pri/GCA_017976375.1_bCucCan1.pri_genomic.CHR.fa")

#load in vcf
vcf = readVcf(paste0("vcfs/",chr,'.SNP.DP3-AC1-MQ40.vcf.gz'), genome = "cuckoo")

#annotate
if (chr == 'chr_MT') {
  code = getGeneticCode("SGC1")
  coding_predictions = predictCoding(query = vcf, subject = txdb, seqSource = fasta_seq, genetic.code= code) #for chr_MT
  dat = as_tibble(coding_predictions) %>% dplyr::select(seqnames,start,REF,varAllele,CONSEQUENCE,GENEID) 
} else {
  code = getGeneticCode("SGC0")
  coding_predictions = predictCoding(query = vcf, subject = txdb, seqSource = fasta_seq, genetic.code= code,) #for other chrs
  dat = as_tibble(coding_predictions) %>% dplyr::select(seqnames,start,REF,varAllele,CONSEQUENCE,GENEID) 

}

write.table(dat %>% unique %>% arrange(seqnames,start),file=paste0('dnds/Annotated_Variants_',chr,'__2024MAR01.txt'),quote=F,sep='\t',row.names=F)