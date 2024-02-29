library(tidyverse)
library(ape)
library(pegas)
library(adegenet)
library(mmod)
library(poppr)

comps = read_tsv('PairwiseComparisons.list',col_names = F)

#save results here 
fdM = NULL

#for chr_MT
for (comp in comps$X1) {
  cat('Working on comparison: ',comp,'\n')
  
  #p1 and p2
  p1 = gsub('__.*','',comp)
  p2 = gsub('.*__','',comp)
  
  seqs1 = read.dna(paste0('fastas/',p1,'.fa'), format = "fasta")
  seqs2 = read.dna(paste0('fastas/',p2,'.fa'), format = "fasta")
  
  combined_seqs = rbind(seqs1,seqs2)
  pops = c(rep(p1,nrow(seqs1)),rep(p2,nrow(seqs2)))
  
  #import to genind, then convert to genclone 
  binary_matrix <- as.matrix(combined_seqs)
  genind_obj <- DNAbin2genind(combined_seqs, pop = pops)
  
  gc = as.genclone(genind_obj)

  #add strata (populations)
  stratadf = data.frame(ind = indNames(genind_obj), stratum = pops)
  strata(gc) = stratadf
  res = poppr.amova(gc,hier=~stratum,clonecorrect = FALSE)
  resd = data.frame(Group = paste0(p1,'__',p2),PHI = res$statphi[[1]])
  
  fdM = rbind(fdM,resd)
}

###for W
fd = NULL
for (comp in comps$X1) {
  cat('Working on comparison: ',comp,'\n')
  
  p1 = gsub('__.*','',comp)
  p2 = gsub('.*__','',comp)
  
  seqs1 = read.dna(paste0('fastas/',p1,'.W.fa'), format = "fasta")
  seqs2 = read.dna(paste0('fastas/',p2,'.W.fa'), format = "fasta")
  
  combined_seqs = rbind(seqs1,seqs2)
  pops = c(rep(p1,nrow(seqs1)),rep(p2,nrow(seqs2)))
  
  binary_matrix <- as.matrix(combined_seqs)
  genind_obj <- DNAbin2genind(combined_seqs, pop = pops)

  gc = as.genclone(genind_obj)
  stratadf = data.frame(ind = indNames(genind_obj), stratum = pops)
  strata(gc) = stratadf
  res = poppr.amova(gc,hier=~stratum,clonecorrect = FALSE)
  resd = data.frame(Group = paste0(p1,'__',p2),PHI = res$statphi[[1]])
  
  fd = rbind(fd,resd)
}

dats = cbind(fd %>% dplyr::rename(PHI_W = PHI),fdM %>% 
        select(-Group)) 
#dats %>% ggplot(aes(x=PHI_W,y=PHI))+
  # geom_smooth(method='lm')+
  # geom_point()+
  # theme_bw()

write_tsv(dats,file='PHI_chr_MT-W_2024FEB28.txt')