library(tidyverse)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
md = md %>% drop_na(Egg) %>% mutate(Group = paste0(Egg,'_',Hap))
keep = md %>% count(Group) %>% filter(n >= 4)
mds = md %>% filter(Group %in% keep$Group)
mds %>% count(Group)

#save metadata, and a file .list for each pop listing the individuals 
write.table(mds,file='relatives_hunt/Metadata_Relatives_Hunt-Mirror.txt',quote=F,sep='\t',row.names=F)
for (group in unique(mds$Group)){
  d = mds %>% filter(Group == group)
  write.table(d$ID,file=paste0('relatives_hunt/pops/',group,'.list'),quote=F,sep='\t',row.names=F,col.names=F)
}
groups = mds %>% select(Group) %>% unique %>% pull(Group)
pairwise_combinations <- combn(groups, 2)

# Convert the combinations into a dataframe
pairwise_combinations_df <- data.frame(
  Group1 = pairwise_combinations[1,],
  Group2 = pairwise_combinations[2,]
)

#and add a phylo category
phylo = pairwise_combinations_df %>% mutate(Group = paste0(Group1,'__',Group2))

phylo = phylo %>% mutate(
  #for each comparison, indicate whether it is ancient comparison (e.g. blue vs noneblue) or intra-clade
  Depth = ifelse(grepl('W1|W2|W3',Group1) & !grepl('W1|W2|W3',Group2),'Ancient',
                 ifelse(grepl('W1|W2|W3',Group2) & !grepl('W1|W2|W3',Group1),'Ancient','Contemporary')),
  Comparison = ifelse(grepl('W1|W2',Group) & Depth == 'Ancient','Blue',
                      ifelse(Depth == 'Ancient','Reversion',
                             ifelse(Group == 'W2__Pphoenicurus__W1__Pphoenicurus' | Group == 'W7__Aarundinaceus__W5__Aarundinaceus','Control',
                                    ifelse((grepl('W1|W2|W3',Group1) & grepl('W1|W2|W3',Group2)) & Depth == 'Contemporary','Reversion','Diversification')))),
  Phylo = paste0(Depth,'_',Comparison))
phylo %>% count(Phylo)

write.table(phylo,'relatives_hunt/Comparisons_MetadataMirror.txt',quote=F,sep='\t',row.names=F)
write.table(mds %>% select(ID,Group),'relatives_hunt/Relatives_HuntMirror.pop',quote=F,sep='\t',row.names=F,col.names=F)

dfs = phylo %>% mutate(Group = paste0(Group1,'@',Group2))
write.table(dfs$Group,file='relatives_hunt/Pairwise_ComparisonsMirror.list',quote=F,sep='\t',row.names=F,col.names=F)
