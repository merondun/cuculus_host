#### Determine egg shift parsimony using binary classifications 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/tree_comparison_mtauto/with_cp')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(ape)
library(caper)
library(ggtree)
library(treeio)
library(phytools)
library(tidyverse)
library(viridis)
library(pheatmap)
library(ggpubr)
library(meRo)
library(data.table)
library(tangler)
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(magick)

# Read in trees and metadata 
md <- read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') %>% filter(Analysis_PopulationGenetics == 1) %>% drop_na(Egg)
md$dummy <- 1 # Assign fake dummy var so that all connections are plotted 

### mtDNA tree
m = read.iqtree('ml_trees/chr_MT_n89.min4.phy.varsites.phy.contree')
m1 <- root(as.phylo(m),outgroup = '387_CP_MBW_RUS_F')
mt_tree <- drop.tip(m1,c('387_CP_MBW_RUS_F','386_CP_MBW_RUS_M'))

### AUTOSOME tree 
a = read.iqtree('ml_trees/autos_LD_n89.min4.phy.varsites.phy.contree')
a1 <- root(as.phylo(a),outgroup = '387_CP_MBW_RUS_F')
a_tree <- drop.tip(a1,c('387_CP_MBW_RUS_F','386_CP_MBW_RUS_M'))

# Store results
results <- list()
for (egg in unique(md$Egg)) { 
  for (tree in c('mt_tree','a_tree')) { 
    
    cat('Working on egg type: ',egg,' for tree: ',tree,'\n')
    
    # Change egg to binary trait, only target egg is 1 all else is 0 
    md_mod <- md %>% mutate(Egg = ifelse(Egg == egg,egg,'E0'))
    
    # Generate base tree 
    targ_tree <- get(tree)
    ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md_mod
    
    #grab only egg
    phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
    egg_mat <- as.matrix(phenos %>% select(Egg))
    phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)
    
    #inspect tree
    ggtree(targ_tree, layout = "rectangular",branch.length='none') %<+% md_mod + 
      geom_tippoint(aes(fill=Egg),pch=21,size=2)+
      scale_fill_brewer(palette='Set2')
    
    #Plot probabilities 
    t2 <- multi2di(targ_tree)
    
    # Quick ape method 
    fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")
    
    # Determine ancestral state likelihood number of shifts using MCMC
    sim <- make.simmap(t2, phenotypes, model="ER", nsim=100,Q='mcmc')
    
    # Output stats on the number of shifts across the 100 replicates 
    vec <- data.frame(shifts = countSimmap(sim)$Tr[,1])
    shifts <- vec %>% mutate(Egg = egg, Tree = tree)
    
    results[[paste0(egg,'_',tree)]] <- shifts
    
  }
}

full_results <- rbindlist(results)

write.table(full_results,file='Results_Binary_EggComparison_2024JULY16.txt',quote=F,sep='\t',row.names=F)

# Read in 
full_results <- read_tsv('Results_Binary_EggComparison_2024JULY16.txt')

# Add egg colors
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
egglevs = md %>% filter(Analysis_Mantel == 1) %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% mutate(col = viridis(11, option='turbo'))
full_results$Egg = factor(full_results$Egg,levels=egglevs$Egg)

# Each egg / tree has 100 bootstrapps, so add an identifier for each
fr <- full_results %>% 
  group_by(Egg, Tree) %>% 
  mutate(rep = row_number()) 

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut <- fr %>% 
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = mt_tree - a_tree)

# Calculate summary stats incl. 95% CIs
confs <- mt_vs_aut %>% group_by(Egg) %>% sum_stats(shifts)

plot_mtaut <- mt_vs_aut %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.15))+
  ggdist::stat_halfeye(adjust = 1,width = .3,.width = 0,justification = -.2, point_colour = NA,alpha = 0.95,normalize='groups')+
  #geom_boxplot(width = .15,outlier.shape = NA, alpha = 0.3,position=position_nudge(x=0.15)) +
  #gghalves::geom_half_point(aes(col=Egg),side='l',range_scale = .4,alpha = .3)+
  scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$col,breaks=egglevs$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut

png('../../figures/20240717_mtDNA-vs-Autosomal-shifts-MCMC.png',units='in',res=300,height=4,width=4)
pdf('../../figures/20240717_mtDNA-vs-Autosomal-shifts-MCMC.pdf',height=4,width=4)
plot_mtaut
dev.off()

#### Plot tree discordance ####
# This will also estimate transitions from each egg type to other eggs,  but doesn't seem as reliable as binary above approach

full_tree_results <- list()
full_egg_results <- list()

for (tree in c('mt_tree','a_tree')) {
  
  cat('Running full ML reconstruction on tree: ',tree,'\n')
  
  # Generate base tree 
  targ_tree <- get(tree)
  ggt <- ggtree(targ_tree, layout = "circular",branch.length='none') %<+% md 
  
  #grab only egg
  phenos <- as.data.frame(ggt$data %>% filter(isTip == TRUE))
  egg_mat <- as.matrix(phenos %>% select(Egg))
  phenotypes <- setNames(as.vector(egg_mat[, 1]), targ_tree$tip.label)
  
  #Plot probabilities 
  t2 <- multi2di(targ_tree)
  
  # Quick ape method 
  fitER <- ape::ace(phenotypes,t2,model="ER",type="discrete")
  
  # Determine ancestral state likelihood number of shifts using MCMC
  simfull <- make.simmap(t2, phenotypes, model="ER", nsim=100,Q='mcmc')
  
  # Output stats on the number of shifts across the 100 replicates 
  vec <- data.frame(shifts = countSimmap(simfull)$Tr[,1])
  shifts <- vec %>% mutate(Tree = tree, MLloglik = fitER$loglik, MLest = fitER$rates, MLse = fitER$se) 
  
  # Extract number straight from the MCMC approach
  mat <- countSimmap(simfull)$Tr
  col_names <- colnames(mat)
  
  # Calculate average transitions for each Egg
  calculate_statistics <- function(mat, prefix) {
    # Identify columns that contain the specific egg prefix
    relevant_cols <- grep(paste0("^", prefix, ",E|E,", prefix, "$"), col_names, value = TRUE)
    
    # Sum the values in these columns
    total_transitions <- rowSums(mat[, relevant_cols])
    
    # Calculate statistics
    stats <- list(
      Mean = mean(total_transitions),
      Median = median(total_transitions),
      SD = sd(total_transitions)
    )
    
    return(stats)
  }
  
  # Apply to each E group
  E_groups <- paste("E", 1:11, sep="")  # Adjust based on your specific E groups
  averages <- sapply(E_groups, function(E) calculate_statistics(mat, E))
  
  egg_shifts <- as.data.frame(t(averages)) %>% mutate(Egg = rownames(.), Tree = tree)
  rownames(egg_shifts) <- NULL
  
  full_tree_results[[tree]] <- shifts
  full_egg_results[[tree]] <- egg_shifts
  
  # Extract nodes and the proportions for pies
  nodes <- data.frame(
    node=1:t2$Nnode+Ntip(t2),
    fitER$lik.anc)
  
  ## cols parameter indicate which columns store stats
  pies <- nodepie(nodes, cols=2:12,outline.color='black',outline.size = 0.1)
  pies <- lapply(pies, function(g) g+scale_fill_manual(values = egglevs$col,breaks=egglevs$Egg))
  
  t3 <- full_join(t2, data.frame(label = names(phenotypes), stat = phenotypes ), by = 'label')
  tp <- ggtree(t3,layout='rectangular',branch.length = 'none') %<+% md
  tp$data$dummy <- 1
  tp_final <- tp + geom_inset(pies, width = .09, height = .09)+ 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg) 
  assign(paste0(tree,'_pies'),tp_final)
  
}

# Grab the model results
full_search_res <- rbindlist(full_tree_results)
full_search_res %>% group_by(Tree,MLloglik,MLest,MLse) %>% sum_stats(shifts)

# Also bind the egg data 
full_search_eggs <- rbindlist(full_egg_results)
full_search_eggs
full_search_eggs$Mean <- unlist(full_search_eggs$Mean)
full_search_eggs$Median <- unlist(full_search_eggs$Median)
full_search_eggs$SD <- unlist(full_search_eggs$SD)

write.table(full_search_res,file='Full_Search_Results_20240806.txt',quote=F,sep='\t',row.names=F)
write.table(full_search_eggs,file='Full_Search_EggResults_20240806.txt',quote=F,sep='\t',row.names=F)

full_search_res <- read_tsv('Full_Search_Results_20240806.txt')

# Plot all connections 
discord <- simple.tanglegram(tree1=mt_tree_pies,tree2=a_tree_pies,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('../../figures/20240717_TreeCompare-mtDNA-Auto-NodePiesML.pdf',height=5,width=7)
discord
dev.off()

# Swap 
discord2 <- simple.tanglegram(tree1=a_tree_pies,tree2=mt_tree_pies,column=dummy,value=1,t2_pad=2,l_color='black',tiplab=F)

pdf('../../figures/20240717_TreeCompare-Auto-mtDNA-NodePiesML.pdf',height=5,width=7)
discord2
dev.off()

# And as before, add the shift data 
frf <- full_search_res %>% select(Tree,shifts) %>% 
  group_by(Tree) %>% 
  mutate(rep = row_number()) 

# And for each bootstrap compare mtDNA vs autosomal shifts 
mt_vs_aut_full <- frf %>%
  pivot_wider(names_from = Tree,values_from = shifts) %>% 
  mutate(shifts = mt_tree - a_tree)

# Calculate summary stats incl. 95% CIs
confs_full <- mt_vs_aut_full %>% sum_stats(shifts)

## For binding with the overall reconstruction! ~~ load in data above FIRST ~~
mt_vs_autboth <- rbind(mt_vs_aut,mt_vs_aut_full %>% mutate(Egg = 'Multiclass'))
confs_all <- rbind(confs,confs_full %>% mutate(Egg = 'Multiclass'))
eggboth <- rbind(egglevs,data.frame(Egg='Multiclass',ord=12,col='grey80')) %>% arrange(desc(ord))
eggboth$Egg <- factor(eggboth$Egg,levels=eggboth$Egg)
mt_vs_autboth$Egg <- factor(mt_vs_autboth$Egg,levels=eggboth$Egg)
confs_all$Egg <- factor(confs_all$Egg,levels=eggboth$Egg)

# Plot
plot_mtaut <- mt_vs_autboth %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_errorbar(data=confs_all,width=0.2,col='black',aes(x=Egg,ymin=conf_low,ymax=conf_high),inherit.aes=FALSE,position=position_nudge(x=-0.15))+
  geom_boxplot(width = .15,outlier.shape = NA, alpha = 0.9,lwd=0.25,position=position_nudge(x=0.15)) +
  #ggdist::stat_halfeye(width = .3,.width = 0,justification = -.2, point_colour = NA,alpha = 0.95,normalize='groups')+
  scale_fill_manual(values=eggboth$col,breaks=eggboth$Egg)+
  scale_color_manual(values=eggboth$col,breaks=eggboth$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip(ylim=c(-15,15))
plot_mtaut

png('../../figures/20240807_mtDNA-vs-Autosomal-shifts-MCMC-Multiclass.png',units='in',res=300,height=4,width=4)
pdf('../../figures/20240807_mtDNA-vs-Autosomal-shifts-MCMC-Multiclass-Boxes.pdf',height=4.5,width=4)
plot_mtaut
dev.off()


# Plot full search egg shifts
eggdf <- as.data.frame(full_search_eggs)
mt_vs_aut <- eggdf %>% 
  select(Mean,Tree,Egg) %>% 
  pivot_wider(names_from = Tree,values_from = Mean) %>% 
  mutate(mt_tree = map_dbl(mt_tree, unlist),
         a_tree = map_dbl(a_tree, unlist),
         shifts = mt_tree - a_tree)
mt_vs_aut$Egg = factor(mt_vs_aut$Egg,levels=egglevs$Egg)


plot_mtaut_alt <- mt_vs_aut %>% ungroup %>% 
  ggplot(aes(x=Egg,y=shifts,fill=Egg))+
  geom_hline(yintercept=0,lty=2)+
  geom_hline(yintercept=c(-1,1),lty=3,col='darkgray')+
  geom_point(pch=21,size=3)+
  scale_fill_manual(values=egglevs$col,breaks=egglevs$Egg)+
  scale_color_manual(values=egglevs$col,breaks=egglevs$Egg)+
  theme_bw()+
  ylab('')+xlab('')+
  coord_flip()
plot_mtaut_alt

#inspect tree
for (egg in unique(egglevs$Egg)) { 
  
  cat('Making individual plot for egg: ',egg,'\n')
  md_inspect <- md %>% mutate(Egg = ifelse(Egg == egg,'1','0'))
  t1 <- ggtree(mt_tree, layout = "rectangular",branch.length='none') %<+% md_inspect + 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=rainbow(3))
  t2 <- ggtree(a_tree, layout = "rectangular",branch.length='none') %<+% md_inspect + 
    geom_tippoint(aes(fill=Egg),pch=21,size=1)+
    scale_fill_manual(values=rainbow(3))
  ts <- simple.tanglegram(tree1=t1,tree2=t2,column=Egg,value=1,t2_pad=2,l_color='black',tiplab=F) + ggtitle(egg)
  assign(paste0('p',egg),ts)
  
}

full_plots <- ggarrange(pE1,pE2,pE3,pE4,pE5,pE6,pE7,pE8,pE9,pE10,pE11,common.legend = TRUE,nrow=4,ncol=3)

pdf('../../figures/20240722_mtDNA-vs-Auto_Individual_Eggs.pdf',height=7,width=5)
full_plots
dev.off()
