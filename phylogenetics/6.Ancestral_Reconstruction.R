#### Plot Trees, Assign Haplogroups 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
.libPaths('~/mambaforge/envs/R/lib/R/library')
library(magick)
library(ggtree)
library(treeio)
library(dplyr)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(phytools)

#import metadata 
md = read.table('~/merondun/cuculus_host/Metadata_Host.txt',sep='\t',comment.char = '',header=TRUE)

###### CIRCULAR mtDNA TREE  #######
#plot tree 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_MT_MNLR.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))
gg = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md
plotTree(iqtr,type='fan',ftype='i')

#grab only egg
phenos = as.data.frame(gg$data %>% filter(isTip == TRUE))
egglevs = phenos %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord) %>% mutate(col = viridis(10, option='turbo'))
plotted_cols = egglevs %>% arrange(Egg)
phenos$Egg = factor(phenos$Egg,levels=egglevs$E)
pal = setNames(plotted_cols$col, plotted_cols$Egg)
palette(pal)

egg = as.matrix(phenos %>% select(Egg))
phenotypes <- setNames(as.vector(egg[, 1]), iqtr$tip.label)

#Plot probabilities 
iqtr_binary <- multi2di(iqtr)
fitER<-ace(phenotypes,iqtr_binary,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)

#plot probabilities
pdf('figures/Ancestral_Probabilities_2024MAR21.pdf',height=10,width=8)
plotTree(iqtr_binary,type="phylogram",fsize=0.8,ftype="i")
nodelabels(node=1:iqtr_binary$Nnode+Ntip(iqtr_binary),
           pie=fitER$lik.anc,piecol=plotted_cols$col,cex=0.5)
#tiplabels(pie=to.matrix(phenotypes,sort(unique(phenotypes))),piecol=plotted_cols$col,cex=0.15)
legend(legend=egglevs$Egg,pch=22,fill=egglevs$col,"bottomleft")
dev.off()

#mcmc approach for marginal states 
tree_simmap <- make.simmap(iqtr_binary, phenotypes, model="ER", nsim=100,Q='mcmc')
pd<-summary(tree_simmap,plot=FALSE)

pdf('figures/Ancestral_MarginalStates_2024MAR21.pdf',height=10,width=8)
#plot(pd,type="phylogram",fsize=0.6,ftype="i")

#without tips
plotd = pd
plotd$tips = NA #set the tips to NA 
plot(plotd,fsize=0.6,ftype="i")

legend(legend=egglevs$Egg,pch=22,fill=egglevs$col,"topleft")
dev.off()

pdf('figures/Ancestral_Probs_vs_MarginalStates_2024MAR21.pdf',height=4.5,width=4)
plot(fitER$lik.anc,pd$ace[1:151,],xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)
dev.off()

#and for W chromosome 
iqtree = read.iqtree('beast_dating/variant_only/ml_trees/chr_W_MNLR.SNP.DP3-AC2-MQ40.min4.phy.contree')
iqtr = midpoint.root(as.phylo(iqtree))
gg = ggtree(iqtr, layout = "circular",branch.length='none') %<+% md
plotTree(iqtr,type='fan',ftype='i')

#grab only egg
phenos = as.data.frame(gg$data %>% filter(isTip == TRUE))
#use the mtDNA color scheme from above 
w_cols = plotted_cols %>% filter(Egg != 'E4') %>% arrange(ord)
w_egg = w_cols %>% arrange(Egg)
w_pal = setNames(w_egg$col, w_egg$Egg)
palette(w_pal)

egg = as.matrix(phenos %>% select(Egg))
phenotypes <- setNames(as.vector(egg[, 1]), iqtr$tip.label)

#Plot probabilities 
iqtr_binary <- multi2di(iqtr)
fitER<-ace(phenotypes,iqtr_binary,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)

#plot probabilities
pdf('figures/Ancestral_Probabilities_WChromosome_2024MAR21.pdf',height=10,width=8)
plotTree(iqtr_binary,type="phylogram",fsize=0.8,ftype="i")
nodelabels(node=1:iqtr_binary$Nnode+Ntip(iqtr_binary),
           pie=fitER$lik.anc,piecol=w_egg$col,cex=0.5)
legend(legend=w_cols$Egg,pch=22,fill=w_cols$col,"bottomleft")
dev.off()

#mcmc approach for marginal states 
tree_simmap <- make.simmap(iqtr_binary, phenotypes, model="ER", nsim=100)
pd<-summary(tree_simmap,plot=FALSE)

pdf('figures/Ancestral_MarginalStates_WChromosome_2024MAR21.pdf',height=10,width=8)
plotd = pd
plotd$tips = NA #set the tips to NA 
plot(plotd,fsize=0.6,ftype="i")
legend(legend=w_cols$Egg,pch=22,fill=w_cols$col,"topleft")
dev.off()

pdf('figures/Ancestral_Probs_vs_MarginalStates_WChromosome_2024MAR21.pdf',height=4.5,width=4)
plot(fitER$lik.anc,pd$ace[1:87,],xlab="marginal ancestral states",
     ylab="posterior probabilities from stochastic mapping")
lines(c(0,1),c(0,1),lty="dashed",col="red",lwd=2)
dev.off()

