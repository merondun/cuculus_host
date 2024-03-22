#### Plot Distance ~ FST / PHIst
setwd('~/merondun/cuculus_host/correlations_geography_genetics/')
.libPaths('~/mambaforge/envs/r/lib/R/library')
library(tidyverse)
library(viridis)
library(geosphere)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(meRo)
library(RColorBrewer)

#Read in metadata
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

d = read_tsv('Pairwise_GeographicDistance_Km_2024FEB28.txt')
f = read_tsv('Pairwise_FST_DistanceGroups_Autosomes_2024FEB28.txt')
p = read_tsv('Pairwise_PHIst_DistanceGroups_W-MT_2024FEB28.txt')

#correlation between W/MT
p %>% ggplot(aes(x=PHI_W,y=PHI))+geom_point()+theme_bw()
cor.test(p$PHI_W,p$PHI)
# 
# Pearson's product-moment correlation
# 
# data:  p$PHI_W and p$PHI
# t = 8.8912, df = 43, p-value = 2.695e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.6693965 0.8884474
# sample estimates:
#       cor 
# 0.8047956 

#prep and merge frames 
p = p %>% dplyr::rename(chr_W = PHI_W, chr_MT = PHI) %>% pivot_longer(!Group) %>% 
  dplyr::rename(FST = value, chr = name) %>% mutate(FST = pmax(0, pmin(1, FST)))
names(f) = c('chr','start','end','snps','FST','Group')
z_comparison = f %>% filter(chr == 'chr_Z') %>% mutate(Subset = ifelse(grepl('__M$',Group),'Males_Only','All_Samples'),Group = gsub('__M$','',Group)) %>% 
  select(-snps) %>% pivot_wider(names_from = Subset,values_from=FST) %>% na.omit 
z_comparison %>% ggplot(aes(x=All_Samples,y=Males_Only))+geom_point()+theme_bw()
cor.test(z_comparison$All_Samples,z_comparison$Males_Only)
# Pearson's product-moment correlation
# 
# data:  z_comparison$All_Samples and z_comparison$Males_Only
# t = 356.72, df = 34922, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8835307 0.8880480
# sample estimates:
#       cor 
# 0.8858103 

zc = f %>% filter(grepl('__M$',Group)) #only grab males for Z 
amc = f %>% filter(chr != 'chr_Z') #exclude the male only comparison for autosomes 
f2 = rbind(zc,amc) %>% mutate(FST = pmax(0, pmin(1, FST))) %>% select(!c(start,end,snps)) %>% mutate(Group = gsub('__M','',Group))
pf = rbind(p,f2)

#merge with geographic distance 
df = left_join(pf,d) 

#assign AvZ 
df = df %>% mutate(AvZ = ifelse(chr == 'chr_Z','Z',ifelse(chr == 'chr_MT','mtDNA',ifelse(chr == 'chr_W','W','Autosome'))))

dfs = df %>% group_by(AvZ,Group,P1,P2,Distance_km) %>% 
  sum_stats(FST)
dfs %>% ggplot(aes(x=log(Distance_km),y=log(mean),col=AvZ))+
  geom_point()+
  geom_smooth()+
  theme_bw()

# Spearman's rho and cor test within each AvZ level and gather coefficients
cors = dfs %>% group_by(AvZ) %>%
  summarize(
    rho = cor(mean, Distance_km,method='spearman'),
    p_value = cor.test(mean, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))

cols = brewer.pal(4,'Dark2')
dfs$AvZ = factor(dfs$AvZ,levels=c('Autosome','Z','W','mtDNA'))

#lmm, account for P1 and P2 
model_summaries <- dfs %>%
  group_by(AvZ) %>%
  do({
    model <- lmer(log(pmax(mean, 0.005)) ~ log(pmax(Distance_km, 0.005)) + (1|P1) + (1|P2), data = .)
    summary_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  }) %>%
  dplyr::bind_rows()
#output fixed effect of distance, bonferroni correction 
model_summaries %>% filter(grepl('Distance',term)) %>% select(-effect,-term) %>% ungroup %>% 
  mutate(padj = p.adjust(p.value,method='bonferroni'),
         signif = ifelse(padj < 0.05,'*','n.s.'))
# A tibble: 4 × 10
# AvZ      estimate std.error statistic    df  p.value conf.low conf.high     padj signif
# <fct>       <dbl>     <dbl>     <dbl> <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr> 
#   1 Autosome    1.22     0.0614     19.9   30.3 6.01e-19   1.10       1.35  2.41e-18 *     
#   2 Z           1.09     0.0747     14.6   34.1 3.07e-16   0.941      1.24  1.23e-15 *     
#   3 W           0.614    0.331       1.86  41.4 7.05e- 2  -0.0538     1.28  2.82e- 1 n.s.  
# 4 mtDNA       0.499    0.183       2.73  42.3 9.30e- 3   0.130      0.868 3.72e- 2 *   
  

# Create the plot
pp1 = ggplot(dfs, aes(x = log(pmax(Distance_km,0.005)), y = log(pmax(mean,0.005)), color = AvZ,shape=AvZ)) +
  geom_point(size=1) + #0.25 for main plot 
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +  #0.5 for main plot 
  # geom_text(data = cors[1,], aes(x = Inf, y = -Inf, label = label), vjust = -6, hjust = 1.73, col = cols[1]) +
  # geom_text(data = cors[2,], aes(x = Inf, y = -Inf, label = label), vjust = -4.5, hjust = 1.73, col = cols[2]) +
  # geom_text(data = cors[3,], aes(x = Inf, y = -Inf, label = label), vjust = -3, hjust = 1.1, col = cols[3]) +
  # geom_text(data = cors[4,], aes(x = Inf, y = -Inf, label = label), vjust = -1.5, hjust = 2.02, col = cols[4]) +
  labs(title = "Relationship between Mean and Distance_km",
       x = "log(Geographic Distance)",
       y = "log(FST | ΦST)",
       color = "AvZ") +
  scale_color_brewer(palette='Dark2')+
  facet_wrap(.~AvZ,scales='free',nrow=4,ncol=1)+ #for dual plot sensitivity
  scale_shape_manual(values=c(15,16,17,3))+
  theme_bw(base_size=6) + theme(legend.position = 'none')
pp1

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations_logFST-logDIST_2024FEB29.pdf',height=1.5,width=1.5)
pp1
dev.off()

# Spearman's rho and cor test within each AvZ level and gather coefficients
dfs = dfs %>% mutate(onephi = mean / (1 - mean))
cors2 = dfs %>% group_by(AvZ) %>%
  summarize(
    rho = cor(onephi, Distance_km,method='spearman'),
    p_value = cor.test(onephi, Distance_km,method='spearman')$p.value) %>% 
  mutate(padj = p.adjust(p_value,method='bonferroni'),
         signif = case_when(
           padj < 0.05 ~ "*",
           TRUE ~ " (n.s.)"),
         label = paste0('rs: ',round(rho, 2), signif))

pp2 = ggplot(dfs, aes(x = Distance_km, y = pmin(onephi,1), color = AvZ,shape=AvZ)) +
  geom_point(size=1) +
  geom_smooth(method = "lm", se = TRUE,lwd=0.75) +
  # geom_text(data = cors[1,], aes(x = Inf, y = -Inf, label = label), vjust = -6, hjust = 1.73, col = cols[1]) + 
  # geom_text(data = cors[2,], aes(x = Inf, y = -Inf, label = label), vjust = -4.5, hjust = 1.73, col = cols[2]) + 
  # geom_text(data = cors[3,], aes(x = Inf, y = -Inf, label = label), vjust = -3, hjust = 1.1, col = cols[3]) + 
  # geom_text(data = cors[4,], aes(x = Inf, y = -Inf, label = label), vjust = -1.5, hjust = 2.02, col = cols[4]) + 
  labs(title = "Relationship between FST | ΦST and Distance",
       x = "Geographic Distance (km)",
       y = "FST | ΦST / 1 - FST | ΦST)",
       color = "AvZ") +
  scale_color_manual('Compartment',values=brewer.pal(4,'Dark2'))+
  scale_shape_manual('Compartment',values=c(16,17,15,18))+
  theme_bw(base_size=6)+
  facet_wrap(.~AvZ,scales='free',nrow=4,ncol=1)+
  theme(legend.position='top')
pp2
pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations_1-logFST-logDIST_2024FEB29.pdf',height=1.5,width=1.5)
pp
dev.off()

write_tsv(cors,file='~/merondun/cuculus_host/correlations_geography_genetics/Correlations_logFST-logDIST_cortest-results.txt')

pdf('~/merondun/cuculus_host/correlations_geography_genetics/Correlations-BothMethods_2024FEB29.pdf',height=6,width=4)
ggarrange(pp1,pp2,common.legend = TRUE,nrow=1,ncol=2)
dev.off()
