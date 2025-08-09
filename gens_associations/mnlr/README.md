# Multinomial logistic regression

Create a MNLR with egg as the response variable and Geography + Autosomal K + Haplogroups as the predictors. In short:

* Only retain response variables where there are at least n=2 observations
* Downsample all response classes so that all classes have n=2 observations
* Fit 7 multinomial logistic regression models, each with n=100 bootstraps using all combinations of predictors
* Extract AUC, and use the model to predict response variable on the full dataset again (too small for unseen data prediction)
* Repeat the above procedure 100 times so that different downsampled observations are included 
* Determine which classes are predicted correctly (% correct) from the confusion matrix on real / predicted responses across bootstraps

Run model:

```R
#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

# Egg associations, MNLR, Model Selection
#### Find associations between MT/W haplotypes and features
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

sp = args[1]

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

#### Count proportions first, count proportions for Egg 
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>%
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>%
  group_by(Egg,name,value) %>%
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

# Bind them together
ap = rbind(ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely
# Plot proportions
ord <- ap %>% select(name, value) %>%
  distinct() %>%
  mutate(ord = as.numeric(gsub("[^0-9.]", "", value))) %>%
  arrange(name, ord)
ap$value <- factor(ap$value,levels=ord$value)
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
ap$Response <- factor(ap$Response,levels=egglev$Egg)
app = ap %>%
  arrange(value) %>%
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_Proportions.pdf'),app,height=3,width=7,dpi=300)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR
vars = 'Egg'

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# Change punctuation e.g. 'A. pal' to A_pal' if necessary
md_cv = md_egg %>% mutate(Egg = gsub('\\. ','_',Egg))

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
md_subbed = md_cv %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

for (rep in seq(1,100,1)){  # Create 100 replicate models
  for (var in vars) { counter = counter + 1;

  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)

  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))

  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)

  # First MNLR on combinations
  formula_1 = as.formula(paste(var, "~ Haplogroup + Ancestry + Geography"))
  m1 = train(formula_1, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_2 = as.formula(paste(var, "~ Haplogroup + Geography"))
  m2 = train(formula_2, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_3 = as.formula(paste(var, "~ Haplogroup "))
  m3 = train(formula_3, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_4 = as.formula(paste(var, "~ Haplogroup + Ancestry"))
  m4 = train(formula_4, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_5 = as.formula(paste(var, "~ Ancestry"))
  m5 = train(formula_5, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_6 = as.formula(paste(var, "~ Ancestry + Geography"))
  m6 = train(formula_6, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  formula_7 = as.formula(paste(var, "~ Geography"))
  m7 = train(formula_7, data = mdi, method = "multinom", trControl = ctrl, metric = "AUC", trace = FALSE)

  models = c('m1','m2','m3','m4','m5','m6','m7')

  # Extract model fit
  for (model in models) {
    # Output model fit from confusion matrix
    mo = get(model)

    # Get AIC
    final_model = mo$finalModel;
    AIC = AIC(final_model)

    # Save the model results
    dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat,dat_best)

    # Also save training confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)

    # Predict against real data
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)

    conf_new = confusionMatrix(new$predicted, new$reference)
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
    new_preds = rbind(new_preds,conf_real)
    rm(conf_real,dat,dat_best)

  } # Exit model loop
  } # Exit response variable loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250330_Model_Selection_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250330_ConfusionMatrix_Boot-2Obs-100Reps_',sp,'.txt'),quote=F,sep='\t',row.names=F)


```

### Plot

```bash
# Plot Egg associations, MNLR, Model Selection 
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

sp = 'CC'
set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_FullGensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

#### Count proportions first, count proportions for Egg
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>%
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>%
  group_by(Egg,name,value) %>%
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

# Bind them together
ap = rbind(ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely
# Plot proportions
ord <- ap %>% select(name, value) %>%
  distinct() %>%
  mutate(ord = as.numeric(gsub("[^0-9.]", "", value))) %>%
  arrange(name, ord)
ap$value <- factor(ap$value,levels=ord$value)
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
ap$Response <- factor(ap$Response,levels=egglev$Egg)
app = ap %>%
  arrange(value) %>%
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_Proportions.pdf'),app,height=3,width=7,dpi=300)

# Read in saved data 
adat = read_tsv(paste0('20250330_Model_Selection_Boot-2Obs-100Reps_',sp,'.txt'))
conf_mats = read_tsv(paste0('20250330_ConfusionMatrix_Boot-2Obs-100Reps_',sp,'.txt'))

# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('A+G+M','G+M','M','A+M','A','A+G','G'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('A+G+M','A+G','A+M','G+M','A','G','M'))
auc_plot = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
auc_plot

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection_AUC.pdf'),
       auc_plot,height=3,width=7,dpi=300)


# Summarize AUC across the core 3 models 
auc_plot_input <- model_dat %>%
  # %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

# Plot
cols <- brewer.pal(3,'Set2')[c(1,2,3)]
#model_dat$Label <- factor(model_dat$Label,levels=c('A','G','M'))

auc_summary_plot <- auc_plot_input %>% 
  filter(Label == 'A' | Label == 'G' | Label == 'M') %>% 
  ggplot(aes(y=Label,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.5,1))+
  theme(legend.position='top')
auc_summary_plot

ggsave(paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection95CI_AUC.pdf'),
       auc_summary_plot,height=2,width=1.5,dpi=300)

write.table(auc_plot_input,file=paste0('~/symlinks/host/figures/20250330_',sp,'_MNLR_ModelSelection_AUCResults.txt'),quote=F,sep='\t',row.names=F)

#order full, single plot, make sure the 3 variables are in order 
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=egglev$Egg),
                                 Reference = factor(Reference,levels=egglev$Egg))

### Plot how the addition of haplogroup improves predictions show A+G (m6) vs A+G+M (m1)
auc_vals <- adat %>% group_by(Model) %>% sum_stats(AUC)

### only geography (G; model 7) 
lab <- auc_vals %>% filter(Model == 'm7') %>% mutate(label = paste0('G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_geo = conf_mats %>% 
  filter(Model == 'm7') %>%  # (G ONLY) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
geo_plot = repredictions_geo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
geo_plot


### ancestry + geography (A+G; m6)
lab <- auc_vals %>% filter(Model == 'm6') %>% mutate(label = paste0('A+G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  # plot model 6 (A+G) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap (A+G+M, m1)
lab <- auc_vals %>% filter(Model == 'm1') %>% mutate(label = paste0('A+G+M: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  # plot model 1 (A+G+M) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
ggsave(paste0('~/symlinks/host/figures/20250330_MNLR_ConfusionMatrix-Repredictions-',sp,'_M1vsM6.pdf'),
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       #dpi=300,height=3,width=1) # optatus 
       dpi=300,height=3.5,width=1.5) # canorus

```

### Plot females-only

```R
# Egg associations, MNLR, Model Selection FEMALES 
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rf25/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 
sp='CC'
sp='CO'

# Read in saved data 
adat = read_tsv(paste0('20250729_Model_Selection_Boot-2Obs-FEMALES_',sp,'.txt'))
conf_mats = read_tsv(paste0('20250729_ConfusionMatrix_Boot-2Obs-FEMALES_',sp,'.txt'))

# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('A+G+M','G+M','M','A+M','A','A+G','G'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('A+G+M','A+G','A+M','G+M','A','G','M'))
auc_plot = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
auc_plot
model_dat

ggsave(paste0('~/symlinks/host/figures/20250729_',sp,'_MNLR_ModelSelection_AUC-FEMALES.pdf'),
       auc_plot,height=3,width=7,dpi=300)

#order full, single plot, make sure the 3 variables are in order 
egglev <- md %>% filter(Analysis_GensAssociations==1) %>%  select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=egglev$Egg),
                                 Reference = factor(Reference,levels=egglev$Egg))

### Plot how the addition of haplogroup improves predictions show A+G (m6) vs A+G+M (m1)
auc_vals <- adat %>% group_by(Model) %>% sum_stats(AUC)

### ancestry + geography (A+G; m6)
lab <- auc_vals %>% filter(Model == 'm6') %>% mutate(label = paste0('A+G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  # plot model 6 (A+G) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap (A+G+M, m1)
lab <- auc_vals %>% filter(Model == 'm1') %>% mutate(label = paste0('A+G+M: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  # plot model 1 (A+G+M) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
ggsave(paste0('~/symlinks/host/figures/20250729_MNLR_ConfusionMatrix-Repredictions-',sp,'_FEMALES.pdf'),
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       dpi=300,height=3,width=1.25) # optatus 
       #dpi=300,height=3.5,width=1.5) # canorus

```



### Sensitivity: Egg Exclusions

```bash
# Egg associations, MNLR, Model Selection: EXCLUDE common eggs and ancient haplogroup eggs 
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

sp = 'CC'

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# For one CC fork, also drop the most abundant egg types (E1 and E6) to see how the results are impacted
md_egg <- md_egg %>% filter(!Egg %in% c('ECC1','ECC6'))

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg we have the same individuals with representative sampling
md_subbed = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)
var = 'Egg'
for (rep in seq(1,100,1)){  # Create 10 replicate models
  counter = counter + 1;
  
  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)
  
  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))
  
  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)
  
  # Function to safely train models
  safe_train <- function(formula, data, method, trControl, metric) {
    tryCatch({
      train(formula, data = data, method = method, trControl = trControl, metric = metric, trace = FALSE)
    }, error = function(e) {
      message(paste("Model failed:", as.character(formula), "Error:", e$message))
      return(NULL)
    })
  }
  
  # Training models inside loop
  m1 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m2 = safe_train(as.formula(paste(var, "~ Haplogroup + Geography")), mdi, "multinom", ctrl, "AUC")
  m3 = safe_train(as.formula(paste(var, "~ Haplogroup")), mdi, "multinom", ctrl, "AUC")
  m4 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry")), mdi, "multinom", ctrl, "AUC")
  m5 = safe_train(as.formula(paste(var, "~ Ancestry")), mdi, "multinom", ctrl, "AUC")
  m6 = safe_train(as.formula(paste(var, "~ Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m7 = safe_train(as.formula(paste(var, "~ Geography")), mdi, "multinom", ctrl, "AUC")
  
  models = list(m1, m2, m3, m4, m5, m6, m7)
  model_names = c("m1", "m2", "m3", "m4", "m5", "m6", "m7")
  
  for (i in seq_along(models)) {
    mo = models[[i]]
    
    if (is.null(mo)) {
      message(paste("Skipping", model_names[i], "due to failure"))
      next
    }
    
    # Proceed with extracting metrics if model trained successfully
    final_model = mo$finalModel
    AIC = AIC(final_model)
    dat = data.frame(Model = model_names[i], Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, 
                     AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, 
                     AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat, dat_best)
    
    # Confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
    conf_new = confusionMatrix(new$predicted, new$reference)
    
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model_names[i], Iteration = counter, Variable = var, 
             Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,
             logloss=dat_best$logLoss, Accuracy = dat_best$Accuracy, 
             AccuracySD=dat_best$AccuracySD)
    
    new_preds = rbind(new_preds, conf_real)
    
  } # Exit model loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250320_Model_Selection_Boot-2Obs_',sp,'-NoE1E6.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250320_ConfusionMatrix_Boot-2Obs_',sp,'-NoE1E6.txt'),quote=F,sep='\t',row.names=F)





######## SECOND
set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# For one CC fork, also drop the most abundant egg types (E1 and E6) to see how the results are impacted
md_egg <- md_egg %>% filter(!Haplogroup %in% c('MCC1','MCC2','MCC3'))

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)

##### Model Selection #####
#assess covariate importance with model selection, using MNLR

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0

# Filter at the very base level, ensuring that across egg we have the same individuals with representative sampling
md_subbed = md_egg %>% mutate_at(c('Egg','Haplogroup','Ancestry','Geography'),as.factor)
var = 'Egg'
for (rep in seq(1,100,1)){  # Create 10 replicate models
  counter = counter + 1;
  
  # Ensure that we have adequate levels, only a sanity since it is already filtered
  retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% pull(var)
  length(retained)
  # Number of samples to subsample
  subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= 2) %>% summarize(min = min(n)) %>% pull(min)
  
  cat('Downsampling to n = ',subsamp,', requiring min = ',2,' for variable: ',var,', ','replicate: ',rep,'\n')
  # Subsampling
  mdi = md_subbed %>%
    filter(!!sym(var) %in% retained) %>%
    group_by(!!sym(var)) %>%
    sample_n(min(n(), subsamp),replace = TRUE) %>%
    ungroup()
  mdi %>% count(!!sym(var))
  
  # Ensure the factor levels are dropped which aren't in the dataset after filtering
  mdi = droplevels(mdi)
  
  # Function to safely train models
  safe_train <- function(formula, data, method, trControl, metric) {
    tryCatch({
      train(formula, data = data, method = method, trControl = trControl, metric = metric, trace = FALSE)
    }, error = function(e) {
      message(paste("Model failed:", as.character(formula), "Error:", e$message))
      return(NULL)
    })
  }
  
  # Training models inside loop
  m1 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m2 = safe_train(as.formula(paste(var, "~ Haplogroup + Geography")), mdi, "multinom", ctrl, "AUC")
  m3 = safe_train(as.formula(paste(var, "~ Haplogroup")), mdi, "multinom", ctrl, "AUC")
  m4 = safe_train(as.formula(paste(var, "~ Haplogroup + Ancestry")), mdi, "multinom", ctrl, "AUC")
  m5 = safe_train(as.formula(paste(var, "~ Ancestry")), mdi, "multinom", ctrl, "AUC")
  m6 = safe_train(as.formula(paste(var, "~ Ancestry + Geography")), mdi, "multinom", ctrl, "AUC")
  m7 = safe_train(as.formula(paste(var, "~ Geography")), mdi, "multinom", ctrl, "AUC")
  
  models = list(m1, m2, m3, m4, m5, m6, m7)
  model_names = c("m1", "m2", "m3", "m4", "m5", "m6", "m7")
  
  for (i in seq_along(models)) {
    mo = models[[i]]
    
    if (is.null(mo)) {
      message(paste("Skipping", model_names[i], "due to failure"))
      next
    }
    
    # Proceed with extracting metrics if model trained successfully
    final_model = mo$finalModel
    AIC = AIC(final_model)
    dat = data.frame(Model = model_names[i], Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = 2,
                     decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, 
                     AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, 
                     AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
    
    dat_best = dat %>% slice_min(logLoss)
    adat = rbind(adat, dat_best)
    
    # Confusion matrix
    pred = mo$pred
    pred_best = pred %>% filter(decay == dat_best$decay)
    predicted_classes = predict(mo, newdata = mdi, type = "raw")
    new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
    conf_new = confusionMatrix(new$predicted, new$reference)
    
    conf_real = as.data.frame(conf_new$table) %>% mutate(
      Prediction = gsub('_','\\. ',Prediction),
      Reference = gsub('_','\\. ',Reference)) %>%
      mutate(Model = model_names[i], Iteration = counter, Variable = var, 
             Subsample = subsamp, MinObs = 2, AUC=dat_best$AUC,
             logloss=dat_best$logLoss, Accuracy = dat_best$Accuracy, 
             AccuracySD=dat_best$AccuracySD)
    
    new_preds = rbind(new_preds, conf_real)
    
  } # Exit model loop
} # Exit iteration loop

# Write the no-blue data
write.table(adat,paste0('20250320_Model_Selection_Boot-2Obs_',sp,'-NoM1M2M3.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250320_ConfusionMatrix_Boot-2Obs_',sp,'-NoM1M2M3.txt'),quote=F,sep='\t',row.names=F)

```

### Plot Sensitivity

Plot the sensitivity analyses above on egg exclusions:

```R
# Egg associations, MNLR, Model Selection: PLOT SENSITIVITY
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
.libPaths('~/mambaforge/envs/rfs/lib/R/library')
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(vegan)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggdist)
library(nnet)
library(purrr)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

sp = 'CC'

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1) %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

# Read in data from sensitivities
# for No M1/M2/M3
adat = read_tsv('20250320_Model_Selection_Boot-2Obs_CC-NoM1M2M3.txt')
conf_mats = read_tsv('20250320_ConfusionMatrix_Boot-2Obs_CC-NoM1M2M3.txt')

# for No E1/E6
adat = read_tsv('20250320_Model_Selection_Boot-2Obs_CC-NoE1E6.txt')
conf_mats = read_tsv('20250320_ConfusionMatrix_Boot-2Obs_CC-NoE1E6.txt')


# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('A+G+M','G+M','M','A+M','A','A+G','G'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('A+G+M','A+G','A+M','G+M','A','G','M'))
auc_plot = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
auc_plot

# For no M1/M2/M3
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection_AUC-NoM1M2M3.pdf',
       auc_plot,height=3,width=7,dpi=300)

# For no E1/E6
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection_AUC-NoE1E6.pdf',
       auc_plot,height=3,width=7,dpi=300)


# Summarize AUC across the core 3 models 
auc_plot_input <- model_dat %>%
  # %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

# Plot
cols <- brewer.pal(3,'Set2')[c(1,2,3)]
#model_dat$Label <- factor(model_dat$Label,levels=c('A','G','M'))

auc_summary_plot <- auc_plot_input %>% 
  filter(Label == 'A' | Label == 'G' | Label == 'M') %>% 
  ggplot(aes(y=Label,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.5,1))+
  theme(legend.position='top')
auc_summary_plot

# for no M1/M2/M3
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection95CI_AUC-NoM1M2M3.pdf',
       auc_summary_plot,height=2,width=1.5,dpi=300)

# for no E1/E6
ggsave('~/symlinks/host/figures/20250326_CC_MNLR_ModelSelection95CI_AUC-NoE1E6.pdf',
       auc_summary_plot,height=2,width=1.5,dpi=300)

#order full, single plot, make sure the 3 variables are in order 
egglev <- md %>% select(Egg,EggOrder) %>% distinct %>% arrange(EggOrder)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=egglev$Egg),
                                 Reference = factor(Reference,levels=egglev$Egg))

### Plot how the addition of haplogroup improves predictions show A+G (m6) vs A+G+M (m1)
auc_vals <- adat %>% group_by(Model) %>% sum_stats(AUC)


### ancestry + geography (A+G; m6)
lab <- auc_vals %>% filter(Model == 'm6') %>% mutate(label = paste0('A+G: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  # plot model 6 (A+G) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap (A+G+M, m1)
lab <- auc_vals %>% filter(Model == 'm1') %>% mutate(label = paste0('A+G+M: ',round(conf_low,2),' - ',round(conf_high,2)))
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  # plot model 1 (A+G+M) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  ggtitle(lab$label)+
  facet_wrap(lab$label~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
# For no M1/M2/M3
ggsave('~/symlinks/host/figures/20250320_MNLR_ConfusionMatrix-Repredictions-CC_M1vsM6-NoM1M2M3.pdf',
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       dpi=300,height=3.5,width=1.5) # canorus

# For no E1/E6
ggsave('~/symlinks/host/figures/20250320_MNLR_ConfusionMatrix-Repredictions-CC_M1vsM6-NoE1E6.pdf',
       ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE),
       dpi=300,height=3.5,width=1.5) # canorus


```



