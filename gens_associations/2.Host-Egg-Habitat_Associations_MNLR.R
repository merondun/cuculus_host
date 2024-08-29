# Egg, Host, Habitat associations, MNLR, Model Selection
#### Find associations between MT/W haplotypes and features 
setwd('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host')
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

set.seed(123)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(Analysis_PopulationGenetics == 1) %>%
  drop_na(Egg) %>%
  select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist)

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

#only retain Host where we have at least 2 cuckoos 
md_egg = md_egg %>% group_by(Host) %>% mutate(TotalHost = n()) %>% ungroup %>% group_by(Habitat) %>% mutate(TotalHabitat = n()) %>% ungroup %>% group_by(Egg) %>% mutate(TotalEgg = n())  %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalHost >= minobs & TotalHabitat >= minobs & TotalEgg >= minobs) 
md_egg %>% filter(!ID %in% md$ID)
#write.table(md$ID,file='randomforest/Samples_Retained_Unrelated_MNLR_2024MAR14.txt',quote=F,sep='\t',row.names=F,col.names=F)

# If you want to exclude the blue clades W1, W2, W3! 
md_egg <-  md_egg %>% filter(!grepl('W1|W2|W3',Haplogroup))

#### Count proportions first, count proportions for host and habitat and egg  
hp = md_egg %>% group_by(Host) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Host,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
tp = md_egg %>% group_by(Habitat) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Habitat,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
ep = md_egg %>% group_by(Egg) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Egg,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
egglev = ep %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord)
ep$Egg = factor(ep$Egg,levels=egglev$Egg)

# Bind them together 
ap = rbind(hp %>% ungroup %>% mutate(Response = Host, variable = 'Host') %>% select(-Host), 
           tp %>% ungroup %>% mutate(Response = Habitat, variable = 'Habitat') %>% select(-Habitat),
           ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

# Just for ordering the covariates nicely 
ord = ap %>% ungroup %>% select(name,value) %>% unique %>% mutate(ord = as.numeric(str_sub(value,2))) %>% group_by(name) %>% arrange(ord)
ap$value = factor(ap$value,levels=ord$value)

# Plot proportions
app = ap %>% 
  ggplot(aes(y=Response,x=value,fill=Proportion,label=Parts))+
  facet_grid(variable~name,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(legend.position = 'top')
app

#pdf('figures/Proportions_HostHabitat_2024MAR14.pdf',height=6.5,width=6)
app
#dev.off()

##### Model Selection #####
#assess covariate importance with model selection, using MNLR
vars = c('Habitat','Host','Egg')
vars = 'Egg'

# With bootstrapping
ctrl = trainControl(method = "boot",   # Use bootstrapping
                    number = 100,      # Replicates 
                    summaryFunction = multiClassSummary,  #caret now attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

# Change punctuation e.g. 'A. pal' to A_pal'
md_cv = md_egg %>% mutate(Host = gsub('\\. ','_',Host))

# For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0
for (minnum in c(2)) {  # Loop through 2, 3, or 4 minimum observations per category
  
  # Filter at the very base level, ensuring that across egg / host / habitat we have the same individuals with representative sampling
  md_subbed = md_cv %>% filter(TotalHost >= minnum & TotalHabitat >= minnum & TotalEgg >= minnum) %>%
    mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)
  
  for (rep in seq(1,10,1)){  # Create 10 replicate models
    for (var in vars) { counter = counter + 1;
    
    # Ensure that we have adequate levels, only a sanity since it is already filtered
    retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% pull(var)
    length(retained)
    # Number of samples to subsample
    subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(min = min(n)) %>% pull(min)
    
    cat('Downsampling to n = ',subsamp,', requiring min = ',minnum,' for variable: ',var,', ','replicate: ',rep,'\n')
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
      dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = minnum,
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
        mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = minnum, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
      new_preds = rbind(new_preds,conf_real)
      rm(conf_real,dat,dat_best)
      
    } # Exit model loop
    } # Exit response variable loop
  } # Exit iteration loop
} # Exit minimum samples loop

# Write the no-blue data
write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_NoW1W2W3_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)
write.table(new_preds,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_NoW1W2W3_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)

write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)
write.table(new_preds,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_Boot-2Obs_2024AUG06.txt',quote=F,sep='\t',row.names=F)

# Or start here and read in saved data 
adat = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Boot-2Obs_2024AUG06.txt')
conf_mats = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix_Boot-2Obs_2024AUG06.txt')

# Plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('W+K+B','W+B','W','W+K','K','K+B','B'))
model_dat = adat %>% left_join(.,leg)

# Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('W+K+B','W+B','W+K','K+B','B','W','K'))
app = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
app 

pdf('figures/2024AUG06_MNLR_ModelSelection_AUC.pdf',height=4,width=6)
app
dev.off()

cols <- brewer.pal(3,'Set2')[c(1,2,3)]
model_dat$Label <- factor(model_dat$Label,levels=c('K','B','W'))
model_dat$Variable <- factor(model_dat$Variable,levels=c('Host','Habitat','Egg'))
auc_plot_input <- model_dat %>%
  filter(MinObs == 2) %>% 
  filter(Label == 'W' | Label == 'K' | Label == 'B') %>% 
  group_by(Label,Variable) %>% 
  sum_stats(AUC)

auc_plot <- auc_plot_input %>% 
  ggplot(aes(y=Variable,fill=Label,x=mean,xmin=conf_low,xmax=conf_high))+
  geom_bar(stat='identity',position=position_dodge(width=0.9))+
  scale_fill_manual(values=cols)+
  xlab('')+ylab('AUC (Mean & 95% CI)')+
  geom_errorbar(width=0.25,position=position_dodge(width=0.9))+
  theme_bw(base_size=6)+
  coord_cartesian(xlim=c(0.7,1))+
  theme(legend.position='top')
auc_plot

pdf('figures/2024AUG06_MNLR_ModelSelection_AUC.pdf',height=2,width=1.5)
auc_plot
dev.off()

write.table(auc_plot_input,file='figures/20240806_AUC_Results_Boot-2Obs.txt',quote=F,sep='\t',row.names=F)

#order full, single plot, make sure the 3 variables are in order 
egglevs = conf_mats %>% filter(Variable == 'Egg') %>% select(Prediction,Variable) %>% mutate(ord = as.numeric(gsub('E','',Prediction)),Prediction) %>% unique %>% arrange(Variable,ord) %>% select(-ord)
hostlevs = conf_mats %>% filter(Variable == 'Host') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
hablevs = conf_mats %>% filter(Variable == 'Habitat') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
all_levs = rbind(hostlevs,egglevs,hablevs)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=all_levs$Prediction),
                                 Reference = factor(Reference,levels=all_levs$Prediction))

### Plot how the addition of haplogroup improves predictions show K+B (m6) vs W+K+B (m1)
# ancestry + geogrpahy
repredictions_ancgeo = conf_mats %>% 
  filter(Model == 'm6') %>%  #main figure, plot model 2 (W+B) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)

#plot 
ancgeo_plot = repredictions_ancgeo %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  #geom_point(size=0.2)+
  facet_wrap(Variable~.,scales='free')+
  scale_color_continuous(low='white',high='black')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeo_plot

# ancestry + geogrpahy + hap 
repredictions_ancgeohap = conf_mats %>% 
  filter(Model == 'm1') %>%  #main figure, plot model 2 (W+B) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)

#plot 
ancgeohap_plot = repredictions_ancgeohap %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  #geom_point(size=0.2)+
  scale_color_continuous(low='white',high='black')+
  facet_wrap(Variable~.,scales='free')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
ancgeohap_plot

pdf('figures/2024AUG08_MNLR_ConfusionMatrix-Repredictions_M1vsM6.pdf',height=4.5,width=5)
ggarrange(ancgeo_plot,ancgeohap_plot,nrow = 2,common.legend = TRUE)
dev.off()

