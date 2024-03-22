#Egg, Host, Habitat associations, MNLR, Model Selection
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
# This initial section assigns ancestry 'K' values for related individuals according to their relatives, since they weren't included in the popgen analysis.
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')

#only grab samples with known egg (hash out filter for related individuals)
md = md %>%
  filter(Analysis_PopulationGenetics == 1) %>%
  drop_na(Egg) %>%
  select(ID = ID, Host = HostParentShort, Egg, Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist)

# #add similar K for relatives which weren't included in pop gen analyses
# unknowns = md %>% filter(is.na(Ancestry))
# rels = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/relatedness/chr_10.CC.rel')
# rels = rels %>% select(IDA = INDV1, IDB = INDV2, phi = RELATEDNESS_PHI)
# relk = left_join(rels,md %>% select(IDA=ID,KA=Ancestry)) %>% left_join(.,md %>% select(IDB=ID,KB=Ancestry)) %>% filter(IDA != IDB)
# assigned = NULL
# for (unk in unique(unknowns$ID)) {
#   targ = relk %>% filter(IDA == unk | IDB == unk) %>% mutate(K = ifelse(is.na(KA),KB,KA)) %>% drop_na(K) %>% slice_max(phi) %>% select(K) %>% unique %>% mutate(ID = unk)
#   assigned = rbind(assigned,targ)
#   cat('Sample : ',unk,' is ancestry: ',targ$K,'\n')
# }
# md_egg = left_join(md,assigned) %>% mutate(Ancestry = ifelse(is.na(Ancestry),K,Ancestry)) %>% select(-K)
# write.table(md_egg,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/MNLR_Input_Egg_RelatedK_2024MAR16.txt',quote=F,sep='\t',row.names=F)

#md_egg = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/MNLR_Input_Egg_RelatedK_2024MAR16.txt')
md_egg = md 

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

#only retain Host where we have at least 2 cuckoos 
md_egg = md_egg %>% group_by(Host) %>% mutate(TotalHost = n()) %>% ungroup %>% group_by(Habitat) %>% mutate(TotalHabitat = n()) %>% ungroup %>% group_by(Egg) %>% mutate(TotalEgg = n())  %>% ungroup
minobs=2
md = md_egg %>% filter(TotalHost >= minobs & TotalHabitat >= minobs & TotalEgg >= minobs) 
md_egg %>% filter(!ID %in% md$ID)
#write.table(md$ID,file='randomforest/Samples_Retained_Unrelated_MNLR_2024MAR14.txt',quote=F,sep='\t',row.names=F,col.names=F)

#### Count proportions first, count proportions for host and habitat and egg  
hp = md %>% group_by(Host) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Host,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
tp = md %>% group_by(Habitat) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Habitat,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
ep = md %>% group_by(Egg) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Egg,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
egglev = ep %>% ungroup %>% select(Egg) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Egg))) %>% arrange(ord)
ep$Egg = factor(ep$Egg,levels=egglev$Egg)

#bind them together 
ap = rbind(hp %>% ungroup %>% mutate(Response = Host, variable = 'Host') %>% select(-Host), 
           tp %>% ungroup %>% mutate(Response = Habitat, variable = 'Habitat') %>% select(-Habitat),
           ep %>% ungroup %>% mutate(Response = Egg, variable = 'Egg') %>% select(-Egg))

#just for ordering the covariates nicely 
ord = ap %>% ungroup %>% select(name,value) %>% unique %>% mutate(ord = as.numeric(str_sub(value,2))) %>% group_by(name) %>% arrange(ord)
ap$value = factor(ap$value,levels=ord$value)

#plot 
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

#set up cross validation
ctrl = trainControl(method = "cv",   #Use cross-validation
                    number = 5,      #number of folds
                    summaryFunction = multiClassSummary,  #caret now attempts to stratify levels
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

#change 'A. pal' to A_pal'
md_cv = md %>% mutate(Host = gsub('\\. ','_',Host))

#For determining which predictors improve model:
adat = NULL; conf_dat = NULL; class_dat = NULL; new_preds = NULL; counter = 0
for (minnum in c(2)) {  #cycle through 2, 3, or 4 minimum observations per category
  
  #Filter at the very base level, ensuring that across egg / host / habitat we have the same individuals with representative sampling
  md_subbed = md_cv %>% filter(TotalHost >= minnum & TotalHabitat >= minnum & TotalEgg >= minnum) %>%
    mutate_at(c('Host','Egg','Habitat','Haplogroup','Ancestry','Geography'),as.factor)
  for (rep in seq(1,100,1)){  #create 100 replicates
    for (var in vars) { counter = counter + 1;
    
    #ensure that we have adequate levels, only a sanity since it is already filtered
    retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% pull(var)
    length(retained)
    #number of samples to subsample
    subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(min = min(n)) %>% pull(min)
    
    cat('Downsampling to n = ',subsamp,', requiring min = ',minnum,' for variable: ',var,', ','replicate: ',rep,'\n')
    #subsampling
    mdi = md_subbed %>%
      filter(!!sym(var) %in% retained) %>%
      group_by(!!sym(var)) %>%
      sample_n(min(n(), subsamp),replace = TRUE) %>%
      ungroup()
    mdi %>% count(!!sym(var))
    
    #ensure the factor levels are dropped which aren't in the dataset after filtering
    mdi = droplevels(mdi)
    
    #First MNLR on combinations
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
    for (model in models) {
      #output model fit from confusion matrix
      mo = get(model)
      
      #get AIC
      final_model = mo$finalModel;
      AIC = AIC(final_model)
      
      #save the model results
      dat = data.frame(Model = model, Iteration = rep, Variable = var, Subsampled = subsamp, MinObs = minnum,
                       decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, AIC = AIC, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
      dat_best = dat %>% slice_min(logLoss)
      adat = rbind(adat,dat_best)
      
      #also save training confusion matrix
      pred = mo$pred
      pred_best = pred %>% filter(decay == dat_best$decay)
      
      #and predict against real data
      predicted_classes = predict(mo, newdata = mdi, type = "raw")
      new = mdi %>% select(reference = !!sym(var)) %>% mutate(predicted = predicted_classes)
      
      conf_new = confusionMatrix(new$predicted, new$reference)
      conf_real = as.data.frame(conf_new$table) %>% mutate(
        Prediction = gsub('_','\\. ',Prediction),
        Reference = gsub('_','\\. ',Reference)) %>%
        mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = minnum, AUC=dat_best$AUC,logloss=dat_best$logLoss,Accuracy = dat_best$Accuracy,AccuracySD=dat_best$AccuracySD)
      new_preds = rbind(new_preds,conf_real)
      rm(conf_real,dat,dat_best)
      
     } #exit model loop
    } #exit response variable loop
  } #exit iteration loop
} #exit minimum samples loop

write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Unrelated-SIMPLE_2024MAR17.txt',quote=F,sep='\t',row.names=F)
write.table(new_preds,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix-Overall_Unrelated-SIMPLE_2024MAR17.txt',quote=F,sep='\t',row.names=F)

adat = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_Unrelated-SIMPLE_2024MAR17.txt')
conf_mats = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/ConfusionMatrix-Overall_Unrelated_2024MAR17.txt')

#plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('W+K+B','W+B','W','W+K','K','K+B','B'))
model_dat = adat %>% left_join(.,leg)

#Plot AUC for each model across all the iterations
model_dat$Label = factor(model_dat$Label, levels=c('W+K+B','W+B','W+K','K+B','B','W','K'))
app = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlab('AUC')+ylab('Model Covariates')+
  theme_bw(base_size=6)+theme(legend.position = 'top',strip.text.y.right = element_text(angle = 0))
app 

pdf('figures/MNLR_ModelSelection_Accuracy-Sensitivity_Unrelated-2024MAR17.pdf',height=4,width=6)
app
dev.off()

#plot boxes, AUC
auc_plot = model_dat %>%
  group_by(Variable,MinObs,Label,Iteration) %>% 
  ungroup %>% 
  ggplot(aes(x=Variable,fill=Label,y=AUC))+
  geom_boxplot(width=0.75,outlier.size = 0.5)+
  scale_fill_manual(values=c(brewer.pal(4,'Greys'),brewer.pal(3,'Set2')[c(2,3,1)]))+
  theme_bw(base_size=6)+
  theme(legend.position='top')
auc_plot

pdf('figures/MNLR_ModelSelection_Accuracy-Overall_Unrelated_2024MAR17.pdf',height=2,width=3.5)
auc_plot
dev.off()

#order full, single plot, make sure the 3 variables are in order 
egglevs = conf_mats %>% filter(Variable == 'Egg') %>% select(Prediction,Variable) %>% mutate(ord = as.numeric(gsub('E','',Prediction)),Prediction) %>% unique %>% arrange(Variable,ord) %>% select(-ord)
hostlevs = conf_mats %>% filter(Variable == 'Host') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
hablevs = conf_mats %>% filter(Variable == 'Habitat') %>% select(Prediction,Variable) %>% unique %>% arrange(Variable,Prediction) 
all_levs = rbind(hostlevs,egglevs,hablevs)
conf_mats = conf_mats %>% mutate(Prediction = factor(Prediction,levels=all_levs$Prediction),
                                 Reference = factor(Reference,levels=all_levs$Prediction))

#and also plot the from the re-predictions on the subsampled data 
repredictions_dat = conf_mats %>% 
  filter(Model == 'm2') %>%  #main figure, plot model 2 (W+B) 
  group_by(Variable,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq), #within each of the 100 replicates, calculate the frequency of correct classifications
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  sum_stats(Proportion)

#plot 
repredictions_plot = repredictions_dat %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=mean,col=sd))+
  geom_tile(col=NA)+
  geom_point(size=0.2)+
  scale_color_continuous(low='white',high='black')+
  scale_fill_continuous(low='white',high='darkblue')+  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
repredictions_plot

pdf('figures/MNLR_ConfusionMatrix-Repredictions-Unrelated_2024MAR17.pdf',height=3.5,width=3.5)
repredictions_plot
dev.off()
