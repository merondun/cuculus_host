#Host & Environment associations, MNLR, Model Selection
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
# # This initial section assigns ancestry 'K' values for related individuals according to their relatives, since they weren't included in the popgen analysis.
# md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt')
# 
# #only grab samples with known egg
# md = md %>% drop_na(Egg) %>%
#   select(ID = ID, Host = HostParentShort, Egg, Environment = Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist)
# unknowns = md %>% filter(is.na(Ancestry))
# 
# #add similar K for relatives which weren't included in pop gen analyses
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
# write.table(md_egg,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/MNLR_Input_Egg_RelatedK_2024MAR05.txt',quote=F,sep='\t',row.names=F)

md_egg = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/MNLR_Input_Egg_RelatedK_2024MAR05.txt')

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md_egg = md_egg %>% mutate_at(c('Host','Egg','Environment','Haplogroup','Ancestry','Geography'),as.factor)

# Define the variables and covariates
variables <- c("Egg", "Host", "Environment")
covariates <- c("Haplogroup", "Ancestry", "Geography")

# Create combinations of variables and covariates
combinations <- expand.grid(Variable = variables, Covariate = covariates, stringsAsFactors = FALSE)

# Function to perform Chi-square test and return a data frame
perform_chi_square <- function(var, cov) {
  # Filter the data for relevant observations
  filtered_md <- md_egg %>%
    select(!!sym(var), !!sym(cov)) %>%
    na.omit()  # Remove NA values which chisq.test cannot handle
  
  # Perform Chi-square test
  test_result <- chisq.test(table(filtered_md[[1]], filtered_md[[2]]))
  
  # Return a data frame with the results
  tibble(
    Variable = var,
    Covariate = cov,
    ChiSq = test_result$statistic,
    p = test_result$p.value
  )
}

# Apply the function to each combination and bind the results into a single data frame
results_chi = map2_df(combinations$Variable, combinations$Covariate, perform_chi_square)
results_chi
results_chi %>% ggplot(aes(x=Variable,col=Covariate,y=ChiSq))+
  geom_point()+
  theme_bw()

#only retain Host where we have at least 2 cuckoos 
md_egg = md_egg %>% group_by(Host) %>% mutate(TotalHost = n()) %>% ungroup %>% group_by(Environment) %>% mutate(TotalEnvironment = n()) %>% ungroup %>% group_by(Egg) %>% mutate(TotalEgg = n())  %>% ungroup
minobs=2
md = md_egg %>% filter(TotalHost >= minobs & TotalEnvironment >= minobs & TotalEgg >= minobs) 
md_egg %>% filter(!ID %in% md$ID)
write.table(md$ID,file='randomforest/Samples_Retained_Related_MNLR_2024MAR05.txt',quote=F,sep='\t',row.names=F,col.names=F)

##### Model Selection ##### 
#assess covariate importance with model selection, using MNLR 
vars = c('Environment','Host','Egg')
sampling = c('Median','Minimum')

#set up cross validation 
ctrl = trainControl(method = "cv",   #Use cross-validation
                    number = 5,      #number of folds
                    summaryFunction = multiClassSummary,  #caret now attempts to stratify levels 
                    classProbs = TRUE, #For classification models, to save class probabilities
                    savePredictions = TRUE) #save predictions for each fold

#change 'A. pal' to A_pal' 
md_cv = md %>% mutate(Host = gsub('\\. ','_',Host),
                      Egg = gsub('\\. ','_',Egg))

#For determining which predictors improve model:
adat = NULL; conf_dat = NULL; counter = 0
for (minnum in c(2,5)) { 
  #Filter at the very base level, ensuring that across egg / host / habitat we have the same individuals with representative sampling 
  md_subbed = md_cv %>% filter(TotalHost > minnum & TotalEnvironment >= minnum & TotalEgg >= minnum) 
  for (rep in seq(1,100,1)){ 
    for (type in sampling) {
      for (var in vars) { counter = counter + 1;
      
      retained = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% pull(var)
      
      if (type == 'Median') {
        subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(median = median(n)) %>% pull(median)
      } else {
        subsamp = md_subbed %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(min = min(n)) %>% pull(min)
      }
      
      cat('Downsampling to n = ',subsamp,' for variable: ',var,', and sampling with: ',type,' ','replicate: ',rep,'\n')
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
      m1 = train(formula_1, data = mdi, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_2 = as.formula(paste(var, "~ Haplogroup + Geography"))
      m2 = train(formula_2, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_3 = as.formula(paste(var, "~ Haplogroup "))
      m3 = train(formula_3, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_4 = as.formula(paste(var, "~ Haplogroup + Ancestry"))
      m4 = train(formula_4, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_5 = as.formula(paste(var, "~ Ancestry"))
      m5 = train(formula_5, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_6 = as.formula(paste(var, "~ Ancestry + Geography"))
      m6 = train(formula_6, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      formula_7 = as.formula(paste(var, "~ Geography"))
      m7 = train(formula_7, data = md_cv, method = "multinom", trControl = ctrl, metric = "Accuracy", trace = FALSE)
      
      models = c('m1','m2','m3','m4','m5','m6','m7')
      for (model in models) {
        #output model fit from confusion matrix
        mo = get(model)
        pred = mo$pred
        
        #confusion Matrix for the predictions across all folds
        conf_matrix = confusionMatrix(pred$pred, pred$obs)
        conf = as.data.frame(conf_matrix$table) %>% mutate(
          Prediction = gsub('_','\\. ',Prediction),
          Reference = gsub('_','\\. ',Reference))
        
        print(mo)
        
        #save the model results 
        dat = data.frame(Model = model, Iteration = rep, Variable = var, Type = type, Subsampled = subsamp, MinObs = minnum,
                         decay = mo$results$decay, logLoss = mo$results$logLoss, logLossSD = mo$results$logLossSD, AUC = mo$results$AUC, Accuracy = mo$results$Accuracy, AccuracySD = mo$results$AccuracySD, Kappa = mo$results$Kappa)
        adat = rbind(adat,dat)
        
        #also save confusion matrix
        conf_dat = rbind(conf_dat,conf %>% mutate(Model = model, Iteration = counter, Variable = var, Subsample = subsamp, MinObs = minnum, Type = type))
        
      }
      }
    }
  }
}

write.table(adat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_AUC_WithReplace-3ResponseRel_Output_2024MAR05.txt',quote=F,sep='\t',row.names=F)
write.table(conf_dat,'/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_AUC_WithReplace-3ResponseRel_Confusion_Matrix_2024MAR05.txt',quote=F,sep='\t',row.names=F)

adat = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_AUC_WithReplace-3ResponseRel_Output_2024MAR05.txt')
conf_dat = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_AUC_WithReplace-3ResponseRel_Confusion_Matrix_2024MAR05.txt')

#re-name egg levels
eggs = read_tsv('/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/manyhost_hunt/Egg_Lookup.txt')

nonegg = conf_dat %>% filter(Variable != 'Egg')
egg = conf_dat %>% filter(Variable == 'Egg')
egg_rename = left_join(egg,eggs %>% dplyr::rename(Prediction=Egg)) %>% 
	select(-Prediction) %>% dplyr::rename(Prediction = Eggtype) %>% 
	left_join(.,eggs %>% dplyr::rename(Reference=Egg)) %>% 
	select(-Reference) %>% dplyr::rename(Reference = Eggtype)

conf_dat = rbind(nonegg,egg_rename)
write.table(conf_dat,file='/dss/dsslegfs01/pr53da/pr53da-dss-0021/projects/2021__Cuckoo_Resequencing/vcfs/all_samples-2022_11/host/randomforest/Model_Selection_AUC_WithReplace-3ResponseRel_Confusion_Matrix_2024MAR05.txt',quote=F,sep='\t',row.names=F)

#plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('W+K+D','W+D','W','W+K','K','K+D','D'))
model_dat = adat %>% 
  left_join(.,leg) %>% 
  mutate(Facet = paste0(Type,'\nn = ',MinObs)) 

#prep for plots, showing sensitivity of AUC
auc_summary = model_dat %>% 
  group_by(Variable,Facet,Label,Iteration) %>% 
  slice_min(logLoss,n=1) %>% ungroup %>% 
  group_by(Variable,Facet,Label) %>% 
  sum_stats(AUC)  %>%
  mutate(lab = paste0(round(conf_low,2),' - ',round(conf_high,2)))
model_dat$Label = factor(model_dat$Label, levels=c('W+K+D','W+D','W+K','K+D','D','W','K'))
app = model_dat %>% ggplot(aes(y=Label,x=AUC,fill=Variable))+
  geom_boxplot()+
  #geom_text(data=auc_summary,aes(y=Label,x=0.55,label=lab),
  #          position=position_dodge(width=1),size=3,col='black')+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlim(c(0.5,1))+xlab('Accuracy')+ylab('Model Covariates')+
  facet_grid(Facet~Variable,scales='free',space='free')+
  theme_bw()+theme(legend.position = 'none',strip.text.y.right = element_text(angle = 0))
app 

pdf('figures/MNLR_ModelSelection_AUC_WithReplace-3Response_2024MAR12.pdf',height=7,width=5)
app
dev.off()

#final plot, overall across all models
auc_overall = model_dat %>% 
  group_by(Variable,Facet,Label,Iteration) %>% 
  filter(MinObs == 2 & Type == 'Median') %>% 
  slice_min(logLoss,n=1) %>% ungroup %>% 
  group_by(Variable,Label) %>% 
  sum_stats(AUC)  %>%
  mutate(lab = paste0(sprintf('%.2f',mean),' +/- ',sprintf('%.2f',sd)))
modp = model_dat %>%
  group_by(Variable,Facet,Label,Iteration) %>% 
  filter(MinObs == 2 & Type == 'Median') %>% 
  slice_min(logLoss,n=1) %>% ungroup %>% 
  ggplot(aes(y=Label,x=AUC,fill=Variable))+
  ggdist::stat_halfeye(point_interval=NULL,alpha = 0.8,justification = -0.1,normalize='groups')+
  geom_point(data=auc_overall, aes(x=median,fill=Variable,y=Label),pch=21,stroke=0.35,size=1)+
  scale_fill_manual(values=rev(viridis(3,option='cividis')))+
  theme_bw(base_size=7)+
  coord_cartesian(xlim=c(0.75,1))+
  theme(legend.position='top')
modp

#adjust = .25,width = .5,.width = 0,trim=TRUE,justification = 0, point_colour = NA,
#just summary, small size 
modp = auc_overall %>% 
  ggplot(aes(y=Label,x=mean,xmin=mean-sd,xmax=pmin(mean+sd,1),col=Variable))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(position=position_dodge(width=0.5),width=0.25)+https://ood-2.ai.lrz.de/rnode/cpu-011.ai.lrz.de/8984/graphics/39b3f2b0-1a75-4ed2-9ee1-aa7db9b6fc78.png
  xlab('ROC AUC (Mean +/- SD')+ylab('')+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  #facet_grid(.~Variable,scales='free')+
  theme_bw(base_size=7)+
  xlim(c(0.5,1.0))+
  theme(legend.position='top')
modp

pdf('figures/MNLR_ModelSelection_AUC-Overall_WithReplacement-3Response_2024MAR12.pdf',height=2,width=2.75)
modp
dev.off()

#plot confusion matrix
conf_plot_dat = conf_dat %>% filter(Model == 'm2') %>% 
  filter(MinObs == 2 & Type == 'Median') %>% 
  group_by(Variable,Type,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq),
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  summarize(MedianFreq = median(Proportion),
            MeanFreq = mean(Proportion))  
conf_plotg = conf_plot_dat %>% 
  filter(Variable == 'Host') %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=MedianFreq))+
  geom_tile()+
  scale_fill_continuous(low='white',high='darkblue')+
  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
conf_plotg

eggdat = conf_plot_dat %>% filter(Variable == 'Egg')
egglevs = eggdat %>% ungroup %>% select(Prediction) %>% unique %>% mutate(ord = as.numeric(gsub('E','',Prediction))) %>% arrange(ord) 
eggdat$Reference = factor(eggdat$Reference,levels=egglevs$Prediction)
eggdat$Prediction = factor(eggdat$Prediction,levels=egglevs$Prediction)
conf_plote = eggdat %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=MedianFreq))+
  geom_tile()+
  scale_fill_continuous(low='white',high='darkblue')+
  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
conf_plote

conf_ploth = conf_plot_dat %>% 
  filter(Variable == 'Environment') %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=MedianFreq))+
  geom_tile()+
  scale_fill_continuous(low='white',high='darkblue')+
  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
conf_ploth

pdf('figures/MNLR_ConfusionMatrix-AUC_Overall_WithReplacement-HOST_2024MAR12.pdf',height=2.5,width=2)
conf_plotg
dev.off()

pdf('figures/MNLR_ConfusionMatrix-AUC_Overall_WithReplacement-ENV_2024MAR12.pdf',height=1.5,width=1)
conf_ploth
dev.off()

pdf('figures/MNLR_ConfusionMatrix-AUC_Overall_WithReplacement-EGG_2024MAR12.pdf',height=2.25,width=1.75)
conf_plote
dev.off()

#what is overall predictability?
model_dat %>% 
  group_by(Variable,Facet,Label,Iteration) %>% 
  filter(MinObs == 2 & Type == 'Median') %>% 
  slice_min(logLoss,n=1) %>% ungroup %>% group_by(Variable,Label) %>% sum_stats(Accuracy) %>% data.frame

