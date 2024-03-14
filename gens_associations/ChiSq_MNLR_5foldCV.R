#Gens & Habitat associations, MNLR, Model Selection
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
library(nnet)
library(broom)
library(meRo) #devtools::install_github('merondun/meRo',force=TRUE,upgrade='never')
library(caret)

set.seed(123)
md = read_tsv('~/merondun/cuculus_host/Metadata_Host.txt') 
md = md %>% filter(Analysis_Associations == 1) %>% 
  select(Gens, Habitat, Haplogroup = Hap, Ancestry = AncestryK5, Geography = KDist) 

##### Initialize, summary stats on raw data #####
#ensure they are all factor variables
md = md %>% mutate_at(c('Gens','Habitat','Haplogroup','Ancestry','Geography'),as.factor)

#only retain gens where we have at least 2 cuckoos 
md = md %>% group_by(Gens) %>% mutate(TotalGens = n()) %>% ungroup %>% group_by(Habitat) %>% mutate(TotalHabitat = n())
minobs=1
md = md %>% filter(TotalGens > minobs & TotalHabitat > minobs)

#raw proportions first, count proportions for Gens and habitat 
hp = md %>% group_by(Gens) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Gens,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique
tp = md %>% group_by(Habitat) %>% mutate(Total = n()) %>% 
  pivot_longer(c(Ancestry,Geography,Haplogroup)) %>% ungroup %>% 
  group_by(Habitat,name,value) %>% 
  summarize(Parts = n(),
            Proportion = Parts / Total,
            Percent = paste0(round( Parts / Total,3)*100,'% (',Total,')')) %>% unique

#bind them together for visualization 
ap = rbind(hp %>% ungroup %>% mutate(Response = Gens, variable = 'Gens') %>% select(-Gens), tp %>% ungroup %>% mutate(Response = Habitat, variable = 'Habitat') %>% select(-Habitat))
#just for ordering the covariates nicely 
ord = ap %>% ungroup %>% select(name,value) %>% unique %>% mutate(ord = as.numeric(str_sub(value,2))) %>% group_by(name) %>% arrange(ord)
ap$value = factor(ap$value,levels=ord$value)

#plot the raw proportions 
app = ap %>% 
  ggplot(aes(x=Response,y=value,fill=Proportion,label=Parts))+
  facet_grid(name~variable,scales='free',space='free')+
  scale_fill_gradient('Proportion Observed',low='yellow',high='red')+
  geom_tile()+
  geom_text(size=1)+
  ylab('')+xlab('')+
  theme_bw(base_size=6)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),
        legend.position = 'right')
app

pdf('figures/Proportions_GensHabitat_2024FEB29_min3.pdf',height=3,width=7.5)
app
dev.off()

# Chi-square tests
chi_square_haplogroup = chisq.test(table(md$Gens, md$Haplogroup))
chi_square_ancestry = chisq.test(table(md$Gens, md$Ancestry))
chi_square_geography = chisq.test(table(md$Gens, md$Geography))

data.frame(Haplogroup=chi_square_haplogroup$statistic,Ancestry=chi_square_ancestry$statistic,Geography=chi_square_geography$statistic)
# Haplogroup Ancestry Geography
# X-squared   403.5243 301.5551  480.0857

##### Model Selection ##### 
#assess covariate importance with model selection, using MNLR 
vars = c('Habitat','Gens')
sampling = c('Median','Minimum')

#set up cross validation 
ctrl = trainControl(method = "cv",   # Use cross-validation
                     number = 5,      # Number of folds
                     summaryFunction = multiClassSummary,  # Change if you need a different summary metric
                     classProbs = TRUE, # For classification models, to save class probabilities
                     savePredictions = TRUE) # To save predictions for each fold

#change 'A. pal' to A_pal' 
md_cv = md %>% mutate(Gens = gsub('\\. ','_',Gens))

#For determining which predictors improve model:
adat = NULL; conf_dat = NULL; counter = 0
for (rep in seq(1,100,1)){ 
  for (minnum in c(2,3,5)) { 
    for (type in sampling) {
      for (var in vars) { counter = counter + 1;
      
      retained = md_cv %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% pull(var)
      
      if (type == 'Median') {
        subsamp = md_cv %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(median = median(n)) %>% pull(median)
      } else {
        subsamp = md_cv %>% ungroup %>% count(!!sym(var)) %>% filter(n >= minnum) %>% summarize(min = min(n)) %>% pull(min)
      }
      
      cat('Downsampling to n = ',subsamp,' for variable: ',var,', and sampling with: ',type,' ','replicate: ',rep,'\n')
      #subsampling 
      mdi = md_cv %>%
        filter(!!sym(var) %in% retained) %>% 
        group_by(!!sym(var)) %>%
        sample_n(min(n(), subsamp)) %>%
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

write.table(adat,'~/merondun/cuculus_host/gens_associations/Model_Selection_Output_2024FEB29.txt',quote=F,sep='\t',row.names=F)
write.table(conf_dat,'~/merondun/cuculus_host/gens_associations/Model_Selection_Confusion_Matrix_2024FEB29.txt',quote=F,sep='\t',row.names=F)

adat = read_tsv('~/merondun/cuculus_host/gens_associations/Model_Selection_Output_2024FEB29.txt')
conf_dat = read_tsv('~/merondun/cuculus_host/gens_associations/Model_Selection_Confusion_Matrix_2024FEB29.txt')

#plot results 
leg = data.frame(Model = c('m1','m2','m3','m4','m5','m6','m7'),
                 Label = c('W+K+D','W+D','W','W+K','K','K+D','D'))
model_dat = adat %>% 
  pivot_longer(!c(Model,Iteration,Variable,Type, Subsampled, MinObs,AIC)) %>% 
  left_join(.,leg) %>% 
  mutate(Facet = paste0(Type,'\nn = ',MinObs)) 

#prep for plots, showing sensitivity of accuracy 
adatp = model_dat %>% filter(name == 'Accuracy') 
#summary, labels 
acc_sum = adatp %>% 
  group_by(Variable,Facet,Label,name) %>% sum_stats(value)  %>%
  mutate(lab = paste0(round(conf_low,2),' - ',round(conf_high,2)))
adatp$Label = factor(adatp$Label, levels=c('W+K+D','W+D','W+K','K+D','D','W','K'))
app = adatp %>% ggplot(aes(y=Label,x=value,fill=name))+
  geom_boxplot()+
  geom_text(data=acc_sum,aes(y=Label,x=-Inf,label=lab),
            position=position_dodge(width=1),hjust=-0.2,size=3,col='black')+
  scale_fill_manual(values=brewer.pal(3,'RdBu'))+
  xlim(c(0,1))+xlab('Accuracy')+ylab('Model Covariates')+
  facet_grid(Facet~Variable,scales='free',space='free')+
  theme_bw()+theme(legend.position = 'none',strip.text.y.right = element_text(angle = 0))
app 

pdf('figures/MNLR_ModelSelection_ACCURACY_2024FEB29.pdf',height=7,width=5)
app
dev.off()

#and AIC
aic_dat = model_dat %>% filter(name == 'Accuracy') %>% select(-name,-value, value = AIC) %>% unique
aic_dat$Label = factor(aic_dat$Label, levels=c('W+K+D','W+D','W+K','K+D','D','W','K'))
aic_sum = aic_dat %>% 
  group_by(Variable,Facet,Label) %>% sum_stats(value)  %>%
  mutate(lab = paste0(round(conf_low,0),' - ',round(conf_high,)),
         xpos = ifelse(Variable == 'Habitat',250, 500))
aip = aic_dat %>% ggplot(aes(y=Label,x=as.numeric(value)))+
  geom_boxplot(fill='cadetblue3')+
  geom_text(data=aic_sum,aes(y=Label,x=xpos,label=lab),
            position=position_dodge(width=1),hjust=1,size=3,col='black')+
  xlab('AIC')+ylab('Model Covariates')+
  facet_grid(Facet~Variable,scales='free')+
  theme_bw()+theme(legend.position = 'none',strip.text.y.right = element_text(angle = 0))
aip 

pdf('figures/MNLR_ModelSelection_AIC_2024FEB29.pdf',height=7,width=5)
aip
dev.off()

#final plot 
modp = rbind(adatp %>% select(-name,-AIC) %>% mutate(metric = 'Accuracy'),aic_dat %>% mutate(metric = 'AIC')) %>% 
  #filter(MinObs == 2 & Type == 'Minimum') %>%
  ggplot(aes(y=Label,x=as.numeric(value),fill=metric))+
  geom_boxplot(outlier.size = 1)+xlab('')+ylab('')+
  scale_fill_manual(values=brewer.pal(3,'Set2'))+
  facet_grid(Variable~metric,scales='free')+
  theme_bw(base_size=8)+
  theme(legend.position='top')
modp

pdf('~/merondun/cuculus_host/gens_associations/MNLR_ModelSelection_Acc-AIC-Overall_2024FEB29.pdf',height=2.75,width=3)
modp
dev.off()

#plot confusion matrix
conf_plot_dat = conf_dat %>% filter(Model == 'm1') %>% 
  group_by(Variable,Type,MinObs,Iteration,Reference) %>% 
  mutate(TotalIteration = sum(Freq),
         Proportion = Freq / TotalIteration) %>% ungroup %>% 
  group_by(Prediction,Reference,Variable) %>% 
  summarize(MedianFreq = median(Proportion),
            MeanFreq = mean(Proportion))  
conf_plot = conf_plot_dat %>% 
  ggplot(aes(x=Prediction,y=Reference,fill=MedianFreq))+
  geom_tile()+
  facet_wrap(Variable~.,scales='free',nrow=2,ncol=1)+
  scale_fill_continuous(low='white',high='darkblue')+
  theme_bw(base_size=8)+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1),legend.position='top')
conf_plot

pdf('~/merondun/cuculus_host/gens_associations/MNLR_ConfusionMatrix-Overall_2024FEB29.pdf',height=4.5,width=2.75)
conf_plot
