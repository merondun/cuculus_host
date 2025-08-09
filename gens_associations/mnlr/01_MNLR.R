#!/usr/bin/env rscript
args <- commandArgs(trailingOnly = TRUE)

# Egg associations, MNLR, Model Selection c optatus
#### Find associations between MT/W haplotypes and features 
setwd('~/EvoBioWolf/CUCKOO_gentes/gens_associations/mnlr/')
#.libPaths('~/mambaforge/envs/rfs/lib/R/library')
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
library(caret)

sp = args[1]

set.seed(123)
md = read_tsv('~/EvoBioWolf/CUCKOO_gentes/Metadata_Host.txt') 

#only grab samples with known egg (hash out filter for related individuals)
md_egg = md %>%
  filter(SpeciesShort == sp & Analysis_GensAssociations == 1 & Sex == 'F') %>%
  select(ID = ID, Egg, Haplogroup, Geography = GeographicGroup,Ancestry = AncestryA5)
md_egg %>% count(Egg)

md_egg = md_egg %>% group_by(Egg) %>% mutate(TotalEgg = n()) %>% ungroup
minobs=2
md_egg = md_egg %>% filter(TotalEgg >= minobs)
md_egg %>% count(Egg)

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

# Filter at the very base level, ensuring that across egg / Egg / habitat we have the same individuals with representative sampling
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
write.table(adat,paste0('20250729_Model_Selection_Boot-2Obs-FEMALES_',sp,'.txt'),quote=F,sep='\t',row.names=F)
write.table(new_preds,paste0('20250729_ConfusionMatrix_Boot-2Obs-FEMALES_',sp,'.txt'),quote=F,sep='\t',row.names=F)


