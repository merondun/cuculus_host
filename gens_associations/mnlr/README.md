# Multinomial logistic regression

Create a MNLR with egg or host or habitat as the response variable and Geography + Autosomal K + Haplogroups as the predictors. In short:

* Only retain response variables where there are at least n=2 observations
* Downsample all response classes so that all classes have n=2 observations
* Fit 7 multinomial logistic regression models, each with n=100 bootstraps using all combinations of predictors
* Extract AUC, and use the model to predict response variable on the full dataset again (too small for unseen data prediction)
* Repeat the above procedure 100 times so that different downsampled observations are included 
* Determine which classes are predicted correctly (% correct) from the confusion matrix on real / predicted responses across bootstraps

Files:

`Model_Selection_Boot-2Obs_2024AUG06.txt`: Model selection outputs for different covariates.

`ConfusionMatrix_Boot-2Obs_2024AUG06.txt`: Reference and predictions (e.g. confusion matrix) from model outputs. 

Files with `NoW1W2W3` are the same as above, except excluding the ancient blue egg W1/W2/W3 haplogroup split. 