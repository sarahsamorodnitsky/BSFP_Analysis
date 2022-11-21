# Analysis Code for ``Bayesian Simultaneous Factorization and Prediction Using Multi-Omic Data"

The files in this repository contain the code to run all analyses used in the article titled ``Bayesian Simultaneous Factorization and Prediction Using Multi-Omic Data." 

* bsfp.R contains all functions needed to run the BSFP model. 
* bpmf_validation.R contains the validation study. Results from validation simulations can be found in the validation_results folder. 
* bpmf_simulation.R contains the simulation study to compare models on how they recovery latent variation structure and predict an outcome. 
* bpmf_imputation.R contains the simulation study to compare different single imputation approaches on how they impute missing values (entrywise missing, blockwise missing, and missing-not-at-random).
* hiv_copd_application_V2.R contains the code to run the data application, where the BSFP model is fit on proteomic and metabolomic data to predict percent predicted forced expiratory volume in 1-second (FEV1pp). 
* hiv_copd_alignment_sensitivity.R evaluates the sensitivity of the alignment to choice of pivot. 
