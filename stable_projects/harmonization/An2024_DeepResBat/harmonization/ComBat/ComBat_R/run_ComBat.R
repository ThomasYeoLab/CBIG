# @Description:
# R function for performing harmonization
# @Author: anlijuncn@gmail.com
# @Version: 1.0
# Written by Lijun An and CBIG under MIT license:
# https://github.com/ThomasYeoLab/CBIG/blob/master/LICENSE.md


ComBatHarm <- function(working_dir, fold_input_path, fold_output_path, train_val_file, test_file, isRef, isSaveTrain){

  ###################### Set Environment ######################
  setwd(working_dir);
  source('neuroCombat/R/imports.R');
  source('neuroCombat/R/neuroCombat.R');
  source('neuroCombat/R/neuroCombat_helpers.R');
  library(matrixStats);
  ###################### Training & Validation ######################
  train_roi_csv = paste(train_val_file, "_ROI.csv",sep = "")
  train_roi = read.csv(file.path(fold_input_path, train_roi_csv))
  train_roi = as.matrix(train_roi)
  train_site_csv = paste(train_val_file, "_site.csv",sep = "")
  train_site = read.csv(file.path(fold_input_path, train_site_csv))
  train_site = as.matrix(train_site)
  train_cov_csv = paste(train_val_file, "_covariates.csv",sep = "")
  train_cov = read.csv(file.path(fold_input_path, train_cov_csv))
  train_cov = as.matrix(train_cov)
  # run ComBat
  if (isRef == 1){
    train_harm <- neuroCombat(dat=train_roi, batch=train_site, mod=train_cov, ref.batch = 0)}
  else{
    train_harm <- neuroCombat(dat=train_roi, batch=train_site, mod=train_cov, ref.batch = NULL)}
  
  # save harmonized training & validation data
  if (isSaveTrain == 1){
    train_harm_ROI = train_harm$dat.combat
    harm_train_csv = paste("harm_", train_val_file, "_ROI.csv", sep = "")
    write.csv(train_harm_ROI, file=file.path(fold_output_path, harm_train_csv))}
  # get estimator
  estimator = train_harm$estimates
  ###################### Test ######################
  test_roi_csv = paste(test_file, "_ROI.csv",sep = "")
  test_roi = read.csv(file.path(fold_input_path, test_roi_csv))
  test_roi = as.matrix(test_roi)
  test_site_csv = paste(test_file, "_site.csv",sep = "")
  test_site = read.csv(file.path(fold_input_path, test_site_csv))
  test_site = as.matrix(test_site)
  test_cov_csv = paste(test_file, "_covariates.csv",sep = "")
  test_cov = read.csv(file.path(fold_input_path, test_cov_csv))
  test_cov = as.matrix(test_cov)
  # run ComBat using estimator from training & validation data
  test_harm <- neuroCombatFromTraining(test_roi, test_site, estimates = estimator, mod = test_cov)
  # save harmonized Tetset data
  harm_test_csv = paste("harm_", test_file, "_ROI.csv", sep = "") 
  write.csv(test_harm, file=file.path(fold_output_path, harm_test_csv))
}