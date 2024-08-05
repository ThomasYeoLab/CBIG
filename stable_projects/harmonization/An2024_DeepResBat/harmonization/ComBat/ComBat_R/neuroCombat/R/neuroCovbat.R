# Author: Jean-Philippe Fortin, fortin946@gmail.com
# Date: July 14 2020
# Projet repo: github.com/Jfortin1/ComBatHarmonization
# This is a modification of the ComBat function code from the sva package that can be found at
# https://bioconductor.org/packages/release/bioc/html/sva.html 
# The original code is under the Artistic License 2.0.
# The present code is under the MIT license
# If using this code, make sure you agree and accept this license. 

#' @export
neuroCovbat <- function(dat, 
  batch, 
  mod=NULL, 
  percent.var = 0.95,
  n.pc = NULL,
  eb=TRUE, 
  parametric=TRUE,
  mean.only=FALSE,
  ref.batch=NULL,
  std.var = TRUE,
  score.eb = FALSE, 
  score.parametric = TRUE,
  verbose=FALSE
){
  dat <- as.matrix(dat)
  .checkConstantRows(dat)
  .checkNARows(dat)
   ## Check for missing values
  hasNAs <- any(is.na(dat))
  if (hasNAs & verbose){
    cat(paste0("[Covbat] Found ", sum(is.na(dat)), " missing data values. \n"))
  }
  if(mean.only){
      if (verbose) cat("[Covbat] Performing ComBat with mean only\n")
  }
  if (eb){
      if (verbose) cat("[Covbat] Performing ComBat with empirical Bayes\n")
  } else {
      if (verbose) cat("[Covbat] Performing ComBat without empirical Bayes (L/S model)\n")
  }

  ##################### Getting design ############################
  dataDict <- getDataDict(batch, mod, verbose=verbose, mean.only=mean.only, ref.batch=ref.batch)
  design <- dataDict[["design"]]
  ####################################################################


  ##################### Data standardization #######################
  if (verbose) cat('[Covbat] Standardizing Data across features\n')
  stdObjects <- getStandardizedData(dat=dat, 
    dataDict=dataDict,
    design=design,
    hasNAs=hasNAs
  )
  s.data <- stdObjects[["s.data"]]
  ####################################################################



  ##################### Getting L/S estimates #######################
  if (verbose) cat("[Covbat] Fitting L/S model and finding priors\n")
  naiveEstimators <- getNaiveEstimators(s.data=s.data,
      dataDict=dataDict, 
      hasNAs=hasNAs,
      mean.only=mean.only
  )
  ####################################################################


  ######################### Getting final estimators ####################
  if (eb){
      if (parametric){
        if (verbose) cat("[Covbat] Finding parametric adjustments\n")}else{
        if (verbose) cat("[Covbat] Finding non-parametric adjustments\n")
      }
      estimators <- getEbEstimators(naiveEstimators=naiveEstimators, 
          s.data=s.data, 
          dataDict=dataDict,
          parametric=parametric,
          mean.only=mean.only
      )
  } else {
      estimators <- getNonEbEstimators(naiveEstimators=naiveEstimators, dataDict=dataDict)
  }
  ####################################################################


  ######################### Correct data #############################
  if (verbose) cat("[Covbat] Adjusting the Data with ComBat\n")
  combatdata <- getCorrectedData(dat=dat,
      s.data=s.data,
      dataDict=dataDict,
      estimators=estimators,
      naiveEstimators=naiveEstimators,
      stdObjects=stdObjects,
      eb=eb
  )
  ####################################################################


  # List of estimates:
  estimates <- list(gamma.hat=naiveEstimators[["gamma.hat"]], 
    delta.hat=naiveEstimators[["delta.hat"]], 
    gamma.star=estimators[["gamma.star"]],
    delta.star=estimators[["delta.star"]], 
    gamma.bar=estimators[["gamma.bar"]], 
    t2=estimators[["t2"]], 
    a.prior=estimators[["a.prior"]], 
    b.prior=estimators[["b.prior"]], 
    stand.mean=stdObjects[["stand.mean"]], 
    mod.mean=stdObjects[["mod.mean"]], 
    var.pooled=stdObjects[["var.pooled"]],
    beta.hat=stdObjects[["beta.hat"]],
    mod=mod, 
    batch=batch, 
    ref.batch=ref.batch, 
    eb=eb, 
    parametric=parametric, 
    mean.only=mean.only
  )

  ######################### CovBat #############################
  # Reference:
  # https://github.com/andy1764/CovBat_Harmonization/blob/master/R/R/covbat.R
  if (verbose) cat("[Covbat] Adjusting the Data with CovBat\n")
  bayesdata <- combatdata - estimates$stand.mean - estimates$mod.mean
  if (verbose) cat("[Covbat] Repalcing NaNs after ComBat with 0\n")
  bayesdata[is.na(bayesdata)] <- 0
  
  # perform PCA on normalized bayesdata
  x_pc <- prcomp(t(bayesdata), center = TRUE, scale. = std.var)
  npc <- which(cumsum(x_pc$sdev^2/sum(x_pc$sdev^2)) > percent.var)[1]
  if (!is.null(n.pc)) {npc <- n.pc}
  scores <- x_pc$x[,1:npc]
  # perform ComBat again on first npc principal components
  # ComBat without covariates to remove site effect in score mean/variance
  scores_combat <- neuroCombat(t(scores), batch,
    eb=score.eb, parametric=score.parametric)
  full_scores <- x_pc$x
  full_scores[, 1:npc] <- t(scores_combat$dat.combat)
  # Project scores back into observation space
  if (std.var) {
    x.covbat <- t(full_scores %*% t(x_pc$rotation)) * 
      matrix(x_pc$scale, dim(bayesdata)[1], dim(bayesdata)[2]) + 
      matrix(x_pc$center, dim(bayesdata)[1], dim(bayesdata)[2])
  } else {
    x.covbat <- t(full_scores %*% t(x_pc$rotation)) + 
      matrix(x_pc$center, dim(bayesdata)[1], dim(bayesdata)[2])
  }
  # add back
  x.covbat <- x.covbat + estimates$stand.mean + estimates$mod.mean
  
  return(list(
    dat.covbat=x.covbat, dat.pca_est=x_pc, dat.covcomb_est=scores_combat$estimates,
    dat.npc=npc,
    dat.combat=combatdata, dat.estimates=estimates))
}



                       
#' @export
neuroCovbatFromTraining <- function(
    dat, batch, estimates, pca_est, npc, covcomb_est, mod=NULL, verbose=FALSE){
  new.levels <- unique(batch)
  missing.levels <- new.levels[!new.levels %in% estimates$batch]
  if (length(missing.levels)!=0){
    stop(paste0("The batches ", missing.levels, " are not part of the training dataset\n"))
  }

  # Step 0: standardize data
  var.pooled <- estimates$var.pooled
  stand.mean <- estimates$stand.mean[,1]
  gamma.star <- estimates$gamma.star
  delta.star <- estimates$delta.star
  B.hat <- estimates$beta.hat
  n.array  <- ncol(dat)
  
  if (!is.null(estimates$mod)){
    if (verbose){
      cat("[CovbatFromTraining] Using mean imputation to account for previous covariates adjustment in training dataset\n")
    }
  }
  if (is.null(mod)){
    stand.mean <- stand.mean+rowMeans(mod.mean)
  } else {
    # stop("Including covariates for ComBat correction on a new dataset is not supported yet\n")
    ################################## Modification by Lijun ##################################
    # 1. get design matrix
    dataDict <- getDataDict(batch, mod, verbose=verbose, mean.only=FALSE, ref.batch=estimates$ref.batch)
    design <- dataDict[["design"]]
    n.batch=dataDict$n.batch
    batches <- dataDict$batches
    batch.design <- dataDict$batch.design
    n.batches <- dataDict$n.batches
    # 2. calculate mod.mean
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    mod.mean <- t(tmp%*%B.hat)
    stand.mean <- stand.mean + mod.mean
    ################################## End Modification ##################################
  }
  # Step 1: standardize data
  bayesdata <- (dat-stand.mean)/sqrt(var.pooled)
  # Step 2: remove estimates
  # gamma     <- t(gamma.star[batch,,drop=FALSE])
  # delta     <- t(delta.star[batch,,drop=FALSE])
  # bayesdata <- (bayesdata-gamma)/sqrt(delta)
  ################################## Modification by Lijun ##################################
  j <- 1
  for (i in batches){
    top <- bayesdata[,i]-t(batch.design[i,]%*%gamma.star)
    bottom <- tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    bayesdata[,i] <- top/bottom
    j <- j+1
  }
  ################################## End Modification ##################################
  # Step 3: transforming to original scale
  # bayesdata <- bayesdata*sqrt(var.pooled) + stand.mean 
  bayesdata <- bayesdata*sqrt(var.pooled)
  if (verbose) cat("[CovbatFromTraining] Repalcing NaNs after ComBat with 0\n")
  bayesdata[is.na(bayesdata)] <- 0
  # 1. perform PCA to get top npc PCs
  full_scores <- predict(pca_est, newdata=t(bayesdata))
  # 2. Harmonize top npc PCs using covcomb_est
  scores <- full_scores[, 1:npc]
  scores_combat <- (scores - covcomb_est$stand.mean[,1]) / sqrt(covcomb_est$var.pooled)
  sample_range <- 1:dim(batch)[1]
  for (j in sample_range){
    batch_j <- batch[j] + 1
    scores_combat[j,] <- (scores_combat[j,] - covcomb_est$gamma.star[batch_j,]) / sqrt(covcomb_est$delta.star[batch_j,])
  }
  scores_combat <- scores_combat * sqrt(covcomb_est$var.pooled) + covcomb_est$stand.mean[,1]
  full_scores2 <- predict(pca_est, newdata=t(bayesdata))
  full_scores2[, 1:npc] <- t(scores_combat)
  # 3. Project back to observation space
  x.covbat <- t(full_scores2 %*% t(pca_est$rotation)) * 
      matrix(pca_est$scale, dim(bayesdata)[1], dim(bayesdata)[2]) + 
      matrix(pca_est$center, dim(bayesdata)[1], dim(bayesdata)[2])
  # 4. Add back stand.mean
  bayesdata <- bayesdata + stand.mean 
  return(bayesdata)
}

