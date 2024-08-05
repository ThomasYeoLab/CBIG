# Author: Jean-Philippe Fortin, fortin946@gmail.com
# Date: July 14 2020
# Projet repo: github.com/Jfortin1/ComBatHarmonization
# This is a modification of the ComBat function code from the sva package that can be found at
# https://bioconductor.org/packages/release/bioc/html/sva.html 
# The original code is under the Artistic License 2.0.
# The present code is under the MIT license
# If using this code, make sure you agree and accept this license. 

#' @export
neuroCombat <- function(dat, 
  batch, 
  mod=NULL, 
  eb=TRUE, 
  parametric=TRUE,
  mean.only=FALSE,
  ref.batch=NULL,
  verbose=FALSE
){
  dat <- as.matrix(dat)
  .checkConstantRows(dat)
  .checkNARows(dat)
   ## Check for missing values
  hasNAs <- any(is.na(dat))
  if (hasNAs & verbose){
    cat(paste0("[neuroCombat] Found ", sum(is.na(dat)), " missing data values. \n"))
  }
  if(mean.only){
      if (verbose) cat("[neuroCombat] Performing ComBat with mean only\n")
  }
  if (eb){
      if (verbose) cat("[neuroCombat] Performing ComBat with empirical Bayes\n")
  } else {
      if (verbose) cat("[neuroCombat] Performing ComBat without empirical Bayes (L/S model)\n")
  }

  ##################### Getting design ############################
  dataDict <- getDataDict(batch, mod, verbose=verbose, mean.only=mean.only, ref.batch=ref.batch)
  design <- dataDict[["design"]]
  ####################################################################


  ##################### Data standardization #######################
  if (verbose) cat('[neuroCombat] Standardizing Data across features\n')
  stdObjects <- getStandardizedData(dat=dat, 
    dataDict=dataDict,
    design=design,
    hasNAs=hasNAs
  )
  s.data <- stdObjects[["s.data"]]
  ####################################################################



  ##################### Getting L/S estimates #######################
  if (verbose) cat("[neuroCombat] Fitting L/S model and finding priors\n")
  naiveEstimators <- getNaiveEstimators(s.data=s.data,
      dataDict=dataDict, 
      hasNAs=hasNAs,
      mean.only=mean.only
  )
  ####################################################################


  ######################### Getting final estimators ####################
  if (eb){
      if (parametric){
        if (verbose) cat("[neuroCombat] Finding parametric adjustments\n")}else{
        if (verbose) cat("[neuroCombat] Finding non-parametric adjustments\n")
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
  if (verbose) cat("[neuroCombat] Adjusting the Data\n")
  bayesdata <- getCorrectedData(dat=dat,
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

  return(list(dat.combat=bayesdata, estimates=estimates))
}



                       
#' @export
neuroCombatFromTraining <- function(dat, batch, estimates, mod=NULL,
verbose=FALSE){
  # cat("[neuroCombatFromTraining] In development ...\n")
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
      cat("[neuroCombatFromTraining] Using mean imputation to account for previous covariates adjustment in training dataset\n")
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
  bayesdata <- bayesdata*sqrt(var.pooled) + stand.mean
  return(bayesdata)
}

