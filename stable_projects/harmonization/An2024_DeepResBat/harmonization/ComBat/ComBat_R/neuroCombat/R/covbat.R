# Author: Andrew Chen, andrewac@pennmedicine.upenn.edu
# Adapted from code by Jean-Philippe Fortin, fortin946@gmail.com, available at
# https://github.com/Jfortin1/ComBatHarmonization
# The original and present code is under the Artistic License 2.0.
# If using this code, make sure you agree and accept this license.

# Correcting Covariance Batch Effects: CovBat
# Performs original ComBat to residualize and harmonize observations followed by
# PCA step to harmonize covariance across sites

#' Correcting Covariance Batch Effects: CovBat
#' 
#' \code{covbat} harmonizes covariance of observations across sites. It
#' first applies ComBat to harmonize mean and variance then adjusts variance of
#' PC scores to harmonize covariance.
#'
#' @param dat A \emph{p x n} matrix (or object coercible by
#'   \link[base]{as.matrix} to a numeric matrix) of observations where \emph{p}
#'   is the number of features and \emph{n} is the number of subjects.
#' @param bat Factor (or object coercible by \link[base]{as.factor} to a 
#'    factor) designating batch IDs.
#' @param mod Optional design matrix of covariates to preserve, usually from 
#'    \link[stats]{model.matrix}.
#' @param percent.var Numeric. The number of harmonized principal component 
#'    scores is selected to explain this proportion of the variance.
#' @param n.pc Optional numeric. If specified, this number of principal
#'    component scores is harmonized. Overrides \code{percent.var}.
#' @param train Optional logical vector specifying subset of observations used
#'    to estimate coefficients.
#' @param mean.only If \code{TRUE}, ComBat step does not harmonize variance.
#' @param std.var If \code{TRUE}, scales variances to be equal to 1 before PCA.
#' @param resid Whether to leave intercept and covariates regressed out.
#' @param eb If \code{TRUE}, uses ComBat model with empirical Bayes for mean
#'   and variance harmonization.
#' @param parametric If \code{TRUE}, uses parametric updates in ComBat step.
#' @param score.eb If \code{TRUE}, uses ComBat model with empirical Bayes for
#'   harmonization of scores.
#' @param score.parametric If \code{TRUE}, uses parametric ComBat is used
#'   for harmonization of scores.
#' @param verbose Whether additional details are printed to console.
#'
#' @return \code{covbat} returns a list containing the following components:
#' \item{dat.covbat}{Harmonized data as a matrix with same dimensions as 
#'    \code{x}}
#' \item{combat.out}{List output of \link[CovBat]{combat} from the ComBat step}
#' \item{combat.scores}{List output of \link[CovBat]{combat} from the score
#'    harmonization step}
#' \item{n.pc}{Number of principal components harmonized}
#' \item{prcomp.out}{PCA results after ComBat, output of \link[stats]{prcomp}}
#' @export
#' 
#' @examples
covbat <- function(dat, bat, mod = NULL, percent.var = 0.95, n.pc = NULL,
                   train = NULL, mean.only = FALSE, std.var = TRUE, 
                   resid = FALSE, eb = TRUE, parametric = TRUE,
                   score.eb = FALSE, score.parametric = TRUE, verbose = FALSE)
{
  dat <- as.matrix(dat)
  
  .checkConstantRows <- function(dat){
    sds <- rowSds(dat)
    ns <- sum(sds==0)
    if (ns>0){
      message <- paste0(ns, " rows (features) were found to be constant 
                        across samples. Please remove these rows before 
                        running ComBat.")
      stop(message)
    }
  }
  .checkConstantRows(dat)
  if (eb){
    if (verbose) cat("[combat] Performing ComBat with empirical Bayes\n")
  } else {
    if (verbose) cat("[combat] Performing ComBat without empirical Bayes 
                     (L/S model)\n")
  }
  # make batch a factor and make a set of indicators for batch
  batch <- as.factor(bat)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])
  
  # Number of covariates or covariate levels
  if (verbose) cat("[combat] Adjusting for",ncol(design)-ncol(batchmod),
                   'covariate(s) or covariate level(s)\n')
  
  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n.batch+1)){
      stop("[combat] The covariate is confounded with batch. Remove the 
           covariate and rerun ComBat.")
    }
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('The covariates are confounded. Please remove one or more of the 
             covariates so the design is not confounded.')
      } else {
        stop("At least one covariate is confounded with batch. Please remove 
             confounded covariates and rerun ComBat.")
      }
    }
  }
  
  ## Standardize Data across features
  if (verbose) cat('[combat] Standardizing Data across features\n')
  
  # Estimate coefficients using training set if specified, otherwise use full
  if (!is.null(train)) {
    design_tr <- design[train,]
    
    B.hat1 <- solve(crossprod(design_tr))
    B.hat1 <- tcrossprod(B.hat1, design_tr)
    B.hat <- tcrossprod(B.hat1, dat[,train])
  } else {
    B.hat1 <- solve(crossprod(design))
    B.hat1 <- tcrossprod(B.hat1, design)
    B.hat <- tcrossprod(B.hat1, dat)
  }
  
  # Standardization Model
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  
  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }	
  s.data <- (dat-stand.mean)/(tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  
  ##Get regression batch effect parameters
  if (eb){
    if (verbose) cat("[combat] Fitting L/S model and finding priors\n")
  } else {
    if (verbose) cat("[combat] Fitting L/S model\n")
  }
  batch.design <- design[,1:n.batch]
  gamma.hat <- tcrossprod(solve(crossprod(batch.design, batch.design)), 
                          batch.design)
  gamma.hat <- tcrossprod(gamma.hat, s.data)
  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,rowVars(s.data[,i], na.rm=TRUE)) # fixed error
  }
  
  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (eb){
    ##Find Priors
    #gamma.bar <- apply(gamma.hat, 1, mean)
    #t2 <- apply(gamma.hat, 1, var)
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apriorMat(delta.hat)
    b.prior <- bpriorMat(delta.hat)
    
    ##Find EB batch adjustments
    if (parametric){
      if (verbose) cat("[combat] Finding parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- it.sol(s.data[,batches[[i]]], gamma.hat[i,], delta.hat[i,],
                       gamma.bar[i], t2[i],a.prior[i], b.prior[i])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    } else {
      if (verbose) cat("[combat] Finding non-parametric adjustments\n")
      for (i in 1:n.batch){
        temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),gamma.hat[i,], 
                           delta.hat[i,])
        gamma.star <- rbind(gamma.star,temp[1,])
        delta.star <- rbind(delta.star,temp[2,])
      }
    }
    
  } 
  
  if (mean.only) {
    delta.star <- array(1, dim = dim(delta.star))
  }
  
  ### Normalize the Data ###
  if (verbose) cat("[combat] Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    if (eb){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/
        tcrossprod(sqrt(delta.star[j,]), rep(1,n.batches[j]))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/
        tcrossprod(sqrt(delta.hat[j,]), rep(1,n.batches[j]))
    }
    j <- j+1
  }
  
  # save ComBat dataset
  comdata <- (bayesdata*(tcrossprod(sqrt(var.pooled), rep(1,n.array)))) +
    stand.mean
  
  bayesdata <- bayesdata * (tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  
  x_pc <- prcomp(t(bayesdata), center = TRUE, scale. = std.var)
  
  # Subset scores based on percent of variance explained
  npc <- which(cumsum(x_pc$sdev^2/sum(x_pc$sdev^2)) > percent.var)[1]
  if (!is.null(n.pc)) {npc <- n.pc}
  # print(npc)
  scores <- x_pc$x[,1:npc]
  
  # ComBat without covariates to remove site effect in score mean/variance
  scores_com <- combat(t(scores), bat, eb = score.eb, 
                              parametric = score.parametric)
  full_scores <- x_pc$x
  full_scores[,1:npc] <- t(scores_com$dat.combat)
  
  # Project scores back into observation space
  if (std.var) {
    x.covbat <- t(full_scores %*% t(x_pc$rotation)) * 
      matrix(x_pc$scale, dim(bayesdata)[1], dim(bayesdata)[2]) + 
      matrix(x_pc$center, dim(bayesdata)[1], dim(bayesdata)[2])
  } else {
    x.covbat <- t(full_scores %*% t(x_pc$rotation)) + 
      matrix(x_pc$center, dim(bayesdata)[1], dim(bayesdata)[2])
  }
  if (resid == FALSE) {
    x.covbat <- x.covbat + stand.mean
  }
  
  return(list(dat.covbat = x.covbat, 
              combat.out = list(dat.combat=comdata, s.data=s.data, 
                                gamma.hat=gamma.hat, delta.hat=delta.hat,
                                gamma.star=gamma.star, delta.star=delta.star,
                                gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, 
                                b.prior=b.prior, batch=batch, mod=mod,
                                stand.mean=stand.mean, 
                                stand.sd=sqrt(var.pooled)[,1],
                                B.hat=B.hat),
              combat.scores = scores_com,
              npc=npc, x.pc = x_pc))
}