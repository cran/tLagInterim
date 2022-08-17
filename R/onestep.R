#' Compute the One-Step AIPW1 and AIPW2 Estimators
#'
#' Compute the one-step AIPW1 and AIPW2 estimators at an interim analysis at
#'   time t; also compute the standard error for the IPW estimator.  The AIPW1 
#'   estimator uses only baseline covariate information to improve precision; 
#'   the AIPW2 estimator also includes time-dependent covariate information.
#'
#' @param b.data A data frame containing the basic observed data on the n
#'   enrolled subjects at the time of an interim analysis at time t, with
#'   columns for the subject identifier, subjID, treatment indicator a, u, 
#'   delta, and outcome Y, where Y is = 0 if not available.
#'
#' @param x.data A data frame containing the baseline covariates and subject
#'   identifier, subjID. Must contain same subjIDs as provided in b.data.
#' 
#' @param t.data A data frame containing the time dependent covariate
#'   information in the form specified by the user. Must use the same subjIDs
#'   as provided in b.data.
#'
#' @param ipwObj An object of class "IPW_*" containing the initial IPW fit.
#'
#' @param f The user-defined function creating the f basis functions (can be NULL)
#'
#' @param h The user-defined function creating the h basis functions (can be NULL)
#'
#' @returns A list object containing
#'
#' \item{beta.aipw1}{The AIPW1 estimate.}
#'
#' \item{se.beta.aipw1}{Standard error for the AIPW1 estimator.}
#'
#' \item{beta.aipw2}{The AIPW2 estimate.}
#'
#' \item{se.beta.aipw2}{Standard error for the AIPW2 estimator.}
#'
#' \item{se.beta.ipw}{Standard error for the IPW estimator.}
#'
#' @include augment.R miscfunc.R
#'
#' @noRd
#' @importFrom stats lm fitted
#' @keywords internal
.onestep <- function(b.data,
                     x.data,
                     t.data,
                     ipwObj,
                     khat,
                     f,
                     h, ...) {
  
  # number of subjects
  n <- nrow(x = b.data)
  
  scorevec <- .scorevec(ipwObj = ipwObj, data = b.data, khat = khat)

  # Get the design matrices for the augmentation
  aug_results <- .augment(b.data = b.data,
                          x.data = x.data,
                          t.data = t.data,
                          scorevec = scorevec,
                          f = f,
                          h = h, ...)

  # Augmented influence function
  scorevec <- scorevec + aug_results$adj
  
  # Get the IPW SE
  IPW <- list("beta" = ipwObj$beta, 
              "se" = drop(x = sqrt(x = crossprod(x = scorevec)) / n))

  # If there are baseline covariates, regress the influence function
  # on the design matrix for only baseline covariate basis functions
  # (AIPW1)
  if (!is.null(x = x.data)) {
    AIPW1 <- .onestep_beta_se(scorevec = scorevec,
                              designmat = aug_results$designmat[,1L:aug_results$fnum],
                              ipw_beta = ipwObj$beta)
  } else {
    AIPW1 <- NULL
  }
  
  # If there are time-dependent covariates, regress the influence
  # function on the full design matrix to get the AIPW (AIPW2); if
  # there are no baseline covariates, the design matrix has columns
  # only for time-dependent covariates
  if (!is.null(x = t.data)) {
    AIPW2 <- .onestep_beta_se(scorevec = scorevec,
                              designmat = aug_results$designmat,
                              ipw_beta = ipwObj$beta)
    
  } else {
    AIPW2 <- NULL
  }
  
  # Return results
  return( list("IPW" = IPW, "AIPW1" = AIPW1, "AIPW2" = AIPW2) )
}

.onestep_beta_se <- function(scorevec, designmat, ipw_beta) {
  
  n <- nrow(x = designmat)
  
  reg_Obj <- stats::lm(scorevec ~ designmat - 1)
  pred_values <- stats::fitted(object = reg_Obj)
  
  # One-step update and sandwich SE
  beta <- ipw_beta - sum(pred_values) / n
  se <- sqrt(x = crossprod(x = scorevec - pred_values)) / n
  
  return( list("beta" = beta, "se" = drop(x = se)) )
}