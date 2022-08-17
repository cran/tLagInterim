#' Compute the Initial IPWCC estimator betahat_init
#'
#' Compute the initial IPW estimator at an interim analysis at time t
#'
#' @param response An S3 object of class 'Mean', "RiskDiff", "RiskRatio",
#'   or "OddsRatio" indicating the type of treatment effect.
#'
#' @param data A data frame containing the basic observed data
#'   on the n enrolled subjects at the time of an interim analysis
#'   at time t, with columns for the subject identifier, subjID, the 
#'   treatment indicator a, u, delta, and outcome Y, where Y is = 0 if not 
#'   available.
#'
#' @returns A list object of S3 class IPW
#' 
#' \item{alpha}{The IPW estimate(s) of nuisance parameters in the given model.}
#'
#' \item{beta}{The IPW estimate of the treatment effect in the given model.}
#'
#' \item{respmat}{Non-zero only for categorical outcomes. A matrix of indicator 
#'   functions.}
#'
#' @noRd
#' @keywords internal
.ipw <- function(response, ...) {
  UseMethod(".ipw", response)
}

.inflTerm <- function(ipwObj, ...) {
  UseMethod(".inflTerm", ipwObj)
}

.scorevec <- function(ipwObj, ...) {
  UseMethod(".scorevec", ipwObj)
}