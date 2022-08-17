#' Compute Effective Sample Size
#'
#' Compute the effective sample size at the current interim analysis
#'   time, which can be used to compute the proportion of information
#'   and thus stopping boundaries (using appropriate software such as
#'   ldbounds) for fixed-sample-based monitoring.
#'
#' @param b.data A data frame containing the basic observed data on the n
#'   enrolled subjects at the time of an interim analysis at time t, with
#'   columns with headers 
#'   * "subjID" (unique subject identifiers),
#'   * "a" (treatment indicator),
#'   * "u" (minimum of time lag or censoring time), 
#'   * "delta" (time lag/censoring indicator), and 
#'   * "Y" (the outcome if it is available, = 0 if not). 
#'
#' @param x.data A data frame whose columns are baseline covariates, which is
#'   input to the user-specified function f (see example) to create the M+1
#'   baseline basis functions f_0, f_1, ..., f_M, where f_0 = 1 for all 
#'   subjects; f_0 must be created in the function f.  
#'
#' @param ipwObj An object inheriting from class "IPW_*" containing the initial
#'   IPW fit.
#'
#' @param beta The estimate of beta for which effective sample size is required.
#'
#' @param sebeta The standard error of beta.
#'
#' @param type Indicator if this is an IPW or AIPW estimator; type = 0
#'   is IPW, type = 1 is AIPW1 or AIPW2.
#' 
#' @param f The user-defined function creating the f basis functions
#'     (can be NULL).
#'
#' @returns The effective sample size.
#'
#' @noRd 
#' @importFrom stats lm fitted
#' @include miscfunc.R
#' @keywords internal
.nEss <- function(ipwObj,
                  b.data,
                  x.data,
                  beta,
                  se.beta,
                  khat,
                  type,
                  f, ...) {
  
  n <- nrow(x = b.data)
  
  infl_term <- .inflTerm(ipwObj = ipwObj, 
                         data = b.data, 
                         beta = beta, 
                         khat = khat)

  # If type=0, this is IPW, and this term is all we need.  If
  # type=1, then we need to include the term involving (A-pi) psi
  # f(X).  

  if (type == 1L) {
    
    Ybeta <- infl_term
    
    fbasis <- .f(f = f, x.data = x.data, ...)
    
    designmat <- fbasis * (b.data$a - mean(x = b.data$a))
    
    wts <- b.data$delta / khat

    regObj <- stats::lm(Ybeta ~ designmat-1, weights = wts)
    
    pred_values <- stats::fitted(object = regObj)
    
    infl_term <- infl_term - pred_values
  }

  # Variance of the full data influence function
  varm <- mean(x = b.data$delta * (infl_term^2) / khat)

  # Effective sample size
  n_ess <- varm / (se.beta^2)

  return( n_ess )
}