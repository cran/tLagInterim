#'  Create Martingale Adjustment Term and Augmentation Terms
#'
#'  Creates the martingale contributions for each treatment group and
#'   place in list; also returns the adjustment term.  The function is
#'   called by other functions that use these quantities.
#'
#' @noRd
#' @param b.data A data frame containing the basic observed data
#'   on the subset of enrolled subjects that received the same treatment, with 
#'   columns for the subject identifier, subjID, 
#'   treatment indicator a, u, delta, and outcome Y, where Y is = 0 if not 
#'   available.
#'
#' @param x.data A data frame containing the baseline covariates and subject
#'   identifier, subjID. Must contain same subjIDs as provided in b.data.
#' 
#' @param t.data A data frame containing the time dependent covariate
#'   information in the form specified by the user. Must use the same subjIDs
#'   as provided in b.data.
#'
#' @param scorevec A numeric vector. The inverse-weighted full-data influence
#'   function for the treatment effect parameter estimator.
#' 
#' @param h The user-defined function creating the h basis functions
#'   (can be NULL). The first 4 inputs to this function must be: `b.data`, 
#'   `x.data`, `t.data`, and `times`.
#'
#' @returns A list object containing
#'
#' \item{infl}{The h-basis function "covariate" contributions for the one-step 
#'   estimator.}
#'
#' \item{adj}{The "adjustment term" contribution to the "dependent variable" 
#'   for the one-step estimator.}
#'
#' \item{L}{The number of h basis functions for each treatment group.}
#' 
#' @include miscfunc.R
#' @keywords internal
.infl <- function(b.data,
                  x.data,
                  t.data,
                  scorevec,
                  h, ...) {

  n_subset <- nrow(x = b.data)
    
  #  censoring times for this trt
  times <- b.data$u[b.data$delta == 0L]
  n_times <- length(x = times)
    
  if (n_times == 0L) {
    stop("no censored cases for tx = ", b.data$a[1L], "; verify data", 
         call. = FALSE)
  }

  # At risk and censoring counting processes
  # {n_subset x n_times}
  dNt <- outer(X = b.data$u, Y = times, FUN = "==") * {b.data$delta == 0L}
  # {n_times}
  dNtsum <- colSums(x = dNt)
    
  # {n_subset x n_times}
  Yt <- outer(X = b.data$u, Y = times, FUN = ">=")
  # {n_times}
  Ysum <- colSums(x = Yt) 

  # Get martingale increment for scorevec 
  # {n_subset x n_times}
  mu_khat <- matrix(data = colSums(x = scorevec * Yt) / Ysum,
                    nrow = n_subset,
                    ncol = n_times,
                    byrow = TRUE)
    
  # {n_subset x n_subset_0}
  dMC <- dNt - t( t(x = Yt) * dNtsum / Ysum)

  # {n_subset}
  adj <- rowSums(x = mu_khat * dMC)
    
  # Get basis functions h for this treatment if t.data is
  # provided
    
  if (is.null(x = t.data)) {
    return( list("infl" = 0.0, "adj" = adj, "L" = 0L) )
  }
      
  hbasis <- .h(h = h,
               b.data = b.data,
               x.data = x.data,
               t.data = t.data, 
               times = times, ...)
  
  L <- dim(x = hbasis)[3L]
  
  infl <- matrix(data = 0.0, nrow = n_subset, ncol = L)
  for (l in 1L:L) {
    thisL <- hbasis[, , l]
    lYsum <- colSums(x = Yt * thisL)
    l_lbar <- t(x = t(x = thisL) - lYsum / Ysum)
    infl[, l] <- rowSums(x = l_lbar * dMC)
  }
  
  return( list("infl" = infl, "adj" = adj, "L" = L) )
}
