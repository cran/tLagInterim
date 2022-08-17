#####################################################################

## Get the design matrices associated with the specified basis
## functions to construct the augmentation terms for each treatment

#####################################################################

#' Construct the design matrices for the "regressions" needed to
#' construct the one-step AIPW1 and AIPW2 estimators
#'
#' @param b.data A data frame containing the basic observed data
#'   on the n enrolled subjects at the time of an interim analysis
#'   at time t, with columns for the subject identify, subjID, 
#'   treatment indicator a, u, delta, and outcome Y, where Y is = 0 if not 
#'   available.
#'
#' @param x.data A data frame containing the baseline covariates and subject
#'   identifier, subjID. Must contain same subjIDs as provided in b.data.
#' 
#' @param t.data A data frame containing the time dependent covariate
#'   information in the form specified by the user. Must use the same subjIDs
#'   as provided in b.data, but does not need to include all.
#'
#' @param scorevec the inverse-weighted full-data influence function
#'   for the treatment effect parameter estimator.
#' 
#' @param f The user-defined function creating the f basis functions
#'   (can be NULL). Function takes 1 input, the baseline covariate
#'   data.frame. Function must return a matrix of dimension (now(x.data) x fnum)
#'
#' @param h The user-defined function creating the h basis functions
#'   (can be NULL). Function takes 4 inputs, `b.data`, `x.data`, `t.data`, and
#'   `times`.
#'
#' @returns A list object
#'
#' \item{designmat}{The design matrix for the regression for AIPW2;
#' that for AIPW1 is a subset of the columns.}
#'
#' \item{adj}{A vector containing the adjustment terms for each subject.}
#'
#' \item{fnum}{The number of columns of `designmat` that correspond to
#' the f basis functions.}
#'
#' \item{L}{The number of h basis functions for each treatment group.}
#'
#' @include martingale.R miscfunc.R
#' @noRd
#' 
.augment <- function(b.data,
                     x.data,
                     t.data,
                     scorevec,
                     f,
                     h, ...) { 

  n <- nrow(x = b.data)
    
  # Get the design matrix for the first augmentation term if x.data
  # is provided
  if (!is.null(x = x.data)) {
    
    fbasis <- .f(f = f, x.data = x.data, ...)
    fnum <- ncol(x = fbasis)
    message(fnum, " f basis functions")
    
    designmat <- fbasis * (b.data$a - mean(x = b.data$a))

  } else {
    
    fbasis <- NULL
    fnum <- 0L
    designmat <- NULL
    
  }
  
  adj <- numeric(length = n)
  Ls <- integer(length = 2L)
  
  for (a in 0L:1L) {
    
    subset_b <- b.data$a == a
    subset_x <- x.data$subjID %in% b.data$subjID[subset_b]
    subset_t <- t.data$subjID %in% b.data$subjID[subset_b]
          
    # Get the martingale increments for the adjustment and
    # time-dependent augmentation terms
    haug <- .infl(b.data = b.data[subset_b, , drop = FALSE],
                  x.data = x.data[subset_x, , drop = FALSE],
                  t.data = t.data[subset_t, , drop = FALSE],
                  scorevec = scorevec[subset_b],
                  h = h, ...)
    Ls[a+1L] <- haug$L

    # Columns for the second augmentation term if t.data is provided
    if (!is.null(x = t.data)) {  
      tmp_matrix <- matrix(data = 0.0, nrow = n, ncol = haug$L)
      tmp_matrix[subset_b, ] <- haug$infl
      designmat <- cbind(designmat, tmp_matrix)
    }
       
    # Adjustment terms
    adj[subset_b] <- haug$adj
  }
  
  if (Ls[2L] != Ls[1L]) {
    stop("h() must turn the same number of basis functions for all treatment ",
         "categories.", call. = FALSE)
  } else {
    if (Ls[1L] > 0L) {
      message(Ls[1L], " h basis functions")
    }
  }
  
  if (is.null(x = designmat)) designmat <- 0.0

  # Return design matrix and adjustment vector; also numbers of
  # basis functions used; if fnum=L=0, there are no covariates of
  # either type
  return( list("designmat" = designmat,
               "adj" = adj,
               "fnum" = fnum,
               "L" = Ls[1L]) )
}
