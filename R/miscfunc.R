#####################################################################

##  Convenience functions

#####################################################################

# expit and logistic functions

#' @noRd
#' @importFrom stats plogis
#' @keywords internal
.expit <- function(x){  stats::plogis(q = x) }

#' @noRd
#' @importFrom stats qlogis
#' @keywords internal
.logit <- function(x){ stats::qlogis(p = x) }

#' Internal f to allow for possibility of additional inputs and to perform
#'   basic checks on the returned object.
#' @noRd
#' @keywords internal
.f <- function(f, x.data, ...) {
  
  args <- list("x.data" = x.data)

  fbasis <- tryCatch(expr = do.call(what = f, args = args),
                     condition = function(e) {
                       stop("`f()` must execute without warnings or errors; ",
                            "the following were returned\n",
                            e$message, call. = FALSE)
                     })
  
  if (!is.matrix(x = fbasis) || nrow(x = fbasis) != nrow(x.data)) {
    stop("`f` must be a function that returns an {nSubject x (nBasis + 1)} ",
         "matrix of basis functions.", 
         call. = FALSE)
  }
  
  cmf <- colMeans(x = fbasis)
  tst <- cmf > {1.0 - 1e-8} & cmf < {1.0 + 1e-8}
  
  if (!any(tst)) {
    stop("a column of 1s must be included in the matrix returned by `f()`",
         call. = FALSE)
  }
  
  return( fbasis )
  
}

#' Internal h to allow for possibility of additional inputs and to perform
#'   basic checks on the returned object.
#' @noRd
#' @keywords internal
.h <- function(h, b.data, x.data, t.data, times, ...) {
  
  args <- list("b.data" = b.data, 
               "x.data" = x.data,
               "t.data" = t.data,
               "times" = times)
  
  hbasis <- tryCatch(expr = do.call(what = h, args = args),
                     condition = function(e) {
                       stop("`h()` must execute without warnings or errors; ",
                            "the following were returned\n",
                            e$message, call. = FALSE)
                     })
  
  if (!is.array(x = hbasis) || length(x = dim(x = hbasis)) != 3L) {
    stop("`h()` must return an array of dimension {nData x nTimes x nBasis} ", 
         call. = FALSE)
  }
  
  h_dim <- dim(x = hbasis)
  if ({h_dim[1L] != nrow(x = b.data)} || {h_dim[2L] != length(x = times)}) {
    stop("`h()` must return an array of dimension {nData x nTimes x nBasis} ", 
         call. = FALSE)
  }
  
  return( hbasis )
  
}