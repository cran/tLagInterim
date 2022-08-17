#'
#' Print results from a tLagInterim object
#'
#' @param x A tLagInterimObj object, returned by tLagInterim().
#' @param ... Ignored.
#'
#' @returns No return value, called to display key results.
#' 
#' @name print
#' @examples
#' data(tLagIntCat)
#' # f basis functions#'
#' f <- function(x.data) {
#'        return( as.matrix(x = cbind(1.0, x.data)) )
#'      }
#'
#' # h basis functions#'
#' h <- function(b.data, x.data, t.data, times) {
#'
#'   # Number of basis functions L 
#'   # (note that the number of basis functions does not and cannot depend
#'   #  on the treatment group; `h` is called internally multiple times -- each
#'   #  call is for a single treatment group.)
#'   L <- 2
#'   
#'   # Number of subjects in data
#'   n_subjects <- nrow(x = b.data)
#'   
#'   # Number of time points
#'   n_times <- length(x = times)
#'    
#'   # Initialize array of basis functions for this trt
#'   h.basis <- array(data = 0.0, dim = c(n_subjects, n_times, L))
#'        
#'   # Indicator of still being in hospital at any censoring time
#'   lindicator <- outer(X = t.data$lu, Y = times, "<=") * {t.data$ldelta == 2L}
#'   h.basis[, , 1L] <- lindicator
#'      
#'   # Time from leaving hospital to obstime for those known to
#'   # leave hospital at each censoring time
#'   h.basis[, , 2L] <- {times - t.data$lu} * lindicator
#'           
#'   return( h.basis )
#' }
#'
#' # fit with only baseline covariates provided, categorical outcome, user-specified f, h
#' res <- tLagInterim(b.data = b.data.cat,
#'                    x.data = x.data.cat,
#'                    t.data = NULL,
#'                    outcome = "categorical",
#'                    f = f, 
#'                    h = h)
#' print(res)
#' 
#' @method print tLagInterimObj
#' @export
#' @importFrom R.utils printf
print.tLagInterimObj <- function(x, ...) {

  headerLog = "   %6.4f  %6.4f  %10.4f  %10.4f  %6.4f %7.4f %7.4f\n"
  headerLogTitle = "   %6s  %6s  %10s  %10s  %6s %8s %7s\n"

  cat("\n\n",
      x$nt," subjects in this analysis; ",
      round(x = 100*x$cens, digits = 1L), "% censored\n",
      sep = "")  

  for (j in 3L:length(x = x)) {
    
    this.est <- attr(x = x, which = "estimator")[j-2L]
    
    cat("\n\n", this.est, " estimator\n\n", sep = "")
    
    R.utils::printf(headerLogTitle,
                   "beta", "se", "lower .95", "upper .95", "T", "n_ess", "info")
    R.utils::printf(headerLog,
                    x[[ this.est ]]["beta"],
                    x[[ this.est ]]["se.beta"],
                    x[[ this.est ]]["lower.95"],
                    x[[ this.est ]]["upper.95"],
                    x[[ this.est ]]["T"],
                    x[[ this.est ]]["ness"],
                    x[[ this.est ]]["info"])
  }    
}