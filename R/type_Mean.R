#' @noRd
#' @include ipw.R
#' @keywords internal
.ipw.Mean <- function(response, data, khat) {
    
  tmp <- data$delta / khat
    
  # Estimates in closed form
  bterm <- tmp * data$a
  aterm <- tmp * (1.0 - data$a)
    
  alpha <- sum(aterm * data$Y) / sum(aterm) 
  beta <- sum(bterm * data$Y) / sum(bterm) - alpha
    
  result <- list("alpha" = alpha, "beta"  = beta)
  
  class(x = result) <- "IPW_Mean"
  
  return( result )
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.inflTerm.IPW_Mean <- function(ipwObj, data, beta = NULL, ...) {
  
  alpha <- ipwObj$alpha
  beta <- ifelse(test = is.null(x = beta), yes = ipwObj$beta, no = beta)
  
  abar <- mean(x = data$a)
  
  infl_term <- data$a * {data$Y - {beta + alpha}} / abar -
               {1.0 - data$a} * {data$Y - alpha} / {1.0 - abar}
  
  return( infl_term )
  
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.scorevec.IPW_Mean <- function(ipwObj, data, khat, beta = NULL, ...) {
  
  inflterm <- .inflTerm(ipwObj = ipwObj, data = data, beta = beta)
  
  scorevec <- inflterm * data$delta / khat
  
  return( scorevec )
}
