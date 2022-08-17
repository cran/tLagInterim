#' @noRd
#' @include ipw.R type_Mean.R
#' @keywords internal
.ipw.RiskDiff <- function(response, data, khat) {
  
  # Estimates in closed form
  bterm <- data$delta * data$a / khat
  aterm <- data$delta / khat - bterm
  
  alpha <- sum(aterm * data$Y) / sum(aterm) 
  beta <- sum(bterm * data$Y) / sum(bterm) - alpha
  
  result <- list("alpha" = alpha, "beta"  = beta)
  
  class(x = result) <- "IPW_RiskDiff"
  
  return( result )
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.inflTerm.IPW_RiskDiff <- function(ipwObj, data, beta = NULL, ...) {
  
  # Use IPW alpha
  alpha <- ipwObj$alpha
  beta <- ifelse(test = is.null(x = beta), yes = ipwObj$beta, no = beta)
  
  abar <- mean(x = data$a)
  
  # {n}
  infl_term <- data$a * {data$Y - {beta + alpha}} / abar -
               {1.0 - data$a} * {data$Y - alpha} / {1.0 - abar}
  
  return( infl_term )
  
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.scorevec.IPW_RiskDiff <- .scorevec.IPW_Mean