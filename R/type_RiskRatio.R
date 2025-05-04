#' @noRd
#' @include ipw.R type_Mean.R
#' @keywords internal
.ipw.RiskRatio <- function(response, data, khat, ...) {
  
  # Estimates in closed form
  bterm <- data$delta * data$a / khat
  aterm <- data$delta / khat - bterm
  
  alpha <- log(x = sum(aterm * data$Y) / sum(aterm))
  beta <- log(x = sum(bterm * data$Y) / sum(bterm)) - alpha
  
  result <- list("alpha" = alpha, "beta"  = beta)
  
  class(x = result) <- "IPW_RiskRatio"
  
  return( result )
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.inflTerm.IPW_RiskRatio <- function(ipwObj, data, beta = NULL, ...) {
  
  alpha <- ipwObj$alpha
  beta <- ifelse(test = is.null(x = beta), yes = ipwObj$beta, no = beta)
  
  ea <- exp(x = alpha)
  eb <- exp(x = beta)
  
  abar <- mean(x = data$a)
  
  infl_term <- data$a * {data$Y - ea * eb} / {ea * eb * abar} -
               {1.0 - data$a} * {data$Y - ea} / {ea * {1.0 - abar}}
  
  return( infl_term )
  
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.scorevec.IPW_RiskRatio <- .scorevec.IPW_Mean
