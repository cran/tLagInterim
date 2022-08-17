#' Treatment-specific Kaplan-Meier Estimates of the Censoring
#' Probabilities
#'
#' Kaplan-Meier estimates of the censoring probabilities for each
#'   subject in the data set, estimated separately for each treatment.
#'   Is used to calculate inverse probability of censoring weights.  
#'
#' @param data A data frame containing the basic observed data
#'   on the n enrolled subjects at the time of an interim analysis
#'   at time t, with columns for treatment indicator a, u, delta,
#'   and outcome Y, where Y is = 0 if not available
#'
#' @returns A vector of length the number of rows in b.data containing
#'   the Kaplan-Meier estimate of the censoring distribution at U(t)
#'   for each subject
#' 
#' @import survival
#' @noRd
.km <- function(data) {

  # Vector to return KM estimated probabilities for each subject,
  # estimated separately by treatment
  Khat <- numeric(length = nrow(x = data))
    
  # Loop through each treatment
  for (a in 0L:1L) {
    
    grp <- data$a == a
    
    ua <- data$u[grp]

    # Get ranks of ua to reordered estimated survival probabilities
    rank.ua <- rank(x = ua, ties.method = "first")

    # Kaplan-Meier estimate, get estimated survival probabilities
    # for each subject
    Khat[grp] <- summary(object = survfit(formula = Surv(u, 1L-delta) ~ 1,
                                          data = data,
                                          subset = grp), 
                         times = ua)$surv[rank.ua]

  }
    
  return( Khat )
}    
