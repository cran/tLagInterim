#' Group Sequential Methods for Interim Monitoring of Randomized
#'   Clinical Trials with Time-lagged Outcome
#'
#' Implements methods for estimation of treatment effect parameters to
#'   support interim monitoring of clinical trials in which the outcome is
#'   ascertained after a time lag, so that not all subjects enrolled at
#'   the time of an interim analysis will have the outcome available.
#'   The methods take advantage of all available data to increase the
#'   precision of the analysis and thus lead to potentially earlier
#'   stopping.
#'
#' The data at the time of the desired interim analysis at time "t"
#'   must be input in one required and two optional data frames. The
#'   required data frame contains the basic information on treatment
#'   assignment, whether or not the outcome is available and the time
#'   lag, and, if available, the outcome itself. The first optional data
#'   frame contains baseline covariate information.  The second optional
#'   data frame must contain information relevant to constructing
#'   time-dependent covariates, and its form is specified by the user;
#'   an example is provided.
#'
#' Three types of outcome are supported: (1) continuous, (2) binary,
#'   and (3) ordered categorical.  For a continuous outcome, the treatment
#'   effect parameter is the difference in treatment means.  For
#'   a categorical outcome, the treatment effect parameter is the log odds
#'   ratio under an assumed proportional odds model.  For a binary
#'   outcome, the treatment effect parameter can be one of (a) the risk
#'   difference (equivalent to the difference in treatment means), (b)
#'   logarithm of the risk ratio (log relative risk), or (c) log odds
#'   ratio.
#'
#' If the outcome is ordered categorical, the categories must be ordered such
#'   that the outcomes are "worse" as one progresses from the base level
#'   to the final level. 
#'
#' If the outcome is binary and its levels are not coded as 0, 1, the coding
#'   for the base level must be provided as input.  The outcome will be recast
#'   internally as 0, 1. The underlying models for each type of treatment 
#'   effect are models for the probability that Y = 1.  There must be at least 
#'   one subject with available outcome equal to each of 0 and 1.
#'
#' The basic analysis data frame b.data must contain the following 
#'   variables for each subject:
#' \describe{
#' \item{subjID}{An identifier unique to each subject.}
#' \item{a}{The treatment assignment indicator; treatments must be binary.}
#' \item{u}{The time lag T at which the outcome was ascertained, if it was 
#'          ascertained, or the censoring time on the scale of subject time.}
#' \item{delta}{The indicator of T <= C, so that the outcome is
#'              observed if delta = 1}
#' \item{Y}{The outcome if it is available (delta=1); otherwise Y should
#'          be set equal to zero (delta=0); thus, Y = delta times outcome}
#'}
#' 
#' Each column of the baseline covariate data frame x.data should be a
#'   baseline covariate.  Data must contain a `subjID` column that contains
#'   the same subject identifiers as used in `b.data`.
#' 
#' The time-dependent data frame must contain the information used to
#'   construct time-dependent covariates in a format that is input into
#'   the user-specified function h() that constructs the
#'   basis functions. As this data frame is only used to construct the h basis 
#'   functions, the format and contents are, for the most part, entirely up to 
#'   the user. The notable exception is that it must contain a `subjID` 
#'   column that contains the same subject identifiers as used in `b.data`.
#'   
#' The function `h` is called multiple times internally -- each call is for
#'   a single treatment group. The function is provided only the data for the 
#'   specific treatment group under consideration, e.g., when estimating the L
#'   basis functions for a = 0, the b.data, x.data, and t.data passed to h() 
#'   contain only the rows for subjects in the a = 0 treatment arm; further,
#'   the nt censoring times are only those for this subset of subjects.
#'
#' The returned object contains the information needed to conduct any
#'   desired interim analysis (information-based or fixed-sample-based)
#'   for efficacy or futility using standard interim analysis software
#'   that assumes the test statistic has independent increments, such as
#'   the R package ldbounds.
#'
#' @param b.data A data frame containing the basic observed data on the n
#'   enrolled subjects at the time of an interim analysis at time t, with
#'   columns with headers 
#'   * "subjID" (unique subject identifiers),
#'   *  "a" (treatment indicator),
#'   *  "u" (minimum of time lag or censoring time), 
#'   *  "delta" (time lag/censoring indicator), and 
#'   *  "Y" (the outcome if it is available, = 0 if not). 
#'
#' @param x.data A data frame whose columns are baseline covariates, which is
#'   input to the user-specified function f (see example) to create the M+1
#'   baseline basis functions f_0, f_1, ..., f_M, where f_0 = 1 for all
#'   subjects; f_0 must be created in the function f.  If not provided or NULL,
#'   the AIPW2 estimator will be computed if t.data and h are provided; 
#'   otherwise, only the IPW estimator will be computed. Must contain a column 
#'   with header "subjID" containing the unique subject identifiers.
#'
#' @param t.data A data frame in whatever form the user specifies containing 
#'   the time-dependent covariate information, which is input along with x.data 
#'   to the user-specified function h (see example) to create the time-dependent 
#'   basis functions h_l, l=1, ..., L.  These basis functions can involve both 
#'   baseline and time-dependent covariates.  If not provided or NULL, the IPW
#'   and AIPW1 estimators will be computed if x.data and f are provided; 
#'   otherwise, only the IPW estimator will be computed. Must contain a column 
#'   with header "subjID" containing the unique subject identifiers.
#'
#' @param outcome Choices are "continuous", "binary", or "categorical".  If 
#'   outcome = "categorical", for each category there must be at least one 
#'   subject with available outcome.  If outcome = "binary", there must be at 
#'   least one subject for each level. If outcome is not specified as one
#'   of "continuous", "binary", or "categorical", an error will be generated.
#'
#' @param trteff If outcome = "binary", must be provided; trteff = "risk.diff" 
#'   for risk difference, trteff = "risk.ratio" for the logarithm of the risk 
#'   ratio (log relative risk), and trteff = "odds.ratio" for the log odds ratio.  
#'   If outcome = "binary" but trteff is not provided, an error will be generated.  
#'   Ignored if outcome = "continuous" or "categorical."
#'
#' @param f A user-specified function taking the data frame x.data as input, 
#'   which returns an (n x M+1) matrix whose first column is all ones and 
#'   remaining columns are the M user-specified basis functions f_1, ..., f_M 
#'   for each subject (see example below).  If `x.data` is not provided, `f` is
#'   ignored.
#'
#' @param h A user-specified function taking the data frames b.data, x.data, 
#'   t.data and a vector of censoring times as input. This function must return
#'   an array of dimension n x nt x L, where n = number of rows of the passed
#'   input b.data, and nt = number of censoring times passed as input,  so that 
#'   the (i,j,l) element of h.basis is the value of the lth basis function h_l 
#'   at the jth censoring time for the ith subject (see example below). If t.data 
#'   is not provided, h is ignored. See Details for further information.
#'
#' @param baseTx Type depends on class of treatment data. Treatment will be
#'   converted to 0/1 internally, this input specifies the value of b.data$a
#'   that is the base (control) value.
#'     
#' @param baseY Used only for binary outcomes. Type depends on class of outcome
#'   data. Outcome will be converted to 0/1 internally, this input specifies 
#'   the value of b.data$Y that is the base (0) value.
#'     
#' @param ... Ignored.
#' 
#' @returns An S3 object of class tLagInterim containing a list of
#'   variable length depending on which estimators can be computed
#'   given the inputs.  The elements of the list have the following
#'   names:
#'
#' \item{nt}{The number of subjects enrolled at the time of the
#' interim analysis.}
#'
#' \item{cens}{The proportion of these subjects for whom the outcome
#'     is not available (i.e., the time lag is censored).}
#'
#' \item{IPW}{A data frame containing the IPW estimate of the
#'    treatment effect parameter, its standard error, a 95% Wald
#'    confidence interval for the treatment effect, the corresponding
#'    Wald test statistic, the effective sample size n_ESS(t) (for
#'    fixed-sample-based monitoring), and the information Inf(t) =
#'    1/(standard error)^2 (for information-based monitoring).}
#'
#' \item{AIPW1}{If x.data and f are provided, a data frame
#'     containing the same information as for the IPW estimator 
#'     for the AIPW1 estimator that incorporates baseline covariate
#'     information only.}
#'
#' \item{AIPW2}{If either (i) x.data and f are not provided
#'     and t.data and h are, or (ii) both x.data and f and t.data and
#'     h are provided, a data frame containing the same information as
#'     for the IPW estimator for the AIPW2 estimator that
#'     incorporates time-dependent covariate information (alone or in
#'     addition to baseline covariate information).}
#'
#' The S3 object has an additional attribute, "estimators", giving a
#' description of which estimators are computed.
#'
#' @include KM.R ipw.R onestep.R ness.R verifyInputs.R
#'
#' @export
#' 
#' @references Tsiatis AA and Davidian M, 
#'   Group sequential methods for interim monitoring of
#'   randomized clinical trials with time-lagged outcome. 
#'   https://arxiv.org/abs/2204.10739.
#'
#' 
#' @examples
#' 
#' # Baseline and time-dependent covariates provided, categorical outcome
#' data(tLagIntCat)
#'
#' # f (basis functions for main effects when x contains continuous and
#' # binary (0/1) covariates); a user-specified function could also
#' # include dummies for categorical covariates, interaction terms,
#' # functions of covariates, etc.
#'
#' f <- function(x.data) {
#'    f.basis <- cbind(1.0, data.matrix(frame = x.data))
#'    return( f.basis )
#' }
#'
#' # h as for the first two simulation scenarios in the paper
#' # (categorical outcome), where t.data has columns "lu" = time to
#' # leaving hospital, death, or censoring, which ever first, and
#' # "ldelta" = 0 (censored), 1 (death), 2 (left hospital).  The basis
#' # functions could also include baseline covariates, although that
#' # is not the case here. 
#' 
#' h <- function(b.data, x.data, t.data, times) {
#'
#'   # Number of basis functions L 
#'   # (note that the number of basis functions does not and cannot depend
#'   #  on the treatment group; `h` is called internally multiple times -- each
#'   #  call is for a single treatment group.)
#'   L <- 2
#'   
#'   # Number of subjects in the provided data
#'   n_data <- nrow(x = b.data)
#'   
#'   # Number of censoring times provided
#'   n_times <- length(x = times)
#'    
#'   # Initialize array of basis functions
#'   h.basis <- array(data = 0.0, dim = c(n_data, n_times, L))
#'        
#'   # Indicator of still being in hospital at any censoring time
#'   lindicator <- outer(X = t.data$lu, Y = times, "<=") * {t.data$ldelta == 2L}
#'   h.basis[, , 1L] <- lindicator
#'      
#'   obstime <- max(b.data$u)
#'   
#'   # Time from leaving hospital to obstime for those known to
#'   # leave hospital at each censoring time
#'   h.basis[, , 2L] <- {obstime - t.data$lu} * lindicator
#'           
#'   #  Return the basis functions
#'   return( h.basis )
#' }
#'
#' # Compute all of IPW, AIPW1, AIPW2
#' tLagInterim(b.data = b.data.cat,
#'             x.data = x.data.cat,
#'             t.data = t.data.cat,
#'             outcome = "categorical",
#'             f = f,
#'             h = h)
#'
#' # Compute IPW, AIPW1 only (no time-dependent covariates)
#' tLagInterim(b.data = b.data.cat,
#'             x.data = x.data.cat,
#'             t.data = NULL,
#'             outcome = "categorical",
#'             f = f,
#'             h = NULL)
#'
#' # Baseline and time-dependent covariates provided, binary outcome, risk ratio
#' data(tLagIntBin)
#' 
#' # Compute all of IPW, AIPW1, AIPW2
#' tLagInterim(b.data = b.data.bin,
#'             x.data = x.data.bin,
#'             t.data = t.data.bin,
#'             outcome = "binary",
#'             trteff = "risk.ratio",
#'             f = f,
#'             h = h)
#'
#' # Compute IPW, AIPW2 only (no baseline covariates)
#' tLagInterim(b.data = b.data.bin,
#'             x.data = NULL,
#'             t.data = t.data.bin,
#'             outcome = "binary",
#'             trteff = "risk.ratio",
#'             f = NULL,
#'             h = h)
#' 
#' 
#' # Baseline and time-dependent covariates provided, continuous outcome
#' data(tLagIntCont)
#' 
#' # h as for the third simulation scenario in the paper (continuous
#' # outcome), where t.data has 5 columns corresponding to the 5
#' # intended times at which longitudinal measures of the outcome are
#' # ascertained, and the last observed measure is carried forward to
#' # all future times if it is not available
#'
#' h <- function(b.data, x.data, t.data, times) {
#'    
#'   # Number of basis functions L 
#'   # (note that the number of basis functions does not and cannot depend
#'   #  on the treatment group; `h` is called internally multiple times -- each
#'   #  call is for a single treatment group.)
#'   L <- 1L
#'
#'   # Number of subjects in provided data
#'   n_data <- nrow(x = b.data)
#'   
#'   # Number of censoring times provided
#'   n_times <- length(x = times)
#'  
#'   ti <- c(0,4,12,24,52)
#'    
#'   # Initialize array of basis functions
#'   h.basis <- array(data = 0.0, dim = c(n_data, n_times, L))
#'        
#'   # last value at each censoring time
#'   # dropping 1st column as it contains subject ids.
#'   h.basis[, , 1L] <- t(apply(X = t.data[,-1L],
#'                              MARGIN = 1L,
#'                              FUN = function(u) { 
#'                                      u[findInterval(x = times, vec = ti)] 
#'                                    }))
#'
#'   #  Return the basis functions
#'   return( h.basis )
#' }
#'
#' # Compute all of IPW, AIPW1, AIPW2
#' tLagInterim(b.data = b.data.cont,
#'             x.data = x.data.cont,
#'             t.data = t.data.cont,
#'             outcome = "continuous",
#'             f = f,
#'             h = h)
tLagInterim <- function(b.data,
                        x.data = NULL,
                        t.data = NULL,
                        outcome = c("continuous", "binary", "categorical"),
                        trteff = c("risk.diff", "risk.ratio", "odds.ratio"),
                        ..., 
                        f = NULL,
                        h = NULL, 
                        baseTx = 0L,
                        baseY = 0L) {
  
  ## check b.data

  b.data <- .verifyDataFrame(df = b.data,
                             required = c("a", "u", "delta", "Y", "subjID"),
                             input = "b.data")

  # 'a' must be binary integer -- converted to 0/1 if non-integer
  b.data$a <- .verifyInteger01(x = b.data$a, 
                               base = baseTx, 
                               x_name = "`a`",
                               base_name = "`baseTx`")
  
  # 'delta' must be 0/1 -- do not convert if non-integer
  if (!is.numeric(x = b.data$delta)) {
    stop("`delta` must be integer 0/1;\n\t", 
         "provided a ", class(b.data$delta), " object",
         call. = FALSE)
  } else if (is.integer(x = b.data$delta)) {
    if (!all(b.data$delta %in% c(0L, 1L))) {
      stop("`delta` must be integer 0/1;\n\t", 
           "provided ", paste(sort(unique(b.data$delta)), collapse = ", "),
           call. = FALSE)
    }
  } else {
    tmp <- as.integer(x = round(x = b.data$delta, digits = 0L))
    if (!isTRUE(x = all.equal(target = tmp, current = b.data$delta))) {
      stop("`delta` must be integer 0/1;\n\t", 
           "provided ", paste(sort(unique(b.data$delta)), collapse = ", "),
           call. = FALSE)
    }
    b.data$delta <- tmp
  }

  # 'u' must be non-negative
  if (any(b.data$u < 0.0)) {
    stop("`u` must be non-negative", call. = FALSE)
  }
  
  ## check x.data
  
  x.data <- .verifyDataFrame(df = x.data,
                             required = c("subjID"),
                             input = "x.data")
  
  if (!is.null(x = x.data)) {
    # x.data and b.data must reference same subject set
    if (!all(x.data$subjID %in% b.data$subjID) ||
        !all(b.data$subjID %in% x.data$subjID)) {
      stop('subject ids of b.data and x.data do not agree', call. = FALSE)
    }
  }
  
  ## check t.data
  
  t.data <- .verifyDataFrame(df = t.data,
                             required = c("subjID"),
                             input = "t.data")
  
  if (!is.null(x = t.data)) {
    # t.data and b.data must reference same subject set
    if (!all(t.data$subjID %in% b.data$subjID) ||
        !all(b.data$subjID %in% t.data$subjID)) {
      stop('subject ids if b.data and t.data do not agree')
    }
  }
  
  # Identify the estimators to be calculated
  
  if (is.null(x = x.data) & is.null(x = t.data)) {
    
    message("no covariates provided; only the IPW estimator will be computed")
    estimators <- c("IPW")
    types <- 0L
    
  } else if (!is.null(x = x.data) & is.null(x = t.data)) {

    message("no time-dependent covariates provided; ",
            "only IPW and AIPW1 estimators will be computed")
    
    .verifyUserFunctions(func = f,
                         required = "x.data",
                         input = "f",
                         estimator = "AIPW1")
    
    estimators <- c("IPW", "AIPW1")
    types <- c(0L, 1L)
    
  } else if (is.null(x = x.data) & !is.null(x = t.data)) {

    message("no baseline covariates provided; ",
            "only IPW and AIPW2 estimators will be computed")
    
    .verifyUserFunctions(func = h,
                         required = c("b.data", "x.data", "t.data", "times"),
                         input = "h",
                         estimator = "AIPW2")
    
    estimators <- c("IPW", "AIPW2")
    types <- c(0L, 0L)
    
  } else {
    
    .verifyUserFunctions(func = f,
                         required = "x.data",
                         input = "f",
                         estimator = "AIPW1")
    
    .verifyUserFunctions(func = h,
                         required = c("b.data", "x.data", "t.data", "times"),
                         input = "h",
                         estimator = "AIPW2")
    
    estimators <- c("IPW", "AIPW1", "AIPW2")
    types <- c(0L, 1L, 1L)
    
  }
       
  # Check if outcome is one of the three supported types
  outcome <- match.arg(arg = outcome)
  
  if (outcome == "binary") {
    trteff <- match.arg(arg = trteff)
  }
  
  # establish S3 class structure for internal procedures
  
  response <- b.data$Y
  eff_class <- switch(EXPR = outcome,
                      "continuous" = "Mean",
                      "binary" = switch(EXPR = trteff,
                                        "risk.diff" = "RiskDiff",
                                        "risk.ratio" = "RiskRatio",
                                        "odds.ratio" = "OddsRatio"),
                      "categorical" = "OddsRatio")
  class(x = response) <- c(outcome, eff_class, class(x = response))

  b.data$Y <- .checkOutcome(outcome = response, 
                            delta = b.data$delta,
                            base = baseY)
  
  # Number of subjects 
  n <- nrow(x = b.data)
  
  # Kaplan-Meier estimates
  khat <- .km(data = b.data)

  # Get the IPW estimate (SE gotten from onestep) 
  ipwObj <- .ipw(response = response, data = b.data, khat = khat)

  # Get AIPW estimate(s) (if applicable) and IPW SE and effective sample sizes
  onestepObj <- .onestep(b.data = b.data, 
                         x.data = x.data, 
                         t.data = t.data, 
                         ipwObj = ipwObj, 
                         khat = khat, 
                         f = f, 
                         h = h)

  ca <- stats::qnorm(p = 0.975)
  
  .summarize <- function(beta, se, type) {
  
    # Get the effective sample size for the IPW estimator and the information (1/SE^2)
    ness <- .nEss(b.data = b.data,
                  x.data = x.data,
                  ipwObj = ipwObj,
                  beta = beta,
                  se.beta = se,
                  type = type,
                  f = f, 
                  khat = khat)
  
    return( data.frame("beta" = beta,
                       "se.beta" = se,
                       "lower.95" = beta - ca * se,
                       "upper.95" = beta + ca * se,
                       "T" = beta / se,
                       "ness" = ness,
                       "info" = 1.0 / se^2) )
  }
  
  beta <- c(onestepObj$IPW$beta, onestepObj$AIPW1$beta, onestepObj$AIPW2$beta)
  se_beta <- c(onestepObj$IPW$se, onestepObj$AIPW1$se, onestepObj$AIPW2$se)

  results <- mapply(FUN = .summarize, beta, se_beta, types)
  colnames(results) <- estimators

  to_return <- list("nt" = n, "cens" = mean(b.data$delta == 0L))
  for (i in seq_along(estimators)) {
    to_return[[ estimators[i] ]] <- unlist(results[,i])
  }

  class(to_return) <- "tLagInterimObj"

  # estimators attribute for printing
  attr(x = to_return, which = "estimators") <- estimators

  return( to_return )

} 