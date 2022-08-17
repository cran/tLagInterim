#' Toy Dataset With a Continuous Outcome For Illustration
#'
#' These data are provided for the purposes of illustrating the use of
#' the software when the outcome of interest is continuous.
#' Though these data were generated to mimic conditions of a clinical trial,
#' they should not be interpreted as representing true clinical trial data.
#' 
#' @usage data(tLagIntCont)
#' 
#' @format The dataset provides three data.frames: `b.data.cont` containing the 
#'   basic observed data on 245 enrolled subjects at the time of an interim 
#'   analysis, with columns with headers 
#'   * "subjID" - unique subject identifiers,
#'   *  "u" - minimum of time lag or censoring time, 
#'   *  "delta" - time lag/censoring indicator, and 
#'   *  "Y" - the outcome if it is available, = 0 if not. 
#'   *  "a" - treatment indicator;
#'   
#'   `x.data.cont` contains the baseline covariates for the 245 subjects.
#'   * "subjID" - unique subject identifiers,
#'   * "X1" - a continuous covariate;
#'   
#'   and `t.data.cont` contains time-dependent covariate information comprising
#'   "subjID" and 6 measurements of a single continuous covariate measured 
#'   at times (t1 = 0, t2 = 4, t3 = 12, t4 = 24, t5 = 52)
#'   
#'
#' @name tLagIntCont
#' @aliases b.data.cont x.data.cont t.data.cont
#' @keywords datasets
NULL

#' Toy Dataset With a Binary Outcome For Illustration
#'
#' These data are provided for the purposes of illustrating the use of
#' the software when the outcome of interest is binary.
#' Though these data were generated to mimic conditions of a clinical trial,
#' they should not be interpreted as representing true clinical trial data.
#'
#' @usage data(tLagIntBin)
#' 
#' @format Each dataset provides three data.frames: `b.data.bin` containing the 
#'   basic observed data on 722 enrolled subjects at the time of an interim 
#'   analysis at time 52 week, with columns with headers 
#'   * "subjID" - unique subject identifiers,
#'   *  "u" - minimum of time lag or censoring time, 
#'   *  "delta" - time lag/censoring indicator, and 
#'   *  "Y" - the outcome if it is available, = 0 if not. 
#'   *  "a" - treatment indicator;
#'   
#'   `x.data.bin` contains the baseline covariates for the 722 subjects.
#'   * "subjID" - unique subject identifiers,
#'   * "X1" - a continuous covariate;
#'   
#'   and `t.data.bin` contains time-dependent covariate information comprising
#'   * "subjID" - unique subject identifiers,
#'   * "lu" - time to leaving hospital, death, or censoring
#'   * "ldelta" (0), death (1), or left hosp (2)
#'   
#'
#' @name tLagIntBin
#' @aliases b.data.bin x.data.bin t.data.bin
#' @keywords datasets
NULL

#' Toy Dataset With a Categorical Outcome For Illustration
#'
#' These data are provided for the purposes of illustrating the use of
#' the software when the outcome of interest is categorical.
#' Though these data were generated to mimic conditions of a clinical trial,
#' they should not be interpreted as representing true clinical trial data.
#'
#' @usage data(tLagIntCat)
#' 
#' @format Each dataset provides three data.frames: `b.data.cat` containing the 
#'   basic observed data on 477 enrolled subjects at the time of an interim 
#'   analysis, with columns with headers 
#'   * "subjID" - unique subject identifiers,
#'   *  "u" - minimum of time lag or censoring time, 
#'   *  "delta" - time lag/censoring indicator, and 
#'   *  "Y" - the outcome if it is available, = 0 if not. 
#'   *  "a" - treatment indicator;
#'   
#'   `x.data.cat` contains the baseline covariates for the 477 subjects.
#'   * "subjID" - unique subject identifiers,
#'   * "X1" - a continuous covariate;
#'   
#'   and `t.data.cat` contains time-dependent covariate information comprising
#'   * "subjID" - unique subject identifiers,
#'   * "lu" - time to leaving hospital, death, or censoring
#'   * "ldelta" (0), death (1), or left hosp (2)
#'   
#'
#' @name tLagIntCat
#' @aliases b.data.cat x.data.cat t.data.cat
#' @keywords datasets
NULL