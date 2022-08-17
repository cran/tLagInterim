#' Verify Vector is Binary Integer
#' 
#' Verifies that x is binary and converts to integer 0/1.
#' 
#' @noRd
#' @param x A vector.
#' @param base A scalar. The value of x to be taken as the base level.
#' @param x_name A character. The variable name.
#' @keywords internal
.verifyInteger01 <- function(x, ...) {
  UseMethod(generic = ".verifyInteger01", object = x)
}

#' @noRd
#' @importFrom stats na.omit
#' @keywords internal
.verifyInteger01.integer <- function(x, base, x_name, base_name) {
  
  levs <- sort(x = unique(x = stats::na.omit(object = x)))
  
  switch(as.character(x = length(x = levs)),
         "0" = stop("`", x_name, "` is all NAs", call. = FALSE),
         "1" = stop("`", x_name, "` is single valued", call. = FALSE),
         "2" = NULL,
         stop("`", x_name, "` is not binary", call. = FALSE))
  
  if (is.null(x = base)) {
    stop("`", base_name, "` = NULL` is not a valid input", call. = FALSE)
  }

  base <- as.integer(x = base)
  
  if (!is.finite(x = base)) {
    stop("`", base_name, " = ", base, "` is not a valid input", call. = FALSE)
  }
  
  if (!any(x == base)) {
    stop("`", base_name, " =  ", base, "` is not present in `", x_name, "`", 
         call. = FALSE)
  }
  
  tmp_x <- integer(length = length(x = x))
  tmp_x[x == base] <- 0L
  tmp_x[x != base] <- 1L
  tmp_x[is.na(x = x)] <- NA
    
  return( tmp_x )
}

#' @noRd
#' @keywords internal
.verifyInteger01.numeric <- function(x, x_name, ...) {
  
  tmp <- as.integer(x = round(x = x, digits = 0L))
    
  if (isTRUE(x = all.equal(target = x, current = tmp))) {
    return( .verifyInteger01(x = tmp, x_name = x_name, ...) )
  } else {
    stop("`", x_name, "` must be an integer", call. = FALSE)
  }
  
}

#' @noRd
#' @keywords internal
.verifyInteger01.factor <- function(x, base, x_name, base_name) {
  
  levs <- levels(x = x)
  base_level <- which(levs == base)

  if (length(base_level) == 0L) {
    stop("`", base_name, " = ", base, "` is not present in `", x_name, "`", 
         call. = FALSE)
  }
  
  x_unclassed <- unclass(x = x)
  attributes(x = x_unclassed) <- NULL
  
  return( .verifyInteger01(x = x_unclassed,
                           base = base_level,
                           x_name = x_name,
                           base_name = base_name) )
}
      
#' @noRd
#' @keywords internal
.verifyInteger01.character <- function(x, ...) {
  return( .verifyInteger01(x = as.factor(x = x), ...) )
}
    
#' @noRd
#' @keywords internal
.verifyInteger01.logical <- function(x, ...) {
  return( .verifyInteger01(x = as.integer(x = x), ...) )
}

#' @noRd
#' @keywords internal
.verifyDataFrame <- function(df, required, input) {
  
  if (is.null(x = df)) return( NULL )
  
  # must a data frame 
  if (!is.data.frame(x = df)) { 
    if (is.matrix(x = df)) {
      df <- as.data.frame(x = df)
    } else {
      stop("`", input, "` must be a data.frame", call. = FALSE)
    }
  }
  
  # must be complete
  if (any(is.na(x = df))) {
    stop('all data.frames must be complete;\n\t', 
         "`", input, "` contains NAs", call. = FALSE)
  }
  
  # must contain required column headers
  if (!all(required %in% colnames(x = df))) {
    stop("`", input, "` must contain column headers ", 
         paste(required, collapse = ", "),
         call. = FALSE)
  }
  
  return( df )
  
}

#' @noRd
#' @keywords internal
.verifyUserFunctions <- function(func, required, input, estimator) {
  
  if (is.null(x = func)) return( NULL )
  
  # must a function
  if (!is.function(x = func)) { 
    stop("function `", input, "` must be provided to compute ",
         estimator, " estimator",
         call. = FALSE)
  }
  
  frmls <- names(x = formals(fun = func))
  if (!all(required %in% frmls)) {
    stop("function `", input, "` must have input argument(s) ",
         paste(required, collapse = ", "), 
         call. = FALSE)
  }
  
  return( NULL )
  
}

.checkOutcome <- function(outcome, ...) {
  UseMethod(".checkOutcome")
}

#' @noRd
#' @keywords internal
.checkOutcome.continuous <- function(outcome, ...) {
  
  if (!is.numeric(x = outcome)) {
    stop("`Y` must be numeric", call. = FALSE)
  }
  
  if (length(x = unique(x = outcome)) <= 3L) {
    warning("`outcome = \"continuous\"`, but Y` takes <= 3 unique values",
            call. = FALSE)
  }
  # strip off outcome type and treatment effect classes
  class(outcome) <- class(outcome)[-c(1L,2L)]
  
  return( outcome )
  
}

#' @noRd
#' @keywords internal
.checkOutcome.binary <- function(outcome, delta, base, ...) {
  
  n_values <- length(x = unique(x = outcome[delta == 1L]))

  if (n_values != 2L) {
    stop("`outcome = \"binary\"`, but Y` does not contain 2 unique values",
         call. = FALSE)
  }
  
  # convert delta = 1 cases to 0/1
  new_outcome <- .verifyInteger01(x = outcome[delta == 1L], 
                                  base = base, 
                                  x_name = "Y", 
                                  base_name = "baseY")
  
  result <- integer(length = length(x = outcome))
  result[delta == 1L] <- new_outcome
  
  return( result )
  
}

#' @noRd
#' @keywords internal
.checkOutcome.categorical <- function(outcome, delta, ...) {
  outcome[delta == 0L] <- NA
  
  if (!is.factor(x = outcome) || !is.ordered(x = outcome)) {
    outcome <- factor(outcome, ordered = TRUE)
  }
  n_values <- nlevels(x = outcome)
  if (n_values > 10L) {
    warning("identified more than 10 outcome categories", call. = FALSE)
  }

  # strip off outcome type and treatment effect classes
  class(outcome) <- class(outcome)[-c(1L,2L)]
  
  return( outcome )
}
