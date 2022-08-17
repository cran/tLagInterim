#' @noRd
#' @importFrom stats na.omit
#' @include ipw.R miscfunc.R
#' @keywords internal
.ipw.OddsRatio <- function(response, data, khat) {
  
  # number of categories
  cats <- levels(x = data$Y)
  if (is.null(cats)) {
    n_cats <- length(unique(data$Y))
  } else {
    n_cats <- nlevels(x = data$Y)
  }
  n_alpha <- n_cats - 1L

  # Response indicator matrix inverse weighted; note that censored
  # subjects are zeroed out by delta {n x n_alpha}
  if (n_cats == 2L) {
    resp_mat <- matrix(data = data$Y, ncol = 1L)
    alpha <- .logit(x = mean(x = data$Y[data$a == 0L]))
  } else {
    resp_mat <- outer(X = data$Y, Y = cats[-n_cats], FUN = "<=")
    resp_mat[is.na(resp_mat)] <- FALSE
    alpha <- .logit(x = colMeans(x = resp_mat))
  }
  

  # Get starting values for proportional odds model parameters
  beta <- 0.0
  
  # Initialize gradient matrix
  grad <- matrix(data = 0.0, nrow = n_cats, ncol = n_cats)
  
  # Start Newton-Raphson iteration
  
  delta_khat <- data$delta / khat
  
  score <- rep(x = Inf, times = n_cats)
  tol <- 1e-5
  imax <- 20L
  iter <- 1L
  
  while (max(abs(x = score)) > tol && iter < imax) { 
    
    # {n_aplha}
    p_A0 <- .expit(x = alpha)
    p_A1 <- .expit(x = alpha + beta)
    
    # exipt_{ij} = expit(alpha_j + A_i beta)
    # {n x n_alpha}
    p <- t(x = outer(X = p_A1 - p_A0, Y = data$a, FUN = "*") + p_A0)
    
    # E{M_alpha_{j}} = sum_i Delta_i / K(U_i,A_i) (R_j - expit_{ij})
    # {n_alpha}
    M_alpha <- colSums(x = {resp_mat - p} * delta_khat)
    
    # E{M_beta = sum_i Delta_i / K(U_i,A_i) sum_j (R_j - expit_{ij})
    M_beta <- sum(data$a * {resp_mat - p} * delta_khat)
    
    # E{M}
    # {n_alpha+1L}
    score <- c(M_alpha, M_beta)  
    
    # Derivatives
    
    # d/dalpha M_alpha
    # {n x n_alpha}
    pm1d <- delta_khat * p * {1.0 - p}
    
    # E{d M / d(alpha, beta)}
    diag(x = grad) <- c(colSums(x = pm1d), sum(pm1d * data$a))
    
    # E{d/dalpha_j M_beta}
    # E{d/dbeta_k M_alpha}
    # {1 x n_alpha}
    tmp <- drop(t(x = data$a) %*% pm1d)
    
    grad[{1L:n_alpha}, n_alpha+1L] <- tmp
    grad[n_alpha+1L, {1L:n_alpha}] <- tmp

    # parameter update
    incr <- tryCatch(expr = solve(a = grad, b = score),
                     error = function(e) {
                       stop("unable to invert gradient in Newton-Raphson",
                            e$message, call. = FALSE)
                     })
    
    alpha <- alpha + incr[1L:n_alpha]
    beta <- beta + incr[n_alpha + 1L]
  }
  
  if (iter >= imax) warning("IPW: Newton-Raphson did not converge", call. = FALSE)
  
  result <- list("alpha" = alpha,
                 "beta"  = beta,
                 "respmat" = resp_mat)
  
  class(x = result) <- "IPW_OddsRatio"
  
  return( result )
  
}

#' @noRd
#' @include ipw.R miscfunc.R
#' @keywords internal
.inflTerm.IPW_OddsRatio <- function(ipwObj, data, khat, beta = NULL, ...) {
  
  # number of enrolled subjects
  n <- nrow(x = data)
  
  # Use IPW alpha
  alpha <- ipwObj$alpha
  beta <- ifelse(test = is.null(x = beta), yes = ipwObj$beta, no = beta)
  
  # number of categories, and -1; assumes all categories are represented in data
  n_cats <- nlevels(x = data$Y)
  n_alpha <- n_cats - 1L
  
  # {n_alpha x n}
  respmat0 <- t(x = ipwObj$respmat)
  
  # {n_alpha}
  p0 <-.expit(x = alpha)
  p00 <- p0 * {1.0 - p0}
  
  # {n_alpha}
  p1 <- .expit(x = alpha + beta)
  p11 <- p1 * {1.0 - p1}

  # {1}
  abar <- mean(x = data$a)

  # {n_alpha}
  pbar <- abar * p11 + {1.0 - abar} * p00

  # {n_alpha}
  grad <- abar * {1.0 - abar} * p11 * p00 / pbar

  # {n x n_alpha}
  infl_term <- data$a * t(x = {(respmat0 - p1) * {1.0 - abar} * p00 / pbar}) -
               {1.0 - data$a} * t(x = {(respmat0 - p0) * abar * p11 / pbar})

  return( rowSums(x = infl_term / sum(grad)) )
}

#' @noRd
#' @include ipw.R
#' @keywords internal
.scorevec.IPW_OddsRatio <- .scorevec.IPW_Mean