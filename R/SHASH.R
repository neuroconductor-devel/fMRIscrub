#' SHASH to normal data transformation
#' 
#' Transform SHASH-distributed data to normal-distributed data.
#'
#' @param x Numeric vector of data to transform.
#' @param mu Parameter that modulates the mean of \code{x}.
#' @param sigma Parameter that modulates the variance of \code{x}.
#'  Must be greater than zero. This parameter is on the logarithm scale. 
#' @param nu Parameter that modulates the skewness of \code{x}.
#' @param tau Parameter that modulates the tailweight of \code{x}.
#'  Must be greater than zero. This parameter is on the logarithm scale.
#' @param inverse Transform normal data to SHASH instead? Default: \code{FALSE}.
#'
#' @return The transformed data.
#' @importFrom fMRItools is_1 is_posNum
#'
SHASH_to_normal <- function(x, mu, sigma, nu, tau, inverse = FALSE){
  stopifnot(is.numeric(x))
  stopifnot(is_1(mu, "numeric"))
  stopifnot(is_posNum(sigma))
  stopifnot(is_1(nu, "numeric"))
  stopifnot(is_posNum(tau))
  stopifnot(is_1(inverse, "logical"))

  sigma <- exp(sigma)
  tau <- exp(tau)

  xtrans <- if (inverse) {
    # normal to SHASH
    (sigma * tau * sinh((asinh(x) + nu ) / tau)) + mu
  } else {
    # SHASH to normal
    sinh((tau * asinh((x - mu)/ (sigma*tau))) - nu)
  }
}

#' Robust outlier detection based on SHASH distribution
#' 
#' A robust outlier detection based on modeling the data as coming from a SHASH 
#'  distribution.
#'
#' @param x The numeric vector in which to detect outliers.
#' @param maxit The maximum number of iterations. Default: \code{10}. 
#'
#' @return A \code{"SHASH_out"} object, i.e. a list with components
#' \describe{
#'  \item{out_idx}{Indices of the detected outliers.}
#'  \item{last_iter}{Last iteration number.}
#'  \item{converged}{Logical indicating whether the convergence criteria was satisfied or not.}
#' }
#' 
#' @importFrom gamlss gamlssML coefAll
#'
SHASH_out <- function(x, maxit = 10){
  nL <- length(x)
  weight_new <- rep(TRUE, nL) # TRUE if not an outlier
  iter <- 0
  success <- 0

  repeat {
    iter <- iter + 1
    weight_old <- weight_new

    # Transform the data.
    mod <- gamlss::gamlssML(
      X~1, 
      family = "SHASHo2", 
      maxit = 10000, 
      weight = as.numeric(weight_new)
    )
    est <- gamlss::coefAll(mod) 
    x_norm <- SHASH_to_normal(
      data = x, 
      mu = est$mu, sigma = est$sigma, nu = est$nu, tau = est$tau, 
      inverse = FALSE
    )
    x_norm_med <- median(x_norm)
    
    # Detect outliers.
    MAD = (1.4826) * median(abs(x_norm - x_norm_med))
    lim_left = x_norm_med - 4 * MAD
    lim_right = x_norm_med + 4 * MAD
    weight_new <- (x_norm > lim_left) & (x_norm < lim_right) # TRUE for non-outliers

    # Check convergence.
    if (all.equal(weight_old, weight_new)) {
      success <- 1
      break
    } else if (iter > maxit) {
      break
    }
  }

  # Return results.
  list(
    out_idx = which(!weight_new),
    last_iter = iter,
    converged = success
  )
}