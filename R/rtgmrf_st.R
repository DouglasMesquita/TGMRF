#' @title Random observations of a Transformed Gaussian Markov Random Field
#'
#' @description Use this function to simulate a data of a spatio-temporal TGMRF.
#'
#' @param rowid Number of lines of a lattice. (Used if neigh = NULL)
#' @param colid Number of columns of a lattice. (Used if neigh = NULL)
#' @param n_var Number of variables in problem (for multivariate problems)
#' @param neigh A object of class nb
#' @param n Vector with size of binomial trials
#' @param rho_s Spatial dependence parameter
#' @param rho_t Temporal dependence parameter
#' @param rho_st Spatio-temporal dependence parameter
#' @param betas Coeficients
#' @param nu Dispersion parameter for gamma models
#' @param tau Vector of precisions of variables
#' @param type_data Depends of family.
#' 'lognormal', 'lognormal-precision', 'gamma-shape', 'gamma-scale', 'weibull-shape', 'weibull-scale' for Poisson family
#' 'beta-logit', 'beta-probit', 'beta-alpha', 'beta-beta' for Binomial family
#' @param family 'poisson' or 'binary'
#' @param seed A seed to reproduce the reuslts
#'
#' @return y Response variable
#' @return X Covariates matrix
#' @return neigh Neighborhood structure
#' @return Q Covariance matrix
#' @return mu Means (TGMRF)
#' @return eps Errors (GMRF)

rtgmrf_st <- function(rowid = 10, colid = 10, n_var = 2, X = NULL,
                      neigh = NULL, n_trials = NULL, E = NULL,
                      rho_s = 0.9,  rho_t = 0,  rho_st = 0,
                      betas = c(-0.1, 0.3, 0.8), intercept = T, nu = 2,
                      tau = 1,
                      type_data = 'gamma-shape', family = 'poisson', mat_type = 'car',
                      seed = 1){

  set.seed(seed)

  Wt <- abs(outer(1:n_var, 1:n_var, "-")) == 1

  if(is.null(neigh)){
    N <- rowid*colid*n_var
    n_reg <- rowid*colid
    neigh <- cell2nb(rowid, colid)
    Ws <- nb2mat(neigh, style='B')
    coord <- expand.grid(y.coord = 1:colid, x.coord = 1:rowid)
  } else{
    N <- length(neigh)*n_var
    n_reg <- length(neigh)
    Ws <- nb2mat(neigh, style='B')
    coord <- NULL
  }

  P <- length(betas)
  if(is.null(intercept)){
    if(is.null(X)){
      X <- scale(matrix(rnorm(n = P*N), ncol = P))
    }
    colnames(X) <- paste0('X', 1:P)
  } else{
    if(is.null(X)){
      X <- scale(matrix(rnorm(n = P*N), ncol = P))
    }

    X <- cbind(1, X)
    betas <- c(intercept, betas)

    colnames(X) <- c("(Intercept)", paste0('X', 1:P))
  }

  Q <- buildQ(Ws = Ws, Wt = Wt, tau = tau, rho_s = rho_s, rho_t = rho_t, rho_st = rho_st)
  sigma <- solve(Q)

  eps <- as.vector(rmvnorm(1, rep(0, N), sigma))
  eps <- eps - mean(eps)

  Xbeta <- as.vector(X%*%betas)

  if(family == 'poisson'){
    if(type_data == 'gamma-scale') mu <- qgamma(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), shape = nu, scale = exp(Xbeta)/nu)
    if(type_data == 'gamma-shape') mu <- qgamma(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), shape = exp(Xbeta)*nu, scale = 1/nu)
    if(type_data == 'gamma-precision') mu <- qgamma(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), shape = exp(Xbeta)^2*nu, scale = 1/(nu*exp(Xbeta)))
    if(type_data == 'lognormal') mu <- qlnorm(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), meanlog = Xbeta, sdlog = sqrt(diag(sigma)/nu))
    if(type_data == 'lognormal-precision') mu <- qlnorm(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), meanlog = Xbeta, sdlog = sqrt(1/nu))
    if(type_data == 'weibull-scale') mu <- qweibull(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), shape = 1/nu, scale = exp(Xbeta))
    if(type_data == 'weibull-shape') mu <- qweibull(pnorm(eps, rep(0, N), sd = sqrt(diag(sigma))), shape = exp(Xbeta), scale = 1/nu)

    y <- rpois(n = N, lambda = mu)
  }

  if(family == 'binary'){
    stop("It is still not ready!")
    if(type_data == 'beta-logit'){
      par <- exp(Xbeta)/(exp(Xbeta) + 1)
      mu <- qbeta(pnorm(eps, rep(0, N), sd = 1), nu*par, nu*(1-par))
    }
    if(type_data == 'beta-probit'){
      par <- pnorm(Xbeta)
      mu <- qbeta(pnorm(eps, rep(0, N), sd = 1), nu*par, nu*(1-par))
    }
    if(type_data == 'beta-alpha'){
      mu <- qbeta(pnorm(eps, rep(0, N), sd = 1), exp(Xbeta), nu)
    }
    if(type_data == 'beta-beta'){
      mu <- qbeta(pnorm(eps, rep(0, N), sd = 1), nu, exp(-Xbeta))
    }

    if(is.null(n_trials)) n_trials <- 1

    y <- rbinom(n = N, size = n_trials, prob = mu)
  }

  if(!is.null(intercept)){
    X <- X[,-1, drop = FALSE]
  }

  return(list(y = y, X = X, reg = rep(1:n_reg, each = n_var), var = rep(1:n_var, times = n_reg),
              neigh = neigh, coord = coord,
              Q = Q, mu = mu, eps = eps))
}
