#' @title Spatio-temporal Transformed Gaussian Markov Random Field
#'
#' @description Uses this function to fit spatio-temporal data using non-guassian copulas.
#'
#' @param y Response variable
#' @param X Matrix or data.frame of covariates
#' @param n_reg Number of regions in study
#' @param n_var Number of variables in study
#' @param beta Vector of the initial values to coeficients vector
#' @param nu A Initial value to variance of copula
#' @param rho_s Spatial dependence parameter
#' @param rho_t Temporal dependence parameter
#' @param rho_st Spatio-temporal dependence parameter
#' @param tau Vector of precisions of variables
#' @param family 'Poisson' or 'Binary'
#' @param type Depends of family.
#' 'log-normal', 'gamma-shape', 'gamma-scale', 'weibull-shape', 'weibull-scale' for Poisson family
#' 'beta-logit', 'beta-probit', 'beta-alpha', 'beta-beta' for Binomial family
#' @param mat_type car or leroux
#' @param method 'arms' or 'metropolis'
#' @param nsim Number of MCMC iterations
#' @param burnin number of discards iterations
#' @param thin Lag to collect the observations
#' @param E Offsets for Poisson data
#' @param n Number of trials in Binomial data
#' @param neigh Neighborhood structure of nb class or a neighborhood matrix
#' @param prior_param List of priors parameters. Gaussian for beta vector and IGamma for nu. Entrys are: 'nu' and 'beta'
#' @param range Range to use during the MCMC. Entrys are: 'rho_s', 'rho_t' and 'rho_st'
#' @param MCMC_config Parameters to tunning MCMC
#' @param fix_rho A list informing if any rho (rho_s, rho_t, rho_st) is fixed
#' @param range A list informing a range for sampling rho (rho_s, rho_t, rho_st). Each entre must be a vector in the interval (-1, 1).
#'
#' @return beta A matrix with samples of beta
#' @return nu A vector with samples of nu
#' @return rho A vector with samples of rho
#'
#' @examples s_tgmrf(1,1,1,1)
#'
#' @export

st_tgmrf <- function(y, X, n_reg, n_var,
                     beta, nu, eps, mu,
                     rho_s, rho_t, rho_st,
                     tau,
                     family, type, mat_type, method,
                     nsim, burnin, thin,
                     E, n, neigh,
                     prior_param,
                     MCMC_config,
                     fix_rho, range,
                     verbose,
                     c_beta, c_eps, c_mu, c_nu, c_rho){

  P <- length(beta)
  N <- length(y)

  mean_beta <- prior_param$beta$mean
  if(!(length(mean_beta) %in% c(1,P))) stop('mean_beta must length 1 or the number of covariates')
  if(length(mean_beta) == 1) mean_beta <- rep(mean_beta, P)
  tau_beta <- prior_param$beta$precision
  if(!(length(tau_beta) %in% c(1,P))) stop('tau_beta must length 1 or the number of covariates')
  if(length(tau_beta) == 1) tau_beta <- rep(tau_beta, P)

  eta_nu <- prior_param$nu$shape
  psi_nu <- prior_param$nu$rate
  if(is.null(eps)){
    eps <- rep(0, N)
  }
  if(is.null(mu)){
    mu <- qgamma(p = 0.5, shape = rep(1, N), rate = 1)
  }

  fix_rho_s <- fix_rho$rho_s
  fix_rho_t <- fix_rho$rho_t
  fix_rho_st <- fix_rho$rho_st

  range_rho_s <- range$rho_s
  range_rho_t <- range$rho_t
  range_rho_st <- range$rho_st

  ninit <- MCMC_config$arms$ninit
  maxpoint <- MCMC_config$arms$maxpoint

  var_eps <- MCMC_config$metropolis$var_eps
  if(is.matrix(var_eps)){
    if(any(!(dim(var_eps) == c(N, N)))) stop(sprintf('var_eps must be a %sx%s matrix', N, N))
  } else{
    if(!(length(var_eps) %in% c(1, N))) stop('var_eps must length 1 or the number of elements')
    if(length(var_eps) == 1) var_eps <- rep(var_eps, N)
    var_eps <- diag(var_eps, ncol = N)
  }

  var_log_mu <- MCMC_config$metropolis$var_log_mu
  if(is.matrix(var_log_mu)){
    if(any(!(dim(var_log_mu) == c(N, N)))) stop(sprintf('var_log_mu must be a %sx%s matrix', N, N))
  } else{
    if(!(length(var_log_mu) %in% c(1, N))) stop('var_log_mu must length 1 or the number of elements')
    if(length(var_log_mu) == 1) var_log_mu <- rep(var_log_mu, N)
    var_log_mu <- diag(var_log_mu, ncol = N)
  }

  var_beta <- MCMC_config$metropolis$var_beta
  if(is.matrix(var_beta)){
    if(any(!(dim(var_beta) == c(P, P)))) stop(sprintf('var_beta must be a %sx%s matrix', P, P))
  } else{
    if(!(length(var_beta) %in% c(1,P))) stop('var_beta must length 1 or the number of covariates')
    if(length(var_beta) == 1) var_beta <- rep(var_beta, P)
    var_beta <- diag(var_beta, ncol = P)
  }
  var_log_nu <- MCMC_config$metropolis$var_log_nu
  var_rho <- MCMC_config$metropolis$var_rho
  if(is.matrix(var_rho)){
    if(any(!(dim(var_rho) == c(3, 3)))) stop('var_rho must be a 3x3 matrix')
  } else{
    if(!(length(var_rho) %in% c(1, 3))) stop('var_rho must length 1 or 3')
    if(length(var_rho) == 1) var_rho <- rep(var_rho, 3)
    var_rho <- diag(var_rho, ncol = 3)
  }

  if(ninit < 3 | !is.numeric(ninit))
    stop("ninit must be a numeric > 2")
  if(maxpoint < 1 | !is.numeric(maxpoint))
    stop("maxpoint must be a numeric > 0")
  if(var_log_nu <= 0)
    stop("var_log_nu must be a numeric > 0")
  if(rho_s < range_rho_s[1] | rho_s > range_rho_s[2])
    stop("rho_s must be in the interval defined in range$rho_s")
  if(rho_t < range_rho_t[1] | rho_t > range_rho_t[2])
    stop("rho_t must be in the interval defined in range$rho_t")
  if(rho_st < range_rho_st[1] | rho_st > range_rho_st[2])
    stop("rho_st must be in the interval defined in range$rho_st")

  if(family == "poisson"){
    if(type == "gamma-scale") type_num <- 1
    else if(type == "gamma-shape") type_num <- 2
    else if(type == "gamma-precision") type_num <- 3
    else if(type == "log-normal") type_num <- 4
    else if(type == "weibull-shape") type_num <- 5
    else if(type == "weibull-scale") type_num <- 6
    else stop("Family for the TGMRF not defined\n")

    out <- mcmc_poisson_st(y = y, X = X, E = E,
                           n_reg = n_reg, n_var = n_var, neigh = neigh,
                           beta = beta, nu = nu, eps = eps, mu = mu,
                           rho_s = rho_s, rho_t = rho_t, rho_st = rho_st, tau = tau,
                           nsim = nsim, burnin = burnin, thin = thin,
                           type = type, type_num = type_num, mat_type = mat_type, method = method,
                           mean_beta = mean_beta, tau_beta = tau_beta,
                           eta_nu = eta_nu, psi_nu = psi_nu,
                           ninit, maxpoint,
                           var_beta = var_beta, var_eps = var_eps, var_log_mu = var_log_mu, var_log_nu = var_log_nu, var_rho = var_rho,
                           fix_rho_s = fix_rho_s, fix_rho_t = fix_rho_t, fix_rho_st = fix_rho_st,
                           range_rho_s = range_rho_s, range_rho_t = range_rho_t, range_rho_st = range_rho_st,
                           verbose = verbose,
                           c_beta = c_beta, c_eps = c_eps, c_mu = c_mu, c_nu = c_nu, c_rho = c_rho)
  } else{
    stop("It is still not ready")

    if(length(n) == 1) n <- rep(n,length(y))
    else if(length(n) != length(y)) stop("n and y not compatible")
    if(type == "beta-logit") type_num <- 0
    else if(type == "beta-probit") type_num <- 1
    else if(type == "beta-alpha") type_num <- 2
    else if(type == "beta-beta") type_num <- 3
    else stop("Family for the TGMRF not defined\n")

    # out <- mcmc_binary_st(y = y, X = X, n = n,
    #                       n_reg = n_reg, n_var = n_var, neigh = neigh,
    #                       beta = beta, nu = nu, eps = eps,
    #                       rho_s = rho_s, rho_t = rho_t, rho_st = rho_st, tau = tau,
    #                       nsim = nsim, burnin = burnin, thin = thin, nsamp = nsamp,
    #                       type = type, type_num = type_num, mat_type = mat_type, method = method,
    #                       range_nu = range_nu,
    #                       mean_beta = mean_beta, tau_beta = tau_beta,
    #                       eta_nu = eta_nu, psi_nu = psi_nu,
    #                       ninit, maxpoint,
    #                       var_beta = var_beta, var_eps = var_eps, var_log_nu = var_log_nu, var_rho = var_rho,
    #                       verbose = verbose)

  }
  return(out)
}
