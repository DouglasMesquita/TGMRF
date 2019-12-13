#' @title MCMC poisson for spatio-temporal data
#'
#' @description MCMC for Poisson data using TGMRF
#'
#' @param y Response variable
#' @param X Matrix or data.frame of covariates
#' @param E Offsets for Poisson data
#' @param n_reg Number of regions in study
#' @param n_var Number of variables in study
#' @param neigh Neighborhood structure of nb class or a neighborhood matrix
#' @param beta Vector of the initial values to coeficients vector
#' @param nu A Initial value to variance of copula
#' @param eps Initial values to eps parameters
#' @param rho_s A Initial value to spatial dependence
#' @param rho_t A Initial value to temporal dependence
#' @param rho_st A Initial value to spatio-temporal dependence
#' @param tau Precision parameter vector of MCAR
#' @param nsim Number of MCMC iterations
#' @param burnin number of discards iterations
#' @param thin Lag to collect the observations
#' @param type 'lognormal', 'lognormal-precision', 'gamma-shape', 'gamma-scale', 'weibull-shape', 'weibull-scale'
#' @param type_num A numeric that represent type in .C functions
#' @param mat_type car or leroux
#' @param mean_beta Prior mean of beta vector
#' @param tau_beta Prior precision of beta vector
#' @param eta_nu Prior Shape of nu parameter
#' @param psi_nu Prior Scale of nu parameter
#' @param method arms or metropolis
#' @param ninit Number of initial points in ARMS method
#' @param maxpoint Maximum number of evaluation a envelope in each iteration
#' @param var_beta Variance of beta proposal (metropolis)
#' @param var_eps Variance of eps proposal (metropolis)
#' @param var_log_mu Variance of log(mu) proposal (metropolis)
#' @param var_rho Variance of rho proposal (metropolis)
#' @param fix_rho_s Is rho_s fixed?
#' @param fix_rho_t Is rho_t fixed?
#' @param fix_rho_st Is rho_st fixed?
#' @param range_rho_s Range to sample rho_s
#' @param range_rho_t Range to sample rho_t
#' @param range_rho_st Range to sample rho_st
#'
#' @return eps A matrix with eps samples
#' @return mu A matrix with intensities samples
#' @return beta A matrix with beta samples
#' @return nu A vector with samples of nu
#' @return rho A vector with samples of rho
#' @return model_info A data_frame with some models' informations
#' @return fit measures A data_frame with some fit measures

mcmc_poisson_st <- function(y, X, E,
                            n_reg, n_var, neigh,
                            beta, nu, eps, mu,
                            rho_s, rho_t, rho_st,
                            tau,
                            nsim, burnin, thin,
                            type, type_num, mat_type, method,
                            tau_beta, mean_beta,
                            eta_nu, psi_nu,
                            ninit, maxpoint,
                            var_beta, var_eps, var_log_mu, var_log_nu, var_rho,
                            fix_rho_s, fix_rho_t, fix_rho_st,
                            range_rho_s, range_rho_t, range_rho_st,
                            verbose,
                            c_beta, c_eps, c_mu, c_nu, c_rho,
                            conj_beta){

  N <- length(eps)
  P <- length(beta)

  Wt <- abs(outer(1:n_var, 1:n_var, "-")) == 1

  if(class(neigh) == 'nb'){
    Ws <- nb2mat(neigh, style = "B")
  } else{
    Ws <- neigh
  }

  ind_n <- diag(nrow = n_reg, ncol = n_reg)
  ind_n_var <- diag(nrow = n_var, ncol = n_var)

  M <- rowSums(kronecker(ind_n, Wt) + kronecker(Ws, Wt) + kronecker(Ws, ind_n_var))

  if(method == 'arms'){
    method_num = 0
  } else{
    method_num = 1
  }

  other <- list(nsim = nsim, nthin = thin, nsample = nsim,
                tau_beta = tau_beta,
                eta_nu = eta_nu, psi_nu = psi_nu)

  foo <- poimcar_cpp(nsim = nsim, burnin = burnin, thin = thin,
                     eps = eps, mu = mu, beta = beta, nu = nu,
                     rho_s = rho_s, rho_t = rho_t, rho_st = rho_st,
                     X = X, y = y, E = E,
                     Ws = Ws, Wt = Wt,
                     N = N, P = P,
                     mean_beta = mean_beta, tau_beta = tau_beta,
                     eta_nu = eta_nu, psi_nu = psi_nu,
                     fix_rho_s = fix_rho_s, fix_rho_t = fix_rho_t, fix_rho_st = fix_rho_st,
                     range_rho_s = range_rho_s, range_rho_t = range_rho_t, range_rho_st = range_rho_st,
                     type = type_num,
                     var_beta_met = var_beta, var_eps_met = var_eps, var_log_mu_met = var_log_mu, var_log_nu_met = var_log_nu, var_rho_met = var_rho,
                     verbose = verbose,
                     c_beta = c_beta, c_eps = c_eps, c_mu = c_mu, c_nu = c_nu, c_rho = c_rho,
                     conj_beta = conj_beta)

  out <- list(eps = foo$eps,
              mu = foo$mu,
              beta = foo$beta,
              nu = foo$nu,
              rho = foo$rho,
              model_info = other)

  if(is.matrix(out$beta)){
    colnames(out$beta) <- colnames(X)
  }
  colnames(out$rho) <- c("rho_s", "rho_t", "rho_st")

  out$fit_measures <- MEASURES(y = y, lambda = out$mu)

  class(out) <- 'tgmrf'
  return(out)
}
