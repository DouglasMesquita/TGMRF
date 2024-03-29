% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{poimcar_cpp}
\alias{poimcar_cpp}
\title{POIMCAR}
\usage{
poimcar_cpp(nsim, burnin, thin, eps, mu, beta, nu, rho_s, rho_t, rho_st, X,
  y, E, Ws, Wt, N, P, mean_beta, tau_beta, eta_nu, psi_nu, fix_rho_s,
  fix_rho_t, fix_rho_st, range_rho_s, range_rho_t, range_rho_st, type,
  var_beta_met, var_eps_met, var_log_mu_met, var_rho_met, var_log_nu_met,
  verbose, c_beta, c_eps, c_mu, c_nu, c_rho, conj_beta)
}
\arguments{
\item{nsim}{MCMC size}

\item{X}{Covariate matrix}

\item{y}{Response}

\item{E}{Offset}

\item{N}{Dimension of the observations}

\item{P}{Dimension of the covariates}

\item{mean_beta}{Mean a priori to beta vector}

\item{tau_beta}{Variance a priori to beta vector}

\item{eta_nu}{Shape a priori to nu}

\item{psi_nu}{Rate a priori to nu}

\item{type}{TGMRF type (1 to 6)}

\item{var_beta_met}{Variance of beta proposal}

\item{var_eps_met}{Variance of eps proposal}

\item{var_rho_met}{Variance of rho proposal}

\item{var_log_nu_met}{Variance of log(nu) proposal}

\item{M}{Number of neighbors in each area}

\item{W}{Matrix with the neighborhood structure}

\item{rangeRho}{Range to sample rho by using ARMS}

\item{method}{ARMS (0) or Metropolis (1)}

\item{ninit}{Number of initial points in ARMS}

\item{maxpoint}{Maximum number of evaluation in each ARMS iteration}

\item{tau}{Vector of tau parameters to construct Q}
}
\description{
Multivariate Poisson regression with CAR covariance structure
}
