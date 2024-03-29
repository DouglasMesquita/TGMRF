% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtgmrf.R
\name{rtgmrf}
\alias{rtgmrf}
\title{Random observations of a Transformed Gaussian Markov Random Field}
\usage{
rtgmrf(rowid = 10, colid = 10, X = NULL, n_var = 1, neigh = NULL,
  n_trials = NULL, E = NULL, rho_s = 0.9, rho_t = 0, rho_st = 0,
  betas = c(-0.1, 0.3, 0.8), intercept = NULL, nu = 2, tau = 1,
  type_data = "gamma-shape", family = "poisson", mat_type = "car",
  seed = 1)
}
\arguments{
\item{rowid}{Number of lines of a lattice. (Used if neigh = NULL)}

\item{colid}{Number of columns of a lattice. (Used if neigh = NULL)}

\item{neigh}{A object of class nb}

\item{betas}{Coeficients}

\item{intercept}{Should use intercept? NULL or coefficiente}

\item{nu}{Dispersion parameter for gamma models}

\item{type_data}{Depends of family.
'lognormal', 'lognormal-precision', 'gamma-shape', 'gamma-scale', 'weibull-shape', 'weibull-scale' for Poisson family
'beta-logit', 'beta-probit', 'beta-alpha', 'beta-beta' for Binomial family}

\item{family}{'poisson' or 'binary'}

\item{seed}{A seed to reproduce the reuslts}

\item{n}{Vector with size of binomial trials}

\item{rho}{Dependence parameter}
}
\value{
y Response variable

X Covariates matrix

neigh Neighborhood structure

Q Covariance matrix

mu Means (TGMRF)

eps Errors (GMRF)
}
\description{
Use this function to simulate a data of a TGMRF.
}
\examples{
#--- library(mvtnorm)
#--- library(spdep)

data(nenia)
coord <- as.data.frame(cbind(x = nenia$x, y = nenia$y))
coordinates(coord) <- c("x","y")
dst <- 3
neigh <- dnearneigh(coord, 0, dst)

data_poisson <- rtgmrf_s(neigh = neigh,
                       rho_s = 0.9, betas = c(1, 0.9, -0.1), nu = 2,
                       type_data = 'gamma-shape', family = 'poisson',
                       seed = 1)

data_binary <- rtgmrf_s(rowid = 10, colid = 10,
                      rho_s = 0.9, betas = c(1, 0.9), nu = 5,
                      type_data = 'gamma-shape', family = 'poisson',
                      seed = 1)
}
