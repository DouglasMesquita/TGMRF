#' @title Random observations of a Transformed Gaussian Markov Random Field
#'
#' @description Use this function to simulate a data of a TGMRF.
#'
#' @param rowid Number of lines of a lattice. (Used if neigh = NULL)
#' @param colid Number of columns of a lattice. (Used if neigh = NULL)
#' @param neigh A object of class nb
#' @param n Vector with size of binomial trials
#' @param rho Dependence parameter
#' @param betas Coeficients
#' @param intercept Should use intercept? NULL or coefficiente
#' @param nu Dispersion parameter for gamma models
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
#'
#' @importFrom spdep cell2nb
#' @importFrom spdep nb2mat
#' @importFrom mvtnorm rmvnorm
#'
#' @export
#'
#' @examples #--- library(mvtnorm)
#' #--- library(spdep)
#'
#' data(nenia)
#' coord <- as.data.frame(cbind(x = nenia$x, y = nenia$y))
#' coordinates(coord) <- c("x","y")
#' dst <- 3
#' neigh <- dnearneigh(coord, 0, dst)
#'
#' data_poisson <- rtgmrf_s(neigh = neigh,
#'                        rho_s = 0.9, betas = c(1, 0.9, -0.1), nu = 2,
#'                        type_data = 'gamma-shape', family = 'poisson',
#'                        seed = 1)
#'
#' data_binary <- rtgmrf_s(rowid = 10, colid = 10,
#'                       rho_s = 0.9, betas = c(1, 0.9), nu = 5,
#'                       type_data = 'gamma-shape', family = 'poisson',
#'                       seed = 1)

rtgmrf <- function(rowid = 10, colid = 10, X = NULL, n_var = 1,
                   neigh = NULL, n_trials = NULL, E = NULL,
                   rho_s = 0.9, rho_t = 0, rho_st = 0,
                   betas = c(-0.1, 0.3, 0.8), intercept = NULL, nu = 2, tau = 1,
                   type_data = 'gamma-shape', family = 'poisson', mat_type = 'car',
                   seed = 1){

  if(!('mvtnorm' %in% loadedNamespaces())){
    stop('The package mvtnorm was not found \n')
  }
  if(!('spdep' %in% loadedNamespaces())){
    stop('The package spdep was not found \n')
  }
  if((family != "poisson") && (family != "binary")){
    stop("The only families accepted are: poisson or binary \n")
  }
  if(!is.null(neigh) & class(neigh) != "nb"){
    stop("The neighbors musta belong to a nb class\n")
  }
  if(family == 'poisson' & !(type_data %in% c('gamma-shape', 'gamma-scale', 'gamma-precision',
                                              'lognormal', 'lognormal-precision', 'weibull-scale', 'weibull-shape'))){
    stop("The only models accepted for poisson family are: 'gamma-shape', 'gamma-scale', 'gamma-precision',
                                              'lognormal', 'lognormal-precision', 'weibull-scale' or 'weibull-shape' \n ")
  }
  if(family == 'binary' & !(type_data %in% c('beta-logit', 'beta-probit', 'beta-alpha', 'beta-beta'))){
    stop("The only models accepted for binary family are: 'beta-logit', 'beta-probit',
         'beta-alpha'or 'beta-beta' \n ")
  }
  if(is.null(n_var)){
    # data <- rtgmrf_s(rowid = rowid, colid = colid, X = X, neigh = neigh, n_trials = n_trials, E = E,
    #                  rho_s = rho_s, betas = betas, intercept = intercept, nu = nu,
    #                  type_data = type_data, family = family, mat_type = mat_type,
    #                  seed = seed)

  }else{
    data <- rtgmrf_st(rowid = rowid, colid = colid, X = X, n_var = n_var, neigh = neigh, n_trials = n_trials, E = E,
                      rho_s = rho_s, rho_t = rho_t, rho_st = rho_st,
                      betas = betas, intercept = intercept, nu = nu,
                      tau = tau,
                      type_data = type_data, family = family, mat_type = mat_type,
                      seed = seed)
  }
  return(data)
}
