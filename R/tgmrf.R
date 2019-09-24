#' @title Transformed Gaussian Markov Random Field
#'
#' @useDynLib TGMRF
#'
#' @description Use this function to fit spatial data using non-guassian copulas.
#'
#' @param data A data.frame with all variables used to fit a model
#' @param formula Formula to fit the model
#' @param neigh Neighborhood structure of nb class or a neighborhood matrix
#' @param scale Logical parameter to scale covariates
#' @param spatial_var Variable name that indicates the column which represents the spatial variable
#' @param neigh_order Order that the spatial_var appear on data
#' @param group_var Variable name that indicates the column which represents the variable variable
#' @param change_vars Variables that change in the time. We'll return n_var parameters for this variable. One for each n_var.
#' @param beta Vector of the initial values to coeficients vector
#' @param nu A Initial value to variance of copula
#' @param rho_s Spatial dependence parameter
#' @param rho_t Temporal dependence parameter
#' @param rho_st Spatio-temporal dependence parameter
#' @param family 'Poisson' or 'Binary'
#' @param type Depends of family.
#' 'lognormal', 'lognormal-precision', 'gamma-shape', 'gamma-scale', 'weibull-shape', 'weibull-scale' for Poisson family
#' 'beta-logit', 'beta-probit', 'beta-alpha', 'beta-beta' for Binomial family
#' @param mat_type car or leroux
#' @param method metropolis or arms
#' @param nsim Number of MCMC iterations
#' @param burnin number of discards iterations
#' @param thin Lag to collect the observations
#' @param E Offsets for Poisson data
#' @param n Number of trials in Binomial data
#' @param prior_param List of priors parameters. Gaussian for beta vector and Gamma for nu. Entrys are: 'nu' (list with shape and scale) and 'beta' (list with mean and precision)
#' @param MCMC_config Parameters to tunning MCMC. Entrys are: 'arms' (list with ninit and maxpoint) and 'metropolis' (list with var_beta, var_eps, var_rho and var_log_nu)
#' @param fix_rho A list informing if any rho (rho_s, rho_t, rho_st) is fixed
#' @param range A list informing a range for sampling rho (rho_s, rho_t, rho_st). Each entre must be a vector in the interval (-1, 1).
#'
#' @return beta A matrix with samples of beta
#' @return nu A vector with samples of nu
#' @return rho A vector with samples of rho
#'
#' @examples ##-- Pacotes ----
#' library(TGMRF)
#' library(tmvtnorm)
#' library(truncnorm)
#' library(OpenMPController)
#'
#' #'#'-- Configurações ----
#' rowid <- 8
#' colid <- 5
#' n_vars <- 2
#' N <- rowid*colid*n_vars
#'
#' betas <- 0.5
#' intercept <- -.5
#' nu <- 0.3
#'
#' P <- length(betas)
#' X <- scale(matrix(c(rnorm(n = rowid*colid*n_vars, mean = 0)), ncol = P))
#'
#' type_data <- 'lognormal'
#' family <- 'poisson'
#' seed <- 123456
#'
#' rho_s <- 0.95*2.25
#' rho_t <- 0.00001
#' rho_st <- 0.00001
#'
#' #'#'-- Dados ----
#' data_poisson <- rtgmrf(rowid = rowid, colid = colid, X = X, n_var = n_vars,
#'                        rho_s = rho_s, rho_t = rho_t, rho_st = rho_st,
#'                        betas = betas, nu = nu, intercept = intercept,
#'                        type_data = type_data, family = family,
#'                        seed = seed)
#'
#' hist(data_poisson$y, breaks = 20)
#'
#' bd <- data.frame(y = data_poisson$y,
#'                  regiao = data_poisson$reg, grupo = data_poisson$var,
#'                  X1 = data_poisson$X,
#'                  eps = data_poisson$eps,
#'                  check.names = F)
#' neigh <- data_poisson$neigh
#'
#' #'#'-- Configurações ----
#' nsim <- 20000
#' burnin <- 0
#' thin <- 1
#' formula <- paste0('y ~ X1')
#'
#' #'#'-- Metropolis rcpp ----
#' vars_beta <- diag(c(0.03, 0.02), nrow = P+1, ncol = P+1)
#' vars_eps <- diag(0.0001, nrow = N, ncol = N)
#' var_log_nu <- 0.85
#' var_rho <- diag(c(0.1, 0.3, 0.05), ncol = 3, nrow = 3)
#'
#' omp_set_num_threads(n = 1)
#' type_model <- type_data
#'
#' c_beta <- 1#'(2.38^2)/P
#' c_eps <- 1#'(2.38^2)/N
#' c_nu <- 1#'(2.38^2)
#' c_rho <- 1#'(2.38^2)/3
#'
#' system.time(
#'   out_metcpp <- tgmrf(data = bd, formula = formula,
#'                       spatial_var = "regiao", group_var = "grupo",
#'                       beta = c(intercept, betas), eps = data_poisson$eps,
#'                       nu = nu,
#'                       rho_s = 0.95, rho_t = 0, rho_st = 0,
#'                       family = family, type = type_model, mat_type = "car", method = "metropolis",
#'                       nsim = nsim, burnin = burnin, thin = thin,
#'                       prior_param = list("nu" = list("shape" = 0.2, "rate" = 0.2)),
#'                       E = NULL, neigh = neigh,
#'                       fix_rho = list("rho_s" = FALSE, "rho_t" = FALSE, "rho_st" = FALSE),
#'                       MCMC_config = list('metropolis' = list('var_beta' = vars_beta,
#'                                                              'var_eps' = vars_eps,
#'                                                              'var_log_nu' = var_log_nu,
#'                                                              'var_rho' = var_rho)),
#'                       range = list("rho_s" = c(-1, 1), "rho_t" = c(-1, 1), "rho_st" = c(-1, 1)),
#'                       scale = FALSE, verbose = TRUE,
#'                       c_beta = c_beta, c_eps = c_eps, c_nu = c_nu, c_rho = c_rho)
#' )
#'
#' summary(out_metcpp)
#'
#' @import RcppArmadillo
#' @import tmvtnorm
#' @importFrom spdep nb2mat
#'
#' @export

tgmrf <- function(data, formula, neigh, scale = T,
                  spatial_var, neigh_order = NULL,
                  group_var = NULL,  change_vars = NULL,
                  beta = NULL, nu = 1, eps = NULL, mu = NULL,
                  rho_s = 0, rho_t = 0, rho_st = 0,
                  tau = NULL,
                  family = "poisson", type = "gamma-shape", mat_type = "car", method = "metropolis",
                  nsim = 1000, burnin = 0, thin = 1,
                  E = NULL, n = 1,
                  prior_param = NULL, MCMC_config = NULL, fix_rho = NULL,
                  range = list("rho_s" = c(-1, 1), "rho_t" = c(-1, 1), "rho_st" = c(-1, 1)),
                  verbose = FALSE,
                  c_beta = NULL, c_eps = NULL, c_mu = NULL, c_nu = NULL, c_rho = NULL){

  ##-- Joining lists
  prior_param_dft <- list("nu" = list(shape = 0.1, rate = 0.1),
                          "beta" = list(mean = 0, precision = 0.01))
  prior_param <- appendList(prior_param_dft, prior_param)

  if(is.null(group_var)){
    MCMC_config_dft <- list('arms' = list('ninit' = 5, 'maxpoint' = 100),
                            'metropolis' = list('var_beta' = 0.1, 'var_eps' = 0.1, 'var_log_mu' = 0.1,
                                                'var_log_nu' = 0.1,
                                                'var_rho' = 0.1))
    MCMC_config <- appendList(MCMC_config_dft, MCMC_config)

    ##--
    fix_rho_dft <- list(rho_s = F)
    fix_rho <- appendList(fix_rho_dft, fix_rho)

    range_dft <- list("rho_s" = c(-1, 1))
    range <- appendList(range_dft, range)
  } else{
    MCMC_config_dft <- list('arms' = list('ninit' = 5, 'maxpoint' = 100),
                            'metropolis' = list('var_beta' = 0.1, 'var_eps' = 0.1, 'var_log_mu' = 0.1,
                                                'var_log_nu' = 0.1,
                                                'var_rho' = c(0.1, 0.1, 0.1)))
    MCMC_config <- appendList(MCMC_config_dft, MCMC_config)

    ##--
    fix_rho_dft <- list(rho_s = F, rho_t = F, rho_st = F)
    fix_rho <- appendList(fix_rho_dft, fix_rho)

    range_dft <- list("rho_s" = c(-1, 1), "rho_t" = c(-1, 1), "rho_st" = c(-1, 1))
    range <- appendList(range_dft, range)
  }

  ##-- Checking parameters
  if(any(abs(unlist(range)) > 1)) stop("range should be a list of range vectors (rho_s, rho_t, rho_st) in the range (-1, 1)")

  if(!is.data.frame(data)) stop("data must be a data.frame")
  if(!(class(neigh) %in% c("nb", "matrix"))) stop("neigh must be a object of class nb or a matrix")

  N <- nrow(data)

  if(!is.null(group_var)){
    n_var <- nrow(unique(data[group_var]))
    n_reg <- N/n_var
    if(is.null(tau)){
      tau <- diag(1, n_var)
    }

    data <- data[order(data[, spatial_var], data[, group_var]),]

  } else{
    data <- data[order(data[, spatial_var]),]
  }

  call_tgmrf <- match.call()

  model_fr <- match.call(expand.dots = FALSE)
  match_strings <- match(c("formula", "data"), names(model_fr), 0L)
  model_fr <- model_fr[c(1L, match_strings)]
  model_fr[[1L]] <- quote(stats::model.frame)
  model_fr <- eval(model_fr, parent.frame())
  model_types <- attr(model_fr, "terms")
  y <- model.response(model_fr, "numeric")
  X <- model.matrix(model_types, model_fr)

  if(scale){
    if("(Intercept)" %in% colnames(X) | "(Intercept)" %in% names(X)){
      pos <- which(colnames(X) == "(Intercept)")
      X <- cbind("(Intercept)" = 1, scale(X[,-pos]))
    } else{
      X <- scale(X)
    }
  }

  if(!is.null(change_vars)){
    n_var_t <- length(change_vars)
    X_t <- X[, change_vars]
    X_t <- kronecker(X_t, matrix(1, ncol = n_var))
    mAux <- kronecker(matrix(1, nrow = n_reg, ncol = n_var_t), diag(1, ncol = n_var, nrow = n_var))
    X_t <- X_t*mAux
    colnames(X_t) <- paste(rep(change_vars, each = n_var), rep(1:n_var, n_var_t), sep = "_")
    pos_var_change <- which(colnames(X) %in% change_vars)
    X <- X[, -pos_var_change]
    X <- cbind(X, X_t)
  }

  P <- ncol(X)

  if(is.null(c_beta)) c_beta <- (2.38^2)/P
  if(is.null(c_eps)) c_eps <- 2.38^2
  if(is.null(c_mu)) c_mu <- 2.38^2
  if(is.null(c_nu)) c_nu <- 2.38^2
  if(is.null(c_rho)) c_rho <- (2.38^2)/3

  if(is.null(beta)) beta <- rep(0, P)
  if(is.null(E)) E <- rep(1, N)
  if(is.null(n)) n <- rep(1, N)

  initial_time <- Sys.time()
  if(is.null(group_var)){
    stop("It is still not ready")

    if((family != "poisson") & (family != "binary")){
      stop("The only families accepted are: poisson or binary\n")
    }
    if(class(neigh) != "nb" & class(neigh) != "matrix"){
      stop("The neighbors musta belong to a nb class\n")
    }
    if(any(!is.numeric(unlist(prior_param)))){
      stop("Prior list must be numeric! \n")
    }
    if(method != 'arms' & method != 'metropolis'){
      stop("The only methods accepted are: arms or metropolis\n")
    }

    # out <- s_tgmrf(y = y, X = X,
    #                beta = beta, nu = nu, rho = rho_s,
    #                family = family, type = type, mat_type = mat_type, method = method,
    #                nsim, burnin, thin, nsamp = nsamp,
    #                E = E, n = n, neigh = neigh,
    #                prior_param = prior_param,
    #                MCMC_config = MCMC_config,
    #                verbose = verbose)

  } else{
    out <- st_tgmrf(y = y, X = X, n_reg = n_reg, n_var = n_var,
                    beta = beta, nu = nu, eps = eps, mu = mu,
                    rho_s = rho_s, rho_t = rho_t, rho_st = rho_st,
                    tau = tau,
                    family = family, type = type, mat_type = mat_type, method = method,
                    nsim = nsim, burnin = burnin, thin = thin,
                    E = E, n = n, neigh = neigh,
                    prior_param = prior_param,
                    MCMC_config = MCMC_config,
                    fix_rho = fix_rho, range = range,
                    verbose = verbose,
                    c_beta = c_beta, c_eps = c_eps, c_mu = c_mu, c_nu = c_nu, c_rho = c_rho)
  }
  final_time <- Sys.time()
  out$time_elapsed <- final_time - initial_time
  out$call <- call_tgmrf

  return(out)
}
