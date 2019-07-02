#' @title Conditional Predictive Ordinate for Poisson data
#'
#' @param y Vector of response variable
#' @param lambda Vector of intensities for each region where y was avaliable
#'
#' @return cpo Conditional Predictive Ordinate
#'

CPO_Poisson <- function(y, lambda){
  niter <- length(lambda)
  px <- dpois(x = y, lambda = lambda)
  px <- ifelse(px < .Machine$double.eps, .Machine$double.eps, px)
  cpo <- niter/sum(1/px, na.rm = T)
  return(cpo)
}

#' @title Log Pseudo Marginal Likelihood
#'
#' @param y Vector of response variable
#' @param lambda Matrix of intensities for each region where y was avaliable in each MCMC iteration
#'
#' @return LPML -2 * Log Pseudo Marginal Likelihood measure

LPML <- function(y, lambda){
  aux <- data.frame(y, t(lambda))
  CPO <- apply(X = aux, MARGIN = 1,
               FUN = function(x) CPO_Poisson(y = x[1], lambda = x[-1]))
  LPML <- sum(log(CPO), na.rm = T)
  return(-2*LPML)
}

#' @title Deviance Information Criterion
#'
#' @param y Vector of response variable
#' @param lambda Matrix of intensities for each region where y was avaliable in each MCMC iteration
#'
#' @return DIC Deviance Information Criterion measure

DIC <- function(y, lambda){
  niter <- nrow(lambda)
  thetaBayes <- colMeans(lambda)
  log_px <- dpois(x = y, lambda = thetaBayes, log = TRUE)
  log_px <- sum(log_px, na.rm = T)
  lambda <- cbind(y, t(lambda))
  esp_log_px <- apply(lambda, MARGIN = 1,
                      FUN = function(x) dpois(x = x[1], lambda = x[-1], log = TRUE))
  esp_log_px <- sum(esp_log_px, na.rm = T)/niter

  pDIC <- 2*(log_px - esp_log_px)
  DIC <- -2*log_px + 2*pDIC
  return(DIC)
}

#' @title Widely Applicable Information Criterion
#'
#' @param y Vector of response variable
#' @param lambda Matrix of intensities for each region where y was avaliable in each MCMC iteration
#'
#' @return WAIC_1 Widely Applicable Information Criterion measure

WAIC <- function(y, lambda){
  aux <- data.frame(y, t(lambda))
  px <- apply(aux, MARGIN = 1,
              FUN = function(x) dpois(x = x[1], lambda = x[-1]))
  px <- ifelse(px < .Machine$double.eps, .Machine$double.eps, px)
  log_px <- log(px)
  lppd <- sum(log(colMeans(px, na.rm = T)))
  mean_log_px <- colMeans(log_px, na.rm = T)
  log_mean_px <- log(colMeans(px, na.rm = T))
  pWAIC <- 2*sum(log_mean_px - mean_log_px)
  WAIC <- -2*(lppd - pWAIC)
  return(WAIC)
}

#' @title Pseudo R square
#'
#' @param y Vector of response variable
#' @param lambda Matrix of intensities for each region where y was avaliable in each MCMC iteration
#'
#' @return R2 Pseudo R2 based on raw residuals or deviance residuals

PSEUDO_R2 <- function(y, lambda, residual = "raw"){
  lambda_hat <- colMeans(lambda)

  raw_res <- y - lambda_hat
  y_bar <- mean(y)

  if(residual == "deviance"){
    r2_num <- y*log(y/lambda_hat)
    r2_den <- y*log(y/y_bar)

    r2 <- 1 - sum(r2_num, na.rm = T)/sum(r2_den, na.rm = T)
  } else{
    var_y <- sum((y - y_bar)^2)
    r2 <- 1 - sum(raw_res^2)/var_y
  }

  # beta_hat <- colMeans(beta)
  # eps_hat <- colMeans(eps)
  # lambda_hat <- colMeans(lambda)
  #
  # sig_f <- var(apply(X, MARGIN = 1, function(x) beta_hat%*%x))
  # sig_eps <- var(eps_hat)
  # sig_res <- var(y - lambda_hat)
  #
  # r2_marg <- sig_f/(sig_f + sig_eps + sig_res)
  # r2_cond <- (sig_f + sig_eps)/(sig_f + sig_eps + sig_res)

  return(r2)
}

#' @title LPML, DIC and WAIC measures
#'
#' @param y Vector of response variable
#' @param lambda Matrix of intensities for each region where y was avaliable in each MCMC iteration
#'
#' @return measures A data.frame with LMPL, DIC and WAIC measures

MEASURES <- function(y, lambda){
  lpml <- LPML(y, lambda)
  dic <- DIC(y, lambda)
  waic <- WAIC(y, lambda)
  #r2 <- PSEUDO_R2(y, lambda)
  #measures <- data.frame('DIC' = dic, '-2*LPML' = lpml, 'WAIC' = waic, 'R-squared' = r2, check.names = F)
  measures <- data.frame('DIC' = dic, '-2*LPML' = lpml, 'WAIC' = waic, check.names = F)
  return(measures = measures)
}

#' @title Ergodic mean
#'
#' @description Return the ergodic mean of a vector
#'
#' @param x Vector to restore the ergodic mean
#'
#' @return erg A vector with the ergodic mean

erg_mean <- function(x){
  erg <- cumsum(x)/(1:length(x))
  return(erg)
}

#' @title Plot MCMC
#'
#' @description Plot MCMC chains in 3D
#'
#' @param chain Vector to restore the ergodic mean
#'
#' @return erg A vector with the ergodic mean
#'
plotMCMC <- function(chain, xlab = 'x', ylab = 'Iteration', zlab = 'Density',
                     alpha = 0.5, col = grey(0.2),
                     bty = 'b2', phi = 20, theta = 50,
                     ticktype = "detailed",
                     lwd_density = 3, lwd_chain = 1, ...){

  #--- Graphical parameters
  n <- length(chain)
  alpha.lines = 0.5
  expand = 0.4
  par(mar=c(2, 2, 2, 2))
  ytext <-  -0.15*n

  zPlot <- density(chain)

  box3D(x0 = min(zPlot$"x"), y0 = 0, z0 = 0,
        x1 = max(zPlot$"x"), y1 = n, z1 = 1.1*max(zPlot$"y"),
        bty = "b2",
        xlab = xlab, ylab = ylab, zlab = zlab,
        expand = expand,
        phi = phi, theta = theta, alpha = 0,
        ticktype = ticktype,
        ...)

  #--- Density
  scatter3D(x = zPlot$"x", y = rep(n, length(zPlot$x)), z = zPlot$"y", type = "l",
            add = T, ylim = c(1, n), col = col, lwd = lwd_density, alpha = alpha.lines)
  #--- Chain
  scatter3D(x = chain, y = 1:n, z = rep(0, n), type = "l",
            col = col, lwd = lwd_chain, alpha = alpha, add = T)

  par(mar=c(5, 4, 4, 2) + 0.1)

}

#' @title Build a dependence matrix
#'
#' @description Build a multivariate dependence structure
#'
#' @param S A neighborhood matrix
#' @param tau Vector of precisions of variables
#' @param rho_s Spatial dependence parameter
#' @param rho_t Temporal dependence parameter
#' @param rho_st Spatio-temporal dependence parameter
#'
#' @return Q a dependence structure

buildQ <- function(Ws, Wt, tau, rho_s = 0.99, rho_t = 0.8, rho_st = 0){
  n <- nrow(Ws)
  n_var <- nrow(Wt)

  ind_n <- diag(nrow = n, ncol = n)
  ind_nvar <- diag(nrow = n_var, ncol = n_var)

  out_diag <- rho_s*kronecker(Ws, ind_nvar) + rho_t*kronecker(ind_n, Wt) + rho_st*kronecker(Ws, Wt)
  if(any(c(rho_s, rho_t, rho_st) != 0)){
    D <- diag(colSums(abs(out_diag) > 0))
  } else{
    D <- diag(colSums(kronecker(ind_n, Wt) + kronecker(Ws, Wt) + kronecker(Ws, ind_nvar)))
  }

  Q <- (D - out_diag)*(1/tau)
  return(Q)
}

#' @title Maximum range dependence parameters
#'
#' @description Find the maximum range of the dependence parameters (comparing with the parameters passed)
#'
#' @param S A neighborhood matrix
#' @param tau Vector of precisions of variables
#' @param rho_s Spatial dependence parameter
#' @param rho_t Temporal dependence parameter
#' @param rho_st Spatio-temporal dependence parameter
#'
#' @return Q a dependence structure

max_range <- function(Ws, Wt, rho_s = NULL, rho_t = NULL, rho_st = NULL){
  n <- nrow(Ws)
  n_vars <- nrow(Wt)

  ind_n <- diag(nrow = n, ncol = n)
  ind_nvar <- diag(nrow = n_vars, ncol = n_vars)

  kron_s <- kronecker(Ws, ind_nvar)
  kron_t <- kronecker(ind_n, Wt)
  kron_st <- kronecker(Ws, Wt)

  if(is.null(rho_s)){rho_s <- NA; rho_s_use <- 0} else{rho_s_use <- rho_s}
  if(is.null(rho_t)){rho_t <- NA; rho_t_use <- 0} else{rho_t_use <- rho_t}
  if(is.null(rho_st)){rho_st <- NA; rho_st_use <- 0} else{rho_st_use <- rho_st}

  Q_max <- kron_s + kron_t + kron_st
  D <- rowSums(Q_max)

  D_s <- rowSums(kron_s)
  D_t <- rowSums(kron_t)
  D_st <- rowSums(kron_st)

  D_s_busy <- rowSums(kron_t)*abs(rho_t_use) + rowSums(kron_st)*abs(rho_st_use)
  D_t_busy <- rowSums(kron_s)*abs(rho_s_use) + rowSums(kron_st)*abs(rho_st_use)
  D_st_busy <- rowSums(kron_s)*abs(rho_s_use) + rowSums(kron_t)*abs(rho_t_use)

  max_s <- min(abs(D-D_s_busy)/D_s)
  if(abs(rho_s_use) > max_s) stop(sprintf("rho_s must be smaller than %s", max_s))
  max_t <- min(abs(D-D_t_busy)/D_t)
  if(abs(rho_t_use) > max_t) stop(sprintf("rho_t must be smaller than %s", max_t))
  max_st <- min(abs(D-D_st_busy)/D_st)
  if(abs(rho_st_use) > max_st) stop(sprintf("rho_st must be smaller than %s", max_st))

  df_out <- data.frame(rho = c("rho_s", "rho_t", "rho_st"),
                       value = c(rho_s, rho_t, rho_st),
                       max_value = c(max_s, max_t, max_st),
                       stringsAsFactors = FALSE)

  df_out <- subset(df_out, is.na(value))
  df_out$value <- NULL

  return(df_out)
}

#' @title Append two lists
#'
#' @description Get commom parameters in two list and generate one append list
#'
#' @param x List base
#' @param y Second list
#'
#' @return x
#'
appendList <- function (x, y){
  xnames <- names(x)
  for (v in names(y)) {
    if(v %in% xnames && is.list(x[[v]]) && is.list(y[[v]])){
      x[[v]] <- appendList(x[[v]], y[[v]])
    } else{
      if(!is.null(y[[v]])){
        x[[v]] <- y[[v]]
      }
    }
  }
  return(x)
}
