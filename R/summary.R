#' Summary for tgmrf_model class (coda is necessary)
#'
#' @param object tgmrf_model object to summarise
#'
#' @method summary tgmrf
#'
#' @export

summary.tgmrf <- function(object, bw = 0.1, HPD = T){

  ans <- list()
  ans$'Call' <- object$call

  object$beta <- as.matrix(object$beta)

  if(HPD){
    hpd_interval <- apply(X = object$beta, MARGIN = 2, function(x) coda::HPDinterval(coda::as.mcmc(x)))

    coef_info <- data.frame('mean' = colMeans(object$beta),
                            'median' = apply(X = object$beta, MARGIN = 2, median),
                            'mode' = apply(X = object$beta, MARGIN = 2, mode, bw = bw),
                            'std_error' = apply(X = object$beta, MARGIN = 2, sd),
                            'lower_95' = hpd_interval[1,],
                            'upper_95' = hpd_interval[2,])
  } else{
    coef_info <- data.frame('mean' = colMeans(object$beta),
                            'median' = apply(X = object$beta, MARGIN = 2, median),
                            'mode' = apply(X = object$beta, MARGIN = 2, mode, bw = bw),
                            'std_error' = apply(X = object$beta, MARGIN = 2, sd),
                            'lower_95' = apply(X = object$beta, MARGIN = 2, function(x) quantile(x = x, probs = 0.025)),
                            'upper_95' = apply(X = object$beta, MARGIN = 2, function(x) quantile(x = x, probs = 0.975)))
  }

  ans$'Coeficients' <- coef_info

  if(!is.null(object$call$group_var)){
    if(HPD){
      hpd_rho <- apply(X = object$rho, MARGIN = 2, function(x) coda::HPDinterval(coda::as.mcmc(x)))
      hpd_nu <- coda::HPDinterval(coda::as.mcmc(object$nu))

      other_info <- data.frame('mean' = c(colMeans(object$rho), mean(object$nu)),
                               'median' = c(apply(X = object$rho, MARGIN = 2, median), median(object$nu)),
                               'mode' = c(apply(X = object$rho, MARGIN = 2, mode, bw = bw), mode(object$nu)),
                               'std_error' = c(apply(X = object$rho, MARGIN = 2, sd), sd(object$nu)),
                               'lower_95' = c(hpd_rho[1,],
                                              hpd_nu[1]),
                               'upper_95' = c(hpd_rho[2, ],
                                              hpd_nu[2]))
    } else{
      other_info <- data.frame('mean' = c(colMeans(object$rho), mean(object$nu)),
                               'median' = c(apply(X = object$rho, MARGIN = 2, median), median(object$nu)),
                               'mode' = c(apply(X = object$rho, MARGIN = 2, mode, bw = bw), mode(object$nu)),
                               'std_error' = c(apply(X = object$rho, MARGIN = 2, sd), sd(object$nu)),
                               'lower_95' = c(apply(X = object$rho, MARGIN = 2, function(x) quantile(x = x, probs = 0.025)),
                                              quantile(x = object$nu, probs = 0.025)),
                               'upper_95' = c(apply(X = object$rho, MARGIN = 2, function(x) quantile(x = x, probs = 0.975)),
                                              quantile(x = object$nu, probs = 0.975)))
    }
    rownames(other_info) <- c("Spatial parameter", "Temporal parameter", "Spatio-temporal parameter", "Dispersion parameter")
    ans$'Other parameters' <- other_info
  } else{
    if(HPD){
      hpd_rho <- coda::HPDinterval(coda::as.mcmc(object$rho))
      hpd_nu <- coda::HPDinterval(coda::as.mcmc(object$nu))

      other_info <- data.frame('mean' = c(mean(object$rho), mean(object$nu)),
                               'median' = c(median(object$rho), median(object$nu)),
                               'mode' = c(mode(object$rho, bw = bw), mode(object$nu, bw = bw)),
                               'std_error' = c(sd(object$rho), sd(object$nu)),
                               'lower_95' = c(hpd_rho[1],
                                              hpd_nu[1]),
                               'upper_95' = c(hpd_rho[2],
                                              hpd_nu[2]))
    } else{
      other_info <- data.frame('mean' = c(mean(object$rho), mean(object$nu)),
                               'median' = c(median(object$rho), median(object$nu)),
                               'mode' = c(mode(object$rho, bw = bw), mode(object$nu, bw = bw)),
                               'std_error' = c(sd(object$rho), sd(object$nu)),
                               'lower_95' = c(quantile(x = object$rho, probs = 0.025),
                                              quantile(x = object$nu, probs = 0.025)),
                               'upper_95' = c(quantile(x = object$rho, probs = 0.975),
                                              quantile(x = object$nu, probs = 0.975)))
    }
    rownames(other_info) <- c("Spatial parameter", "Dispersion parameter")
    ans$'Other parameters' <- other_info
  }

  ans$'Fit measures' <- object$fit_measures
  ans
}
