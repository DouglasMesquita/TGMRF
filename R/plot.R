#' Plot for tgmrf_model class
#'
#' @param object tgmrf_model object to plot
#'
#' @method plot tgmrf
#'
#' @export

plot.tgmrf <- function(object, series = T, series3D = F, erg = F, boxplot = F, coord = NULL,
                       xlim = NULL, ylim = NULL,
                       col = '#4682b460', pch = 20,
                       main = 'Spatial effects', bty = "l", ask = T, ...){

  temporal <- ncol(object$rho) == 3

  if(temporal){
    data <- data.frame(object$beta, object$rho, nu = object$nu)
  } else{
    data <- data.frame(object$beta, rho = object$rho, nu = object$nu)
  }

  if(names(data)[1] == 'X.Intercept.'){
    names(data)[1] <- 'Intercept'
  }
  P <- ncol(object$beta)
  N <- nrow(object$beta)

  ##-- To plot the MCMC chains ----
  if(series){
    ##-- Coefficients ----
    for(i in 1:P){
      xlab <- 'Iteration'
      ylab <- colnames(object$beta)[i]
      plot(x = 1:N, y = object$beta[,i],
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = xlab, ylab = ylab,
           bty = bty, ...)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dependence parameters ----
    if(temporal){
      ##-- + Spatial parameter ----
      plot(x = 1:N, y = data$rho_s,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatial parameter',
           bty = bty, ...)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Temporal parameter ----
      plot(x = 1:N, y = data$rho_t,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Temporal parameter',
           bty = bty, ...)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Temporal parameter ----
      plot(x = 1:N, y = data$rho_st,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatio-temporal parameter',
           bty = bty, ...)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    } else{
      ##-- + Spatial parameter ----
      plot(x = 1:N, y = data$rho,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatial parameter',
           bty = bty, ...)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dispersion parameter ----
    plot(x = 1:N, y = data$nu,
         col = 'grey20', type = 'l', cex = 1.5,
         xlab = 'Iteration', ylab = 'Dispersion parameter',
         bty = bty, ...)
  }
  ##-- To plot the MCMC chains in 3D ----
  if(series3D){
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    if('plot3D' %in% installed.packages()){
      library(plot3D)
    } else{
      stop('The package plot3D was not found \n')
    }
    ##-- Coefficients ----
    for(i in 1:P){
      ylab <- 'Iteration'
      zlab <- 'Density'
      xlab <- colnames(object$beta)[i]
      plotMCMC(chain = object$beta[,i], xlab = xlab, ylab = ylab, zlab = zlab)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dependence parameters ----
    if(temporal){
      ##-- + Spatial parameter ----
      plotMCMC(chain = data$rho_s, xlab = 'Spatial parameter', ylab = 'Iteration', zlab = 'Density')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Temporal parameter ----
      plotMCMC(chain = data$rho_t, xlab = 'Temporal parameter', ylab = 'Iteration', zlab = 'Density')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Spatio-temporal parameter ----
      plotMCMC(chain = data$rho_st, xlab = 'Spatio-temporal parameter', ylab = 'Iteration', zlab = 'Density')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    } else{
      ##-- + Spatial parameter ----
      plotMCMC(chain = data$rho, xlab = 'Spatial parameter', ylab = 'Iteration', zlab = 'Density')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dispersion parameter ----
    plotMCMC(chain = data$nu, xlab = 'Dispersion parameter', ylab = 'Iteration', zlab = 'Density')
  }
  ##-- To plot the ergodic mean ----
  if(erg){
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    ##-- Coefficients ----
    for(i in 1:P){
      em <- erg_mean(object$beta[,i])
      xlab <- 'Iteration'
      ylab <- colnames(object$beta)[i]
      mean_point <- mean(object$beta[,i])
      plot(x = 1:N, y = em,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = xlab, ylab = ylab, main = 'Ergodic mean',
           bty = bty, ...)
      abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dependence parameters ----
    if(temporal){
      ##-- + Spatial parameter ----
      em <- erg_mean(data$rho_s)
      mean_point <- mean(data$rho_s)
      plot(x = 1:N, y = em,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatial parameter', main = 'Ergodic mean',
           bty = bty, ...)
      abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Spatial parameter ----
      em <- erg_mean(data$rho_t)
      mean_point <- mean(data$rho_t)
      plot(x = 1:N, y = em,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Temporal parameter', main = 'Ergodic mean',
           bty = bty, ...)
      abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
      ##-- + Spatial parameter ----
      em <- erg_mean(data$rho_st)
      mean_point <- mean(data$rho_st)
      plot(x = 1:N, y = em,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatio-temporal parameter', main = 'Ergodic mean',
           bty = bty, ...)
      abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    } else{
      ##-- + Spatial parameter ----
      em <- erg_mean(data$rho)
      mean_point <- mean(data$rho)
      plot(x = 1:N, y = em,
           col = 'grey20', type = 'l', cex = 1.5,
           xlab = 'Iteration', ylab = 'Spatial parameter', main = 'Ergodic mean',
           bty = bty, ...)
      abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dispersion parameter ----
    em <- erg_mean(data$nu)
    mean_point <- mean(data$nu)
    plot(x = 1:N, y = em,
         col = 'grey20', type = 'l', cex = 1.5,
         xlab = 'Iteration', ylab = 'Dispersion parameter', main = 'Ergodic mean',
         bty = bty, ...)
    abline(h = mean_point, col = 'grey30', lwd = 1.5, lty = 'dashed')
  }
  ##-- To plot the boxplot ----
  if(boxplot){
    ##-- Coefficients ----
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    boxplot(object$beta, pch = 20, bty = bty)
    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    ##-- Dependence parameters ----
    if(temporal){
      ##-- + Spatial parameter ----
      boxplot(data[c("rho_s", "rho_t", "rho_st")], pch = 20, xlab = "Dependence parameters", bty = bty)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    } else{
      ##-- + Spatial parameter ----
      boxplot(data$rho, pch = 20, xlab = "Spatial parameter", bty = bty)
      if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    }
    ##-- Dispersion parameter ----
    boxplot(data$nu, pch = 20, xlab = "Dispersion parameter", bty = bty)
  }
  ##-- To plot the spatial effects in a simple map ----
  if(!is.null(coord)){

    if(ask) invisible(readline(prompt = "Press [enter] to see the next plot..."))
    if(class(coord) == 'SpatialPoints'){
      coord <- coord@coords
    }

    x_coord <- coord[, 1]
    y_coord <- coord[, 2]
    cex_scale <- colMeans(x$mu)*4/max(colMeans(x$mu)) + 1
    cnst <- 0.08

    if(is.null(xlim)){
      xlim <- range(x_coord)
      xlim <- c(xlim[1] - cnst*abs(xlim[2]), xlim[2] + cnst*abs(xlim[2]))
    }
    if(is.null(ylim)){
      ylim <- range(y_coord)
      ylim <- c(ylim[1] - cnst*abs(ylim[2]), ylim[2] + cnst*abs(ylim[2]))
    }

    plot.new()
    plot.window(xlim = xlim, ylim = ylim, asp = 1)

    points(x = x_coord, y = y_coord, xlim = xlim, ylim = ylim,
           cex = cex_scale, pch = pch, col = col, ...)

    polygon(x = c(xlim[1], xlim[2], xlim[2], xlim[1]),
            y = c(ylim[1], ylim[1], ylim[2], ylim[2]),
            col = 'transparent', border = 'grey30')

    title(main = main)
  }
}
