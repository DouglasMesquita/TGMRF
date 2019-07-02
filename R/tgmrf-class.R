#' An object of 'tgmrf' class
#'
#' @slot eps A matrix with errors samples in each region
#' @slot mu A matrix with means samples in each region
#' @slot beta A matrix with coefficients samples
#' @slot nu A vector with dispersion parameter sample
#' @slot rho A vector with spatial parameter sample
#' @slot model_info A list with other informations about the model
#' @slot fit_measures A data.frame with fit measures
#' @slot time_elapsed A numeric with time elapsed running
#' @slot call A call that result in this object
#'
#' @export

tgmrf_class <- setClass(Class = 'tgmrf',
                        representation(eps = 'matrix',
                                       mu = 'matrix',
                                       beta = 'matrix',
                                       nu = 'numeric',
                                       rho = 'numeric',
                                       model_info = 'list',
                                       fit_measures = 'numeric',
                                       time_elapsed = 'numeric',
                                       call = 'call'))
