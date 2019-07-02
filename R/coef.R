#' Coefficients for tgmrf_model class
#'
#' @param object tgmrf_model object to restore the coefficients
#'
#' @method coef tgmrf
#'
#' @export

coef.tgmrf <- function(object){
  coefficients <- colMeans(object$beta)
  return(coefficients)
}
