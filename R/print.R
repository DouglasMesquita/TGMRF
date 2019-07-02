#' Print for tgmrf_model class
#'
#' @param object tgmrf_model object to print
#'
#' @method print tgmrf
#'
#' @export

print.tgmrf <- function(object){
  cat('Transformed Gaussian Markov Random Field: \n')
  cat('\n')
  cat('Call:')
  cat('\n')
  print(object$call)
  cat('\n')
  cat('Coefficients:')
  cat('\n')
  print(coef(object))
  cat('\n')
  cat('Fit measures: \n')
  print(object$fit_measures)
}
