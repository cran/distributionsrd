#' @importFrom flexmix stepFlexmix FLXMCdist1 FLXMRglm
#' @importFrom methods new
#' @importFrom modeltools parameters
#' @import stats
#' @importFrom Rdpack reprompt

utils::globalVariables(c(
  "dist1", "dist2", "dist3", "coeff1", "coeff2", "coeff3",
  "alpha1", "alpha2", "c1", "c2", "n", "truncation", "lowertrunc", "uppertrunc"
))
