#' Vuong's closeness test
#'
#' Likelihood ratio test for model selection using the Kullback-Leibler information criterion \insertCite{vuong1989likelihood}{distributionsrd}
#'
#' @param x,y vector of log-likelihoods
#' @param np.x,np.y Number of paremeters respectively
#' @param corr type of correction for parameters, defaults to none.
#'
#' @return returns data frame with test statistic, p-value and character vector indicating the test outcome.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' x <- rlnorm(1e4, meanlog = -0.5, sdlog = 0.5)
#' pareto_fit <- combdist.mle(x = x, dist = "pareto")
#' pareto_loglike <- dcombdist(x = x, dist = "pareto", coeff = pareto_fit$coefficients, log = TRUE)
#' lnorm_fit <- combdist.mle(x = x, dist = "lnorm")
#' lnorm_loglike <- dcombdist(x = x, dist = "lnorm", coeff = lnorm_fit$coefficients, log = TRUE)
#'
#' llr_vuong(x = pareto_loglike, y = lnorm_loglike, np.x = pareto_fit$np, np.y = lnorm_fit$np)
#'
#' # BIC type parameter correction
#' llr_vuong(x = pareto_loglike, y = lnorm_loglike, np.x = pareto_fit$np, np.y = lnorm_fit$np,
#' corr = "BIC")
#'
#' # AIC type parameter correction
#' llr_vuong(x = pareto_loglike, y = lnorm_loglike, np.x = pareto_fit$np, np.y = lnorm_fit$np,
#' corr = "AIC")
#' @export

llr_vuong <- function(x, y, np.x, np.y, corr = c("none", "BIC", "AIC")) {
  corr <- match.arg(corr)

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) != length(y)) {
    stop("Lenght of vectors should equal")
  }
  n <- length(x)

  if (corr == "BIC") {
    K <- (np.x - np.y) * log(n) / 2
  } else if (corr == "AIC") {
    K <- np.x - np.y
  } else {
    K <- 0
  }

  lr <- sum(x - y) - K
  omega.2 <- c((n - 1) / n * var(x - y))
  vuongval <- (1 / sqrt(n)) * lr / sqrt(omega.2)
  pval <- pnorm(vuongval, lower.tail = FALSE)
  pval.comp <- pnorm(vuongval, lower.tail = TRUE)
  decision <- ifelse((pval <= 0.05), "x",
    ifelse((pval.comp <= 0.05), "y", "I")
  )

  return_list <- data.frame(testval = vuongval, pval = pval, decision = decision)
}
