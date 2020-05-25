#' Combined distributions
#'
#' Density, distribution function, quantile function, raw moments and random generation for combined (empirical, single, composite and finite mixture) truncated or complete distributions.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter
#' @param dist character vector denoting the distribution(s).
#' @param prior Numeric vector of prior coefficients, defaults to single vector with value one.
#' @param coeff list of parameters for the distribution(s).
#' @param lowertrunc,uppertrunc lowertrunc- and uppertrunc truncation points, defaults to 0 and Inf respectively
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \leq x]} \eqn{(E[x^r|X \leq y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#' @param compress Logical indicating whether return values from individual densities of finite mixtures should be gathered or not, defaults to TRUE.
#'
#' @return dcombdist gives the density, pcombdist gives the distribution function, qcombdist gives the quantile function, mcombdist gives the rth moment of the distribution and rcombdist generates random deviates.
#'
#'  The length of the result is determined by n for rcombdist, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#'
#' \donttest{
#' # Load necessary tools
#' data("fit_US_cities")
#' library(tidyverse)
#' x <- rcombdist(
#'   n = 25359, dist = "lnorm",
#'   prior = subset(fit_US_cities, (dist == "lnorm" & components == 5))$prior[[1]],
#'   coeff = subset(fit_US_cities, (dist == "lnorm" & components == 5))$coefficients[[1]]
#' ) # Generate data from one of the fitted functions
#'
#' # Evaluate functioning of dcomdist by calculating log likelihood for all distributions
#' loglike <- fit_US_cities %>%
#'   group_by(dist, components, np, n) %>%
#'   do(loglike = sum(dcombdist(dist = .[["dist"]], x = sort(x), prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]], log = TRUE))) %>%
#'   unnest(cols = loglike) %>%
#'   mutate(NLL = -loglike, AIC = 2 * np - 2 * (loglike), BIC = log(n) * np - 2 * (loglike)) %>%
#'   arrange(NLL)
#'
#' # Evaluate functioning of mcombdist and pcombdist by calculating NMAD
#' #(equivalent to the Kolmogorov-Smirnov test statistic for the zeroth moment
#' #of the distribution) for all distributions
#' nmad <- fit_US_cities %>%
#'   group_by(dist, components, np, n) %>%
#'   do(
#'     KS = max(abs(pempirical(q = sort(x), data = x) - pcombdist(dist = .[["dist"]],
#'     q = sort(x), prior = .[["prior"]][[1]], coeff = .[["coefficients"]][[1]]))),
#'     nmad_0 = nmad_test(r = 0, dist = .[["dist"]], x = sort(x), prior = .[["prior"]][[1]],
#'     coeff = .[["coefficients"]][[1]], stat = "max"),
#'     nmad_1 = nmad_test(r = 1, dist = .[["dist"]], x = sort(x), prior = .[["prior"]][[1]],
#'     coeff = .[["coefficients"]][[1]], stat = "max")
#'   ) %>%
#'   unnest(cols = c(KS, nmad_0, nmad_1)) %>%
#'   arrange(nmad_0)
#'
#' # Evaluate functioning of qcombdist pcombdist by calculating NMAD (equivalent to the Kolmogorov-
#' #Smirnov test statistic for the zeroth moment of the distribution) for all distributions
#' test <- fit_US_cities %>%
#'   group_by(dist, components, np, n) %>%
#'   do(out = qcombdist(pcombdist(2, dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]], log.p = TRUE),
#'     dist = .[["dist"]], prior = .[["prior"]][[1]], coeff = .[["coefficients"]][[1]],
#'     log.p = TRUE
#'   )) %>%
#'   unnest(cols = c(out))
#' }
#'
#' @name combdist

#' @rdname combdist
#' @export

dcombdist <- function(x, dist, prior = c(1), coeff, log = FALSE, compress = TRUE, lowertrunc = 0, uppertrunc = Inf) {
  dist <- unlist(strsplit(dist, split = "_"))

  if (length(grep("trunc", dist)) > 0) {

    # Single Truncated distribution
    if (length(prior) == 1) {

      # Single distribution
      d <- dtruncdist(x = x, dist = dist, coeff = coeff, log = log, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    } else {

      # Finite Mixture

      if (log) {
        coeff <- as.matrix(coeff)
        d <- sapply(1:length(prior), function(i, coeff) {
          coeff.int <- coeff[, i]
          names(coeff.int) <- rownames(coeff)
          return((log(prior[i]) + dtruncdist(x = x, coeff = coeff.int, log = TRUE, lowertrunc = lowertrunc, uppertrunc = uppertrunc)))
        }, coeff = coeff)

        colnames(d) <- names(prior)

        if (compress) {
          d <- log(rowSums(exp(d)))
        }
      } else {
        coeff <- as.matrix(coeff)
        d <- sapply(1:length(prior), function(i, coeff, log) {
          coeff.int <- coeff[, i]
          names(coeff.int) <- rownames(coeff)
          return(prior[i] * dtruncdist(x = x, coeff = coeff.int, log = log, lowertrunc = lowertrunc, uppertrunc = uppertrunc))
        }, coeff = coeff, log = log)

        colnames(d) <- names(prior)

        if (compress) {
          d <- rowSums(d)
        }
      }
    }
  } else if (is.null(coeff)) {

    # Empirical distribution
    d <- do.call(paste0("d", dist), c(x = x, list(truncation = truncation), list(x = x), log = log))
  } else if (length(dist) > 1) {

    # Composite distribution
    d <- dcomposite(x = x, dist = dist, coeff = coeff, log = log)
  } else if (length(prior) == 1) {

    # Single distribution
    d <- do.call(paste0("d", dist), c(list(x = x), coeff, log = log))
  } else {


    # Finite Mixture

    if (log) {
      coeff <- as.matrix(coeff)
      d <- sapply(1:length(prior), function(i, coeff) {
        coeff.int <- coeff[, i]
        names(coeff.int) <- rownames(coeff)
        return((log(prior[i]) + do.call(paste0("d", dist), c(list(x = x), coeff.int, log = TRUE))))
      }, coeff = coeff)

      colnames(d) <- names(prior)

      if (compress) {
        d <- log(rowSums(exp(d)))
      }
    } else {
      coeff <- as.matrix(coeff)
      d <- sapply(1:length(prior), function(i, coeff, log) {
        coeff.int <- coeff[, i]
        names(coeff.int) <- rownames(coeff)
        return(prior[i] * do.call(paste0("d", dist), c(list(x = x), coeff.int, log = FALSE)))
      }, coeff = coeff, log = log)

      colnames(d) <- names(prior)

      if (compress) {
        d <- rowSums(d)
      }
    }
  }

  return(d)
}

#' @rdname combdist
#' @export

pcombdist <- function(q, dist, prior = 1, coeff, log.p = FALSE, lower.tail = TRUE, compress = TRUE, lowertrunc = NULL, uppertrunc = NULL) {
  dist <- unlist(strsplit(dist, split = "_"))

  if (length(grep("trunc", dist)) > 0) {

    # Single Truncated distribution
    if (length(prior) == 1) {

      # Single distribution
      p <- ptruncdist(q = q, dist = dist, coeff = coeff, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    } else {
      if (log.p) {
        # Finite mixture
        coeff <- as.matrix(coeff)
        p <- sapply(1:length(prior), function(i, coeff, log.p, lower.tail, lowertrunc, uppertrunc) {
          coeff.int <- coeff[, i]
          names(coeff.int) <- rownames(coeff)
          return(log(prior[i]) + ptruncdist(q = q, coeff = coeff.int, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc))
        }, coeff = coeff, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
        p <- matrix(p, nrow = length(q), ncol = length(prior))
        colnames(p) <- names(prior)
        if (compress & length(prior) > 1) {
          p <- log(rowSums(exp(p)))
        }
      } else {

        # Finite mixture
        coeff <- as.matrix(coeff)
        p <- sapply(1:length(prior), function(i, coeff, log.p, lower.tail, lowertrunc, uppertrunc) {
          coeff.int <- coeff[, i]
          names(coeff.int) <- rownames(coeff)
          return(prior[i] * ptruncdist(q = q, coeff = coeff.int, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc))
        }, coeff = coeff, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
        p <- matrix(p, nrow = length(q), ncol = length(prior))
        colnames(p) <- names(prior)
        if (compress & length(prior) > 1) {
          p <- rowSums(p)
        }
      }
    }
  } else if (length(dist) > 1) {

    # Composite distribution
    p <- pcomposite(q = q, dist = dist, coeff = coeff, log.p = log.p, lower.tail = lower.tail)
  } else if (length(prior) == 1) {

    # Single distribution
    p <- do.call(paste0("p", dist), c(list(q = q), coeff, list(log.p = log.p), list(lower.tail = lower.tail)))
  } else {
    if (log.p) {
      # Finite mixture
      coeff <- as.matrix(coeff)
      p <- sapply(1:length(prior), function(i, coeff, log.p, lower.tail) {
        coeff.int <- coeff[, i]
        names(coeff.int) <- rownames(coeff)
        return(log(prior[i]) + do.call(paste0("p", dist), c(list(q = q), coeff.int, list(log.p = log.p), list(lower.tail = lower.tail))))
      }, coeff = coeff, log.p = log.p, lower.tail = lower.tail)
      p <- matrix(p, nrow = length(q), ncol = length(prior))
      colnames(p) <- names(prior)
      if (compress & length(prior) > 1) {
        p <- log(rowSums(exp(p)))
      }
    } else {
      # Finite mixture
      coeff <- as.matrix(coeff)
      p <- sapply(1:length(prior), function(i, coeff, log.p, lower.tail) {
        coeff.int <- coeff[, i]
        names(coeff.int) <- rownames(coeff)
        return(prior[i] * do.call(paste0("p", dist), c(list(q = q), coeff.int, list(log.p = log.p), list(lower.tail = lower.tail))))
      }, coeff = coeff, log.p = log.p, lower.tail = lower.tail)
      p <- matrix(p, nrow = length(q), ncol = length(prior))
      colnames(p) <- names(prior)
      if (compress & length(prior) > 1) {
        p <- rowSums(p)
      }
    }
  }

  return(p)
}

#' @rdname combdist
#' @export

qcombdist <- function(p, dist, prior, coeff, log.p = FALSE, lower.tail = TRUE) {
  dist <- unlist(strsplit(dist, split = "_"))

  if (length(grep("trunc", dist)) > 0) {
    if (length(prior) == 1 & !is.matrix(coeff)) {

      # Single truncated distribution
      q <- qtruncdist(p = p, dist = dist, coeff = coeff, log.p = log.p, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    } else {
      p <- if (log.p) {
        exp(p)
      } else {
        p
      }
      p <- if (lower.tail) {
        p
      } else {
        1 - p
      }

      # Finite Mixture of truncate distributions
      qf <- function(p) {
        uniroot(function(q, dist, prior, coeff, p, lowertrunc, uppertrunc) {
          pcombdist(q = q, dist = dist, prior = prior, coeff = coeff, lowertrunc = lowertrunc, uppertrunc = uppertrunc) - p
        }, interval = c(0, 100), extendInt = "yes", tol = 1e-10, dist = dist, prior = prior, coeff = coeff, p = p, lowertrunc = lowertrunc, uppertrunc = uppertrunc)$root
      }
      qf <- Vectorize(qf)
      q <- qf(p)
    }
  } else if (is.null(coeff)) {

    # Empirical distribution
    q <- do.call(paste0("q", dist), c(list(p = p), list(coeff), lower.tail = lower.tail, log.p = log.p))
  } else if (length(dist) > 1) {

    # Composite distribution
    q <- qcomposite(p = p, dist = dist, coeff = coeff, log.p = log.p, lower.tail = lower.tail)
  } else if (length(prior) == 1) {

    # Single distribution
    q <- do.call(paste0("q", dist), c(list(p = p), coeff, log.p = log.p, lower.tail = lower.tail))
  } else {
    p <- if (log.p) {
      exp(p)
    } else {
      p
    }
    p <- if (lower.tail) {
      p
    } else {
      1 - p
    }

    qf <- function(p) {
      uniroot(function(q, dist, prior, coeff, p) {
        pcombdist(q = q, dist = dist, prior = prior, coeff = coeff) - p
      }, interval = c(0, 100), extendInt = "yes", tol = 1e-10, dist = dist, prior = prior, coeff = coeff, p = p)$root
    }
    qf <- Vectorize(qf)
    q <- qf(p)
  }

  return(q)
}

#' @rdname combdist
#' @export

mcombdist <- function(r, truncation = NULL, dist, prior = 1, coeff, lower.tail = TRUE, compress = TRUE, uppertrunc = 0, lowertrunc = Inf) {
  dist <- unlist(strsplit(dist, split = "_"))

  if (length(grep("trunc", dist)) > 0) {
    if (length(prior) == 1 & !is.matrix(coeff)) {

      # Single truncated distribution
      m <- mtruncdist(r = r, truncation = truncation, coeff = coeff, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    } else {

      # Finite Mixture truncated
      coeff <- as.matrix(coeff)
      m <- sapply(1:length(prior), function(i, coeff, truncation, r, lower.tail, lowertrunc, uppertrunc) {
        coeff.int <- coeff[, i]
        names(coeff.int) <- rownames(coeff)
        return(prior[i] * mtruncdist(r = r, truncation = truncation, coeff = coeff.int, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc))
      }, coeff = coeff, truncation = truncation, r = r, lower.tail = lower.tail, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
      m <- matrix(m, nrow = length(truncation), ncol = length(prior))
      colnames(m) <- names(prior)
      if (compress & length(prior) > 1) {
        m <- rowSums(m)
      }
    }
  } else if (length(dist) > 1) {

    # Composite distribution
    m <- mcomposite(r = r, truncation = truncation, dist = dist, coeff = coeff, lower.tail = lower.tail)
  } else if (length(prior) == 1 & !is.matrix(coeff)) {

    # Single distribution
    m <- do.call(paste0("m", dist), c(list(r = r), list(truncation = truncation), coeff, list(lower.tail = lower.tail)))
  } else {


    # Finite Mixture
    coeff <- as.matrix(coeff)
    m <- sapply(1:length(prior), function(i, coeff, truncation, r, lower.tail) {
      coeff.int <- coeff[, i]
      names(coeff.int) <- rownames(coeff)
      return(prior[i] * do.call(paste0("m", dist), c(list(r = r), list(truncation = truncation), coeff.int, list(lower.tail = lower.tail))))
    }, coeff = coeff, truncation = truncation, r = r, lower.tail = lower.tail)
    m <- matrix(m, nrow = length(truncation), ncol = length(prior))
    colnames(m) <- names(prior)
    if (compress & length(prior) > 1) {
      m <- rowSums(m)
    }
  }

  return(m)
}


#' @rdname combdist
#' @export

rcombdist <- function(n, dist, prior, coeff, uppertrunc = NULL, lowertrunc = NULL) {
  dist <- unlist(strsplit(dist, split = "_"))

  if (length(grep("trunc", dist)) > 0) {
    if (!is.matrix(coeff)) {
      x <- rtruncdist(n = n, dist = dist, coeff = coeff, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    } else {
      components <- sample(1:length(prior), prob = prior, size = n, replace = TRUE)

      coeff.list <- lapply(1:nrow(coeff), function(i) {
        coeff[i, components]
      })
      names(coeff.list) <- rownames(coeff)

      x <- rtruncdist(n = n, dist = dist, coeff = coeff.list, lowertrunc = lowertrunc, uppertrunc = uppertrunc)
    }
  } else if (length(dist) > 1) {

    # Composite distribution
    x <- rcomposite(n = n, dist = dist, coeff = coeff)
  } else if (!is.matrix(coeff)) {
    x <- do.call(paste0("r", dist), c(list(n = n), coeff))
  } else {
    components <- sample(1:length(prior), prob = prior, size = n, replace = TRUE)

    coeff.list <- lapply(1:nrow(coeff), function(i) {
      coeff[i, components]
    })
    names(coeff.list) <- rownames(coeff)

    x <- do.call(paste0("r", dist), c(n = n, coeff.list))
  }

  return(x)
}

#' Combined coefficients of power-law transformed combined distribution
#'
#' Coefficients of a power-law transformed combined distribution
#' @param dist character vector denoting the distribution(s).
#' @param prior Numeric vector of prior coefficients, defaults to single vector with value one.
#' @param coeff list of parameters for the distribution(s).
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#' @param nested logical indicating whether results should be returned in a nested list or flat list, defaults to FALSE.
#'
#'
#' @return Returns a nested or flat list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#' @examples
#'
#' \donttest{
#'
#' # Load necessary tools
#' data("fit_US_cities")
#' library(tidyverse)
#'
#'
#' ## Comparing probabilites of power-law transformed transformed variables
#' prob <- fit_US_cities %>%
#'   filter(!(dist %in% c(
#'     "exp", "invpareto_exp_pareto", "exp_pareto", "invpareto_exp",
#'     "gamma", "invpareto_gamma_pareto", "gamma_pareto", "invpareto_gamma"
#'   ))) %>%
#'   group_by(dist, components, np, n) %>%
#'   do(prob = pcombdist(q = 1.1, dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]])) %>%
#'   unnest(cols = c(prob))
#' fit_US_cities_plt <- fit_US_cities %>%
#'   filter(!(dist %in% c(
#'     "exp", "invpareto_exp_pareto", "exp_pareto", "invpareto_exp",
#'     "gamma", "invpareto_gamma_pareto", "gamma_pareto", "invpareto_gamma"
#'   ))) %>%
#'   group_by(dist, components, np, n, convergence) %>%
#'   do(results = as_tibble(combdist_plt(dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]], a = 2, b = 0.5, nested = TRUE))) %>%
#'   unnest(cols = c(results))
#' prob$prob_plt <- fit_US_cities_plt %>%
#'   group_by(dist, components, np, n) %>%
#'   do(prob_plt = pcombdist(q = 2 * 1.1^0.5, dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]])) %>%
#'   unnest(cols = c(prob_plt)) %>%
#'   .$prob_plt
#' prob <- prob %>%
#'   mutate(check = abs(prob - prob_plt))
#'
#' prob <- fit_US_cities %>%
#'   filter(!(dist %in% c(
#'     "exp", "invpareto_exp_pareto", "exp_pareto", "invpareto_exp",
#'     "gamma", "invpareto_gamma_pareto", "gamma_pareto", "invpareto_gamma"
#'   ))) %>%
#'   group_by(dist, components, np, n) %>%
#'   do(prob = pcombdist(q = 2 * 1.1^0.5, dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]])) %>%
#'   unnest(cols = c(prob))
#' fit_US_cities_plt <- fit_US_cities %>%
#'   filter(!(dist %in% c(
#'     "exp", "invpareto_exp_pareto", "exp_pareto", "invpareto_exp",
#'     "gamma", "invpareto_gamma_pareto", "gamma_pareto", "invpareto_gamma"
#'   ))) %>%
#'   group_by(dist, components, np, n, convergence) %>%
#'   do(results = as_tibble(combdist_plt(dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]], a = 2, b = 0.5, nested = TRUE, inv = TRUE))) %>%
#'   unnest(cols = c(results))
#' prob$prob_plt <- fit_US_cities_plt %>%
#'   group_by(dist, components, np, n) %>%
#'   do(prob_plt = pcombdist(q = 1.1, dist = .[["dist"]], prior = .[["prior"]][[1]],
#'   coeff = .[["coefficients"]][[1]])) %>%
#'   unnest(cols = c(prob_plt)) %>%
#'   .$prob_plt
#' prob <- prob %>%
#'   mutate(check = abs(prob - prob_plt))
#' }
#'
#' @export

combdist_plt <- function(dist, prior = NULL, coeff, a = 1, b = 1, inv = FALSE, nested = FALSE) {
  dist <- unlist(strsplit(dist, split = "_"))
  dist <- gsub("trunc", "", dist)

  if (length(dist) > 1) {

    # Composite distribution
    newcoeff <- composite_plt(dist = dist, coeff = coeff, a = a, b = b, inv = inv)$coefficients
  } else if (!is.matrix(coeff)) {

    # Single distribution
    newcoeff <- do.call(paste0(dist, "_plt"), c(coeff, list(a = a), list(b = b), list(inv = inv)))$coefficients
  } else {

    # Finite Mixture
    newcoeff <- sapply(1:ncol(coeff), function(i, coeff, a, b, inv) {
      coeff.int <- coeff[, i]
      names(coeff.int) <- rownames(coeff)
      return(do.call(paste0(dist, "_plt"), c(coeff.int, list(a = a), list(b = b), list(inv = inv)))$coefficients)
    }, coeff = coeff, a = a, b = b, inv = inv)
  }

  if (is.null(prior)) {
    return_list <- list(coefficients = newcoeff)
  } else {
    if (nested) {
      return_list <- list(prior = list(prior), coefficients = list(newcoeff))
    } else {
      return_list <- list(prior = (prior), coefficients = (newcoeff))
    }
  }
  return(return_list)
}

#' Combined distributions MLE
#'
#' Maximum Likelihood estimation for combined ( single, composite and finite mixture) truncated or complete distributions.
#' @param x data vector
#' @param start named numeric vector holding the starting values for the coefficients.
#' @param lower,upper Lower and upper bounds to the estimated coefficients, defaults to -Inf and Inf respectively.
#' @param components number of components for a mixture distribution.
#' @param steps number of steps taken in stepflexmix, defaults to 1.
#' @param dist character vector denoting the distribution(s).
#' @param lowertrunc,uppertrunc lowertrunc- and uppertrunc truncation points, defaults to 0 and Inf respectively
#' @param nested logical indicating whether results should be returned in a nested list or a flat list  form, defaults to FALSE.
#' @param ... Additional arguments.
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{dist}{Character vector denoting the distributions, separated by an underscore}
#' \item{components}{Nr. of combined distributions}
#' \item{prior}{Weights assigned to the respective component distributions}
#' \item{coefficients}{Named vector of coefficients}
#' \item{convergence}{logical indicator of convergence}
#' \item{n}{Length of the fitted data vector}
#' \item{np}{Nr. of coefficients}
#' }
#'
#' @examples
#'
#' \donttest{
#' x <- rdoubleparetolognormal(1e3)
#' combdist.mle(x = x, dist = "doubleparetolognormal") # Double-Pareto Lognormal
#' combdist.mle(x = x, components = 2, dist = "lnorm", steps = 20) # FMM with 2 components
#' combdist.mle( x = x, dist = c("invpareto", "lnorm", "pareto"),
#' start = c(coeff1.k = 1, coeff2.meanlog = mean(log(x)), coeff2.sdlog = sd(log(x)), coeff3.k = 1),
#' lower = c(1e-10, -Inf, 1e-10, 1e-10), upper = c(Inf, Inf, Inf, Inf), nested = TRUE)
#' # composite distribution
#' }
#'
#' @export

combdist.mle <- function(x, dist, start = NULL, lower = NULL, upper = NULL, components = 1, nested = FALSE, steps = 1, lowertrunc = 0, uppertrunc = Inf,  ...) {
  dist <- unlist(strsplit(dist, split = "_"))

  # estimation
  if (length(dist) > 1) {

    # Composite distrbution
    output <- composite.mle(x = x, dist = dist, start = start, ...)

    dist <- paste0(dist, collapse = "_")
    components <- length(dist)
    coeff <- output$coefficients
    convergence <- output$convergence
    np <- length(output$coefficients)
    prior <- 1
  } else if (dist %in% c("pareto", "invpareto", "leftparetolognormal", "doubleparetolognormal", "rightparetolognormal")) {
    output <- do.call(paste0(dist, ".mle"), c(list(x = x), lower = lower, upper = upper, ...))

    if (dist %in% c("leftparetolognormal", "rightparetolognormal")) {
      components <- 2
    } else if (dist %in% c("pareto", "invpareto")) {
      components <- 1
    } else {
      components <- 3
    }

    coeff <- output$coefficients
    convergence <- output$convergence
    np <- length(output$coefficients)
    prior <- 1
  } else {

    # Finite Mixtures
    if (dist %in% c("lnorm", "exp", "gamma", "weibull", "burr", "frechet", "lnormtrunc")) {

      if (!is.null(lowertrunc)) {
        x <- x[(x >= lowertrunc & x <= uppertrunc)]
        init.cluster <- kmeans(x, components)$cluster
        output <- try(flexmix::stepFlexmix(x ~ 1, k = components, cluster = init.cluster, model = FLXMCdist1(dist = dist, lowertrunc = lowertrunc, uppertrunc = uppertrunc), nrep = steps, control = list(minprior = 0, iter.max = 1e4), ...))
      } else {
        output <- try(flexmix::initFlexmix(x ~ 1, k = components, init = list(name = "sem.em"), model = FLXMCdist1(dist = dist), nrep = steps, control = list(minprior = 0, iter.max = 1e4), ...))
      }

      prior <- output@prior
      coeff <- parameters(output)

      if (dist == "exp") {
        coeff <- matrix(coeff, ncol = components)
        rownames(coeff) <- "rate"
      }

      convergence <- output@converged
      np <- output@df
    } else if (dist %in% c("norm")) {

      # Initialization
      if (is.null(start)) {
        start <- kmeans(x, components)$cluster
      }

      output <- try(stepFlexmix(x ~ 1, k = components, cluster = start, model = FLXMRglm(family = "gaussian"), control = list(minprior = 0, iter.max = 1e4), nrep = steps))

      prior <- output@prior
      pars.tab <- parameters(output)
      rownames(pars.tab) <- c("mean", "sd")

      convergence <- output@converged
      np <- output@df
    }
  }

  # Recover results
  if (inherits(output, "try-error")) {
    fit <- list(dist = dist, components = components, pars = NA, np = 0, n = n)
  }

  return_list <- list(dist = dist, components = components, prior = prior, coefficients = coeff, np = np, n = length(x), convergence = convergence)
  if (nested) {
    return_list <- list(dist = dist, components = components, prior = list(prior), coefficients = list(coeff), np = np, n = length(x), convergence = convergence)
  }

  return(return_list)
}
