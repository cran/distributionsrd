#' The two- or three-composite distribution
#'
#' Density, distribution function, quantile function, raw moments and random generation for the two- or three-composite distribution.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param r rth raw moment of the Pareto distribution
#' @param truncation lower truncation parameter
#' @param dist character vector denoting the distribution of the first-, second- (and third) component respectively. If only two components are provided, the distribution reduces to the two-component distribution.
#' @param coeff named numeric vector holding the coefficients of the first-, second- (and third) component, predeced by coeff1., coeff2. (and  coeff3.), respectively. Coefficients for the last component do not have to be provided for the two-component distribution and will be disregarded.
#' @param startc starting values for the lower and upper cutoff, defaults to c(1,1).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities (moments) are \eqn{P[X \leq x]} \eqn{(E[x^r|X \leq y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}
#'
#' @details These derivations are based on the two-composite distribution proposed by \insertCite{bakar2015modelling}{distributionsrd}. Probability Distribution Function:
#'
#'  \deqn{f(x) =  \{ \begin{array}{lrl}
#'  \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \frac{m_1(x)}{M_1(c_1)}  & { if} & 0<x \leq c_1 \\
#'  \frac{1}{1 + \alpha_1 + \alpha_2} \frac{m_2(x)}{M_2(c_2) - M_2(c_1)} & { if} & c_1<x \leq c_2 \\
#'  \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{m_3(x)}{1-M_3(c_2)}  & { if} & c_{2}<x < \infty \\
#'  \end{array} .  }
#'
#'  Cumulative Distribution Function:
#'
#'  \deqn{  \{
#'  \begin{array}{lrl}
#'  \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \frac{M_1(x)}{M_1(c_1)} & { if} & 0<x \leq c_1 \\
#'  \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} + \frac{1}{1 + \alpha_1 + \alpha_2}\frac{M_2(x) - M_2(c_1)}{M_2(c_2) - M_2(c_1)} & { if} & c_1<x \leq c_2 \\
#'  \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{M_3(x) - M_3(c_2)}{1-M_3(c_2)} & { if} & c_{2}<x < \infty \\
#'  \end{array}
#'  .
#'  }
#'
#'  Quantile function
#'
#'  \deqn{
#'  Q(p) =  \{
#'  \begin{array}{lrl}
#'  Q_1( \frac{1 + \alpha_1 + \alpha_2}{\alpha_1} p  M_1(c_1) )  & { if} & 0<x \leq \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \\
#'  Q_2[((p - \frac{\alpha_1}{1 + \alpha_1 + \alpha_2})(1 + \alpha_1 + \alpha_2)(M_2(c_2) - M_2(c_1))) + M_2(c_1)] & { if} & \frac{\alpha_1}{1 + \alpha_1 + \alpha_2}<x \leq \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2} \\
#'  Q_3[((p - \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2})(\frac{1 + \alpha_1 + \alpha_2}{\alpha_2})(1 - M_3(c_2))) + M_3(c_2)] & { if} & \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2}<x < \infty \\
#'  \end{array}
#'  .
#'  }
#'
#'  The lower y-bounded r-th raw moment of the distribution equals
#'
#'  \deqn{
#'  \mu_y^r =
#'  \{
#'  \begin{array}{lrl}
#'  \frac{\alpha_1}{1 + \alpha_1 + \alpha_2}\frac{ (\mu_1)_y^r - (\mu_1)_{c_1}^r}{M_1(c_1)} +  \frac{1}{1 + \alpha_1 + \alpha_2}\frac{ (\mu_2)_{c_1}^r - (\mu_2)_{c_2}^r }{M_2(c_2) - M_2(c_1)} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_y^r}{1-M_3(c_2)} & { if} & 0< y \leq c_2 \\
#'  \frac{1}{1 + \alpha_1 + \alpha_2} \frac{(\mu_2)_y^r - (\mu_2)_{c_2}^r }{M_2(c_2) - M_2(c_1)} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_{c_2}^r}{1-M_3(c_2)} & { if} & c_1< y \leq c_2\\
#'  \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_y^r}{1-M_3(c_2)}  & { if} & c_2< y < \infty \\
#'  \end{array}
#'  .
#'  }
#'
#' @return dcomposite returns the density, pcomposite the distribution function, qcomposite the quantile function, mcomposite the rth moment of the distribution and rcomposite generates random deviates.
#'
#'  The length of the result is determined by n for rcomposite, and is the maximum of the lengths of the numerical arguments for the other functions.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#'
#' #' ## Three-component distribution
#' dist <- c("invpareto", "lnorm", "pareto")
#' coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1.5, coeff1.k = 1.5)
#'
#' # Compare density with the Double-Pareto Lognormal distribution
#' plot(x = seq(0, 5, length.out = 1e3), y = dcomposite(x = seq(0, 5, length.out = 1e3),
#' dist = dist, coeff = coeff))
#' lines(x = seq(0, 5, length.out = 1e3), y = ddoubleparetolognormal(x = seq(0, 5, length.out = 1e3)))
#'
#' # Demonstration of log functionality for probability and quantile function
#' qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist,
#' coeff = coeff, log.p = TRUE)
#'
#' # The zeroth truncated moment is equivalent to the probability function
#' pcomposite(2, dist = dist, coeff = coeff)
#' mcomposite(truncation = 2, dist = dist, coeff = coeff)
#'
#' # The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 3, coeff1.k = 1.5)
#' x <- rcomposite(1e5, dist = dist, coeff = coeff)
#'
#' mean(x)
#' mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff)
#'
#' ## Two-component distribution
#' dist <- c("lnorm", "pareto")
#' coeff <- coeff <- c(coeff2.k = 1.5, coeff1.meanlog = -0.5, coeff1.sdlog = 0.5)
#'
#' # Compare density with the Right-Pareto Lognormal distribution
#' plot(x = seq(0, 5, length.out = 1e3), y = dcomposite(x = seq(0, 5, length.out = 1e3),
#' dist = dist, coeff = coeff))
#' lines(x = seq(0, 5, length.out = 1e3), y = drightparetolognormal(x = seq(0, 5, length.out = 1e3)))
#'
#' # Demonstration of log functionality for probability and quantile function
#' qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist,
#' coeff = coeff, log.p = TRUE)
#'
#' # The zeroth truncated moment is equivalent to the probability function
#' pcomposite(2, dist = dist, coeff = coeff)
#' mcomposite(truncation = 2, dist = dist, coeff = coeff)
#'
#' # The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#' #for large enough samples.
#' coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 3)
#' x <- rcomposite(1e5, dist = dist, coeff = coeff)
#'
#' mean(x)
#' mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff)
#'
#' sum(x[x > quantile(x, 0.1)]) / length(x)
#' mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff)
#' @name composite

NULL

#' @rdname composite
#' @export

dcomposite <- function(x, dist, coeff, startc = c(1, 1), log = FALSE) {
  list2env(coeffcomposite(dist = dist, coeff = coeff, startc), envir = environment())

  # Call functions
  fun1 <- function(d, z, logit = FALSE) {
    do.call(paste0(d, dist1), c(list(z), coeff1, list(log = logit)))
  }
  fun2 <- function(d, z, logit = FALSE) {
    do.call(paste0(d, dist2), c(list(z), coeff2, list(log = logit)))
  }
  fun3 <- function(d, z, logit = FALSE) {
    do.call(paste0(d, dist3), c(list(z), coeff3, list(log = logit)))
  }


  if (log) {
    # Calculate density
    d <- rep(NA, length(x))
    d[(x <= c1)] <- log(alpha1 / (1 + alpha1 + alpha2)) + fun1("d", x[(x <= c1)], logit = TRUE) - log(fun1("p", c1))
    d[(c1 < x & x <= c2)] <- log(1 / (1 + alpha1 + alpha2)) + fun2("d", x[(c1 < x & x <= c2)], logit = TRUE) - log((fun2("p", c2) - fun2("p", c1)))
    d[x > c2] <- log(alpha2 / (1 + alpha1 + alpha2)) + fun3("d", x[(x > c2)], logit = TRUE) - log(ifelse(alpha2 == 0, 1, (1 - fun3("p", c2))))
  } else {
    # Calculate density
    d <- rep(NA, length(x))
    d[(x <= c1)] <- alpha1 / (1 + alpha1 + alpha2) * fun1("d", x[(x <= c1)]) / fun1("p", c1)
    d[(c1 < x & x <= c2)] <- 1 / (1 + alpha1 + alpha2) * fun2("d", x[(c1 < x & x <= c2)]) / (fun2("p", c2) - fun2("p", c1))
    d[x > c2] <- alpha2 / (1 + alpha1 + alpha2) * fun3("d", x[(x > c2)]) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))
  }

  return(d)
}

#' @rdname composite
#' @export

pcomposite <- function(q, dist, coeff, startc = c(1, 1), log.p = FALSE, lower.tail = TRUE) {
  list2env(coeffcomposite(dist = dist, coeff = coeff, startc), envir = environment())

  # Call functions
  fun1 <- function(d, z) {
    do.call(paste0(d, dist1), c(list(z), coeff1))
  }
  fun2 <- function(d, z) {
    do.call(paste0(d, dist2), c(list(z), coeff2))
  }
  fun3 <- function(d, z) {
    do.call(paste0(d, dist3), c(list(z), coeff3))
  }

  # Calculate probability
  p <- rep(NA, length(q))
  p[(q <= c1)] <- alpha1 / (1 + alpha1 + alpha2) * fun1("p", q[(q <= c1)]) / fun1("p", c1)
  p[(c1 < q & q <= c2)] <- alpha1 / (1 + alpha1 + alpha2) + 1 / (1 + alpha1 + alpha2) * (fun2("p", q[(c1 < q & q <= c2)]) - fun2("p", c1)) / (fun2("p", c2) - fun2("p", c1))
  p[q > c2] <- (alpha1 + 1) / (1 + alpha1 + alpha2) + alpha2 / (1 + alpha1 + alpha2) * (fun3("p", q[q > c2]) - fun3("p", c2)) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))

  p <- if (lower.tail) {
    p
  } else {
    1 - p
  }
  p <- if (log.p) {
    log(p)
  } else {
    p
  }
  return(p)
}

#' @rdname composite
#' @export

qcomposite <- function(p, dist, coeff, startc = c(1, 1), log.p = FALSE, lower.tail = TRUE) {
  list2env(coeffcomposite(dist = dist, coeff = coeff, startc), envir = environment())

  # Call functions
  fun1 <- function(d, z) {
    do.call(paste0(d, dist1), c(list(z), coeff1))
  }
  fun2 <- function(d, z) {
    do.call(paste0(d, dist2), c(list(z), coeff2))
  }
  fun3 <- function(d, z) {
    do.call(paste0(d, dist3), c(list(z), coeff3))
  }

  p <- if (log.p) {
    exp(p)
  } else {
    p
  }
  # Calculate quantiles
  p <- if (lower.tail) {
    p
  } else {
    1 - p
  }


  q <- rep(NA, length(q))
  q[(p <= (alpha1 / (1 + alpha1 + alpha2)))] <- fun1("q", ((1 + alpha1 + alpha2) / alpha1 * p[(p <= (alpha1 / (1 + alpha1 + alpha2)))] * fun1("p", c1)))
  q[((alpha1 / (1 + alpha1 + alpha2)) < p & p <= ((1 + alpha1) / (1 + alpha1 + alpha2)))] <- fun2("q", ((p[((alpha1 / (1 + alpha1 + alpha2)) < p & p <= ((1 + alpha1) / (1 + alpha1 + alpha2)))] -
    (alpha1 / (1 + alpha1 + alpha2))) * (1 + alpha1 + alpha2) * (fun2("p", c2) - fun2("p", c1)) +
    fun2("p", c1)))
  q[(p > ((1 + alpha1) / (1 + alpha1 + alpha2)))] <- fun3("q", ((p[(p > ((1 + alpha1) / (1 + alpha1 + alpha2)))] -
    ((1 + alpha1) / (1 + alpha1 + alpha2))) * (1 + alpha1 + alpha2) / alpha2 * (1 - fun3("p", c2)) +
    fun3("p", c2)))
  # For some obscure reason, qpareto does not return the correct value if p=1 for the assigned values
  q[p == 1] <- fun3("q", 1)

  return(q)
}

#' @rdname composite
#' @export

mcomposite <- function(r = 0, truncation = 0, dist, coeff, startc = c(1, 1), lower.tail = TRUE) {
  list2env(coeffcomposite(dist = dist, coeff = coeff, startc), envir = environment())

  # Call functions
  fun1 <- function(d, z) {
    do.call(paste0(d, dist1), c(list(z), coeff1))
  }
  fun2 <- function(d, z) {
    do.call(paste0(d, dist2), c(list(z), coeff2))
  }
  fun3 <- function(d, z) {
    do.call(paste0(d, dist3), c(list(z), coeff3))
  }

  # Call functions
  mfun1 <- function(trunc) {
    do.call(paste0("m", dist1), c(list(r = r), coeff1, list(truncation = trunc, lower.tail = FALSE)))
  }
  mfun2 <- function(trunc) {
    do.call(paste0("m", dist2), c(list(r = r), coeff2, list(truncation = trunc, lower.tail = FALSE)))
  }
  mfun3 <- function(trunc) {
    do.call(paste0("m", dist3), c(list(r = r), coeff3, list(truncation = trunc, lower.tail = FALSE)))
  }

  # Calculate raw moment
  m <- rep(NA, length(truncation))
  m[(truncation <= c1)] <- alpha1 / (1 + alpha1 + alpha2) * (mfun1(truncation[(truncation <= c1)]) - mfun1(c1)) / fun1("p", c1) +
    1 / (1 + alpha1 + alpha2) * (mfun2(c1) - mfun2(c2)) / (fun2("p", c2) - fun2("p", c1)) +
    alpha2 / (1 + alpha1 + alpha2) * mfun3(c2) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))
  m[(c1 < truncation & truncation <= c2)] <- 1 / (1 + alpha1 + alpha2) * (mfun2(truncation[(c1 < truncation & truncation <= c2)]) - mfun2(c2)) / (fun2("p", c2) - fun2("p", c1)) +
    alpha2 / (1 + alpha1 + alpha2) * mfun3(c2) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))
  m[truncation > c2] <- alpha2 / (1 + alpha1 + alpha2) * mfun3(truncation[(truncation > c2)]) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))

  notrunc <- alpha1 / (1 + alpha1 + alpha2) * (mfun1(0) - mfun1(c1)) / fun1("p", c1) +
    1 / (1 + alpha1 + alpha2) * (mfun2(c1) - mfun2(c2)) / (fun2("p", c2) - fun2("p", c1)) +
    alpha2 / (1 + alpha1 + alpha2) * mfun3(c2) / ifelse(alpha2 == 0, 1, (1 - fun3("p", c2)))

  if (lower.tail) {
    m <- notrunc - m
  }

  return(m)
}

#' @rdname composite
#' @export

rcomposite <- function(n, dist, coeff, startc = c(1, 1)) {
  r <- qcomposite(runif(n), dist = dist, coeff = coeff, startc = startc)
  return(r)
}

#' Composite coefficients after power-law transformation
#'
#' Coefficients of a power-law transformed composite distribution
#' @param dist character vector denoting the distribution of the first-, second- (and third) component respectively. If only two components are provided, the distribution reduces to the two-component distribution.
#' @param coeff named numeric vector holding the coefficients of the first-, second- (and third) component, predeced by coeff1., coeff2. (and  coeff3.), respectively. Coefficients for the last component do not have to be provided for the two-component distribution and will be disregarded.
#' @param a,b constant and power of power-law transformation, defaults to 1 and 1 respectively.
#' @param inv logical indicating whether coefficients of the outcome variable of the power-law transformation should be returned (FALSE) or whether coefficients of the input variable being power-law transformed should be returned (TRUE). Defaults to FALSE.
#'
#' @return Returns a named list containing
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' }
#'
#'  ## Comparing probabilites of power-law transformed transformed variables
#'  dist <- c("invpareto", "lnorm", "pareto")
#'  coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1.5, coeff1.k = 1.5)
#'
#'  pcomposite(3,dist=dist,coeff=coeff)
#'  newcoeff = composite_plt(dist=dist,coeff=coeff,a=5,b=7)$coefficients
#'  pcomposite(5*3^7,dist=dist,coeff=newcoeff)
#'
#'  pcomposite(5*0.9^3,dist=dist,coeff=coeff)
#'  newcoeff = composite_plt(dist=dist,coeff=coeff,a=5,b=3,inv=TRUE)$coefficients
#'  pcomposite(0.9,dist=dist,coeff=newcoeff)
#'
#' @export

composite_plt <- function(dist, coeff, a = 1, b = 1, inv = FALSE) {
  list2env(coeffcomposite(dist = dist, coeff = coeff), envir = environment())

  # Call functions
  newcoeff1 <- do.call(paste0(dist1, "_plt"), c(coeff1, list(a = a), list(b = b), list(inv = inv)))$coefficients
  newcoeff2 <- do.call(paste0(dist2, "_plt"), c(coeff2, list(a = a), list(b = b), list(inv = inv)))$coefficients
  newcoeff3 <- do.call(paste0(dist3, "_plt"), c(coeff3, list(a = a), list(b = b), list(inv = inv)))$coefficients

  if (dist1 == "invpareto") {
    newcoeff1 <- newcoeff1["k"]
  }
  if (dist2 == "pareto") {
    newcoeff2 <- newcoeff2["k"]
  }
  if (dist3 == "pareto") {
    newcoeff3 <- newcoeff3["k"]
  }

  if (length(dist) == 2) {
    newcoeff <- c(newcoeff1, newcoeff2)
    names(newcoeff) <- c(paste0("coeff1.", names(newcoeff1)), paste0("coeff2.", names(newcoeff2)))
  } else {
    newcoeff <- c(newcoeff1, newcoeff2, newcoeff3)
    names(newcoeff) <- c(paste0("coeff1.", names(newcoeff1)), paste0("coeff2.", names(newcoeff2)), paste0("coeff3.", names(newcoeff3)))
  }

  return_list <- list(coefficients = newcoeff)
  return(return_list)
}


#' Parametrise two-/three- composite distribution
#'
#' Determines the weights and cutoffs of the three-composite distribution numerically applying te continuity- and differentiability condition.
#' @param dist character vector denoting the distribution of the first-, second- (and third) component respectively. If only two components are provided, the distribution reduces to the two-component distribution.
#' @param coeff named numeric vector holding the coefficients of the first-, second- (and third) component, predeced by coeff1., coeff2. (and  coeff3.), respectively. Coefficients for the last component do not have to be provided for the two-component distribution and will be disregarded.
#' @param startc starting values for the lower and upper cutoff, defaults to c(1,1).
#'
#' @details The continuity condition implies
#'
#'  \deqn{ \alpha_1 = \frac{m_2(c_1) M_1(c_1)}{m_1(c_1)[M_2(c_2) - M_2(c_1)]}, \qquad \alpha_2 = \frac{m_2(c_2) [1 - M_3(c_2)]}{m_3(c_2) [M_2(c_2) - M_2(c_1)]} }
#'
#'  The differentiability condition implies
#'
#'  \deqn{ \frac{d}{dc_1} ln[\frac{m_1(c_1)}{m_2(c_1)}] = 0, \qquad \frac{d}{dc_2} ln[\frac{m_2(c_2)}{m_3(c_2)}] = 0 }
#'
#' @return Returns a named list containing a the separate distributions and their respective coefficients, as well as the cutoffs and weights of the composite distribution
#'
#'
#' @examples
#'
#' # Three-composite distribution
#' dist <- c("invpareto", "lnorm", "pareto")
#' coeff <- c(coeff1.k = 1, coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1)
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))
#'
#' # Two-composite distribution
#' dist <- c("lnorm", "pareto")
#' coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 1.5)
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))
#'
#' dist <- c("invpareto", "lnorm")
#' coeff <- c(coeff1.k = 1.5, coeff2.meanlog = 2, coeff2.sdlog = 0.5)
#' coeffcomposite(dist = dist, coeff = coeff, startc = c(1, 1))
#' #'
#' @export

coeffcomposite <- function(dist, coeff, startc = c(1, 1)) {

  # Call functions
  fun1 <- function(d, z) {
    do.call(paste0(d, dist1), c(list(z), coeff1))
  }
  fun2 <- function(d, z) {
    do.call(paste0(d, dist2), c(list(z), coeff2))
  }
  fun3 <- function(d, z) {
    do.call(paste0(d, dist3), c(list(z), coeff3))
  }

  for (i in 1:length(dist)) {
    # Distribution
    assign(paste0("dist", i), dist[i])

    # Coefficients
    assign("int", coeff[grep(paste0("coeff", i, "."), names(coeff))])
    names(int) <- gsub(paste0("coeff", i, "."), "", names(int))
    assign(paste0("coeff", i), int)
  }

  # lapply(1:length(dist), function(i) {
  #   # Distribution
  #   assign(paste0("dist", i), dist[i], envir = parent.frame(n = 2))
  #
  #   # Coefficients
  #   assign("int", coeff[grep(paste0("coeff", i, "."), names(coeff))])
  #   names(int) <- gsub(paste0("coeff", i, "."), "", names(int))
  #   assign(paste0("coeff", i), int, envir = parent.frame(n = 2))
  # })

  if (length(dist) == 3) {

    # Determine cuttofs from differentiability condition
    # Check whether cutoffs are part of the parameter space
    if (dist1 == "invpareto") {
      coeff1 <- c(coeff1, xmax = 1e50)
    }
    if (dist3 == "pareto") {
      coeff3 <- c(coeff3, xmin = 1e-50)
    }

    start <- start
    ratio1 <- function(x) {
      -log(fun2("d", x) / fun1("d", x))
    }
    ratio2 <- function(x) {
      -log(fun2("d", x) / fun3("d", x))
    }

    c1 <- suppressWarnings(nlminb(startc[1], ratio1, lower = 1e-10, upper = Inf, control = list(rel.tol = 1e-15))$par)
    c2 <- suppressWarnings(nlminb(startc[2], ratio2, lower = 1e-10, upper = Inf, control = list(rel.tol = 1e-15))$par)

    # Update cutoffs that are part of the parameter space
    if (dist1 == "invpareto") {
      coeff1["xmax"] <- c1
    }
    if (dist3 == "pareto") {
      coeff3["xmin"] <- c2
    }

    # Determine mixing weights from continuity condition
    alpha1 <- (fun2("d", c1) * (fun1("p", c1))) / (fun1("d", c1) * (fun2("p", c2) - fun2("p", c1)))
    alpha2 <- (fun2("d", c2) * (1 - fun3("p", c2))) / (fun3("d", c2) * (fun2("p", c2) - fun2("p", c1)))

    # Warning message
    if (c2 == Inf) warning("Upper cutoff (c2) has reached maxmium value of Inf. Try fitting distribution without right-side composite.", immediate. = TRUE, call. = FALSE)
  } else {

    # Replacement distribution, has no influence on results
    dist3 <- "lnorm"
    coeff3 <- c(meanlog = 0, sdlog = 0)
    alpha2 <- 0
    c2 <- Inf

    # Determine cuttofs from differentiability condition
    # Check whether cutoffs are part of the parameter space
    if (dist1 == "invpareto") {
      coeff1 <- c(coeff1, xmax = 1e50)
    }
    if (dist2 == "pareto") {
      coeff2 <- c(coeff2, xmin = 1e-50)
    }

    start <- start

    # Make sure we are looking for a minimum
    # Can use abs here because ratio should never go below one, so ratio can never switch sign (either its negative or positive)
    ratio1 <- function(x) {
      -abs(log(fun2("d", x) / fun1("d", x)))
    }

    c1 <- suppressWarnings(nlminb(startc[1], ratio1, lower = 1e-10, upper = Inf, control = list(rel.tol = 1e-15))$par)

    # Update cutoffs that are part of the parameter space
    if (dist1 == "invpareto") {
      coeff1["xmax"] <- c1
    }
    if (dist2 == "pareto") {
      coeff2["xmin"] <- c1
    }

    # Determine mixing weights from continuity condition
    alpha1 <- (fun2("d", c1) * (fun1("p", c1))) / (fun1("d", c1) * (fun2("p", c2) - fun2("p", c1)))
  }

  # Warning message
  if (c1 == 1e-10) warning("Lower cutoff (c1) has reached minimum value of 1e-10. Try fitting distribution without left-side composite.", immediate. = TRUE, call. = FALSE)

  return_list <- list(dist1 = dist1, dist2 = dist2, dist3 = dist3, coeff1 = coeff1, coeff2 = coeff2, coeff3 = coeff3, c1 = c1, c2 = c2, alpha1 = alpha1, alpha2 = alpha2)
  return(return_list)
}


#' Composite MLE
#'
#' Maximum likelihood estimation of the parameters of the two-/three- composite distribution
#' @param x data vector
#' @param dist character vector denoting the distribution of the first-, second- (and third) component respectively. If only two components are provided, the distribution reduces to the two-component distribution.
#' @param start named numeric vector holding the coefficients of the first-, second- (and third) component, predeced by coeff1., coeff2. (and  coeff3.), respectively. Coefficients for the last component do not have to be provided for the two-component distribution and will be disregarded.
#' @param lower,upper Lower and upper bounds to the estimated coefficients, defaults to -Inf and Inf respectively.
#'
#' @return Returns a named list containing a
#' \describe{
#' \item{coefficients}{Named vector of coefficients}
#' \item{convergence}{logical indicator of convergence}
#' \item{cutoffs}{Cutoffs of the composite distribution}
#' \item{n}{Length of the fitted data vector}
#' \item{np}{Nr. of coefficients}
#' \item{components}{Nr. of components}
#' }
#'
#' @examples
#'
#' dist <- c("invpareto", "lnorm", "pareto")
#' coeff <- c(
#'   coeff1.k = 1.5, coeff2.meanlog = -0.5,
#'   coeff2.sdlog = 0.5, coeff3.k = 1.5
#' )
#' lower <- c(1e-10, -Inf, 1e-10, 1e-10)
#' upper <- c(Inf, Inf, Inf, Inf)
#' x <- rcomposite(1e3, dist = dist, coeff = coeff)
#' \donttest{
#' composite.mle(x = x, dist = dist, start = coeff + 0.2, lower = lower, upper = upper)
#' }
#' #'
#' @export

composite.mle <- function(x, dist, start, lower = NULL, upper = NULL) {
  if (is.null(lower)) {
    lower <- rep(-Inf, length(start))
  }
  if (is.null(upper)) {
    upper <- rep(Inf, length(start))
  }

  mle <- function(coeff, dist, x, names) {
    names(coeff) <- names

    d <- dcomposite(x = x, dist = dist, coeff = coeff, log = TRUE)

    # Update cutoff starting values
    nlogl <- -sum(d)
    if (is.nan(nlogl)) {
      nlogl <- 1e20
    }
    if (is.infinite(nlogl)) {
      nlogl <- 1e20
    }

    return(nlogl)
  }

  optim.out <- suppressWarnings(nlminb(start = as.numeric(start), objective = mle, lower = lower, upper = upper, dist = dist, x = x, names = names(start)))
  coeff <- optim.out$par
  names(coeff) <- names(start)

  cutoffs <- coeffcomposite(dist = dist, coeff = coeff)[c("c1", "c2")]
  mode(cutoffs) <- "numeric"

  convergence <- ifelse(optim.out$convergence == 0, TRUE, FALSE)

  return_list <- list(coefficients = coeff, convergence = convergence, cutoffs = cutoffs, np = length(coeff), n = length(x), components = length(dist))
  return(return_list)
}
