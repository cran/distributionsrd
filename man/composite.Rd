% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composite.R
\name{composite}
\alias{composite}
\alias{dcomposite}
\alias{pcomposite}
\alias{qcomposite}
\alias{mcomposite}
\alias{rcomposite}
\title{The two- or three-composite distribution}
\usage{
dcomposite(x, dist, coeff, startc = c(1, 1), log = FALSE)

pcomposite(q, dist, coeff, startc = c(1, 1), log.p = FALSE, lower.tail = TRUE)

qcomposite(p, dist, coeff, startc = c(1, 1), log.p = FALSE, lower.tail = TRUE)

mcomposite(
  r = 0,
  truncation = 0,
  dist,
  coeff,
  startc = c(1, 1),
  lower.tail = TRUE
)

rcomposite(n, dist, coeff, startc = c(1, 1))
}
\arguments{
\item{x, q}{vector of quantiles}

\item{dist}{character vector denoting the distribution of the first-, second- (and third) component respectively. If only two components are provided, the distribution reduces to the two-component distribution.}

\item{coeff}{named numeric vector holding the coefficients of the first-, second- (and third) component, predeced by coeff1., coeff2. (and  coeff3.), respectively. Coefficients for the last component do not have to be provided for the two-component distribution and will be disregarded.}

\item{startc}{starting values for the lower and upper cutoff, defaults to c(1,1).}

\item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; if TRUE (default), probabilities (moments) are \eqn{P[X \leq x]} \eqn{(E[x^r|X \leq y])}, otherwise, \eqn{P[X > x]} \eqn{(E[x^r|X > y])}}

\item{p}{vector of probabilities}

\item{r}{rth raw moment of the Pareto distribution}

\item{truncation}{lower truncation parameter}

\item{n}{number of observations}
}
\value{
dcomposite returns the density, pcomposite the distribution function, qcomposite the quantile function, mcomposite the rth moment of the distribution and rcomposite generates random deviates.

 The length of the result is determined by n for rcomposite, and is the maximum of the lengths of the numerical arguments for the other functions.
}
\description{
Density, distribution function, quantile function, raw moments and random generation for the two- or three-composite distribution.
}
\details{
These derivations are based on the two-composite distribution proposed by \insertCite{bakar2015modelling}{distributionsrd}. Probability Distribution Function:

 \deqn{f(x) =  \{ \begin{array}{lrl}
 \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \frac{m_1(x)}{M_1(c_1)}  & { if} & 0<x \leq c_1 \\
 \frac{1}{1 + \alpha_1 + \alpha_2} \frac{m_2(x)}{M_2(c_2) - M_2(c_1)} & { if} & c_1<x \leq c_2 \\
 \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{m_3(x)}{1-M_3(c_2)}  & { if} & c_{2}<x < \infty \\
 \end{array} .  }

 Cumulative Distribution Function:

 \deqn{  \{
 \begin{array}{lrl}
 \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \frac{M_1(x)}{M_1(c_1)} & { if} & 0<x \leq c_1 \\
 \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} + \frac{1}{1 + \alpha_1 + \alpha_2}\frac{M_2(x) - M_2(c_1)}{M_2(c_2) - M_2(c_1)} & { if} & c_1<x \leq c_2 \\
 \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{M_3(x) - M_3(c_2)}{1-M_3(c_2)} & { if} & c_{2}<x < \infty \\
 \end{array}
 .
 }

 Quantile function

 \deqn{
 Q(p) =  \{
 \begin{array}{lrl}
 Q_1( \frac{1 + \alpha_1 + \alpha_2}{\alpha_1} p  M_1(c_1) )  & { if} & 0<x \leq \frac{\alpha_1}{1 + \alpha_1 + \alpha_2} \\
 Q_2[((p - \frac{\alpha_1}{1 + \alpha_1 + \alpha_2})(1 + \alpha_1 + \alpha_2)(M_2(c_2) - M_2(c_1))) + M_2(c_1)] & { if} & \frac{\alpha_1}{1 + \alpha_1 + \alpha_2}<x \leq \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2} \\
 Q_3[((p - \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2})(\frac{1 + \alpha_1 + \alpha_2}{\alpha_2})(1 - M_3(c_2))) + M_3(c_2)] & { if} & \frac{1+\alpha_1}{1 + \alpha_1 + \alpha_2}<x < \infty \\
 \end{array}
 .
 }

 The lower y-bounded r-th raw moment of the distribution equals

 \deqn{
 \mu_y^r =
 \{
 \begin{array}{lrl}
 \frac{\alpha_1}{1 + \alpha_1 + \alpha_2}\frac{ (\mu_1)_y^r - (\mu_1)_{c_1}^r}{M_1(c_1)} +  \frac{1}{1 + \alpha_1 + \alpha_2}\frac{ (\mu_2)_{c_1}^r - (\mu_2)_{c_2}^r }{M_2(c_2) - M_2(c_1)} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_y^r}{1-M_3(c_2)} & { if} & 0< y \leq c_2 \\
 \frac{1}{1 + \alpha_1 + \alpha_2} \frac{(\mu_2)_y^r - (\mu_2)_{c_2}^r }{M_2(c_2) - M_2(c_1)} + \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_{c_2}^r}{1-M_3(c_2)} & { if} & c_1< y \leq c_2\\
 \frac{\alpha_2}{1 + \alpha_1 + \alpha_2} \frac{(\mu_3)_y^r}{1-M_3(c_2)}  & { if} & c_2< y < \infty \\
 \end{array}
 .
 }
}
\examples{

#' ## Three-component distribution
dist <- c("invpareto", "lnorm", "pareto")
coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 1.5, coeff1.k = 1.5)

# Compare density with the Double-Pareto Lognormal distribution
plot(x = seq(0, 5, length.out = 1e3), y = dcomposite(x = seq(0, 5, length.out = 1e3),
dist = dist, coeff = coeff))
lines(x = seq(0, 5, length.out = 1e3), y = ddoubleparetolognormal(x = seq(0, 5, length.out = 1e3)))

# Demonstration of log functionality for probability and quantile function
qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist,
coeff = coeff, log.p = TRUE)

# The zeroth truncated moment is equivalent to the probability function
pcomposite(2, dist = dist, coeff = coeff)
mcomposite(truncation = 2, dist = dist, coeff = coeff)

# The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#for large enough samples.
coeff <- c(coeff2.meanlog = -0.5, coeff2.sdlog = 0.5, coeff3.k = 3, coeff1.k = 1.5)
x <- rcomposite(1e5, dist = dist, coeff = coeff)

mean(x)
mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff)

sum(x[x > quantile(x, 0.1)]) / length(x)
mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff)

## Two-component distribution
dist <- c("lnorm", "pareto")
coeff <- coeff <- c(coeff2.k = 1.5, coeff1.meanlog = -0.5, coeff1.sdlog = 0.5)

# Compare density with the Right-Pareto Lognormal distribution
plot(x = seq(0, 5, length.out = 1e3), y = dcomposite(x = seq(0, 5, length.out = 1e3),
dist = dist, coeff = coeff))
lines(x = seq(0, 5, length.out = 1e3), y = drightparetolognormal(x = seq(0, 5, length.out = 1e3)))

# Demonstration of log functionality for probability and quantile function
qcomposite(pcomposite(2, dist = dist, coeff = coeff, log.p = TRUE), dist = dist,
coeff = coeff, log.p = TRUE)

# The zeroth truncated moment is equivalent to the probability function
pcomposite(2, dist = dist, coeff = coeff)
mcomposite(truncation = 2, dist = dist, coeff = coeff)

# The (truncated) first moment is equivalent to the mean of a (truncated) random sample,
#for large enough samples.
coeff <- c(coeff1.meanlog = -0.5, coeff1.sdlog = 0.5, coeff2.k = 3)
x <- rcomposite(1e5, dist = dist, coeff = coeff)

mean(x)
mcomposite(r = 1, lower.tail = FALSE, dist = dist, coeff = coeff)

sum(x[x > quantile(x, 0.1)]) / length(x)
mcomposite(r = 1, truncation = quantile(x, 0.1), lower.tail = FALSE, dist = dist, coeff = coeff)
}
\references{
\insertAllCited{}
}
