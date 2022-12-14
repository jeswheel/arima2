% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arima2.R
\name{arima2}
\alias{arima2}
\title{ARIMA Modelling of Time Series}
\usage{
arima2(
  x,
  order = c(0L, 0L, 0L),
  seasonal = list(order = c(0L, 0L, 0L), period = NA),
  xreg = NULL,
  include.mean = TRUE,
  nrestart = 10,
  transform.pars = TRUE,
  fixed = NULL,
  init = NULL,
  method = c("CSS-ML", "ML", "CSS"),
  n.cond,
  SSinit = c("Gardner1980", "Rossignol2011"),
  optim.method = "BFGS",
  optim.control = list(),
  kappa = 1e+06
)
}
\arguments{
\item{x}{a univariate time series}

\item{order}{A specification of the non-seasonal part of the ARIMA
    model: the three integer components \eqn{(p, d, q)} are the AR order, the
    degree of differencing, and the MA order.}

\item{seasonal}{A specification of the seasonal part of the ARIMA
    model, plus the period (which defaults to \code{frequency(x)}).
    This may be a list with components \code{order} and
    \code{period}, or just a numeric vector of length 3 which
    specifies the seasonal \code{order}.  In the latter case the
    default period is used.}

\item{xreg}{Optionally, a vector or matrix of external regressors,
    which must have the same number of rows as \code{x}.}

\item{include.mean}{Should the ARMA model include a mean/intercept term?  The
    default is \code{TRUE} for undifferenced series, and it is ignored
    for ARIMA models with differencing.}

\item{nrestart}{Number of random restarts to use in fitting the model. If
nrestart = 0L, then the function is equivalent to \code{\link[stats:arima]{stats::arima()}}.}

\item{transform.pars}{logical; if true, the AR parameters are
    transformed to ensure that they remain in the region of
    stationarity.  Not used for \code{method = "CSS"}.  For
    \code{method = "ML"}, it has been advantageous to set
    \code{transform.pars = FALSE} in some cases, see also \code{fixed}.}

\item{fixed}{optional numeric vector of the same length as the total
    number of coefficients to be estimated.  It should be of the form
    \deqn{(\phi_1, \ldots, \phi_p, \theta_1, \ldots, \theta_q,
      \Phi_1, \ldots, \Phi_P, \Theta_1, \ldots, \Theta_Q, \mu),
    }{(phi_1, ..., phi_p, theta_1, ..., theta_p,Phi_1,..., Phi_P, Theta_1,...,Theta_Q, mu),}
    where \eqn{\phi_i}{phi_i} are the AR coefficients,
    \eqn{\theta_i}{theta_i} are the MA coefficients,
    \eqn{\Phi_i}{Phi_i} are the seasonal AR coefficients,
    \eqn{\Theta_i}{Theta_i} are the seasonal MA coefficients and
    \eqn{\mu}{mu} is the intercept term.  Note that the \eqn{\mu}{mu}
    entry is required if and only if \code{include.mean} is \code{TRUE}.
    In particular it should not be present if the model is an ARIMA
    model with differencing.

    The entries of the \code{fixed} vector should consist of the
    values at which the user wishes to \dQuote{fix} the corresponding
    coefficient, or \code{NA} if that coefficient should \emph{not} be
    fixed, but estimated.

    The argument \code{transform.pars} will be set to \code{FALSE} if any
    AR parameters are fixed.  A warning will be given if \code{transform.pars}
    is set to (or left at its default) \code{TRUE}.  It may be wise to set
    \code{transform.pars = FALSE} even when fixing MA parameters,
    especially at values that cause the model to be nearly non-invertible.
  }

\item{init}{optional numeric vector of initial parameter
    values.  Missing values will be filled in, by zeroes except for
    regression coefficients.  Values already specified in \code{fixed}
    will be ignored.}

\item{method}{fitting method: maximum likelihood or minimize conditional
sum-of-squares. The default (unless there are missing values) is to use
conditional-sum-of-squares to find starting values, then maximum
likelihood. Note that the random restart algorithm is only implemented
for method = "CSS-ML" or method = "ML", as method = "CSS" is not expected
to maximize the model likelihood.}

\item{n.cond}{only used if fitting by conditional-sum-of-squares: the
    number of initial observations to ignore.  It will be ignored if
    less than the maximum lag of an AR term.}

\item{SSinit}{a string specifying the algorithm to compute the
    state-space initialization of the likelihood; see
    \code{\link[stats]{KalmanLike}} for details.   Can be abbreviated.}

\item{optim.method}{The value passed as the \code{method} argument to
    \code{\link[stats]{optim}}.}

\item{optim.control}{List of control parameters for \code{\link[stats]{optim}}.}

\item{kappa}{the prior variance (as a multiple of the innovations
    variance) for the past observations in a differenced model.  Do not
    reduce this.}
}
\value{
A list of class \code{"Arima"} with components:

  \item{coef}{a vector of AR, MA and regression coefficients, which can
    be extracted by the \code{\link[stats]{coef}} method.}

  \item{sigma2}{the MLE of the innovations variance.}

  \item{var.coef}{the estimated variance matrix of the coefficients
    \code{coef}, which can be extracted by the \code{\link[stats]{vcov}} method.}

  \item{loglik}{the maximized log-likelihood (of the differenced data),
    or the approximation to it used.}

  \item{arma}{A compact form of the specification, as a vector giving
    the number of AR, MA, seasonal AR and seasonal MA coefficients,
    plus the period and the number of non-seasonal and seasonal
    differences.}

  \item{aic}{the AIC value corresponding to the log-likelihood. Only
    valid for \code{method = "ML"} fits.}

  \item{residuals}{the fitted innovations.}

  \item{call}{the matched call.}

  \item{series}{the name of the series \code{x}.}

  \item{code}{the convergence value returned by \code{\link[stats]{optim}}.}

  \item{n.cond}{the number of initial observations not used in the fitting.}

  \item{nobs}{the number of \dQuote{used} observations for the fitting,
    can also be extracted via \code{\link[stats]{nobs}()} and is used by
    \code{\link[stats]{BIC}}.}

  \item{model}{A list representing the Kalman Filter used in the
    fitting.  See \code{\link[stats]{KalmanLike}}.}
}
\description{
Fit an ARIMA model to a univariate time series. This function builds on
the ARIMA model fitting approach used in \code{\link[stats:arima]{stats::arima()}} by fitting
model parameters via a random restart algorithm.
}
\examples{
arima2(LakeHuron, order = c(2, 0, 1), nrestart = 10)
}
