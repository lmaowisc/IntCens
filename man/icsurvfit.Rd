% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{icsurvfit}
\alias{icsurvfit}
\title{Fit Interval-Censored Model}
\usage{
icsurvfit(L, R, Z, model = "NP", maxiter = 500)
}
\arguments{
\item{L}{Left endpoints (may contain \code{NA}, \code{-Inf}).}

\item{R}{Right endpoints (may contain \code{NA}, \code{Inf}).}

\item{Z}{A matrix or data frame of covariates.}

\item{model}{One of \code{"NP"}, \code{"PO"}, or \code{"PH"}. Defaults to \code{"NP"}
for nonparametric (no covariates).}

\item{maxiter}{Maximum number of iterations for the fitting algorithm.}
}
\value{
A list (class \code{"icsurvfit"}) with components:
\describe{
  \item{\code{beta}}{Estimated regression coefficients (if \code{model="PO"} or \code{"PH"}).}
  \item{\code{y}}{Baseline function estimates at each time index.}
  \item{\code{t}}{Ordered time points at which \code{y} is evaluated.}
  \item{\code{var}}{Estimated variance-covariance matrix of \code{beta} (if applicable).}
  \item{\code{j}}{Number of iterations used.}
  \item{\code{model}}{The specified model type.}
  \item{\code{conv}}{Logical indicating whether the algorithm converged.}
  \item{\code{varnames}}{Covariate names.}
  \item{\code{call}}{Matched call.}
}
}
\description{
Fits an interval-censored model of various types (\code{NP}, \code{PO}, or \code{PH})
by iterating over baseline and regression parameter updates.
}
\details{
\describe{
  \item{\code{NP}}{Nonparametric (no covariate effects).}
  \item{\code{PO}}{Proportional odds.}
  \item{\code{PH}}{Proportional hazards.}
}
Internally, this function sets up the data, identifies which observations are
left-, right-, or interval-censored, and then calls various helper routines
(e.g., \code{\link{profile1}}, \code{\link{NR2}}) to perform the iterative fitting.
}
\examples{
# Simple example
data(bcos)

n <- nrow(bcos) # sample size
L <- bcos$left + rnorm(n, 0, 0.0001) # add random noise
R <- bcos$right + rnorm(n, 0, 0.0001)
Z <- as.numeric(bcos$treatment == "RadChem")

# Nonparametric
obj_np <- icsurvfit(L, R)

# Cox model (PH)
obj_ph <- icsurvfit(L, R, Z, model = "PH")
obj_ph
#> Call:
#>   icsurvfit(L = L, R = R, Z = Z, model = "PH")
#>
#> NPMLE of proportional hazards for interval-censored data:
#>
#>   ICM algorithm converges in 20 iterations.
#>
#> Maximum Likelihood Estimates for Regression parameters:
#>
#>   Estimate  StdErr z.value  p.value
#> [1,]  0.78883 0.29344  2.6882 0.007183 **
#>   ---
#>   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Proportional odds model (PO)
obj_po <- icsurvfit(L, R, Z, model = "PO")
obj_po
#> Call:
#>   icsurvfit(L = L, R = R, Z = Z, model = "PO")
#>
#> NPMLE of proportional odds for interval-censored data:
#>
#>   ICM algorithm converges in 20 iterations.
#>
#> Maximum Likelihood Estimates for Regression parameters:
#>
#>   Estimate  StdErr z.value p.value
#> [1,]  0.89529 0.41137  2.1764 0.02953 *
#>   ---
#>   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
