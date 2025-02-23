############################
# Helper Functions (internal)
############################

#' Compute the CDF Value for epn*y (Internal)
#'
#' Computes \eqn{1 - H(epn \times y)}, where \code{H} is a user-supplied function.
#'
#' @param epn Numeric value, typically \eqn{\exp(Z \beta)} or similar.
#' @param y Numeric scalar or vector.
#' @param H A function that takes one argument (numeric) and returns a numeric value.
#'
#' @return A numeric value, \eqn{1 - H(epn \times y)}.
#' @keywords internal
FdG.f <- function(epn, y, H) {
  return(1 - H(epn * y))
}

#' First Derivative Eta1 (Internal)
#'
#' Computes \eqn{-epn \times H1(epn \times y)}, typically used in partial derivatives of log-likelihood.
#'
#' @param epn Numeric value, e.g., \eqn{\exp(Z \beta)}.
#' @param y Numeric scalar or vector.
#' @param H1 A function. The first derivative of \code{H}.
#'
#' @return A numeric value.
#' @keywords internal
eta1G.f <- function(epn, y, H1) {
  return(-epn * H1(epn * y))
}

#' Second Derivative Eta2 (Internal)
#'
#' Computes \eqn{-epn^2 \times H2(epn \times y)}.
#'
#' @param epn Numeric value, e.g., \eqn{\exp(Z \beta)}.
#' @param y Numeric scalar or vector.
#' @param H2 A function. The second derivative of \code{H}.
#'
#' @return A numeric value.
#' @keywords internal
eta2G.f <- function(epn, y, H2) {
  return(-(epn^2) * H2(epn * y))
}

#' Zeta1 Calculation (Internal)
#'
#' Computes \eqn{Fd^{-2} \times \eta1^2 - Fd^{-1} \times \eta2}, used for partial derivatives.
#'
#' @param Fd Numeric value, typically \code{FdG.f(epn, y, H)}.
#' @param eta1 Numeric value, typically \code{eta1G.f(epn, y, H1)}.
#' @param eta2 Numeric value, typically \code{eta2G.f(epn, y, H2)}.
#'
#' @return A numeric value for the zeta1 term.
#' @keywords internal
zeta1.f <- function(Fd, eta1, eta2) {
  return(Fd^(-2) * eta1^2 - Fd^(-1) * eta2)
}

#' Zeta21 Calculation (Internal)
#'
#' Computes \eqn{(Fd2 - Fd1)^{-2} \times \eta1^2 + (Fd2 - Fd1)^{-1} \times \eta2}.
#'
#' @param Fd1 Numeric value (lower bound CDF).
#' @param Fd2 Numeric value (upper bound CDF).
#' @param eta1 Numeric, e.g., \code{eta1G.f(epn, x, H1)}.
#' @param eta2 Numeric, e.g., \code{eta2G.f(epn, x, H2)}.
#'
#' @return A numeric value for the zeta21 term.
#' @keywords internal
zeta21.f <- function(Fd1, Fd2, eta1, eta2) {
  return((Fd2 - Fd1)^(-2) * eta1^2 + (Fd2 - Fd1)^(-1) * eta2)
}

#' Zeta22 Calculation (Internal)
#'
#' Computes \eqn{(Fd2 - Fd1)^{-2} \times \eta1^2 - (Fd2 - Fd1)^{-1} \times \eta2}.
#'
#' @inheritParams zeta21.f
#'
#' @return A numeric value for the zeta22 term.
#' @keywords internal
zeta22.f <- function(Fd1, Fd2, eta1, eta2) {
  return((Fd2 - Fd1)^(-2) * eta1^2 - (Fd2 - Fd1)^(-1) * eta2)
}

#' Zeta3CR Calculation (Internal)
#'
#' Computes \eqn{(ext - Fd)^{-2} \times \eta1^2 + (ext - Fd)^{-1} \times \eta2}.
#'
#' @param Fd Numeric value.
#' @param eta1 Numeric value.
#' @param eta2 Numeric value.
#' @param ext Numeric scalar, an external parameter or bound (e.g., 1).
#'
#' @return A numeric value for the zeta3 term.
#' @keywords internal
zeta3CR.f <- function(Fd, eta1, eta2, ext) {
  return((ext - Fd)^(-2) * eta1^2 + (ext - Fd)^(-1) * eta2)
}

#' First Derivative of Beta (xi1) (Internal)
#'
#' Computes \eqn{-z \times epn \times y \times H1(epn \times y)}, used in partial derivative w.r.t. \code{beta}.
#'
#' @param epn Numeric, e.g. \eqn{\exp(Z \beta)}.
#' @param z Numeric or vector, representing the covariate(s).
#' @param y Numeric, typically \eqn{y[j]} in iterative algorithms.
#' @param H1 A function, the first derivative of \code{H}.
#'
#' @return A numeric value for \code{xi1}.
#' @keywords internal
xi1G.f <- function(epn, z, y, H1) {
  return(-z * epn * y * H1(epn * y))
}

#' Second Derivative of Beta (xi2) (Internal)
#'
#' Computes \eqn{-z z^T \times epn \times y \times (H1(epn \times y) + y \times epn \times H2(epn \times y))}.
#'
#' @param epn Numeric, e.g. \eqn{\exp(Z \beta)}.
#' @param z Numeric or vector, representing the covariate(s) (possibly multiple dimensions).
#' @param y Numeric, typically \eqn{y[j]}.
#' @param H1,H2 Functions, first and second derivatives of \code{H}.
#'
#' @return A matrix or numeric value, the second derivative term.
#' @keywords internal
xi2G.f <- function(epn, z, y, H1, H2) {
  return(-z %*% t(z) * epn * y * (H1(epn * y) + y * epn * H2(epn * y)))
}


############################
# Transformation Functions
############################

#' Log-Transform Functions Constructor (Internal)
#'
#' Creates a list of transformation functions and their derivatives based on the
#' parameter \code{gamma}. Used to transform (and invert) data in estimation routines.
#'
#' @param gamma Numeric scalar. If \code{gamma = 0}, the transform is identity; otherwise,
#'   a log-based transform is used.
#'
#' @return A list of 5 functions:
#'   \describe{
#'     \item{\code{transform[[1]]}}{The main transform.}
#'     \item{\code{transform[[2]]}}{The first derivative of the transform.}
#'     \item{\code{transform[[3]]}}{The second derivative.}
#'     \item{\code{transform[[4]]}}{The third derivative.}
#'     \item{\code{transform[[5]]}}{The fourth derivative.}
#'   }
#' @keywords internal
trans_log <- function(gamma) {
  transform <- list()
  if (gamma == 0) {
    transform[[1]] <- function(x) {
      return(x)
    }
    transform[[2]] <- function(x) {
      return(1)
    }
    transform[[3]] <- function(x) {
      return(0)
    }
    transform[[4]] <- transform[[3]]
    transform[[5]] <- transform[[3]]
  } else {
    transform[[1]] <- function(x) {
      return(1 / gamma * log(1 + gamma * x))
    }
    transform[[2]] <- function(x) {
      return(1 / (1 + gamma * x))
    }
    transform[[3]] <- function(x) {
      return(-gamma / (1 + gamma * x)^2)
    }
    transform[[4]] <- function(x) {
      return(2 * gamma^2 / (1 + gamma * x)^3)
    }
    transform[[5]] <- function(x) {
      return(-6 * gamma^3 / (1 + gamma * x)^4)
    }
  }
  return(transform)
}

#' Inverse Log-Transform Function Constructor (Internal)
#'
#' Creates an inverse function for the transformations constructed by \code{\link{trans_log}}.
#'
#' @param gamma Numeric scalar. If 0, the inverse transform is the identity function.
#'
#' @return A single function, \code{Ginv}, which computes the inverse transform.
#' @keywords internal
trans_inv <- function(gamma) {
  if (gamma == 0) {
    Ginv <- function(y) {
      return(y)
    }
  } else {
    Ginv <- function(y) {
      return(1 / gamma * (exp(gamma * y) - 1))
    }
  }
  return(Ginv)
}


#' Baseline Hazard Function Constructor (Internal)
#'
#' Given functions \code{G}, \code{G1}, and \code{G2}, returns a list of baseline
#' functions \code{H}, \code{H1}, and \code{H2} used in survival or hazard calculations.
#'
#' @param G A function \code{G(x)}, e.g., a cumulative hazard transform.
#' @param G1 A function, the first derivative of \code{G}.
#' @param G2 A function, the second derivative of \code{G}.
#'
#' @return A list with three functions:
#'   \describe{
#'     \item{\code{H(x)}}{Computes \eqn{\exp(-G(x))}.}
#'     \item{\code{H1(x)}}{First derivative of \code{H(x)}.}
#'     \item{\code{H2(x)}}{Second derivative of \code{H(x)}.}
#'   }
#' @keywords internal
H.func <- function(G, G1, G2) {
  H <- function(x) {
    return(exp(-G(x)))
  }

  H1 <- function(x) {
    return(-H(x) * G1(x))
  }

  H2 <- function(x) {
    return(-H1(x) * G1(x) - H(x) * G2(x))
  }

  return(list(H = H, H1 = H1, H2 = H2))
}
