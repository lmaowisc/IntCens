#' Compute Convex Minorant
#'
#' Computes a piecewise-linear convex minorant of two numeric vectors \code{G}
#' and \code{Q} by iteratively identifying minimum slopes and assigning them to
#' the output vector \code{y}. This method is sometimes used in isotonic or
#' shape-constrained regression contexts.
#'
#' @param G Numeric vector of \eqn{x}-coordinates (must be the same length as \code{Q}).
#' @param Q Numeric vector of \eqn{y}-coordinates (must be the same length as \code{G}).
#'
#' @return A numeric vector \code{y} of length \code{length(G)}, giving the
#'   slopes or values associated with the convex minorant.
#'
#' @details
#' Internally, the function calculates slopes for segments between points
#' in \code{G} and picks the minimum slope iteratively, ensuring the resulting
#' piecewise-linear function is convex.
#'
#'
CM <- function(G, Q) {
  m <- length(G) + 1
  y <- rep(NA, m - 1)
  Q <- c(0, Q)
  G <- c(0, G)
  j <- 1
  while (j < m) {
    slopes <- rep(NA, m - j)
    for (k in (j + 1):m) {
      slopes[k - j] <- (Q[k] - Q[j]) / (G[k] - G[j])
    }
    j1 <- which.min(slopes) + j
    y[j:(j1 - 1)] <- min(slopes)
    j <- j1
  }
  return(y)
}


#' Order and Rank Interval-Censored Data
#'
#' Processes interval-censored data (and corresponding indicators) to produce
#' an ordered set of exam times and flags indicating censoring status. It also
#' adjusts certain \code{delta} and \code{gamma} values based on constraints
#' so that the first observation is always an event time, and the last
#' observation is always right-censored.
#'
#' @param delta Numeric (0/1). Event indicator for left endpoint \code{U}.
#' @param gamma Numeric (0/1). Event indicator for right endpoint \code{V}.
#' @param U Numeric vector of left endpoints of censoring intervals (may contain \code{NA}).
#' @param V Numeric vector of right endpoints of censoring intervals (may contain \code{NA}).
#' @param Z A matrix (or data frame) of covariates, with one row per subject.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{t}}{A condensed, ordered set of exam times.}
#'   \item{\code{index}}{A matrix of dimension \code{n x 2}, giving the ranks (indices)
#'        for \code{U} and \code{V} after ordering.}
#'   \item{\code{data}}{A matrix (or data frame) containing the updated
#'        \code{delta}, \code{gamma}, \code{U}, \code{V}, and \code{Z} values.}
#' }
#'
#' @details
#' This function ensures that the earliest recorded time in \code{t} corresponds
#' to a left endpoint that is truly an event, while the latest recorded time
#' is right-censored. It also removes observations that are no longer valid
#' after these adjustments. The vector \code{flag} is used internally to mark
#' observation types (1 = event, 2 = left-censored, 3 = both, 4 = right-censored).
#'
#' @examples
#' # (Not run) Example usage (requires suitable data)
#' # delta <- c(1, 0, 0)
#' # gamma <- c(0, 1, 0)
#' # U <- c(2, 2, 3)
#' # V <- c(NA, 4, 5)
#' # Z <- matrix(rnorm(9), nrow = 3)
#' # order_rank(delta, gamma, U, V, Z)
#'
order_rank <- function(delta, gamma, U, V, Z) {
  n <- length(delta)

  # condense the exam times
  t <- NULL
  flag <- NULL
  obs <- NULL

  for (i in 1:n) {
    if (delta[i] == 1) {
      t <- c(t, U[i])
      flag <- c(flag, 1)
      obs <- c(obs, i)
    } else {
      if (gamma[i] == 1) {
        t <- c(t, V[i])
        flag <- c(flag, 3)
        obs <- c(obs, i)
        t <- c(t, U[i])
        flag <- c(flag, 2)
        obs <- c(obs, i)
      } else {
        t <- c(t, V[i])
        flag <- c(flag, 4)
        obs <- c(obs, i)
      }
    }
  }
  o <- order(t)
  t <- t[o]
  flag <- flag[o]
  obs <- obs[o]

  # Adjust to ensure T_(1) is an event
  i <- 1
  del <- NULL
  while (flag[i] != 1) {
    j <- obs[i]
    if (flag[i] == 2) {
      delta[j] <- 1
      gamma[j] <- 0
      U[j] <- V[j]
      V[j] <- NA
      flag[i] <- 1
    } else {
      if (flag[i] == 4) {
        del <- c(del, j)
        flag[i] <- NA
        t[i] <- NA
      }
    }
    i <- i + 1
  }

  # Adjust to ensure T_(m) is right-censored
  i <- length(flag)
  while (flag[i] != 4) {
    j <- obs[i]
    if (flag[i] == 3) {
      gamma[j] <- 0
      V[j] <- U[j]
      U[j] <- NA
      flag[i] <- 4
    } else {
      if (flag[i] == 1) {
        del <- c(del, j)
        flag[i] <- NA
        t[i] <- NA
      }
    }
    i <- i - 1
  }

  # Remove invalid entries
  if (!is.null(del)) {
    U <- U[-del]
    V <- V[-del]
    Z <- Z[-del, ]
    delta <- delta[-del]
    gamma <- gamma[-del]
  }
  flag <- flag[!is.na(flag)]
  t <- t[!is.na(t)]
  flag <- as.numeric(flag == 1 | flag == 3)

  # Aggregate
  m <- length(t)
  n <- length(delta)
  key <- rep(NA, m)
  key[1] <- 1
  aggr <- rep(NA, m)
  aggr[1] <- 1
  for (i in 2:m) {
    jump <- as.numeric(flag[i - 1] == 0 & flag[i] == 1)
    key[i] <- key[i - 1] + jump
    aggr[i] <- jump
  }

  index <- matrix(0, n, 2)
  for (i in 1:n) {
    if (delta[i] == 1) {
      index[i, 1] <- key[sum(U[i] >= t - 1e-10)]
    } else {
      if (gamma[i] == 1) {
        index[i, 1] <- key[sum(U[i] >= t - 1e-10)]
        index[i, 2] <- key[sum(V[i] >= t - 1e-10)]
      } else {
        index[i, 2] <- key[sum(V[i] >= t - 1e-10)]
      }
    }
  }

  t <- t[aggr == 1]

  return(list(
    t = t,
    index = index,
    data = cbind(delta, gamma, U, V, Z)
  ))
}


#' Compute Score Components for Interval-Censored Data
#'
#' Computes partial derivatives or score components for interval-censored data,
#' given current estimates of baseline functions and a regression parameter
#' vector \code{beta}. This function is typically part of an iterative algorithm
#' for semiparametric estimation.
#'
#' @param Hf A list containing baseline functions or related objects (e.g., \code{H}, \code{H1}, \code{H2}).
#' @param beta A numeric vector of regression coefficients for the covariates in \code{Z}.
#' @param y A numeric vector (often the current estimate of the cumulative baseline
#'   hazard or a similar function) of length \code{p}.
#' @param index An \code{n x 2} matrix of indices referencing positions in \code{y} for each subject.
#' @param delta Numeric (0/1) event indicator for the left endpoint \code{U}.
#' @param gamma Numeric (0/1) event indicator for the right endpoint \code{V}.
#' @param Z A matrix or data frame of covariates, with \code{n} rows (one per subject).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{G}}{A numeric vector (length \code{p}) of cumulative increments related to \code{gamma}.}
#'   \item{\code{Q}}{A numeric vector (length \code{p}) of combined increments (e.g., \code{dW + y[i] * dG}).}
#' }
#'
#' @details
#' This function uses auxiliary routines (e.g., \code{FdG.f()}, \code{eta1G.f()},
#' \code{eta2G.f()}, \code{zeta1.f()}, etc.) to compute partial derivatives of
#' the log-likelihood or estimating equations with respect to \code{y}.
#'
#' @examples
#' # (Not run) Example with placeholder data:
#' # Hf <- list(H = NULL, H1 = NULL, H2 = NULL)
#' # beta <- c(0.1, -0.2)
#' # y <- rep(0.05, 5)
#' # index <- matrix(c(1, 2, 2, 3), nrow = 2, byrow = TRUE)
#' # delta <- c(1, 0)
#' # gamma <- c(0, 1)
#' # Z <- matrix(c(1, 0, 0, 1), nrow = 2)
#' # CSD.f2(Hf, beta, y, index, delta, gamma, Z)
#'
CSD.f2 <- function(Hf, beta, y, index, delta, gamma, Z) {

  H <- Hf[[1]]
  H1 <- Hf[[2]]
  H2 <- Hf[[3]]

  p <- length(y)
  n <- length(delta)
  dW <- rep(0, p)
  dG <- rep(0, p)
  dQ <- rep(0, p)
  W <- rep(0, p)
  G <- rep(0, p)
  Q <- rep(0, p)

  epn <- exp(Z %*% beta)

  for (i in 1:n) {

    # 1) Left endpoint contribution
    if (index[i, 1] != 0) {
      if (delta[i] == 1) {
        j <- index[i, 1]
        x <- y[j]
        Fd <- FdG.f(epn[i], x, H)
        if (Fd <= 1e-6) Fd <- 1e-6

        eta1 <- eta1G.f(epn[i], x, H1)
        eta2 <- eta2G.f(epn[i], x, H2)
        zeta1 <- zeta1.f(Fd, eta1, eta2)

        dW[j] <- dW[j] + eta1 / Fd
        dG[j] <- dG[j] + zeta1
      } else {
        if (gamma[i] == 1) {
          j1 <- index[i, 1]
          j2 <- index[i, 2]
          x1 <- y[j1]
          x2 <- y[j2]
          Fd1 <- FdG.f(epn[i], x1, H)
          Fd2 <- FdG.f(epn[i], x2, H)
          if (Fd2 <= Fd1 + 1e-6) Fd2 <- Fd1 + 1e-6

          eta1 <- eta1G.f(epn[i], x1, H1)
          eta2 <- eta2G.f(epn[i], x1, H2)
          zeta2 <- zeta21.f(Fd1, Fd2, eta1, eta2)
          dW[j1] <- dW[j1] - (Fd2 - Fd1)^(-1) * eta1
          dG[j1] <- dG[j1] + zeta2
        }
      }
    }

    # 2) Right endpoint contribution
    if (index[i, 2] != 0) {
      if (gamma[i] == 1) {
        j1 <- index[i, 1]
        j2 <- index[i, 2]
        x1 <- y[j1]
        x2 <- y[j2]
        Fd1 <- FdG.f(epn[i], x1, H)
        Fd2 <- FdG.f(epn[i], x2, H)
        if (Fd2 <= Fd1 + 1e-6) Fd2 <- Fd1 + 1e-6

        eta1 <- eta1G.f(epn[i], x2, H1)
        eta2 <- eta2G.f(epn[i], x2, H2)
        zeta2 <- zeta22.f(Fd1, Fd2, eta1, eta2)
        dW[j2] <- dW[j2] + (Fd2 - Fd1)^(-1) * eta1
        dG[j2] <- dG[j2] + zeta2
      } else {
        if (gamma[i] == 0) {
          j <- index[i, 2]
          x <- y[j]
          Fd <- FdG.f(epn[i], x, H)
          eta1 <- eta1G.f(epn[i], x, H1)
          eta2 <- eta2G.f(epn[i], x, H2)
          if (Fd >= 1 - 1e-6) Fd <- 1 - 1e-6

          zeta3 <- zeta3CR.f(Fd, eta1, eta2, 1)

          dW[j] <- dW[j] - eta1 / (1 - Fd)
          dG[j] <- dG[j] + zeta3
        }
      }
    }
  }

  # Combine increments
  for (i in 1:p) {
    if (is.na(dW[i])) dW[i] <- 0
    dQ[i] <- dW[i] + y[i] * dG[i]
    G[i] <- sum(dG[1:i])
    Q[i] <- sum(dQ[1:i])
    if (Q[i] < 0) Q[i] <- 0
  }

  return(list(G = G, Q = Q))
}
