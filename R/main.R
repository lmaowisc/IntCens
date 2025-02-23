##########################
# Internal Helper Functions
##########################

#' Profile Update for y (Internal)
#'
#' Performs an iterative update of the \code{y} vector given current estimates
#' of the baseline function \code{Hf} and regression parameter \code{beta}.
#' Uses the function \code{\link{CSD.f2}} to compute partial derivatives,
#' then applies the convex minorant function \code{\link{CM}}.
#'
#' @param Hf A list of functions \code{(H, H1, H2)} defining the baseline hazard (or similar).
#' @param beta Numeric vector of regression coefficients.
#' @param yini Initial guess for \code{y}.
#' @param delta,gamma Numeric 0/1 vectors indicating left and interval censoring status.
#' @param index A 2-column integer matrix mapping each subject's interval to indices in \code{y}.
#' @param Z A matrix or data frame of covariates.
#' @param eps Convergence tolerance for the iterative algorithm.
#' @param maxiter Maximum number of iterations to perform.
#'
#' @return A numeric vector \code{y} that is the updated estimate after
#'   the iterative procedure.
#' @keywords internal
profile1 <- function(Hf, beta, yini, delta, gamma, index, Z, eps = 1e-4, maxiter = 100) {
  y <- yini
  yp <- 0

  i <- 1
  while (sum(abs(y - yp)) > eps & i <= maxiter) {
    yp <- y
    csd <- CSD.f2(Hf = Hf, beta = beta, y = y, index = index,
                  delta = delta, gamma = gamma, Z = Z)
    G <- csd$G
    Q <- csd$Q
    y <- CM(G = G, Q = Q)
    i <- i + 1
  }

  # Optional message if non-convergence
  if (i == (maxiter + 1)) {
    # message("profile1 did not converge.")
  }

  return(y)
}


#' Newton-Raphson Update for Beta (Internal)
#'
#' Implements a Newton-Raphson step to update the regression parameter \code{beta}
#' using the partial derivatives of the log-likelihood. Internally calls helper
#' functions for computing the score and information matrix components (e.g. \code{\link{xi1G.f}}).
#'
#' @param Hf A list of functions \code{(H, H1, H2)} defining the baseline hazard (or similar).
#' @param beta Numeric vector of regression coefficients (initial/current guess).
#' @param y Numeric vector of baseline-related parameters (often updated by \code{\link{profile1}}).
#' @param delta,gamma Numeric 0/1 vectors indicating left and interval censoring status.
#' @param index A 2-column integer matrix mapping each subject's interval to indices in \code{y}.
#' @param Z A matrix or data frame of covariates.
#' @param eps Convergence tolerance for updating \code{beta}.
#' @param maxiter Maximum number of Newton-Raphson iterations.
#'
#' @return Updated \code{beta} as a numeric vector.
#' @keywords internal
NR2 <- function(Hf, beta, y, delta, gamma, index, Z, eps = 1e-4, maxiter = 20) {

  n <- length(delta)
  p <- length(beta)
  H <- Hf[[1]]
  H1 <- Hf[[2]]
  H2 <- Hf[[3]]

  betap <- rep(0, length(beta)) - 1
  k <- 1

  while (sum(abs(beta - betap)) > eps & k <= maxiter) {
    score <- rep(0, p)
    info <- matrix(0, p, p)
    epn <- exp(Z %*% beta)

    for (i in 1:n) {
      if (delta[i] == 1) {
        j <- index[i, 1]
        x <- y[j]
        Fd <- FdG.f(epn[i], x, H)
        if (Fd <= 1e-10) Fd <- 1e-10

        xi1 <- xi1G.f(epn[i], Z[i, ], x, H1)
        xi2 <- xi2G.f(epn[i], Z[i, ], x, H1, H2)

        score <- score + xi1 / Fd
        info <- info + Fd^(-2) * xi1 %*% t(xi1) - xi2 / Fd

      } else {
        if (gamma[i] == 1) {
          j1 <- index[i, 1]
          j2 <- index[i, 2]
          x1 <- y[j1]
          x2 <- y[j2]
          Fd1 <- FdG.f(epn[i], x1, H)
          Fd2 <- FdG.f(epn[i], x2, H)
          if (Fd2 <= Fd1 + 1e-10) Fd2 <- Fd1 + 1e-10

          xi11 <- xi1G.f(epn[i], Z[i, ], x1, H1)
          xi12 <- xi1G.f(epn[i], Z[i, ], x2, H1)

          xi21 <- xi2G.f(epn[i], Z[i, ], x1, H1, H2)
          xi22 <- xi2G.f(epn[i], Z[i, ], x2, H1, H2)

          score <- score + (Fd2 - Fd1)^(-1) * (xi12 - xi11)
          info <- info + (Fd2 - Fd1)^(-2) * (xi12 - xi11) %*% t(xi12 - xi11) -
            (Fd2 - Fd1)^(-1) * (xi22 - xi21)

        } else {
          j <- index[i, 2]
          x <- y[j]
          Fd <- FdG.f(epn[i], x, H)

          xi1 <- xi1G.f(epn[i], Z[i, ], x, H1)
          xi2 <- xi2G.f(epn[i], Z[i, ], x, H1, H2)

          if (Fd >= 1 - 1e-10) Fd <- 1 - 1e-10

          score <- score - xi1 / (1 - Fd)
          info <- info + (1 - Fd)^(-2) * xi1 %*% t(xi1) + (1 - Fd)^(-1) * xi2
        }
      }
    }

    betap <- beta
    beta <- beta + solve(info) %*% score
    k <- k + 1
  }

  return(beta)
}


#' Log-Likelihood Calculation (Internal)
#'
#' Computes the log-likelihood for interval-censored data given a current estimate
#' of the parameter vector \code{beta} and baseline functions in \code{Hf}.
#'
#' @param H A function representing the transformed baseline function \eqn{H(x) = \exp(-G(x))}.
#' @param beta Numeric vector of regression coefficients.
#' @param y Numeric vector of baseline-related parameters.
#' @param delta,gamma Numeric 0/1 vectors indicating left and interval censoring status.
#' @param index A 2-column integer matrix mapping each subject's interval to indices in \code{y}.
#' @param Z A matrix or data frame of covariates.
#'
#' @return The numeric value of the log-likelihood.
#' @keywords internal
lik.f2 <- function(H, beta, y, delta, gamma, index, Z) {
  n <- length(delta)
  lik <- 0
  epn <- exp(Z %*% beta)

  for (i in 1:n) {
    if (delta[i] == 1) {
      j <- index[i, 1]
      x <- y[j]
      Fd <- FdG.f(epn[i], x, H)
      if (Fd <= 1e-6) Fd <- 1e-6
      lik <- lik + log(Fd)
    } else {
      if (gamma[i] == 1) {
        j1 <- index[i, 1]
        j2 <- index[i, 2]
        x1 <- y[j1]
        x2 <- y[j2]
        Fd1 <- FdG.f(epn[i], x1, H)
        Fd2 <- FdG.f(epn[i], x2, H)
        if (Fd2 <= Fd1 + 1e-6) Fd2 <- Fd1 + 1e-6
        lik <- lik + log(Fd2 - Fd1)
      } else {
        j <- index[i, 2]
        x <- y[j]
        Fd <- FdG.f(epn[i], x, H)
        if (Fd >= 1) Fd <- 1 - 1e-10
        lik <- lik + log(1 - Fd)
      }
    }
  }
  return(lik)
}


#' Higher-Level Profile Update (Internal)
#'
#' Repeatedly calls \code{\link{profile1}} to update \code{y}, converging to a solution
#' based on the current \code{beta} and baseline hazard functions \code{Hf}.
#'
#' @inheritParams profile1
#'
#' @return The final log-likelihood (numeric) after convergence.
#' @keywords internal
profile2 <- function(Hf, y, beta, delta, gamma, index, Z, eps = 1e-4, maxiter = 100) {
  yp <- 0
  j <- 1

  while (sum(abs(y - yp)) > eps & j <= maxiter) {
    yp <- y
    y <- profile1(Hf = Hf, beta = beta, yini = y,
                  delta = delta, gamma = gamma, index = index, Z = Z,
                  eps = 1e-3, maxiter = 100)
    j <- j + 1
  }

  lik <- lik.f2(Hf[[1]], beta, y, delta, gamma, index, Z)
  return(lik)
}


#' Quadratic Approximation for Profile Likelihood (Internal)
#'
#' Approximates the information matrix via finite differences by perturbing
#' \code{beta} in each dimension, re-profiling over \code{y}, and computing
#' the difference in the log-likelihood.
#'
#' @param n Sample size, used to scale perturbations \code{ptb = 1/sqrt(n)}.
#' @param Hf A list of functions \code{(H, H1, H2)} defining the baseline hazard.
#' @param beta Numeric vector of regression coefficients.
#' @param y Numeric vector of baseline parameters (updated by \code{\link{profile1}}).
#' @inheritParams profile1
#'
#' @return A \code{p x p} matrix, an approximation to the negative Hessian (information).
#' @keywords internal
quad.profile2 <- function(n, Hf, beta, y, delta, gamma, index, Z,
                          eps = 1e-4, maxiter = 100) {
  p <- length(beta)
  ptb <- 1 / sqrt(n)

  # Likelihood at beta
  y <- profile1(Hf = Hf, yini = y, beta = beta,
                delta = delta, gamma = gamma, index = index, Z = Z,
                eps = eps, maxiter = maxiter)
  lik0 <- lik.f2(Hf[[1]], beta, y, delta, gamma, index, Z)

  info <- matrix(NA, p, p)

  # Diagonal entries
  for (i in 1:p) {
    beta.tmp <- beta
    beta.tmp[i] <- beta.tmp[i] + ptb
    y.tmp <- profile1(Hf = Hf, yini = y, beta = beta.tmp,
                      delta = delta, gamma = gamma, index = index, Z = Z,
                      eps = eps, maxiter = maxiter)
    lik.tmp <- lik.f2(Hf[[1]], beta.tmp, y.tmp, delta, gamma, index, Z)
    info[i, i] <- -2 * (lik.tmp - lik0) / ptb^2
  }

  # Off-diagonal entries
  if (p > 1) {
    for (i in 2:p) {
      for (j in 1:(i - 1)) {
        beta.tmp <- beta
        beta.tmp[i] <- beta.tmp[i] + ptb
        beta.tmp[j] <- beta.tmp[j] + ptb

        y.tmp <- profile1(Hf = Hf, yini = y, beta = beta.tmp,
                          delta = delta, gamma = gamma, index = index, Z = Z,
                          eps = eps, maxiter = maxiter)
        lik.tmp <- lik.f2(Hf[[1]], beta.tmp, y.tmp, delta, gamma, index, Z)

        val <- (-2 * (lik.tmp - lik0) / ptb^2 - info[i, i] - info[j, j]) / 2
        info[i, j] <- val
        info[j, i] <- val
      }
    }
  }

  return(info)
}


#' Internal Composite Function for Interval Censoring (Internal)
#'
#' Provides an internal composite function that runs a full iterative procedure
#' for fitting an interval-censored model, supporting different transformations.
#' Typically called by the user-facing \code{\link{icsurvfit}} function.
#'
#' @param delta,gamma,U,V Numeric vectors representing censoring indicators and intervals.
#' @param Z Matrix or data frame of covariates.
#' @param model Character string, one of \code{"NP"}, \code{"PO"}, or \code{"PH"}.
#' @param maxiter Maximum number of iterations.
#'
#' @return A list with elements \code{beta}, \code{y}, \code{t}, \code{var}, and convergence info.
#' @keywords internal
#' @export
ICSurvICM <- function(delta, gamma, U, V, Z, model, maxiter=500){
  if (model!="PO"&model!="PH"){

    G.l=trans_log(0)
    Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

    #dummy data
    n=length(delta)
    Z=as.matrix(rep(0,n))

    obj=order_rank(delta,gamma,U,V,Z)

    data=obj$data
    t=obj$t
    index=obj$index
    delta=data[,1]
    gamma=data[,2]
    U=data[,3]
    V=data[,4]
    Z=as.matrix(data[,5:ncol(data)])


    mp=length(t)
    y=(1:mp)/mp



    yp=0


    beta=rep(0,ncol(Z))
    betap=beta-1

    eps=1e-3

    j=1

    while(sum(abs(y-yp))>eps&j<=maxiter){
      yp=y

      y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)

      j=j+1
    }
    if (j==maxiter+1){
      print(conv=F)
    }else{
      conv=T
      var=NULL
    }
  }else{
    if (model=="PO"){
      r=1
    }else{
      r=0
    }

    G.l=trans_log(r)
    Hf=H.func(G.l[[1]],G.l[[2]],G.l[[3]])

    #   delta=data[,1]
    #   gamma=data[,2]
    #   U=data[,3]
    #   V=data[,4]
    #   Z=as.matrix(data[,5:ncol(data)])

    Z=as.matrix(Z)
    obj=order_rank(delta,gamma,U,V,Z)

    data=obj$data
    t=obj$t
    index=obj$index
    delta=data[,1]
    gamma=data[,2]
    U=data[,3]
    V=data[,4]
    Z=as.matrix(data[,5:ncol(data)])
    n=length(delta)

    mp=length(t)
    y=(1:mp)/mp



    yp=0


    beta=rep(0,ncol(Z))
    betap=beta-1

    eps=1e-4

    j=1

    while(sum(abs(beta-betap))>eps&j<=maxiter){
      yp=y

      y=profile1(Hf=Hf,beta=beta,yini=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-3,maxiter=200)
      betap=beta
      beta=NR2(Hf=Hf,beta=beta,y=y,delta=delta,gamma=gamma,index=index,Z=Z,eps=1e-4,maxiter=20)
      j=j+1
    }
    if (j==maxiter+1){
      print("conv=F")
    }else{
      conv=T
      var=solve(quad.profile2(n,Hf,beta,y,delta,gamma,index,Z,eps=1e-3,maxiter=maxiter))
    }
  }
  obj=list(beta=beta,y=y,t=t,var=var,j=j,model=model,conv=conv,varnames=colnames(Z),call=match.call())
  class(obj)="ICSurv"
  return(obj)
}


##########################
# Main Exported Function
##########################

#' Fit Interval-Censored Model
#'
#' Fits an interval-censored model of various types (\code{NP}, \code{PO}, or \code{PH})
#' by iterating over baseline and regression parameter updates.
#'
#' @param L Left endpoints (may contain \code{NA}, \code{-Inf}).
#' @param R Right endpoints (may contain \code{NA}, \code{Inf}).
#' @param Z A matrix or data frame of covariates.
#' @param model One of \code{"NP"}, \code{"PO"}, or \code{"PH"}. Defaults to \code{"NP"}
#'   for nonparametric (no covariates).
#' @param maxiter Maximum number of iterations for the fitting algorithm.
#'
#' @details
#' \describe{
#'   \item{\code{NP}}{Nonparametric (no covariate effects).}
#'   \item{\code{PO}}{Proportional odds.}
#'   \item{\code{PH}}{Proportional hazards.}
#' }
#' Internally, this function sets up the data, identifies which observations are
#' left-, right-, or interval-censored, and then calls various helper routines
#' (e.g., \code{\link{profile1}}, \code{\link{NR2}}) to perform the iterative fitting.
#'
#' @return A list (class \code{"icsurvfit"}) with components:
#' \describe{
#'   \item{\code{beta}}{Estimated regression coefficients (if \code{model="PO"} or \code{"PH"}).}
#'   \item{\code{y}}{Baseline function estimates at each time index.}
#'   \item{\code{t}}{Ordered time points at which \code{y} is evaluated.}
#'   \item{\code{var}}{Estimated variance-covariance matrix of \code{beta} (if applicable).}
#'   \item{\code{j}}{Number of iterations used.}
#'   \item{\code{model}}{The specified model type.}
#'   \item{\code{conv}}{Logical indicating whether the algorithm converged.}
#'   \item{\code{varnames}}{Covariate names.}
#'   \item{\code{call}}{Matched call.}
#' }
#'
#' @examples
#' # Simple example
#' data(bcos)
#'
#' n <- nrow(bcos) # sample size
#' L <- bcos$left + rnorm(n, 0, 0.0001) # add random noise
#' R <- bcos$right + rnorm(n, 0, 0.0001)
#' Z <- as.numeric(bcos$treatment == "RadChem")
#'
#' # Nonparametric
#' obj_np <- icsurvfit(L, R)
#'
#' # Cox model (PH)
#' obj_ph <- icsurvfit(L, R, Z, model = "PH")
#' obj_ph
#' #> Call:
#' #>   icsurvfit(L = L, R = R, Z = Z, model = "PH")
#' #>
#' #> NPMLE of proportional hazards for interval-censored data:
#' #>
#' #>   ICM algorithm converges in 20 iterations.
#' #>
#' #> Maximum Likelihood Estimates for Regression parameters:
#' #>
#' #>   Estimate  StdErr z.value  p.value
#' #> [1,]  0.78883 0.29344  2.6882 0.007183 **
#' #>   ---
#' #>   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#'
#' # Proportional odds model (PO)
#' obj_po <- icsurvfit(L, R, Z, model = "PO")
#' obj_po
#' #> Call:
#' #>   icsurvfit(L = L, R = R, Z = Z, model = "PO")
#' #>
#' #> NPMLE of proportional odds for interval-censored data:
#' #>
#' #>   ICM algorithm converges in 20 iterations.
#' #>
#' #> Maximum Likelihood Estimates for Regression parameters:
#' #>
#' #>   Estimate  StdErr z.value p.value
#' #> [1,]  0.89529 0.41137  2.1764 0.02953 *
#' #>   ---
#' #>   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#' @export
icsurvfit <- function(L, R, Z, model = "NP", maxiter = 500) {
  # ... (the body of your function exactly as provided) ...
  # Keep all lines the same, except for doc block above
  # (Shown for completeness below)

  # delta, gamma, U, V
  # Convert L & R
  L[is.na(L)] <- 0
  L[L <= 0] <- 0
  L[L == -Inf] <- 0
  R[is.na(R)] <- Inf

  if (any(L > R)){
    stop("L must be less than or equal to R")
  }

  fl <- ifelse(L == 0, 1, ifelse(L > 0 & R < Inf, 2, 3))
  delta <- (fl == 1) + 0
  gamma <- (fl == 2) + 0
  U <- ifelse(fl == 3, NA, ifelse(fl == 2, L, R))
  V <- ifelse(fl == 1, NA, ifelse(fl == 2, R, L))

  if (model != "PO" & model != "PH") {
    G.l <- trans_log(0)
    Hf <- H.func(G.l[[1]], G.l[[2]], G.l[[3]])

    n <- length(delta)
    Z <- as.matrix(rep(0, n))

    obj <- order_rank(delta, gamma, U, V, Z)

    data <- obj$data
    t <- obj$t
    index <- obj$index
    delta <- data[, 1]
    gamma <- data[, 2]
    U <- data[, 3]
    V <- data[, 4]
    Z <- as.matrix(data[, 5:ncol(data)])

    mp <- length(t)
    y <- (1:mp) / mp

    yp <- 0
    beta <- rep(0, ncol(Z))
    betap <- beta - 1
    eps <- 1e-3
    j <- 1

    while (sum(abs(y - yp)) > eps & j <= maxiter) {
      yp <- y
      y <- profile1(Hf = Hf, beta = beta, yini = y,
                    delta = delta, gamma = gamma, index = index,
                    Z = Z, eps = 1e-3, maxiter = 200)
      j <- j + 1
    }
    if (j == maxiter + 1) {
      print(conv = FALSE)
    } else {
      conv <- TRUE
      var <- NULL
    }

  } else {
    if (model == "PO") {
      r <- 1
    } else {
      r <- 0
    }
    G.l <- trans_log(r)
    Hf <- H.func(G.l[[1]], G.l[[2]], G.l[[3]])

    Z <- as.matrix(Z)
    obj <- order_rank(delta, gamma, U, V, Z)

    data <- obj$data
    t <- obj$t
    index <- obj$index
    delta <- data[, 1]
    gamma <- data[, 2]
    U <- data[, 3]
    V <- data[, 4]
    Z <- as.matrix(data[, 5:ncol(data)])
    n <- length(delta)

    mp <- length(t)
    y <- (1:mp) / mp

    yp <- 0
    beta <- rep(0, ncol(Z))
    betap <- beta - 1
    eps <- 1e-4
    j <- 1

    while (sum(abs(beta - betap)) > eps & j <= maxiter) {
      yp <- y
      y <- profile1(Hf = Hf, beta = beta, yini = y,
                    delta = delta, gamma = gamma, index = index,
                    Z = Z, eps = 1e-3, maxiter = 200)
      betap <- beta
      beta <- NR2(Hf = Hf, beta = beta, y = y,
                  delta = delta, gamma = gamma, index = index,
                  Z = Z, eps = 1e-4, maxiter = 20)
      j <- j + 1
    }

    if (j == maxiter + 1) {
      print("conv = F")
      conv <- FALSE
    } else {
      conv <- TRUE
      var <- solve(quad.profile2(n, Hf, beta, y, delta, gamma, index, Z,
                                 eps = 1e-3, maxiter = maxiter))
    }
  }

  obj <- list(
    beta = beta,
    y = y,
    t = t,
    var = var,
    j = j,
    model = model,
    conv = conv,
    varnames = colnames(Z),
    call = match.call()
  )
  class(obj) <- "icsurvfit"
  return(obj)
}
