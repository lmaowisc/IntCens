#' Print Method for icsurvfit Objects
#'
#' Prints a summary of the fitted model, including the call, convergence status,
#' estimated survival function (for \code{model="NP"}), or estimated regression
#' coefficients (for \code{model="PO"} or \code{"PH"}).
#'
#' @param object An object of class \code{"icsurvfit"} returned by
#'   \code{\link{icsurvfit}}.
#' @param ... Further arguments passed to or from other methods. Not currently used.
#' @importFrom stats pnorm printCoefmat
#' @method print icsurvfit
#' @export
print.icsurvfit <- function(object, ...) {
  cat("\nCall:\n")
  print(object$call)
  cat("\n")
  if (object$conv) {
    beta <- object$beta
    var <- object$var
    y <- object$y
    t <- object$t
    model <- object$model
    varnames <- object$varnames

    if (model != "PO" & model != "PH") {
      result <- cbind(t, exp(-y))
      rownames(result) <- 1:length(y)
      colnames(result) <- c("Time", "Survival rate")

      cat("NPMLE for interval-censored data:\n")
      cat("ICM algorithm converges in", object$j, "iterations.\n\n")

      cat("Estimated survival function\n")
      cat("---------------------------\n")
      print(result)
      cat("\n")
    } else {
      p <- length(beta)
      se <- sqrt(diag(var))
      pv <- 2 * (1 - pnorm(abs(beta / se)))

      if (model == "PO") {
        mod <- "proportional odds"
      } else {
        mod <- "proportional hazards"
      }

      cat("NPMLE of", mod, "for interval-censored data:\n\n")
      cat("ICM algorithm converges in", object$j, "iterations.\n\n")

      table <- cbind(
        Estimate = as.vector(beta),
        StdErr   = se,
        z.value  = as.vector(beta / se),
        p.value  = as.vector(pv)
      )
      rownames(table) <- varnames

      cat("Maximum Likelihood Estimates for Regression parameters:\n\n")
      printCoefmat(table, P.values = TRUE, has.Pvalue = TRUE)
      cat("\n")
    }
  } else {
    cat("Convergence failure.\n")
  }
}


#' Plot Method for icsurvfit Objects
#'
#' Plots the estimated survival function from an \code{icsurvfit} object.
#' If \code{model="NP"}, the nonparametric estimate is plotted. For
#' \code{model="PO"} or \code{"PH"}, the survival curve for a specific covariate
#' value \code{z} can be overlaid. By default, it plots a step function.
#'
#' @param object An object of class \code{"icsurvfit"} returned by
#'   \code{\link{icsurvfit}}.
#' @param z Numeric (vector of) covariate values for which to plot the survival
#'   curve. If \code{NULL} (default) and \code{model} is nonparametric,
#'   plots the baseline curve.
#' @param xlab,ylab Axis labels.
#' @param lty Line type. Default is \code{1}.
#' @param frame.plot Logical indicating whether to draw a box around the plot.
#' @param add Logical indicating whether to add the curve to an existing plot.
#'   If \code{FALSE} (default), a new plot is created.
#' @param ylim Numeric vector of length 2 specifying the \eqn{y}-axis range.
#' @param ... Additional graphical parameters passed to \code{\link{plot}} or
#'   \code{\link{lines}}.
#' @importFrom stats stepfun
#' @importFrom graphics lines
#' @method plot icsurvfit
#' @export
plot.icsurvfit <- function(object, z = NULL, xlab = "Time", ylab = "Survival rate",
                           lty = 1, frame.plot = FALSE, add = FALSE,
                           ylim = c(0, 1), ...) {

  beta <- object$beta
  var <- object$var
  y <- object$y
  t <- object$t
  model <- object$model
  varnames <- object$varnames

  if (model != "PO" & model != "PH") {
    # Nonparametric case
    if (!add) {
      plot(stepfun(t, c(1, exp(-y))),
           do.points = FALSE, ylim = ylim,
           lty = lty, xlab = xlab, ylab = ylab,
           frame.plot = frame.plot, ...
      )
    } else {
      lines(stepfun(t, c(1, exp(-y))),
            do.points = FALSE, lty = lty, ...
      )
    }
  } else {
    # Parametric (Proportional Odds or Proportional Hazards)
    if (is.null(z)) z <- 0  # or some default if user doesn't supply

    if (model == "PO") {
      # PO => S(t | z) = (1 + exp(sum(beta*z)) * y)^(-1)
      curveVals <- (1 + exp(sum(beta * z)) * y)^(-1)
      if (!add) {
        plot(stepfun(t, c(1, curveVals)), do.points = FALSE, ylim = ylim,
             lty = lty, xlab = xlab, ylab = ylab,
             frame.plot = frame.plot, ...
        )
      } else {
        lines(stepfun(t, c(1, curveVals)),
              do.points = FALSE, lty = lty, ...
        )
      }
    } else {
      # PH => S(t | z) = exp(- exp(sum(beta*z)) * y)
      curveVals <- exp(-exp(sum(beta * z)) * y)
      if (!add) {
        plot(stepfun(t, c(1, curveVals)), do.points = FALSE, ylim = ylim,
             lty = lty, xlab = xlab, ylab = ylab,
             frame.plot = frame.plot, ...
        )
      } else {
        lines(stepfun(t, c(1, curveVals)),
              do.points = FALSE, lty = lty, ...
        )
      }
    }
  }
}


