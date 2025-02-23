#' Breast Cosmesis Data
#'
#' A dataset containing interval-censored observations and a treatment indicator
#' (radiation vs. radiation plus chemotherapy). The variables are as follows:
#'
#' @format A data frame with 94 rows and 3 variables:
#' \describe{
#'   \item{left}{Left endpoints of the censoring intervals (numeric).}
#'   \item{right}{Right endpoints of the censoring intervals (numeric).}
#'   \item{treatment}{A factor with levels \code{"Rad"} (radiation only)
#'       or \code{"RadChem"} (radiation + chemotherapy).}
#' }
#'
#' @docType data
#' @usage data(bcos)
#' @keywords datasets
#'
#' @details
#' These data are often used to illustrate methods for interval-censored survival analysis.
#' Observations with \code{Inf} in the right endpoint represent right-censored data.
#'
#' @examples
#' data(bcos)
#' head(bcos)
"bcos"
