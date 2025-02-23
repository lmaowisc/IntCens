

# Preambles ---------------------------------------------------------------

# install.packages("devtools")
# install.packages("usethis")
# install.packages("rmarkdown")
# install.packages("fs")

library(devtools)
library(roxygen2)
library(knitr)
library(rmarkdown)
library(usethis)
# library(fs)


use_git_config(
  user.name = "lmaowisc",
  user.email = "lmao@biostat.wisc.edu"
)
use_git()

# usethis:::use_devtools()


# Create, edit, and load code ---------------------------------------------
## create R file in "./R"
use_r("wr_helpers.R")
## Test your function in the new package
### devtools::load_all()
### Ctrl+Shift+L
load_all()



# Check package -----------------------------------------------------------
## 3 types of messages
## • ERRORs: Severe problems - always fix.
## • WARNINGs: Problems that you should fix, and must fix if you’re planning to
## submit to CRAN.
## • NOTEs: Mild problems or, in a few cases, just an observation.
## • When submitting to CRAN, try to eliminate all NOTEs.
check()



# Licenses ----------------------------------------------------------------

use_ccby_license()


# The DESCRIPTION file ----------------------------------------------------

# Type: Package
# Package: poset
# Title: Analysis of Partially Ordered Data
# Version: 1.0
# Author: Lu Mao
# Maintainer: Lu Mao <lmao@biostat.wisc.edu>
#   Description: Win ratio
# License: CC BY 4.0
# URL: https://sites.google.com/view/lmaowisc/
#   Depends:
#   R (>= 3.10)
# Suggests: knitr, rmarkdown
# VignetteBuilder:
#   knitr
# Config/testthat/edition: 3
# Encoding: UTF-8
# LazyData: true
# RoxygenNote: 7.3.1

# Commit changes to git ---------------------------------------------------

# Prerequisites:
# • GitHub account
# • create_github_token() - follow instructions
# • gitcreds::gitcreds_set() - paste PAT
# • git_sitrep() - verify
# use_github()
# - push content to new repository on GitHub

# usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)

edit_r_environ()



###### GitHub #######
create_github_token()
gitcreds::gitcreds_set()
git_sitrep()

use_github()
####################

# Documentation -----------------------------------------------------------

# RStudio: Code > Insert Roxygen Skeleton
## Ctrl+Alt+Shift+R
# Special comments (#') above function
#   definition in R/*.R

#### ----------- an example ---
#' Multiplicative win ratio (MWR) regression analysis
#'
#' @description Fit a multiplicative win ratio (MWR) regression model to
#' partially ordered outcome against covariates
#' @return An object of class \code{MWRreg} with the following components:
#' \item{beta}{A vector of estimated regression coefficients.}
#' @seealso \code{\link{wprod}}
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases MWRreg
#' @keywords MWRreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' Under revision.
#' @examples
#' set.seed(12345)
#' n <- 200
#' Z <- cbind(rnorm(n),rnorm(n))
#' \dontrun{
#'   use_git()
#' }
####

### Steps ###
# Go to function definition: Ctrl+.(type function name)
# • Cursor in function definition
# • Insert roxygen skeleton (Ctrl+Alt+Shift+R)
# • Complete the roxygen fields
# • document() (Ctrl+Shift+D) - create .rd files in ./man
# • ?myfunction

document()

# document() updates the NAMESPACE file with directives from Roxygen
# comments in your R code about export() and import()


# Package-level documentation ---------------------------------------------

use_package_doc()
#> ✔ Writing 'R/mypackage-package.R'
#> • Modify ‘R/mypackage-package.R’
document()
# ?poset

check() # again



# Install package to your library -----------------------------------------

install()

# build documentation pdf
# If you're in the package root
build_manual(pkg = ".", path = ".")

use_pkgdown_github_pages()
use_pkgdown_github_pages()
build_readme()

# Run once to configure your package to use pkgdown
# usethis::use_pkgdown()
# pkgdown::build_site()

# usethis::use_pkgdown_github_pages()

?icsurvfit
# Work space  ------------------------------------------------------------
library(IntCens)
## edit and test code ======================================
# Simple example
data(bcos)
n <- nrow(bcos) # sample size
L <- bcos$left + rnorm(n,0,0.0001) # add random noise
R <- bcos$right + rnorm(n,0,0.0001)
Z <- as.numeric(bcos$treatment=="RadChem")

# Nonparametric
obj_np <- icsurvfit(L, R)
# Cox model
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
# Propotional odds model
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



## BMA HIV study
## Bangkok Metropolitan Administration HIV Study
data <- read.table("bam.txt")


# install the IntCens package from local file
# install.packages("IntCens_0.1.0.tar.gz",
#                  repos = NULL,
#                  type = "source")
# library(IntCens)

## BMA data
head(data)

# get the response data for ICSurvICM()
delta <- data$delta
gamma <- data$gamma
n <- nrow(data)
U <- data$U
V <- data$V

L <- ifelse(delta == 1, 0, ifelse(gamma == 1, U, V))
R <- ifelse(delta == 1, U, ifelse(gamma == 1, V, Inf))

PH_fit <- icsurvfit(L, R, Z = data[,5:9], model = "PH")
PH_fit$beta
PH_fit$var
plot(PH_fit)

PO_fit <- icsurvfit(L, R, Z = data[,5:9], model = "PO")


print.ICSurv(obj.PH)

document()

## Read and document data

bcos <- read.table("bcos.txt")
bcos
# Save and document data bcos
# Now save it to the /data directory as an .rda file
usethis::use_data(bcos, overwrite = TRUE)



n <- nrow(bcos)

set.seed(123)
L <- bcos$left + rnorm(n,0,0.0001)
R <- bcos$right + rnorm(n,0,0.0001)

Z <- as.numeric(bcos$treatment=="RadChem")

obj_np <- icsurvfit(L, R, Z)

obj_ph <- icsurvfit(L, R, Z, model = "PH", maxiter = 100)
obj_po <- icsurvfit(L, R, Z, model = "PO", maxiter = 100)

#organize the data into desired format
delta=rep(0,n)
gamma=rep(0,n)
U=rep(0,n)
V=rep(0,n)
for (i in 1:n){
  if (bcos$left[i]==0){
    U[i]=bcos$right[i]
    delta[i]=1
    gamma[i]=0
  }
  else{

    if (bcos$right[i]==Inf){
      V[i]=bcos$left[i]
      delta[i]=0
      gamma[i]=0
    }
    else{
      U[i]=bcos$left[i]
      V[i]=bcos$right[i]
      delta[i]=0
      gamma[i]=1
    }
  }
}
U=U+rnorm(n,0,0.0001)
V=V+rnorm(n,0,0.0001)


obj.PH <- ICSurvICM(delta, gamma, U, V, Z, model = "PH")

obj.PH$beta


?ICSurvICM
