

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



use_pkgdown_github_pages()
use_pkgdown_github_pages()
build_readme()

# Run once to configure your package to use pkgdown
usethis::use_pkgdown()
pkgdown::build_site()

usethis::use_pkgdown_github_pages()


# Work space  ------------------------------------------------------------
library(WR)
## edit and test code ======================================
### test wrtest()

data("non_ischemic")
df <- non_ischemic
id <- df$ID
time <- df$time
status <- df$status
Z <- df[, 4:ncol(df)]


# get the number of rows and number of covariates.
nr <- nrow(non_ischemic)
p <- ncol(non_ischemic)-3

# extract ID, time, status and covariates matrix Z from the data.
# note that: ID, time and status should be column vector.
# covariatesZ should be (nr, p) matrix.
id <- non_ischemic[,"ID"]
time <- non_ischemic[,"time"]
status <- non_ischemic[,"status"]
Z <- as.matrix(non_ischemic[, 4:(3+p)],nr,p)

system.time({
obj <- pwreg(id, time, status, Z, strata = NULL, fixedL = TRUE, eps = 1e-4,
       maxiter = 50)
})

library(dplyr)
library(tidyr)
library(stringr)



data("non_ischemic")
df <- non_ischemic
id <- df$ID
time <- df$time
status <- df$status
Z <- df[, 4:ncol(df)]


# library(microbenchmark)
system.time({
obj <- pwreg1(id, time, status, Z, eps = 1e-8)
})
# obj$beta
# user  system elapsed
# 0.44    0.09    0.54
obj

z1 <- non_ischemic[1, 4:ncol(non_ischemic)]
z2 <- z1
z1$trt_ab = 1
preds <- predict(obj, z1 = z1, z2 = z2, contrast = TRUE)

preds
# # A tibble: 357 × 5
# time win_prob win_prob_se loss_prob loss_prob_se
# <dbl>    <dbl>       <dbl>     <dbl>        <dbl>
# 1     0  0.00629     0.00319   0.00763      0.00387
# 2     2  0.00755     0.00355   0.00916      0.00428
# 3     5  0.00881     0.00396   0.0107       0.00481
# 4     7  0.0101      0.00432   0.0122       0.00528
# 5     8  0.0126      0.00493   0.0153       0.00606
# 6     9  0.0138      0.00532   0.0168       0.00657
# ?predict.pwreg1




str(x)

x <- obj



insert_comment <- function() {
  rstudioapi::insertText(text = "#> ")
}

x$pw_win |> count(component)


devtools::install()


?residuals
?pwreg1


obj$wfun
# obj$resids_n
obj_resids <- residuals(obj)

install.packages("reprex")
y <- 1:3
y

library(ggplot2)

obj_resids |>
  ggplot(aes(x = factor(id), y = cook_d)) +
  geom_segment(aes(x = factor(id), xend = factor(id), y = 0, yend = cook_d)) +
  geom_point() +
  geom_text(data = obj_resids |> filter(cook_d > 0.02),
            aes(label = id), vjust = - 1) +
  theme_minimal()

# transform
arctanh <- function(x) {
  0.5 * log((1 + x) / (1 - x))
}




sum(obj_resids$M)

obj_resids |>
  ggplot(aes(x = factor(id), y = M)) +
  geom_segment(aes(x = factor(id), xend = factor(id), y = 0, yend = M)) +
  geom_point() +
  geom_text(data = obj_resids |> filter(abs(M) > 0.45), aes(label = id, vjust = - sign(M))) +
  theme_minimal()

obj_resids |>
  ggplot(aes(x = factor(id), y = hii)) +
  geom_segment(aes(x = factor(id), xend = factor(id), y = 0, yend =  hii)) +
  geom_point() +
  geom_text(data = obj_resids |> filter(hii > 0.1), aes(label = id, vjust = - 1)) +
  theme_minimal()


summary(non_ischemic$bipllvef)

beta <- obj$beta
se <- sqrt(diag(obj$Var))
obj$iter
Zn <- obj$Zn

p <- length(beta)
n <- nrow(Zn)

cbind(beta, se)
# beta          se
# trt_ab          0.19280012 0.124426769
# age            -0.01281768 0.005637639
# sex            -0.15887407 0.127243976
# Black.vs.White -0.30261101 0.143782446
# Other.vs.White -0.35805309 0.337188274
# bmi            -0.01792220 0.009612621
# bipllvef        0.02153551 0.008511336
# hyperten       -0.02757744 0.143281955
# COPD           -0.40365707 0.203150789
# diabetes        0.06664939 0.139771040
# acei           -0.10443828 0.154637988
# betab          -0.53240924 0.324031930
# smokecurr      -0.05821292 0.165329855

# Old results from pwreg()
# Total number of pairs: 101475
# Wins-losses on death:  7644 (7.5%)
# Wins-losses on non-fatal event:  78387 (77.2%)
# Indeterminate pairs 15444 (15.2%)
#
# Newton-Raphson algorithm converged in 5 iterations.
#
# Overall test: chisq test with 13 degrees of freedom;
# Wald statistic 24.9 with p-value 0.02392931
#
# Estimates for Regression parameters:
#   Estimate         se z.value p.value
# trt_ab          0.1906687  0.1264658  1.5077 0.13164
# age            -0.0128306  0.0057285 -2.2398 0.02510 *
#   sex            -0.1552923  0.1294198 -1.1999 0.23017
# Black.vs.White -0.3026335  0.1461330 -2.0709 0.03836 *
#   Other.vs.White -0.3565390  0.3424360 -1.0412 0.29779
# bmi            -0.0181310  0.0097582 -1.8580 0.06316 .
# bipllvef        0.0214905  0.0086449  2.4859 0.01292 *
#   hyperten       -0.0318291  0.1456217 -0.2186 0.82698
# COPD           -0.4023069  0.2066821 -1.9465 0.05159 .
# diabetes        0.0703990  0.1419998  0.4958 0.62006
# acei           -0.1068201  0.1571317 -0.6798 0.49662
# betab          -0.5344979  0.3289319 -1.6250 0.10417
# smokecurr      -0.0602350  0.1682826 -0.3579 0.72039
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# install.package("remotes")   #In case you have not installed it.
remotes::install_github("GuangchuangYu/hexSticker")




# Create hexsticker -------------------------------------------------------
library(tidyverse)
# Example data
data <- tibble(
  y = rep(1:5, each = 2),
  x = c(
    -5, 5,
    -4, 4,
    -3, 3,
    -2, 2,
    -1, 1
  )

)

# levels(data$Covariate)


# Create the forest plot
ggplot(data, aes(x = y, y = x)) +
  geom_bar(stat = "identity", fill = "grey", width = 0.925) +
  coord_flip() +
  # plot a transparent line segment at x = 0
  geom_segment(aes(y = 0, yend = 0, x = 0.525, xend = 5.475), color = "#2C3E50", linewidth = 0.5, linetype = 1) +
  theme_void()


ggsave("fp.png", width = 2, height = 1, dpi = 300)


# library(hexSticker)
sticker(
   "fp.png", # Base image file path
  package = "WR",          # Package name to display on the sticker
  p_size = 20,             # Font size of the package name
  p_color = "#ffffff",     # Color of the package name
  s_x = 1, s_y = 0.8,      # Position of the base image
  s_width = 0.7,           # Width of the base image
  h_fill = "#2C3E50",      # Background color of the sticker
  h_color = "#2C3E50",     # Border color of the sticker
  filename = "logo.png", # Output file name
  url = "https://lmaowisc.github.io/WR/", # URL to link to
  # adjust url size
  u_size = 4,             # Font size of the URL
)




# devtools::install_github("lmaowisc/WR")
library(WR)
library(tidyverse) # For data wrangling and visualization
# Generate data in long format
n <- 10
df_y <- tibble(
  id = c(1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 9, 9, 10),
  time = c(0.2, 2, 2.3, 1, 4.5, 3.6, 2.2, 3.2, 5,
           2.1, 3.8, 1.2, 4.5, 2.3),
  # status = 1 for death, 2 for hospitalization, 0 for censoring
  status = c(2, 1, 0, 2, 0, 1, 2, 1, 1, 0, 1, 2, 1, 0)
)
# Covariates
set.seed(123)
Z <- tibble(
  id = 1:n,
  Z1 = rnorm(n),
  Z2 = rbinom(n, 1, 0.5)
)
# Merge by id
df <- left_join(df_y, Z, by = "id")
# Fit the model
obj <- pwreg1(df$id, df$time, df$status, df[, 4:ncol(df)])
# Calculate residuals
resids <- residuals(obj)
resids
#># A tibble: 10 × 10
#> id delta R Lambda M r hii cook_d Z1 Z2
#><dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 0 1 0.324 -0.324 -0.324 0.179 0.143 0.161 -0.037
#> 2 0.444 0.444 0.161 0.284 0.639 0.107 0.056 -0.116 0.088
#> 3 0.444 1 0.651 -0.206 -0.206 0.206 0.072 -0.244 -0.096
#> 4 0.222 0.667 0.286 -0.063 -0.095 0.066 0.002 -0.006 -0.027
#> 5 0.111 0.889 0.378 -0.267 -0.301 0.061 0.025 -0.05 -0.052
#> 6 0.667 0.667 0.462 0.205 0.307 0.254 0.099 0.333 0.058
#> 7 0.333 0.333 0.142 0.191 0.574 0.047 0.01 0.048 0.079
#> 8 0.333 0.667 0.152 0.181 0.272 0.426 0.219 -0.167 -0.013
#> 9 0.556 1 0.695 -0.139 -0.139 0.337 0.076 0.133 0.088
#> 10 0.444 0.444 0.305 0.139 0.313 0.316 0.067 -0.092 -0.088


# Install and load WR
devtools::install_github("lmaowisc/WR")
library(WR)
library(tidyverse) # For data wrangling and visualization
# Generate data in long format
n <- 10
df_y <- tibble(
  id = c(1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 8, 9, 9, 10),
  time = c(0.2, 2, 2.3, 1, 4.5, 3.6, 2.2, 3.2, 5, 2.1, 3.8, 1.2, 4.5, 2.3),
  # status = 1 for death, 2 for hospitalization, 0 for censoring
  status = c(2, 1, 0, 2, 0, 1, 2, 1, 1, 0, 1, 2, 1, 0)
)

# Covariates
set.seed(123)
Z <- tibble(
  id = 1:n,
  Z1 = rnorm(n),
  Z2 = rbinom(n, 1, 0.5)
)
# Merge by id
df <- left_join(df_y, Z, by = "id")
# Fit the model
obj <- pwreg1(df$id, df$time, df$status, df[, 4:ncol(df)])
# Calculate predictions
# Specify covariate vectors for comparison
z1 <- c(1, 0)
z2 <- c(0, 0)
# Predict win-loss probabilities, along with
# win ratio, win odds, and net benefit
preds <- predict(obj, z1 = z1, z2 = z2, contrast = TRUE)
preds[1:9] # Basic output
preds[c(1, 10:19)] # Additional output by contrast = TRUE
