#' @title Simulated data.
#'
#' @description \code{sim_data} is used to generate simulate estimated with sparse linear. Assume Y = -5D + 25var1_square + u and D = -5var1_square.
#'
#' @param nunber_sample number of sample in simulated data.
#' @param number_covariate number of covariates.
#'
#' @return A dataframe containing follows:
#' \item{vari}{ ith control}
#' \item{vari_square}{ squated ith control}
#' \item{Y}{ Outcome correlated with var1}
#' \item{D}{ Treatment correlated with var2}
#'
#' @export
#'
#' @examples
#' sim <- sim_data(nunber_sample = 1000, number_covariate = 5)
#'
#' @references No.
#


sim_data <- function(nunber_sample,
                     number_covariate){
  require(tidyverse)
  # Parameter ----
  N <- nunber_sample
  L <- number_covariate
  # data
  X <-
    matrix(runif(N*L,-5,5), N, L) |>
    as.data.frame()
  colnames(X) <- str_c("var",1:L)
  X2 <- X^2
  colnames(X2) <- str_c(colnames(X),"_square")
  df <-
    X |>
    bind_cols(X2) |>
    mutate(D = -5*var1_square + rnorm(N,0,10),
           Y = 5*D + 25*var1_square + rnorm(N,0,100)
           )
  df
    }
