#' @title Simulated prediction quality with sparse data
#'
#' @description \code{sim_pred_trace_n} is used to \code{sim_pred_trace_n} is
#' used to simulate prediction accuracy (out-of-sample MSE) with sparse linear data..
#'
#' @param seq Integer sequance from minimum to maximum sample size.
#' @param number_covariate Number of covariates.
#' @param core Number of core to use function.
#' @param oracle Logic Include oracle OLS or not
#' @param LASSO Logic Include LASSO or not
#' @param wrong Logic Include Wrong OLS or not
#' @param over Logic Include over-parametric OLS or not
#'
#' @return A dataframe containing follows:
#' \item{est}{ Out-of-sample MSE.}
#' \item{sample}{ Sample size in simulation data.}
#' \item{method}{ Prediction method.}
#'
#' @export
#'
#' @examples
#' sim <- sim_pred_trace_n(seq = c(100:150); number_covariate = 10; core = 1)
#'
#' @references No.
#

sim_pred_trace_n <- function(seq,
                             number_covariate,
                             core = 1,
                             oracle = FALSE,
                             LASSO = FALSE,
                             wrong = FALSE,
                             over = FALSE){
  require(tidyverse)
  require(gamlr)
  require(furrr)
  # Parameter ----
  L <- number_covariate
  # Individual estimater ----
  pred_oracle <- function(X,D,Y,X_new,D_new,Y_new,oracle){
    result <- NA
    if (oracle == TRUE){train <- bind_cols(D = D,X)
    test <- bind_cols(D = D_new, X_new)
    Y.hat <-
      lm(Y ~ 0 + D + var1_square,
       data = train) |>
      predict(test)
    result <- mean((5*D_new + 25*X_new$var1_square - Y.hat)^2)}
    result
  }
  pred_LASSO <- function(X,D,Y,X_new,D_new,Y_new,LASSO){
    result <- NA
    if (LASSO == TRUE){x <- bind_cols(D = D,X) |> as.matrix()
    x.test <- bind_cols(D = D_new,X_new) |> as.matrix()
    Y.hat <-
      gamlr(x = x,
          y = Y) |>
      predict(x.test)
    result <- mean((5*D_new + 25*X_new$var1_square - Y.hat)^2)}
    result
  }
  pred_wrong_simple <- function(X,D,Y,X_new,D_new,Y_new,wrong){
    result <- NA
    if (wrong == TRUE){train <- bind_cols(D = D,X)
    test <- bind_cols(D = D_new, X_new)
    Y.hat <-
      lm(Y ~ 0 + D + var1,
       data = train) |>
      predict(test)
    result <- mean((5*D_new + 25*X_new$var1_square - Y.hat)^2)}
    result
  }
  pred_over_parameter <- function(X,D,Y,X_new,D_new,Y_new,over){
    result <- NA
    if (over == TRUE){train <- bind_cols(D = D,X)
    test <- bind_cols(D = D_new, X_new)
    Y.hat <-
      lm(Y ~ .,
       data = train) |>
      predict(test)
    result <- mean((5*D_new + 25*X_new$var1_square - Y.hat)^2)}
    result
  }
  # Aggregation ----
  repeat_est <- function(N){
    oracle = oracle
    LASSO = LASSO
    wrong = wrong
    over = over
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
             Y = 5*D + 25*var1_square + rnorm(N,0,100))
    X <-
      matrix(runif(1000*L,-5,5), 1000, L) |>
      as.data.frame()
    colnames(X) <- str_c("var",1:L)
    X2 <- X^2
    colnames(X2) <- str_c(colnames(X),"_square")
    df_new <-
      X |>
      bind_cols(X2) |>
      mutate(D = -5*var1_square + rnorm(1000,0,10),
             Y = 5*D + 25*var1_square + rnorm(1000,0,100))
    Y <- df$Y
    D <- df$D
    X <- select(df,-Y,-D)
    Y_new <- df_new$Y
    D_new <- df_new$D
    X_new <- select(df_new,-Y,-D)
    tibble(est = c(pred_oracle(X,D,Y,X_new,D_new,Y_new,oracle),
                   pred_over_parameter(X,D,Y,X_new,D_new,Y_new,over),
                   pred_wrong_simple(X,D,Y,X_new,D_new,Y_new,wrong),
                   pred_LASSO(X,D,Y,X_new,D_new,Y_new,LASSO)),
           sample = c(N,N,N,N),
           method = c("correct parametric",
                      "over-parameteric",
                      "almost-correct parametric",
                      "LASSO"))
  }
  # Repeat----
  plan(multisession, workers = core)
  future_map_dfr(seq,
                 repeat_est,
                 .progress = TRUE,
                 .options = furrr_options(seed = NULL))
  }
