#' @title Simulated targeting estimation with sparse data.
#'
#' @description \code{sim_tau} is used to simulate estimated targeting
#' parameter with sparse linear data given N.
#'
#' @param nunber_sample number of sample in simulated data.
#' @param number_covariate number of covariates.
#' @param core number of core to use function.
#'
#' @return A dataframe containing follows:
#' \item{est}{ Out-of-sample MSE.}
#' \item{sample}{ Sample size in simulation data.}
#' \item{method}{ Prediction method.}
#'
#' @export
#'
#' @examples
#' sim <- sim_fix_n(number_sample = 100; number_covariate = 10; core = 1)
#'
#' @references No.
#


sim_tau <- function(nunber_sample,
                      number_covariate,
                      number_interation,
                      oracle = FALSE,
                      LASSO = FALSE,
                      RL_forest= FALSE,
                      RL = FALSE,
                      unpenalty_LASSO = FALSE,
                      wrong_simple = FALSE,
                      over_parameter = FALSE,
                    adaptive_LASSO = FALSE,
                    MCP_LASSO = FALSE,
                      core = 1){
  require(tidyverse)
  require(gamlr)
  require(furrr)
  require(glmnet)
  require(ranger)
  require(ncvreg)
  # Parameter ----
  N <- nunber_sample
  L <- number_covariate
  b <- number_interation
  # Individual estimater ----
    est_oracle <- function(X,D,Y){
    x <- X$var1_square
    lm(Y ~ 0 + D + x) |>
      coef() %>%
      .[1]
    }
  est_LASSO <- function(X,D,Y){
    x <- bind_cols(D = D,X) |> as.matrix()
    fit <- gamlr(x = x,
                 y = Y)
    coef(fit)[2]
  }
  est_adaptive_LASSO <- function(X,D,Y){
    x <- bind_cols(D = D,X) |> as.matrix()
    temp <- lm(Y ~ x) |> coef()
    weights <- 1/abs(temp[-1])
    cv <- cv.glmnet(x = x,
                    y = Y,
                    penalty.factor = 1/weights)
    fit <- glmnet(x = x,
                  y = Y,
                  penalty.factor = 1/weights,
                  lambda = cv$lambda.min)
    coef(fit)[2]
  }
  est_MCP_LASSO <- function(X,D,Y){
    x <- bind_cols(D = D,X) |> as.matrix()
    fit <- cv.ncvreg(X = x,
                    y = Y)
    coef(fit)[2]
  }
  est_RL_forest <- function(X,D,Y){
    x <- X |> as.matrix()
    fit <- ranger(y = Y,
                  x = X)
    Y.hat <- fit$predictions
    fit <- ranger(y = D,
                  x = X)
    D.hat <- fit$predictions
    Y.oht <- Y - Y.hat
    D.oht <- D - D.hat
    fit <- lm(Y.oht ~ 0 + D.oht)
    coef(fit)
  }
  est_RL <- function(X,D,Y){
    x <- X |> as.matrix()
    fit <- doubleML(x = x,
                    d = D,
                    y = Y)
    coef(fit)
  }
  est_unpenalty_LASSO <- function(X,D,Y){
    x <- bind_cols(D = D,X) |> as.matrix()
    fit <- gamlr(x = x,
                 y = Y,
                 free = "D")
    coef(fit)[2]
  }
  est_wrong_simple <- function(X,D,Y){
    x <- X$var1
    lm(Y ~ 0 + D + x) |>
      coef() %>%
      .[1]
  }
  est_over_parameter <- function(X,D,Y){
    x <- X |> as.matrix()
    lm(Y ~ D + x) |>
      coef() %>%
      .[2]
  }
  # Aggregation ----
  repeat_est <- function(i){
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
    Y <- df$Y
    D <- df$D
    X <- select(df,-Y,-D)
    tibble(est = c(if (oracle) {est_oracle(X,D,Y)} else {NA},
                   if (over_parameter) {est_over_parameter(X,D,Y)} else{NA},
                   if (wrong_simple) {est_wrong_simple(X,D,Y)} else{NA},
                   if (LASSO) {est_LASSO(X,D,Y)} else{NA},
                   if (unpenalty_LASSO) {est_unpenalty_LASSO(X,D,Y)} else{NA},
                   if (RL) {est_RL(X,D,Y)} else{NA},
                   if (RL_forest) {est_RL_forest(X,D,Y)} else{NA},
                   if (adaptive_LASSO) {est_adaptive_LASSO(X,D,Y)} else{NA},
                   if (MCP_LASSO) {est_MCP_LASSO(X,D,Y)} else{NA}
                   ),
           sample = c(N,N,N,N,N,N,N,N,N),
           method = c("correct parametric",
                      "over-parameteric",
                      "almost-correct parametric",
                      "LASSO",
                      "LASSO without penalty",
                      "R_learner",
                      "R_learner with Forest",
                      "Adaptive LASSO",
                      "MCP LASSO")
           )
  }
  # Repeat----
  plan(multisession, workers = core)
  future_map_dfr(1:b,
                 repeat_est,
                 .progress = TRUE,
                 .options = furrr_options(seed = NULL))
}
