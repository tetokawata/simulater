% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_pred_trace_n.R
\name{sim_pred_trace_n}
\alias{sim_pred_trace_n}
\title{Simulated prediction quality with sparse data}
\usage{
sim_pred_trace_n(
  seq,
  number_covariate,
  core = 1,
  oracle = FALSE,
  LASSO = FALSE,
  wrong = FALSE,
  over = FALSE
)
}
\arguments{
\item{seq}{Integer sequance from minimum to maximum sample size.}

\item{number_covariate}{Number of covariates.}

\item{core}{Number of core to use function.}

\item{oracle}{Logic Include oracle OLS or not}

\item{LASSO}{Logic Include LASSO or not}

\item{wrong}{Logic Include Wrong OLS or not}

\item{over}{Logic Include over-parametric OLS or not}
}
\value{
A dataframe containing follows:
\item{est}{ Out-of-sample MSE.}
\item{sample}{ Sample size in simulation data.}
\item{method}{ Prediction method.}
}
\description{
\code{sim_pred_trace_n} is used to \code{sim_pred_trace_n} is
used to simulate prediction accuracy (out-of-sample MSE) with sparse linear data..
}
\examples{
sim <- sim_pred_trace_n(seq = c(100:150); number_covariate = 10; core = 1)

}
\references{
No.
}
