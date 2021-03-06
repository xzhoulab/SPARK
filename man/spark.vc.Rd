% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcestimating.R
\name{spark.vc}
\alias{spark.vc}
\title{Fitting the count-based spatial model to estimate the parameters}
\usage{
spark.vc(
  object,
  covariates = NULL,
  lib_size = NULL,
  fit.maxiter = 500,
  fit.tol = 1e-05,
  fit.model = "poisson",
  num_core = 1,
  verbose = FALSE
)
}
\arguments{
\item{object}{SPARK object}

\item{covariates}{The covariates in experiments, i.e. confounding factors/batch effect}

\item{lib_size}{The read depth for each cell}

\item{fit.maxiter}{Iteration}

\item{fit.tol}{Tolerance}

\item{fit.model}{The model to be either "poisson" or "gaussian"}

\item{num_core}{The number of core used when fitting the model}

\item{verbose}{Output fitting information}
}
\description{
Fitting the count-based spatial model to estimate the parameters
}
