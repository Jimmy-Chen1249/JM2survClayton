#' JM2survClayton: Joint Modeling of a Longitudinal Outcome and Two Correlated Survival Endpoints via Clayton Copula
#'
#' Implements EM-based joint models linking a longitudinal outcome with two correlated
#' time-to-event endpoints through shared random effects. Dependence between the two
#' survival endpoints is modeled using a Clayton survival copula.
#'
#' The package supports:
#' \itemize{
#'   \item Gaussian and Poisson longitudinal submodels;
#'   \item One- or two-dimensional random effects;
#'   \item User-supplied time-varying covariate functions for the survival association.
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib JM2survClayton, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL