#' Summarize a fitted joint model
#'
#' @param object An object of class `mjoint_fit`.
#' @param ... Unused.
#'
#' @return A list invisibly.
#' @export
summary.mjoint_fit <- function(object, ...) {
  stopifnot(inherits(object, "mjoint_fit"))
  
  th <- object$theta
  
  cat("=================================================\n")
  cat("Joint model fit (Clayton copula)\n")
  cat("Longitudinal type:", object$y_type, "\n")
  cat("EM iterations:", object$iter, "\n")
  cat("-------------------------------------------------\n")
  
  cat("alpha:\n")
  print(th$alpha)
  
  if (!is.null(th$sigma2)) {
    cat("\nsigma2:\n")
    print(th$sigma2)
  }
  
  if (!is.null(th$D2)) {
    cat("\nRandom-effects covariance (D2):\n")
    print(th$D2)
  }
  if (!is.null(th$Sigma.b)) {
    cat("\nRandom-effects covariance (Sigma.b):\n")
    print(th$Sigma.b)
  }
  
  cat("\nSurvival 1 gamma:\n")
  print(th$gamma.1)
  cat("\nSurvival 2 gamma:\n")
  print(th$gamma.2)
  
  cat("\nClayton rho:\n")
  print(th$rho)
  
  cat("=================================================\n")
  
  invisible(object)
}