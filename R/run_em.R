#' Run EM iterations for the joint model
#'
#' Internal wrapper running the EM algorithm until convergence.
#' This version supports:
#' - Gaussian longitudinal: theta has sigma2
#' - Poisson longitudinal: theta has no sigma2
#' - Random effects covariance stored as either `D2` or `Sigma.b`
#'
#' @param theta Initial parameter list.
#' @param em_step Function(theta, nMC) -> list(theta=theta.new).
#' @param nMC Initial Monte Carlo sample size.
#' @param tol Convergence tolerance (absolute/relative combined).
#' @param iter_max Max iterations.
#' @param has_sigma2 Logical; TRUE for Gaussian, FALSE for Poisson.
#' @param verbose Logical; print iteration info if TRUE.
#'
#' @return list(theta=final_theta, iter=iterations_used)
#' @keywords internal
#' @noRd
run_em <- function(theta,
                   em_step,
                   nMC,
                   tol,
                   iter_max,
                   has_sigma2 = TRUE,
                   verbose = FALSE) {
  
  tolerance <- Inf
  iter <- 0
  
  ## helper: extract covariance field name
  get_cov <- function(th) {
    if (!is.null(th$Sigma.b)) return(th$Sigma.b)
    if (!is.null(th$D2)) return(th$D2)
    stop("theta must contain either Sigma.b or D2.")
  }
  set_cov <- function(th, val) {
    if (!is.null(th$Sigma.b)) th$Sigma.b <- val
    if (!is.null(th$D2)) th$D2 <- val
    th
  }
  
  while (tolerance >= tol && iter < iter_max) {
    
    ## old
    alpha.old   <- theta$alpha
    cov.old     <- get_cov(theta)
    gamma.1.old <- theta$gamma.1
    gamma.2.old <- theta$gamma.2
    haz.1.old   <- theta$haz.1
    haz.2.old   <- theta$haz.2
    rho.old     <- theta$rho
    if (has_sigma2) sigma2.old <- theta$sigma2
    
    ## one EM step
    out <- em_step(theta, nMC)
    theta.new <- out$theta
    
    ## new
    alpha.new   <- theta.new$alpha
    cov.new     <- get_cov(theta.new)
    gamma.1.new <- theta.new$gamma.1
    gamma.2.new <- theta.new$gamma.2
    haz.1.new   <- theta.new$haz.1
    haz.2.new   <- theta.new$haz.2
    rho.new     <- theta.new$rho
    if (has_sigma2) sigma2.new <- theta.new$sigma2
    
    ## flatten for convergence check
    if (has_sigma2) {
      par_old <- unlist(list(alpha.old, sigma2.old, as.vector(cov.old), gamma.1.old, gamma.2.old, rho.old))
      par_new <- unlist(list(alpha.new, sigma2.new, as.vector(cov.new), gamma.1.new, gamma.2.new, rho.new))
    } else {
      par_old <- unlist(list(alpha.old, as.vector(cov.old), gamma.1.old, gamma.2.old, rho.old))
      par_new <- unlist(list(alpha.new, as.vector(cov.new), gamma.1.new, gamma.2.new, rho.new))
    }
    
    abssub <- abs(par_new - par_old)
    rbssub <- abssub / pmax(abs(par_old), .Machine$double.eps)
    
    tolerance <- min(max(abssub), max(rbssub))
    
    iter <- iter + 1
    theta <- theta.new
    
    if (isTRUE(verbose)) {
      message(sprintf("[EM] iter=%d, nMC=%d, tol=%.3e", iter, nMC, tolerance))
    }
    
    ## increase MC size gradually (cap)
    nMC <- min(nMC + floor(nMC / 3), 2000)
  }
  
  if (tolerance < tol && isTRUE(verbose)) message("EM converged.")
  list(theta = theta, iter = iter)
}