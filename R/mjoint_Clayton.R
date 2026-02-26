#' Fit a joint model with two correlated survival endpoints via Clayton copula
#'
#' Fits a joint model that links a longitudinal outcome with **two correlated**
#' time-to-event endpoints via **shared random effects**. Dependence between the
#' two event times is modeled by a **Clayton survival copula**. Estimation is
#' performed via a Monte Carlo EM algorithm (MCEM).
#'
#' The longitudinal submodel can be:
#' \itemize{
#'   \item a Gaussian linear mixed model (LMM), or
#'   \item a Poisson generalized linear mixed model (GLMM).
#' }
#'
#' The two survival submodels are Cox-type models and can include
#' user-specified time-varying covariates via \code{timevars1} and \code{timevars2}.
#'
#' @param Mixeffect Mixed-effects model formula for the longitudinal process,
#'   e.g., \code{y ~ time + X1 + X2 + (1 | id)} or \code{y ~ time + X1 + X2 + (1 + time | id)}.
#' @param Longdat Longitudinal \code{data.frame} containing at least \code{id}, \code{time},
#'   the response \code{y}, and covariates used in \code{Mixeffect}.
#' @param Surdat1 Survival \code{data.frame} for endpoint 1. Must contain \code{id},
#'   event time \code{surtime}, censoring indicator \code{cens}, baseline covariates in \code{formSurv1},
#'   and any columns referenced by \code{timevars1}.
#' @param Surdat2 Survival \code{data.frame} for endpoint 2 (same format as \code{Surdat1}).
#' @param formSurv1 Cox formula for endpoint 1, e.g. \code{Surv(surtime, cens) ~ W}.
#' @param formSurv2 Cox formula for endpoint 2.
#' @param timevars1 Named list \code{list(fixed = FUN, rand = FUN)} for endpoint 1.
#'   \code{fixed(t, index, Sur = ..., ...)} returns a matrix of time-varying covariates
#'   entering the fixed part of the survival model (dimension \code{p1}).
#'   \code{rand(t, index, ...)} returns the random-effect design used in the association
#'   term (dimension matches the random-effect dimension in \code{Mixeffect}).
#' @param timevars2 Named list \code{list(fixed = FUN, rand = FUN)} for endpoint 2.
#' @param y_type Character string specifying the longitudinal distribution:
#'   \code{"gaussian"} or \code{"poisson"}.
#' @param rho0 Initial value for the Clayton copula parameter \eqn{\rho > 0}.
#'   Default is \code{0.5}.
#' @param nMC Initial Monte Carlo sample size used in the E-step.
#'   The algorithm may adaptively increase it.
#' @param iter_max Maximum number of EM iterations.
#' @param tol Convergence tolerance for the EM algorithm.
#' @param rho.region Search interval for \eqn{\rho} used by \code{\link[stats]{optimize}}.
#' @param eps_prob Numerical guard used in likelihood evaluations to avoid \code{log(0)}
#'   or division by zero.
#' @param verbose Logical; if \code{TRUE}, print progress information.
#'
#' @return An object of class \code{"mjoint_fit"} containing parameter estimates and
#'   model components. Use \code{\link{summary.mjoint_fit}} to summarize.
#'
#' @examples
#' \dontrun{
#' ## ------------------------------------------------------------
#' ## Example datasets shipped with the package
#' ## ------------------------------------------------------------
#' ## Each dataset is a list called `blocks` with:
#' ##   blocks$Longdat, blocks$Surdat1, blocks$Surdat2
#' ##
#' ## Available datasets:
#' ##   Gaussian_Clayton_One
#' ##   Poisson_Clayton_One
#' ##   Gaussian_Clayton_Two
#' ##   Poisson_Clayton_Two
#'
#' ## ------------------------------------------------------------
#' ## Helper: build time-varying covariates (fixed) and Z_sur (rand)
#' ## - fixed uses B1,B2,B3 columns in Surdat
#' ## - rand depends on random-effect dimension:
#' ##     re_dim=1 -> 1
#' ##     re_dim=2 -> (1, t)
#' ## ------------------------------------------------------------
#' make_timevars <- function(re_dim = 1) {
#'   stopifnot(re_dim %in% c(1, 2))
#'
#'   fixed_fun <- function(t, index, Sur = NULL, ...) {
#'     if (is.null(Sur)) stop("fixed_fun: Sur is NULL (build_surv_block must pass Sur=).")
#'     B1 <- Sur$B1[index]; B2 <- Sur$B2[index]; B3 <- Sur$B3[index]
#'     z  <- B1 * (t <= B3) + B2 * (t > B3)
#'     matrix(z, ncol = 1)
#'   }
#'
#'   rand_fun <- function(t, index, ...) {
#'     if (re_dim == 1) matrix(1, nrow = length(t), ncol = 1) else cbind(1, t)
#'   }
#'
#'   list(fixed = fixed_fun, rand = rand_fun)
#' }
#'
#' formSurv1 <- survival::Surv(surtime, cens) ~ W
#' formSurv2 <- survival::Surv(surtime, cens) ~ W
#'
#' ## ------------------------------------------------------------
#' ## Case 1g: Gaussian longitudinal + 1D random effects
#' ## ------------------------------------------------------------
#' data(Gaussian_Clayton_One)
#' Longdat  <- Gaussian_Clayton_One$Longdat
#' Surdat1  <- Gaussian_Clayton_One$Surdat1
#' Surdat2  <- Gaussian_Clayton_One$Surdat2
#'
#' Mixeffect <- y ~ time + X1 + X2 + (1 | id)
#' timevars1 <- make_timevars(re_dim = 1)
#' timevars2 <- make_timevars(re_dim = 1)
#'
#' set.seed(1)
#' fit_1g <- mjoint_Clayton(
#'   Mixeffect = Mixeffect,
#'   Longdat   = Longdat,
#'   Surdat1   = Surdat1,
#'   Surdat2   = Surdat2,
#'   formSurv1 = formSurv1,
#'   formSurv2 = formSurv2,
#'   timevars1 = timevars1,
#'   timevars2 = timevars2,
#'   y_type    = "gaussian",
#'   rho0      = 0.5,
#'   nMC       = 100,
#'   iter_max  = 200,
#'   tol       = 5e-3,
#'   rho.region = c(1e-4, 30),
#'   verbose   = TRUE
#' )
#' summary(fit_1g)
#'
#' ## ------------------------------------------------------------
#' ## Case 1p: Poisson longitudinal + 1D random effects
#' ## ------------------------------------------------------------
#' data(Poisson_Clayton_One)
#' Longdat  <- Poisson_Clayton_One$Longdat
#' Surdat1  <- Poisson_Clayton_One$Surdat1
#' Surdat2  <- Poisson_Clayton_One$Surdat2
#'
#' Mixeffect <- y ~ time + X1 + X2 + (1 | id)
#' timevars1 <- make_timevars(re_dim = 1)
#' timevars2 <- make_timevars(re_dim = 1)
#'
#' set.seed(1)
#' fit_1p <- mjoint_Clayton(
#'   Mixeffect = Mixeffect,
#'   Longdat   = Longdat,
#'   Surdat1   = Surdat1,
#'   Surdat2   = Surdat2,
#'   formSurv1 = formSurv1,
#'   formSurv2 = formSurv2,
#'   timevars1 = timevars1,
#'   timevars2 = timevars2,
#'   y_type    = "poisson",
#'   rho0      = 0.5,
#'   nMC       = 100,
#'   iter_max  = 200,
#'   tol       = 5e-3,
#'   rho.region = c(1e-4, 30),
#'   verbose   = TRUE
#' )
#' summary(fit_1p)
#'
#' ## ------------------------------------------------------------
#' ## Case 2g: Gaussian longitudinal + 2D random effects
#' ## ------------------------------------------------------------
#' data(Gaussian_Clayton_Two)
#' Longdat  <- Gaussian_Clayton_Two$Longdat
#' Surdat1  <- Gaussian_Clayton_Two$Surdat1
#' Surdat2  <- Gaussian_Clayton_Two$Surdat2
#'
#' Mixeffect <- y ~ time + X1 + X2 + (1 + time | id)
#' timevars1 <- make_timevars(re_dim = 2)
#' timevars2 <- make_timevars(re_dim = 2)
#'
#' set.seed(1)
#' fit_2g <- mjoint_Clayton(
#'   Mixeffect = Mixeffect,
#'   Longdat   = Longdat,
#'   Surdat1   = Surdat1,
#'   Surdat2   = Surdat2,
#'   formSurv1 = formSurv1,
#'   formSurv2 = formSurv2,
#'   timevars1 = timevars1,
#'   timevars2 = timevars2,
#'   y_type    = "gaussian",
#'   rho0      = 0.5,
#'   nMC       = 100,
#'   iter_max  = 200,
#'   tol       = 5e-3,
#'   rho.region = c(1e-4, 30),
#'   verbose   = TRUE
#' )
#' summary(fit_2g)
#'
#' ## ------------------------------------------------------------
#' ## Case 2p: Poisson longitudinal + 2D random effects
#' ## ------------------------------------------------------------
#' data(Poisson_Clayton_Two)
#' Longdat  <- Poisson_Clayton_Two$Longdat
#' Surdat1  <- Poisson_Clayton_Two$Surdat1
#' Surdat2  <- Poisson_Clayton_Two$Surdat2
#'
#' Mixeffect <- y ~ time + X1 + X2 + (1 + time | id)
#' timevars1 <- make_timevars(re_dim = 2)
#' timevars2 <- make_timevars(re_dim = 2)
#'
#' set.seed(1)
#' fit_2p <- mjoint_Clayton(
#'   Mixeffect = Mixeffect,
#'   Longdat   = Longdat,
#'   Surdat1   = Surdat1,
#'   Surdat2   = Surdat2,
#'   formSurv1 = formSurv1,
#'   formSurv2 = formSurv2,
#'   timevars1 = timevars1,
#'   timevars2 = timevars2,
#'   y_type    = "poisson",
#'   rho0      = 0.5,
#'   nMC       = 100,
#'   iter_max  = 200,
#'   tol       = 5e-3,
#'   rho.region = c(1e-4, 30),
#'   verbose   = TRUE
#' )
#' summary(fit_2p)
#' }
#'
#' @export
mjoint_Clayton <- function(Mixeffect,
                           Longdat, Surdat1, Surdat2,
                           formSurv1, formSurv2,
                           timevars1, timevars2,
                           y_type = c("gaussian", "poisson"),
                           rho0 = 0.5,
                           nMC = 100,
                           iter_max = 500,
                           tol = 5e-3,
                           rho.region = c(1e-4, 30),
                           eps_prob = 1e-100,
                           verbose = FALSE) {

  y_type <- match.arg(as.character(y_type), c("gaussian", "poisson"))

  ## =========================
  ## 1) longitudinal initialization
  ## =========================
  if (identical(y_type, "gaussian")) {
    lfit <- lme4::lmer(formula = Mixeffect, data = Longdat)
    alpha0 <- lme4::fixef(lfit)
    sigma20 <- stats::sigma(lfit)^2
    D20 <- matrix(unlist(summary(lfit)$varcor),
                  nrow = sqrt(length(unlist(summary(lfit)$varcor))))
    em_fun <- stepEM_Clayton_Gaussian
    has_sigma2 <- TRUE
  } else {
    lfit <- lme4::glmer(formula = Mixeffect, data = Longdat, family = stats::poisson())
    alpha0 <- summary(lfit)$coefficients[, 1]
    sigma20 <- NULL
    D20 <- matrix(unlist(summary(lfit)$varcor),
                  nrow = sqrt(length(unlist(summary(lfit)$varcor))))
    em_fun <- stepEM_Clayton_Poisson
    has_sigma2 <- FALSE
  }

  ## build longitudinal blocks for EM (same as你原逻辑)
  mf_long <- model.frame(lfit)
  y_vec <- model.response(mf_long)
  X_mat <- lme4::getME(lfit, "X")
  id_vec <- mf_long$id

  ni <- as.numeric(table(id_vec))
  y.list <- split(as.numeric(y_vec), id_vec)

  X.df <- data.frame(id = id_vec, X_mat, check.names = FALSE)
  X.list <- by(X.df, X.df$id, function(u) as.matrix(u[, -1, drop = FALSE]))
  n <- length(ni)

  ## random design matrix Z (from random effect terms)
  rand_term <- lme4::findbars(formula(lfit))[[1]]
  rand_rhs <- rand_term[[2]]
  rand_form <- as.formula(paste("~", deparse(rand_rhs)))
  Z_all <- model.matrix(rand_form, data = mf_long)

  Z_df <- data.frame(id = id_vec, Z_all, check.names = FALSE)
  Z.lon <- by(Z_df, Z_df$id, function(u) as.matrix(u[, -1, drop = FALSE]))
  Z.lont <- lapply(Z.lon, t)

  XtX <- lapply(X.list, crossprod)
  XtX.inv <- solve(Reduce("+", XtX))
  Xty <- mapply(function(x, y) crossprod(x, y), x = X.list, y = y.list, SIMPLIFY = FALSE)
  XtZ <- mapply(function(x, z) crossprod(x, z), x = X.list, z = Z.lon, SIMPLIFY = FALSE)

  l <- list(
    ni = ni,
    y = y.list,
    X = X.list,
    Z.lon = Z.lon,
    Z.lont = Z.lont,
    XtX.inv = XtX.inv,
    Xty = Xty,
    XtZ = XtZ
  )

  ## =========================
  ## 2) survival blocks (endpoint 1 & 2)
  ## =========================
  surv1 <- build_surv_block(
    Sur = Surdat1,
    formSurv = formSurv1,
    timevars = timevars1,
    data = Longdat,
    lfit = lfit,
    n = n,
    init_fun = initsSurv_init
  )

  t.1 <- with(surv1$t_list, list(
    W.1 = W,
    Sur1 = Sur,
    Sur1.1 = Sur.1,
    Sur1.1.list = Sur.1.list,
    q.1 = q,
    nev.1 = nev,
    nev.uniq.1 = nev.uq,
    tj.1 = tj,
    maxT.1 = maxT
  ))

  s.1 <- with(surv1$s_list, list(
    Zdat.sur.1 = Zdat.sur,
    time_sur_list_1 = time_sur_list,
    Z.1 = Z,
    Zt.1 = Zt,
    Z.sur.1 = Z.sur,
    Zt.sur.1 = Zt.sur,
    IW.1 = IW
  ))

  su1.para <- surv1$su_para

  surv2 <- build_surv_block(
    Sur = Surdat2,
    formSurv = formSurv2,
    timevars = timevars2,
    data = Longdat,
    lfit = lfit,
    n = n,
    init_fun = initsSurv_init
  )

  t.2 <- with(surv2$t_list, list(
    W.2 = W,
    Sur2 = Sur,
    Sur2.1 = Sur.1,
    Sur2.1.list = Sur.1.list,
    q.2 = q,
    nev.2 = nev,
    nev.uniq.2 = nev.uq,
    tj.2 = tj,
    maxT.2 = maxT
  ))

  s.2 <- with(surv2$s_list, list(
    Zdat.sur.2 = Zdat.sur,
    time_sur_list_2 = time_sur_list,
    Z.2 = Z,
    Zt.2 = Zt,
    Z.sur.2 = Z.sur,
    Zt.sur.2 = Zt.sur,
    IW.2 = IW
  ))

  su2.para <- surv2$su_para

  ## =========================
  ## 3) initial theta
  ## =========================
  theta <- list(
    alpha = alpha0,
    sigma2 = if (has_sigma2) sigma20 else NULL,
    D2 = D20,
    gamma.1 = su1.para$gamma,
    haz.1 = su1.para$haz,
    gamma.2 = su2.para$gamma,
    haz.2 = su2.para$haz,
    rho = rho0
  )

  ## =========================
  ## 4) EM driver
  ## =========================
  em_step <- function(theta, nMC) {
    em_fun(
      theta = theta,
      l = l,
      t.1 = t.1, s.1 = s.1,
      t.2 = t.2, s.2 = s.2,
      nMC = nMC,
      rho.region = rho.region,
      eps_prob = eps_prob,
      verbose = verbose
    )
  }

  em_out <- run_em(
    theta = theta,
    em_step = em_step,
    nMC = nMC,
    tol = tol,
    iter_max = iter_max,
    has_sigma2 = has_sigma2,
    verbose = verbose
  )

  theta <- em_out$theta
  iter <- em_out$iter

  out <- list(
    theta = theta,
    iter = iter,
    long = list(alpha = theta$alpha,
                sigma2 = if (!is.null(theta$sigma2)) theta$sigma2 else NA_real_,
                D2 = theta$D2),
    surv1 = list(gamma = theta$gamma.1, haz = theta$haz.1),
    surv2 = list(gamma = theta$gamma.2, haz = theta$haz.2),
    rho = theta$rho,
    y_type = y_type
  )
  class(out) <- "mjoint_fit"
  out
}
