#' Build survival blocks used by the EM algorithm
#'
#' This is an internal helper that converts one survival endpoint into the
#' list-structure required by the EM step:
#' - baseline covariates (W) from a Cox model,
#' - event-time grid tj and subject-specific grid time list,
#' - time-varying fixed covariates Z^*(t) (entering gamma_{k2}),
#' - time-varying random-effect design Z_sur(t) (entering beta_k),
#' - IW matrices (diagonal selector matrices),
#' - initial values (gamma, baseline hazard) via an `init_fun`.
#'
#' IMPORTANT:
#' - This function does NOT assume any specific form of Z_sur(t).
#'   The user provides it through `timevars$rand`.
#' - It supports arbitrary fixed dimension p (timevars$fixed can return p>1 columns),
#'   and arbitrary random-effect dimension r (timevars$rand returns r columns).
#'
#' @param Sur A data.frame for one endpoint, containing:
#'   id column (default "id"), survival time column (default "surtime"),
#'   event indicator column (default "cens"), plus any baseline covariates in `formSurv`.
#'   It may also contain B1/B2/B3 etc. that are used inside `timevars$fixed`.
#' @param formSurv A survival formula `Surv(time,event) ~ covariates` for `survival::coxph`.
#' @param timevars A named list with two functions:
#'   - fixed(t, index, Sur = ..., ...) -> matrix (length(t) x p)
#'   - rand(t, index, Sur = ..., ...)  -> matrix (length(t) x r)
#'   Both functions are evaluated per subject on that subject's grid.
#' @param data Longitudinal data (passed to init_fun if needed).
#' @param lfit Fitted longitudinal mixed model object (passed to init_fun if needed).
#' @param n Number of subjects.
#' @param id_name Subject id column name in `Sur`. Default "id".
#' @param time_name_sur Survival time column name in `Sur`. Default "surtime".
#' @param event_name_sur Event indicator column name in `Sur`. Default "cens".
#' @param jitter_ties Logical; if TRUE, jitter tied survival times slightly.
#' @param jitter_eps Std dev of jitter noise.
#' @param init_fun Function that computes initial `(gamma, haz)` for this endpoint.
#'   It must accept `...` and ignore unused arguments safely.
#'   See `initsSurv_init()` in `R/initsSurv.R`.
#' @param ... Extra arguments forwarded to `timevars$fixed`, `timevars$rand`,
#'   and `init_fun`.
#'
#' @return A list with components:
#' - t_list: W, Sur, Sur.1, Sur.1.list, q, nev, nev.uq, tj, maxT, timevars
#' - s_list: Zdat.sur, time_sur_list, Z, Zt, Z.sur, Zt.sur, IW
#' - su_para: list(gamma, haz)
#'
#' @keywords internal
#' @noRd
build_surv_block <- function(
    Sur,
    formSurv,
    timevars,
    data,
    lfit,
    n,
    id_name = "id",
    time_name_sur = "surtime",
    event_name_sur = "cens",
    jitter_ties = TRUE,
    jitter_eps = 1e-6,
    init_fun = NULL,
    ...
) {
  ## ---------- basic checks ----------
  stopifnot(is.data.frame(Sur))
  stopifnot(inherits(formSurv, "formula"))
  
  if (!is.list(timevars) || is.null(timevars$fixed) || is.null(timevars$rand)) {
    stop("timevars must be a named list: list(fixed = FUN, rand = FUN).")
  }
  if (!is.function(timevars$fixed) || !is.function(timevars$rand)) {
    stop("timevars$fixed and timevars$rand must be functions.")
  }
  
  if (!id_name %in% names(Sur)) stop("Sur must contain column: ", id_name)
  if (!time_name_sur %in% names(Sur)) stop("Sur must contain column: ", time_name_sur)
  if (!event_name_sur %in% names(Sur)) stop("Sur must contain column: ", event_name_sur)
  
  ## enforce ordering by id (important for consistent indexing)
  Sur <- Sur[order(Sur[[id_name]]), , drop = FALSE]
  
  ## ---------- optional jitter on tied times ----------
  if (isTRUE(jitter_ties)) {
    tt <- as.numeric(Sur[[time_name_sur]])
    dup <- duplicated(tt) | duplicated(tt, fromLast = TRUE)
    if (any(dup)) {
      tt[dup] <- tt[dup] + stats::rnorm(sum(dup), mean = 0, sd = jitter_eps)
    }
    tt <- pmax(tt, .Machine$double.eps)
    Sur[[time_name_sur]] <- tt
  }
  
  ## ---------- Cox fit to get baseline covariates W ----------
  sfit <- survival::coxph(formSurv, data = Sur, x = TRUE)
  q <- ncol(sfit$x)
  
  sfit.st <- survival::survfit(sfit)
  tj  <- sfit.st$time[sfit.st$n.event > 0]
  nev <- sfit.st$n.event[sfit.st$n.event > 0]
  nev.uq <- length(tj)
  
  ## Sur.1: one row per subject
  Sur.1 <- data.frame(
    id    = Sur[[id_name]],
    sfit$x,
    T     = sfit$y[, 1],
    delta = sfit$y[, 2],
    check.names = FALSE
  )
  if (q > 0) names(Sur.1)[2:(q + 1)] <- colnames(sfit$x)
  
  ## how many grid points tj <= T_i
  Sur.1$tj.ind <- vapply(Sur.1$T, function(Ti) sum(tj <= Ti), numeric(1))
  
  ## list per subject
  Sur.1.list <- split(Sur.1, Sur.1$id)
  
  ## W list: each element is numeric vector length q
  W <- lapply(Sur.1.list, function(df) {
    if (q == 0) return(numeric(0))
    as.numeric(df[1, 2:(q + 1), drop = TRUE])
  })
  
  ## subject with largest observed EVENT time (fallback if no event)
  if (any(Sur.1$delta == 1)) {
    maxT <- which.max(ifelse(Sur.1$delta == 1, Sur.1$T, -Inf))
  } else {
    maxT <- which.max(Sur.1$T)
  }
  
  ## ---------- build subject-specific time grids ----------
  ## Each subject i uses tj[1:nj_i], where nj_i = max(tj.ind_i, 1)
  nj_i <- pmax(Sur.1$tj.ind, 1L)
  
  Zdat.sur <- data.frame(
    id = rep(Sur.1$id, nj_i),
    time = unlist(lapply(seq_len(nrow(Sur.1)), function(i) tj[seq_len(nj_i[i])])),
    row.names = NULL
  )
  
  ## per subject time vector
  time_sur_list <- split(Zdat.sur$time, Zdat.sur$id)
  
  ## ---------- build Z(t) and Z_sur(t) from user functions ----------
  Z <- vector("list", length = n)
  Z.sur <- vector("list", length = n)
  
  for (i in seq_len(n)) {
    id_i <- as.character(i)
    tvec <- time_sur_list[[id_i]]
    if (is.null(tvec) || length(tvec) == 0) {
      ## should not happen, but be defensive
      tvec <- tj[1]
      time_sur_list[[id_i]] <- tvec
    }
    
    ## IMPORTANT: pass Sur=Sur so user functions can access B1/B2/B3 etc.
    Zf <- timevars$fixed(tvec, i, Sur = Sur, ...)
    Zr <- timevars$rand(tvec, i, Sur = Sur, ...)
    
    Zf <- as.matrix(Zf)
    Zr <- as.matrix(Zr)
    
    if (nrow(Zf) != length(tvec)) {
      stop(sprintf("timevars$fixed must return %d rows (id=%d).", length(tvec), i))
    }
    if (nrow(Zr) != length(tvec)) {
      stop(sprintf("timevars$rand must return %d rows (id=%d).", length(tvec), i))
    }
    
    Z[[i]] <- Zf
    Z.sur[[i]] <- Zr
  }
  names(Z) <- as.character(seq_len(n))
  names(Z.sur) <- as.character(seq_len(n))
  
  Zt <- lapply(Z, t)
  Zt.sur <- lapply(Z.sur, t)
  
  ## ---------- IW (diagonal matrices) ----------
  IW <- lapply(seq_len(n), function(i) diag(length(time_sur_list[[as.character(i)]])))
  names(IW) <- as.character(seq_len(n))
  
  ## ---------- initial values: gamma + baseline hazard increments ----------
  ## default init: Cox with baseline covariates only
  if (is.null(init_fun)) {
    init_fun <- function(Sur, formSurv, ...) {
      fit0 <- survival::coxph(formSurv, data = Sur, x = TRUE)
      list(gamma = stats::coef(fit0),
           haz   = survival::coxph.detail(fit0)$hazard)
    }
  }
  
  ## The init_fun MUST tolerate extra arguments via `...`.
  su.para <- init_fun(
    Sur = Sur,
    formSurv = formSurv,
    timevars = timevars,
    Sur.1 = Sur.1,
    time_sur_list = time_sur_list,
    data = data,
    lfit = lfit,
    q = q,
    ...
  )
  
  ## ---------- return ----------
  list(
    t_list = list(
      W          = W,
      Sur        = Sur,
      Sur.1      = Sur.1,
      Sur.1.list = Sur.1.list,
      q          = q,
      nev        = nev,
      nev.uq     = nev.uq,
      tj         = tj,
      maxT       = maxT,
      timevars   = timevars
    ),
    s_list = list(
      Zdat.sur      = Zdat.sur,
      time_sur_list = time_sur_list,
      Z             = Z,
      Zt            = Zt,
      Z.sur         = Z.sur,
      Zt.sur        = Zt.sur,
      IW            = IW
    ),
    su_para = su.para
  )
}