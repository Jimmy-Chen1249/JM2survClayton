#' Initial values for the survival submodel (Cox + AG approximation)
#'
#' This file provides an initialization routine for the survival part of the joint model.
#' The function `initsSurv_init()` is designed to be passed as `init_fun` into
#' `build_surv_block()`.
#'
#' Strategy:
#' - Construct an Andersen-Gill (start, stop, status) dataset using longitudinal visit times.
#' - Add baseline covariates W (from Cox design), time-varying fixed covariate Z^*(t),
#'   and a proxy of the shared random effect from the fitted longitudinal model (lfit).
#' - Fit a Cox model on this AG dataset to get initial gamma and baseline hazard.
#'
#' IMPORTANT:
#' - This routine does NOT assume any specific Z_sur(t) form. It only needs Z^*(t) (fixed part).
#' - For stability, we add a small epsilon to the terminal time.
#'
#' @keywords internal
#' @noRd
NULL

## -------------------------------------------------------------------
## helper: safe jitter for tied survival times (optional external use)
## -------------------------------------------------------------------
#' @keywords internal
#' @noRd
jitter_surtime <- function(surtime, eps = 1e-4) {
  surtime <- as.numeric(surtime)
  dup_idx <- duplicated(surtime) | duplicated(surtime, fromLast = TRUE)
  jitter_idx <- dup_idx & !is.na(surtime)
  if (any(jitter_idx)) {
    surtime[jitter_idx] <- surtime[jitter_idx] + stats::rnorm(sum(jitter_idx), 0, eps)
  }
  surtime[!is.na(surtime)] <- pmax(surtime[!is.na(surtime)], .Machine$double.eps)
  surtime
}

## -------------------------------------------------------------------
## Main init_fun: compatible with build_surv_block(...)
## -------------------------------------------------------------------
#' @keywords internal
#' @noRd
initsSurv_init <- function(
    Sur,
    formSurv,
    timevars,
    Sur.1,
    time_sur_list,
    data,
    lfit,
    q,
    eps = 1e-6,
    ...
) {
  ## This signature MUST tolerate extra args via `...`.
  
  ## Proxy for shared random effects from the longitudinal fit:
  ## b_proxy = (Xbeta + Zb) - Xbeta
  XbetaZb <- stats::predict(lfit, re.form = NULL)
  Xbeta   <- stats::predict(lfit, re.form = NA)
  b_proxy <- as.numeric(XbetaZb - Xbeta)
  
  ## align b_proxy with `data$id`
  b_df <- data.frame(id = data$id, b = b_proxy, row.names = NULL)
  
  ## Build AG data by subject
  ## We rely on longitudinal visit times in `data$time` by default.
  ## If your longitudinal time variable has different name, pass it via ... and use it here.
  timeVar <- if ("time" %in% names(data)) "time" else stop("Longitudinal data must contain 'time' column for init.")
  
  dataAG_list <- by(data, data$id, FUN = function(u) {
    sid <- u$id[1]
    
    Ti <- Sur.1[Sur.1$id == sid, "T"] + eps
    di <- Sur.1[Sur.1$id == sid, "delta"]
    
    start_all <- u[[timeVar]]
    keep <- which(start_all < Ti)
    if (length(keep) == 0) return(NULL)
    
    id.t <- max(keep)
    
    start <- start_all[1:id.t]
    stop  <- c(start_all[2:id.t], Ti)
    status <- rep(0, id.t)
    status[id.t] <- di
    
    ## baseline W (q columns): take from Sur.1 row
    if (q > 0) {
      Wmat <- Sur.1[Sur.1$id == sid, 2:(q + 1), drop = FALSE]
      Wmat <- as.matrix(Wmat)
      Wmat <- Wmat[rep(1, id.t), , drop = FALSE]
      colnames(Wmat) <- names(Sur.1)[2:(q + 1)]
    }
    
    ## time-varying fixed covariate Z^*(t): use timevars$fixed
    Zf <- timevars$fixed(start, sid, Sur = Sur, ...)
    Zf <- as.matrix(Zf)
    
    ## proxy of shared random effect: match rows 1:id.t
    b1 <- b_df$b[b_df$id == sid]
    b1 <- b1[1:id.t]
    
    ## assemble
    if (q > 0) {
      out <- data.frame(
        id = rep(sid, id.t),
        start = start,
        stop = stop,
        status = status,
        Wmat,
        Zf,
        beta_1 = b1,
        row.names = NULL,
        check.names = FALSE
      )
    } else {
      out <- data.frame(
        id = rep(sid, id.t),
        start = start,
        stop = stop,
        status = status,
        Zf,
        beta_1 = b1,
        row.names = NULL,
        check.names = FALSE
      )
    }
    out
  })
  
  dataAG <- do.call(rbind, dataAG_list)
  if (is.null(dataAG) || nrow(dataAG) == 0) stop("initsSurv_init: AG data is empty.")
  
  ## Cox formula: baseline W + Z^*(t) + beta_1 proxy
  z_names <- setdiff(colnames(dataAG), c("id","start","stop","status","beta_1"))
  rhs <- c(z_names, "beta_1")
  form_full <- stats::as.formula(paste0("survival::Surv(start, stop, status) ~ ", paste(rhs, collapse = " + ")))
  
  fit <- survival::coxph(form_full, data = dataAG, singular.ok = TRUE)
  
  gamma <- stats::coef(fit)
  haz <- survival::coxph.detail(fit)$hazard
  
  list(gamma = gamma, haz = haz)
}