#' One EM step for Clayton joint model (Gaussian longitudinal)
#'
#' Internal MCEM step for:
#' - Gaussian longitudinal submodel,
#' - two Cox-type survival submodels,
#' - Clayton survival copula dependence.
#'
#' This function is called by `mjoint_Clayton()` via `run_em()`.
#'
#' Notes:
#' - Supports random-effect dimension r >= 1.
#' - Supports fixed time-varying dimension p >= 1.
#' - Uses general Z_sur(t) provided by `build_surv_block()` (from `timevars$rand`).
#'
#' @param theta Current parameter list.
#' @param l Longitudinal block list.
#' @param t.1,s.1 Survival blocks for endpoint 1.
#' @param t.2,s.2 Survival blocks for endpoint 2.
#' @param nMC Monte Carlo sample size.
#' @param rho.region Optimization interval for rho.
#' @param eps_prob Numerical guard (e.g., for log/exp).
#' @param verbose Print debug info if TRUE.
#'
#' @return list(theta = updated_theta)
#' @keywords internal
#' @noRd
stepEM_Clayton_Gaussian <- function(theta, l, t.1, s.1, t.2, s.2,
                                    nMC,
                                    rho.region = c(1e-4, 30),
                                    eps_prob = 1e-300,
                                    verbose = FALSE) {
  ## >>> keep your existing body, but ensure hazard/rho update call:
  ## haz.hat.1 <- update_haz_breslow_Zsur(...)
  ## haz.hat.2 <- update_haz_breslow_Zsur(...)
  ## rho update: Q_rho_Zsur(...)
  ## =========================
  ## 0) unpack theta + detect random-effect dimension r
  ## =========================
  alpha  <- theta$alpha
  sigma2 <- theta$sigma2
  gamma.1 <- theta$gamma.1
  gamma.2 <- theta$gamma.2
  haz.1 <- theta$haz.1
  haz.2 <- theta$haz.2
  rho   <- theta$rho
  
  # D matrix name can be D2 (old one-dim version) or Sigma.b (two-dim version)
  if (!is.null(theta$Sigma.b)) {
    Sigma.b <- theta$Sigma.b
  } else if (!is.null(theta$D2)) {
    Sigma.b <- theta$D2
  } else {
    stop("theta must contain either Sigma.b or D2.")
  }
  r_dim <- ncol(Sigma.b)
  
  ## =========================
  ## 1) unpack longitudinal blocks
  ## =========================
  ni <- l$ni
  X  <- l$X
  y  <- l$y
  XtX.inv <- l$XtX.inv
  Xty <- l$Xty
  XtZ <- l$XtZ
  Z.lon  <- l$Z.lon
  Z.lont <- l$Z.lont
  
  n <- length(ni)
  
  ## =========================
  ## 2) unpack survival blocks (naming consistent with your new build_surv_block pipeline)
  ## =========================
  # --- T1
  W.1 <- t.1$W.1
  Sur1 <- t.1$Sur1
  Sur1.1 <- t.1$Sur1.1
  Sur1.1.list <- t.1$Sur1.1.list
  q.1 <- t.1$q.1
  tj.1 <- t.1$tj.1
  nev.1 <- t.1$nev.1
  nev.uniq.1 <- t.1$nev.uniq.1
  
  Zdat.sur.1 <- s.1$Zdat.sur.1
  Z.sur.1 <- s.1$Z.sur.1
  Zt.sur.1 <- s.1$Zt.sur.1
  Z.1 <- s.1$Z.1
  Zt.1 <- s.1$Zt.1
  IW.1 <- s.1$IW.1
  time_sur_list_1 <- if (!is.null(s.1$time_sur_list_1)) s.1$time_sur_list_1 else split(Zdat.sur.1$time, Zdat.sur.1$id)
  
  # --- T2
  W.2 <- t.2$W.2
  Sur2 <- t.2$Sur2
  Sur2.1 <- t.2$Sur2.1
  Sur2.1.list <- t.2$Sur2.1.list
  q.2 <- t.2$q.2
  tj.2 <- t.2$tj.2
  nev.2 <- t.2$nev.2
  nev.uniq.2 <- t.2$nev.uniq.2
  
  Zdat.sur.2 <- s.2$Zdat.sur.2
  Z.sur.2 <- s.2$Z.sur.2
  Zt.sur.2 <- s.2$Zt.sur.2
  Z.2 <- s.2$Z.2
  Zt.2 <- s.2$Zt.2
  IW.2 <- s.2$IW.2
  time_sur_list_2 <- if (!is.null(s.2$time_sur_list_2)) s.2$time_sur_list_2 else split(Zdat.sur.2$time, Zdat.sur.2$id)
  
  # fixed part dimension p_k is inferred from Z.* columns
  p.1 <- ncol(as.matrix(Z.1[[1]]))
  p.2 <- ncol(as.matrix(Z.2[[1]]))
  
  ## =========================
  ## 3) Monte Carlo sampling b | y
  ## =========================
  Sigmai.inv <- lapply(ni, function(ii) diag(x = rep(1 / sigma2, ii), ncol = ii))
  
  Dinv <- solve(Sigma.b)
  Ai <- mapply(
    FUN = function(zt, s, z) solve((zt %*% s %*% z) + Dinv),
    z = Z.lon, zt = Z.lont, s = Sigmai.inv,
    SIMPLIFY = FALSE
  )
  
  Mi <- mapply(
    function(a, zt, s, yi, xi) as.vector(a %*% (zt %*% s %*% (yi - xi %*% alpha))),
    a = Ai, zt = Z.lont, s = Sigmai.inv, yi = y, xi = X,
    SIMPLIFY = FALSE
  )
  
  Zq <- randtoolbox::sobol(nMC, dim = r_dim, normal = TRUE, scrambling = 1)
  bi.y <- mapply(
    function(m, a) {
      C <- chol(a)
      matrix(rep(m, nMC), nrow = nMC, byrow = TRUE) + (Zq %*% C)
    },
    m = Mi, a = Ai,
    SIMPLIFY = FALSE
  )
  names(bi.y) <- names(Ai)
  
  ## =========================
  ## 4) Cox pieces under current gamma
  ## =========================
  # W^T gamma_{k1}
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  
  # exp{ Z^*(t) gamma_{k2} } : Z^*(t) is p_k dim
  expZ.1 <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1 + 1):(q.1 + p.1)])), z = Z.1, SIMPLIFY = FALSE)
  expZ.2 <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2 + 1):(q.2 + p.2)])), z = Z.2, SIMPLIFY = FALSE)
  
  # association parameter(s) beta_k : K=1 so scalar, located after q+p
  beta1 <- as.numeric(gamma.1[-(1:(q.1 + p.1))])
  beta2 <- as.numeric(gamma.2[-(1:(q.2 + p.2))])
  
  # expWArma scaling matrix (r_dim x r_dim)
  gam.1.scale <- diag(rep(beta1, r_dim), r_dim, r_dim)
  gam.2.scale <- diag(rep(beta2, r_dim), r_dim, r_dim)
  
  # build IZ.sur lists for expWArma
  IZ.sur.1 <- mapply(function(x, y) t(x %*% y), x = IW.1, y = Z.sur.1, SIMPLIFY = FALSE)
  IZ.sur.2 <- mapply(function(x, y) t(x %*% y), x = IW.2, y = Z.sur.2, SIMPLIFY = FALSE)
  
  expZ.surb.1 <- expWArma(IZ.sur.1, bi.y, gam.1.scale, Sur1.1.list)
  expZ.surb.2 <- expWArma(IZ.sur.2, bi.y, gam.2.scale, Sur2.1.list)
  
  ## =========================
  ## 5) posterior weights pb.yt via joint survival likelihood
  ## =========================
  logf_sur <- mapply(function(w1, z1, zsur1, h1, w2, z2, zsur2, h2) {
    
    H1 <- as.vector(t(t(zsur1) * z1) %*% haz.1[1:ncol(zsur1)]) * exp(w1)
    H2 <- as.vector(t(t(zsur2) * z2) %*% haz.2[1:ncol(zsur2)]) * exp(w2)
    
    if (h1$delta == 1 && h2$delta == 1) {
      (log(haz.1[ncol(zsur1)]) + w1 + log(z1[ncol(zsur1)]) + log(zsur1[, ncol(zsur1)])) - H1 +
        (log(haz.2[ncol(zsur2)]) + w2 + log(z2[ncol(zsur2)]) + log(zsur2[, ncol(zsur2)])) - H2 +
        log((rho + 1) * (exp(-H1) * exp(-H2)) ^ (-(rho + 1)) *
              (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(2 * rho + 1) / rho))
    } else if (h1$delta == 1 && h2$delta == 0) {
      (log(haz.1[ncol(zsur1)]) + w1 + log(z1[ncol(zsur1)]) + log(zsur1[, ncol(zsur1)])) - H1 +
        log(exp(-H1)^(-(rho + 1)) * (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(rho + 1) / rho))
    } else if (h1$delta == 0 && h2$delta == 1) {
      (log(haz.2[ncol(zsur2)]) + w2 + log(z2[ncol(zsur2)]) + log(zsur2[, ncol(zsur2)])) - H2 +
        log(exp(-H2)^(-(rho + 1)) * (exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1) ^ (-(rho + 1) / rho))
    } else {
      log((exp(-H1)^(-rho) + exp(-H2)^(-rho) - 1)^(-1 / rho))
    }
  },
  w1 = Wtgam.1, z1 = expZ.1, zsur1 = expZ.surb.1, h1 = Sur1.1.list,
  w2 = Wtgam.2, z2 = expZ.2, zsur2 = expZ.surb.2, h2 = Sur2.1.list,
  SIMPLIFY = FALSE)
  
  f_pdf <- lapply(logf_sur, exp)
  den <- lapply(f_pdf, mean)
  
  for (k in which(unlist(den) == 0)) {
    f_pdf[[k]] <- rep(1, length(f_pdf[[1]]))
    den[[k]] <- 1
  }
  
  pb.yt <- mapply(function(f, d) {
    if (d == 1 || is.nan(d) || is.na(d) || is.infinite(d)) rep(1, length(f)) else f / d
  }, f = f_pdf, d = den, SIMPLIFY = FALSE)
  
  ## =========================
  ## 6) E-step: Eb, EbbT
  ## =========================
  Eb <- mapply(function(b, pb) colMeans(b * pb), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  EbbT <- mapply(function(b, pb) crossprod(b, (b * pb)) / nrow(b), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  
  expvstargam.1 <- mapply(function(w1, z1, zsur1) exp(w1) * t(t(zsur1) * z1),
                          w1 = Wtgam.1, z1 = expZ.1, zsur1 = expZ.surb.1, SIMPLIFY = FALSE)
  expvstargam.2 <- mapply(function(w2, z2, zsur2) exp(w2) * t(t(zsur2) * z2),
                          w2 = Wtgam.2, z2 = expZ.2, zsur2 = expZ.surb.2, SIMPLIFY = FALSE)
  
  ## =========================
  ## 7) hazard update (E-step)
  ## =========================
  haz.hat.1 <- update_haz_breslow_Zsur(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    time_sur_list = s.1$time_sur_list_1,
    Zsur_list = s.1$Z.sur.1,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta1
  )
  
  haz.hat.2 <- update_haz_breslow_Zsur(
    tj = tj.2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    time_sur_list = s.2$time_sur_list_2,
    Zsur_list = s.2$Z.sur.2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  ## =========================
  ## 8) gamma update (E-step) — safe
  ## =========================
  g1 <- gammaUpdate_fixDelta_general(
    bi.y, Zt.sur.1, expvstargam.1, pb.yt, haz.hat.1,
    W.1, Zt.1, Sur1.1.list,
    K = 1, q = q.1, qq = p.1, nev = nev.uniq.1, jcount = nev.1
  )
  gDelta.1 <- g1$gDelta
  
  g2 <- gammaUpdate_fixDelta_general(
    bi.y, Zt.sur.2, expvstargam.2, pb.yt, haz.hat.2,
    W.2, Zt.2, Sur2.1.list,
    K = 1, q = q.2, qq = p.2, nev = nev.uniq.2, jcount = nev.2
  )
  gDelta.2 <- g2$gDelta
  ## =========================
  ## 9) M-step updates
  ## =========================
  Sigma.b.new <- Reduce("+", EbbT) / n
  rownames(Sigma.b.new) <- colnames(Sigma.b.new) <- rownames(Sigma.b)
  
  rr <- mapply(function(x1, x2, b) x1 - (x2 %*% b), x1 = Xty, x2 = XtZ, b = Eb)
  rr.sum <- rowSums(rr)
  alpha.new <- as.vector(XtX.inv %*% rr.sum)
  names(alpha.new) <- names(alpha)
  
  SSq <- mapply(function(yi, xi, zi, b, b2) {
    residFixed <- (yi - xi %*% alpha.new)
    t(residFixed) %*% (residFixed - 2 * (zi %*% b)) + sum(diag(crossprod(zi) %*% b2))
  }, yi = y, xi = X, zi = Z.lon, b = Eb, b2 = EbbT)
  
  sigma2.new <- sum(SSq) / sum(unlist(ni))
  
  gamma.1.new <- gamma.1 + gDelta.1
  gamma.2.new <- gamma.2 + gDelta.2
  
  ## =========================
  ## 10) update haz under NEW gamma, then update rho
  ## =========================
  gamma.1 <- gamma.1.new; gamma.2 <- gamma.2.new
  
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  expZ.1  <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)])), z = Z.1, SIMPLIFY = FALSE)
  beta1   <- as.numeric(gamma.1[-(1:(q.1+p.1))])
  
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  expZ.2  <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)])), z = Z.2, SIMPLIFY = FALSE)
  beta2   <- as.numeric(gamma.2[-(1:(q.2+p.2))])
  
  haz.1.new <- update_haz_breslow_Zsur(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ = expZ.1,
    time_sur_list = s.1$time_sur_list_1,
    Zsur_list = s.1$Z.sur.1,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta1
  )
  
  haz.2.new <- update_haz_breslow_Zsur(
    tj = tj.2,
    Sur1.1 = Sur2.1,
    Wtgam = Wtgam.2,
    expZ = expZ.2,
    time_sur_list = s.2$time_sur_list_2,
    Zsur_list = s.2$Z.sur.2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  # rho update
  log_rho <- function(rho) {
    Q_rho_Zsur(
      rho = rho,
      Sur1.1 = Sur1.1, Sur2.1 = Sur2.1,
      tj.1 = tj.1, tj.2 = tj.2,
      haz.1 = haz.1.new, haz.2 = haz.2.new,
      Wtgam.1 = Wtgam.1, Wtgam.2 = Wtgam.2,
      expZ.1 = expZ.1, expZ.2 = expZ.2,
      time_sur_list_1 = time_sur_list_1,
      time_sur_list_2 = time_sur_list_2,
      Zsur_list_1 = Z.sur.1,
      Zsur_list_2 = Z.sur.2,
      b_y = bi.y, pb_y = pb.yt,
      beta1 = beta1, beta2 = beta2
    )
  }
  rho.new <- tryCatch(optimize(log_rho, rho.region)$minimum, error = function(e) rho)
  
  ## =========================
  ## 11) return theta.new (keep original field names)
  ## =========================
  theta.new <- theta
  theta.new$alpha <- alpha.new
  theta.new$sigma2 <- sigma2.new
  theta.new$gamma.1 <- gamma.1.new
  theta.new$gamma.2 <- gamma.2.new
  theta.new$haz.1 <- haz.1.new
  theta.new$haz.2 <- haz.2.new
  theta.new$rho <- rho.new
  
  if (!is.null(theta$Sigma.b)) theta.new$Sigma.b <- Sigma.b.new
  if (!is.null(theta$D2))      theta.new$D2 <- Sigma.b.new
  
  list(theta = theta.new)
}