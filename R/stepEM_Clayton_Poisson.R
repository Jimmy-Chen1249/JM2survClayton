#' One EM step for Clayton joint model (Poisson longitudinal)
#'
#' Internal MCEM step for:
#' - Poisson longitudinal submodel (GLMM),
#' - two Cox-type survival submodels,
#' - Clayton survival copula dependence.
#'
#' @inheritParams stepEM_Clayton_Gaussian
#' @return list(theta = updated_theta)
#' @keywords internal
#' @noRd
stepEM_Clayton_Poisson <- function(theta, l, t.1, s.1, t.2, s.2,
                                   nMC,
                                   rho.region = c(1e-4, 30),
                                   eps_prob = 1e-300,
                                   verbose = FALSE) {
  ## =========================
  ## 0) unpack theta + detect random-effect dimension r
  ## =========================
  alpha   <- theta$alpha
  gamma.1 <- theta$gamma.1
  gamma.2 <- theta$gamma.2
  haz.1   <- theta$haz.1
  haz.2   <- theta$haz.2
  rho     <- theta$rho
  
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
  ni      <- l$ni
  X       <- l$X
  y       <- l$y
  XtX.inv <- l$XtX.inv
  Xty     <- l$Xty
  XtZ     <- l$XtZ
  Z.lon   <- l$Z.lon
  Z.lont  <- l$Z.lont
  n       <- length(ni)
  
  ## =========================
  ## 2) unpack survival blocks
  ## =========================
  ## --- T1
  W.1         <- t.1$W.1
  Sur1        <- t.1$Sur1
  Sur1.1      <- t.1$Sur1.1
  Sur1.1.list <- t.1$Sur1.1.list
  q.1         <- t.1$q.1
  tj.1        <- t.1$tj.1
  nev.1       <- t.1$nev.1
  nev.uniq.1  <- t.1$nev.uniq.1
  
  Zdat.sur.1 <- s.1$Zdat.sur.1
  Z.sur.1    <- s.1$Z.sur.1
  Zt.sur.1   <- s.1$Zt.sur.1
  Z.1        <- s.1$Z.1
  Zt.1       <- s.1$Zt.1
  IW.1       <- s.1$IW.1
  time_sur_list_1 <- if (!is.null(s.1$time_sur_list_1)) s.1$time_sur_list_1 else split(Zdat.sur.1$time, Zdat.sur.1$id)
  
  ## --- T2
  W.2         <- t.2$W.2
  Sur2        <- t.2$Sur2
  Sur2.1      <- t.2$Sur2.1
  Sur2.1.list <- t.2$Sur2.1.list
  q.2         <- t.2$q.2
  tj.2        <- t.2$tj.2
  nev.2       <- t.2$nev.2
  nev.uniq.2  <- t.2$nev.uniq.2
  
  Zdat.sur.2 <- s.2$Zdat.sur.2
  Z.sur.2    <- s.2$Z.sur.2
  Zt.sur.2   <- s.2$Zt.sur.2
  Z.2        <- s.2$Z.2
  Zt.2       <- s.2$Zt.2
  IW.2       <- s.2$IW.2
  time_sur_list_2 <- if (!is.null(s.2$time_sur_list_2)) s.2$time_sur_list_2 else split(Zdat.sur.2$time, Zdat.sur.2$id)
  
  ## fixed part dimension p_k inferred from Z.* columns
  p.1 <- ncol(as.matrix(Z.1[[1]]))
  p.2 <- ncol(as.matrix(Z.2[[1]]))
  
  ## =========================
  ## 3) Monte Carlo sampling b | y  (Poisson)
  ## =========================
  make_pd <- function(S) {
    S <- 0.5 * (S + t(S))
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 1e-8) S <- S + diag(1e-6 - pmin(0, min(ev)) + 1e-8, ncol(S))
    S
  }
  Sigma.prop <- make_pd(Sigma.b)
  
  ## Xbeta_i = X_i %*% alpha  (vector length n_i)
  Xbeta_list <- lapply(X, function(Xi) drop(Xi %*% alpha))
  
  tune_sd0 <- function(Xbeta_i, Zi, yi, Sigma_i, u0,
                       sd0 = 1.0, B_tune = 400, max_tune = 8,
                       ar_lo = 0.15, ar_hi = 0.40) {
    for (kk in seq_len(max_tune)) {
      out  <- uSamplerPoissonCpp_n_acc(Xbeta_i, Zi, u0, yi, Sigma_i, B = B_tune, sd0 = sd0)
      ar   <- out$acc_rate
      u0   <- out$sample[nrow(out$sample), ]
      if (!is.finite(ar)) break
      if (ar < ar_lo) { sd0 <- 0.8 * sd0; next }
      if (ar > ar_hi) { sd0 <- 1.2 * sd0; next }
      break
    }
    list(sd0 = sd0, u0 = u0)
  }
  
  bi.y <- mapply(function(yi, Xbeta_i, Zi) {
    q  <- ncol(Zi)
    u0 <- mvtnorm::rmvnorm(1, rep(0, q), Sigma.prop)[1, ]
    tun <- tune_sd0(Xbeta_i, Zi, yi, Sigma.prop, u0, sd0 = 1.0, B_tune = 400)
    out <- uSamplerPoissonCpp_n_acc(Xbeta_i, Zi, tun$u0, yi, Sigma.prop, B = nMC, sd0 = tun$sd0)
    out$sample
  }, yi = y, Xbeta_i = Xbeta_list, Zi = Z.lon, SIMPLIFY = FALSE)
  
  ## =========================
  ## 4) Cox pieces under current gamma
  ## =========================
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  
  expZ.1 <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1 + 1):(q.1 + p.1)])), z = Z.1, SIMPLIFY = FALSE)
  expZ.2 <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2 + 1):(q.2 + p.2)])), z = Z.2, SIMPLIFY = FALSE)
  
  beta1 <- as.numeric(gamma.1[-(1:(q.1 + p.1))])
  beta2 <- as.numeric(gamma.2[-(1:(q.2 + p.2))])
  
  gam.1.scale <- diag(rep(beta1, r_dim), r_dim, r_dim)
  gam.2.scale <- diag(rep(beta2, r_dim), r_dim, r_dim)
  
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
  Eb   <- mapply(function(b, pb) colMeans(b * pb), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  EbbT <- mapply(function(b, pb) crossprod(b, (b * pb)) / nrow(b), b = bi.y, pb = pb.yt, SIMPLIFY = FALSE)
  
  expvstargam.1 <- mapply(function(w1, z1, zsur1) exp(w1) * t(t(zsur1) * z1),
                          w1 = Wtgam.1, z1 = expZ.1, zsur1 = expZ.surb.1, SIMPLIFY = FALSE)
  expvstargam.2 <- mapply(function(w2, z2, zsur2) exp(w2) * t(t(zsur2) * z2),
                          w2 = Wtgam.2, z2 = expZ.2, zsur2 = expZ.surb.2, SIMPLIFY = FALSE)
  
  ## =========================
  ## 7) hazard update (E-step) — general Zsur
  ## =========================
  haz.hat.1 <- update_haz_breslow_Zsur(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ  = expZ.1,
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
    expZ  = expZ.2,
    time_sur_list = s.2$time_sur_list_2,
    Zsur_list = s.2$Z.sur.2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  ## =========================
  ## 8) gamma update (E-step) — safe (general)
  ## =========================
  g1 <- gammaUpdate_fixDelta_general(
    bi.y, Zt.sur.1, expvstargam.1, pb.yt, haz.hat.1,
    W.1, Zt.1, Sur1.1.list,
    K = 1, q = q.1, qq = p.1, nev = nev.uniq.1, jcount = nev.1
  )
  gDelta.1 <- as.vector(g1$gDelta)
  
  g2 <- gammaUpdate_fixDelta_general(
    bi.y, Zt.sur.2, expvstargam.2, pb.yt, haz.hat.2,
    W.2, Zt.2, Sur2.1.list,
    K = 1, q = q.2, qq = p.2, nev = nev.uniq.2, jcount = nev.2
  )
  gDelta.2 <- as.vector(g2$gDelta)
  
  ## =========================
  ## 9) M-step updates
  ## =========================
  ## Sigma.b
  Sigma.b.new <- Reduce("+", EbbT) / n
  rownames(Sigma.b.new) <- colnames(Sigma.b.new) <- rownames(Sigma.b)
  
  ## alpha (Poisson NR)
  Eexpb <- mapply(function(b, pb, xi, zi) {
    exp(xi %*% alpha) * apply(exp(b %*% t(zi)) * pb, 2, mean)
  }, b = bi.y, pb = pb.yt, xi = X, zi = Z.lon, SIMPLIFY = FALSE)
  
  ScoreAlpha <- mapply(function(yi, xi, eexpb) {
    colSums(yi * xi - xi * as.numeric(eexpb))
  }, yi = y, xi = X, eexpb = Eexpb, SIMPLIFY = FALSE)
  
  dScoreAlpha <- mapply(function(xi, eexpb) {
    Reduce("+", lapply(seq_len(nrow(xi)), function(i) -xi[i, ] %*% t(xi[i, ]) * as.numeric(eexpb)[i]))
  }, xi = X, eexpb = Eexpb, SIMPLIFY = FALSE)
  
  alpha.new <- as.vector(alpha - MASS::ginv(Reduce("+", dScoreAlpha)) %*% Reduce("+", ScoreAlpha))
  names(alpha.new) <- names(alpha)
  
  gamma.1.new <- gamma.1 + gDelta.1
  gamma.2.new <- gamma.2 + gDelta.2
  
  ## =========================
  ## 10) update haz under NEW gamma, then update rho
  ## =========================
  gamma.1 <- gamma.1.new
  gamma.2 <- gamma.2.new
  
  Wtgam.1 <- mapply(function(w) as.numeric(w %*% gamma.1[1:q.1]), w = W.1, SIMPLIFY = FALSE)
  expZ.1  <- mapply(function(z) exp(as.numeric(z %*% gamma.1[(q.1+1):(q.1+p.1)])), z = Z.1, SIMPLIFY = FALSE)
  beta1   <- as.numeric(gamma.1[-(1:(q.1 + p.1))])
  
  Wtgam.2 <- mapply(function(w) as.numeric(w %*% gamma.2[1:q.2]), w = W.2, SIMPLIFY = FALSE)
  expZ.2  <- mapply(function(z) exp(as.numeric(z %*% gamma.2[(q.2+1):(q.2+p.2)])), z = Z.2, SIMPLIFY = FALSE)
  beta2   <- as.numeric(gamma.2[-(1:(q.2 + p.2))])
  
  haz.1.new <- update_haz_breslow_Zsur(
    tj = tj.1,
    Sur1.1 = Sur1.1,
    Wtgam = Wtgam.1,
    expZ  = expZ.1,
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
    expZ  = expZ.2,
    time_sur_list = s.2$time_sur_list_2,
    Zsur_list = s.2$Z.sur.2,
    b_y = bi.y,
    pb_y = pb.yt,
    beta_shared = beta2
  )
  
  ## rho update (Clayton copula) — general Zsur
  log_rho <- function(rho_tmp) {
    Q_rho_Zsur(
      rho = rho_tmp,
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
  theta.new$alpha   <- alpha.new
  theta.new$gamma.1 <- gamma.1.new
  theta.new$gamma.2 <- gamma.2.new
  theta.new$haz.1   <- haz.1.new
  theta.new$haz.2   <- haz.2.new
  theta.new$rho     <- rho.new
  
  if (!is.null(theta$Sigma.b)) theta.new$Sigma.b <- Sigma.b.new
  if (!is.null(theta$D2))      theta.new$D2      <- Sigma.b.new
  
  if (verbose) {
    cat("stepEM_Clayton_Poisson: rho=", rho.new, "\n")
  }
  
  list(theta = theta.new)
}