## =========================================================
## Unified Clayton joint-data generator (1D/2D RE, Gaussian/Poisson longitudinal)
## Integrates: Gdata_Clay, Gdata_Clay_P, Gdata_Clay_Two, Gdata_Clay_Two_P
## =========================================================
Gen_Joint_Clayton <- function(
    n = 200,
    seed = 1,
    MR = NULL,                 # target marginal censoring rate (vector length 2 or scalar)
    CR = c(4, 6),              # censoring upper bounds if MR is NULL
    rho0 = 0.5,                # Clayton parameter
    re_dim = 2,                # 1 or 2
    y_type = c("gaussian", "poisson"),
    ## ---- true parameters (defaults = what you gave) ----
    k1 = 0.1, k2 = 5,
    alpha0 = matrix(c(0.5, -1, 1, 1), nrow = 4),         # (Intercept, time, X1, X2)
    sigma0e = 0.5,
    Sigma_b0 = matrix(c(0.5, 0.2, 0.2, 0.5), 2, 2),      # for re_dim=2
    D0 = 1,                                              # for re_dim=1
    gamma011 = 0.3, gamma012 = -0.3, beta01 = 0.4,
    gamma021 = 0.4, gamma022 =  0.3, beta02 = 0.3,
    ## ---- longitudinal visit grid ----
    t_max = 2, mgrid = 20
) {
  
  y_type <- match.arg(y_type)
  
  ## ========= helper: one draw given CR =========
  gen_once <- function(CR_use) {
    
    set.seed(seed)
    
    ## -----------------------------
    ## 1) Clayton copula uniforms
    ## -----------------------------
    U1 <- runif(n)
    V1 <- runif(n)
    U2 <- (U1^(-rho0) * V1^(-rho0/(rho0 + 1)) - U1^(-rho0) + 1)^(-1/rho0)
    Tstar <- cbind(U1, U2)  # both in (0,1)
    
    ## -----------------------------
    ## 2) baseline cov + random effects
    ## -----------------------------
    W <- as.matrix(runif(n, 0, 1))
    
    if (re_dim == 1) {
      b <- rnorm(n, 0, D0)   # n x 1 (vector)
    } else if (re_dim == 2) {
      b <- MASS::mvrnorm(n, rep(0, ncol(Sigma_b0)), Sigma_b0)  # n x 2
    } else {
      stop("re_dim must be 1 or 2.")
    }
    
    ## piecewise covariate B(t):  Bi1 I(t<=Bi3) + Bi2 I(t>Bi3)
    B11 <- sample(c(-1, 1), n, replace = TRUE)
    B12 <- sample(c(-1, 1), n, replace = TRUE)
    B13 <- runif(n, 0, 1.75)
    
    B21 <- sample(c(-1, 1), n, replace = TRUE)
    B22 <- sample(c(-1, 1), n, replace = TRUE)
    B23 <- runif(n, 0, 1.75)
    
    ## -----------------------------
    ## 3) generate T1, T2 (piecewise closed form)
    ##    NOTE: matches your four files' algebra exactly (up to notation)
    ## -----------------------------
    if (re_dim == 1) {
      
      ## Use -log(1 - pnorm(Tstar)) as in your 1D code (Clay / Clay_P)
      E1 <- -log(1 - pnorm(Tstar)[, 1])
      E2 <- -log(1 - pnorm(Tstar)[, 2])
      
      ## T1
      Temp1 <- log(1 + k2 * E1 / (k1 * exp(W %*% gamma011 + beta01 * b + gamma012 * B11))) / k2
      Temp1[is.na(Temp1)] <- 0
      T1 <- Temp1 * (Temp1 <= B13)
      
      Temp2 <- log((k2 * E1 / (k1 * exp(W %*% gamma011 + beta01 * b)) -
                      exp(k2 * B13) * (exp(B11 * gamma012) - exp(B12 * gamma012)) +
                      exp(B11 * gamma012)) / exp(B12 * gamma012)) / k2
      Temp2[is.na(Temp2)] <- 0
      T1 <- T1 + Temp2 * (Temp2 > B13)
      
      ## T2
      Temp1 <- log(1 + k2 * E2 / (k1 * exp(W %*% gamma021 + beta02 * b + gamma022 * B21))) / k2
      Temp1[is.na(Temp1)] <- 0
      T2 <- Temp1 * (Temp1 <= B23)
      
      Temp2 <- log((k2 * E2 / (k1 * exp(W %*% gamma021 + beta02 * b)) -
                      exp(k2 * B23) * (exp(B21 * gamma022) - exp(B22 * gamma022)) +
                      exp(B21 * gamma022)) / exp(B22 * gamma022)) / k2
      Temp2[is.na(Temp2)] <- 0
      T2 <- T2 + Temp2 * (Temp2 > B23)
      
      T1[is.nan(T1)] <- 10
      T2[is.nan(T2)] <- 10
      
      T_orig <- cbind(T1, T2)
      
    } else {
      ## re_dim == 2: matches your Gdata_Clay_Two.R version
      ## NOTE: you used -log(1 - Tstar)[,k] in the latest “Two” code you pasted
      ## (i.e. treating Tstar itself as uniform). We keep that to match your current code.
      E1 <- -log(1 - Tstar[, 1])
      E2 <- -log(1 - Tstar[, 2])
      
      B11.t <- exp(W %*% gamma011 + beta01 * b[, 1] + gamma012 * B11)
      B12.t <- exp(W %*% gamma011 + beta01 * b[, 1] + gamma012 * B12)
      
      den1 <- (k2 + b[, 2] * beta01)
      Temp1 <- log(1 + den1 * E1 / (k1 * B11.t)) / den1
      Temp1[is.na(Temp1)] <- 0
      T1 <- Temp1 * (Temp1 <= B13)
      
      Temp2 <- log((den1 * E1 + k1 * B11.t +
                      exp(k2 * B13 + b[, 2] * beta01 * B13) * (k1 * B12.t - k1 * B11.t)) /
                     (k1 * B12.t)) / den1
      Temp2[is.na(Temp2)] <- 0
      T1 <- T1 + Temp2 * (Temp2 > B13)
      
      B21.t <- exp(W %*% gamma021 + beta02 * b[, 1] + gamma022 * B21)
      B22.t <- exp(W %*% gamma021 + beta02 * b[, 1] + gamma022 * B22)
      
      den2 <- (k2 + b[, 2] * beta02)
      Temp1 <- log(1 + den2 * E2 / (k1 * B21.t)) / den2
      Temp1[is.na(Temp1)] <- 0
      T2 <- Temp1 * (Temp1 <= B23)
      
      Temp2 <- log((den2 * E2 + k1 * B21.t +
                      exp(k2 * B23 + b[, 2] * beta02 * B23) * (k1 * B22.t - k1 * B21.t)) /
                     (k1 * B22.t)) / den2
      Temp2[is.na(Temp2)] <- 0
      T2 <- T2 + Temp2 * (Temp2 > B23)
      
      T1[is.nan(T1)] <- mean(T1, na.rm = TRUE)
      T2[is.nan(T2)] <- mean(T2, na.rm = TRUE)
      
      T_orig <- cbind(T1, T2)
    }
    
    ## -----------------------------
    ## 4) censoring
    ## -----------------------------
    C1 <- runif(n, 0, CR_use[1])
    C2 <- runif(n, 0, CR_use[2])
    delta <- (T_orig <= cbind(C1, C2)) * 1
    
    T_cen <- cbind(pmin(T_orig[, 1], C1), pmin(T_orig[, 2], C2))
    T_max <- apply(T_cen, 1, max)
    
    ## -----------------------------
    ## 5) longitudinal grid + X + y
    ## -----------------------------
    t_cond <- seq(0, t_max, length.out = mgrid)
    t_list <- lapply(1:n, function(i) t_cond[t_cond < T_max[i]])
    ni <- vapply(t_list, length, integer(1))
    
    X3 <- rbinom(n, 1, 0.5)
    X4 <- runif(n, 0, 2)
    
    X <- cbind(
      1,
      unlist(t_list),
      rep(X3, ni),
      rep(X4, ni)
    )
    
    if (y_type == "gaussian") {
      e <- unlist(lapply(1:n, function(i) rnorm(ni[i], 0, sigma0e)))
      
      if (re_dim == 1) {
        y <- as.numeric(X %*% alpha0 + rep(b, ni) + e)
      } else {
        y <- as.numeric(X %*% alpha0 + rep(b[, 1], ni) + unlist(t_list) * rep(b[, 2], ni) + e)
      }
      
    } else { # poisson
      if (re_dim == 1) {
        mu <- as.numeric(exp(X %*% alpha0 + rep(b, ni)))
      } else {
        mu <- as.numeric(exp(X %*% alpha0 + rep(b[, 1], ni) + unlist(t_list) * rep(b[, 2], ni)))
      }
      y <- rpois(length(mu), mu)
    }
    
    ## -----------------------------
    ## 6) pack
    ## -----------------------------
    obse.no <- unlist(lapply(1:n, function(i) seq_len(ni[i])))
    id <- rep(seq_len(n), ni)
    
    out <- list(
      id = id,
      ni = ni,
      obse.no = obse.no,
      X = X,
      y = y,
      T.cen = T_cen,
      W = W,
      delta = delta,
      X_al = cbind(X3, X4),
      t = t_list,
      B1 = cbind(B11, B12, B13),
      B2 = cbind(B21, B22, B23)
    )
    
    if (re_dim == 1) {
      out$b.long <- rep(b, ni)
      out$b.surv <- b
    } else {
      out$b.long1 <- rep(b[, 1], ni)
      out$b.long2 <- rep(b[, 2], ni)
      out$b.surv  <- b
    }
    
    out
  }
  
  ## ========= if MR is provided, pick CR by grid search (like Cen_Clay*) =========
  if (!is.null(MR)) {
    if (length(MR) == 1) MR <- rep(MR, 2)
    
    CR_grid <- seq(1, 3, length.out = 20)
    CRR <- cbind(rep(CR_grid, each = length(CR_grid)),
                 rep(CR_grid, length(CR_grid)))
    
    mr <- matrix(0, nrow(CRR), 2)
    for (i in seq_len(nrow(CRR))) {
      dat_tmp <- gen_once(CRR[i, ])
      mr[i, ] <- 1 - colMeans(dat_tmp$delta)
    }
    i_opt <- which.min(rowMeans(abs(mr - matrix(MR, nrow(CRR), 2, byrow = TRUE))))
    return(gen_once(CRR[i_opt, ]))
  } else {
    return(gen_once(CR))
  }
}