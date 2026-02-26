#' Baseline hazard update and rho objective (Clayton survival copula)
#'
#' This file contains two core utilities used by the EM algorithm:
#' 1) Breslow-type baseline hazard increments update with posterior weights.
#' 2) Q-function in rho for the Clayton survival copula, using posterior weights.
#'
#' Both functions are fully general w.r.t. the random-effect design Z_sur(t):
#' - re_dim = 1: Z_sur(t) is (1)
#' - re_dim = 2: Z_sur(t) can be (1, t)
#' - etc.
#'
#' @keywords internal
#' @noRd
NULL

#' Breslow update for baseline hazard increments with general Z_sur(t)
#'
#' @param tj Global event time grid (vector).
#' @param Sur1.1 One-row-per-subject data.frame: id, T, delta, tj.ind, ...
#' @param Wtgam List: exp part W_i^T gamma_{k1} is provided as numeric scalar per id.
#' @param expZ List: exp(Z^*(t)^T gamma_{k2}) at subject grid times.
#' @param time_sur_list List: subject grid times (increasing).
#' @param Zsur_list List: subject Z_sur(t) matrices (n_i x r).
#' @param b_y List: MC draws of b, each is nMC x r.
#' @param pb_y List: posterior weights, each is length nMC (mean=1 recommended).
#' @param beta_shared Scalar association parameter beta_k (for K=1 case).
#' @param eps Numerical guard for denominator.
#'
#' @return haz vector, length length(tj).
#' @keywords internal
#' @noRd
update_haz_breslow_Zsur <- function(
    tj,
    Sur1.1,
    Wtgam,
    expZ,
    time_sur_list,
    Zsur_list,
    b_y,
    pb_y,
    beta_shared,
    eps = 1e-12
) {
  J <- length(tj)
  haz <- numeric(J)
  
  ## fast & stable: idx = max{m: time_i[m] <= uu}
  find_idx <- function(time_i, uu) {
    idx <- findInterval(uu, time_i)  # 0..length(time_i)
    if (idx < 1) idx <- 1
    idx
  }
  
  for (j in seq_len(J)) {
    uu <- tj[j]
    
    ## numerator: #events at uu
    dNj <- sum(Sur1.1$delta == 1 & Sur1.1$T == uu)
    
    ## risk set: T_i >= uu
    risk <- which(Sur1.1$T >= uu)
    if (length(risk) == 0) {
      haz[j] <- 0
      next
    }
    
    denom <- 0
    
    for (ii in risk) {
      id_i <- as.character(Sur1.1$id[ii])
      
      time_i <- time_sur_list[[id_i]]
      if (is.null(time_i) || length(time_i) == 0) next
      idx_i <- find_idx(time_i, uu)
      
      ## scalar exp(W^T gamma) * exp(Z^*(uu)^T gamma)
      sc_ij <- exp(Wtgam[[id_i]]) * expZ[[id_i]][idx_i]
      
      ## MC expectation of exp(beta * Zsur(uu) %*% b)
      bmat <- as.matrix(b_y[[id_i]])          # nMC x r
      pb   <- as.numeric(pb_y[[id_i]])        # length nMC
      Zsur_i <- as.matrix(Zsur_list[[id_i]])  # n_i x r
      
      zrow <- Zsur_i[idx_i, , drop = FALSE]   # 1 x r
      lin  <- as.numeric(bmat %*% t(zrow))    # nMC
      
      E_ij <- mean(exp(beta_shared * lin) * pb)
      
      denom <- denom + sc_ij * E_ij
    }
    
    denom <- max(denom, eps)
    haz[j] <- dNj / denom
  }
  
  haz
}

#' Q-function in rho for Clayton survival copula with general Z_sur(t)
#'
#' This returns the NEGATIVE expected complete-data log-likelihood (w.r.t rho),
#' so it can be passed directly to `optimize()` (which minimizes).
#'
#' @param rho Clayton parameter (>0).
#' @param Sur1.1,Sur2.1 One-row-per-subject survival summaries for the two endpoints.
#' @param tj.1,tj.2 Global event time grids.
#' @param haz.1,haz.2 Baseline hazard increments (length >= length(tj)).
#' @param Wtgam.1,Wtgam.2 Lists: W^T gamma_{k1} per subject.
#' @param expZ.1,expZ.2 Lists: exp(Z^*(t)^T gamma_{k2}) per subject grid.
#' @param time_sur_list_1,time_sur_list_2 Lists: subject grid times.
#' @param Zsur_list_1,Zsur_list_2 Lists: subject Z_sur(t) matrices.
#' @param b_y,pb_y Lists: MC b draws and posterior weights.
#' @param beta1,beta2 Scalars association parameters (K=1 case).
#'
#' @return scalar (negative Q).
#' @keywords internal
#' @noRd
Q_rho_Zsur <- function(
    rho,
    Sur1.1, Sur2.1,
    tj.1, tj.2,
    haz.1, haz.2,
    Wtgam.1, Wtgam.2,
    expZ.1, expZ.2,
    time_sur_list_1, time_sur_list_2,
    Zsur_list_1, Zsur_list_2,
    b_y, pb_y,
    beta1, beta2
) {
  n <- nrow(Sur1.1)
  
  ## Clayton survival copula density c(u,v) with u=S1, v=S2
  clayton_c <- function(u, v, rho) {
    (1 + rho) * (u * v)^(-(1 + rho)) * (u^(-rho) + v^(-rho) - 1)^(-(2 + 1/rho))
  }
  
  out_i <- numeric(n)
  
  for (i in seq_len(n)) {
    
    T1 <- Sur1.1$T[i]; d1 <- Sur1.1$delta[i]
    T2 <- Sur2.1$T[i]; d2 <- Sur2.1$delta[i]
    
    id_i <- as.character(Sur1.1$id[i])
    
    t1_i <- time_sur_list_1[[id_i]]
    t2_i <- time_sur_list_2[[id_i]]
    
    ## find the last grid index <= T
    j1 <- max(which(t1_i <= T1)); if (!is.finite(j1)) j1 <- 1L
    j2 <- max(which(t2_i <= T2)); if (!is.finite(j2)) j2 <- 1L
    
    eta1 <- exp(Wtgam.1[[id_i]])
    eta2 <- exp(Wtgam.2[[id_i]])
    
    ez1 <- expZ.1[[id_i]]
    ez2 <- expZ.2[[id_i]]
    
    bmat <- as.matrix(b_y[[id_i]])       # nMC x r
    pb   <- as.numeric(pb_y[[id_i]])     # length nMC
    
    ## Zsur rows at grid
    Zsur1 <- as.matrix(Zsur_list_1[[id_i]])  # n1 x r
    Zsur2 <- as.matrix(Zsur_list_2[[id_i]])  # n2 x r
    
    ## cumulative hazards Hk(Tk|b)
    lin1 <- beta1 * as.numeric(bmat %*% t(Zsur1[1:j1, , drop = FALSE]))  # nMC * j1 (vectorized)
    lin1 <- matrix(lin1, nrow = nrow(bmat), byrow = FALSE)
    
    lin2 <- beta2 * as.numeric(bmat %*% t(Zsur2[1:j2, , drop = FALSE]))
    lin2 <- matrix(lin2, nrow = nrow(bmat), byrow = FALSE)
    
    H1m <- as.numeric(exp(lin1) %*% (haz.1[1:j1] * (eta1 * ez1[1:j1])))
    H2m <- as.numeric(exp(lin2) %*% (haz.2[1:j2] * (eta2 * ez2[1:j2])))
    
    S1m <- exp(-H1m)
    S2m <- exp(-H2m)
    
    ## hazard at event time (grid j1/j2)
    lin1_T <- beta1 * as.numeric(bmat %*% t(Zsur1[j1, , drop = FALSE]))
    lin2_T <- beta2 * as.numeric(bmat %*% t(Zsur2[j2, , drop = FALSE]))
    
    lam1_T <- haz.1[j1] * eta1 * ez1[j1] * exp(lin1_T)
    lam2_T <- haz.2[j2] * eta2 * ez2[j2] * exp(lin2_T)
    
    ## Clayton survival copula C(u,v)
    Cuv <- (S1m^(-rho) + S2m^(-rho) - 1)^(-1/rho)
    
    ## partial derivatives
    dC_du <- Cuv^(rho + 1) * S1m^(-rho - 1)
    dC_dv <- Cuv^(rho + 1) * S2m^(-rho - 1)
    
    cden <- clayton_c(S1m, S2m, rho)
    
    loglik_m <- if (d1 == 1 && d2 == 1) {
      log(pmax(lam1_T, 1e-300)) + log(pmax(lam2_T, 1e-300)) + log(pmax(cden, 1e-300))
    } else if (d1 == 1 && d2 == 0) {
      log(pmax(lam1_T, 1e-300)) + log(pmax(dC_du, 1e-300))
    } else if (d1 == 0 && d2 == 1) {
      log(pmax(lam2_T, 1e-300)) + log(pmax(dC_dv, 1e-300))
    } else {
      log(pmax(Cuv, 1e-300))
    }
    
    ## posterior-weighted expectation
    out_i[i] <- mean(loglik_m * pb)
  }
  
  ## optimize() minimizes -> return negative
  -sum(out_i)
}