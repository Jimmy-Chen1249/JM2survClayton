// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// =============================== 工具函数（稳健版） ===============================
  
  // 多元正态 N(0, Sigma) 的 log 密度（数值稳健：chol + solve），ok = false 表示 Sigma 非 SPD
static double ldmvnorm0(const arma::vec& x, const arma::mat& Sigma, bool& ok) {
  ok = true;
  if (Sigma.n_rows != Sigma.n_cols || Sigma.n_rows != (int)x.n_elem) { ok=false; return R_NegInf; }
  arma::mat L;
  bool chol_ok = arma::chol(L, Sigma, "lower");
  if (!chol_ok) { ok=false; return R_NegInf; }
  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y = arma::solve(arma::trimatl(L), x);  // L y = x
  const double quad  = arma::dot(y, y);            // x^T Sigma^{-1} x
  const double k = (double)x.n_elem;
  return -0.5 * (k*std::log(2.0*M_PI) + logdet + quad);
}

// Poisson 对数似然（向量化）：ll = y^T eta - sum(exp(eta))
// eta = Xbeta + Z u；对 eta 截断以防 exp 溢出
static double ll_pois(const arma::vec& Xbeta, const arma::mat& Z,
                      const arma::vec& u, const arma::vec& y) {
  arma::vec eta = Xbeta + Z * u;
  // double 精度下 exp(>709) 溢出，留余量
  eta = arma::clamp(eta, -700.0, 700.0);
  return arma::dot(y, eta) - arma::sum(arma::exp(eta));
}

// ---------------------- 老接口需要的 min0（保留以兼容） ----------------------
  // [[Rcpp::export]]
double min0(double a, double b) { return (a < b) ? a : b; }

// =============================== 对数似然（新旧接口） ===============================
  
  // 新：已预先计算 Xbeta 的版本（推荐）
// [[Rcpp::export]]
double loglikelihoodPoisson_preX(const arma::vec& Xbeta, const arma::mat& Z,
                                 const arma::vec& u,    const arma::vec& y,
                                 const arma::mat& Sigma) {
  bool ok = true;
  double val = ll_pois(Xbeta, Z, u, y) + ldmvnorm0(u, Sigma, ok);
  return ok ? val : R_NegInf;
}

// 旧：仍接受 beta 与 X，内部计算 Xbeta，保持向后兼容
// [[Rcpp::export]]
double loglikelihoodPoissonCpp_n(const arma::vec& beta,  const arma::mat& Sigma,
                                 const arma::vec& u,     const arma::vec& y,
                                 const arma::mat& X,     const arma::mat& Z) {
  arma::vec Xbeta = X * beta;
  return loglikelihoodPoisson_preX(Xbeta, Z, u, y, Sigma);
}

// 接受概率（旧接口，保留以兼容；内部走新实现）
// [[Rcpp::export]]
double logAcceptPoisson_n(const arma::vec& beta, const arma::mat& Sigma,
                          const arma::vec& ucur, const arma::vec& uprop,
                          const arma::vec& y,    const arma::mat& X,
                          const arma::mat& Z) {
  arma::vec Xbeta = X * beta;
  const double ll_c = loglikelihoodPoisson_preX(Xbeta, Z, ucur, y, Sigma);
  const double ll_p = loglikelihoodPoisson_preX(Xbeta, Z, uprop, y, Sigma);
  return min0(0.0, ll_p - ll_c);
}

// =============================== 采样器（相关性随机游走） ===============================
  
  // 新：返回采样矩阵 + 接受率（用于自适应调参）
// 提议：u' = u + sd0 * L z，L 为 chol(Sigma)
// [[Rcpp::export]]
List uSamplerPoissonCpp_n_acc(const arma::vec& Xbeta, const arma::mat& Z,
                              const arma::vec& u0,    const arma::vec& y,
                              const arma::mat& Sigma, int B, double sd0) {
  Rcpp::RNGScope scope;

  const int k = u0.n_rows;
  arma::mat L;
  if (!arma::chol(L, Sigma, "lower")) {
    stop("Sigma is not SPD in proposal (chol failed).");
  }

  arma::mat usample(B, k, arma::fill::zeros);
  arma::vec ucur = u0, uprop(k, arma::fill::zeros);
  usample.row(0) = ucur.t();

  double ll_cur = loglikelihoodPoisson_preX(Xbeta, Z, ucur, y, Sigma);
  int acc = 0;

  for (int i = 1; i < B; ++i) {
    arma::vec z = arma::randn<arma::vec>(k);
    uprop = ucur + sd0 * (L * z);
    double ll_prop = loglikelihoodPoisson_preX(Xbeta, Z, uprop, y, Sigma);

    const double loga = ll_prop - ll_cur; // 对称提议
    if (std::log(R::runif(0.0, 1.0)) < loga) {
      ucur  = uprop;
      ll_cur = ll_prop;
      acc++;
    }
    usample.row(i) = ucur.t();
  }

  const double acc_rate = (B > 1) ? (double)acc / (double)(B - 1) : NA_REAL;
  return List::create(_["sample"] = usample, _["acc_rate"] = acc_rate);
}

// 旧：保持原签名与返回值（仅返回采样矩阵），内部同样使用相关性提议
// [[Rcpp::export]]
arma::mat uSamplerPoissonCpp_n(const arma::vec& beta,  const arma::mat& Sigma,
                               const arma::vec& u0,    const arma::vec& y,
                               const arma::mat& X,     const arma::mat& Z,
                               int B, double sd0) {
  Rcpp::RNGScope scope;

  const int k = u0.n_rows;
  arma::mat L;
  if (!arma::chol(L, Sigma, "lower")) {
    stop("Sigma is not SPD in proposal (chol failed).");
  }

  arma::vec Xbeta = X * beta;

  arma::mat usample(B, k, arma::fill::zeros);
  arma::vec ucur = u0, uprop(k, arma::fill::zeros);
  usample.row(0) = ucur.t();

  double ll_cur = loglikelihoodPoisson_preX(Xbeta, Z, ucur, y, Sigma);

  for (int i = 1; i < B; ++i) {
    arma::vec z = arma::randn<arma::vec>(k);
    uprop = ucur + sd0 * (L * z);
    double ll_prop = loglikelihoodPoisson_preX(Xbeta, Z, uprop, y, Sigma);

    const double loga = ll_prop - ll_cur; // 对称提议
    if (std::log(R::runif(0.0, 1.0)) < loga) {
      ucur  = uprop;
      ll_cur = ll_prop;
    }
    usample.row(i) = ucur.t();
  }
  return usample;
}

// =============================== 可选：多元 t 的 log 密度（若你别处需要） ===============================
// [[Rcpp::export]]
double ldmt(const arma::vec& x, double df, const arma::mat& Sigma) {
  if (df <= 0.0) return R_NegInf;
  arma::mat L;
  if (!arma::chol(L, Sigma, "lower")) return R_NegInf;
  const int k = x.n_elem;
  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y = arma::solve(arma::trimatl(L), x);
  const double quad = arma::dot(y, y);
  return R::lgammafn(0.5 * (df + k)) - R::lgammafn(0.5 * df)
       - 0.5 * (k * std::log(df) + k * std::log(M_PI) + logdet)
       - 0.5 * (df + k) * std::log1p(quad / df);
}