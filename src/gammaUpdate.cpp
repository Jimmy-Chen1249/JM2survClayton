#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List gammaUpdate_fixDelta_general(const Rcpp::List& b_,   // list: nMC x r
                                 const Rcpp::List& z_,   // list: r x (nj*K)  (or transposed)
                                 const Rcpp::List& w_,   // list: nMC x (nj*K)
                                 const Rcpp::List& pb_,  // list: nMC weights
                                 const arma::vec& haz,   // length >= nj
                                 const Rcpp::List& v_,   // list: q vector
                                 const Rcpp::List& u_,   // list: nj x qq (or compatible)
                                 const Rcpp::List& h_,   // list DF with delta, tj.ind
                                 const int& K,
                                 const int& q,
                                 const int& qq,
                                 const int& nev,
                                 const arma::vec& jcount) {

  const int P = q + qq + K;
  const int nsubj = w_.size();

  arma::mat Si      = arma::zeros<arma::mat>(P, nsubj); // integral score part (no delta)
  arma::mat Evstari = arma::zeros<arma::mat>(P, nsubj); // event part (has delta)
  arma::vec S       = arma::zeros<arma::vec>(P);

  arma::mat I      = arma::zeros<arma::mat>(P, P);      // info (no delta)
  arma::mat Gammaj = arma::zeros<arma::mat>(P, nev);    // correction building blocks (no delta)
  arma::mat Gamma  = arma::zeros<arma::mat>(P, P);

  auto coerce_u_to_mat = [&](SEXP uu, int nj) -> arma::mat {
    // Goal: return nj x qq
    if (Rf_isMatrix(uu)) {
      arma::mat um = Rcpp::as<arma::mat>(uu);

      // common bad case: 1 x nj -> transpose
      if ((int)um.n_rows == 1 && (int)um.n_cols == nj) um = um.t();

      // common bad case: qq x nj -> transpose
      if ((int)um.n_rows == qq && (int)um.n_cols == nj) um = um.t();

      // if qq==1 and um is nj x something, take first col
      if (qq == 1 && (int)um.n_rows == nj && (int)um.n_cols != 1) {
        um = um.col(0);
      }

      // now must be nj x qq
      if ((int)um.n_rows != nj) stop("u_ has incompatible shape: expected nj rows.");
      if ((int)um.n_cols != qq) stop("u_ has incompatible shape: expected qq columns.");
      return um;
    } else {
      // vector: treat as nj x 1 (qq must be 1)
      arma::vec uv = Rcpp::as<arma::vec>(uu);
      if ((int)uv.n_elem != nj) stop("u_ vector has incompatible length: expected nj.");
      if (qq != 1) stop("u_ provided as vector, but qq != 1. Provide u as nj x qq matrix.");
      return arma::mat(uv); // nj x 1
    }
  };

  // loop subjects
  for (int i = 0; i < nsubj; i++) {

    Rcpp::checkUserInterrupt();

    arma::mat b  = Rcpp::as<arma::mat>(b_[i]);   // nMC x r
    arma::mat z  = Rcpp::as<arma::mat>(z_[i]);   // r x (nj*K) expected (or transposed)
    arma::mat w  = Rcpp::as<arma::mat>(w_[i]);   // nMC x (nj*K)
    arma::vec pb = Rcpp::as<arma::vec>(pb_[i]);  // nMC
    arma::vec v  = Rcpp::as<arma::vec>(v_[i]);   // q
    arma::mat u  = coerce_u_to_mat(u_[i], (int)(w.n_cols / K));

    Rcpp::DataFrame h = Rcpp::as<Rcpp::DataFrame>(h_[i]);
    int tj_ind = h["tj.ind"];
    if (tj_ind == 0) continue;

    int delta = h["delta"];

    // infer nj
    int ncols_w = (int)w.n_cols;
    if (K <= 0) stop("K must be positive.");
    if (ncols_w % K != 0) stop("w_ n_cols must be divisible by K (expect nj*K columns).");
    int nj = ncols_w / K;
    if ((int)haz.n_elem < nj) stop("haz length < nj.");

    // ensure u has nj rows
    if ((int)u.n_rows != nj) stop("u_ rows != nj after coercion.");

    // infer r and fix z orientation if needed
    int r = (int)b.n_cols;
    if ((int)z.n_rows != r && (int)z.n_cols == r) {
      z = z.t();
    }
    if ((int)z.n_rows != r) stop("z_ has incompatible shape: expected r rows.");

    if ((int)z.n_cols != nj * K) stop("z_ must have nj*K columns.");

    // ---------- core ----------
    arma::mat Ii_int = arma::zeros<arma::mat>(P, P);

    arma::mat bzt   = b * z;       // nMC x (nj*K)
    arma::mat bztev = bzt % w;     // nMC x (nj*K)

    // Eexpvj(t): use FIRST block (cols 0..nj-1) as common baseline grid
    arma::rowvec Ew_all = arma::mean(w.each_col() % pb, 0);  // 1 x (nj*K)
    arma::rowvec Ew     = Ew_all.cols(0, nj - 1);            // 1 x nj

    arma::rowvec hazj   = arma::trans(haz.subvec(0, nj - 1));  // 1 x nj
    arma::rowvec Eexpvj = Ew % hazj;                            // 1 x nj

    double Eexpv = arma::as_scalar(arma::sum(Eexpvj));

    // uEexpv: length qq; u2Eexpv: qq x qq
    arma::rowvec uEexpv(qq, arma::fill::zeros);
    arma::mat    u2Eexpv(qq, qq, arma::fill::zeros);

    for (int t = 0; t < nj; t++) {
      arma::rowvec ut = u.row(t);         // 1 x qq
      double et = Eexpvj(t);
      uEexpv += et * ut;
      u2Eexpv += et * (ut.t() * ut);      // qq x qq
    }

    // ---------- info baseline blocks ----------
    // v-v
    if (q > 0) {
      Ii_int.submat(0, 0, q - 1, q - 1) = (v * v.t()) * Eexpv;
    }

    // v-u and u-v
    if (q > 0 && qq > 0) {
      for (int j = 0; j < q; j++) {
        for (int jj = 0; jj < qq; jj++) {
          Ii_int(j, q + jj) = v(j) * uEexpv(jj);
          Ii_int(q + jj, j) = Ii_int(j, q + jj);
        }
      }
    }

    // u-u
    if (qq > 0) {
      Ii_int.submat(q, q, q + qq - 1, q + qq - 1) = u2Eexpv;
    }

    // ---------- K blocks ----------
    // outj: 1 x (nj*K) = mean(bztev*pb) % (haz repeated K times)
    arma::rowvec haz_rep = arma::repmat(hazj, 1, K); // 1 x (nj*K)
    arma::rowvec Ebztev  = arma::mean(bztev.each_col() % pb, 0);
    arma::rowvec outj    = Ebztev % haz_rep;

    // info diagonal for K blocks: mean((bzt % bztev)*pb) % haz_rep
    arma::mat bzt2ev = (bzt % bztev); // nMC x (nj*K)
    arma::rowvec Ebzt2ev = arma::mean(bzt2ev.each_col() % pb, 0);
    arma::rowvec IiKdiag_all = Ebzt2ev % haz_rep;

    arma::rowvec Eb = arma::mean(b.each_col() % pb, 0); // 1 x r

    for (int k = 0; k < K; k++) {

      int c1 = nj * k;
      int c2 = nj * (k + 1) - 1;

      // event part: delta * (Zb)(T_i)  (take last grid col for this block)
      Evstari(q + qq + k, i) =
        delta * arma::as_scalar(z.col(c2).t() * Eb.t());

      // integral score: sum_t outj(block)
      Si(q + qq + k, i) = arma::sum(outj.subvec(c1, c2));

      // info diag
      Ii_int(q + qq + k, q + qq + k) =
        arma::sum(IiKdiag_all.subvec(c1, c2));

      // off-diag between k and k2
      for (int k2 = k + 1; k2 < K; k2++) {
        int d1 = nj * k2;
        int d2 = nj * (k2 + 1) - 1;

        arma::mat cross = (bztev.cols(c1, c2) % bzt.cols(d1, d2)); // nMC x nj
        arma::rowvec Ecross = arma::mean(cross.each_col() % pb, 0); // 1 x nj
        double val = arma::as_scalar(arma::sum(Ecross % hazj));

        Ii_int(q + qq + k, q + qq + k2) = val;
        Ii_int(q + qq + k2, q + qq + k) = val;
      }

      // cross with baseline v
      if (q > 0) {
        for (int j = 0; j < q; j++) {
          Ii_int(j, q + qq + k) = Si(q + qq + k, i) * v(j);
          Ii_int(q + qq + k, j) = Ii_int(j, q + qq + k);
        }
      }

      // cross with u (qq dims): sum_t outj(t)*u(t,jj)
      if (qq > 0) {
        arma::rowvec out_block = outj.subvec(c1, c2); // 1 x nj
        for (int jj = 0; jj < qq; jj++) {
          arma::rowvec ucol = u.col(jj).t(); // 1 x nj
          double uk = arma::as_scalar(arma::sum(out_block % ucol));
          Ii_int(q + jj, q + qq + k) = uk;
          Ii_int(q + qq + k, q + jj) = uk;
        }
      }

      // Gammaj rows for k (NO delta): for each t add outj(c1+t)
      for (int t = 0; t < nj; t++) {
        if (t < nev) {
          Gammaj(q + qq + k, t) += outj(c1 + t);
        }
      }
    }

    // ---------- baseline parts: Evstari/Si/Gammaj ----------
    if (q > 0) {
      Evstari.submat(0, i, q - 1, i) = delta * v;
      Si.submat(0, i, q - 1, i)      = v * Eexpv;

      for (int t = 0; t < nj; t++) {
        if (t < nev) {
          Gammaj.submat(0, t, q - 1, t) += v * Eexpvj(t);
        }
      }
    }

    if (qq > 0) {
      // event part: delta * u(T_i)  (take last grid row)
      Evstari.submat(q, i, q + qq - 1, i) = delta * u.row(nj - 1).t();

      // integral score: uEexpv
      Si.submat(q, i, q + qq - 1, i) = uEexpv.t();

      for (int t = 0; t < nj; t++) {
        if (t < nev) {
          Gammaj.submat(q, t, q + qq - 1, t) += u.row(t).t() * Eexpvj(t);
        }
      }
    }

    // subject contribution (delta only in event part)
    S += (Evstari.col(i) - Si.col(i));
    I += Ii_int;
  }

  // build Gamma correction
  for (int t = 0; t < nev; t++) {
    if (t < (int)jcount.n_elem && jcount(t) > 0) {
      Gamma += (Gammaj.col(t) * Gammaj.col(t).t()) / jcount(t);
    }
  }

  arma::vec gDelta;
  arma::mat A = I - Gamma;

  // solve with fallback
  bool ok = arma::solve(gDelta, A, S, arma::solve_opts::fast + arma::solve_opts::no_approx);
  if (!ok) {
    // fallback: ridge + approx
    arma::mat Ar = A + 1e-8 * arma::eye(P, P);
    ok = arma::solve(gDelta, Ar, S, arma::solve_opts::fast);
    if (!ok) {
      // last resort
      gDelta = arma::pinv(Ar) * S;
    }
  }

  return List::create(
    Named("gDelta") = gDelta,
    Named("I")      = I,
    Named("Gamma")  = Gamma,
    Named("S")      = S
  );
}