// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

//' @title max_range_cpp
//'
//' @description Max range to sample rho using rcpp
//'
//' @param Ws Spatial adjacency matrix
//' @param Wt Temporal adjacency matrix
//' @param rho_s Spatial parameter
//' @param rho_t Temporal parameter
//' @param rho_st Spatio-temporal parameter
//'
// [[Rcpp::export]]
Rcpp::List max_range_cpp(arma::mat Ws, arma::mat Wt, double rho_s, double rho_t, double rho_st){

  int i, n, n_var;
  double max_s, max_t, max_st;

  n = Ws.n_rows;
  n_var = Wt.n_rows;

  arma::mat ind_n = arma::eye<arma::mat>(n, n);
  arma::mat ind_nvar = arma::eye<arma::mat>(n_var, n_var);

  arma::mat kron_s(n*n_var, n*n_var);
  arma::mat kron_t(n*n_var, n*n_var);
  arma::mat kron_st(n*n_var, n*n_var);
  kron_s.zeros();
  kron_t.zeros();
  kron_st.zeros();

  kron_s = arma::kron(Ws, ind_nvar);
  kron_t = arma::kron(ind_n, Wt);
  kron_st = arma::kron(Ws, Wt);

  arma::mat Q_max = kron_s + kron_t + kron_st;
  arma::mat D = arma::cumsum(Q_max, 0);

  arma::mat D_s = arma::cumsum(kron_s, 0);
  arma::mat D_t = arma::cumsum(kron_t, 0);
  arma::mat D_st = arma::cumsum(kron_st, 0);

  arma::mat D_s_aux = arma::cumsum(kron_s, 0);
  arma::mat D_t_aux = arma::cumsum(kron_t, 0);
  arma::mat D_st_aux = arma::cumsum(kron_st, 0);

  arma::vec D_vec(n*n_var), D_s_vec(n*n_var), D_t_vec(n*n_var), D_st_vec(n*n_var);

  for(i = 0; i < n*n_var; i++){
    D_vec(i) = D(n*n_var - 1, i);
    D_s_vec(i) = D_s_aux(n*n_var - 1, i);
    D_t_vec(i) = D_t_aux(n*n_var - 1, i);
    D_st_vec(i) = D_st_aux(n*n_var - 1, i);
  }

  arma::vec D_s_busy(n*n_var), D_t_busy(n*n_var), D_st_busy(n*n_var), D_tot_busy(n*n_var);

  D_s_busy = D_t_vec*::fabs(rho_t) + D_st_vec*::fabs(rho_st);
  D_t_busy = D_s_vec*::fabs(rho_s) + D_st_vec*::fabs(rho_st);
  D_st_busy = D_s_vec*::fabs(rho_s) + D_t_vec*::fabs(rho_t);
  D_tot_busy = D_s_vec*::fabs(rho_s) + D_t_vec*::fabs(rho_t) + D_st_vec*::fabs(rho_st);

  arma::vec max_s_vec(n*n_var), max_t_vec(n*n_var), max_st_vec(n*n_var), max_tot_vec(n*n_var);

  for(i = 0; i < n*n_var; i++){
    max_s_vec(i) = (D_vec(i) - D_s_busy(i))/D_s_vec(i);
    max_t_vec(i) = (D_vec(i) - D_t_busy(i))/D_t_vec(i);
    max_st_vec(i) = (D_vec(i) - D_st_busy(i))/D_st_vec(i);
    max_tot_vec(i) = (D_vec(i) - D_tot_busy(i));
  }

  max_s = min(max_s_vec);
  max_t = min(max_t_vec);
  max_st = min(max_st_vec);
  double tot = min(max_tot_vec);

  Rcpp::List params = Rcpp::List::create(Rcpp::Named("rho_s") = max_s,
                                         Rcpp::Named("rho_t") = max_t,
                                         Rcpp::Named("rho_st") = max_st,
                                         Rcpp::Named("tot") = tot);

  return params;
}

//' @title buildQC_cpp
//'
//' @description CAR covariance matrix
//'
//' @param M A vector containing the number of neighbors
//' @param W A matrix containing the neighborhood structure
//' @param rho The spatial dependence parameter
//'
// [[Rcpp::export]]
arma::mat buildQC_cpp(arma::vec M, arma::mat W, arma::vec rho){
  int N = M.n_elem;
  int i = 0, j = 0;

  arma::mat Q(N, N);
  Q.zeros();

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      Q(i, j) = (rho(0) * W(i, j));

  for(i = 0; i < N; i++)
    Q(i, i) = M(i);

  return Q;
}

//' @title buildQL_cpp
//'
//' @description Leroux covariance matrix
//'
//' @param R TO DO
//' @param rho Dependence parameter
//'
// [[Rcpp::export]]
arma::mat buildQL_cpp(arma::mat R, arma::vec rho){
  int N = R.n_rows;
  int i = 0, j = 0;

  arma::mat Q(N, N);
  Q.zeros();

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      Q(i, j) = ((1-rho(0)) + Q(i, j));

  for(i = 0; i < N; i++)
    Q(i, i) = (rho(0) * R(i));

  return Q;
}

//' @title buildQST_cpp
//'
//' @description Build a Spatio-temporal covariance matrix
//'
//' @param Ws Spatial neighborhod matrix
//' @param Wt Temporal neighborhod matrix
//' @param D A vector containing the number of neighbors on diagonal
//' @param rho_s Spatial dependence parameter
//' @param rho_t Temporal dependence parameter
//' @param rho_st Spatio-temporal dependence parameter
//'
// [[Rcpp::export]]
void buildQST_cpp(arma::mat& Q, arma::mat& Ws, arma::mat& Wt, double& rho_s, double& rho_t, double& rho_st){
  int i = 0, nReg = Ws.n_rows, nVar = Wt.n_rows;

  arma::mat Inreg = arma::eye<arma::mat>(nReg, nReg);
  arma::mat Invar = arma::eye<arma::mat>(nVar, nVar);

  arma::mat kronAux1(nVar*nReg, nVar*nReg);
  arma::mat kronAux2(nVar*nReg, nVar*nReg);
  arma::mat kronAux3(nVar*nReg, nVar*nReg);
  kronAux1.zeros();
  kronAux2.zeros();
  kronAux3.zeros();

  kronAux1 = arma::kron(Ws, Invar);
  kronAux2 = arma::kron(Inreg, Wt);
  kronAux3 = arma::kron(Ws, Wt);

  Q.zeros();
  Q = -(rho_s * kronAux1 + rho_t * kronAux2 + rho_st * kronAux3);

  arma::mat D = arma::cumsum(kronAux1 + kronAux2 + kronAux3, 0);

  for(i = 0; i < nReg*nVar; i++)
    Q(i, i) = D(nVar*nReg - 1, i);
}

//' @title rmvnorm_cpp
//'
//' @description Random values from a multivariate normal
//'
//' @param n number of samples
//' @param mean vector of means
//' @param sigma covariance matrix
//'
// [[Rcpp::export]]
arma::vec rmvnorm_cpp(arma::mat sigma){
  arma::vec result;
  int ncols = sigma.n_cols;

  arma::mat Y = arma::randn(1, ncols);
  Y = Y * chol(sigma);
  result = Y.row(0).t();

  return result;
}

//' @title rtruncnorm_cpp
//'
//' @description To sample from a multivariate truncated normal
//'
//' @param a Lower bounds vector
//' @param b Upper bounds vector
//' @param mean Vector of means
//' @param var Covariance matrix
//
// [[Rcpp::export]]
double rtruncnorm_cpp(double a, double b, double mean, double var){
  Rcpp::Environment truncnorm("package:truncnorm");
  Rcpp::Function rtruncnorm = truncnorm["rtruncnorm"];

  double res = Rcpp::as<double>(rtruncnorm(1L, a, b, mean, sqrt(var)));

  return res;
}

//' @title rtmvnorm_cpp
//'
//' @description To sample from a multivariate truncated normal
//'
//' @param a Lower bounds vector
//' @param b Upper bounds vector
//' @param mean Vector of means
//' @param var Covariance matrix
//
// [[Rcpp::export]]
arma::vec rtmvnorm_cpp(arma::vec a, arma::vec b, arma::vec mean, arma::mat var){
  int i = 0, j = 0;

  Rcpp::Environment tmvtnorm("package:tmvtnorm");
  Rcpp::Function rtmvnorm = tmvtnorm["rtmvnorm"];

  arma::mat D = arma::eye<arma::mat>(var.n_rows, var.n_cols);

  NumericVector a_nv = NumericVector(a.begin(), a.end());
  NumericVector b_nv = NumericVector(b.begin(), b.end());
  NumericVector mean_nv = NumericVector(mean.begin(), mean.end());
  NumericMatrix var_nv(var.n_rows, var.n_cols);

  for(i = 0; i < var.n_rows; i++) {
    for(j = 0; j < var.n_cols; j++) {
      var_nv(i, j) = var(i, j);
    }
  }

  arma::vec res = as<arma::vec>(rtmvnorm(1L, mean_nv, var_nv, a_nv, b_nv, D, R_NilValue, "gibbs"));

  return res;
}

//' @title dtruncnorm_cpp
//'
//' @description To get the density of a truncated normal
//'
//' @param a Lower bound vector
//' @param b Upper bound vector
//' @param mean Normal expectation vector
//' @param var Normal variance matrix
//' @param x Point to evaluate the density
//
// [[Rcpp::export]]
double dtruncnorm_cpp(double x, double a, double b, double mean, double var, int l){
  Rcpp::Environment truncnorm("package:truncnorm");
  Rcpp::Function dtruncnorm = truncnorm["dtruncnorm"];

  double res = Rcpp::as<double>(dtruncnorm(x, a, b, mean, sqrt(var)));

  if(l) res = log(res);

  return res;
}

//' @title dtmvnorm_cpp
//'
//' @description To get the density of a multivariate truncated normal
//'
//' @param a Lower bound vector
//' @param b Upper bound vector
//' @param mean Normal expectation vector
//' @param var Normal variance matrix
//' @param x Point to evaluate the density
//
// [[Rcpp::export]]
double dtmvnorm_cpp(arma::vec x, arma::vec a, arma::vec b, arma::vec mean, arma::mat var){
  int i = 0, j = 0;

  Rcpp::Environment tmvtnorm("package:tmvtnorm");
  Rcpp::Function dtmvnorm = tmvtnorm["dtmvnorm"];

  NumericVector x_nv = NumericVector(x.begin(), x.end());
  NumericVector a_nv = NumericVector(a.begin(), a.end());
  NumericVector b_nv = NumericVector(b.begin(), b.end());
  NumericVector mean_nv = NumericVector(mean.begin(), mean.end());
  NumericMatrix var_nv(var.n_rows, var.n_cols);

  for(i = 0; i < var.n_cols; i++) {
    for(j = 0; j < var.n_cols; j++) {
      var_nv(i, j) = var(i, j);
    }
  }

  double res = as<double>(dtmvnorm(x_nv, mean_nv, var_nv, a_nv, b_nv, 1L));

  return res;
}

//' @title dens_beta_cpp
//'
//' @description To usage in MCMC estimation
//'
// [[Rcpp::export]]
double dens_beta_cpp(arma::vec x_beta, arma::vec x_eps, double x_nu, arma::mat sigma,
                     Rcpp::List& params){
  int i;
  double out = 0.0, mu = 0.0;

  /// Arguments
  int N = params["N"];
  int P = params["P"];
  int type = params["type"];

  arma::vec y = params["y"];
  arma::vec e = params["e"];
  arma::mat X = params["X"];

  arma::vec tau_beta = params["tau_beta"];
  arma::vec mean_beta = params["mean_beta"];

  /// Auxiliar vector
  arma::vec Xbet(N);

  /// Updating Xbeta
  Xbet = X*x_beta;

  /// Variable
  for(i = 0; i < N; i++){
    if(type == 1) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), x_nu, exp(Xbet(i))/x_nu, 1, 1);
    if(type == 2) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*x_nu, 1/x_nu, 1, 1);
    if(type == 3) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*exp(Xbet(i))*x_nu, 1/(x_nu*exp(Xbet(i))), 1, 1);
    if(type == 4) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(sigma(i, i)/x_nu), 1, 1);
    if(type == 5) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(1/x_nu), 1, 1);
    if(type == 6) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), 1/x_nu, exp(Xbet(i)), 1, 1);
    if(type == 7) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i)), 1/x_nu, 1, 1);

    out = out - e(i)*mu + y(i)*log(mu);
  }

  out = out -
    arma::as_scalar((0.01/2.0)*(arma::trans(x_beta)*x_beta)); /// prior beta

  return out;
}

void metropolis_beta_cpp(arma::vec beta_prev, arma::vec eps_prev, double nu_prev, arma::mat sigma,
                         arma::mat& varBeta, arma::vec& beta_samp, double c_beta,
                         Rcpp::List& params){

  double prevloglik, loglik, aux, ratio, u;

  int N = params["N"];
  int P = params["P"];

  arma::vec beta_proposal(P);

  // beta
  beta_proposal = rmvnorm_cpp(varBeta*c_beta) + beta_prev;

  // likelihood
  prevloglik = dens_beta_cpp(beta_prev, eps_prev, nu_prev, sigma, params);
  loglik = dens_beta_cpp(beta_proposal, eps_prev, nu_prev, sigma, params);

  aux = loglik - prevloglik;

  if(aux > 0){
    ratio = 1.0;
  } else {
    ratio = exp(aux);
  }

  u = Rcpp::as<double>(runif(1, 0, 1));

  if(u < ratio){
    beta_samp = beta_proposal;
  } else {
    beta_samp = beta_prev;
  }
}

//' @title dens_eps_cpp
//'
//' @description To usage in MCMC estimation
//'
// [[Rcpp::export]]
double dens_eps_cpp(arma::vec x_beta, arma::vec& x_eps, arma::vec& x_mu, double x_nu,
                    arma::mat& Q, arma::mat& sigma, arma::sp_mat& Qsparse, int i,
                    Rcpp::List& params){
  double out = 0.0, mu, logdet, sign;
  double jaceps = 0.0;

  /// Arguments
  int N = params["N"];
  int P = params["P"];
  int type = params["type"];

  arma::vec y = params["y"];
  arma::vec e = params["e"];
  arma::mat X = params["X"];

  /// Auxiliar vector
  arma::vec Xbet(N);

  /// Updating Xbeta
  Xbet = X*x_beta;

  /// Mu
  if(type == 1) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), x_nu, exp(Xbet(i))/x_nu, 1, 1);
  if(type == 2) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*x_nu, 1/x_nu, 1, 1);
  if(type == 3) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*exp(Xbet(i))*x_nu, 1/(x_nu*exp(Xbet(i))), 1, 1);
  if(type == 4) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(sigma(i, i)/x_nu), 1, 1);
  if(type == 5) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(1/x_nu), 1, 1);
  if(type == 6) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), 1/x_nu, exp(Xbet(i)), 1, 1);
  if(type == 7) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i)), 1/x_nu, 1, 1);

  x_mu(i) = mu;
  out = out - e(i)*mu + y(i)*log(mu);

  out = out -
    arma::as_scalar(0.5*arma::trans(x_eps)*Qsparse*x_eps);               /// eps

  return out;
}

void metropolis_eps_cpp(arma::vec beta_prev, arma::vec eps_prev, arma::vec mu_prev, double nu_prev,
                        arma::mat& Q_prev, arma::mat& sigma_prev, arma::sp_mat& Qsparse_prev, int i,
                        arma::mat& varEps, arma::vec& eps_samp, arma::vec& mu_samp, double c_eps,
                        Rcpp::List& params){

  double prevloglik, loglik, prevdens, dens, aux, ratio, u;

  int N = params["N"];
  int P = params["P"];

  double eps_proposal;

  arma::vec eps(N);
  arma::vec mu(N);

  eps_proposal = Rcpp::as<double>(rnorm(1, eps_prev(i), sqrt(varEps(i, i)*c_eps)));

  eps = eps_prev;
  mu = mu_prev;

  eps(i) = eps_proposal;
  eps(N-1) = -sum(eps.subvec(0, (N-2)));

  // likelihood
  prevloglik = dens_eps_cpp(beta_prev, eps_prev, mu, nu_prev, Q_prev, sigma_prev, Qsparse_prev, i, params);
  loglik = dens_eps_cpp(beta_prev, eps, mu, nu_prev, Q_prev, sigma_prev, Qsparse_prev, i, params);

  aux = loglik - prevloglik;

  if(aux > 0){
    ratio = 1.0;
  } else {
    ratio = exp(aux);
  }

  u = Rcpp::as<double>(runif(1, 0, 1));

  if(u < ratio){
    eps_samp(i) = eps(i);
    eps_samp(N-1) = eps(N-1);
    mu_samp(i) = mu(i);
  } else {
    eps_samp(i) = eps_prev(i);
    mu_samp(i) = mu_prev(i);
  }
}

//' @title dens_cpp
//'
//' @description To usage in MCMC estimation
//'
// [[Rcpp::export]]
double dens_rho_nu_cpp(arma::vec x_beta, arma::vec x_eps, arma::vec x_mu, arma::vec x_rho, double x_nu,
                       arma::mat& Q, arma::mat& sigma, arma::sp_mat& Qsparse, Rcpp::List& params){
  int i;
  double out = 0.0, mu, logdet, sign;
  double jaceps = 0.0;

  /// Arguments
  int N = params["N"];
  int P = params["P"];
  int type = params["type"];

  arma::vec y = params["y"];
  arma::vec e = params["e"];
  arma::mat X = params["X"];

  arma::vec tau_beta = params["tau_beta"];
  arma::vec mean_beta = params["mean_beta"];

  double eta_nu = params["eta_nu"];
  double psi_nu = params["psi_nu"];

  /// Auxiliar vector
  arma::vec Xbet(N);

  /// Updating Xbeta
  Xbet = X*x_beta;

  /// Variable
  x_nu = exp(x_nu);

  for(i = 0; i < N; i++){
    if(type == 1) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), x_nu, exp(Xbet(i))/x_nu, 1, 1);
    if(type == 2) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*x_nu, 1/x_nu, 1, 1);
    if(type == 3) mu = R::qgamma(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i))*exp(Xbet(i))*x_nu, 1/(x_nu*exp(Xbet(i))), 1, 1);
    if(type == 4) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(sigma(i, i)/x_nu), 1, 1);
    if(type == 5) mu = R::qlnorm(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), Xbet(i), sqrt(1/x_nu), 1, 1);
    if(type == 6) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), 1/x_nu, exp(Xbet(i)), 1, 1);
    if(type == 7) mu = R::qweibull(R::pnorm(x_eps(i), 0, sqrt(sigma(i, i)), 1, 1), exp(Xbet(i)), 1/x_nu, 1, 1);

    x_mu(i) = mu;
    out = out - e(i)*x_mu(i) + y(i)*log(x_mu(i));
  }

  arma::log_det(logdet, sign, Q);

  out = out +
    0.5*logdet -                                                          /// rho
    arma::as_scalar(0.5*arma::trans(x_eps)*Qsparse*x_eps) -               /// eps|rho
    (eta_nu - 1)*log(x_nu) - psi_nu*x_nu +                                /// prior nu
    log(x_nu);                                                            /// Jacobian nu

  return out;
}

void metropolis_rho_nu_cpp(arma::vec beta_prev, arma::vec eps_prev, arma::vec mu_prev, arma::vec rho_prev, double nu_prev,
                           arma::mat& Q_prev, arma::mat& sigma_prev, arma::sp_mat& Qsparse_prev,
                           arma::mat& varRho, double& varLogNu,
                           arma::vec& max_rho, arma::vec& range_prev, arma::vec range_rho_s, arma::vec range_rho_t, arma::vec range_rho_st,
                           int fix_rho_s, int fix_rho_t, int fix_rho_st,
                           double c_rho, double c_nu,
                           arma::vec& rho_samp, double& nu_samp,
                           Rcpp::List& params){

  double prevloglik, loglik, aux, ratio, u;
  arma::vec prevdens(3), dens(3);
  arma::vec range(3);
  Rcpp::List range_list;

  int N = params["N"];
  int P = params["P"];

  double nu_proposal;
  arma::vec rho_proposal(3);

  arma::mat Ws = params["Ws"];
  arma::mat Wt = params["Wt"];

  // rho_s
  if(!fix_rho_s){
    range_list = max_range_cpp(Ws, Wt, 0, 0, 0);
    range(0) = range_list["rho_s"];
    rho_proposal(0) = rtruncnorm_cpp(-range(0), range(0), rho_prev(0), varRho(0, 0)*c_rho);

    prevdens(0) = dtruncnorm_cpp(rho_prev(0), -range_prev(0), range_prev(0), rho_proposal(0), varRho(0, 0)*c_rho, 1L);
    dens(0) = dtruncnorm_cpp(rho_proposal(0), -range(0), range(0), rho_prev(0), varRho(0, 0)*c_rho, 1L);
  } else{
    rho_proposal(0) = rho_prev(0);

    prevdens(0) = 0;
    dens(0) = 0;
  }

  // rho_t
  if(!fix_rho_t){
    range_list = max_range_cpp(Ws, Wt, rho_proposal(0), 0, 0);
    range(1) = range_list["rho_t"];
    rho_proposal(1) = rtruncnorm_cpp(-range(1), range(1), rho_prev(1), varRho(1, 1)*c_rho);

    prevdens(1) = dtruncnorm_cpp(rho_prev(1), -range_prev(1), range_prev(1), rho_proposal(1), varRho(1, 1)*c_rho, 1L);
    dens(1) = dtruncnorm_cpp(rho_proposal(1), -range(1), range(1), rho_prev(1), varRho(1, 1)*c_rho, 1L);
  } else{
    rho_proposal(1) = rho_prev(1);

    prevdens(1) = 0;
    dens(1) = 0;
  }

  // rho_st
  if(!fix_rho_st){
    range_list = max_range_cpp(Ws, Wt, rho_proposal(0), rho_proposal(1), 0);
    range(2) = range_list["rho_st"];
    rho_proposal(2) = rtruncnorm_cpp(-range(2), range(2), rho_prev(2), varRho(2, 2)*c_rho);

    prevdens(2) = dtruncnorm_cpp(rho_prev(2), -range_prev(2), range_prev(2), rho_proposal(2), varRho(2, 2)*c_rho, 1L);
    dens(2) = dtruncnorm_cpp(rho_proposal(2), -range(2), range(2), rho_prev(2), varRho(2, 2)*c_rho, 1L);
  } else{
    rho_proposal(2) = rho_prev(2);

    prevdens(2) = 0;
    dens(2) = 0;
  }

  // rho's proposal is not symetric
  // double a1, a2;
  // a1 = dtmvnorm_cpp(rho_prev, -1*max_rho, max_rho, rho_proposal, varRho*c_rho);
  // a2 = dtmvnorm_cpp(rho_proposal, -1*max_rho, max_rho, rho_prev, varRho*c_rho);
  //
  // printf("%f - %f \n", a1, a2);

  /// Auxiliar matrix
  arma::mat sigma(N, N);
  arma::mat Q(N, N);

  sigma.zeros();
  Q.zeros();

  /// Building matrix
  buildQST_cpp(Q, Ws, Wt, rho_proposal(0), rho_proposal(1), rho_proposal(2));

  //Solve Q
  sigma = arma::inv_sympd(Q);
  arma::sp_mat Qsparse(Q);

  // nu
  nu_proposal = Rcpp::as<double>(rnorm(1, nu_prev, sqrt(varLogNu*c_nu)));

  // likelihood
  prevloglik = dens_rho_nu_cpp(beta_prev, eps_prev, mu_prev, rho_prev, nu_prev, Q_prev, sigma_prev, Qsparse_prev, params);
  loglik = dens_rho_nu_cpp(beta_prev, eps_prev, mu_prev, rho_proposal, nu_proposal, Q, sigma, Qsparse, params);

  // aux = loglik - prevloglik + prevdens - dens;
  aux = loglik - prevloglik + prevdens(0) + prevdens(1) + prevdens(2) - dens(0) - dens(1) - dens(2);

  if(aux > 0){
    ratio = 1.0;
  } else {
    ratio = exp(aux);
  }

  u = Rcpp::as<double>(runif(1, 0, 1));

  if(u < ratio){
    nu_samp = nu_proposal;
    rho_samp = rho_proposal;
    Q_prev = Q;
    sigma_prev = sigma;
    Qsparse_prev = Qsparse;
    range_prev = range;
  } else {
    nu_samp = nu_prev;
    rho_samp = rho_prev;
  }
}

//' @title POIMCAR
//'
//' @description Multivariate Poisson regression with CAR covariance structure
//'
//' @param nsim MCMC size
//' @param X Covariate matrix
//' @param y Response
//' @param E Offset
//' @param M Number of neighbors in each area
//' @param W Matrix with the neighborhood structure
//' @param N Dimension of the observations
//' @param P Dimension of the covariates
//' @param mean_beta Mean a priori to beta vector
//' @param tau_beta Variance a priori to beta vector
//' @param eta_nu Shape a priori to nu
//' @param psi_nu Rate a priori to nu
//' @param rangeRho Range to sample rho by using ARMS
//' @param type TGMRF type (1 to 6)
//' @param method ARMS (0) or Metropolis (1)
//' @param ninit Number of initial points in ARMS
//' @param maxpoint Maximum number of evaluation in each ARMS iteration
//' @param var_beta_met Variance of beta proposal
//' @param var_eps_met Variance of eps proposal
//' @param var_rho_met Variance of rho proposal
//' @param var_log_nu_met Variance of log(nu) proposal
//' @param tau Vector of tau parameters to construct Q
//'
// [[Rcpp::export]]
Rcpp::List poimcar_cpp(int nsim, int burnin, int thin,
                       arma::vec eps, arma::vec mu, arma::vec beta, double nu,
                       double rho_s, double rho_t, double rho_st,
                       arma::mat X, arma::vec y, arma::vec E,
                       arma::mat Ws, arma::mat Wt,
                       int N, int P,
                       arma::vec mean_beta, arma::vec tau_beta, double eta_nu, double psi_nu,
                       int fix_rho_s, int fix_rho_t, int fix_rho_st,
                       arma::vec range_rho_s, arma::vec range_rho_t, arma::vec range_rho_st,
                       int type,
                       arma::mat var_beta_met, arma::mat var_eps_met, arma::mat var_log_mu_met, arma::mat var_rho_met, double var_log_nu_met,
                       int verbose,
                       double c_beta, double c_eps, double c_mu, double c_nu, double c_rho){

  int i = 0, j = 0, k = 0, l = 0, ntot = 0, nt = 1;
  double Snu, Snuaux;
  arma::vec Seps(N), Smu(N), Sbeta(P), Srho(3), Xbeta(N);
  arma::vec max_rho(3), range(3);
  Rcpp::List range_aux;
  arma::mat Q(N, N), sigma(N, N);

  Q.zeros();
  sigma.zeros();

  ntot = nsim*thin + burnin;

  int stop_point = 100, limit = ntot, buffer = 100;
  int pct = ntot/10, cont = 1, pos = 1;

  arma::mat beta_out(nsim, P), eps_out(nsim, N), mu_out(nsim, N), rho_out(nsim, 3);
  arma::vec nu_out(nsim);
  arma::mat beta_aux(buffer, P), eps_aux(buffer, N), mu_aux(buffer, N), rho_aux(buffer, 3);
  arma::vec nu_aux(buffer);

  Rcpp::List params = Rcpp::List::create(Rcpp::Named("N") = N, Rcpp::Named("P") = P,
                                         Rcpp::Named("y") = y, Rcpp::Named("X") = X,
                                         Rcpp::Named("e") = E, Rcpp::Named("type") = type,
                                         Rcpp::Named("Ws") = Ws, Rcpp::Named("Wt") = Wt,
                                         Rcpp::Named("tau_beta") = tau_beta, Rcpp::Named("mean_beta") = mean_beta,
                                         Rcpp::Named("eta_nu") = eta_nu, Rcpp::Named("psi_nu") = psi_nu);

  /// Set the inital values to the model
  Snu = nu;
  Srho(0) = rho_s;
  Srho(1) = rho_t;
  Srho(2) = rho_st;

  Seps = eps;
  Smu = mu;
  Sbeta = beta;

  range_aux = max_range_cpp(Ws, Wt, 0, 0, 0);
  max_rho(0) = range_aux["rho_s"];
  max_rho(1) = range_aux["rho_t"];
  max_rho(2) = range_aux["rho_st"];
  range = max_rho;

  if(verbose){
    printf("Running...\n");
  }

  /// Calculating Q and sigma
  buildQST_cpp(Q, Ws, Wt, Srho(0), Srho(1), Srho(2));
  sigma = arma::inv_sympd(Q);
  arma::sp_mat Qsparse(Q);

  // Loop
  for(i = 0; i < ntot; i++){
    /// Sampling rho and nu together
    Snuaux = log(Snu);
    metropolis_rho_nu_cpp(Sbeta, Seps, Smu, Srho, Snuaux,
                          Q, sigma, Qsparse,
                          var_rho_met, var_log_nu_met,
                          max_rho, range, range_rho_s, range_rho_t, range_rho_st,
                          fix_rho_s, fix_rho_t, fix_rho_st,
                          c_rho, c_nu,
                          Srho, Snuaux,
                          params);
    Snu = exp(Snuaux);

    /// Sampling beta
    metropolis_beta_cpp(Sbeta, Seps, Snu,
                        sigma,
                        var_beta_met,
                        Sbeta, c_beta,
                        params);

    /// Sampling eps
    for(j = 0; j < (N-1); j++){
      metropolis_eps_cpp(Sbeta, Seps, Smu, Snu,
                         Q, sigma, Qsparse, j,
                         var_eps_met,
                         Seps, Smu, c_eps,
                         params);
    }

    /// Saving chain (+ 1 because vector start in position 0)
    if((i + 1) > burnin & (i + 1)%thin == 0){
      beta_out.row(k) = Sbeta.t();
      eps_out.row(k) = Seps.t();
      mu_out.row(k) = Smu.t();
      rho_out.row(k) = Srho.t();
      nu_out(k) = Snu;

      k += 1;
    }

    /// Auxiliar to calculate variances
    if(l > buffer - 1) l = 0;

    beta_aux.row(l) = Sbeta.t();
    eps_aux.row(l) = Seps.t();
    mu_aux.row(l) = Smu.t();
    rho_aux.row(l) = Srho.t();
    nu_aux(l) = Snuaux;
    l += 1;

    /// Recalculating variances
    if(i == stop_point & i < limit){
      nt = nt + 1;

      var_beta_met = var_beta_met*(nt - 1)/nt + cov(beta_aux.submat(0, 0, buffer-1, P-1))/nt;
      var_eps_met = var_eps_met*(nt - 1)/nt + cov(eps_aux.submat(0, 0, buffer-1, N-1))/nt;
      var_rho_met = var_rho_met*(nt - 1)/nt + cov(rho_aux.submat(0, 0, buffer-1, 2))/nt;

      for(j = 0; j < P; j++) var_beta_met(j, j) += 0.00001;
      for(j = 0; j < N; j++) var_eps_met(j, j) += 0.00001;
      for(j = 0; j < 3; j++) var_rho_met(j, j) += 0.00001;

      // var_rho_met(0, 1) = 0;
      // var_rho_met(0, 2) = 0;
      // var_rho_met(1, 0) = 0;
      // var_rho_met(1, 2) = 0;
      // var_rho_met(2, 0) = 0;
      // var_rho_met(2, 1) = 0;

      var_log_nu_met = var_log_nu_met*(nt - 1)/nt + arma::var(nu_aux.subvec(0, buffer-1))/nt + 0.00001;

      stop_point = stop_point + buffer;
    }

    /// %
    if(verbose){
      pos += 1;
      if(i == cont*pct){
        printf("%d%% completed... \n", cont*10);
        cont += 1;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("beta") = beta_out,
                            Rcpp::Named("eps") = eps_out,
                            Rcpp::Named("mu") = mu_out,
                            Rcpp::Named("nu") = nu_out,
                            Rcpp::Named("rho") = rho_out);
}
