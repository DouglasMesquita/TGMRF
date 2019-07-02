// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// max_range_cpp
Rcpp::List max_range_cpp(arma::mat Ws, arma::mat Wt, double rho_s, double rho_t, double rho_st);
RcppExport SEXP _TGMRF_max_range_cpp(SEXP WsSEXP, SEXP WtSEXP, SEXP rho_sSEXP, SEXP rho_tSEXP, SEXP rho_stSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ws(WsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Wt(WtSEXP);
    Rcpp::traits::input_parameter< double >::type rho_s(rho_sSEXP);
    Rcpp::traits::input_parameter< double >::type rho_t(rho_tSEXP);
    Rcpp::traits::input_parameter< double >::type rho_st(rho_stSEXP);
    rcpp_result_gen = Rcpp::wrap(max_range_cpp(Ws, Wt, rho_s, rho_t, rho_st));
    return rcpp_result_gen;
END_RCPP
}
// cov_cpp
void cov_cpp(arma::mat A, arma::mat& B);
RcppExport SEXP _TGMRF_cov_cpp(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type B(BSEXP);
    cov_cpp(A, B);
    return R_NilValue;
END_RCPP
}
// buildQC_cpp
arma::mat buildQC_cpp(arma::vec M, arma::mat W, arma::vec rho);
RcppExport SEXP _TGMRF_buildQC_cpp(SEXP MSEXP, SEXP WSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(buildQC_cpp(M, W, rho));
    return rcpp_result_gen;
END_RCPP
}
// buildQL_cpp
arma::mat buildQL_cpp(arma::mat R, arma::vec rho);
RcppExport SEXP _TGMRF_buildQL_cpp(SEXP RSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(buildQL_cpp(R, rho));
    return rcpp_result_gen;
END_RCPP
}
// buildQST_cpp
void buildQST_cpp(arma::mat& Q, arma::mat& Ws, arma::mat& Wt, double& rho_s, double& rho_t, double& rho_st);
RcppExport SEXP _TGMRF_buildQST_cpp(SEXP QSEXP, SEXP WsSEXP, SEXP WtSEXP, SEXP rho_sSEXP, SEXP rho_tSEXP, SEXP rho_stSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Ws(WsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Wt(WtSEXP);
    Rcpp::traits::input_parameter< double& >::type rho_s(rho_sSEXP);
    Rcpp::traits::input_parameter< double& >::type rho_t(rho_tSEXP);
    Rcpp::traits::input_parameter< double& >::type rho_st(rho_stSEXP);
    buildQST_cpp(Q, Ws, Wt, rho_s, rho_t, rho_st);
    return R_NilValue;
END_RCPP
}
// rmvnorm_cpp
arma::vec rmvnorm_cpp(arma::mat sigma);
RcppExport SEXP _TGMRF_rmvnorm_cpp(SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm_cpp(sigma));
    return rcpp_result_gen;
END_RCPP
}
// rtruncnorm_cpp
double rtruncnorm_cpp(double a, double b, double mean, double var);
RcppExport SEXP _TGMRF_rtruncnorm_cpp(SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(rtruncnorm_cpp(a, b, mean, var));
    return rcpp_result_gen;
END_RCPP
}
// rtmvnorm_cpp
arma::vec rtmvnorm_cpp(arma::vec a, arma::vec b, arma::vec mean, arma::mat var);
RcppExport SEXP _TGMRF_rtmvnorm_cpp(SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(rtmvnorm_cpp(a, b, mean, var));
    return rcpp_result_gen;
END_RCPP
}
// dtruncnorm
double dtruncnorm(double x, double a, double b, double mean, double var);
RcppExport SEXP _TGMRF_dtruncnorm(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(dtruncnorm(x, a, b, mean, var));
    return rcpp_result_gen;
END_RCPP
}
// dtruncnorm_cpp
double dtruncnorm_cpp(double x, double a, double b, double mean, double var);
RcppExport SEXP _TGMRF_dtruncnorm_cpp(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(dtruncnorm_cpp(x, a, b, mean, var));
    return rcpp_result_gen;
END_RCPP
}
// dtmvnorm_cpp
double dtmvnorm_cpp(arma::vec x, arma::vec a, arma::vec b, arma::vec mean, arma::mat var);
RcppExport SEXP _TGMRF_dtmvnorm_cpp(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(dtmvnorm_cpp(x, a, b, mean, var));
    return rcpp_result_gen;
END_RCPP
}
// dens_eps_pois_cpp
double dens_eps_pois_cpp(double x, Rcpp::List& params, int& ele, arma::vec& eps, double& nu, arma::sp_mat& Q, arma::mat dsigma, arma::vec& Xbeta);
RcppExport SEXP _TGMRF_dens_eps_pois_cpp(SEXP xSEXP, SEXP paramsSEXP, SEXP eleSEXP, SEXP epsSEXP, SEXP nuSEXP, SEXP QSEXP, SEXP dsigmaSEXP, SEXP XbetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int& >::type ele(eleSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dsigma(dsigmaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Xbeta(XbetaSEXP);
    rcpp_result_gen = Rcpp::wrap(dens_eps_pois_cpp(x, params, ele, eps, nu, Q, dsigma, Xbeta));
    return rcpp_result_gen;
END_RCPP
}
// dens_cpp
double dens_cpp(arma::vec x_beta, arma::vec x_eps, arma::vec x_rho, double x_nu, Rcpp::List& params, arma::vec& mu, arma::mat& cor_mat);
RcppExport SEXP _TGMRF_dens_cpp(SEXP x_betaSEXP, SEXP x_epsSEXP, SEXP x_rhoSEXP, SEXP x_nuSEXP, SEXP paramsSEXP, SEXP muSEXP, SEXP cor_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x_beta(x_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_eps(x_epsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_rho(x_rhoSEXP);
    Rcpp::traits::input_parameter< double >::type x_nu(x_nuSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type cor_mat(cor_matSEXP);
    rcpp_result_gen = Rcpp::wrap(dens_cpp(x_beta, x_eps, x_rho, x_nu, params, mu, cor_mat));
    return rcpp_result_gen;
END_RCPP
}
// poimcar_cpp
Rcpp::List poimcar_cpp(int nsim, int burnin, int thin, arma::vec eps, arma::vec mu, arma::vec beta, double nu, double rho_s, double rho_t, double rho_st, arma::mat X, arma::vec y, arma::vec E, arma::mat Ws, arma::mat Wt, int N, int P, arma::vec mean_beta, arma::vec tau_beta, double eta_nu, double psi_nu, int fix_rho_s, int fix_rho_t, int fix_rho_st, arma::vec range_rho_s, arma::vec range_rho_t, arma::vec range_rho_st, int type, arma::mat var_beta_met, arma::mat var_eps_met, arma::mat var_rho_met, double var_log_nu_met, int verbose, double c_beta, double c_eps, double c_nu, double c_rho);
RcppExport SEXP _TGMRF_poimcar_cpp(SEXP nsimSEXP, SEXP burninSEXP, SEXP thinSEXP, SEXP epsSEXP, SEXP muSEXP, SEXP betaSEXP, SEXP nuSEXP, SEXP rho_sSEXP, SEXP rho_tSEXP, SEXP rho_stSEXP, SEXP XSEXP, SEXP ySEXP, SEXP ESEXP, SEXP WsSEXP, SEXP WtSEXP, SEXP NSEXP, SEXP PSEXP, SEXP mean_betaSEXP, SEXP tau_betaSEXP, SEXP eta_nuSEXP, SEXP psi_nuSEXP, SEXP fix_rho_sSEXP, SEXP fix_rho_tSEXP, SEXP fix_rho_stSEXP, SEXP range_rho_sSEXP, SEXP range_rho_tSEXP, SEXP range_rho_stSEXP, SEXP typeSEXP, SEXP var_beta_metSEXP, SEXP var_eps_metSEXP, SEXP var_rho_metSEXP, SEXP var_log_nu_metSEXP, SEXP verboseSEXP, SEXP c_betaSEXP, SEXP c_epsSEXP, SEXP c_nuSEXP, SEXP c_rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< double >::type rho_s(rho_sSEXP);
    Rcpp::traits::input_parameter< double >::type rho_t(rho_tSEXP);
    Rcpp::traits::input_parameter< double >::type rho_st(rho_stSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Ws(WsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Wt(WtSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mean_beta(mean_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau_beta(tau_betaSEXP);
    Rcpp::traits::input_parameter< double >::type eta_nu(eta_nuSEXP);
    Rcpp::traits::input_parameter< double >::type psi_nu(psi_nuSEXP);
    Rcpp::traits::input_parameter< int >::type fix_rho_s(fix_rho_sSEXP);
    Rcpp::traits::input_parameter< int >::type fix_rho_t(fix_rho_tSEXP);
    Rcpp::traits::input_parameter< int >::type fix_rho_st(fix_rho_stSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type range_rho_s(range_rho_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type range_rho_t(range_rho_tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type range_rho_st(range_rho_stSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_beta_met(var_beta_metSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_eps_met(var_eps_metSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_rho_met(var_rho_metSEXP);
    Rcpp::traits::input_parameter< double >::type var_log_nu_met(var_log_nu_metSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type c_beta(c_betaSEXP);
    Rcpp::traits::input_parameter< double >::type c_eps(c_epsSEXP);
    Rcpp::traits::input_parameter< double >::type c_nu(c_nuSEXP);
    Rcpp::traits::input_parameter< double >::type c_rho(c_rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(poimcar_cpp(nsim, burnin, thin, eps, mu, beta, nu, rho_s, rho_t, rho_st, X, y, E, Ws, Wt, N, P, mean_beta, tau_beta, eta_nu, psi_nu, fix_rho_s, fix_rho_t, fix_rho_st, range_rho_s, range_rho_t, range_rho_st, type, var_beta_met, var_eps_met, var_rho_met, var_log_nu_met, verbose, c_beta, c_eps, c_nu, c_rho));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TGMRF_max_range_cpp", (DL_FUNC) &_TGMRF_max_range_cpp, 5},
    {"_TGMRF_cov_cpp", (DL_FUNC) &_TGMRF_cov_cpp, 2},
    {"_TGMRF_buildQC_cpp", (DL_FUNC) &_TGMRF_buildQC_cpp, 3},
    {"_TGMRF_buildQL_cpp", (DL_FUNC) &_TGMRF_buildQL_cpp, 2},
    {"_TGMRF_buildQST_cpp", (DL_FUNC) &_TGMRF_buildQST_cpp, 6},
    {"_TGMRF_rmvnorm_cpp", (DL_FUNC) &_TGMRF_rmvnorm_cpp, 1},
    {"_TGMRF_rtruncnorm_cpp", (DL_FUNC) &_TGMRF_rtruncnorm_cpp, 4},
    {"_TGMRF_rtmvnorm_cpp", (DL_FUNC) &_TGMRF_rtmvnorm_cpp, 4},
    {"_TGMRF_dtruncnorm", (DL_FUNC) &_TGMRF_dtruncnorm, 5},
    {"_TGMRF_dtruncnorm_cpp", (DL_FUNC) &_TGMRF_dtruncnorm_cpp, 5},
    {"_TGMRF_dtmvnorm_cpp", (DL_FUNC) &_TGMRF_dtmvnorm_cpp, 5},
    {"_TGMRF_dens_eps_pois_cpp", (DL_FUNC) &_TGMRF_dens_eps_pois_cpp, 8},
    {"_TGMRF_dens_cpp", (DL_FUNC) &_TGMRF_dens_cpp, 7},
    {"_TGMRF_poimcar_cpp", (DL_FUNC) &_TGMRF_poimcar_cpp, 37},
    {NULL, NULL, 0}
};

RcppExport void R_init_TGMRF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
