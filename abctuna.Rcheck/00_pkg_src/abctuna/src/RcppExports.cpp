// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// initpdyn
RcppExport SEXP initpdyn(SEXP dm_, SEXP srec_, SEXP psi_, SEXP M_, SEXP mata_, SEXP wta_, SEXP sela_, SEXP hinit_);
RcppExport SEXP _abctuna_initpdyn(SEXP dm_SEXP, SEXP srec_SEXP, SEXP psi_SEXP, SEXP M_SEXP, SEXP mata_SEXP, SEXP wta_SEXP, SEXP sela_SEXP, SEXP hinit_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dm_(dm_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type srec_(srec_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type psi_(psi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type M_(M_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type mata_(mata_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type wta_(wta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sela_(sela_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type hinit_(hinit_SEXP);
    rcpp_result_gen = Rcpp::wrap(initpdyn(dm_, srec_, psi_, M_, mata_, wta_, sela_, hinit_));
    return rcpp_result_gen;
END_RCPP
}
// msypdyn
RcppExport SEXP msypdyn(SEXP dm_, SEXP srec_, SEXP R0_, SEXP hh_, SEXP psi_, SEXP M_, SEXP mata_, SEXP wta_, SEXP sela_, SEXP hinit_);
RcppExport SEXP _abctuna_msypdyn(SEXP dm_SEXP, SEXP srec_SEXP, SEXP R0_SEXP, SEXP hh_SEXP, SEXP psi_SEXP, SEXP M_SEXP, SEXP mata_SEXP, SEXP wta_SEXP, SEXP sela_SEXP, SEXP hinit_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dm_(dm_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type srec_(srec_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type R0_(R0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type hh_(hh_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type psi_(psi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type M_(M_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type mata_(mata_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type wta_(wta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sela_(sela_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type hinit_(hinit_SEXP);
    rcpp_result_gen = Rcpp::wrap(msypdyn(dm_, srec_, R0_, hh_, psi_, M_, mata_, wta_, sela_, hinit_));
    return rcpp_result_gen;
END_RCPP
}
// pdyn
RcppExport SEXP pdyn(SEXP dm_, SEXP srec_, SEXP R0_, SEXP hh_, SEXP psi_, SEXP epsr_, SEXP spr0_, SEXP M_, SEXP mata_, SEXP wta_, SEXP sela_, SEXP Ninit_, SEXP Cb_);
RcppExport SEXP _abctuna_pdyn(SEXP dm_SEXP, SEXP srec_SEXP, SEXP R0_SEXP, SEXP hh_SEXP, SEXP psi_SEXP, SEXP epsr_SEXP, SEXP spr0_SEXP, SEXP M_SEXP, SEXP mata_SEXP, SEXP wta_SEXP, SEXP sela_SEXP, SEXP Ninit_SEXP, SEXP Cb_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dm_(dm_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type srec_(srec_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type R0_(R0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type hh_(hh_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type psi_(psi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type epsr_(epsr_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type spr0_(spr0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type M_(M_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type mata_(mata_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type wta_(wta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sela_(sela_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ninit_(Ninit_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Cb_(Cb_SEXP);
    rcpp_result_gen = Rcpp::wrap(pdyn(dm_, srec_, R0_, hh_, psi_, epsr_, spr0_, M_, mata_, wta_, sela_, Ninit_, Cb_));
    return rcpp_result_gen;
END_RCPP
}
// pdynlfcpue
RcppExport SEXP pdynlfcpue(SEXP dm_, SEXP srec_, SEXP R0_, SEXP hh_, SEXP psi_, SEXP epsr_, SEXP spr0_, SEXP M_, SEXP mata_, SEXP wta_, SEXP sela_, SEXP Ninit_, SEXP Cb_, SEXP pla_, SEXP fref_);
RcppExport SEXP _abctuna_pdynlfcpue(SEXP dm_SEXP, SEXP srec_SEXP, SEXP R0_SEXP, SEXP hh_SEXP, SEXP psi_SEXP, SEXP epsr_SEXP, SEXP spr0_SEXP, SEXP M_SEXP, SEXP mata_SEXP, SEXP wta_SEXP, SEXP sela_SEXP, SEXP Ninit_SEXP, SEXP Cb_SEXP, SEXP pla_SEXP, SEXP fref_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type dm_(dm_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type srec_(srec_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type R0_(R0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type hh_(hh_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type psi_(psi_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type epsr_(epsr_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type spr0_(spr0_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type M_(M_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type mata_(mata_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type wta_(wta_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type sela_(sela_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Ninit_(Ninit_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Cb_(Cb_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type pla_(pla_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type fref_(fref_SEXP);
    rcpp_result_gen = Rcpp::wrap(pdynlfcpue(dm_, srec_, R0_, hh_, psi_, epsr_, spr0_, M_, mata_, wta_, sela_, Ninit_, Cb_, pla_, fref_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_abctuna_initpdyn", (DL_FUNC) &_abctuna_initpdyn, 8},
    {"_abctuna_msypdyn", (DL_FUNC) &_abctuna_msypdyn, 10},
    {"_abctuna_pdyn", (DL_FUNC) &_abctuna_pdyn, 13},
    {"_abctuna_pdynlfcpue", (DL_FUNC) &_abctuna_pdynlfcpue, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_abctuna(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}