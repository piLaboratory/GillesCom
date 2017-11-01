// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// create_community
void create_community(arma::vec abundance, arma::mat interaction, arma::vec K, arma::vec d0, arma::vec b, arma::vec m, double save_int);
RcppExport SEXP _GillesCom_create_community(SEXP abundanceSEXP, SEXP interactionSEXP, SEXP KSEXP, SEXP d0SEXP, SEXP bSEXP, SEXP mSEXP, SEXP save_intSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type abundance(abundanceSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type interaction(interactionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type save_int(save_intSEXP);
    create_community(abundance, interaction, K, d0, b, m, save_int);
    return R_NilValue;
END_RCPP
}
// load_community
void load_community(arma::vec abundance, arma::mat trajectories, arma::mat interaction, arma::vec K, arma::vec d0, arma::vec b, arma::vec m, double time, double save_int);
RcppExport SEXP _GillesCom_load_community(SEXP abundanceSEXP, SEXP trajectoriesSEXP, SEXP interactionSEXP, SEXP KSEXP, SEXP d0SEXP, SEXP bSEXP, SEXP mSEXP, SEXP timeSEXP, SEXP save_intSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type abundance(abundanceSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type interaction(interactionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type save_int(save_intSEXP);
    load_community(abundance, trajectories, interaction, K, d0, b, m, time, save_int);
    return R_NilValue;
END_RCPP
}
// abundance
arma::vec abundance();
RcppExport SEXP _GillesCom_abundance() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(abundance());
    return rcpp_result_gen;
END_RCPP
}
// K
arma::vec K();
RcppExport SEXP _GillesCom_K() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(K());
    return rcpp_result_gen;
END_RCPP
}
// d0
arma::vec d0();
RcppExport SEXP _GillesCom_d0() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(d0());
    return rcpp_result_gen;
END_RCPP
}
// birth
arma::vec birth();
RcppExport SEXP _GillesCom_birth() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(birth());
    return rcpp_result_gen;
END_RCPP
}
// migration
arma::vec migration();
RcppExport SEXP _GillesCom_migration() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(migration());
    return rcpp_result_gen;
END_RCPP
}
// elapsed_time
double elapsed_time();
RcppExport SEXP _GillesCom_elapsed_time() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(elapsed_time());
    return rcpp_result_gen;
END_RCPP
}
// save_int
double save_int();
RcppExport SEXP _GillesCom_save_int() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(save_int());
    return rcpp_result_gen;
END_RCPP
}
// get_interaction
arma::mat get_interaction();
RcppExport SEXP _GillesCom_get_interaction() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_interaction());
    return rcpp_result_gen;
END_RCPP
}
// trajectories
arma::mat trajectories();
RcppExport SEXP _GillesCom_trajectories() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(trajectories());
    return rcpp_result_gen;
END_RCPP
}
// Cbdm
void Cbdm(int count);
RcppExport SEXP _GillesCom_Cbdm(SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type count(countSEXP);
    Cbdm(count);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GillesCom_create_community", (DL_FUNC) &_GillesCom_create_community, 7},
    {"_GillesCom_load_community", (DL_FUNC) &_GillesCom_load_community, 9},
    {"_GillesCom_abundance", (DL_FUNC) &_GillesCom_abundance, 0},
    {"_GillesCom_K", (DL_FUNC) &_GillesCom_K, 0},
    {"_GillesCom_d0", (DL_FUNC) &_GillesCom_d0, 0},
    {"_GillesCom_birth", (DL_FUNC) &_GillesCom_birth, 0},
    {"_GillesCom_migration", (DL_FUNC) &_GillesCom_migration, 0},
    {"_GillesCom_elapsed_time", (DL_FUNC) &_GillesCom_elapsed_time, 0},
    {"_GillesCom_save_int", (DL_FUNC) &_GillesCom_save_int, 0},
    {"_GillesCom_get_interaction", (DL_FUNC) &_GillesCom_get_interaction, 0},
    {"_GillesCom_trajectories", (DL_FUNC) &_GillesCom_trajectories, 0},
    {"_GillesCom_Cbdm", (DL_FUNC) &_GillesCom_Cbdm, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_GillesCom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
