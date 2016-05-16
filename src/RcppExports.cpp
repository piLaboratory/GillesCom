// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// create_community
void create_community(arma::vec abundance, arma::mat interaction, arma::vec K, arma::vec d0, arma::vec b, arma::vec m, double save_int);
RcppExport SEXP GillesCom_create_community(SEXP abundanceSEXP, SEXP interactionSEXP, SEXP KSEXP, SEXP d0SEXP, SEXP bSEXP, SEXP mSEXP, SEXP save_intSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
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
void load_community(arma::vec abundance, arma::mat trajectories, arma::mat interaction, arma::vec K, arma::vec d0, arma::vec b, arma::vec m, double time, double save_int, int cycles);
RcppExport SEXP GillesCom_load_community(SEXP abundanceSEXP, SEXP trajectoriesSEXP, SEXP interactionSEXP, SEXP KSEXP, SEXP d0SEXP, SEXP bSEXP, SEXP mSEXP, SEXP timeSEXP, SEXP save_intSEXP, SEXP cyclesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type abundance(abundanceSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type trajectories(trajectoriesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type interaction(interactionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d0(d0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Rcpp::traits::input_parameter< double >::type save_int(save_intSEXP);
    Rcpp::traits::input_parameter< int >::type cycles(cyclesSEXP);
    load_community(abundance, trajectories, interaction, K, d0, b, m, time, save_int, cycles);
    return R_NilValue;
END_RCPP
}
// abundance
arma::vec abundance();
RcppExport SEXP GillesCom_abundance() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(abundance());
    return __result;
END_RCPP
}
// K
arma::vec K();
RcppExport SEXP GillesCom_K() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(K());
    return __result;
END_RCPP
}
// d0
arma::vec d0();
RcppExport SEXP GillesCom_d0() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(d0());
    return __result;
END_RCPP
}
// birth
arma::vec birth();
RcppExport SEXP GillesCom_birth() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(birth());
    return __result;
END_RCPP
}
// migration
arma::vec migration();
RcppExport SEXP GillesCom_migration() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(migration());
    return __result;
END_RCPP
}
// elapsed_time
double elapsed_time();
RcppExport SEXP GillesCom_elapsed_time() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(elapsed_time());
    return __result;
END_RCPP
}
// elapsed_cycles
int elapsed_cycles();
RcppExport SEXP GillesCom_elapsed_cycles() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(elapsed_cycles());
    return __result;
END_RCPP
}
// save_int
double save_int();
RcppExport SEXP GillesCom_save_int() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(save_int());
    return __result;
END_RCPP
}
// get_interaction
arma::mat get_interaction();
RcppExport SEXP GillesCom_get_interaction() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(get_interaction());
    return __result;
END_RCPP
}
// trajectories
arma::mat trajectories();
RcppExport SEXP GillesCom_trajectories() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(trajectories());
    return __result;
END_RCPP
}
// Cbdm
void Cbdm(int count);
RcppExport SEXP GillesCom_Cbdm(SEXP countSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type count(countSEXP);
    Cbdm(count);
    return R_NilValue;
END_RCPP
}
// Tbdm
void Tbdm(double time);
RcppExport SEXP GillesCom_Tbdm(SEXP timeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type time(timeSEXP);
    Tbdm(time);
    return R_NilValue;
END_RCPP
}
