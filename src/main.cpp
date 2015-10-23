#include<iostream>
#include<fstream>
#include<string>
#include<Rcpp.h>

class Community {
  private:
    Rcpp::NumericVector abundance;
    Rcpp::NumericMatrix interaction;
  public:
    Rcpp::NumericVector const get_abundance() {return abundance;}
    int nspecies() {return abundance.size();}
    Community(Rcpp::NumericVector _abundance, Rcpp::NumericMatrix _interaction) {
      abundance = _abundance;
      interaction = _interaction;
    }
};

// Global var???
Community *C = NULL;

// [[Rcpp::export]]
void create_community(Rcpp::NumericVector abundance, Rcpp::NumericMatrix interaction) {
  if (C!=NULL) std::cout << "Warning: overwriting previous Community" << std::endl;
  C = new Community(abundance, interaction);
}

//[[Rcpp::export]]
Rcpp::NumericVector abundance() {
  if (C==NULL) return 0;
  return C->get_abundance();
}
