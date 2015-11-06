#include<iostream>
#include<fstream>
#include<string>
#include<RcppArmadillo.h>

Rcpp::Function warning("warning"); 

// Samples a single integer from 0:(size-1) with probabilities prob
// similar to R function sample(1:length(prob) - 1, 1, prob)
// Notice that prob does not need to be normalized
int Csample (arma::vec prob) {
  double choice = R::runif(0,1);
  prob = prob / sum(prob);
  double cum = 0;
  for (int i = 0; i < prob.n_elem; i++) {
    cum += prob(i);
    if(cum > choice) return i;
  }
  // should never reach
  return prob.n_elem;
}

class Community {
  private:
    // abundance vector: counts how many individuals from each species
    arma::vec abundance;
    // history matrix: one row for the abundances at each specified time
    arma::mat history;
    // interaction matrix
    arma::mat interaction;
    // support capacity
    arma::vec K;
    // death rate when abundance=0
    arma::vec d0;
    // current death rate ("cached" for bdm())
    arma::vec dslope;
    // birth rate
    arma::vec b;
    // migration rate
    arma::vec m;
    // current time, last time and history saving interval
    double time, oldtime, save_int;
  public:
    arma::vec get_abundance() const {return abundance;}
    arma::vec get_K() const {return K;}
    arma::mat get_history() const {return history;}
    double get_time() const {return time;}
    Community(arma::vec _abundance, arma::mat _interaction,
        arma::vec _K, arma::vec _d0, arma::vec _b,
        arma::vec _m, double _save_int) {
      abundance = _abundance;
      history = abundance.t();
      interaction = _interaction;
      K = _K; d0 = _d0; b = _b; m = _m;
      time  = 0; oldtime = 0; save_int = _save_int;
      dslope = (b-d0)/K; //slope of the density-dependent linear relation of death rate to N
    }
    void saveHistory() {
      history.insert_rows(history.n_rows, abundance.t());
    }
    void bdm() {
      // % performs element-wise multiplication
      arma::vec d = d0 + dslope % (abundance.t() * interaction).t();
      for (int i = 0; i < abundance.n_elem; i++) if(abundance(i) == 0) d(i) = 0;
      //Gillespie weights for each specie, which are the sum of their rates
      arma::vec w = (abundance % (b + d)) + m; 
      //sampling which species will suffer the next action, proportionaly to their weights
      int c = Csample(w);
      // Should the selected species gain or lose an individual?
      double choice = R::runif(0,1);
      if ( choice > (b(c)*abundance(c)+m(c)) / (b(c)*abundance(c)+m(c)+d(c)*abundance(c)))
        abundance(c) --;
      else 
        abundance(c) ++;
      // advances the simulation clock
      double elapsed = R::rexp(1.0 / sum(w));
      if (elapsed > save_int) warning("Time elapsed larger than save interval!");
      // only saves history if we have completed a saving period
      if (((int) (time / save_int)) !=  ((int) ((time+elapsed)/save_int)))
        saveHistory();
      time += elapsed;
      return;
    }
};

// Global var???
Community *C = NULL;

// [[Rcpp::export]]
void create_community(arma::vec abundance, arma::mat interaction,
        arma::vec K, arma::vec d0, arma::vec b,
        arma::vec m, double save_int) {
  if (C!=NULL) warning("Warning: overwriting previous Community");
  C = new Community(abundance, interaction, K, d0, b, m, save_int);
}

//[[Rcpp::export]]
arma::vec abundance() {
  if (C==NULL) return arma::vec(1, arma::fill::zeros);
  return C->get_abundance();
}

//[[Rcpp::export]]
arma::vec K() {
  if (C==NULL) return arma::vec(1, arma::fill::zeros);
  return C->get_K();
}

//[[Rcpp::export]]
double time() {
  if (C==NULL) return 0;
  return C->get_time();
}

//[[Rcpp::export]]
arma::mat history() {
  if (C==NULL) return arma::mat(1, 1, arma::fill::zeros);
  return C->get_history();
}

//[[Rcpp::export]]
void bdm(int count = 1) {
  if (C==NULL) return;
  for (int i = 0; i < count; i ++)
    C->bdm();
}
