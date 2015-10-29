#include<iostream>
#include<fstream>
#include<string>
#include<RcppArmadillo.h>

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
    double time;
  public:
    arma::vec get_abundance() const {return abundance;}
    double get_time() const {return time;}
    Community(arma::vec _abundance, arma::mat _interaction,
        arma::vec _K, arma::vec _d0, arma::vec _b,
        arma::vec _m) {
      abundance = _abundance;
      interaction = _interaction;
      K = _K; d0 = _d0; b = _b; m = _m;
      time  = 0;
      dslope = (b-d0)/K; //slope of the density-dependent linear relation of death rate to N
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
      time += R::rexp(1.0 / sum(w));
      return;
    }
};

// Global var???
Community *C = NULL;

// [[Rcpp::export]]
void create_community(arma::vec abundance, arma::mat interaction,
        arma::vec K, arma::vec d0, arma::vec b,
        arma::vec m) {
  if (C!=NULL) std::cout << "Warning: overwriting previous Community" << std::endl;
  C = new Community(abundance, interaction, K, d0, b, m);
}

//[[Rcpp::export]]
arma::vec abundance() {
  if (C==NULL) return arma::vec(1, arma::fill::zeros);
  return C->get_abundance();
}

//[[Rcpp::export]]
double time() {
  if (C==NULL) return 0;
  return C->get_time();
}

//[[Rcpp::export]]
void bdm(int count = 1) {
  if (C==NULL) return;
  for (int i = 0; i < count; i ++)
    C->bdm();
}
