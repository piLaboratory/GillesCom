#include<iostream>
#include<fstream>
#include<string>
#include<RcppArmadillo.h>

// Samples a single integer from 0:(size-1) with probabilities prob
// similar to R function sample(1:N-1, 1, prob)
// Notice that prob does not need to be normalized
// [[Rcpp::export]]
int sample (int size, arma::vec prob) {
  double choice = R::runif(0,1);
  prob = prob / sum(prob);
  double cum = 0;
  for (int i = 0; i < prob.size(); i++) {
    cum += prob(i);
    if(cum > choice) return i;
  }
  // should never reach
  return prob.size();
}

class Community {
  private:
    arma::vec abundance;
    arma::mat interaction;
    arma::vec K;
    arma::vec d0;
    arma::vec b;
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
    }
    arma::vec bdm() {
      arma::vec d = (b-d0)/K; //slope of the density-dependent linear relation of death rate to N
      arma::vec N = trans(abundance.t() * interaction);
      for(int i = 0; i < abundance.n_elem; i++) if(abundance(i) == 0) N(i) = 0;
      d = d0 + d % N; // element wise multiplication
      arma::vec w = abundance % (b + d) + m; //Gillespie weights for each specie, which are the sum of their rates
      int c = sample(abundance.size(), w);//sampling which species will suffer the next action, proportionaly to their weights
      abundance(c) ++;
      
//    d <- (b-d0)/K # 
//    dt <- d0+d*N # death rates of each species
//    w <- N0*(b+dt) + m  # 
//    i <- sample((1:length(N)), size=1, prob=w) ## 
//    N0[i] <- N0[i] + sample(c(1, -1), size=1, prob= c(b[i]*N0[i]+m[i], dt[i]*N0[i])) ## Sampling if the selected species will gain or loss an individual
      return abundance;
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
  if (C==NULL) return 0;
  return C->get_abundance();
}

//[[Rcpp::export]]
double time() {
  if (C==NULL) return 0;
  return C->get_time();
}

//[[Rcpp::export]]
arma::vec bdm() {
  if (C==NULL) return 0;
  return C->bdm();
}
