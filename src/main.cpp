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

// interpolates a function f(time) at instant t
double interpol (arma::vec time, arma::vec f, double t) {
    int i;
    // This is inefficient, should be changed some day
    // maybe we could cache the last "i" used?
    for (i = 0; i < time.n_elem; i++) 
        if (time[i] > t) break;
    if (i == 0) return f[0];
    if (i == time.n_elem) return f[f.n_elem - 1];
    // this formula incurs on cancelation errors if the time step is small!
    // find a better formula ASAP
    return f[i-1] + (f[i] - f[i-1]) * (t - time[i-1]) / (time[i] - time[i-1]); 
}

class Community {
  public:
    // abundance vector: counts how many individuals from each species
    arma::vec abundance;
    // trajectories matrix: one row for the abundances at each specified time
    arma::mat trajectories;
    // interaction matrix
    arma::mat interaction;
    // stochastic effects matrix
    arma::mat stochastic;
    // support capacity
    arma::vec K;
    // death rate when abundance=0
    arma::vec d0;
    // birth rate
    arma::vec b;
    // migration rate
    arma::vec m;
    // current time, last time and trajectories saving interval
    double time, save_int;
    // simulated cycles so far
    int cycles;
  public:
    // overloaded constructor for restoring a saved sim
    // TODO: collapse both constructors, as we really don't need two!
    Community(arma::vec _abundance, arma::mat _trajectories, arma::mat _interaction,
        arma::vec _K, arma::vec _d0, arma::vec _b,
        arma::vec _m, double _time, double _save_int, int _cycles, arma::mat _stochastic) {
      abundance = _abundance;
      trajectories = _trajectories;
      interaction = _interaction;
      K = _K; d0 = _d0; b = _b; m = _m;
      time  = _time; save_int = _save_int; cycles = _cycles; stochastic = _stochastic;
    }
    Community(arma::vec _abundance, arma::mat _interaction,
        arma::vec _K, arma::vec _d0, arma::vec _b,
        arma::vec _m, double _save_int, arma::mat _stochastic) {
      abundance = _abundance;
      trajectories = abundance.t();
      interaction = _interaction;
      K = _K; d0 = _d0; b = _b; m = _m;
      time  = 0; cycles = 0; save_int = _save_int; stochastic = _stochastic;
    }
    Community(arma::vec _abundance, arma::mat _interaction, double _save_int) {
        abundance = _abundance; save_int = _save_int;
        interaction = _interaction;
        time = 0; cycles = 0; trajectories = abundance.t();
    }
    void saveHistory() {
      trajectories.insert_rows(trajectories.n_rows, abundance.t());
    }
    void bdm() {
      double mult; arma::vec instant_b = b;
      // stochastic multiplier for the birth rate:
      if (stochastic.n_elem > 0 ) {
        mult = interpol ( stochastic.col(0), stochastic.col(1), time );
        instant_b = b * mult;
      }
      arma::vec dslope = (instant_b-d0)/K; //slope of the density-dependent linear relation of death rate to N
      // % performs element-wise multiplication
      arma::vec d = d0 + dslope % (abundance.t() * interaction).t();
      for (int i = 0; i < abundance.n_elem; i++) if(abundance(i) == 0) d(i) = 0;
      //Gillespie weights for each specie, which are the sum of their rates
      arma::vec w = (abundance % (instant_b + d)) + m; 
      //sampling which species will suffer the next action, proportionaly to their weights
      int c = Csample(w);
      // Should the selected species gain or lose an individual?
      double choice = R::runif(0,1);
      if ( choice > (instant_b(c)*abundance(c)+m(c)) / (instant_b(c)*abundance(c)+m(c)+d(c)*abundance(c)))
        abundance(c) --;
      else 
        abundance(c) ++;
      // advances the simulation clock
      double elapsed = R::rexp(1.0 / sum(w));
      if (elapsed > save_int) warning("Time elapsed larger than save interval!");
      // only saves trajectories if we have completed a saving period
      if (((int) (time / save_int)) !=  ((int) ((time+elapsed)/save_int)))
        saveHistory();
      time += elapsed; cycles += 1;
      return;
    }
};

//void Cbdm(int count = 1) {
//  if (C==NULL) return;
//  for (int i = 0; i < count; i ++)
//    C->bdm();
//}

//void Tbdm(double time) {
//  if (C==NULL) return;
//  while (C->get_time() < time)
//    C->bdm();
//}

RCPP_MODULE (Community) {
    using namespace Rcpp;
    class_<Community> ("Community")
        .constructor<arma::vec, arma::mat, double> ()
        .field("abundance", &Community::abundance)
        .field("trajectories", &Community::trajectories)
        .field("interaction", &Community::interaction)
        .field("stochastic", &Community::stochastic)
        .field("K", &Community::K)
        .field("d0", &Community::d0)
        .field("b", &Community::b)
        .field("m", &Community::m)
        .field("time", &Community::time)
        .field("save_int", &Community::save_int)
        .field("cycles", &Community::cycles)
        .method("bdm", &Community::bdm)
        ;
}

