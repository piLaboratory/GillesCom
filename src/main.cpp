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
      // Caching the time interval for environmental interpolation
      // As this is only a cache, it does not need file persistence
      int cached_i;
      // interpolates the function environmental(t) at instant "time"
      double interpol_envir () {
          int i;
          for (i = cached_i; i < environmental.n_rows; i++) 
              if (environmental(i,0) > time) break;
          cached_i = i - 1;
          if (i == 0) return environmental.col(1)[0];
          if (i == environmental.n_rows) return environmental(environmental.n_rows - 1,1);
          return environmental(i-1,1) + (environmental(i,1) - environmental(i-1,1))
                 * (time - environmental(i-1,0)) / (environmental(i,0) - environmental(i-1,0)); 
      }
  public:
    // abundance vector: counts how many individuals from each species
    arma::vec abundance;
    // trajectories matrix: one row for the abundances at each specified time
    Rcpp::List trajectories;
    // interaction matrix
    arma::mat interaction;
    // environmental effects matrix
    arma::mat environmental;
    // environmental strength vector
    arma::vec environmental_strength;
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
    // TODO: Learn how to use a constructor with 7+ parameters in the Rcpp module??
    Community(arma::vec _abundance, arma::mat _interaction, double _save_int) {
        cached_i = 0;
        abundance = _abundance; save_int = _save_int;
        interaction = _interaction;
        time = 0; cycles = 0; 
        trajectories = Rcpp::List::create (abundance.t(), 
                arma::vec(1, arma::fill::zeros),
                arma::vec(1, arma::fill::zeros));
    }
    void saveHistory() {
////  This code could use a cleanup, but the performance impact here is negligible
      arma::mat tt = trajectories[0];
      tt.insert_rows(tt.n_rows, abundance.t());
      trajectories[0] = tt;
      arma::vec t = trajectories[1];
      t.insert_rows(t.n_rows, 1);
      t[ t.n_rows - 1 ] = time;
      trajectories[1] = t;
      arma::vec t2 = trajectories[2];
      t2.insert_rows(t2.n_rows, 1);
      t2[ t2.n_rows - 1 ] = cycles;
      trajectories[2] = t2;
    }
    void bdm() {
      double mult; arma::vec instant_K = K;
      // environmental multiplier for K:
      if (environmental.n_elem > 0 ) {
        mult = interpol_envir();
        instant_K = K + ((mult - 1.0) * K) % environmental_strength;
      }
      arma::vec dslope = (b-d0)/instant_K; //slope of the density-dependent linear relation of death rate to N
      // % performs element-wise multiplication
      arma::vec d = d0 + dslope % (abundance.t() * interaction).t();
      for (int i = 0; i < abundance.n_elem; i++) if(abundance(i) == 0) d(i) = 0;
      //Gillespie weights for each specie, which are the sum of their rates
      arma::vec w = (abundance % (b + d)) + m; 
      //sampling which species will suffer the next action, proportionally to their weights
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
      // only saves trajectories if we have completed a saving period
      if (((int) (time / save_int)) !=  ((int) ((time+elapsed)/save_int)))
        saveHistory();
      time += elapsed; cycles += 1;
      return;
    }
    void Cbdm(int count = 1) {
      for (int i = 0; i < count; i ++)
        bdm();
    }
    void Tbdm(double _time) {
        while (time < _time) {
            bdm();
        }
    }
};

RCPP_MODULE (Community) {
    using namespace Rcpp;
    class_<Community> ("Community")
        .constructor<arma::vec, arma::mat, double> ()
        .field("abundance", &Community::abundance)
        .field("trajectories", &Community::trajectories)
        .field("interaction", &Community::interaction)
        .field("environmental", &Community::environmental)
        .field("environmental_strength", &Community::environmental_strength)
        .field("K", &Community::K)
        .field("d0", &Community::d0)
        .field("b", &Community::b)
        .field("m", &Community::m)
        .field("time", &Community::time)
        .field("save_int", &Community::save_int)
        .field("cycles", &Community::cycles)
        .method("Cbdm", &Community::Cbdm)
        .method("Tbdm", &Community::Tbdm)
        ;
}
