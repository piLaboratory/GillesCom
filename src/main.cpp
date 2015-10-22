#include<iostream>
#include<fstream>
#include<string>
#include<Rcpp.h>

class Community {
  private:
    std::vector<int> abundance;
  public:
    std::vector<int> const get_abundance() {return abundance;}
    int get_nspecies() {return abundance.size();}
    Community() {
      
    }
};

// Global var???
double x = 0;

// [[Rcpp::export]]
double test_basic(){
  x++;

	return x;
}
