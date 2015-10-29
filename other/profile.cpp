////////////// CODE FOR PROFILING SRC //////////////
#include "../src/main.cpp"

int main () {
  int N = 100;
  arma::vec abundance (N, arma::fill::ones);
  arma::mat interaction = arma::randu(N, N);
  arma::vec K (N); K.fill(100);
  arma::vec d0(N, arma::fill::zeros);
  arma::vec b(N, arma::fill::ones);
  arma::vec m(N); m.fill(0.05);
  create_community(abundance, interaction, K, d0, b, m);
  bdm(1000000);
  std::cout << time() << std::endl;
  return 0;
}


