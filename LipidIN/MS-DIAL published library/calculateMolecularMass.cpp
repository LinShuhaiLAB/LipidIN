#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculateMolecularMass(CharacterVector formulas) {
  int n = formulas.size();
  NumericVector masses(n);
  
  // Map of element symbols to atomic masses
  std::map<std::string, double> atomicMasses = {
    {"H", 1.007825032},
    {"C", 12},
    {"O", 15.99491462},
    {"N", 14.003074004},
    {"P", 30.973761998},
    {"F", 18.998403163},
    {"S", 31.972071174}
  };
  
  for (int i = 0; i < n; i++) {
    std::string formula = as<std::string>(formulas[i]);
    double mass = 0.0;
    double multiplier = 1.0;
    double digit = 0;
    
    for (size_t j = 0; j < formula.length(); j++) {
      char c = formula[j];
      
      if (std::isalpha(c)) {
        std::string symbol(1, std::toupper(c));
        mass += atomicMasses[symbol] * multiplier;
        digit = atomicMasses[symbol];
        multiplier = 1.0;
      } 
      if (std::isdigit(c)) {
        int multiplier = c - '0';

        char cc = formula[j + 1];
        if (j +1 < formula.length() && std::isdigit(cc)) {
          int multiplierd = formula[j + 1] - '0';
          multiplier = 10.0 * multiplier + multiplierd;
          j++;
        }
        mass += digit * (multiplier - 1);
        multiplier = 1.0;
      }
    }
    
    masses[i] = mass;
  }
  
  return masses;
}
