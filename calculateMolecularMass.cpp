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
      
      if (std::isalpha(c)) {// 判断是否字符
        std::string symbol(1, std::toupper(c));// 将其转大写
        mass += atomicMasses[symbol] * multiplier;
        digit = atomicMasses[symbol];
        multiplier = 1.0;
      } 
      if (std::isdigit(c)) {// 判断是否数字
        int multiplier = c - '0';// 转换为数字
        // 判断下一个是不是也是数字
        char cc = formula[j + 1];
        if (j +1 < formula.length() && std::isdigit(cc)) {
          int multiplierd = formula[j + 1] - '0';// 转换为数字,个位
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