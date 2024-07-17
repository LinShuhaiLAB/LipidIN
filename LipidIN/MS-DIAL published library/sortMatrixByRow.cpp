#include <Rcpp.h>
using namespace Rcpp;


bool compareStrings(const std::string& str1, const std::string& str2) {
  return str1 < str2;
}


CharacterMatrix sortMatrixByRow(CharacterMatrix mat) {
  int rows = mat.nrow();
  int cols = mat.ncol();
  for (int i = 0; i < rows; i++) {
    CharacterVector row = mat(i, _);
    std::vector<std::string> vec(row.begin(), row.end());
    std::sort(vec.begin(), vec.end(), compareStrings);
    for (int j = 0; j < cols; j++) {
      row[j] = vec[j];
    }
    mat(i, _) = row;
  }
  return mat;
}
