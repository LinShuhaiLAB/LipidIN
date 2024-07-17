#include <Rcpp.h>
using namespace Rcpp;


inline double absDiff(double x, double y) {
  return fabs(x - y);
}


DataFrame removeRowsWithinError(DataFrame db, NumericVector a, NumericVector mz, double errorThreshold) {
  int n = db.nrows();

  LogicalVector deleteRows(n, false);
  

  for (int i = 0; i < n; ++i) {
    double mzValue = mz[i];
    

    bool withinError = false;
    for (int j = 0; j < a.size(); ++j) {
      if (absDiff(mzValue, a[j])/mzValue <= errorThreshold) {
        withinError = true;
        break;
      }
    }
    

    if (withinError) {
      deleteRows[i] = true;
    }
  }
  

  DataFrame filteredDb = deleteRows;
  
  return filteredDb;
}
