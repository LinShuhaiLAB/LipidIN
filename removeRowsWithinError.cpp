#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate absolute difference between two values
inline double absDiff(double x, double y) {
  return fabs(x - y);
}

// [[Rcpp::export]]
DataFrame removeRowsWithinError(DataFrame db, NumericVector a, NumericVector mz, double errorThreshold) {
  int n = db.nrows();
  
  // Create a logical vector to mark rows for deletion
  LogicalVector deleteRows(n, false);
  
  // Iterate over each row in the database
  for (int i = 0; i < n; ++i) {
    double mzValue = mz[i];
    
    // Check if the mz value is within the error threshold of any value in 'a'
    bool withinError = false;
    for (int j = 0; j < a.size(); ++j) {
      if (absDiff(mzValue, a[j])/mzValue <= errorThreshold) {
        withinError = true;
        break;
      }
    }
    
    // Mark row for deletion if mz value is within error threshold
    if (withinError) {
      deleteRows[i] = true;
    }
  }
  
  // Filter the database using the logical vector
  DataFrame filteredDb = deleteRows;
  
  return filteredDb;
}
