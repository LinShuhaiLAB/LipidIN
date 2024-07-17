#include <Rcpp.h>
using namespace Rcpp;

IntegerVector rtcluster(const IntegerVector& X, const NumericVector& rt, double error) {
  int n = rt.size();
  IntegerVector clusters(n);
  int clusterCount = 0;
  
  for(int i = 0; i < n; i++) {
    double currentMz = rt[i];
    if(clusters[i] == 0) {
      clusters[i] = ++clusterCount;
    }
    
    for(int j = i + 1; j < n; j++) {
      if(clusters[j] == 0) {
        double diff = std::abs(currentMz - rt[j]);
        if(diff <= error) {
          clusters[j] = clusters[i];
        }
      }
    }
  }
  return clusters;
}



#include <Rcpp.h>
using namespace Rcpp;

NumericVector propscale(const IntegerVector& labels, const NumericVector& values) {
  int n = labels.size();
  NumericVector normalizedValues(n);
  

  IntegerVector uniqueLabels = Rcpp::unique(labels);
  int numLabels = uniqueLabels.size();
  

  for(int i = 0; i < numLabels; i++) {
    int label = uniqueLabels[i];

    NumericVector currentValues;
    for(int j = 0; j < n; j++) {
      if(labels[j] == label) {
        currentValues.push_back(values[j]);
      }
    }
    

    double minValue = Rcpp::min(currentValues);
    double maxValue = Rcpp::max(currentValues);

    for(int j = 0; j < n; j++) {
      if(labels[j] == label) {
        double normalizedValue = (values[j] - minValue) / (maxValue - minValue);
        normalizedValues[j] = normalizedValue;
      }
    }
  }
  
  return normalizedValues;
}

#include <Rcpp.h>
using namespace Rcpp;


IntegerVector rawmzcluster(const IntegerVector& X, const NumericVector& rawmz, double error) {
  int n = rawmz.size();
  IntegerVector clusters(n);
  int clusterCount = 0;
  
  for(int i = 0; i < n; i++) {
    double currentMz = rawmz[i];
    if(clusters[i] == 0) {
      clusters[i] = ++clusterCount;
    }
    
    for(int j = i + 1; j < n; j++) {
      if(clusters[j] == 0) {
        double diff = std::abs(currentMz - rawmz[j])/currentMz;
        if(diff <= error) {
          clusters[j] = clusters[i];
        }
      }
    }
  }
  return clusters;
}




