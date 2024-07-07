#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
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
// [[Rcpp::export]]
NumericVector propscale(const IntegerVector& labels, const NumericVector& values) {
  int n = labels.size();
  NumericVector normalizedValues(n);
  
  // 在labels向量中查找唯一标签
  IntegerVector uniqueLabels = Rcpp::unique(labels);
  int numLabels = uniqueLabels.size();
  
  // 对每个唯一标签进行归一化
  for(int i = 0; i < numLabels; i++) {
    int label = uniqueLabels[i];
    
    // 找到属于当前标签的所有数值
    NumericVector currentValues;
    for(int j = 0; j < n; j++) {
      if(labels[j] == label) {
        currentValues.push_back(values[j]);
      }
    }
    
    // 计算当前数值的最小值和最大值
    double minValue = Rcpp::min(currentValues);
    double maxValue = Rcpp::max(currentValues);
    
    // 归一化当前标签下的数值
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

// [[Rcpp::export]]
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




