#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector block_minority_prop_cpp(IntegerVector house_block,
                                      IntegerVector resident_house,
                                      IntegerVector resident_minority,
                                      int n_block) {
  IntegerVector house_minority(house_block.size(), 0);
  
  int n_res = resident_house.size();
  for (int i = 0; i < n_res; i++) {
    int h = resident_house[i] - 1;
    if (h >= 0 && h < house_minority.size()) {
      house_minority[h] = resident_minority[i];
    }
  }
  
  NumericVector count(n_block, 0.0);
  NumericVector minority_count(n_block, 0.0);
  int n_house = house_block.size();
  for (int i = 0; i < n_house; i++) {
    if (IntegerVector::is_na(house_block[i])) continue;  // ← NAをスキップ
    int b = house_block[i] - 1;
    if (b < 0 || b >= n_block) continue;                 // ← 範囲外も防御
    count[b] += 1;
    if (house_minority[i] == 1) minority_count[b] += 1;
  }
  
  NumericVector prop(n_block);
  for (int i = 0; i < n_block; i++) {
    prop[i] = (count[i] > 0) ? minority_count[i] / count[i] : NA_REAL;
  }
  return prop;
}