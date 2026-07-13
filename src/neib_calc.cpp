#include <Rcpp.h>
using namespace Rcpp;

// -------------------------------------------------------------
// neib_ethnicity_prop
// 各セルの近傍における majority/minority 割合を計算する。
//   neib_ind      : 近傍インデックス行列 (n_cell x max_neib)。
//                   範囲外は -1（0-indexed済み）で渡す。
//   house_minority: 各セルの minority フラグ (0/1)。空室は NA。
// 戻り値: n_cell x 2 行列（列0=majority割合, 列1=minority割合）
// -------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix neib_ethnicity_prop(IntegerMatrix neib_ind,
                                  NumericVector house_minority) {
  int n_cell   = neib_ind.nrow();
  int max_neib = neib_ind.ncol();
  NumericMatrix out(n_cell, 2);
  
  for (int i = 0; i < n_cell; i++) {
    double sum_mino = 0.0;
    int    count    = 0;
    for (int k = 0; k < max_neib; k++) {
      int idx = neib_ind(i, k);
      if (idx < 0) continue;                     // 範囲外セル
      double val = house_minority[idx];
      if (NumericVector::is_na(val)) continue;   // 空室
      sum_mino += val;
      count++;
    }
    if (count > 0) {
      double mino_prop = sum_mino / count;
      out(i, 0) = 1.0 - mino_prop;
      out(i, 1) = mino_prop;
    }
    // count == 0 の場合は 0 のまま（元のNaN->0処理と同じ挙動）
  }
  return out;
}

// -------------------------------------------------------------
// neib_SES_prop
// 各セルの近傍における SES(0-5) 別割合を計算する。
// -------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix neib_SES_prop(IntegerMatrix neib_ind,
                            NumericVector house_SES) {
  int n_cell   = neib_ind.nrow();
  int max_neib = neib_ind.ncol();
  NumericMatrix out(n_cell, 6);
  
  for (int i = 0; i < n_cell; i++) {
    int counts[6] = {0, 0, 0, 0, 0, 0};
    int total = 0;
    for (int k = 0; k < max_neib; k++) {
      int idx = neib_ind(i, k);
      if (idx < 0) continue;
      double val = house_SES[idx];
      if (NumericVector::is_na(val)) continue;
      int ses = (int)val;
      if (ses >= 0 && ses <= 5) { counts[ses]++; total++; }
    }
    if (total > 0) {
      for (int j = 0; j < 6; j++) out(i, j) = (double)counts[j] / total;
    }
  }
  return out;
}