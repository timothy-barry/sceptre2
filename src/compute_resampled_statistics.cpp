#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_resampled_statistics(NumericVector a, NumericVector b, IntegerMatrix synthetic_treatment_idxs) {
  int nrow = synthetic_treatment_idxs.nrow(), ncol = synthetic_treatment_idxs.ncol();
  NumericVector out(ncol);
  double top;
  double bottom;
  int idx;
  for (int j = 0; j < ncol; j++) {
    top = 0;
    bottom = 0;
    for (int i = 0; i < nrow; i++) {
      idx = synthetic_treatment_idxs(i, j);
      top += a[idx];
      bottom += b[idx];
    }
    out[j] = top/sqrt(bottom);
  }
  return out;
}
