#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_ma(NumericVector x, int n_p) {
  int i;
  int n = x.size();
  int nn = n - n_p * 2;
  int nn_p = n_p;
  double ma_tmp = 0;

  NumericVector ans(n - 2 * nn_p);
  NumericVector ma(nn + nn_p + 1);
  NumericVector ma2(nn + 2);
  NumericVector ma3(nn);

  ma_tmp = 0;
  for(i = 0; i < nn_p; ++i) {
    ma_tmp = ma_tmp + x[i];
  }
  ma[0] = ma_tmp / nn_p;
  for(i = nn_p; i < nn + 2 * nn_p; ++i) {
    ma_tmp = ma_tmp - x[i - nn_p] + x[i];
    ma[i - nn_p + 1] = ma_tmp / nn_p;
  }

  ma_tmp = 0;
  for(i = 0; i < nn_p; ++i) {
    ma_tmp = ma_tmp + ma[i];
  }
  ma2[0] = ma_tmp / nn_p;

  for(i = nn_p; i < nn + nn_p + 1; ++i) {
    ma_tmp = ma_tmp - ma[i - nn_p] + ma[i];
    ma2[i - nn_p + 1] = ma_tmp / nn_p;
  }

  ma_tmp = 0;

  for(i = 0; i < 3; ++i) {
    ma_tmp = ma_tmp + ma2[i];
  }
  ans[0] = ma_tmp / 3;

  for(i = 3; i < nn + 2; ++i) {
    ma_tmp = ma_tmp - ma2[i - 3] + ma2[i];
    ans[i - 2] = ma_tmp / 3;
  }

  return ans;
}
