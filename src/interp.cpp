#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector c_interp(NumericVector m, NumericVector fits, NumericVector slopes, NumericVector at) {

  int i, j;
  double u, h, u2, u3;

  int n_at = at.size();
  NumericVector ans(n_at);

  j = 0; // index of leftmost vertex
  for(i = 0; i < n_at; ++i) {
    if(at[i] > m[j + 1]) j++;
    h = (m[j + 1] - m[j]);
    u = (at[i] -  m[j]) / h;
    u2 = u * u;
    u3 = u2 * u;
    ans[i] = (2 * u3 - 3 * u2 + 1) * fits[j] +
             (3 * u2 - 2 * u3)     * fits[j + 1] +
             (u3 - 2 * u2 + u)     * slopes[j] * h +
             (u3 - u2)             * slopes[j + 1] * h;
  }

  return ans;
}

/*** R
if(any(is.nan(fits))) {
  ind <- !is.nan(fits)
  c_interp(m[ind], fits[ind], slopes[ind], at)
} else {
  c_interp(m, fits, slopes, at)
}
*/
