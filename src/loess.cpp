#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List c_loess(
  NumericVector xx,      // time values - should be 1:n unless there are some NAs
  NumericVector yy,      // the corresponding y values
  int degree,            // degree of smoothing
  int span,              // span of smoothing
  NumericVector ww,      // weights
  IntegerVector m,       // points at which to evaluate the smooth
  IntegerVector l_idx,   // index of left starting points
  NumericVector max_dist // distance between nn bounds for each point
) {
  int span2, span3, offset;
  int i, j;
  double r, tmp1, tmp2;

  int n = xx.size();
  int n_m = m.size();

  NumericVector x(span);
  NumericVector w(span);
  NumericVector xw(span);
  NumericVector x2w(span);
  NumericVector x3w(span);

  NumericVector result(n_m);
  NumericVector slopes(n_m);

  // variables for storing determinant intermediate values
  double a, b, c, d, e, a1, b1, c1, a2, b2, c2, det;

  span3 = span;
  if(span > n) {
    span = n;
  }

  span2 = (span - 1) / 2;

  // want to start storing results at index 0, corresponding to the lowest m
  offset = m[0];

  // loop through all values of m
  for(i = 0; i < n_m; i++) {
    a = 0.0;

    // get weights, x, and a
    for(j = 0; j < span; j++) {
      w[j] = 0.0;
      x[j] = xx[l_idx[i] + j] - (double)m[i];

      // r = std::fabs(x[j]);
      r = (x[j] > 0) ? x[j] : -x[j];
      // tricube
      tmp1 = r / max_dist[i];
      // manual multiplication is much faster than pow()
      tmp2 = 1.0 - tmp1 * tmp1 * tmp1;
      w[j] = tmp2 * tmp2 * tmp2;

      // scale by user-defined weights
      w[j] = w[j] * ww[l_idx[i] + j];

      a = a + w[j];
    }

    if(degree == 0) {
        // TODO: make sure denominator is not 0
       a1 = 1 / a;
       for(j = 0; j < span; j++) {
          // l_i[j] = w[j] * a1;
          result[i] = result[i] + w[j] * a1 * yy[l_idx[i] + j];
       }
    } else {
      // get xw, x2w, b, c for degree 1 or 2
      b = 0.0;
      c = 0.0;
      for(j = 0; j < span; j++) {
        xw[j] = x[j] * w[j];
        x2w[j] = x[j] * xw[j];
        b = b + xw[j];
        c = c + x2w[j];
      }
      if(degree == 1) {
        // TODO: make sure denominator is not 0
        det = 1 / (a * c - b * b);
        a1 = c * det;
        b1 = -b * det;
        c1 = a * det;
        for(j=0; j < span; j++) {
          result[i] = result[i] + (w[j] * a1 + xw[j] * b1) * yy[l_idx[i] + j];
          slopes[i] = slopes[i] + (w[j] * b1 + xw[j] * c1) * yy[l_idx[i] + j];
        }
      } else {
        // TODO: make sure degree > 2 cannot be specified (and < 0 for that matter)
        // get x3w, d, and e for degree 2
        d = 0.0;
        e = 0.0;
        for(j = 0; j < span; j++) {
          x3w[j] = x[j] * x2w[j];
          d = d + x3w[j];
          e = e + x3w[j] * x[j];
        }
        a1 = e * c - d * d;
        b1 = c * d - e * b;
        c1 = b * d - c * c;
        a2 = c * d - e * b;
        b2 = e * a - c * c;
        c2 = b * c - d * a;
        // TODO: make sure denominator is not 0
        det = 1 / (a * a1 + b * b1 + c * c1);
        a1 = a1 * det;
        b1 = b1 * det;
        c1 = c1 * det;
        a2 = a2 * det;
        b2 = b2 * det;
        c2 = c2 * det;
        for(j=0; j < span; j++) {
          result[i] = result[i] + (w[j] * a1 + xw[j] * b1 + x2w[j] * c1) * yy[l_idx[i] + j];
          slopes[i] = slopes[i] + (w[j] * a2 + xw[j] * b2 + x2w[j] * c2) * yy[l_idx[i] + j];
        }
      }
    }
  }

  List ret;
  ret["result"] = result;
  ret["slopes"] = slopes;
  return ret;
}
