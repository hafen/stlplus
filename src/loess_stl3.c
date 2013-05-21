#include <R.h>
#include <Rmath.h>

double loess_stl2( double *xx,
                   double *yy,
                   int *nn, 
                   int *ddegree, 
                   int *sspan, 
                   double *ww,
                   int *m, 
                   int *nn_m, 
                   int *ll_idx,
                   int *rr_idx,
                   double *mmax_dist,
                   double *result,
                   double *slopes
                 );
inline double my_abs(double x);

double loess_stl2( double *xx,      // time values - should be 1:n unless there are some NAs
                   double *yy,      // the corresponding y values
                   int *nn,         // the length of x
                   int *ddegree,    // degree of smoothing
                   int *sspan,      // span of smoothing
                   double *ww,
                   int *m,          // points at which to evaluate the smooth
                   int *nn_m,       // number of points at which to evaluate the smooth
                   int *ll_idx,
                   int *rr_idx,
                   double *mmax_dist,
                   double *result,  // vector of length nn_m
                   double *slopes
) {
   int n, n_m, degree, span, span2, span3, idx, offset, ll, rr, ok;
   int i, j, k, ii, jj;
   double *x, *w, *xw, *x2w, *x3w, max_dist, r, tmp1, tmp2;
   int l_idx, r_idx, tmp;

   n = *(nn);
   n_m = *(nn_m);
   degree = *(ddegree);
   span =  *(sspan);

   x = (double *) R_alloc(span, sizeof(double));
   w = (double *) R_alloc(span, sizeof(double));
   xw = (double *) R_alloc(span, sizeof(double));
   x2w = (double *) R_alloc(span, sizeof(double));
   x3w = (double *) R_alloc(span, sizeof(double));
   
   // variables for storing determinant intermediate values
   double a, b, c, d, e, a1, b1, c1, a2, b2, c2, det;
   
   span3 = span;
   if(span > n) span = n;
   span2 = (span - 1)/2;

   // want to start storing results at index 0, corresponding to the lowest m
   offset = *(m);
   
   // Rprintf("%d\n", n_m);
   
   // loop through all values of m
   for(i = 0; i < n_m; i++) {

      l_idx = *(ll_idx + i);
      r_idx = *(rr_idx + i);
      max_dist = *(mmax_dist + i);

      // Rprintf("%d, %d, %d, %f, %f, %f\n", n, span, *(m + i), *(xx + l_idx), *(xx + r_idx), max_dist);
      a = 0.0;
      
      // get weights, x, and a
      for(j = 0; j < span; j++) {
         *(w + j) = 0.0;
         *(x + j) = *(xx + l_idx + j) - (double)*(m + i);
      
         r = my_abs(*(x + j));
         // tricube - hard coded instead of using pow()
         tmp1 = r/max_dist;
         tmp2 = 1.0 - tmp1*tmp1*tmp1;
         *(w + j) = tmp2 * tmp2 * tmp2;
         
         // scale by user-defined weights
         *(w + j) = *(w + j) * *(ww + l_idx + j);                  

         a = a + *(w + j);
      }

      if(degree==0) {
         //TODO: make sure a is not zero
         a1 = 1/a;
         for(j=0; j < span; j++) {
            // *(l_i + j) = *(w + j) * a1;
            *(result + i) = *(result + i) + *(w + j) * a1 * *(yy + l_idx + j);
         }
      } else {
      
         // get xw, x2w, b, c for degree 1 or 2
         b = 0.0;
         c = 0.0;
         for(j = 0; j < span; j++) {
            *(xw + j) = *(x + j) * *(w + j);
            *(x2w + j) = *(x + j) * *(xw + j);
            b = b + *(xw + j);
            c = c + *(x2w + j);
         }
         if(degree==1) {
            // TODO: make sure this is 0
            det = 1/(a * c - b * b);
            a1 = c * det;
            b1 = -b * det;
            c1 = a * det;
            for(j=0; j < span; j++) {
               *(result + i) = *(result + i) + (*(w + j) * a1 + *(xw + j) * b1) * *(yy + l_idx + j);
               *(slopes + i) = *(slopes + i) + (*(w + j) * b1 + *(xw + j) * c1) * *(yy + l_idx + j);
            }
         } else {
            // TODO: make sure degree > 2 cannot be specified (and < 0 for that matter)
            // get x3w, d, and e for degree 2
            d = 0.0;
            e = 0.0;
            for(j = 0; j < span; j++) {
               *(x3w + j) = *(x + j) * *(x2w + j);
               d = d + *(x3w + j);
               e = e + *(x3w + j) * *(x + j);
            }
            a1 = e * c - d * d;
            b1 = c * d - e * b;
            c1 = b * d - c * c;
            a2 = c * d - e * b;
            b2 = e * a - c * c;
            c2 = b * c - d * a;
            // TODO: make sure this is 0
            det = 1/(a*a1 + b*b1 + c*c1);
            a1 = a1 * det;
            b1 = b1 * det;
            c1 = c1 * det;
            a2 = a2 * det;
            b2 = b2 * det;
            c2 = c2 * det;
            for(j=0; j < span; j++) {
               *(result + i) = *(result + i) + (*(w + j) * a1 + *(xw + j) * b1 + *(x2w + j) * c1) * *(yy + l_idx + j);
               *(slopes + i) = *(slopes + i) + (*(w + j) * a2 + *(xw + j) * b2 + *(x2w + j) * c2) * *(yy + l_idx + j);
            }
         }
      }
   }
}      

inline double my_abs(double x)
{
	return (x > 0) ? x : -x;
}



// old code for nearest neighbors
// 
// l_idx = r_idx = *(m + i) - 1;
// if(l_idx < 0) {
//    l_idx = 0; r_idx=0;
// }
// if(r_idx > n-1) {
//    r_idx = n - 1; l_idx = n - 1;         
// }
//    
// for(j = 1; j < span; j++)
// {
//    if(l_idx == 0) {
//       r_idx = r_idx + 1;
//    }
//    else if(r_idx == n-1) {
//       l_idx = l_idx - 1;
//    }
//    else {
//       tmp1 = *(xx + l_idx - 1) - (double)*(m + i);
//       tmp2 = *(xx + r_idx + 1) - (double)*(m + i);
//       if(my_abs(tmp1) > my_abs(tmp2)) {
//          r_idx = r_idx + 1;
//       }
//       else {
//          l_idx = l_idx - 1;
//       }
//       // Rprintf("%f, %f, %d, %d, %f, %f\n", tmp1, tmp2, l_idx, r_idx, *(xx + l_idx), *(xx + r_idx));
//    }
// }
// 
// tmp1 = (double)*(m + i) - *(xx + l_idx);
// tmp2 = *(xx + r_idx) - (double)*(m + i);
// 
// if(tmp1 > tmp2)
//    max_dist = tmp1;
// else
//    max_dist = tmp2;
// 
