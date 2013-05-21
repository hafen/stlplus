#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void ma(double *x, int *n, int *n_p, double *ans) 
{   
   int i;

   int nn = *n - *n_p * 2;
   int nn_p = *n_p;

   double *ma, *ma2, *ma3;
   double ma_tmp;
   
   ma = (double *) malloc((nn + nn_p + 1)* sizeof(double));
   ma2 = (double *) malloc((nn + 2)* sizeof(double));
   ma3 = (double *) malloc(nn * sizeof(double));

   ma_tmp = 0;
   for(i = 0; i < nn_p; ++i)
   {
      ma_tmp = ma_tmp + *(x + i);
   }
   *(ma) = ma_tmp / nn_p;
   for(i = nn_p; i < nn + 2*nn_p; ++i)
   {
      ma_tmp = ma_tmp - *(x + i - nn_p) + *(x + i);
      *(ma + i - nn_p + 1) = ma_tmp / nn_p; 
   }

   ma_tmp = 0;
   for(i = 0; i < nn_p; ++i)
   {
      ma_tmp = ma_tmp + *(ma + i);
   }
   *(ma2) = ma_tmp / nn_p;
   for(i = nn_p; i < nn + nn_p + 1; ++i)
   {
      ma_tmp = ma_tmp - *(ma + i - nn_p) + *(ma + i);
      *(ma2 + i - nn_p + 1) = ma_tmp / nn_p; 
   }

   ma_tmp = 0;
   for(i = 0; i < 3; ++i)
   {
      ma_tmp = ma_tmp + *(ma2 + i);
   }
   *(ans) = ma_tmp / 3;
   for(i = 3; i < nn + 2; ++i)
   {
      ma_tmp = ma_tmp - *(ma2 + i - 3) + *(ma2 + i);
      *(ans + i - 2) = ma_tmp / 3; 
   }
   free(ma);
   free(ma2);
   free(ma3);
}
