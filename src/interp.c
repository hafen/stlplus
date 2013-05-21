# // x_l is location of fit at nodes
# // f_x is fit at nodes
# // f_x_p is slope at fit at nodes
# // nn is length of f_x
# // x is where to compute
# // y is where to store fits
# // NN is length of x
# 

void stl2_interp(int *x_l, double *f_x, double *f_x_p, int *nn, int *x, int *NN, double *y);

void stl2_interp(int *x_l, double *f_x, double *f_x_p, int *nn, int *x, int *NN, double *y)
{
   int i, j, n, N;
   double u, h;

   n = *nn;
   N = *NN;

   j = 0; // index of leftmost vertex
   for(i = 0; i < N; ++i)
   {
      if(*(x+i) > *(x_l+j+1)) j++;
      h = (*(x_l+j+1) - *(x_l+j));
      u = (*(x+i) -  *(x_l+j)) / h;
      *(y+i) = (2*u*u*u - 3*u*u + 1) * *(f_x + j) + 
               (3*u*u - 2*u*u*u)     * *(f_x + j + 1) + 
               (u*u*u - 2*u*u + u)   * *(f_x_p + j)*h + 
               (u*u*u - u*u)         * *(f_x_p + j + 1)*h;
   }
}