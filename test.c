
// test if openMP is running well

#ifdef _OPENMP
  #include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include "util.h"

SEXP test(SEXP Sigma_r, SEXP n_r){

  const int incOne = 1;  int info; int k; 
  int n = INTEGER(n_r)[0]; double *Sigma = REAL(Sigma_r); 
  double a = 0.0; int N = n*n;

  GetRNGstate();

#ifdef _OPENMP
  #pragma omp parallel //num_threads(3)
  {
   #pragma omp for //firstprivate(i) lastprivate(k) 
   for(k=0;k<4;k++){
    int N = n*n; double *Lo = (double *) S_alloc(N, sizeof(double));
    F77_NAME(dcopy)(&N, Sigma, &incOne, Lo, &incOne);
    for(int i = 0; i < n; i++) Lo[i*n+i] += rnorm(10.0, 1.0);
    Rprintf("%4.6f.\n", Lo[0]);
    mtrxInv(Lo, n, info);
    Rprintf("%4.6f...\n", Lo[0]);
    a += Lo[0];
  }
  }
#else
 Rprintf('not used\n')
#endif
  
  Rprintf("sum is %4.6f.\n", a);
  //cout<<a<<endl;
  //F77_NAME(dcopy)(&N, Lo, &incOne, Sigma, &incOne);
  return(R_NilValue);
 }
