
#ifdef _OPENMP
  #include <omp.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <R.h>   
#include <math.h>  
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

SEXP proposed(SEXP X_r, SEXP Y_r, SEXP n_r, SEXP m_r, SEXP p_r, SEXP N_r, SEXP len_r, SEXP e_r, SEXP M_r, SEXP stats_r, SEXP T1_r){
                         
  const int incOne=1; 
  const char ntran='N'; 
  const double one = 1.0; const double zero = 0.0;
  int info; 
  int n = INTEGER(n_r)[0]; int m = INTEGER(m_r)[0]; int p = INTEGER(p_r)[0]; int len = INTEGER(len_r)[0]; int M = INTEGER(M_r)[0]; 
  double *X = REAL(X_r); double *Y = REAL(Y_r); double *N = REAL(N_r);
  double *e = REAL(e_r);
  double *T1 = REAL(T1_r); double *stats = REAL(stats_r);
  int i, q; const int d = 4; 
  int np = n*p, md=m*d; 
  double *T = (double *) S_alloc(len, sizeof(double)); 
  
  GetRNGstate();
  
  #ifdef _OPENMP
   #pragma omp parallel for
   for(i=0;i<len;i++){
   
   double *X0 = (double*) malloc(sizeof(double) * n);
   double *Y0 = (double*) malloc(sizeof(double) * m); 
   double *u = (double*) malloc(sizeof(double) * p);
   double *Z = (double*) malloc(sizeof(double) * m);
   double *Phi1 = (double*) malloc(sizeof(double) * d);
   double counter;
   int j, k;
   
    for(j=0;j<p;j++) u[j] = N[j*len+i];
    F77_NAME(dgemv)(&ntran, &n, &p, &one, X, &n, u, &incOne, &zero, X0, &incOne); // X0 = X*u
    F77_NAME(dgemv)(&ntran, &m, &p, &one, Y, &m, u, &incOne, &zero, Y0, &incOne); // Y0 = Y*u
    
    for(k=0;k<d;k++) Phi1[k] = 0.0; 
    for(j=0;j<m;j++){ 
      counter = 0.0;
      for(k=0;k<n;k++) if(X0[k]<=Y0[j]) counter++;
      Z[j] = counter/n;
      for(k=0;k<d;k++) Phi1[k] += sqrt(2.0)*cos(Z[j]*PI*(k+1));
    }
    
    T[i] = fabs(Phi1[0]);   
    for(k=1;k<d;k++) if(fabs(Phi1[k]) > T[i]) T[i] = fabs(Phi1[k]);
    T[i] /= sqrt(m);
    
    free(X0);
    free(Y0);
    free(u);
    free(Z);
    free(Phi1);
    
    }
  #else
    Rprintf('openMP not used\n')
  #endif
  
  stats[0] = T[0]; 
  for(i=1;i<len;i++) if(T[i] > stats[0]) stats[0] = T[i];
  stats[0] *= sqrt(n/(m+n+0.0));  
  //Rprintf("stats is %4.6f.\n", stats[0]);
  //F77_NAME(dcopy)(&N, Lo, &incOne, Sigma, &incOne);
  
  
#ifdef _OPENMP
  #pragma omp parallel for//parallel //num_threads(4)
  for(q=0;q<M;q++){
   //Rprintf("slave %2d starts\n", q);
   //int me;
   //me = omp_get_thread_num();
   //Rprintf("number of threads: %2d\n", me);
   
   double *X1 = (double*) malloc(sizeof(double) * n);
   double *u1 = (double*) malloc(sizeof(double) * p);
   double *Z1 = (double*) malloc(sizeof(double) * n);
   double *Phi11 = (double*) malloc(sizeof(double) * d);
   double *T0 = (double*) malloc(sizeof(double) * len);
   double counter1;
   int i0, j0, k0;  
   
    // double *e = (double *) S_alloc(n, sizeof(double));
    // for(i=0;i<n;i++) e[i] = rnorm(0.0, 1.0);
    
    for(i0=0;i0<len;i0++){
    for(j0=0;j0<p;j0++) u1[j0] = N[j0*len+i0];
    F77_NAME(dgemv)(&ntran, &n, &p, &one, X, &n, u1, &incOne, &zero, X1, &incOne); // X1 = X*u
    
    for(k0=0;k0<d;k0++) Phi11[k0] = 0.0; 
    for(j0=0;j0<n;j0++){ 
      counter1 = 0.0;
      for(k0=0;k0<n;k0++) if(X1[k0]<=X1[j0]) counter1++;
      Z1[j0] = counter1/n;
      for(k0=0;k0<d;k0++) 
        Phi11[k0] += e[q*n+j0]*sqrt(2.0)*cos(Z1[j0]*PI*(k0+1));
    }
    T0[i0] = fabs(Phi11[0]);   
   for(k0=1;k0<d;k0++) if(fabs(Phi11[k0]) > T0[i0]) T0[i0] = fabs(Phi11[k0]);
    T0[i0] /= sqrt(n);
    
    if(i0==0) T1[q] = T0[i0];
    else if(T0[i0]>T1[q]) T1[q] = T0[i0]; 
    }
    //Rprintf("slave %2d ends\n", q);
    
    free(X1); 
    free(u1); 
    free(Z1);
    free(Phi11);
    free(T0); 
     
   }
#else
  Rprintf('openMP not used\n')
#endif
  
  return(R_NilValue);
 }
