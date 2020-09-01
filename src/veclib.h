#ifndef ___VECLIB_H
#define ___VECLIB_H

#include <math.h>						// powf, ...
#include <R.h>
#include <Rdefines.h>				
#include "memory.h"					// malloc; free 
#include <R_ext/Linpack.h>	// inverse; choleski; determinant
#include "Rmath.h"					// random numbers; distributions


// fast copy of array X to array res, type double
void copy(double *X, double *res, long n)
{
	memcpy(res, X, n * sizeof(double));
}

// maximum absolute value
double maxabs(double *X, long n)
{
	double res = fabs(X[0]);
	for(long i=1; i < n; i++)
		res = fmax(res, fabs(X[i]));
	return res;
}

// maximum absolute value of selected indices
double maxabsInds(double *X, int *inds, long n_inds)
{
	double res = fabs(X[inds[0]]);
	for(long i=1; i < n_inds; i++)
		res = fmax(res, fabs(X[inds[i]]));
	return res;
}

void print(double *X, long k, long m)
{
	Rprintf("%ld x %ld matrix:\n", k, m);
	for(long i=0; i < k; i++) {
		for(long j=0; j < m; j++)
			Rprintf(" %8.3f", X[j*k + i]);
		Rprintf("\n");
	}
}

// copy some columns of a matrix with n rows -- not tested yet!!
void copyCols(double *X, double *res, long n, long *inds, long n_inds)
{
	for(long i=0; i < n_inds; i++)
		memcpy(res + n * sizeof(double) * i, 
					 X + n * sizeof(double) * inds[i],
					 n * sizeof(double));
}

// product of k x m matrix and m x k matrix (XY); only diagonal! (k x 1)
void XYdiag(double *X, double *Y, double *res, long k, long m)
{
	double tmp;
	for(long i=0; i < k; i++) {
		tmp = 0.0;
		for(long ind = 0; ind < m; ind++)
			tmp += X[ind * k + i] * Y[ind + i * m];
		res[i] = tmp;
	}
}


// cross-product of k x m matrix and k x n matrix (X'Y)  ; result is m x n
void XtY(double *X, double *Y, double *res, long k, long m, long n)
{
	long i, j, ind;
	double tmp;
	
	for(i=0; i < m; i++)
		for(j=0; j < n; j++) {
			tmp = 0.0;
			for(ind = 0; ind < k; ind++)
				tmp += X[ind + i*k] * Y[ind + j*k];
			res[i + j*m] = tmp;
		}
}

// X'X (X is k x k)
void XtXsym(double *X, double *res, long *k_l)
{
	long k = (long)*k_l;
	long i, j, ind;
	double tmp;
	
	for(i=0; i < k; i++)
		for(j=i; j < k; j++) {
			tmp = 0.0;
			for(ind = 0; ind < k; ind++)
				tmp += X[ind + i*k] * X[ind + j*k];
			res[i + j*k] = res[j + i*k] = tmp;
		}
}

// X'X (X is k x m)  ;  result is m x m
void XtXasy(double *X, double *res, long k, long m)
{
	long i, j, ind;
	double tmp;
	
	for(i=0; i < m; i++)
		for(j=i; j < m; j++) {
			tmp = 0.0;
			for(ind = 0; ind < k; ind++)
				tmp += X[ind + i*k] * X[ind + j*k];
			res[i + j*m] = res[j + i*m] = tmp;
		}
}

// t(X) (X is k x m)
void trans(double *X, double *res, long k, long m)
{
	for(long i=0; i < k; i++)
		for(long j=0; j < m; j++) {
			res[i*m + j] = X[i + j*k];
		  //Rprintf("trans: %f", res[i*m + j]);
		}
}


// compute inverse and determinant; A_doub is changed
void linpack_inv_det(double *A_doub, long *size, double *logdet)
{
  int i, j , c1, c2, ok, N, n;
  double *A, *det;
  n = (long) *size;                        // type-cast  
  N = n * n; 
  if (NULL == (A = (double *) R_alloc(N, sizeof(double))))
  {
	 error("no memory available\n");
  }
  if (NULL == (det = (double *) R_alloc(2, sizeof(double))))
  {
	 error("no memory available\n");
  }
	
  /* transpose matrix for fortran (reverse order) */
  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      A[j + n * i] = A_doub[j + i * n];
  c1 = n;
  c2 = n;
  
  /* factorize a symmetric matrix */
  F77_NAME(dpofa)(A, &c1, &c2, &ok);
  //Rprintf("ok=%d\n", ok);                                                                                                           
  
  /* compute the determinant and inverse of a certain
	 double precision symmetric positive definite matrix
	 as result A is upper half of inverse */
  ok = 11; // calc both determ+inverse
  F77_NAME(dpodi)(A, &c1, &c2, det, &ok);
  //Rprintf("ok=%d\n", ok);  
  
  // copy result                                                       
  for (i=0; i < n; i++) {
    A_doub[(n+1) * i] = A[(n+1) * i];            // diagonal elements  
    for(j=0; j < i; j++)
      A_doub[i + j*n] = A_doub[j + i*n] = A[i*n + j];
  }
	
  *logdet = log(det[0]) +  M_LN10 * det[1];
}

// compute determinant; A_doub is unchanged
void linpack_det(double *A_doub, long *size, double *logdet)
{
  int i, j , c1, c2, ok, N, n;
  double *A, *det;
  
  n = (long) *size;                        // type-cast                                                                               
  N = n * n; 
  if(NULL == (A = (double *) R_alloc(N, sizeof(double))))
  {
	 error("no memory available\n");
  }
  if(NULL == (det = (double *) R_alloc(2, sizeof(double))))
  {
	 error("no memory available\n");
  }
	
  /* transpose matrix for fortran (reverse order) */
  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      A[j + n * i] = A_doub[j + i * n];
  c1 = n;
  c2 = n;
  
  /* factorize a symmetric matrix */
  F77_NAME(dpofa)(A, &c1, &c2, &ok);
  //Rprintf("ok=%d\n", ok);                                                                                                           
  
  /* compute the determinant and inverse of a certain
	 double precision symmetric positive definite matrix
	 as result A is upper half of inverse */
  ok = 10; // calc determ only
  F77_NAME(dpodi)(A, &c1, &c2, det, &ok);
  //Rprintf("ok=%d\n", ok);  
	
  *logdet = log(det[0]) +  M_LN10 * det[1];
}

// compute inverse; A_doub is changed
void linpack_inv(double *A_doub, long *size)
{
  int i, j , c1, c2, ok, N, n;
  double *A, *det;
  //char uplo = 'L';                // Lower or Upper tri                                                                  
	
  n = (long) *size;                        // type-cast                                                                               
  N = n * n; 
  if(NULL == (A = (double *) R_alloc(N, sizeof(double))))
  {
	 error("no memory available\n");
  }
  if(NULL == (det = (double *) R_alloc(2, sizeof(double))))
  {
	 error("no memory available\n");
  }
	
  /* transpose matrix for fortran (reverse order) */
  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      A[j + n * i] = A_doub[j + i * n];
  c1 = n;
  c2 = n;
  
  /* factorize a symmetric matrix */
  F77_NAME(dpofa)(A, &c1, &c2, &ok);
  //Rprintf("ok=%d\n", ok);                                                                                                           
  
  /* compute the determinant and inverse of a certain
	 double precision symmetric positive definite matrix
	 as result A is upper half of inverse */
  ok = 01; // calc inverse only
  F77_NAME(dpodi)(A, &c1, &c2, det, &ok);
  //Rprintf("ok=%d\n", ok);  
  
  // copy result                                                       
  for (i=0; i < n; i++) {
    A_doub[(n+1) * i] = A[(n+1) * i];            // diagonal elements  
    for(j=0; j < i; j++)
      A_doub[i + j*n] = A_doub[j + i*n] = A[i*n + j];
  }

}

void linpack_choleski(double *A_doub, long *size)
{
  int i, j , c1, c2, job, ok, N, n, *jpvt;
  double *A, *work;
  //char uplo = 'L';                // Lower or Upper tri                                                                  
	
  n = (int) *size;                        // type-cast                                                                               
  N = n * n; 
  if (NULL == (A = (double *) R_alloc(N, sizeof(double))))
  {
	 error("no memory available\n");
  }
  if(NULL == (work = (double *) R_alloc(N, sizeof(double))))
  {
	 error("no memory available\n");
  }
  if (NULL == (jpvt = (int *) R_alloc(n, sizeof(int))))
  {
	 error("no memory available\n");
  }
	
  /* transpose matrix for fortran (reverse order) */
  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      A[j + n * i] = A_doub[j + i * n];
  c1 = n;
  c2 = n;
	job = 0; // if job .eq. 0, no pivoting is done; else yes !!??
  
  /* dchdc computes the cholesky decomposition of a positive definite
	 matrix.  a pivoting option allows the user to estimate the
	 condition of a positive definite matrix or determine the rank
	 of a positive semidefinite matrix. */
  F77_NAME(dchdc)(A, &c1, &c2, work, jpvt, &job, &ok);
  //Rprintf("ok=%d\n", ok);                                                                                                           
	//for(i=0; i < N; i++) Rprintf("%d \t %e \n", i, A[i]);
  
  // copy result                                                       
  for (i=0; i < n; i++) {
    A_doub[(n+1) * i] = A[(n+1) * i];            // diagonal elements  
    for(j=0; j < i; j++) {
			A_doub[j + n * i] = A[j + i * n];		// upper element
			A_doub[i + n * j] = 0.0;						// lower element
		}
  }
	

}


void testRmath()
{
	double res;
	res = R_pow(3.0, 2.0);
	Rprintf("3^2 is %f. \n", res);
	
	// dnorm : x, mean = 0, sd = 1, log = FALSE
	res = dnorm(1.0, 0.0, 1.0, 0);
	Rprintf("normal density on 1 is %f. \n", res);
}

void summe(double *x, long *n, double *res)
{	
	long i;
	*res = 0 ;
	for (i = 0; i < *n; i++)
		*res += x[i];
}


#endif
