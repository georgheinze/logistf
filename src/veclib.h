#ifndef ___VECLIB_H
#define ___VECLIB_H

#define USE_FC_LEN_T

#include <math.h>						// powf, ...
#include <R.h>
#include <Rdefines.h>				
#include "memory.h"					// malloc; free 
//#include <R_ext/RS.h>       // Needed for F77_CALL macro 
#include <R_ext/Lapack.h>	// inverse; choleski; determinant
#include "Rmath.h"					// random numbers; distributions

#ifndef FCONE
# define FCONE
#endif

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
// cross-product of m x k matrix and k x n matrix (XY)  ; result is m x n
void XY(double *X, double *Y, double *res,long k, long m, long n)
{
	long i, j, ind;
	double tmp;
	
	for(i=0; i < m; i++)
		for(j=0; j < n; j++) {
			tmp = 0.0;
			for(ind = 0; ind < m; ind++)
				tmp += X[i + ind*m] * Y[ind + j*k];
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


void lapack_det(double *A_doub, long *size, double *logdet)
{
    int i, j;
    int n = (int) *size;
    int lda = (int) *size;
    int N = n * lda;
    int info;
    char uplo = 'U';  // Use upper triangle
    
    double *A = (double *) R_alloc(N, sizeof(double));
    if (A == NULL) error("no memory available\n");
    
    // Copy matrix to Fortran column-major format
    for (i = 0; i < n; i++)
        for (j = 0; j < lda; j++)
            A[j + n * i] = A_doub[j + i * n];
    
    
    // Cholesky factorization
    F77_CALL(dpotrf)(&uplo, &n, A, &lda, &info FCONE);
    if (info != 0)
        error("LAPACK dpotrf failed: matrix is not positive definite (info = %d)", info);
    
    // Compute log determinant = 2 * sum(log(diagonal elements))
    double sum_log_diag = 0.0;
    for (i = 0; i < n; i++) {
        double diag_val = A[i + i * n];
        if (diag_val <= 0.0)
            error("Non-positive diagonal in Cholesky factor.");
        sum_log_diag += log(diag_val);
    }
    
    *logdet = 2.0 * sum_log_diag;
}

void lapack_inv(double *A_doub, long *size)
{
    int i, j, n, N, info;
    double *A;
    char uplo = 'U';  // LAPACK stores upper or lower Cholesky
    
    n = (int) *size;
    N = n * n;
    
    if ((A = (double *) R_alloc(N, sizeof(double))) == NULL)
        error("no memory available\n");
    
    // Copy to column-major layout
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A[j + n * i] = A_doub[j + i * n];
    
    // Step 1: Cholesky decomposition
    F77_CALL(dpotrf)(&uplo, &n, A, &n, &info FCONE);
    if (info != 0)
        error("LAPACK dpotrf failed: matrix not positive definite (info = %d)", info);
    
    // Step 2: Compute inverse using Cholesky factor
    F77_CALL(dpotri)(&uplo, &n, A, &n, &info FCONE);
    if (info != 0)
        error("LAPACK dpotri failed: inversion failed (info = %d)", info);
    
    // Copy upper triangle to full matrix (symmetric)
    for (i = 0; i < n; i++) {
        A_doub[(n + 1) * i] = A[(n + 1) * i];  // Diagonal
        for (j = 0; j < i; j++) {
            double val = A[j + i * n];           // Upper triangle value
            A_doub[i + j * n] = A_doub[j + i * n] = val;  // Symmetric
        }
    }
}




void testRmath(void)
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
