#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "memory.h"
#include "Rmath.h"
#include "veclib.h"

void logistffit_revised(double *x, int *y, int *n_l, int *k_l, 
                double *weight, double *offset, 
                double *beta,
                int *colfit, int *ncolfit_l, int *firth_l, 
                int *maxit, double *maxstep, int *maxhs,
                double *lconv, double *gconv, double *xconv, double* tau, 
                // output: 
                double *fisher_cov,		// k x k
                double *Ustar,				// k
                double *pi,						// n
                double *Hdiag,				// n
                double *loglik,				// 1
                int *evals,
                int *iter,
                double *convergence		// 3
)
{
  long n = (long)*n_l, k = (long)*k_l, firth = (long)*firth_l, ncolfit = (long)*ncolfit_l;
  double wi, logdet;
  long i, j, halfs;
  double loglik_old, loglik_change = 5.0;
  
  double *xt;
  double *xw2;
  double *xw2t;
  double *tmp;
  double *beta_old;
  double *w;
  double *xw2_reduced_augmented;
  double *xw2_reduced_augmented_t;
  double *fisher_cov_reduced_augmented;
  double *cov_augmented;
  double *delta;
  int *selcol;
  
  // memory allocations
  if (NULL == (xt = (double *) R_alloc(n * k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (xw2 = (double *) R_alloc(k * n, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (xw2t = (double *) R_alloc(n * k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (tmp = (double *) R_alloc(n * k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (beta_old = (double *) R_alloc(k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (w = (double *) R_alloc(n, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (xw2_reduced_augmented = (double *) R_alloc(n*ncolfit, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (xw2_reduced_augmented_t = (double *) R_alloc(n*ncolfit, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (fisher_cov_reduced_augmented = (double *) R_alloc(ncolfit*ncolfit, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (cov_augmented = (double *) R_alloc(k*k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (delta = (double *) R_alloc(k, sizeof(double)))){ error("no memory available\n");}
  if (NULL == (selcol = (int *) R_alloc(ncolfit, sizeof(double)))){ error("no memory available\n");}
  
  //Initialise delta: 
  for(i=0; i < k; i++) {
    delta[i] = 0.0;
  }
  
  // which columns to select based on the columns to fit:
  for(i=0; i < ncolfit; i++){
    selcol[i] = colfit[i] - 1;
  }
    
  // Calculate initial likelihood
  trans(x, xt, n, k);
  XtY(xt, beta, pi, k, n, 1);	//init of pred prob
  for(i = 0; i < n; i++){
    pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));	
  }
  
  //Calculation of hat matrix diagonal; needed for loglik calculation on augmented dataset and in first iteration of main loop
  //-- Calculation of X^T W^(1/2)
  for(i = 0; i < n; i++) { 
    wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i]));
    for(j = 0; j < k; j++)
      xw2[i*k + j] = x[i + j*n] * wi; 
  }
  //-- Transpose: W^(1/2)^T X
  trans(xw2, xw2t, k, n);
  //-- Calculation of XWX
  XtXasy(xw2t, fisher_cov, n, k);
  //-- Invert:
  linpack_det(fisher_cov, &k, &logdet);
  if (logdet < (-200)) {	
    error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
  }
  else {
    linpack_inv(fisher_cov, &k); 
  }
  //-- Calculation of X^T W^(1/2) (X^TWX)^(-1)
  XtY(xw2, fisher_cov, tmp, k, n, k);
  XYdiag(tmp, xw2, Hdiag, n, k);

    *evals = 1, *iter = 0;
  int bStop = 0;
  
  
  // Calculation of loglikelihood using augmented dataset if firth:
  *loglik = 0.0;
  loglik_old = 0.0;
  for(i = 0; i < n; i++){ 
	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
    	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
        	  if(firth){
        	    // weight first replication of dataset with h_i * tau 
        	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
        	    // weight first replication of dataset with h_i * tau with opponent y
        	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
        	  }
      } else {
          warning("fitted probabilities numerically 0 or 1 occurred");
          *loglik = loglik_old;
          bStop = 1;
          break;
      }
    }

  //Calculation of initial U*:
  if(firth){
    for(i=0; i < n; i++){
      w[i] = ((weight[i] *(double)y[i]-pi[i]) + 2 * *tau * Hdiag[i] * (0.5 - pi[i]));
    }
  } else {
    for(i=0; i < n; i++){
      w[i] = weight[i] * ((double)y[i] - pi[i]);
    }
  }
  XtY(x, w, Ustar, n, k, 1);
  
  //Start of iteration: 
  if(*maxit > 0){ // in case of maxit == 0 only evaluate likelihood
    for(;;){
    //--Save iteration values:
      loglik_old = *loglik;
      copy(beta, beta_old, k);

      //--Calculation of (X^TWX)^(-1) using augmented dataset and only columns in selcol (columns to fit: colfit - 1)
      if(ncolfit > 0 && (selcol[0] != -1)) { // selcol[0] == -1 in case of just evaluating likelihood
        //-- Calculation of X^T W^(1/2)
        //---- XW^(1/2)

        for(i = 0; i < n; i++) {
          if(firth){
            wi = sqrt((weight[i] + 2 * Hdiag[i] * *tau) * pi[i] * (1.0 - pi[i])); 
          } else {
            wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i])); 
          }
          for(j = 0; j < ncolfit; j++){
            xw2_reduced_augmented[i*ncolfit + j] =  x[i + selcol[j]*n] * wi;
          }
        }

        //---- W^(1/2)^TX^T
        trans(xw2_reduced_augmented, xw2_reduced_augmented_t, ncolfit, n); 
        XtXasy(xw2_reduced_augmented_t, fisher_cov_reduced_augmented, n, ncolfit);
        linpack_inv(fisher_cov_reduced_augmented, &ncolfit);
        
        //---- Remapping of fisher_cov_(reduced)_augmented
        for(i = 0; i < k*k; i++) { //Initialisation:
          cov_augmented[i] = 0.0;
        }
        for(i=0; i < ncolfit; i++){
          for(j=0; j < ncolfit; j++) {
            cov_augmented[selcol[i] + k*selcol[j]] = fisher_cov_reduced_augmented[i + ncolfit*j];
          }
        }
        
        // Actual computation of delta:
        XtY(cov_augmented, Ustar, delta, k, k, 1);
        
        // Check for maxstep:
        double mx = maxabs(delta, k) / *maxstep;
        if(mx > 1.0){
          for(i=0; i < k; i++) {
            delta[i] /= mx;
          }
        }
      }
      
      //Update coefficient vector beta:
      for(i=0; i < k; i++){
        beta[i] += delta[i];
      }
      
      //Start step-halvings
      for(halfs = 1; halfs <= *maxhs; halfs++) {
        //Calculate loglik:
        //--Update pi:
        XtY(xt, beta, pi, k, n, 1);
        for(i = 0; i < n; i++){
          pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));	
        }
        
        //--Calculation of hat matrix diagonal;
        //-- Calculation of X^T W^(1/2)
        for(i = 0; i < n; i++) { 
          wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i]));
          for(j = 0; j < k; j++)
            xw2[i*k + j] = x[i + j*n] * wi; 
        }
        //-- Transpose: W^(1/2)^T X
        trans(xw2, xw2t, k, n);
        //-- Calculation of XWX
        XtXasy(xw2t, fisher_cov, n, k);
        //-- Invert:
        linpack_inv(fisher_cov, &k); 
        //-- Calculation of X^T W^(1/2) (X^TWX)^(-1)
        XtY(xw2, fisher_cov, tmp, k, n, k);
        XYdiag(tmp, xw2, Hdiag, n, k);
        
        // Calculation of loglikelihood using augmented dataset if firth:
        *loglik = 0.0;
        for(i = 0; i < n; i++){ 
    	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
        	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
            	  if(firth){
            	    // weight first replication of dataset with h_i * tau 
            	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
            	    // weight first replication of dataset with h_i * tau with opponent y
            	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
            	  }
          } else {
              warning("fitted probabilities numerically 0 or 1 occurred");
              *loglik = loglik_old;
              bStop = 1;
              break;
          }
        }
        //Increase evaluation counter
        (*evals)++;
        
        //Convergence check:
        if(*loglik >= (loglik_old - *lconv)){
          break;
        }
        
        //Calculation of U*: (needed as a return value)
        if(firth){
          for(i=0; i < n; i++){
            w[i] = (weight[i] *((double)y[i]-pi[i]) + 2 * *tau * Hdiag[i] * (0.5 - pi[i]));
          }
        } else {
          for(i=0; i < n; i++){
            w[i] = weight[i] * ((double)y[i] - pi[i]);
          }
        }
        XtY(x, w, Ustar, n, k, 1);
        
        //Update beta:
        for(i=0; i < k; i++){
          delta[i] /= 2.0;
          beta[i] -= delta[i];
        }
      }
      
      if(*maxhs == 0){ //if no half stepping: Update pi and compute Hdiag for the next iteration + compute loglik to check for convergence
        //Update predicted prob: 
        XtY(xt, beta, pi, k, n, 1);
        for(i = 0; i < n; i++){
          pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));	
        }
        //Calculation of hat matrix diagonal for next iteration; needed for loglik calculation on augmented dataset
        //If step halfing is activated - Hdiag is computed there
        //-- Calculation of X^T W^(1/2)
        for(i = 0; i < n; i++) { 
          wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i]));
          for(j = 0; j < k; j++)
            xw2[i*k + j] = x[i + j*n] * wi; 
        }
        //-- Transpose: W^(1/2)^T X
        trans(xw2, xw2t, k, n);
        //-- Calculation of XWX
        XtXasy(xw2t, fisher_cov, n, k);
        //-- Invert:
        linpack_det(fisher_cov, &k, &logdet);
        if (logdet < (-200)) {	
          error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
        }
        else {
          linpack_inv(fisher_cov, &k); 
        }
        //-- Calculation of X^T W^(1/2) (X^TWX)^(-1)
        XtY(xw2, fisher_cov, tmp, k, n, k);
        XYdiag(tmp, xw2, Hdiag, n, k);
        // Calculation of loglikelihood using augmented dataset if firth:
        *loglik = 0.0;
        for(i = 0; i < n; i++){ 
    	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
        	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
            	  if(firth){
            	    // weight first replication of dataset with h_i * tau 
            	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
            	    // weight first replication of dataset with h_i * tau with opponent y
            	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
            	  }
          } else {
              warning("fitted probabilities numerically 0 or 1 occurred");
              *loglik = loglik_old;
              bStop = 1;
              break;
          }
        }
        //Increase evaluation counter
        (*evals)++;
        
        //Calculation of U*:
        if(firth){
          for(i=0; i < n; i++){
            w[i] = (weight[i] *((double)y[i]-pi[i]) + 2 * *tau * Hdiag[i] * (0.5 - pi[i]));
          }
        } else {
          for(i=0; i < n; i++){
            w[i] = weight[i] * ((double)y[i] - pi[i]);
          }
        }
        XtY(x, w, Ustar, n, k, 1);
      }
      
      loglik_change = *loglik - loglik_old;
      
      //Check convergence of main loop:
      if((*iter >= *maxit) || (
        (maxabsInds(delta, selcol, ncolfit) <= *xconv) && 
          (maxabsInds(Ustar, selcol, ncolfit) < *gconv) &&
          (loglik_change < *lconv))){
        bStop = 1;
      }
      //if((*iter < *maxit) && (halfs >= *maxhs) && (*maxhs >= 5) && (*loglik < loglik_old)){
      //  bStop = 1;   // stop if half-stepping with at least 5 steps was not successful;
      //}
      if(bStop){
        break;
      }
      //Increase iteration counter
      (*iter)++;
      
    } //End of iterations
    
    if(ncolfit > 0 && (selcol[0] != -1)) {
      copy(cov_augmented, fisher_cov, k*k);
    }
    
    convergence[0] = loglik_change;
    convergence[1] = maxabsInds(Ustar, selcol, ncolfit);
    convergence[2] = maxabsInds(delta, selcol, ncolfit);
  }
}



void logistffit_IRLS(double *x, int *y, int *n_l, int *k_l, 
								double *weight, double *offset, 
								double *beta, // beta is I/O
								int *colfit, int *ncolfit_l, int *firth_l, 
								int *maxit, double *maxstep, int *maxhs,
								double *lconv, double *gconv, double *xconv, double* tau, 
								// output: 
								double *fisher_cov,		// k x k
								double *pi,						// n
								double *Hdiag,				// n
								double *loglik,
								int *evals,
								int *iter,
								double *convergence		// 3
){
	long n = (long)*n_l, k = (long)*k_l, ncolfit = (long)*ncolfit_l, firth = (long)*firth_l;
	double wi, wi_augmented, logdet;
	long i, j;
	// memory allocations

	double *xt;
	double *beta_old;
	double *xw2;
	double *xw2t;
	double *tmpNxK;
	int *selcol;
	double *newresponse; // newresponse of IRLS
	double *delta;
	double *fisher_cov_reduced_augmented;
	double *xw2_reduced_augmented;
	double *xw2t_reduced_augmented;
	double *tmp1_reduced;
	double *covs_full;
	double *xw2_reduced;

	 if (NULL == (xt = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (beta_old = (double *) R_alloc(k ,sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2 = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2t = (double *) R_alloc(n * k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (tmpNxK = (double *) R_alloc(n * k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (selcol = (int *) R_alloc(ncolfit, sizeof(int)))){error("no memory available\n");}
	 if (NULL == (newresponse = (double *) R_alloc(n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (delta = (double *) R_alloc( k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (fisher_cov_reduced_augmented = (double *) R_alloc(ncolfit* ncolfit, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2_reduced_augmented = (double *) R_alloc(ncolfit* n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2t_reduced_augmented = (double *) R_alloc(ncolfit* n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (tmp1_reduced = (double *) R_alloc(ncolfit* n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (covs_full = (double *) R_alloc(k* k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2_reduced = (double *) R_alloc(ncolfit* n, sizeof(double)))){error("no memory available\n");}
	
	
	trans(x, xt, n, k);
	
	// init loglik
	double loglik_old, loglik_change = 5.0;
	
	*evals = 0, *iter = 0;
	int bStop = 0;
	
	//init delta: difference between beta values in iteration i-1 and i
	//and beta
	for(i=0; i < k; i++){
			delta[i] = 0.0;
	 }
	
	for(i=0; i < ncolfit; i++){
	    selcol[i] = colfit[i] - 1;
	 }
	
	//Calculate initial likelihood and Hdiag for first iteration:
	// calculation of pi
	XtY(xt, beta, pi, k, n, 1);
	for(i = 0; i < n; i++){
	  pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));
	} 
	
	// XW^(1/2)
	for(i = 0; i < n; i++) {
	  wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i])); 
	  for(j = 0; j < k; j++){
	    xw2[i*k + j] = x[i + j*n] * wi;
	  }
	}
	
	//Calculation of Hat diag:
	trans(xw2, xw2t, k, n); //W^(1/2)^TX^T
	XtXasy(xw2t, fisher_cov, n, k); //X^TWX
	linpack_det(fisher_cov, &k, &logdet);
    if (logdet < (-200)) {	
        error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
    }
    else {
        linpack_inv(fisher_cov, &k); 
    }
	XtY(xw2, fisher_cov, tmpNxK, k, n, k);
	XYdiag(tmpNxK, xw2, Hdiag, n, k);
	
	// Calculation of loglikelihood using augmented dataset if firth:
	*loglik = 0.0;
	loglik_old = 0.0;
	for(i = 0; i < n; i++){ 
	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
    	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
        	  if(firth){
        	    // weight first replication of dataset with h_i * tau 
        	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
        	    // weight first replication of dataset with h_i * tau with opponent y
        	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
        	  }
      } else {
          warning("fitted probabilities numerically 0 or 1 occurred");
          *loglik = loglik_old;
          bStop = 1;
          break;
      }
	}
	(*evals)++;
	
  //Start IRLS: 
  if(*maxit > 0 && !bStop){
  	for(;;){
      loglik_old = *loglik;
      copy(beta, beta_old, k);
  
      XtY(xt, beta_old, newresponse, k, n, 1);
      for(i=0; i < n; i++){
          wi = pi[i] * (1.0 - pi[i]) * (weight[i] + 2* *tau * Hdiag[i]); //W
        if(firth){
          newresponse[i] += 1/wi*(weight[i] * ((double)y[i]-pi[i]) + 2 * *tau * Hdiag[i] * (0.5 - pi[i]));
        } else {
          newresponse[i] += 1/wi*weight[i]*(((double)y[i]-pi[i]));
        }
      }
      
      // Fisher cov based on augmented dataset and normal X^TW (see iteration formula for beta_new): 
      for(i = 0; i < n; i++) {
            wi_augmented =  pi[i] * (1.0 - pi[i])* (weight[i] + 2* *tau * Hdiag[i]); 
        for(j = 0; j < ncolfit; j++){
          xw2_reduced_augmented[i*ncolfit + j] =  x[i + selcol[j]*n] * sqrt(wi_augmented);
          xw2_reduced[i*ncolfit + j] = x[i + selcol[j]*n] * wi_augmented;
        }
      }
      
      //---- W^(1/2)^T X^T
      trans(xw2_reduced_augmented, xw2t_reduced_augmented, ncolfit, n); 
      XtXasy(xw2t_reduced_augmented, fisher_cov_reduced_augmented, n, ncolfit);
      linpack_det(fisher_cov_reduced_augmented, &ncolfit, &logdet);
      if (logdet < (-200)) {	
        error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
      }
      linpack_inv(fisher_cov_reduced_augmented, &ncolfit);
  
      	
      //(X^TWX)^(-1)X^TW
    	XY(fisher_cov_reduced_augmented, xw2_reduced, tmp1_reduced, ncolfit, ncolfit, n);
      	  
    	double tmp;
    	for(j = 0; j < ncolfit; j++) {
    	   tmp = 0.0;
    	   for(i = 0; i < n; i++){
    	       tmp += tmp1_reduced[j + i*ncolfit]*newresponse[i];
    	   }
    	   beta[selcol[j]] = tmp;
    	}

    	//Calculate likelihood and hdiag for next iteration
    	// calculation of pi
    	XtY(xt, beta, pi, k, n, 1);
    	for(i = 0; i < n; i++){
    	  pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));
    	} 
    	
    	// XW^(1/2)
    	for(i = 0; i < n; i++) {
    	  wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i])); 
    	  for(j = 0; j < k; j++){
    	    xw2[i*k + j] = x[i + j*n] * wi;
    	  }
    	}
  
    	//Calculation of Hat diag:
    	trans(xw2, xw2t, k, n); //W^(1/2)^TX^T
    	XtXasy(xw2t, fisher_cov, n, k); //X^TWX
    	linpack_det(fisher_cov, &k, &logdet);
        if (logdet < (-200)) {	
            error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
        }
        linpack_inv(fisher_cov, &k);
    	
    	XtY(xw2, fisher_cov, tmpNxK, k, n, k);
    	XYdiag(tmpNxK, xw2, Hdiag, n, k);
    	
    	// Calculation of loglikelihood using augmented dataset if firth:
    	*loglik = 0.0;
    	for(i = 0; i < n; i++){ 
    	   if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
    	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
        	  if(firth){
        	    // weight first replication of dataset with h_i * tau 
        	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
        	    // weight first replication of dataset with h_i * tau with opponent y
        	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
        	  }
    	   } else {
    	       warning("fitted probabilities numerically 0 or 1 occurred");
    	       *loglik = loglik_old;
    	       bStop = 1;
    	       break;
    	   }
    	}
    	(*evals)++;
    	
    	loglik_change = *loglik - loglik_old;
    	for(i=0; i < k; i++){
    		delta[i] = beta[i]-beta_old[i];
    	}
    	(*iter)++;
      		
    	if((*iter >= *maxit) || ((maxabsInds(delta, selcol, ncolfit) <= *xconv) && (loglik_change < *lconv)) ) {
    	    bStop = 1;
    	}
      		
    	if(bStop){
    		break;
    	}
    }
}

	
	// return adjusted vcov matrix if not all variables were fitted:
	for(i = 0; i < k*k; i++) {
	  covs_full[i] = 0.0; // init 0
	}
	if(*maxit > 0){
  	for(i=0; i < ncolfit; i++){
  	  for(j=0; j < ncolfit; j++) {
  	    covs_full[selcol[i] + k*selcol[j]] = fisher_cov_reduced_augmented[i + ncolfit*j];
  	  }   
  	}
  	copy(covs_full, fisher_cov, k*k);
  	
  	convergence[0] = loglik_change;
  	convergence[2] = maxabsInds(delta, selcol, ncolfit);
	}
}



// profile likelihood
void logistplfit(double *x, int *y, int *n_l, int *k_l, 
							double *weight, double *offset, 
							double *beta, // beta is I/O (init)
							int *iSel, int *which, double *LL0, int *firth_l, 
							// control parameter:
							int *maxit, double *maxstep, int *maxhs,
							double *lconv, double *xconv, double *tau,
							// output: 
							double *betahist,			// k * maxit
							double *loglik,				// 1
							int *iter,						// 1 
							double *convergence		// 2
							)
{
	long n = (long)*n_l, k = (long)*k_l, firth = (long)*firth_l;
	double wi, logdet;
	long i, j, halfs;
	double loglik_old, lambda, mx, wi_augmented;
	
	int bStop = 0;
	
	// memory allocations
	double *xt;
	double *beta_old;
	double *xw2;
	double *xw2t;
	double *xw2_augmented;
	double *xw2t_augmented;
	double *tmpNxK;
	double *w;
	double *tmpKx1;
	double *Vinv;
	double *cov;
	double *tmp1x1;
	double *delta;
	double *XBeta;	
	double *fisher;
	double *fisher_augmented;
	double *Ustar;
	double *pi;
	double *Hdiag;

	 if (NULL == (xt = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (beta_old = (double *) R_alloc(k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2 = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2t = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2_augmented = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (xw2t_augmented = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (tmpNxK = (double *) R_alloc(k * n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (w = (double *) R_alloc(n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (tmpKx1 = (double *) R_alloc(k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (Vinv = (double *) R_alloc(k * k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (cov = (double *) R_alloc(k * k,sizeof(double)))){error("no memory available\n");}
	 if (NULL == (tmp1x1 = (double *) R_alloc(1, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (delta = (double *) R_alloc(k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (XBeta = (double *) R_alloc(n, sizeof(double)))){error("no memory available\n");}	
	 if (NULL == (fisher = (double *) R_alloc(k * k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (fisher_augmented = (double *) R_alloc(k * k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (Ustar = (double *) R_alloc(k, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (pi = (double *) R_alloc(n, sizeof(double)))){error("no memory available\n");}
	 if (NULL == (Hdiag = (double *) R_alloc(n, sizeof(double)))){error("no memory available\n");}
	
	// init pi
	trans(x, xt, n, k);
	//Calculate initial likelihood and Hdiag for first iteration:
	// calculation of pi
	XtY(xt, beta, pi, k, n, 1);
	for(i = 0; i < n; i++){
	  pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));
	} 
	
	// XW^(1/2)
	for(i = 0; i < n; i++) {
	  wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i])); 
	  for(j = 0; j < k; j++){
	    xw2[i*k + j] = x[i + j*n] * wi;
	  }
	}
	
	//Calculation of Hat diag:
	trans(xw2, xw2t, k, n); //W^(1/2)^TX^T
	XtXasy(xw2t, fisher, n, k); //X^TWX
	linpack_det(fisher, &k, &logdet);
    if (logdet < (-200)) {	
        error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
    }
    else {
        linpack_inv(fisher, &k); 
    }
	XtY(xw2, fisher, tmpNxK, k, n, k);
	XYdiag(tmpNxK, xw2, Hdiag, n, k);
	
	// Calculation of loglikelihood using augmented dataset if firth:
	*loglik = 0.0;
	for(i = 0; i < n; i++){ 
	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
    	      *loglik += y[i] * weight[i] * log(pi[i]) + (1-y[i]) * weight[i] * log(1.0-pi[i]);
        	  if(firth){
        	    // weight first replication of dataset with h_i * tau 
        	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1-y[i]) * Hdiag[i] * *tau * log(1-pi[i]);
        	    // weight first replication of dataset with h_i * tau with opponent y
        	    *loglik += (1-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1-pi[i]);
        	  }
      } else {
          warning("fitted probabilities numerically 0 or 1 occurred");
          break;
      }
	}
	
	// Fisher cov based on augmented dataset and normal X^TW (see iteration formula for beta_new): 
    for(i = 0; i < n; i++) {
        if(firth){
            wi_augmented =  pi[i] * (1.0 - pi[i])* (weight[i] + 2 * *tau * Hdiag[i]); 
        } else {
            wi_augmented = weight[i] * pi[i] * (1.0-pi[i]);
        }
        for(j = 0; j < k; j++){
            xw2_augmented[i*k + j] =  x[i + j*n] * sqrt(wi_augmented);
        }
    }
      
    //---- W^(1/2)^T X^T
    trans(xw2_augmented, xw2t_augmented, k, n);
	XtXasy(xw2t_augmented, fisher_augmented, n, k);
    linpack_det(fisher_augmented, &k, &logdet);
    if (logdet < (-200)) {	
        error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
    }
    linpack_inv(fisher_augmented, &k);

	*iter = 0;
	for(;;) {
		//Calculation of U*:
        if(firth){
          for(i=0; i < n; i++){
            w[i] = weight[i] * ((double)y[i]-pi[i]) + 2 * *tau * Hdiag[i] * (0.5 - pi[i]);
          }
        } else {
          for(i=0; i < n; i++){
            w[i] = weight[i] * ((double)y[i] - pi[i]);
          }
        }
        XtY(x, w, Ustar, n, k, 1);
		
		XtY(Ustar, fisher_augmented, tmpKx1, k, 1, k); 
		XtY(tmpKx1, Ustar, tmp1x1, k, 1, 1); //U*^T (X^TWX)^(-1) u*
		
		double underRoot = (-1.0) * (2.0* (*LL0 - *loglik) - tmp1x1[0]) / fisher_augmented[k*((*iSel)-1) + (*iSel)-1]; 
		lambda = (underRoot < 0.0) ? 0.0 : (double)(*which) * sqrt(underRoot); 
		
		//add lambda to r-th entry in U*: 
		Ustar[(*iSel)-1] += lambda;
		XtY(fisher_augmented, Ustar, delta, k, k, 1);
        
		mx = maxabs(delta, k) / *maxstep;
		if(mx > 1.0) {
		    for(i=0; i < k; i++){
		        delta[i] /= mx;
		    }
		}
		
		for(i=0; i < k; i++){
		    beta[i] += delta[i];
		}
		
		loglik_old = *loglik;
		
		for(halfs = 0;;) {
			// calculation of pi
        	XtY(xt, beta, pi, k, n, 1);
        	for(i = 0; i < n; i++){
        	  pi[i] = 1.0 / (1.0 + exp( - pi[i] - offset[i]));
        	} 
        	
        	// XW^(1/2)
        	for(i = 0; i < n; i++) {
        	  wi = sqrt(weight[i] * pi[i] * (1.0 - pi[i])); 
        	  for(j = 0; j < k; j++){
        	    xw2[i*k + j] = x[i + j*n] * wi;
        	  }
        	}
        	
        	//Calculation of Hat diag:
        	trans(xw2, xw2t, k, n); //W^(1/2)^TX^T
        	XtXasy(xw2t, fisher, n, k); //X^TWX
        	linpack_det(fisher, &k, &logdet);
            if (logdet < (-200)) {
                error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
            }
            else {
                linpack_inv(fisher, &k); 
            }
        	XtY(xw2, fisher, tmpNxK, k, n, k);
        	XYdiag(tmpNxK, xw2, Hdiag, n, k);
        	
			// Calculation of loglikelihood using augmented dataset if firth:
        	*loglik = 0.0;
        	for(i = 0; i < n; i++){ 
        	  if(R_FINITE(log(1.0-pi[i])) && R_FINITE(log(pi[i]))){
            	      *loglik += y[i] * weight[i] * log(pi[i]) + (1.0-y[i]) * weight[i] * log(1.0-pi[i]);
                	  if(firth){
                	    // weight first replication of dataset with h_i * tau 
                	    *loglik += y[i] * Hdiag[i] * *tau * log(pi[i]) + (1.0-y[i]) * Hdiag[i] * *tau * log(1.0-pi[i]);
                	    // weight second replication of dataset with h_i * tau with opponent y
                	    *loglik += (1.0-y[i]) * Hdiag[i] * *tau * log(pi[i]) + y[i] * Hdiag[i] * *tau * log(1.0-pi[i]);
                	  }
              } else {
                  warning("fitted probabilities numerically 0 or 1 occurred");
                  bStop = 1;
                  *loglik = loglik_old;
                  break;
              }
        	}
        	
        	if(bStop){
        	    break;
        	 }
        	
        	// Fisher cov based on augmented dataset if firth
            for(i = 0; i < n; i++) {
                if(firth){
                    wi_augmented =  pi[i] * (1.0 - pi[i])* (weight[i] + 2 * *tau * Hdiag[i]); 
                } else {
                    wi_augmented = weight[i] * pi[i] * (1.0-pi[i]);
                    }
                for(j = 0; j < k; j++){
                    xw2_augmented[i*k + j] =  x[i + j*n] * sqrt(wi_augmented);
                }
            }
              
            //---- W^(1/2)^T X^T
            trans(xw2_augmented, xw2t_augmented, k, n);
        	XtXasy(xw2t_augmented, fisher_augmented, n, k);
            linpack_det(fisher_augmented, &k, &logdet);
            if (logdet < (-200)) {	
                error("Determinant of Fisher information matrix was %lf \n", exp(logdet));
            }
            linpack_inv(fisher_augmented, &k);
        	
        	halfs++;
			
			if((halfs > *maxhs) || ((fabs(*loglik - *LL0) < fabs(loglik_old - *LL0)) && (*loglik > *LL0)))
				break; 
			
			for(i=0; i < k; i++) {
				delta[i] /= 2.0;
				beta[i] -= delta[i];
			}
			
		}
		

		(*iter)++;
		
		for(i=0; i < k; i++){
		    betahist[i * (*maxit) + (*iter) - 1] = beta[i];
		}
		
		if((*iter >= *maxit) || ((fabs(*loglik - *LL0) <= *lconv) && (maxabs(delta, k) < *xconv))){
		    break;
		}
	}
	
	convergence[0] = fabs(*loglik - *LL0);
	convergence[1] = maxabs(delta, k);
}





