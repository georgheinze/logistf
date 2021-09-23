logistpl <- function(x, y, init=NULL, i, LL.0, firth, which = -1, offset=rep(0, length(y)), weight=rep(1,length(y)), plcontrol, modcontrol) {
    n<-nrow(x)
    k<-ncol(x)
    if (is.null(init)) init<-rep(0,k)
    beta<-init
    if (is.null(offset)) offset=rep(0,n)
    if (is.null(weight)) weight=rep(1,n)
    if (missing(plcontrol)) {
        plcontrol<-logistpl.control()
    }    
    if (missing(modcontrol)) {
      modcontrol <- logistf.mod.control()
    }
    tau <- modcontrol$tau
    if (!is.numeric(tau) | length(tau)>1){
      stop("Invalid value for degree of penalization tau: Must be numeric.")
    }
    if(!is.null(modcontrol$terms.fit)) {
      stop("Please call logistpl with all terms in modcontrol$terms.fit.")
    }
      
    maxit<-plcontrol$maxit
    maxstep<-plcontrol$maxstep
    maxhs<-plcontrol$maxhs
    xconv<-plcontrol$xconv
    lconv<-plcontrol$lconv
    firth <- if(firth) 1 else 0
    loglik <- iter <- warning_prob <- 0
    conv <- double(2)
    betahist <- matrix(double(k * maxit), maxit) 
    mode(x) <- mode(weight) <- mode(beta) <- mode(offset) <- mode(LL.0) <- "double"
    mode(y) <- mode(firth) <- mode(n) <- mode(k) <- "integer"
    mode(maxstep) <- mode(lconv) <- mode(xconv) <- mode(loglik) <- mode(tau) <- "double"
    mode(maxit) <- mode(maxhs) <- mode(i) <- mode(which) <- mode(iter) <- mode(warning_prob) <- "integer"
    
    res <- .C("logistplfit", x, y, n, k, weight, offset, beta=beta, i, which, LL.0, firth, maxit, 
    maxstep, maxhs, lconv, xconv, tau, betahist=betahist, loglik=loglik, iter=iter, conv=conv, warning_prob = warning_prob,
    PACKAGE="logistf")
    
    #if(res$iter>=maxit){
    #warning(paste("Maximum number of iterations exceeded. Try to increase the number of iterations or alter step size by passing 'pl.control(maxit=..., maxstep=...)' to parameter plcontrol"))
    #}
    if(res$warning_prob){
      warning("fitted probabilities numerically 0 or 1 occurred")
    }
    
    res <- res[c("beta", "betahist", "loglik", "iter", "conv")]
    res$betahist <- head(res$betahist, res$iter)
    res$beta <- res$beta[i]
    res
}

