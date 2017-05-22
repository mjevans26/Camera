spNmixZ <- function(n, X, M, l, niters, xlims, ylims, tune=c(0.1, 0.1, 3))
{
  #J = traps
  J <- nrow(n)
  #K = occassions
  K <- ncol(n)
  #L = number marked individuals
  I <- nrow(l)
  # initial values
  S <- cbind(runif(M,0,15), runif(M,0,15))
  D <- e2dist1(S, X)
  sigma <-runif(1, .25, 1.25)
  lam0 <- runif(1, .1, 1)
  lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
  w <- c(rep(1,I),rbinom(M-I, 1, .7))
  psi <- runif(1, .2, .8)
  #z = array of possible individuals * traps * occassions
  z <- array(NA, c(M,J,K))
  for(i in 1:I){
    for(k in 1:K){
      z[i,l[i,j],k] <- 1
    }
  }
  for(j in 1:J) {
    for(k in 1:K) {
      z[I+1:M,j,k] <- rmultinom(1, n[j,k] - sum(z[1:I,j,k]), lam[,j]*w)
    }
  }
  # in case the first sigma.cand is rejected
  ll <- sum(dpois(z, lam*w, log=TRUE))
  
  # matrix to hold samples
  out <- matrix(NA, nrow=niters, ncol=4)
  colnames(out) <- c("sigma", "lam0", "psi", "N")
  
  cat("\ninitial values =", c(sigma, lam0, psi, sum(w)), "\n\n")
  
  for(iter in 1:niters) {
    
    if(iter %% 100 ==0) {
      cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
      cat("current =", out[iter-1,], "\n")
      cat("  Acceptance rates\n")
      cat("    S =", Sups/M, "\n")
      cat("    w =", wUps/M, "\n")
    }
    
    # update sigma
    sigma.cand <- rnorm(1, sigma, tune[1])
    if(sigma.cand > 0) {
      lam.cand <- lam0*exp(-(D*D)/(2*sigma.cand*sigma.cand))
      ll<- sum(dpois(z, lam*w, log=TRUE))
      llcand<- sum(dpois(z, lam.cand*w, log=TRUE) )
      if(runif(1)<exp( llcand  - ll ) ){
        ll <- llcand
        lam <- lam.cand
        sigma <- sigma.cand
      }
    }
    
    # update lam0
    lam0.cand <- rnorm(1, lam0, tune[2])
    if(lam0.cand>0) {
      lam.cand <- lam0.cand*exp(-(D*D)/(2*sigma*sigma))
      llcand<- sum(dpois(z, lam.cand*w, log=TRUE) )
      if(runif(1) < exp( llcand - ll ) ) {
        lam0<-lam0.cand
        lam<-lam.cand
        ll <- llcand
      }
    }
    
    # update w
    wUps <- 0
    for(i in 1:M) {
      wcand<-w
      wcand[i] <- if(w[i]==0) 1 else 0
      llcand <- sum(dpois(z, lam*wcand, log=TRUE) )
      prior <- dbinom(w[i], 1, psi, log=TRUE)
      prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
      if(runif(1) < exp((llcand+prior.cand) - (ll+prior))) {
        w <- wcand
        ll <- llcand
        wUps <- wUps+1
      }
    }
    
    # update psi
    psi <- rbeta(1, 1+sum(w), 1+M-sum(w))
    
    # update S
    Sups <- 0
    for(i in 1:M) {
      Scand <- c(rnorm(1, S[i,1], tune[3]),
                 rnorm(1, S[i,2], tune[3]))
      inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
        Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
      if(!inbox)
        next
      dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
      lam.cand <- lam
      lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
      llcand <- sum(dpois(z, lam.cand*w, log=TRUE) )
      if(runif(1)< exp(llcand - ll)) {
        ll <- llcand
        S[i,] <- Scand
        lam <- lam.cand
        D[i,] <- dtmp
        Sups <- Sups+1
      }
    }
    
    # update z
    for(j in 1:J) {
      zip <- lam[,j]*w
      probs <- zip/sum(zip)
      for(k in 1:K)
        z[,j,k] <- rmultinom(1, n[j,k], probs)
    }
    #        }
    
    ll <- sum(dpois(z, lam*w, log=TRUE))
    
    out[iter,] <- c(sigma,lam0,psi,sum(w) )
  }
  last <- list(S=S, lam=lam, w=w)
  list(out=out, last=last)
}

