#' Calculate distance between sample locations and activity centers
#' Faster than using loops
#' @param x matrix of activity center coordinates
#' @param y matrix of sampling location coordinates
e2dist1 <- function (x, y)
{
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

#' Draw posterior samples using MCMC
#'
#' @param n numeric matrix (k x j) of detection counts
#' @param X numeric matrix (2 x j) with X and Y coordinates of sampling locations
#' @param M upper limit of population size
#' @param niters number of MCMC iterations
#' @param xlims vector of maximum and minimum x coordinates of buffered sampling area
#' @param ylims vector of maximum and minimum y coordinates of buffered sampling area
#' @param tune numeric vector standard deviation for parameter resampling
#' @param cov vector of covariate values at sampling locations
#' @examples
#'spNmix_het(n = EN$Y, X = EN$X, M = 1000, niters = 100000, xlims = EN$xlims, ylims = EN14$ylims,
#'            tune = c(0.2, 10, 0.2, 0.2, 5), cov = EN14$forest, monitorS = FALSE)

spNmix_het <- function(n, X, M, niters, xlims, ylims, tune=c(0.2, 10, 0.2, 0.2, 5),cov,
                       monitorS=FALSE)
{
  K <- ncol(n)
  
  # initial values
  S <- cbind(runif(M, xlims[1], xlims[2]),
             runif(M, ylims[1], ylims[2]))
  #D is a matrix of distances between activity centers (rows) and traps (cols)
  D <- e2dist1(S, X)
  sigint<-runif(1,0,2)
  sigbeta<-runif(1,500,2000)
  sigma <-((cov/max(cov))^sigint)*sigbeta
  lamint <- runif(1,-2,2)
  lambeta<-runif(1,2,7)
  lam0 <- exp(lamint + (lambeta*cov))/(1+exp(lamint + (lambeta*cov)))
  lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
  w <- rbinom(M, 1, .5)
  psi <- runif(1, .2, .8)
  lamv.curr <- colSums(lam*w)
  # just in case the first sigma is rejected
  ll <- sum(dpois(n, lamv.curr, log=TRUE))
  
  # matrix to hold samples
  out <- matrix(NA, nrow=niters, ncol=7)
  colnames(out) <- c("lamint", "lambeta", "sigint","sigbeta", "psi", "N", "ll")
  
  Sout <- NULL
  if(monitorS)
    Sout <- array(NA, c(M, 2, niters))
  
  cat("\ninitial values =", c(lamint, lambeta, sigint, sigbeta, psi, sum(w)), "\n\n")
  
  for(iter in 1:niters) {
    
    if(iter %% 100 ==0) {
      cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
      cat("current =", out[iter-1,], llscand, "\n")
      cat("  Acceptance rates\n")
      cat("    S =", Sups/M, "\n")
      cat("    w =", wUps/M, "\n")
      
    }
    
    # update sigbeta/sigint
    sigint.cand<-rnorm(1,sigint,tune[1])
    sigbeta.cand<-rnorm(1,sigbeta,tune[2])
    sigma.cand <-((cov/max(cov))^sigint.cand)*sigbeta.cand
    if(!(0%in%sigma.cand)) {
      lam.cand <- t(lam0*exp(t(-(D*D))/(2*sigma.cand*sigma.cand)))
      lamv.cand <- colSums(lam.cand*w)	    
      #When sigma is too small p becomes zero, which creates -inf log likelihood at sites with a detection.  
      #This replaces zero with minimum other value
      lamv.cand[lamv.cand==0]<-min(lamv.cand[lamv.cand!=0])
      lamv.curr <- colSums(lam*w)
      ll<- sum(dpois(n, lamv.curr, log=TRUE))
      llscand<- sum(dpois(n, lamv.cand, log=TRUE) )
      if(runif(1) < exp( llscand  - ll ) ){
        ll <- llscand
        lamv.curr <- lamv.cand
        lam <- lam.cand
        sigma <- sigma.cand
        sigbeta<-sigbeta.cand
        sigint<-sigint.cand
      }
    }
    
    # update lambeta
    lamint.cand<-rnorm(1, lamint,tune[3])
    lambeta.cand<-rnorm(1,lambeta,tune[4])
    lam0.cand <- exp(lamint.cand + (lambeta.cand*cov))/(1+exp(lamint.cand + (lambeta.cand*cov)))
    if(lam0.cand>0) {
      lam.cand <- t(lam0.cand*exp(t(-(D*D))/(2*sigma*sigma)))
      lamv.cand <- colSums(lam.cand*w)
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      if(runif(1) < exp( llcand - ll ) ) {
        ll <- llcand
        lamv.curr<-lamv.cand
        lam0<-lam0.cand
        lam<-lam.cand
      }
    }
    
    
    # update w, which is a vector of 0 and 1, of length M
    wUps <- 0
    for(i in 1:M) {
      wcand <- w
      wcand[i] <- if(w[i]==0) 1 else 0
      lamv.cand <- colSums(lam*wcand)
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      prior <- dbinom(w[i], 1, psi, log=TRUE)
      prior.cand <- dbinom(wcand[i], 1, psi, log=TRUE)
      if(runif(1) < exp((llcand+prior.cand) - (ll+prior))) {
        w <- wcand
        lamv.curr <- lamv.cand
        ll <- llcand
        wUps <- wUps+1
      }
    }
    
    # update psi
    psi <- rbeta(1, 1+sum(w), 1+M-sum(w))
    
    # update S
    Sups <- 0
    for(i in 1:M) {
      Scand <-c(rnorm(1, S[i,1], tune[5]), rnorm(1, S[i,2], tune[5]))
      inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
        Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
      if(!inbox)
        next
      dtmp <- sqrt( (Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2 )
      lam.cand <- lam
      lam.cand[i,] <-  lam0*exp(-(dtmp*dtmp)/(2*sigma*sigma))
      lamv.cand <- colSums(lam.cand*w)
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      if(runif(1)< exp(llcand - ll)) {
        ll <- llcand
        lamv.curr <- lamv.cand
        S[i,] <- Scand
        lam <- lam.cand
        D[i,] <- dtmp
        Sups <- Sups+1
      }
    }
    out[iter,] <- c(lamint, lambeta, sigint, sigbeta, psi, sum(w), ll )
    if(monitorS)
      Sout[1:sum(w),,iter] <- S[w==1,]
  }
  last <- list(S=S, lam=lam, w=w)
  list(out=out, last=last, Sout=Sout)
}

summary(EN.test.sig$out)

# Histories
plot(EN.test.sig$out[1:110000,1], type="l", ylab="lamint")
plot(EN.test.sig$out[1:110000,2], type="l", ylab="lambeta")
plot(EN.test.sig$out[1:110000,3], type="l", ylab="sigmaint")
plot(EN.test.sig$out[1:110000,4], type="l", ylab="sigmabeta")
plot(EN.test.sig$out[1:110000,5], type="l", ylab="lam0")
plot(EN.test.sig$out[1:110000,6], type="l", ylab="psi")
plot(EN.test.sig$out[1:110000,7], type="l", ylab="N")

# Densities
plot(density(EN.test.sig$out[10001:110000,1]), slab = "lamint", main = "")
plot(density(EN.test.sig$out[10001:110000,2]), slab = "lambeta", main = "")
plot(density(EN.test.sig$out[10001:110000,3]), slab = "sigint", main = "")
plot(density(EN.test.sig$out[10001:110000,4]), xlab="sigmabeta", main="")
plot(density(EN.test.sig$out[10001:110000,5]), xlab="lam0", main="")
plot(density(EN.test.sig$out[10001:110000,6]), xlab="psi", main="")
hist(EN.test.sig$out[1001:11000,7], xlab="N", freq=FALSE, main="")
