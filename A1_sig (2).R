EN.test.sig<-spNmix.sig(n = EN14$y,X = EN14$X, M = 1000, niters = 50000,
                        xlims = EN14$xlims, ylims = EN14$ylims, tune=c(0.2,10,0.15,5), cov=EN14$forest)
#tune is sd for distribution of sigint, sigbeta, distribution of lam0, distribution of S
spNmix.sig <- function(n, X, M, niters, xlims, ylims, tune=c(0.1,10, 0.15, 2),cov,
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
  sigma <-((cov/1)^sigint)*sigbeta
  #sigma <-runif(1, 500, 2000)
  lam0 <- runif(1,.1, 1)
  #lam is a matrix same size as D, of detection probabilities for each individual at each trap
  lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
  w <- rbinom(M, 1, .5)
  psi <- runif(1, .2, .8)
  #lamv.curr is vector of length ntraps of total detection probabilities
  lamv.curr <- colSums(lam*w)
  # just in case the first sigma is rejected
  ll <- sum(dpois(n, lamv.curr, log=TRUE))
  
  # matrix to hold samples
  out <- matrix(NA, nrow=niters, ncol=6)
  colnames(out) <- c("sigint","sigbeta","lam0", "psi", "N", "ll")
  
  Sout <- NULL
  if(monitorS)
    Sout <- array(NA, c(M, 2, niters))
  
  cat("\ninitial values =", c(sigint, sigbeta, lam0, psi, sum(w)), "\n\n")
  
  for(iter in 1:niters) {
    
    if(iter %% 100 ==0) {
      cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
      cat("current =", out[iter-1,], llscand, "\n")
      cat("  Acceptance rates\n")
      cat("    S =", Sups/M, "\n")
      cat("    w =", wUps/M, "\n")
      
    }
    
    # update sigbeta
    sigint.cand<-rnorm(1,sigint,tune[1])
    sigbeta.cand<-rnorm(1,sigbeta,tune[2])
#monmolecular
    sigma.cand <-((cov/1)^sigint.cand)*sigbeta.cand
#log
#    sigma.cand <- exp(log(sigbeta.cand)+(sigint.cand*cov))
    if(!(0%in%sigma.cand)) {
      lam.cand <- t(lam0*exp(t(-(D*D))/(2*sigma.cand*sigma.cand)))
      lamv.cand <- colSums(lam.cand*w)	    
      #when sigma is too small p becomes zero, which creates -inf log likelihood at sites with a detection.  this replaces zero with minimum other value
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
    
    # update lam0
    lam0.cand <- rnorm(1, lam0, tune[3])
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
      Scand <-c(rnorm(1, S[i,1], tune[4]), rnorm(1, S[i,2], tune[4]))
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
    out[iter,] <- c(sigint,sigbeta,lam0,psi,sum(w),ll )
    if(monitorS)
      Sout[1:sum(w),,iter] <- S[w==1,]
  }
  last <- list(S=S, lam=lam, w=w)
  list(out=out, last=last, Sout=Sout)
}

summary(EN.test.sig$out)

# Histories
par(mfrow=c(2,2))
plot(EN.test.sig$out[1:11000,1], type="l", ylab="sigmaint")
plot(EN.test.sig$out[1:11000,2], type="l", ylab="sigmabeta")
plot(EN.test.sig$out[1:11000,3], type="l", ylab="lam0")
plot(EN.test.sig$out[1:11000,4], type="l", ylab="psi")
plot(EN.test.sig$out[1:50000,5], type="l", ylab="N")

# Densities
par(mfrow=c(2,2))
plot(density(EN.test.sig$out[1001:11000,1]), slab = "sigint", main = "")
plot(density(EN.test.sig$out[1001:11000,2]), xlab="sigmabeta", main="")
plot(density(EN.test.sig$out[1001:11000,3]), xlab="lam0", main="")
plot(density(EN.test.sig$out[1001:11000,4]), xlab="psi", main="")
hist(EN.test.sig$out[1001:11000,5], xlab="N", freq=FALSE, main="")

hist(EN.test.sig$out[10001:50000,5]+(0.01*Area), breaks = 20, freq = FALSE, col = rgb(0.1,0.1,0.1,0.5), xlim = c(0,600), xaxt = "n", yaxt = "n", ylim = c(0, 0.01))
hist(EN.10.15.2$out[10001:50000,4], xaxt = "n", yaxt = "n", freq = FALSE, col = rgb(0.8,0.8,0.8,0.5), breaks = 15, add = T)
axis(side = 1, at = seq(0, 900, 100), labels = seq(round(0/Area,2), round(900/Area, 2), round(100/Area, 2)), cex.axis = 1.5)
axis(side = 2, at = seq(0, 0.01, 0.001), labels = seq(0, 0.10, 0.01), cex.axis = 1.5)
box(col = "black")