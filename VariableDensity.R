
library(rgdal)
DEV<-readOGR("D:/MikeGradSchool/PhD/GIS/SMRC", "NEB_DensClass_SPm")
xlims <- DEV@bbox[1,]
ylims <- DEV@bbox[2,]
pts <- SpatialPoints(S, proj4string = DEV@proj4string)
w[,2] <- over(pts, DEV)[,1]

psi <- c(rbeta(1, 1+sum(w[which(w[,2]==1), 1]), 1+M-sum(w[which(w[,2]==1), 1])),
         rbeta(1, 1+sum(w[which(w[,2]==2), 1]), 1+M-sum(w[which(w[,2]==2), 1])),
         rbeta(1, 1+sum(w[which(w[,2]==3), 1]), 1+M-sum(w[which(w[,2]==3), 1]))
)

spNmix_den <- function(n, X, M, niters, xlims, ylims, tune=c(0.2, 10, 0.2, 0.2, 5),cov,
                       mask, monitorS=FALSE)
{
  K <- ncol(n)
  levs <- length(unique(mask@data[,1]))
  ## initial values
  ## generate random points inside mask instead of bounding box
  pts <- spsample(mask, M, type= "random")
  S <- pts@coords
  #S <- cbind(runif(M, xlims[1], xlims[2]),
  #           runif(M, ylims[1], ylims[2]))
  ##D is a matrix of distances between activity centers (rows) and traps (cols)
  D <- e2dist1(S, X)
  sigint<-runif(1,0,2)
  sigbeta<-runif(1,500,2000)
  sigma <-((cov/max(cov))^sigint)*sigbeta
  lamint <- runif(1,-2,2)
  lambeta<-runif(1,2,7)
  lam0 <- exp(lamint + (lambeta*cov))/(1+exp(lamint + (lambeta*cov)))
  lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
  #w <- rbinom(M, 1, .5)
  ## create w as matrix with 2nd column indicating which portion of
  ## density mask individual i falls in.  Starts at 1
  w <- matrix(c(rbinom(M, 1, .5), rep(1,M)), ncol = 2)
  #psi <- runif(1, .2, .8)
  ## psi as vector of inclusion probabilities for individuals with one
  ## value per level of density mask
  psi <- rep(runif(1, .2, .8), levs)
  lamv.curr <- colSums(lam*w[,1])
  # just in case the first sigma is rejected
  ll <- sum(dpois(n, lamv.curr, log=TRUE))
  
  # matrix to hold samples
  out <- matrix(NA, nrow=niters, ncol=6+levs)
  colnames(out) <- c("lamint", "lambeta", "sigint","sigbeta", paste("psi",1:levs), "N", "ll")
  
  Sout <- NULL
  if(monitorS)
    Sout <- array(NA, c(M, 2, niters))
  
  cat("\ninitial values =", c(lamint, lambeta, sigint, sigbeta, psi, sum(w[,1])), "\n\n")
  
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
      lamv.cand <- colSums(lam.cand*w[,1])	    
      #When sigma is too small p becomes zero, which creates -inf log likelihood at sites with a detection.  
      #This replaces zero with minimum other value
      lamv.cand[lamv.cand==0]<-min(lamv.cand[lamv.cand!=0])
      lamv.curr <- colSums(lam*w[,1])
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
      lamv.cand <- colSums(lam.cand*w[,1])
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
      wcand[i,1] <- if(w[i,1]==0) 1 else 0
      lamv.cand <- colSums(lam*wcand[,1])
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      prior <- dbinom(w[i,1], 1, psi[w[i,2]], log=TRUE)
      prior.cand <- dbinom(wcand[i,1], 1, psi[w[i,2]], log=TRUE)
      if(runif(1) < exp((llcand+prior.cand) - (ll+prior))) {
        w <- wcand
        lamv.curr <- lamv.cand
        ll <- llcand
        wUps <- wUps+1
      }
    }
    
    # update psi
    psi <- vapply(1:levs, function(x){
      rbeta(1, 1+sum(w[which(w[,2]==x), 1]), 1+M-sum(w[which(w[,2]==x), 1]))
    },
    FUN.VALUE = 0,
    USE.NAMES = FALSE)
    
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
      lamv.cand <- colSums(lam.cand*w[,1])
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
    out[iter,] <- c(lamint, lambeta, sigint, sigbeta, psi, sum(w[,1]), ll )
    if(monitorS)
      Sout[1:sum(w[,1]),,iter] <- S[w==1,]
  }
  last <- list(S=S, lam=lam, w=w)
  list(out=out, last=last, Sout=Sout)
}