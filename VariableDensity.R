EN14 <- readRDS(file = "input_data.rds")
library(rgdal)
DEV<-readOGR("I:/MikeGradSchool/PhD/GIS/SMRC", "NEB_DensClass_Simple_SPm")
bounds <- readOGR("I:/MikeGradSchool/PhD/GIS/SMRC/Buffers", "NEB10kBuff")
DEV@data[,1] <- as.numeric(DEV@data[,1])

psi <- c(rbeta(1, 1+sum(w[which(w[,2]==1), 1]), 1+M-sum(w[which(w[,2]==1), 1])),
         rbeta(1, 1+sum(w[which(w[,2]==2), 1]), 1+M-sum(w[which(w[,2]==2), 1])),
         rbeta(1, 1+sum(w[which(w[,2]==3), 1]), 1+M-sum(w[which(w[,2]==3), 1]))
)

spNmix_den <- function(n, X, M, niters, tune=c(0.2, 10, 0.2, 0.2, 5),cov,
                       mask, area, monitorS=FALSE)
{
  K <- ncol(n)
  levs <- length(unique(mask@data[,1]))
  xlims <- area@polygons[[1]]@Polygons[[1]]@coords[,1]
  ylims <- area@polygons[[1]]@Polygons[[1]]@coords[,2]
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
  w <- matrix(c(rbinom(M, 1, .5), over(pts, mask)[,1]), ncol = 2)
  #psi <- runif(1, .2, .8)
  ## psi as vector of inclusion probabilities for individuals with one
  ## value per level of density mask
  psi <- runif(levs, .2, .8)
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
  print(iter)
  #print(psi)
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
    if(!(0 %in% sigma.cand)) {
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
    if(!(0 %in% lam0.cand)) {
      lam.cand <- t(lam0.cand*exp(t(-(D*D))/(2*sigma*sigma)))
      lamv.cand <- colSums(lam.cand*w[,1])
      lamv.cand[lamv.cand==0]<-min(lamv.cand[lamv.cand!=0])
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      if(runif(1) < exp( llcand - ll ) ) {
        ll <- llcand
        lamv.curr<-lamv.cand
        lam0<-lam0.cand
        lam<-lam.cand
      }
    }
    
    #print(c(NA %in% w[,1], NA %in% w[,2]))
    # update w, which is a vector of 0 and 1, of length M
    wUps <- 0
    for(i in 1:M) {
      wcand <- w
      wcand[i,1] <- if(w[i,1]==0) 1 else 0
      lamv.cand <- colSums(lam*wcand[,1])
      llcand <- sum(dpois(n, lamv.cand, log=TRUE) )
      prior <- dbinom(w[i,1], 1, psi[w[i,2]], log=TRUE)
      prior.cand <- dbinom(wcand[i,1], 1, psi[w[i,2]], log=TRUE)
      #tryCatch(
      if(runif(1) < exp((llcand + prior.cand) - (ll + prior))) {
        w <- wcand
        lamv.curr <- lamv.cand
        ll <- llcand
        wUps <- wUps+1
      }#, error = return(out))
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
      #inbox <- Scand[1]>=xlims[1] & Scand[1]<=xlims[2] &
        #Scand[2]>=ylims[1] & Scand[2]<=ylims[2]
      ##alternative using sp to check that new activity center is within 
      ## study area
      #if(is.na(over(pts[i], mask)[,1]))
      if(point.in.polygon(Scand[1], Scand[2], xlims, ylims) != 1)
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
    ## assign mask region value to w based on new locations for S
    pts <- SpatialPoints(S, proj4string = mask@proj4string)
    w[,2] <- over(pts, mask)[,1]
    w[which(is.na(w[,2])),2] <- which(psi == min(psi))
    out[iter,] <- c(lamint, lambeta, sigint, sigbeta, psi, sum(w[,1]), ll )
    if(monitorS)
      Sout[1:sum(w[,1]),,iter] <- S[w==1,]
  }
  last <- list(S=S, lam=lam, w=w)
  list(out=out, last=last, Sout=Sout)
}

test <- spNmix_den(EN14$y, EN14$X, 700, 50000, tune=c(0.2, 10, 0.2, 0.2, 5), EN14$forest,
                   DEV, bounds, monitorS=FALSE)