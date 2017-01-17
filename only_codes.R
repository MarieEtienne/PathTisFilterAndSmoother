rm(list=ls())
set.seed(21)
nYear <- 30
biomass <- rep(NA, nYear)
capt <- rep(NA, nYear)
K <- 300 ## carrying capacity
r <- 1.02 ## recruitement rate
q <- 0.5
epsilon <- 10

sigmaBiomass <- 0.3
biomass[1] <- rnorm(1, mean=K, sd=K/10)
capt[1] <- runif(1,0.1,0.5)*biomass[1]
mn_innov = -sigmaBiomass^2/2
for( year in 2:nYear){
  Bk = biomass[year-1]
  #Captures are a random fraction of biomass
  innov <- rnorm(1, mean=mn_innov,sd=sigmaBiomass)
  biomass[year] <- (Bk+r*Bk*(1- Bk/K)-capt[year-1])*exp(innov)
  capt[year] <- runif(1,0.1,0.5)*biomass[year]
}
sqrtD <- sqrt((1+r)^2-4*capt*r/K)
ub <- (-(1+r)-sqrtD)/(-2*r/K)
lb <- (-(1+r)+sqrtD)/(-2*r/K)
bounds <- cbind(lb,ub)
ub+r*ub*(1-ub/K)-capt

plot(1:nYear,biomass, ylim=range(c(0,range(biomass),range(bounds))), 
     main='Evolution de la biomasse',xlab="Année",type="b",pch=20,
     ylab="Masse",bty="n")
lines(1:(nYear),capt,type="b",col="red",pch=20)
legend("topright",pch=20,col=c("black","red"),
       legend = c("Biomasse","Captures"),
       bty="n",horiz = T)
matplot(bounds,pch=20,add=T,col="black",type="l",lty=3)

obs_error <- rnorm(nYear, mean=mn_innov, sd=sigmaBiomass)
abundanceIndex <- q * biomass*exp(obs_error)
lines(1:nYear, abundanceIndex, type='b',col="blue",pch=18,
      ylim=range(c(0,range(abundanceIndex))))
axis(side=4,at=pretty(range(abundanceIndex),n = 2),col = "blue",
     col.ticks = "blue",labels=F ,lwd.ticks=0)
mtext(side=4,text = "Indice d'abondance",col = "blue",at=150,line=1)
mtext(side=4,text = pretty(range(abundanceIndex),n=2),col = "blue",
      at=pretty(range(abundanceIndex),n=2))
# On what conditions the deterministic term of population dynamics --------



r0 <- function(G, theta,method="M",Y0=NA,C0=NA){
  if(method=="M"){
    return(exp(runif(G, min=log(theta$K/5), max=log(5*theta$K))))
  }
  if(method=="P"){
    sqrtD <- sqrt((1+theta$r)^2-4*C0*theta$r/theta$K)
    ub <- (-(1+theta$r)-sqrtD)/(-2*theta$r/theta$K)
    lb <- (-(1+theta$r)+sqrtD)/(-2*theta$r/theta$K)
    X0 <- sapply(1:G,function(i){
      ok <- F
      while(!ok){
        cand <- exp(rnorm(1,log(Y0)-log(theta$q)-theta$mn_innov,
                          theta$sigmaBiomass))
        ok <- (cand < ub) & (cand>lb)
      }
      return(cand)
    })
    return(X0)
  }
}

dr0 <- function(x0, theta,method="M",Y0=NA,C0=NA){
  if(method=="M"){
    return(rep(1/length(x0), length(x0) ) )
  }
  if(method=="P"){
    sqrtD <- sqrt((1+theta$r)^2-4*C0*theta$r/theta$K)
    ub <- (-(1+theta$r)-sqrtD)/(-2*theta$r/theta$K)
    lb <- (-(1+theta$r)+sqrtD)/(-2*theta$r/theta$K)
    mn_ref <- log(Y0)-log(q)-theta$mn_innov
    norm_cst <- (pnorm(log(ub),mn_ref,theta$sigmaBiomass)-
                   pnorm(log(lb),mn_ref,theta$sigmaBiomass))
    dx <- dnorm(log(x0),mn_ref,theta$sigmaBiomass)/norm_cst
    return(dx)
  }
}

rPropag <- function(x,theta, sdPropag=sigmaBiomass,
                    min=NA,
                    method="M",Ynex=NA,
                    Cnex=NA){
  ## min is the minimal acceptable value
  if(method=="M"){
    meanValue <- ifelse(x>min, x-min, 1)
    xprime <- min+exp(rnorm(length(x), 
                            mean=log(meanValue)-theta$sdPropag^2/2,
                            sd=theta$sdPropag))
    return(xprime)
  }
  if(method=="P"){
    sqrtD <- sqrt((1+theta$r)^2-4*Cnex*theta$r/theta$K)
    ub <- (-(1+theta$r)-sqrtD)/(-2*theta$r/theta$K)
    lb <- (-(1+theta$r)+sqrtD)/(-2*theta$r/theta$K)
    v_gau  <- (theta$sdPropag^2*theta$sdObs^2)/(theta$sdPropag^2+theta$sdObs^2)
    obs    <- exp(log(Ynex)-log(q)-theta$mn_innov)
    mn_gau <- v_gau*(x/(theta$sdPropag^2)+obs/(theta$sdObs^2))
    Xnex <- sapply(mn_gau,function(mn){
      ok <- F
      while(!ok){
        cand <- rnorm(1,mn,sqrt(v_gau))
        ok <- cand > lb & cand < ub
      }
      return(cand)
    })
    return(Xnex)
  }
}


dPropag <- function(y, x,theta,sdPropag=sigmaBiomass, min=NA,
                    method="M",Ynex=NA,
                    Cnex=NA){
  if(method=="M"){
    meanValue <- ifelse(x>min, x-min, 1)
    return(dnorm(log(y-min), mean=log(meanValue)-sdPropag^2/2, sd=sdPropag))
  }
  if(method=="P"){
    sqrtD <- sqrt((1+theta$r)^2-4*Cnex*theta$r/theta$K)
    ub <- (-(1+theta$r)-sqrtD)/(-2*theta$r/theta$K)
    lb <- (-(1+theta$r)+sqrtD)/(-2*theta$r/theta$K)
    v_gau = (theta$sdPropag^2*theta$sdObs^2)/(theta$sdPropag^2+theta$sdObs^2)
    obs <- exp(log(Ynex)-log(q)-theta$mn_innov)
    dnex <- sapply(1:length(y),function(i){
      mn = v_gau*(x[i]/(theta$sdPropag^2)+obs/(theta$sdObs^2))
      dnorm(y[i],mn,v_gau)/(pnorm(ub,mn,sqrt(v_gau))-pnorm(lb,mn,sqrt(v_gau)))
    #normalized by the area of possible zone
    })
    return(dnex)
  }
}

mtheta<- function(x,xprime, theta, capt){
  dnorm(log(xprime), 
        mean=log(x+theta$r*x*(1-x/theta$K)- capt)-
          theta$mn_innov, 
        sd=theta$sigmaBiomass)
}

g<-function(x,y,theta){
  dnorm(log(y), mean=log(theta$q*x)-theta$mn_innov,
        sd=theta$sigmaBiomass)
}


sdPropag <- 10
sdObs <- 50
theta=list(q=q, K=K, sigmaBiomass=sigmaBiomass, r=r,
           mn_innov=mn_innov,sdPropag=sdPropag,sdObs=sdObs)
G <- 100
X <- matrix(NA, ncol=nYear, nrow=G)
X[,1] <- r0(G, theta,method = "P",Y0 = abundanceIndex[1],C0 = capt[1] )
weightX <- matrix(NA, ncol=nYear, nrow=G)
weightX[,1] <- dr0(x0=X[,1],theta = theta,method = "P",
                   Y0 = abundanceIndex[1],C0 = capt[1])
weightX[,1] <- weightX[,1]/sum(weightX[,1])
points(rep(1,G),X[,1],cex=10*weightX[,1],col="darkgreen")

method="P"
for(year in 2:nYear){
  X[,year] <- rPropag(x = X[,year-1],theta = theta,
                      Ynex = abundanceIndex[year],
                      Cnex = capt[year],
                      sdPropag = sdPropag,
                      min=capt[year-1],
                      method=method)
  mth <- mtheta(x=X[,year-1],xprime= X[,year],
                theta= theta, capt=capt[year-1])
  gfact <- g(X[,year], abundanceIndex[year], theta)
  pfact <- dPropag(y = X[,year],x= X[,year-1],theta = theta,
                   sdPropag = sdPropag,method =method,
                   Ynex = abundanceIndex[year],Cnex = capt[year],
                   min=capt[year-1])
  weightX[,year]=(weightX[,year-1]
                  *mth*gfact/pfact)
  weightX[is.nan(weightX[,year]),year] <- 0
  weightX[,year] <-   weightX[,year]/sum(  weightX[,year])
  points(rep(year,G),X[,year],cex=10*weightX[,year],col="darkgreen")
}

sdPropag <- 50
sdObs <- 50
theta=list(q=q, K=K, sigmaBiomass=sigmaBiomass, r=r,
           mn_innov=mn_innov,sdPropag=sdPropag,sdObs=sdObs)
G <- 50
X <- matrix(NA, ncol=nYear, nrow=G)
X[,1] <- r0(G, theta,method = "P",Y0 = abundanceIndex[1],C0 = capt[1] )
weightX <- matrix(NA, ncol=nYear, nrow=G)
weightX[,1] <- dr0(x0=X[,1],theta = theta,method = "P",
                   Y0 = abundanceIndex[1],C0 = capt[1])
weightX[,1] <- weightX[,1]/sum(weightX[,1])
points(rep(1,G),X[,1],cex=10*weightX[,1],col="darkgreen")

method="P"
for(year in 2:nYear){
  inds_rspl <- sample(1:G,size=G,replace=T,prob = weightX[,year-1])
  X[,year] <- rPropag(x = X[inds_rspl,year-1],theta = theta,
                      Ynex = abundanceIndex[year],
                      Cnex = capt[year],
                      sdPropag = sdPropag,
                      min=capt[year-1],
                      method=method)
  mth <- mtheta(x=X[inds_rspl,year-1],xprime= X[,year],
                theta= theta, capt=capt[year-1])
  gfact <- g(X[,year], abundanceIndex[year], theta)
  pfact <- dPropag(y = X[,year],x= X[inds_rspl,year-1],theta = theta,
                   sdPropag = sdPropag,method =method,
                   Ynex = abundanceIndex[year],Cnex = capt[year],
                   min=capt[year-1])
  weightX[,year]=(mth*gfact/pfact)
  weightX[,year] <-   weightX[,year]/sum(  weightX[,year])
  points(rep(year,G),X[,year],cex=10*weightX[,year],col="darkgreen")
}

plot(density(X[,5]))
lines(density(X[,5],weights = weightX[,5]),col="red")
points(x=X[,5],y=rep(0,G),pch=20,cex=weightX[,5]*10)
