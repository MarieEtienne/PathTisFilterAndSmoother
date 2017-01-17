## ------------------------------------------------------------------------

# Population Dynamic model simulation ----------------------------------
rm(list=ls())
set.seed(11)
nYear <- 30
biomass <- rep(NA, nYear)
capt <- rep(NA, nYear)
K <- 300 ## carrying capacity
r <- 1.02 ## recruitement rate
q <- 0.5
epsilon <- 10

sigmaBiomass <- 0.3
biomass[1] <- rnorm(1, mean=0.9*K, sd=K/10)
capt[1] <- runif(1,0.1,0.5)*biomass[1]
for( year in 2:nYear){
  Bk = biomass[year-1]
  capt[year-1] <- runif(1,0.1,0.5)*Bk
  #Captures are a random fraction of biomass
  innov <- rnorm(1, mean=-sigmaBiomass^2/2,sd=sigmaBiomass)
  biomass[year] <- (Bk+r*Bk*(1- Bk/K)-capt[year-1])*exp(innov)
  capt[year] <- runif(1,0.1,0.5)*biomass[year]
}
obs_error <- rnorm(nYear, mean=-sigmaBiomass^2/2, sd=sigmaBiomass)
abundanceIndex <- q * biomass*exp(obs_error)
#Calcul de la fenêtre des possible pour les Xs
sqrtD <- sqrt((1+r)^2-4*capt*r/K)
M0s <- (-(1+r)-sqrtD)/(-2*r/K)
m0s <- (-(1+r)+sqrtD)/(-2*r/K)
bounds <- cbind(m0s,M0s)

## ------------------------------------------------------------------------
par(mfrow=c(1,1))
plot(1:nYear,biomass, ylim=range(c(0,range(biomass),range(bounds))),
     main='Evolution de la biomasse',xlab="Année",type="b",pch=20,
     ylab="Masse",bty="n")
lines(1:(nYear),capt,type="b",col="red",pch=20)
legend("topright",pch=20,col=c("black","red"),
       legend = c("Biomasse","Captures"),
       bty="n",horiz = T)

lines(1:nYear, abundanceIndex, type='b',col="blue",pch=18,
      ylim=range(c(0,range(abundanceIndex))))
axis(side=4,at=pretty(range(abundanceIndex),n = 2),col = "blue",
     col.ticks = "blue",labels=F ,lwd.ticks=0)
mtext(side=4,text = "Indice d'abondance",col = "blue",at=150,line=1)
mtext(side=4,text = pretty(range(abundanceIndex),n=2),col = "blue",
      at=pretty(range(abundanceIndex),n=2))
matplot(bounds,lty=2,col="darkgrey",add=T,type="l")

## ------------------------------------------------------------------------
mtheta<- function(x_old,x_new, capt_old,theta){
  r =theta$r; K=theta$K; sigmaBiomass = theta$sigmaBiomass;
  mn_innov =-sigmaBiomass^2/2
  dnorm(log(x_new), 
        mean=log(x_old+r*x_old*(1-x_old/K)- capt_old)-mn_innov, 
        sd=sigmaBiomass)
}
g<-function(x,y,theta){
  q =theta$q;sigmaBiomass=theta$sigmaBiomass;mn_innov =-sigmaBiomass^2/2
  dnorm(log(y), mean=log(q*x)-mn_innov,sd=sigmaBiomass)
}

## ------------------------------------------------------------------------
#Fonction de génération
rf0 <- function(G, theta,Y0=NA,C0=NA){
    #Extraction des paramètres nécessaires
    r =theta$r; sigmaBiomass = theta$sigmaBiomass;q=theta$q;
    mn_innov = -sigmaBiomass^2/2;#Pourrait être défini autrement
    Delta <- (1+r)^2-4*C0*r/K
    sqrtD <- sqrt(Delta)
    M0 <- (-(1+r)-sqrtD)/(-2*r/K)
    m0 <- (-(1+r)+sqrtD)/(-2*r/K)
    X0 <- sapply(1:G,function(i){
      ok <- F
      while(!ok){#On simule jusqu'à tomber dans la zone 
        #acceptable au vu des captures
        cand <- exp(rnorm(1,log(Y0)-log(q)-mn_innov,
                          sigmaBiomass))
        ok <- (cand < M0) & (cand>m0)
      }
      return(cand)
    })
    return(X0)
}
df0 <- function(x0, theta,Y0=NA,C0=NA){
  #Extraction des paramètres nécessaires
  r =theta$r; sigmaBiomass = theta$sigmaBiomass;q=theta$q;
  mn_innov = -sigmaBiomass^2/2;#Pourrait être défini autrement
  Delta <- (1+r)^2-4*C0*r/K
  sqrtD <- sqrt(Delta)
  M0 <- (-(1+r)-sqrtD)/(-2*r/K)
  m0 <- (-(1+r)+sqrtD)/(-2*r/K)
  mn_ref <- log(Y0)-log(q)-mn_innov
  #Définition de la constante de normalisation
  norm_cst <- (pnorm(log(M0),mn_ref,sigmaBiomass)-
                 pnorm(log(m0),mn_ref,sigmaBiomass))
  dx <- dnorm(log(x0),mn_ref,sigmaBiomass)/norm_cst
  return(dx)
}

## ------------------------------------------------------------------------
rPropag <- function(x_old,theta,Y_new=NA,C_new=NA){
  sdRW = theta$sdRW; sigmaBiomass = theta$sigmaBiomass;
  mn_innov = -sigmaBiomass^2/2; q =theta$q;r=theta$r; 
  K = theta$K;
  sqrtD <- sqrt((1+r)^2-4*C_new*r/K)
  M_new <- (-(1+r)-sqrtD)/(-2*r/K)
  m_new <- (-(1+r)+sqrtD)/(-2*r/K)
  v_gau  <- (sdRW^2*sigmaBiomass^2)/(sdRW^2+sigmaBiomass^2)
  mn2    <- log(Y_new)-log(q)-mn_innov
  mn_gau <- v_gau*(log(x_old)/(sdRW^2)+mn2/(sigmaBiomass^2))
  Xnex <- sapply(mn_gau,function(mn){
    ok <- F
    while(!ok){
      cand <- exp(rnorm(1,mn,sqrt(v_gau)))
      ok <- (cand > m_new) & (cand < M_new)
    }
    return(cand)
  })
  return(Xnex)
}
dPropag <- function(x_new, x_old,theta,Y_new=NA,C_new=NA){
  sdRW = theta$sdRW; sigmaBiomass = theta$sigmaBiomass;
  mn_innov = -sigmaBiomass^2/2; q =theta$q;r=theta$r; 
  K = theta$K;
  sqrtD <- sqrt((1+r)^2-4*C_new*r/K)
  M_new <- (-(1+r)-sqrtD)/(-2*r/K)
  m_new <- (-(1+r)+sqrtD)/(-2*r/K)
  v_gau  <- (sdRW^2*sigmaBiomass^2)/(sdRW^2+sigmaBiomass^2)
  mn2    <- log(Y_new)-log(q)-mn_innov
  mn_gau <- v_gau*(log(x_old)/(sdRW^2)+mn2/(sigmaBiomass^2))
  dnex <- sapply(1:length(x_new),function(i){
    mn = mn_gau[i]
    (dnorm(log(x_new[i]),
           mn,
           sqrt(v_gau))
    /(pnorm(log(M_new),mn,sqrt(v_gau))
      -pnorm(log(m_new),mn,sqrt(v_gau))))
    #normalized by the area of possible zone
  })
  return(dnex)
}

## ----echo= -c(1:2)-------------------------------------------------------
#Initialisation
sdRW <- 10
theta=list(q=q, K=K, sigmaBiomass=sigmaBiomass, r=r,sdRW=sdRW)
G <- 1000#50 particules
X <- matrix(NA, ncol=nYear, nrow=G)
X[,1] <- rf0(G=G,theta= theta,
             Y0 = abundanceIndex[1],C0 = capt[1] )
weightX <- matrix(NA, ncol=nYear, nrow=G)
weightX[,1] <- df0(x0=X[,1],theta = theta,
                   Y0 = abundanceIndex[1],C0 = capt[1])
weightX[,1] <- weightX[,1]/sum(weightX[,1])
for(year in 2:nYear){
  X[,year] <- rPropag(x_old = X[,year-1],
                      theta = theta,
                      Y_new = abundanceIndex[year],
                      C_new = capt[year])
  mth <- mtheta(x_old =X[,year-1],
                x_new= X[,year],
                theta= theta,
                capt_old=capt[year-1])
  gfact <- g(X[,year], abundanceIndex[year], theta)
  pfact <- dPropag(x_new = X[,year],x_old= X[,year-1],
                   theta = theta, Y_new = abundanceIndex[year],
                   C_new = capt[year])
  weightX[,year]=weightX[,year-1]*(mth*gfact/pfact)
  weightX[,year] <-   weightX[,year]/sum(  weightX[,year])
}

## ------------------------------------------------------------------------
matplot(t(weightX),type="b",pch=20,xlab="Indice k",ylab="Poids",
        main="Evolution des poids de filtrage\n cas sans resampling")

## ----echo= -c(1:2)-------------------------------------------------------
plot(1:nYear,biomass, ylim=range(c(range(biomass),range(bounds))),
     main='Evolution de la biomasse',xlab="Année",type="n",pch=20,
     ylab="Masse",bty="n")
matplot(bounds,lty=2,col="darkgrey",add=T,type="l")
#Initialisation
sdRW <- 10
theta=list(q=q, K=K, sigmaBiomass=sigmaBiomass, r=r,sdRW=sdRW)
G <- 1000#50 particules
X <- matrix(NA, ncol=nYear, nrow=G)
X[,1] <- rf0(G=G,theta= theta,
             Y0 = abundanceIndex[1],C0 = capt[1] )
weightX <- matrix(NA, ncol=nYear, nrow=G)
weightX[,1] <- df0(x0=X[,1],theta = theta,
                   Y0 = abundanceIndex[1],C0 = capt[1])
weightX[,1] <- weightX[,1]/sum(weightX[,1])
points(rep(1,G),X[,1],cex=weightX[,1]/max(weightX[,1]),
       col="darkgreen",pch=20)

for(year in 2:nYear){
  inds_rspl <- sample(1:G,size=G,replace=T,prob = weightX[,year-1])
  X[,year] <- rPropag(x_old = X[inds_rspl,year-1],
                      theta = theta,
                      Y_new = abundanceIndex[year],
                      C_new = capt[year])
  mth <- mtheta(x_old =X[inds_rspl,year-1],
                x_new= X[,year],
                theta= theta,
                capt_old=capt[year-1])
  gfact <- g(X[,year], abundanceIndex[year], theta)
  pfact <- dPropag(x_new = X[,year],x_old= X[inds_rspl,year-1],
                   theta = theta, Y_new = abundanceIndex[year],
                   C_new = capt[year])
  weightX[,year]=(mth*gfact/pfact)
  weightX[,year] <-   weightX[,year]/sum(  weightX[,year])
  points(rep(year,G), X[,year],col="darkgreen",
         cex=weightX[,year]/max(weightX[,year]),pch=20)
}
lines(biomass,type="b",pch=20)

## ------------------------------------------------------------------------
matplot(t(weightX),type="b",pch=20,xlab="Indice k",ylab="Poids",
        main="Evolution des poids de filtrage\n cas sans resampling")

## ------------------------------------------------------------------------
# par(mfrow=c(4,4))
# sapply(1:nYear,function(i){
#   plot(density(X[,i],weights = weightX[,i]),main=i)
#   abline(v=biomass[i])
# })
# par(mfrow=c(4,4))
# sapply(1:nYear,function(i){
#   plot(X[,i],weightX[,i],main=i)
#   abline(v=biomass[i])
# })
# par(mfrow=c(4,4))
# sapply(1:nYear,function(i){
#   hist(weightX[,i],main=i)
#   abline(v=biomass[i])
# })

