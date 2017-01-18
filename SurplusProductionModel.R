## ----modelsimulation-----------------------------------------------------
# Population Dynamic model simulation ----------------------------------
rm(list=ls()); set.seed(15);par(mfrow=c(1,1))
nYear <- 30
biomass <- rep(NA, nYear);capt <- rep(NA, nYear)
K <- 300 ## Capacité d'accueil
r <- 1.02 ## Taux de recrutement
q <- 0.5 #Capturabilité
sigmaBiomass <- 0.3
biomass[1] <- rnorm(1, mean=0.9*K, sd=K/10)
capt[1] <- runif(1,0.1,0.5)*biomass[1]
for( year in 2:nYear){
 Bk = biomass[year-1]
 capt[year-1] <- runif(1,0.1,0.5)*Bk
 #Les captures sont une fraction aléatoire de la biomasse
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

## ----plotmodel,echo=FALSE,fig.height=5,fig.width=6-----------------------
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

## ----fonction_mg---------------------------------------------------------
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

## ----f0------------------------------------------------------------------
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

## ----Propag--------------------------------------------------------------
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
   (dnorm(log(x_new[i]),mn,sqrt(v_gau))
   /(pnorm(log(M_new),mn,sqrt(v_gau))-pnorm(log(m_new),mn,sqrt(v_gau))))
   #normalized by the area of possible zone
 })
 return(dnex)
}

## ----first_filt----------------------------------------------------------
#Initialisation
sdRW <- 10#Ecart type de la marche aleatoire
theta=list(q=q, K=K, sigmaBiomass=sigmaBiomass,
          r=r,sdRW=sdRW)#Paramètres pour la génération de particules
#et calculs de poids
G <- 30#Nombre de particules
Xs1 <- matrix(NA, ncol=nYear, nrow=G)#Particules
Xs1[,1] <- rf0(G=G,theta= theta,
            Y0 = abundanceIndex[1],C0 = capt[1] )
ws1 <- matrix(NA, ncol=nYear, nrow=G)#Poids
ws1[,1] <- df0(x0=Xs1[,1],theta = theta,
             Y0 = abundanceIndex[1],C0 = capt[1])
ws1[,1] <-  ws1[,1]/sum(ws1[,1])
for(year in 2:nYear){
 Xs1[,year] <- rPropag(x_old = Xs1[,year-1],theta = theta,
                       Y_new = abundanceIndex[year],
                       C_new = capt[year])#Génération de nvelles particules
 mth <- mtheta(x_old =Xs1[,year-1],x_new= Xs1[,year],
               theta= theta,capt_old=capt[year-1])
 #densité de transition du modele
 gfact <- g(Xs1[,year], abundanceIndex[year], theta)#Dnsité d'observation
 pfact <- dPropag(x_new = Xs1[,year],x_old= Xs1[,year-1],
                  theta = theta, Y_new = abundanceIndex[year],
                  C_new = capt[year])#Densité de la loi de proposition
 ws1[,year]=ws1[,year-1]*(mth*gfact/pfact)#Actualisation des poids
 ws1[,year] <-   ws1[,year]/sum(ws1[,year])#Normalisation des poids
}

## ----com_density,eval=FALSE----------------------------------------------
## density(Xs1[,k],weights = ws1[,k])

## ----first_dens_filt,echo=F,fig.height=5,fig.width=6---------------------
k=nYear-5#Indice à représenter
plot(density(Xs1[,k],weights = ws1[,k]),
    main=bquote(paste('Loi de filtrage de '*'X'[.(k)])),
    xlab="Biomasse",ylab="Densité")
abline(v=biomass[k],lty=2)
points(Xs1[,k],rep(0,G) ,col="green",
       cex=ws1[,k]/max(ws1[,k]),pch=20)

## ----first_ws,echo=FALSE,fig.height=5,fig.width=6------------------------
matplot(t(ws1),type="b",pch=20,xlab="Indice k",ylab="Poids",
       main="Evolution des poids de filtrage\n Cas sans resampling",
       lty=1)

## ----dec_filt_run,echo=F-------------------------------------------------
Geneal <- matrix(NA, ncol=nYear, nrow=G)
Geneal[,1] <- 1:G
Xs2 <- matrix(NA, ncol=nYear, nrow=G)
Xs2[,1] <- rf0(G=G,theta= theta,
            Y0 = abundanceIndex[1],C0 = capt[1] )
ws2 <- matrix(NA, ncol=nYear, nrow=G)
ws2[,1] <- df0(x0=Xs2[,1],theta = theta,
                  Y0 = abundanceIndex[1],C0 = capt[1])
ws2[,1] <- ws2[,1]/sum(ws2[,1])
for(year in 2:nYear){
 inds_rspl <- sample(1:G,size=G,replace=T,
                     prob = ws2[,year-1])#Rééchantillonnage
 Gnew <- Geneal[inds_rspl,]
 Gnew[,nYear] <- 1:G
 Gnew[,year-1] <- inds_rspl
 Geneal <- Gnew
 Xs2[,year] <- rPropag(x_old = Xs2[inds_rspl,year-1],
                     theta = theta,
                     Y_new = abundanceIndex[year],
                     C_new = capt[year])
 mth <- mtheta(x_old =Xs2[inds_rspl,year-1],
               x_new= Xs2[,year],
               theta= theta,
               capt_old=capt[year-1])
 gfact <- g(Xs2[,year], abundanceIndex[year], theta)
 pfact <- dPropag(x_new = Xs2[,year],x_old= Xs2[inds_rspl,year-1],
                  theta = theta, Y_new = abundanceIndex[year],
                  C_new = capt[year])
 ws2[,year]=(mth*gfact/pfact)
 ws2[,year] <-   ws2[,year]/sum(ws2[,year])
}

## ----sec_filt_show,eval=FALSE--------------------------------------------
## #Initialisation, identique à précédemment
## for(year in 2:nYear){
##  inds_rspl <- sample(1:G,size=G,replace=T,
##                      prob = ws2[,year-1])#Rééchantillonnage
##  Xs2[,year] <- rPropag(x_old = Xs2[inds_rspl,year-1],
##                      theta = theta,
##                      Y_new = abundanceIndex[year],
##                      C_new = capt[year])
##  mth <- mtheta(x_old =Xs2[inds_rspl,year-1],
##                x_new= Xs2[,year],
##                theta= theta,
##                capt_old=capt[year-1])
##  gfact <- g(Xs2[,year], abundanceIndex[year], theta)
##  pfact <- dPropag(x_new = Xs2[,year],x_old= Xs2[inds_rspl,year-1],
##                   theta = theta, Y_new = abundanceIndex[year],
##                   C_new = capt[year])
##  ws2[,year]=(mth*gfact/pfact)
##  ws2[,year] <-   ws2[,year]/sum(ws2[,year])
## }

## ----sec_dens,echo=F,fig.height=5,fig.width=6----------------------------
plot(density(Xs2[,k],weights = ws2[,k]),
    main=bquote(paste('Loi de filtrage de '*'X'[.(k)])),
    xlab="Biomasse",ylab="Densité")
abline(v=biomass[k],lty=2)
points(Xs2[,k],rep(0,G) ,col="green",
       cex=ws2[,k]/max(ws2[,k]),pch=20)

## ----sec_ws,echo=F,fig.height=5,fig.width=6------------------------------
matplot(t(ws2),type="b",pch=20,xlab="Indice k",ylab="Poids",
       main="Evolution des poids de filtrage\n Cas avec resampling")

## ----sec_filt_plot,echo=F,fig.height=5,fig.width=6,message=FALSE,fig.show=T,results="hide"----
plot(1:nYear,biomass, ylim=range(c(range(biomass),range(bounds))),
    main='Dynamique des particules simulées',
    xlab="Année",type="n",pch=20, ylab="Masse",bty="n")
matplot(bounds,lty=2,col="darkgrey",add=T,type="l")
sapply(1:nYear,function(year){
  points(rep(year,G), Xs2[,year],col="green",
         cex=ws2[,year]/max(ws2[,year]),pch=20)
})
lines(biomass,type="b",pch=20)

## ----traj_first,echo=F,fig.height=5,fig.width=6--------------------------
cols=ws1[,nYear]/max(ws1[,nYear])
matplot(t(Xs1),type="l",lwd=cols,lty=1,
      xlab="Année",ylab="Masse",
      main="Trajectoires particulaires\n Cas sans rééchantillonnage")
rect(par("usr")[1],par("usr")[3],
    par("usr")[2],par("usr")[4],col="lightgray")
matplot(t(Xs1)[,order(cols)],add =T,
       type="l",lty=1,col=gray(1-cols[order(cols)]))

## ----geneal,echo=F,fig.width=6,fig.height=6------------------------------
par(mfrow=c(4,1),mar=c(1,1,1,1))
matplot((nYear-2):nYear,ylab="",xlab="",xaxt="n",yaxt="n",col="black",
       t(Geneal)[(nYear-2):nYear,],lty=1,type="b",pch=20,xlim=c(1,nYear))
matplot((nYear-10):nYear,ylab="",xlab="",xaxt="n",yaxt="n",col="black",
       t(Geneal)[(nYear-10):nYear,],lty=1,type="b",pch=20,xlim=c(1,nYear))
matplot((nYear-20):nYear,ylab="",xlab="",xaxt="n",yaxt="n",col="black",
       t(Geneal)[(nYear-20):nYear,],lty=1,type="b",pch=20,xlim=c(1,nYear))
matplot((1):nYear,ylab="",xlab="",xaxt="n",yaxt="n",col="black",
       t(Geneal)[(1):nYear,],lty=1,type="b",pch=20,xlim=c(1,nYear))

