## clean history
rm(list=ls())

## library
library(MASS)
library(mvtnorm)
library(LaplacesDemon)


## working path
dir_car <- "~/Desktop/Dissertation/ODFL/Data Application Revision/Cardiotoxicity "
## read data
df <- read.csv(paste(dir_car,"final_data.csv", sep = '/'), header = TRUE)
# df$intersection <- df$Sex*df$Dose
## split the data by ID
dnew <- split(df, f = df$ID)
#dnew <- dnew[-c(30,78)]
head(dnew)
## basic information if list
n <- length(dnew) # n = 91
mi <- rep(0, n)
for(i in 1:n){
  mi[i] <- nrow(dnew[[i]])
}
N <- sum(mi) # N = 309
Z <- list()
Z_all <- list()
Y <- list()
G_all <- list()
U_all <- list()
Z_else <- list()
for(i in 1:n){
  Z[[i]] <- matrix(dnew[[i]][,7], nrow = nrow(dnew[[i]]), ncol = 1)
  Y[[i]] <- matrix(dnew[[i]][,2], nrow = nrow(dnew[[i]]), ncol = 1)
  len <- nrow(Z[[i]])
  z <- unlist(Z[[i]])[,1]
  if(len > 1){
    u <- matrix(rep(0, (len-1)), nrow = (len-1), ncol = 1)
    for(j in 1:(len-1)){
      u[j,] <- Z[[i]][j+1, ] - Z[[i]][j,]
    }
  }
  U_all[[i]] <- u
  Z_all[[i]] <- matrix(c(rep(1, len), z), nrow = len, ncol = 2)
  Z_else[[i]] <- as.matrix(dnew[[i]][,c(3:5, 8)], nrow = len, ncol = 4)
  G_all[[i]] <- cbind(Z_all[[i]], Z_else[[i]])
}
Zm <- unlist(Z)
Xm <- min(Zm):max(Zm)
plot(c(1), c(1), type = 'n', ylim = c(-6, 6), xlim = c(min(Xm), max(Xm)), xlab = "time point",
     ylab = "response")
cl <- colors()
for (i in 1:n){
  lines(unlist(Z[[i]])[,1], unlist(Y[[i]]), type = 'b', col = cl[i], lwd = 0.5)
}

## initial values for MCMC samples 
Wt <- rep(1.5, n)
Vt <- list()
Sigt <- list()
rho <- 0.5
sigt <- 1
for(i in 1:n){
  mi <- nrow(Z[[i]])
  Vt[[i]] <- matrix(rep(1.3, (mi-1)), nrow = (mi-1), ncol = 1)
  one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
  cov_y <- (1-rho)*diag(1, mi) + rho*one
  Sigt[[i]] <- sigt*cov_y
}
sigg0 <- 1; sigg1 <- 1; sigg2 <- 1
kappat <- 1
zetat <- 1 
gamma0t <- 1; gamma1t <- 1; gamma2t <- 1
rt <- 2
etat <- 2
zt <- 4

thetat <- matrix(rep(1, 2), nrow = 2, ncol = 1)
bett <- matrix(rep(1, 4), nrow = 4, nco = 1) # regression coefficient for covariates (sex, dose and age)
#Thetat <- matrix(rep(0, 6), nrow = 6, ncol = 1)
niter <- 1000 # samples 
nburn <- 1000 # burning 
aiter <- niter + nburn # all
Wtt <- matrix(rep(0, niter*n), nrow = niter, ncol = n)
kappatt <- rep(0, niter)
gamma0tt <- rep(0, niter); gamma1tt <- rep(0, niter); gamma2tt <- rep(0, niter)
zetatt <- rep(0, niter)
sigtt <- rep(0, niter)
rtt <- rep(0, niter)
etatt <- rep(0, niter)
rhott <- rep(0, niter)
thetatt <- matrix(0, niter, 2)
bettt <- matrix(0, niter, 4) 

for(k in 1:aiter){
  #####################################################
  ################## update Wi ########################
  #####################################################
  
  for (i in 1:n) {
    Wci <- rlnorm(1, meanlog = Wt[i] , sdlog = 1)
    sc <- (etat^2/sigt)*solve(Sigt[[i]])
    mi <- nrow(Z[[i]])
    Wti <- Wt[i]
    
    logWc <- diag(log(Wci), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWc <- t(logWc)
    
    logWt <- diag(log(Wti), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)

    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
    R_i <- as.matrix(Y[[i]][,1], mi, 1) - unlist(Z_all[[i]])%*%thetat - Zi%*%bett # residual matrix 
    
    RSSj <- rep(0, mi-1)
    for(j in 1:(mi-1)){
      RSSj[j] <- Vt[[i]][j,] - gamma0t - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,]
    }
    RSS <- sum(RSSj)
    
    numer <- -log(Wci) - (tlogWc%*%sc%*%logWc - (etat/sigt)*tlogWc%*%solve(Sigt[[i]])%*%R_i)
    -((log(Wci) - 1)^2) - (zetat)*((mi-1)*(gamma1t*Wci)^2 - (gamma1t*Wci)*RSS)
    
    denom <- -log(Wti) - (tlogWt%*%sc%*%logWt - (etat/sigt)*tlogWt%*%solve(Sigt[[i]])%*%R_i)
    -((log(Wti) - 1)^2) - (zetat)*((mi-1)*(gamma1t*Wti)^2 - (gamma1t*Wti)*RSS)
    
    #numer <- -sum(log(Wc[i])) - 0.5*sum(tlogWc%*%sc%*%logWc - 2*(etat/sigt)*tlogWc%*%solve(Sigt[[i]])%*%R_i)
    #- (1/kappat)*(sum(log(Wc[i])^2)) - (0.5/zetat)*sum((mi-1)*(gamma1t*Wc[i])^2 - 2*(gamma1t*Wc[i])*RSS)
    
    #denom <- -sum(log(Wt[i])) - 0.5*sum(tlogWt%*%sc%*%logWt - 2*(etat/sigt)*tlogWt%*%solve(Sigt[[i]])%*%R_i)
    #- (1/kappat)*(sum(log(Wt[i])^2)) - (0.5/zetat)*sum((mi-1)*(gamma1t*Wt[i])^2 - 2*(gamma1t*Wt[i])*RSS)
    
    trant <- dlnorm(Wti, meanlog = Wci, sdlog = 1, log = TRUE)
    tranc <- dlnorm(Wci, meanlog = Wti, sdlog = 1, log = TRUE)
    
    logalpha <- numer - denom + trant - tranc 
    
    logu <- log(runif(1, 0, 1))
    
    if(logu<logalpha){
      Wt[i] <- Wci
    }
  }
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################# update V ##########################
  #####################################################
  Vc <- list()
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    Vci <- mvrnorm(1, Vt[[i]], Sigma = diag(1, mi-1))
    Vci <- as.matrix(Vci, ncol = 1, nrow = mi-1)
    Vc[[i]] <- Vci
  }
  
  for(i in 1:n){
    mi_u <- nrow(Z[[i]]) - 1
    numerV <- 0
    denomV <- 0
    
    for(j in 1:mi_u){
      numerVj <- Vc[[i]][j,]*unlist(U_all[[i]][j,]) - (0.5/zetat)*(Vc[[i]][j,] - gamma0t - gamma1t*Wt[i] - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,])^2 - 
        rt*log(rt + exp(Vc[[i]][j,])) - unlist(U_all[[i]][j,])*log(rt + exp(Vc[[i]][j,]))
      denomVj <- Vt[[i]][j,]*unlist(U_all[[i]][j,]) - (0.5/zetat)*(Vt[[i]][j,] - gamma0t - gamma1t*Wt[i] - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,])^2 - 
        rt*log(rt + exp(Vt[[i]][j,])) - unlist(U_all[[i]][j,])*log(rt + exp(Vt[[i]][j,]))
      
      numerV <- numerV + numerVj
      denomV <- denomV + denomVj
    }
    
    logalphaV <- numerV - denomV
    loguV <- log(runif(1, 0, 1))
    if(loguV<logalphaV){
      Vt[[i]] <- Vc[[i]]
    }
  }
  
   
  
  #####################################################
  ################# update kappa ######################
  #####################################################
  Ak <- 1 # hyper parameter
  Bk <- 1 # hyper parameter
  Akstar <- Ak + n/2
  Bkstar <- Bk + sum((log(Wt[i]))^2)
  kappat <- rinvgamma(1, shape=Akstar, scale = Bkstar)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ############ update theta ###########################
  #####################################################
  DF <- matrix(rep(0, 2), 2, 2)
  EF <- matrix(rep(0, 2), 2, 1)
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
    Hyi <- as.matrix(Y[[i]][,1], mi, 1) - etat*logWt - Zi%*%bett
    DFi <- as.matrix((1/sigt)*t(unlist(Z_all[[i]]))%*%solve(Sigt[[i]])%*%unlist(Z_all[[i]]), 
                     nrow = mi, ncol = mi)
    EFi <- as.matrix((1/sigt)*t(unlist(Z_all[[i]]))%*%solve(Sigt[[i]])%*%Hyi, 
                     nrow = mi, ncol = 1)
    
    DF <- DF + DFi
    EF <- EF + EFi
  }
  S <- solve(DF)
  M <- S%*%EF
  thetat <- matrix(mvrnorm(1, M, S), 2, 1)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ############ update beta ############################
  #####################################################
  DB <- matrix(rep(0, 16), 4, 4)
  EB <- matrix(rep(0, 4), 4, 1)
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
    Hybi <- as.matrix(Y[[i]][,1], mi, 1) - unlist(Z_all[[i]])%*%thetat - etat*logWt
    DBi <- as.matrix((1/sigt)*t(Zi)%*%solve(Sigt[[i]])%*%Zi, 
                     nrow = mi, ncol = mi)
    EBi <- as.matrix((1/sigt)*t(Zi)%*%solve(Sigt[[i]])%*%Hybi, 
                     nrow = mi, ncol = 1)
    
    DB <- DB + DBi
    EB <- EB + EBi
  }
  Sb <- solve(DB)
  Mb <- Sb%*%EB
  bett <- matrix(mvrnorm(1, Mb, Sb), 4, 1)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ############## update sigma2 ########################
  #####################################################
  As <- 1 # hyper parameter
  Bs <- 1 # hyper parameter
  Asstar <- As + N/2
  RSS <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
    Ri <- as.matrix(Y[[i]][,1], mi, 1) - unlist(Z_all[[i]])%*%thetat - etat*logWt - Zi%*%bett
    RSSi <- t(Ri)%*%solve(Sigt[[i]])%*%Ri
    
    RSS <- RSS + RSSi
    
  }
  Bsstar <- 0.5*RSS + Bs
  sigt <- rinvgamma(1, shape = Asstar, scale = Bsstar)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ############### update eat ##########################
  #####################################################
  iDl <- 0
  El <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    #    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    #    tlogWt <- t(logWt)
    ones <- matrix(1,nrow=mi,ncol=1)
    onet <- t(ones)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
    R_i <- as.matrix(Y[[i]][,1], mi, 1) - unlist(Z_all[[i]])%*%thetat - Zi%*%bett
    
    #    iDli <- (1/sigt)*tlogWt%*%solve(Sigt[[i]])%*%logWt
    #    Eli <- (1/sigt)*tlogWt%*%R_i
    iDli <- Wt[i]^2*onet%*%solve(Sigt[[i]])%*%ones
    Eli <- Wt[i]*onet%*%R_i
    
    iDl <- iDl + iDli 
    El <- El + Eli
  }
  Dl <- 1/iDl
  
  etat <- rnorm(1, Dl*El, sqrt(Dl*sigt))
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################ update zeta2 #######################
  #####################################################
  Az <- 100 # hyper parameter
  Bz <- 100 # hyper parameter
  RSSv <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      RSSvi <- (Vt[[i]][j,] - gamma0t - gamma1t*Wt[i] - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,])^2
    }
    RSSv <- RSSv + RSSvi
  }
  Azstar <- (N - n)/2 + Az
  Bzstar <- RSSv/2 + Bz
  zetat <- rinvgamma(1, shape=Azstar, scale=Bzstar)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################# update gamma0 #####################
  #####################################################
  RSS0 <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      RSS0i <- Vt[[i]][j,] - gamma1t*Wt[i] - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,]
    }
    RSS0 <- RSS0 + RSS0i
  }
  iD0 <- (2*n/zetat + 1/sigg0)
  E0 <- RSS0/zetat
  D0 <- 1/iD0
  gamma0t <- rnorm(1, D0*E0, sqrt(D0))
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################# update gamma1 #####################
  #####################################################
  E1 <- 0
  iD1 <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      RSSg1i <- Vt[[i]][j,] - gamma0t - gamma2t*(as.matrix(Y[[i]][,1], mi, 1))[j,]
    }
    E1i <- (Wt[i]*RSSg1i)/zetat
    E1 <- E1 + E1i
    
    iD1i <- (2*(Wt[i]^2))/zetat + 1/sigg1
    iD1 <- iD1 + iD1i
  }
  D1 <- 1/iD1
  gamma1t <- rnorm(1, D1*E1, sqrt(D1))
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################# update gamma2 #####################
  #####################################################
  E2 <- 0
  iD2 <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      RSS2i <- ((as.matrix(Y[[i]][,1], mi, 1))[j,])*(Vt[[i]][j,] - gamma0t - gamma1t*Wt[i])
      iD2i <- (((as.matrix(Y[[i]][,1], mi, 1))[j,])^2)/zetat + 1/sigg2
    }
    E2i <- RSS2i/zetat
    E2 <- E2 + E2i
    
    iD2 <- iD2 + iD2i
  }
  D2 <- 1/iD2
  gamma2t <- rnorm(1, D2*E2, sqrt(D2))
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################ update r ###########################
  #####################################################
  Ar <- 50 # hyper parameter
  Br <- 50 # hyper parameter
  rc <- rgamma(1, rt^2, rate=rt)
  
  numerRj <- rep(0, (mi-1))
  denomRj <- rep(0, (mi-1))
  numerRi <- rep(0, n)
  denomRi <- rep(0, n)
  
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      numerRj[j] <- log(gamma(unlist(U_all[[i]][j,]) + rc)) + (Ar + rc-1)*log(rc) + (-rc*Br)-
        (log(gamma(rc)) + rc*log(rc+exp(Vt[[i]][j,])) + unlist(U_all[[i]][j,])*log(rc+exp(Vt[[i]][j,])))
      denomRj[j] <- log(gamma(unlist(U_all[[i]][j,]) + rt)) + (Ar + rt-1)*log(rt) + (-rt*Br)-
        (log(gamma(rt)) + rt*log(rt+exp(Vt[[i]][j,])) + unlist(U_all[[i]][j,])*log(rt+exp(Vt[[i]][j,])))
      
    }
    numerRi[i] <- sum(numerRj)
    denomRi[i] <- sum(denomRj)
  }
  numerR <- sum(numerRi)
  denomR <- sum(denomRi)
  
  trantcR <- dgamma(rt, shape = rc^2, rate = rt, log = TRUE)
  tranctR <- dgamma(rc, shape = rt^2, rate = rc, log = TRUE)
  
  logalphaR <- numerR-denomR + trantcR - tranctR
  loguR <- log(runif(1, 0, 1))
  if(loguR<logalphaR){
    rt <- rc
  }
  
 
  #####################################################
  ############### update rho (Sig) ####################
  #####################################################
  zc <- rnorm(1, zt, 0.1)
  rhoc <- (exp(2*zc)-1)/(exp(2*zc)+1)
  rhot <- (exp(2*zt)-1)/(exp(2*zt)+1)
  Sigc <- list()
  Sigt <- list()
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
    Sigci <- (1-rhoc)*diag(1, mi) + rhoc*one
    Sigc[[i]] <- Sigci
    Sigti <- (1-rhot)*diag(1, mi) + rhot*one
    Sigt[[i]] <- Sigti
  }
  
  # transition
  logzc <- zc - 2*log(exp(zc)+1)
  logzt <- zt - 2*log(exp(zt)+1)
  numerir <- rep(0, n)
  denomir <- rep(0, n)

  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)

    Rii <- as.matrix(Y[[i]][,1], mi, 1) - unlist(Z_all[[i]])%*%thetat - etat*log(Wt[i]) - Zi%*%bett
    numerir[i] <- (-0.5)*log(det(Sigc[[i]]))-(1/2*sigt)*t(Rii)%*%solve(Sigc[[i]])%*%Rii
    denomir[i] <- (-0.5)*log(det(Sigt[[i]]))-(1/2*sigt)*t(Rii)%*%solve(Sigt[[i]])%*%Rii

  }
  numerr <- sum(numerir, na.rm = TRUE)
  denomr <- sum(denomir, na.rm = TRUE)
  
  logalpharh <- numerr-denomr+logzt - logzc
  logu <- log(runif(1, 0, 1))
  if(logu<logalpharh){
    zt <- zc
  }
  rhot <- (exp(2*zt)-1)/(exp(2*zt)+1)

  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ########### collecte MCMC samples ###################
  #####################################################
  if(k>nburn){
    Wtt[k-nburn,] <- Wt
    kappatt[k-nburn] <- kappat 
    thetatt[k-nburn,] <- thetat
    bettt[k-nburn,] <- bett
    etatt[k-nburn] <- etat
    sigtt[k-nburn] <- sigt
    zetatt[k-nburn] <- zetat
    gamma0tt[k-nburn] <- gamma0t
    gamma1tt[k-nburn] <- gamma1t
    gamma2tt[k-nburn] <- gamma2t
    rtt[k-nburn] <- rt
    rhott[k-nburn] <- rhot
  }
}

W_bp <- matrix(rep(0, n), nrow = 1, ncol = n)
for(i in 1:n){
  W_bp[,i] <- mean(Wtt[,i])
}

(eta_bp <- mean(etatt)) # posterior mean
(sd(etatt)) # standard deviation 
quantile(etatt, probs = c(0.025, 0.975))

gamma1 <- mean(gamma1tt) # posterior mean 
sd(gamma1tt) # standard deviation
######### history Yij-1 gamma2 ##############################################
gamma2 <- mean(gamma2tt) # posterior mean 
quantile(gamma2tt, probs = c(0.025, 0.975))
sd(gamma2tt) # standard deviation

(rho <- mean(rhott)) 
sd(rhott)
quantile(rhott, probs=c(0.025, 0.975))

apply(bettt, 2, mean)
apply(bettt, 2, sd)
apply(bettt, 2, quantile,  probs = c(0.025, 0.975), na.rm = TRUE)

apply(thetatt, 2, mean)
apply(thetatt, 2, sd)
apply(thetatt, 2, quantile, probs=c(0.025, 0.975))
theta <- apply(thetatt, 2, mean)

#############################################################################
################### plot data after regression ##############################
#############################################################################
plot(c(1), c(1), type = 'n', ylim = c(-6, 6), xlim = c(min(Xm), max(Xm)), xlab = "time point",
     ylab = "response")
cl <- colors()
for (i in 1:n){
  lines(unlist(Z[[i]])[,1], unlist(Y[[i]]), type = 'l', col = 'gray', lwd = 0.5)
}
################# add lines #################################################
################# estimated lines ###########################################
Z_cov <- cbind(0.5177994, 0.4919094, 6.467024, 0.2977346)
Wmean <- mean(W_bp)
yhat <- theta[2]*Xm + theta[1] + eta_bp*log(Wmean)
lines(Xm, yhat, lwd = 2, lty = 1, col = 'blue')

cpo <- rep(0, n)
logcpo <- rep(0, n)
for(i in 1:n){
  mi <- nrow(Z[[i]])
  mu_overall <- 0
  Sigm_i <- matrix(rep(0, mi*mi), mi, mi)
  for(kk in 1:1000){
    mu_overall <- unlist(Z[[i]])%*%thetatt[kk,2] + etatt[kk]*log(Wtt[kk,i]) + c(Z_cov%*%bettt[kk,])
    one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
    Sigm_i <- sigtt[kk]*((1-rhott[kk])*diag(1, mi) + rhott[kk]*one)
  }
  
  Yii <- matrix(unlist(Y[[i]]), ncol = mi)
  densi <- dmvnorm(Yii, mean = mu_overall, sigma = Sigm_i)
  
  cpo[i] <- 1/(mean(1/densi))
  logcpo[i] <- log(cpo[i])
}
## small LPML suggests possible outliers, high leverage
## thus, higher LPML suggest better model 
lpml <- mean(logcpo) 

Rst.Car.NB.BP <- list()
Rst.Car.NB.BP[[1]] <- W_bp
Rst.Car.NB.BP[[2]] <- eta_bp
Rst.Car.NB.BP[[3]] <- sd(etatt)
Rst.Car.NB.BP[[4]] <- gamma2
Rst.Car.NB.BP[[5]] <- sd(gamma2tt)
Rst.Car.NB.BP[[6]] <- mean(rhott)
Rst.Car.NB.BP[[7]] <- sd(rhott)
Rst.Car.NB.BP[[8]] <- apply(thetatt, 2, mean)
Rst.Car.NB.BP[[9]] <- apply(thetatt, 2, sd)
Rst.Car.NB.BP[[10]] <- apply(bettt, 2, mean)
Rst.Car.NB.BP[[11]] <- apply(bettt, 2, sd)
Rst.Car.NB.BP[[12]] <- lpml
names(Rst.Car.NB.BP) <- c("mean of W bp", "eta hat", "eta sd", 
                          "gamma2 hat", "sd gamma2", "rho hat", 
                          "sd rho", "theta hat", "theta sd", 
                          "beta hat", "beta sd", "lpml")

save(Rst.Car.NB.BP, file = paste(dir_car, "Rst.Car.NB.BP.RData", sep = '/'))


