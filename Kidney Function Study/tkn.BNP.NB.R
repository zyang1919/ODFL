## clean history
rm(list=ls())

## library
library(MASS)
library(mvtnorm)
library(LaplacesDemon)
getK <- function(x)
{
  # This function returns K matrix, and requires unique design points x
  #
  # Example
  # x <- 1:100
  # K <- getK(x)
  
  n <- length(x)
  dx <- diff(x)
  dx2 <- diff(x,lag=2)
  
  Qt <- matrix(rep(0,(n-2)*n),n-2,n)
  for (k in 1:(n-2)){
    Qt[k,k]   <- 1/dx[k]
    Qt[k,k+1] <- -1/dx[k+1]-1/dx[k]
    Qt[k,k+2] <- 1/dx[k+1]
  }
  
  R  <- matrix(rep(0,(n-2)*(n-2)),n-2,n-2)
  R[1,1] <- dx[1]
  R[1,2] <- dx[2]
  for (k in 2:(n-3)){
    R[k,k-1] <- dx[k]
    R[k,k]   <- 2*dx2[k]
    R[k,k+1] <- dx[k+1]
  }
  R[n-2,n-3] <- dx[n-2]
  R[n-2,n-2] <- dx[n-1]
  R <- R/6
  sR <- eigen(R)
  sD <- sR$values
  sD <- sD[sD>10^(-4)]
  ni <- length(sD)
  sS <- sR$vectors
  invR <- sS[,1:ni]%*%solve(diag(sD))%*%t(sS[, 1:ni])
  
  K <- t(Qt)%*%invR%*%Qt
  
  K
}

## working path
dir_tkn <- "~/Desktop/Dissertation/ODFL/Data Application Revision/kidney data"
## read data
df <- read.csv(paste(dir_tkn,"data_final_1207.csv", sep = '/'), header = TRUE)
#df$time <- rank(df$flw, ties.method = "min")
dnew <- split(df, f = df$id_auto)
head(dnew)
## basic information if list
n <- length(dnew) # n = 73
mi <- rep(0, n)
for(i in 1:n){
  mi[i] <- nrow(dnew[[i]])
}
N <- sum(mi) # N = 863
Z <- list()
Z_all <- list()
Y <- list()
G_all <- list()
U_all <- list()
Z_else <- list()
for(i in 1:n){
  Z[[i]] <- matrix(dnew[[i]][,7], nrow = nrow(dnew[[i]]), ncol = 1)
  Y[[i]] <- matrix(dnew[[i]][,6], nrow = nrow(dnew[[i]]), ncol = 1)
  len <- nrow(Z[[i]])
  z <- unlist(Z[[i]])[,1]
  if(len > 1){
    u <- matrix(rep(0, (len-1)), nrow = (len-1), ncol = 1)
    for(j in 1:(len-1)){
      u[j,] <- Z[[i]][j+1, ] - Z[[i]][j,]
    }
  }
  U_all[[i]] <- u
  Z_all[[i]] <- matrix(c(rep(1, len), z, z^2), nrow = len, ncol = 3)
  Z_else[[i]] <- as.matrix(dnew[[i]][,2:3], nrow = len, ncol = 2)
  G_all[[i]] <- cbind(Z_all[[i]], Z_else[[i]])
}
Zm <- unlist(Z)
Xv <- as.vector(Zm)
Xm <- sort(unique(Xv))

## get K matrix 
K <- getK(Xm)


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

ft <- rep(0, length(Xm)) # BNP only
alp  <- 0.01 # BNP only
taut <- alp/sigt # BNP only
bett <- matrix(rep(1, 2), nrow = 2, nco = 1) #  # regression coefficient for covariates (sex, dose and age)
#Thetat <- matrix(rep(0, 6), nrow = 6, ncol = 1)
niter <- 1000 # samples 
nburn <- 1000 # burning 
aiter <- niter + nburn # all
Wtt <- matrix(rep(0, niter*n), nrow = niter, ncol = n)
kappatt <- rep(0, niter)
gamma0tt <- rep(0, niter)
gamma1tt <- rep(0, niter)
gamma2tt <- rep(0, niter)
zetatt <- rep(0, niter)
sigtt <- rep(0, niter)
rtt <- rep(0, niter)
etatt <- rep(0, niter)
rhott <- rep(0, niter)
ftt <- matrix(0, niter, length(Xm)) # BNP only 
tautt <- rep(0, niter) # BNP only
alptt <- rep(0, niter) # BNP only 
bettt <- matrix(0, niter, 2) 

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
    
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    R_i <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - Zi%*%bett # residual matrix 
    
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
  ############ update f ###############################
  #####################################################
  J <- list() # operator J for f
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    Ji <- matrix(rep(0, mi*length(Xm)), mi, length(Xm))
    for(kk in 1:mi){
      #Ji[kk, Z[[i]][kk,]] <- 1
      Ji[kk,] <- as.numeric(Z[[i]][kk,]==Xm)
    }
    J[[i]] <- Ji
  }
  Dfhalf <- matrix(rep(0, length(Xm)*length(Xm)), nrow =  length(Xm), ncol = length(Xm))
  Ef <- matrix(rep(0, length(Xm)), nrow = length(Xm), ncol= 1)
  for(i in 1:n){
    Dfhalfi <- t(unlist(J[[i]]))%*%solve(Sigt[[i]])%*%unlist(J[[i]])
    Dfhalf <- Dfhalf + Dfhalfi
    
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    Efi <- t(unlist(J[[i]]))%*%solve(Sigt[[i]])%*%(as.matrix(Y[[i]][,1], mi, 1) - etat*logWt - Zi%*%bett)
    Ef <- Ef + Efi
  }
  Df <- Dfhalf + alp*K
  V <- solve(Df)*sigt
  Mu <- solve(Df)%*%Ef
  ft <- matrix(mvrnorm(1, Mu, V), nrow = length(Xm), ncol = 1)  
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ############ update beta ############################
  #####################################################
  DB <- matrix(rep(0, 4), 2, 2)
  EB <- matrix(rep(0, 2), 2, 1)
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    Hybi <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - etat*logWt
    DBi <- as.matrix((1/sigt)*t(Zi)%*%solve(Sigt[[i]])%*%Zi, 
                     nrow = mi, ncol = mi)
    EBi <- as.matrix((1/sigt)*t(Zi)%*%solve(Sigt[[i]])%*%Hybi, 
                     nrow = mi, ncol = 1)
    
    DB <- DB + DBi
    EB <- EB + EBi
  }
  Sb <- solve(DB)
  Mb <- Sb%*%EB
  bett <- matrix(mvrnorm(1, Mb, Sb), 2, 1)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################ update tau #########################
  #####################################################
  At <- 1 # hyper parameter
  Bt <- 1 # hyper parameter 
  lXm <- length(Xm)
  Atstar <- (lXm - 2)/2 + At
  RSSt <- t(ft)%*%K%*%ft
  #for(i in 1:n){
  # RSSti <- t(ft)%*%K%*%ft
  #  RSSt <- RSSti + RSSt
  #}
  Btstar <- RSSt/2 + Bt
  taut <- rgamma(1, Atstar, rate = Btstar)
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
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    Ri <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - etat*logWt - Zi%*%bett
    RSSi <- t(Ri)%*%solve(Sigt[[i]])%*%Ri
    
    RSS <- RSS + RSSi
    
  }
  Bsstar <- 0.5*RSS + Bs
  sigt <- rinvgamma(1, shape = Asstar, scale = Bsstar)
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################ update alpha #######################
  #####################################################
  alp <- taut*sigt
  
  
  
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
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    R_i <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - Zi%*%bett
    
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
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 2)
    
    Rii <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - etat*log(Wt[i]) - Zi%*%bett
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
    ftt[k-nburn,] <- ft # BNP only
    bettt[k-nburn,] <- bett
    tautt[k-nburn] <- taut # BNP only 
    alptt[k-nburn] <- alp # BNP only 
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

W_bnp <- matrix(rep(0, n), nrow = 1, ncol = n)
for(i in 1:n){
  W_bnp[,i] <- mean(Wtt[,i])
}

##################### coefficient of W, eta #################################
(eta_bnp <- mean(etatt)) # posterior mean
(sd(etatt)) # standard deviation 

##################### regression function f #################################
f_bnp <- as.matrix(apply(ftt, 2, mean), length(Xm), 1)

gamma1 <- mean(gamma1tt) # posterior mean 
sd(gamma1tt) # standard deviation
######### history Yij-1 gamma2 ##############################################
gamma2 <- mean(gamma2tt) # posterior mean 
quantile(gamma2tt, probs = c(0.025, 0.975))
sd(gamma2tt) # standard deviation

(rho <- mean(rhott)) 
sd(rhott)
quantile(rhott, probs=c(0.025, 0.975))

beta_bnp <- matrix(apply(bettt, 2, mean), 2, 1)
apply(bettt, 2, sd)
apply(bettt, 2, quantile,  probs = c(0.025, 0.975), na.rm = TRUE)

Z_mean <- cbind(69.07623, 0.5874439)

Wmean <- mean(W_bnp)
yhat <- f_bnp + eta_bnp*log(Wmean) + c(Z_mean%*%beta_bnp)
lines(Xm, yhat, lwd = 2, lty = 1)

cpo <- rep(0, n)
logcpo <- rep(0, n)
for(i in 1:n){
  mi <- nrow(Z[[i]])
  mu_overall <- 0
  Sigm_i <- matrix(rep(0, mi*mi), mi, mi)
  for(kk in 1:1000){
    mu_overall <- ftt[kk,][is.element(Xm, Z[[i]][,1])] + etatt[kk]*log(Wtt[kk,i]) + c(Z_mean%*%bettt[kk,])
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


Rst.tkn.NB.BNP <- list()
Rst.tkn.NB.BNP[[1]] <- W_bnp
Rst.tkn.NB.BNP[[2]] <- eta_bnp
Rst.tkn.NB.BNP[[3]] <- sd(etatt)
Rst.tkn.NB.BNP[[4]] <- gamma2
Rst.tkn.NB.BNP[[5]] <- sd(gamma2tt)
Rst.tkn.NB.BNP[[6]] <- mean(rhott)
Rst.tkn.NB.BNP[[7]] <- sd(rhott)
Rst.tkn.NB.BNP[[8]] <- f_bnp
Rst.tkn.NB.BNP[[9]] <- apply(bettt, 2, mean)
Rst.tkn.NB.BNP[[10]] <- apply(bettt, 2, sd)
Rst.tkn.NB.BNP[[11]] <- lpml
names(Rst.tkn.NB.BNP) <- c("mean of W bnp", "eta hat", "eta sd", 
                           "gamma2 hat", "sd gamma2", "rho hat", 
                           "sd rho", "f hat", 
                           "beta hat", "beta sd", "lpml")

save(Rst.tkn.NB.BNP, file = paste(dir_tkn, "Rst.tkn.NB.BNP.RData", sep = '/'))

