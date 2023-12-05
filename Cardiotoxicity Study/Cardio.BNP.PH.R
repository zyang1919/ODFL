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

gammat <- 1
lambdat <- 0.5

etat <- 2
zt <- 4

ft <- rep(0, length(Xm)) # BNP only
alp  <- 0.01 # BNP only
taut <- alp/sigt # BNP only
bett <- matrix(rep(1, 4), nrow = 4, nco = 1) # regression coefficient for covariates (sex, dose and age)

niter <- 1000 # samples 
nburn <- 1000 # burning 
aiter <- niter + nburn # all
Wtt <- matrix(rep(0, niter*n), nrow = niter, ncol = n)
kappatt <- rep(0, niter)
gammatt <- rep(0, niter)

sigtt <- rep(0, niter)
lambdatt <- rep(0, niter)

etatt <- rep(0, niter)
rhott <- rep(0, niter)
ftt <- matrix(0, niter, length(Xm)) # BNP only 
tautt <- rep(0, niter) # BNP only
alptt <- rep(0, niter) # BNP only 
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
    R_i <- as.matrix(Y[[i]][,1], mi, 1) - ft[is.element(Xm, Z[[i]][,1])] - Zi%*%bett # residual matrix 
    
    expLUWc <- 0 
    for(j in 1:(mi-1)){
      expLUWci <- gammat*(as.matrix(Y[[i]][,1], mi, 1)) - lambdat*(unlist(U_all[[i]][j,]))*Wci*exp(gammat*(as.matrix(Y[[i]][,1], mi, 1)))
      expLUWc <- expLUWc + expLUWci
    }
    expLUWt <- 0 
    for(j in 1:(mi-1)){
      expLUWti <- gammat*(as.matrix(Y[[i]][,1], mi, 1)) - lambdat*(unlist(U_all[[i]][j,]))*Wti*exp(gammat*(as.matrix(Y[[i]][,1], mi, 1)))
      expLUWt <- expLUWt + expLUWti
    }
    
    numer <- -log(Wci) - (tlogWc%*%sc%*%logWc - (etat/sigt)*tlogWc%*%solve(Sigt[[i]])%*%R_i)
    -((log(Wci) - 1)^2) + expLUWc
    
    denom <- -log(Wti) - (tlogWt%*%sc%*%logWt - (etat/sigt)*tlogWt%*%solve(Sigt[[i]])%*%R_i)
    -((log(Wti) - 1)^2) + expLUWt
    
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
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
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
  DB <- matrix(rep(0, 16), 4, 4)
  EB <- matrix(rep(0, 4), 4, 1)
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    tlogWt <- t(logWt)
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
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
  bett <- matrix(mvrnorm(1, Mb, Sb), 4, 1)
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
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
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
    
    Zi <- matrix(sapply(Z_else[[i]], as.numeric), mi, ncol = 4)
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
  ################# update gamma #####################
  #####################################################
  gammac <- rnorm(1, mean = gammat, sd = 1)
  
  loggc <- 0
  loggt <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      loggci <- sum(gammac*(unlist(Y[[i]][j,])) - lambdat*(unlist(U_all[[i]][j,]))*Wt[i]*exp(gammac*unlist(Y[[i]][j,])))
      loggti <- sum(gammat*(unlist(Y[[i]][j,])) - lambdat*(unlist(U_all[[i]][j,]))*Wt[i]*exp(gammat*unlist(Y[[i]][j,])))
    }
    loggc <- loggci + loggc
    loggt <- loggti + loggt
  }
  
  numer <- loggc
  denom <- loggt
  
  logalpha <- numer - denom 
  logu <- log(runif(1, 0, 1))
  if(logu<logalpha){
    gammat <- gammac
  }
  #####################################################
  #####################################################
  #####################################################
  
  #####################################################
  ################ update lambda ######################
  #####################################################
  Al <- 1 # hyper parameter
  Bl <- 1 # hyper parameter
  
  Alstar <- n + Al
  UWe <- 0
  for(i in 1:n){
    mi <- nrow(Z[[i]])
    for(j in 1:(mi-1)){
      UWei <- (unlist(U_all[[i]][j,]))*Wt[i]*exp(gammat*unlist(Y[[i]][j,]))
    }
    UWe <- UWe + UWei
  }
  Blstar <- UWe + Bl
  lambdat <- rgamma(1, Alstar, rate = Blstar)
  
  #####################################################
  #####################################################
  
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
    gammatt[k-nburn] <- gammat
    lambdatt[k-nburn] <- lambdat
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

(gamma <- mean(gammatt)) # posterior mean 
quantile(gammatt, probs = c(0.025, 0.975))
sd(gammatt) # standard deviation

(rho <- mean(rhott)) 
sd(rhott)
quantile(rhott, probs=c(0.025, 0.975))

beta_bnp <- matrix(apply(bettt, 2, mean), 4, 1)
apply(bettt, 2, sd)
apply(bettt, 2, quantile,  probs = c(0.025, 0.975), na.rm = TRUE)

Z_cov <- cbind(0.5177994, 0.4919094, 6.467024, 0.2977346)


cpo <- rep(0, n)
logcpo <- rep(0, n)
for(i in 1:n){
  mi <- nrow(Z[[i]])
  mu_overall <- 0
  Sigm_i <- matrix(rep(0, mi*mi), mi, mi)
  for(kk in 1:1000){
    mu_overall <- ftt[kk,][is.element(Xm, Z[[i]][,1])] + etatt[kk]*log(Wtt[kk,i])
    one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
    Sigm_i <- sigtt[kk]*((1-rhott[kk])*diag(1, mi) + rhott[kk]*one)
  }
  
  Yii <- matrix(unlist(Y[[i]]), ncol = mi)
  densi <- dmvnorm(Yii, mean = mu_overall, sigma = Sigm_i)
  
  cpo[i] <- 1/(mean(1/densi))
  logcpo[i] <- log(cpo[i])
}
lpml <- mean(logcpo)


Rst.Car.PH.BNP <- list()
Rst.Car.PH.BNP[[1]] <- W_bnp
Rst.Car.PH.BNP[[2]] <- eta_bnp
Rst.Car.PH.BNP[[3]] <- sd(etatt)
Rst.Car.PH.BNP[[4]] <- gamma
Rst.Car.PH.BNP[[5]] <- sd(gammatt)
Rst.Car.PH.BNP[[6]] <- mean(rhott)
Rst.Car.PH.BNP[[7]] <- sd(rhott)
Rst.Car.PH.BNP[[8]] <- f_bnp
Rst.Car.PH.BNP[[9]] <- apply(bettt, 2, mean)
Rst.Car.PH.BNP[[10]] <- apply(bettt, 2, sd)
Rst.Car.PH.BNP[[11]] <- lpml
names(Rst.Car.PH.BNP) <- c("mean of W bnp", "eta hat", "eta sd", 
                           "gamma hat", "sd gamma", "rho hat", 
                           "sd rho", "f hat", 
                           "beta hat", "beta sd", "lpml")

save(Rst.Car.PH.BNP, file = paste(dir_car, "Rst.Car.PH.BNP.RData", sep = '/'))



