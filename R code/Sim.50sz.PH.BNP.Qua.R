## ## clean history
rm(list=ls())

## library
library(MASS)
library(mvtnorm)
library(LaplacesDemon)

## path
dir_odfl <- "~/Desktop/Dissertation/ODFL/Simulation Revision"

## read data 
Yall_50 <- get(load(paste(dir_odfl, "Simulated 100 Response (Y) n=50.RData", sep = '/')))
Uall_50 <- get(load(paste(dir_odfl, "Simulated 100 Response (Uall) n=50.RData", sep = '/')))
Zall_50 <- get(load(paste(dir_odfl, "Simulated 100 Response (Zall) n=50.RData", sep = '/')))
Xmall_50 <- get(load(paste(dir_odfl, "Simulated 100 Response (Xm) n=50.RData", sep = '/')))

## get n and m 
n <- 50
mi <- 5
N <- n*mi
#################

getK <- function(x){
  # This function returns K matrix, and requires unique design points x
  
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


## 200 simulation 
nsim <- 200
f_sim <- list()
f_sd_sim <- list()
eta_sim <- rep(0, nsim)
eta_sd <- rep(0, nsim)
W_sim <- list()
gamma_sim <- rep(0, nsim)
gamma_sd_sim <- rep(0, nsim)
lambda_sim <- rep(0, nsim)
lambda_sd_sim <- rep(0, nsim)
rho_sim <- rep(0, nsim)
rho_sd_sim <- rep(0, nsim)
lpml_sim <- rep(0, nsim)
imse_sim <- rep(0, nsim)
eta_mse <- rep(0, nsim)


for(s in 1:nsim){
  ## get data for each simulation 
  Y <- Yall_50[[s]]
  U_all <- Uall_50[[s]]
  Z_all <- Zall_50[[s]]
  Xm <- Xmall_50[[s]]
  K <- getK(Xm)
  X <- matrix(rep(0, mi*n), nrow = mi, ncol = n)
  for(i in 1:n){
    X[,i] <- as.matrix(unlist(Z_all[[i]][, 2]), mi, 1)
  }
  
  
  
  #############################################################################
  ######################### main algorithm ####################################
  #############################################################################
  ## initial values for MCMC samples 
  Wt <- rep(1, n)
  ## mi = 4
  ## Vt <- as.matrix(rbind(rep(1.33, n), rep(1.33, n), rep(1.33, n)), nrow = mi-1, ncol = n)
  ## mi = 5
  #Vt <- as.matrix(rbind(rep(1.5, n), rep(1.5, n), rep(1.5, n), rep(1.5, n)), nrow = mi-1, ncol = n)
  rho <- 0.5
  sigt <- 1
  one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
  cov_y <- (1-rho)*diag(1, mi) + rho*one
  Sigt <- cov_y # initial value for Sig

  kappat <- 1
  
  lambdat <- 0.5
  gammat <- 1
  
  alp  <- 0.01 # BNP only
  taut <- alp/sigt # BNP only

  etat <- 2 # etat <- 0.2
  ft <- rep(0, length(Xm)) # BNP only
  zt <- 1
  
  
  #############################################################################
  #################### set up sample containers ###############################
  #############################################################################
  niter <- 5000 # samples 
  nburn <- 5000 # burning 
  aiter <- niter + nburn # all 
  
  #Vtt <-rep(matrix(0, niter, mi-1), n)
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
  
  for (k in 1:aiter) {
    
    #####################################################
    ################## update Wi ########################
    #####################################################
    sc <- (etat^2/sigt)*solve(Sigt)
    #Wc <- rlnorm(n, meanlog = Wt, sdlog = 1)
    ## proposal from log normal distribution 
    for (i in 1:n) {
      Wci <- rlnorm(1, meanlog = Wt[i] , sdlog = 1)
      Wti <- Wt[i]
      
      #logWc <- diag(log(Wc[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
      logWc <- diag(log(Wci), mi)%*%matrix(1, nrow = mi, ncol = 1)
      tlogWc <- t(logWc)
      
      logWt <- diag(log(Wti), mi)%*%matrix(1, nrow = mi, ncol = 1)
      tlogWt <- t(logWt)
      
      R_i <- as.matrix(Y[,i], mi, 1) - ft[is.element(Xm, X[,i])]
      expLUWc <- 0 
      for(j in 1:(mi-1)){
        expLUWci <- gammat*(as.matrix(Y[,i], mi, 1))[j,] - lambdat*U_all[j,i]*Wci*exp(gammat*(as.matrix(Y[,i], mi, 1))[j,])
        expLUWc <- expLUWc + expLUWci
      }
      
      expLUWt <- 0 
      for(j in 1:(mi-1)){
        expLUWti <- gammat*(as.matrix(Y[,i], mi, 1))[j,] - lambdat*U_all[j,i]*Wti*exp(gammat*(as.matrix(Y[,i], mi, 1))[j,])
        expLUWt <- expLUWt + expLUWti
      }
      
      numer <- -log(Wci) - (tlogWc%*%sc%*%logWc - (etat/sigt)*tlogWc%*%solve(Sigt)%*%R_i)
      -((log(Wci) - 1)^2)  + expLUWc
      
      denom <- -log(Wti) - (tlogWt%*%sc%*%logWt - (etat/sigt)*tlogWt%*%solve(Sigt)%*%R_i)
      -((log(Wti) - 1)^2) + expLUWt
      
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
    Bkstar <- Bk + sum((log(Wt[i]) - 1)^2)
    #ikappat <- rgamma(1, Akstar, rate = Bkstar)
    #kappat <- 1/ikappat
    kappat <- rinvgamma(1, shape=Akstar, scale = Bkstar)
    #####################################################
    #####################################################
    #####################################################
    
    #####################################################
    ############ update f ###############################
    #####################################################
    J <- list() # operator J for f
    for(i in 1:n){
      Ji <- matrix(rep(0, mi*length(Xm)), mi, length(Xm))
      for(kk in 1:mi){
        Ji[kk,] <- as.numeric(X[kk,i]==Xm)
      }
      J[[i]] <- Ji
    }
    Dfhalf <- matrix(rep(0, length(Xm)*length(Xm)), nrow =  length(Xm), ncol = length(Xm))
    Ef <- matrix(rep(0, length(Xm)), nrow = length(Xm), ncol= 1)
    for(i in 1:n){
      Dfhalfi <- t(unlist(J[[i]]))%*%solve(Sigt)%*%unlist(J[[i]])
      Dfhalf <- Dfhalf + Dfhalfi
      
      logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
      Efi <- t(unlist(J[[i]]))%*%solve(Sigt)%*%(Y[,i] - etat*logWt)
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
    ################ update tau #########################
    #####################################################
    At <- 1 # hyper parameter
    Bt <- 1 # hyper parameter 
    Atstar <- (length(Xm) - 2)/2 + At
    RSSt <- t(ft)%*%K%*%ft
    #for(i in 1:n){
    #RSSti <- t(ft)%*%K%*%ft
    #RSSt <- RSSti + RSSt
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
      logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
      tlogWt <- t(logWt)
      
      Ri <- as.matrix(Y[,i], mi, 1) - ft[is.element(Xm, X[,i])] - etat*logWt
      RSSi <- t(Ri)%*%solve(Sigt)%*%Ri
      
      RSS <- RSS + RSSi
      
    }
    Bsstar <- 0.5*RSS + Bs
    #isigt <- rgamma(1, Asstar, rate = Bsstar)
    #sigt <- 1/isigt
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
      ones <- matrix(1,nrow=mi,ncol=1)
      onet <- t(ones)
      
      R_i <- as.matrix(Y[,i], mi, 1) - ft[is.element(Xm, X[,i])]
      iDli <- Wt[i]^2*onet%*%solve(Sigt)%*%ones
      Eli <- Wt[i]*onet%*%R_i
      
      iDl <- iDl + iDli 
      El <- El + Eli
    }
    Dl <- 1/iDl
    etat <- rnorm(1, Dl*El, sqrt(Dl*sigt))
    
    #  Dl <- 0
    #  El <- 0
    #  for(i in 1:n){
    #    logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
    #    tlogWt <- t(logWt)
    
    #    R_i <- as.matrix(Y[,i], mi, 1) - unlist(Z_all[[i]])%*%thetat
    #    Dli <- as.matrix((1/sigt)*tlogWt%*%solve(Sigt)%*%logWt, 1, 1)
    #    Eli <- as.matrix((1/sigt)*tlogWt%*%solve(Sigt)%*%R_i, 1, 1)
    
    #    Dl <- Dl + Dli
    #    El <- El + Eli
    #  }
    #  S_e <- 1/Dl
    #  M_e <- S_e*El
    #  etat <- rnorm(1, M_e, sqrt(S_e))
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
      for(j in 1:(mi-1)){
        loggci <- sum(gammac*Y[j,i] - lambdat*U_all[j,i]*Wt[i]*exp(gammac*Y[j,i]))
        loggti <- sum(gammat*Y[j,i] - lambdat*U_all[j,i]*Wt[i]*exp(gammat*Y[j,i]))
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
      for(j in 1:(mi-1)){
        UWei <- U_all[j,i]*Wt[i]*exp(gammat*Y[j,i])
      }
      UWe <- UWe + UWei
    }
    Blstar <- UWe + Bl
    lambdat <- rgamma(1, Alstar, rate = Blstar)
    
    #####################################################
    #####################################################
    #####################################################
    
    
    
    #####################################################
    ############### update rho (Sig) ####################
    ####################################################
    
    zc <- rnorm(1, zt, 0.1)
    rhoc <- (exp(2*zc)-1)/(exp(2*zc)+1)
    rhot <- (exp(2*zt)-1)/(exp(2*zt)+1)
    one <- matrix(1, mi, mi, byrow = TRUE)
    Sigc <- ((1-rhoc)*diag(1, mi) + rhoc*one)
    Sigt <- ((1-rhot)*diag(1, mi) + rhot*one)
    # transition
    logzc <- zc - 2*log(exp(zc)+1)
    logzt <- zt - 2*log(exp(zt)+1)
    numerr <- 0
    denomr <- 0
    for(i in 1:n){
      logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
      Rii <- as.matrix(Y[,i], mi, 1) - ft[is.element(Xm, X[,i])] - etat*logWt
      
      numerir <- (-0.5)*log(det(Sigc))-(1/2*sigt)*t(Rii)%*%solve(Sigc)%*%Rii
      denomir <- (-0.5)*log(det(Sigt))-(1/2*sigt)*t(Rii)%*%solve(Sigt)%*%Rii
      
      numerr <- numerr+numerir
      denomr <- denomr+denomir
    }
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
      tautt[k-nburn] <- taut # BNP only 
      alptt[k-nburn] <- alp # BNP only 
      etatt[k-nburn] <- etat
      sigtt[k-nburn] <- sigt
      gammatt[k-nburn] <- gammat
      lambdatt[k-nburn] <- lambdat

      rhott[k-nburn] <- rhot
      #Vtt[[k-nburn]] <- Vt
    }
  }
  ## get mean W
  W <- matrix(rep(0, n), nrow = 1, ncol = n)
  for(i in 1:n){
    W[,i] <- mean(Wtt[,i])
  }
  
  ### calculate lpml ####
  cpo <- rep(0, n)
  logcpo <- rep(0, n)
  for(i in 1:n){
    mu_overall <- 0
    Sigm_i <- matrix(rep(0, mi*mi), mi, mi)
    for(kk in 1:5000){
      mu_overall <- ftt[kk,][is.element(Xm, X[,i])] + etatt[kk]*log(Wtt[kk,i])
      Sigm_i <- sigtt[kk]*((1-rhott[kk])*diag(1, mi) + rhott[kk]*one)
    }
    densi <- dmvnorm(Y[,i], mean = mu_overall, sigma = Sigm_i)
    
    cpo[i] <- 1/(mean(1/densi))
    logcpo[i] <- log(cpo[i])
  }
  lpml <- mean(logcpo)
  
  
  f <- as.matrix(apply(ftt, 2, mean), length(Xm), 1)
  ## collecte simulation result 
  f_sim[[s]] <- as.matrix(apply(ftt, 2, mean), length(Xm), 1)
  f_sd_sim[[s]] <- apply(ftt, 2, sd)
  eta_sim[s] <- mean(etatt)
  eta_sd[s] <- sd(etatt)
  W_sim[[s]] <- W
  gamma_sim[s] <- mean(gammatt)
  gamma_sd_sim[s] <- sd(gammatt)
  lambda_sim[s] <- mean(lambdatt)
  lambda_sd_sim[s] <- sd(lambdatt)
  rho_sim[s] <- mean(rhott)
  rho_sd_sim[s] <- sd(rhott)
  lpml_sim[s] <- lpml
  eta_mse[s] <- (mean(etatt) - 1)^2
  
  
  ## calculate IMSE 
  X_design <- Xm
  dx <- diff(X_design)
  nx <- length(X_design)-1
  ft_i <-5 + (-0.01)*X_design + (-0.05)*X_design^2 
  fhat <- f_sim[[s]]
  
  imse_sim[s] <- (ft_i[1:nx] - fhat[1:nx])^2%*%dx
  
}
Rst.50sz.PH.BNP.200sim <- list()
Rst.50sz.PH.BNP.200sim[[1]] <- eta_sim
Rst.50sz.PH.BNP.200sim[[2]] <- eta_sd
Rst.50sz.PH.BNP.200sim[[3]] <- gamma_sim
Rst.50sz.PH.BNP.200sim[[4]] <- gamma_sd_sim
Rst.50sz.PH.BNP.200sim[[5]] <- lambda_sim
Rst.50sz.PH.BNP.200sim[[6]] <- lambda_sd_sim
Rst.50sz.PH.BNP.200sim[[7]] <- rho_sim
Rst.50sz.PH.BNP.200sim[[8]] <- rho_sd_sim
Rst.50sz.PH.BNP.200sim[[9]] <- lpml_sim
Rst.50sz.PH.BNP.200sim[[10]] <- f_sim
Rst.50sz.PH.BNP.200sim[[11]] <- f_sd_sim
Rst.50sz.PH.BNP.200sim[[12]] <- W_sim
Rst.50sz.PH.BNP.200sim[[13]] <- eta_mse
Rst.50sz.PH.BNP.200sim[[14]] <- imse_sim
names(Rst.50sz.PH.BNP.200sim) <- c("eta", "eta_sd", "gamma", "gamma_sd", 
                                  "lambda", "lambda_sd", 
                                  "rho", "rho_sd", 
                                  "lpml_sim", "f_sim", "f_sd", 
                                  "W_sim", "eta mse", "imse sim")
save(Rst.50sz.PH.BNP.200sim, file = paste(dir_odfl, "Rst.50sz.PH.Qua.BNP.200sim.RData", sep = '/'))




