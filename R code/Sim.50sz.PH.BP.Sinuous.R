## ## clean history
rm(list=ls())

## library
library(MASS)
library(mvtnorm)
library(LaplacesDemon)

## path
dir_sinuous <- "~/Desktop/Dissertation/ODFL/Simulation Revision/Data_Sinuous"

## read data 
Yall_50 <- get(load(paste(dir_sinuous, "Sinuous Simulated 200 Response (Y) n=100.RData", sep = '/')))
Uall_50 <- get(load(paste(dir_sinuous, "Sinuous Simulated 200 Response (Uall) n=100.RData", sep = '/')))
Zall_50 <- get(load(paste(dir_sinuous, "Sinuous Simulated 200 Response (Zall) n=100.RData", sep = '/')))
Xmall_50 <- get(load(paste(dir_sinuous, "Sinuous Simulated 200 Response (Xm) n=100.RData", sep = '/')))

## get n and m 
n <- 100
mi <- 5
N <- n*mi
#################

## 100 simulation 
nsim <- 200
theta_sim <- list()
theta_sd_sim <- list()
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
  
  
  
  #############################################################################
  ######################### main algorithm ####################################
  #############################################################################
  ## initial values for MCMC samples 
  Wt <- rep(1, n)
  ## mi = 4
  ## Vt <- as.matrix(rbind(rep(1.33, n), rep(1.33, n), rep(1.33, n)), nrow = mi-1, ncol = n)
  ## mi = 5
  #  Vt <- as.matrix(rbind(rep(1.5, n), rep(1.5, n), rep(1.5, n), rep(1.5, n)), nrow = mi-1, ncol = n)
  rho <- 0.5
  sigt <- 1
  one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
  cov_y <- (1-rho)*diag(1, mi) + rho*one
  Sigt <- cov_y # initial value for Sig
  
  kappat <- 1
  gammat <- 1
  
  etat <- 2 # etat <- 0.2
  thetat <- rep(0, 4) 
  lambdat <- 0.5
  zt <- 1
  
  
  #############################################################################
  #################### set up sample containers ###############################
  #############################################################################
  niter <- 5000 # samples 
  nburn <- 5000 # burning 
  aiter <- niter + nburn # all 
  
  #  Vtt <-rep(matrix(0, niter, mi-1), n)
  Wtt <- matrix(rep(0, niter*n), nrow = niter, ncol = n)
  kappatt <- rep(0, niter)
  gammatt <- rep(0, niter)
  
  sigtt <- rep(0, niter)
  lambdatt <- rep(0, niter)
  etatt <- rep(0, niter)
  rhott <- rep(0, niter)
  thetatt <- matrix(0, niter, 4)
  
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
      
      R_i <- as.matrix(Y[,i], mi, 1) - unlist(Z_all[[i]])%*%thetat
      expLUWc <- 0 
      for(j in 1:(mi-1)){
        expLUWci <- gammat*(as.matrix(Y[,i], mi, 1))[j,] - lambdat*U_all[j,]*Wci[i]*exp(gammat*(as.matrix(Y[,i], mi, 1))[j,])
        expLUWc <- expLUWc + expLUWci
      }
      expLUWt <- 0 
      for(j in 1:(mi-1)){
        expLUWti <- gammat*(as.matrix(Y[,i], mi, 1))[j,] - lambdat*U_all[j,]*Wti[i]*exp(gammat*(as.matrix(Y[,i], mi, 1))[j,])
        expLUWt <- expLUWt + expLUWti
      }
      
      numer <- -log(Wci) - (tlogWc%*%sc%*%logWc - (etat/sigt)*tlogWc%*%solve(Sigt)%*%R_i)
      -((log(Wci) - 1)^2) + expLUWc
      
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
    ############ update theta ###########################
    #####################################################
    DF <- matrix(rep(0, 16), 4, 4)
    EF <- matrix(rep(0, 4), 4, 1)
    for(i in 1:n){
      logWt <- diag(log(Wt[i]), mi)%*%matrix(1, nrow = mi, ncol = 1)
      tlogWt <- t(logWt)
      
      Hyi <- Y[,i] - etat*logWt
      DFi <- as.matrix((1/sigt)*t(unlist(Z_all[[i]]))%*%solve(Sigt)%*%unlist(Z_all[[i]]), 
                       nrow = mi, ncol = mi)
      EFi <- as.matrix((1/sigt)*t(unlist(Z_all[[i]]))%*%solve(Sigt)%*%Hyi, 
                       nrow = mi, ncol = 1)
      
      DF <- DF + DFi
      EF <- EF + EFi
    }
    S <- solve(DF)
    M <- S%*%EF
    thetat <- matrix(mvrnorm(1, M, S), 4, 1)
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
      
      Ri <- as.matrix(Y[,i], mi, 1) - unlist(Z_all[[i]])%*%thetat - etat*logWt
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
    ############### update eat ##########################
    #####################################################
    iDl <- 0
    El <- 0
    for(i in 1:n){
      ones <- matrix(1,nrow=mi,ncol=1)
      onet <- t(ones)
      
      R_i <- as.matrix(Y[,i], mi, 1) - unlist(Z_all[[i]])%*%thetat
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
      Rii <- as.matrix(Y[,i], mi, 1) - unlist(Z_all[[i]])%*%thetat - etat*logWt
      
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
      thetatt[k-nburn,] <- thetat
      etatt[k-nburn] <- etat
      sigtt[k-nburn] <- sigt
      gammatt[k-nburn] <- gammat
      lambdatt[k-nburn] <- lambdat
      rhott[k-nburn] <- rhot
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
      mu_overall <- unlist(Z_all[[i]])%*%thetatt[kk,] + etatt[kk]*log(Wtt[kk,i])
      Sigm_i <- sigtt[kk]*((1-rhott[kk])*diag(1, mi) + rhott[kk]*one)
    }
    densi <- dmvnorm(Y[,i], mean = mu_overall, sigma = Sigm_i)
    
    cpo[i] <- 1/(mean(1/densi))
    logcpo[i] <- log(cpo[i])
  }
  lpml <- mean(logcpo)
  
  
  ## collecte simulation result 
  theta_sim[[s]] <- apply(thetatt, 2, mean)
  theta_sd_sim[[s]] <- apply(thetatt, 2, quantile, probs = c(0.025, 0.975))
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
  Xi <- Xm
  theta.hat <- theta_sim[[s]]
  
  X_design <- rbind(rep(1, length(Xi)), Xi, (Xi)^2, (Xi)^3)
  dx <- diff(Xi)
  nx <- length(Xi)-1
  ytrue <- 5*sin(pi*(Xi-10)/10)/(1+2*((Xi-10)/5)^2*(sign((Xi-10)/5)+1))+4
  yhat <- t(X_design)%*%theta.hat 
  
  imse_sim[s] <- (ytrue[1:nx] - yhat[1:nx])^2%*%dx
  
  
}

Rst.100sz.PH.BP.Sinuous.200sim <- list()
Rst.100sz.PH.BP.Sinuous.200sim[[1]] <- eta_sim
Rst.100sz.PH.BP.Sinuous.200sim[[2]] <- eta_sd
Rst.100sz.PH.BP.Sinuous.200sim[[3]] <- gamma_sim
Rst.100sz.PH.BP.Sinuous.200sim[[4]] <- gamma_sd_sim
Rst.100sz.PH.BP.Sinuous.200sim[[5]] <- lambda_sim
Rst.100sz.PH.BP.Sinuous.200sim[[6]] <- lambda_sd_sim
Rst.100sz.PH.BP.Sinuous.200sim[[7]] <- rho_sim
Rst.100sz.PH.BP.Sinuous.200sim[[8]] <- rho_sd_sim
Rst.100sz.PH.BP.Sinuous.200sim[[9]] <- lpml_sim
Rst.100sz.PH.BP.Sinuous.200sim[[10]] <- theta_sim
Rst.100sz.PH.BP.Sinuous.200sim[[11]] <- theta_sd_sim
Rst.100sz.PH.BP.Sinuous.200sim[[12]] <- W_sim
Rst.100sz.PH.BP.Sinuous.200sim[[13]] <- eta_mse
Rst.100sz.PH.BP.Sinuous.200sim[[14]] <- imse_sim
names(Rst.100sz.PH.BP.Sinuous.200sim) <- c("eta", "eta_sd", "gamma", "gamma_sd", 
                                          "lambda", "lambda_sd", 
                                          "rho", "rho_sd", 
                                          "lpml_sim", "theta_sim", "theta_sd", 
                                          "W_sim", "eta mse", "imse sim")
dir_rst <- "~/Desktop/Dissertation/ODFL/Simulation Revision/Result Sinuous"
save(Rst.100sz.PH.BP.Sinuous.200sim, file = paste(dir_rst, "Rst.100sz.PH.BP.Sinuous.200sim.RData", sep = '/'))








