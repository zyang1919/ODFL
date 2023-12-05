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

#################


## 100 simulation 
nsim <- 200
theta_sim <- list()
theta_sd_sim <- list()
eta_sim <- rep(0, nsim)
eta_sd <- rep(0, nsim)
W_sim <- list()
gamma1_sim <- rep(0, nsim)
gamma1_sd_sim <- rep(0, nsim)
gamma2_sim <- rep(0, nsim)
gamma2_sd_sim <- rep(0, nsim)
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
  Vt <- as.matrix(rbind(rep(1.5, n), rep(1.5, n), rep(1.5, n), rep(1.5, n)), nrow = mi-1, ncol = n)
  rho <- 0.5
  sigt <- 1
  one <- matrix(1, nrow = mi, ncol = mi, byrow = FALSE)
  cov_y <- (1-rho)*diag(1, mi) + rho*one
  Sigt <- cov_y # initial value for Sig
  sigg0 <- 1; sigg1 <- 1; sigg2 <- 1
  kappat <- 1
  zetat <- 1 
  gamma0t <- 1; gamma1t <- 1; gamma2t <- 1
  rt <- 2 # rt <- 1
  etat <- 2 # etat <- 0.2
  thetat <- rep(0, 4) 
  zt <- 1
  
  
  #############################################################################
  #################### set up sample containers ###############################
  #############################################################################
  niter <- 5000 # samples 
  nburn <- 5000 # burning 
  aiter <- niter + nburn # all 
  
  Vtt <-rep(matrix(0, niter, mi-1), n)
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
      RSSj <- rep(0, mi-1)
      for(j in 1:(mi-1)){
        RSSj[j] <- Vt[j,i] - gamma0t - gamma2t*(as.matrix(Y[,i], mi, 1))[j,]
      }
      RSS <- sum(RSSj)
      
      numer <- -log(Wci) - (tlogWc%*%sc%*%logWc - (etat/sigt)*tlogWc%*%solve(Sigt)%*%R_i)
      -((log(Wci) - 1)^2) - (zetat)*((mi-1)*(gamma1t*Wci)^2 - (gamma1t*Wci)*RSS)
      
      denom <- -log(Wti) - (tlogWt%*%sc%*%logWt - (etat/sigt)*tlogWt%*%solve(Sigt)%*%R_i)
      -((log(Wti) - 1)^2) - (zetat)*((mi-1)*(gamma1t*Wti)^2 - (gamma1t*Wti)*RSS)
      
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
    
    Vc <- matrix(rep(0, 4*n), mi-1, n)
    for(j in 1:n){
      Vc[,j] <- mvrnorm(1, Vt[,j], Sigma = diag(1, mi-1))
    }
    
    #numeri <- rep(0, n)
    #denomi <- rep(0, n)
    #numerj <- rep(0, (mi-1))
    #denomj <- rep(0, (mi-1))
    
    #for(i in 1:n){
    #for(j in 1:(mi-1)){
    #  numerj[j] <- Vc[j,i]*U_all[j,i] - (0.5/zetat)*(Vc[j,i] - gamma0t - gamma1t*Wt[i] - gamma2t*Y[j,i])^2 - 
    #    rt*log(rt + exp(Vc[j,i])) - U_all[j,i]*log(rt + exp(Vc[j,i])) 
    #  denomj[j] <- Vt[j,i]*U_all[j,i] - (0.5/zetat)*(Vt[j,i] - gamma0t - gamma1t*Wt[i] - gamma2t*Y[j,i])^2 - 
    #   rt*log(rt + exp(Vt[j,i])) - U_all[j,i]*log(rt + exp(Vt[j,i]))
    #  }
    #numeri <- sum(numerj)
    #denomi <- sum(denomj)
    
    #logalpha <- numeri - denomi
    #logu <- log(runif(1, 0, 1))
    
    #if(logu<logalpha){
    #  Vt[,i] <- Vc[,i]
    #}
    #}
    #  numer <- sum(numeri)
    #  denom <- sum(denomi)
    #  logalpha <- numer - denom
    #  logu <- log(runif(1, 0, 1))
    #  if(logu<logalpha){
    #    Vt <- Vc
    #  }
    
    numerV <- 0
    denomV <- 0
    for(j in 1:(mi-1)){
      numerVj <- Vc[j,]*U_all[j,] - (0.5/zetat)*(Vc[j,] - gamma0t - gamma1t*Wt - gamma2t*Y[j,])^2 - 
        rt*log(rt + exp(Vc[j,])) - U_all[j,]*log(rt + exp(Vc[j,]))
      denomVj <- Vt[j,]*U_all[j,] - (0.5/zetat)*(Vt[j,] - gamma0t - gamma1t*Wt - gamma2t*Y[j,])^2 - 
        rt*log(rt + exp(Vt[j,])) - U_all[j,]*log(rt + exp(Vt[j,]))
      
      numerV <- numerV + numerVj
      denomV <- denomV + denomVj
    }
    logalphaV <- numerV - denomV
    loguV <- log(runif(n, 0, 1))
    
    Vt[,loguV<logalphaV] <- Vc[,loguV<logalphaV]
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
    ################ update zeta2 #######################
    #####################################################
    Az <- 100 # hyper parameter
    Bz <- 100 # hyper parameter
    RSSv <- 0
    for(i in 1:n){
      for(j in 1:(mi-1)){
        RSSvi <- (Vt[j,i] - gamma0t - gamma1t*Wt[i] - gamma2t*Y[j,i])^2
      }
      RSSv <- RSSv + RSSvi
    }
    Azstar <- (n*(mi-1))/2 + Az
    Bzstar <- RSSv/2 + Bz
    #izetat <- rgamma(1, Azstar, rate = Bzstar)
    #zetat <- 1/izetat
    zetat <- rinvgamma(1, shape=Azstar, scale=Bzstar)
    #####################################################
    #####################################################
    #####################################################
    
    #####################################################
    ################# update gamma0 #####################
    #####################################################
    RSS0 <- 0
    for(i in 1:n){
      for(j in 1:(mi-1)){
        RSS0i <- Vt[j,i] - gamma1t*Wt[i] - gamma2t*Y[j,i]
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
      for(j in 1:(mi-1)){
        RSSg1i <- Vt[j,i] - gamma0t - gamma2t*Y[j,i]
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
      for(j in 1:(mi-1)){
        RSS2i <- Y[j,i]*(Vt[j,i] - gamma0t - gamma1t*Wt[i])
        iD2i <- ((Y[j,i])^2)/zetat + 1/sigg2
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
    ####################################################
    Ar <- 50
    Br <- 50
    rc <- rgamma(1, rt^2, rate=rt)
    
    
    numerRj <- rep(0, (mi-1))
    denomRj <- rep(0, (mi-1))
    numerRi <- rep(0, n)
    denomRi <- rep(0, n)
    
    for(i in 1:n){
      for(j in 1:(mi-1)){
        numerRj[j] <- log(gamma(U_all[j,i] + rc)) + (Ar + rc-1)*log(rc) + (-rc*Br)-
          (log(gamma(rc)) + rc*log(rc+exp(Vt[j,i])) + U_all[j,i]*log(rc+exp(Vt[j,i])))
        denomRj[j] <- log(gamma(U_all[j,i] + rt)) + (Ar + rt-1)*log(rt) + (-rt*Br)-
          (log(gamma(rt)) + rt*log(rt+exp(Vt[j,i])) + U_all[j,i]*log(rt+exp(Vt[j,i])))
        
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
      zetatt[k-nburn] <- zetat
      gamma0tt[k-nburn] <- gamma0t
      gamma1tt[k-nburn] <- gamma1t
      gamma2tt[k-nburn] <- gamma2t
      rtt[k-nburn] <- rt
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
  gamma1_sim[s] <- mean(gamma1tt)
  gamma1_sd_sim[s] <- sd(gamma1tt)
  gamma2_sim[s] <- mean(gamma2tt)
  gamma2_sd_sim[s] <- sd(gamma2tt)
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
Rst.100sz.NB.BP.Sinuous.200sim <- list()
Rst.100sz.NB.BP.Sinuous.200sim[[1]] <- eta_sim
Rst.100sz.NB.BP.Sinuous.200sim[[2]] <- eta_sd
Rst.100sz.NB.BP.Sinuous.200sim[[3]] <- gamma1_sim
Rst.100sz.NB.BP.Sinuous.200sim[[4]] <- gamma1_sd_sim
Rst.100sz.NB.BP.Sinuous.200sim[[5]] <- gamma2_sim
Rst.100sz.NB.BP.Sinuous.200sim[[6]] <- gamma2_sd_sim
Rst.100sz.NB.BP.Sinuous.200sim[[7]] <- rho_sim
Rst.100sz.NB.BP.Sinuous.200sim[[8]] <- rho_sd_sim
Rst.100sz.NB.BP.Sinuous.200sim[[9]] <- lpml_sim
Rst.100sz.NB.BP.Sinuous.200sim[[10]] <- theta_sim
Rst.100sz.NB.BP.Sinuous.200sim[[11]] <- theta_sd_sim
Rst.100sz.NB.BP.Sinuous.200sim[[12]] <- W_sim
Rst.100sz.NB.BP.Sinuous.200sim[[13]] <- eta_mse
Rst.100sz.NB.BP.Sinuous.200sim[[14]] <- imse_sim
names(Rst.100sz.NB.BP.Sinuous.200sim) <- c("eta", "eta_sd", "gamma1", "gamma1_sd", 
                                          "gamma2", "gamma2_sd", 
                                          "rho", "rho_sd", 
                                          "lpml_sim", "theta_sim", "theta_sd", 
                                          "W_sim", "eta mse", "imse sim")
dir_rst <- "~/Desktop/Dissertation/ODFL/Simulation Revision/Result Sinuous"
save(Rst.100sz.NB.BP.Sinuous.200sim, file = paste(dir_rst, "Rst.100sz.NB.BP.Sinuous.200sim.RData", sep = '/'))


