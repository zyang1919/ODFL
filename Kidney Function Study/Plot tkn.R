## clean history
rm(list=ls())

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
Xm <- min(Zm):max(Zm)
Xv <- as.vector(Zm)
Xm1 <- sort(unique(Xv))
plot(c(1), c(1), type = 'n', ylim = c(3, 110), xlim = c(min(Xm), max(Xm)), xlab = " ",
     ylab = "")
cl <- colors()
for (i in 1:n){
  lines(unlist(Z[[i]])[,1], unlist(Y[[i]]), type = 'l', col = rgb(0.3, 0.3, 0.3), lwd = 1)
}

Z_mean <- cbind(69.07623, 0.5874439)
cols <- rainbow(3, s=0.5)

## get result 
nb.bp <- get(load(paste(dir_tkn, "Rst.tkn.NB.BP.RData", sep = '/')))
Wmean <- mean(nb.bp$`mean of W bp`)
theta <- nb.bp$`theta hat`
eta_bp <- nb.bp$`eta hat`
bet_bp <- nb.bp$`beta hat`
yhat <- theta[2]*Xm + theta[3]*Xm^2 + theta[1] + eta_bp*log(Wmean) + c(Z_mean%*%bet_bp)
lines(Xm, yhat, lwd = 5, lty = 1, col = "blue")

nb.bnp <- get(load(paste(dir_tkn, "Rst.tkn.NB.BNP.RData", sep = '/')))
Wmean <- mean(nb.bnp$`mean of W bnp`)
f_bnp <- nb.bnp$`f hat`
eta_bnp <- nb.bnp$`eta hat`
beta_bnp <- nb.bnp$`beta hat`
yhat <- f_bnp + eta_bnp*log(Wmean) + c(Z_mean%*%beta_bnp)
lines(Xm1, yhat, lwd = 5, lty = 1, col = "red")
