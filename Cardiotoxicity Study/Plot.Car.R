## clean history
rm(list=ls())

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
n <- length(dnew) 
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

## plot data
plot(c(1), c(1), type = 'n', ylim = c(-6, 6), xlim = c(min(Xm), max(Xm)), xlab = "",
     ylab = "")
for (i in 1:n){
  lines(unlist(Z[[i]])[,1], unlist(Y[[i]]), type = 'l', col = rgb(0.1, 0.1, 0.1), lwd = 0.5)
}
Z_cov <- cbind(0.5177994, 0.4919094, 6.467024, 0.2977346)
cols <- rainbow(3, s=0.5)


## get result 
dir_car <- "~/Desktop/Dissertation/ODFL/Data Application Revision/Cardiotoxicity "

nb.bp <- get(load(paste(dir_car, "Rst.Car.NB.BP.RData", sep = '/')))
Wmean <- mean(nb.bp$`mean of W bp`)
theta <- nb.bp$`theta hat`
eta_bp <- nb.bp$`eta hat`
bet_bp <- nb.bp$`beta hat`
yhat <- theta[2]*Xm + theta[1] + eta_bp*log(Wmean) + c(Z_cov%*%bet_bp)
lines(Xm, yhat, lwd = 5, lty = 1, col = "blue")

nb.bnp <- get(load(paste(dir_car, "Rst.Car.NB.BNP.RData", sep = '/')))
Wmean <- mean(nb.bnp$`mean of W bnp`)
f_bnp <- nb.bnp$`f hat`
eta_bnp <- nb.bnp$`eta hat`
beta_bnp <- nb.bnp$`beta hat`
yhat <- f_bnp + eta_bnp*log(Wmean) + c(Z_cov%*%beta_bnp)
lines(Xm, yhat, lwd = 5, lty = 1, col = "red")

