## ## clean history
rm(list=ls())

## read data 
dir_odfl <- "~/Desktop/Dissertation/ODFL/Simulation Revision/Data"

## data set
Yall_50 <- get(load(paste(dir_odfl, "Simulated 200 Response (Y) n=50.RData", sep = '/')))
Uall_50 <- get(load(paste(dir_odfl, "Simulated 200 Response (Uall) n=50.RData", sep = '/')))
Zall_50 <- get(load(paste(dir_odfl, "Simulated 200 Response (Zall) n=50.RData", sep = '/')))
Xmall_50 <- get(load(paste(dir_odfl, "Simulated 200 Response (Xm) n=50.RData", sep = '/')))

## result 
dir_rst <- "~/Desktop/Dissertation/ODFL/Simulation Revision/Result"
rst <- get(load(paste(dir_rst, "Rst.50sz.NB.BP.Qua.200sim.RData", sep = '/')))
rst.bnp <- get(load(paste(dir_rst, "Rst.50sz.NB.Qua.BNP.200sim.RData", sep = '/')))


## using second data set 
Y <- Yall_50[[2]]
Z_all <- Zall_50[[2]]
Xm <- Xmall_50[[2]]
n <- ncol(Y)
theta.bp <- rst$theta_sim[[2]]
theta.bp <- matrix(theta.bp, nrow = 3)
eta.bp <- rst$eta[2]
W.bp <- rst$W_sim[[2]]
W.bp <- mean(W.bp)
yhat <- theta.bp[1,] + theta.bp[2,]*Xm + theta.bp[3,]*(Xm)^2 + eta.bp*log(W.bp)


## making plot
plot(c(1), c(1), type='n', ylim=c(min(Y), max(Y)), xlim = c(min(Xm), max(Xm)),
     xlab = "", ylab = "")
cols <- rainbow(3, s=0.5)
for(i in 1:n){
  lines(unlist(Z_all[[i]])[,2], Y[,i], type='l', lwd=0.5, col = rgb(0.3,0.3,0.3))
}
lines(Xm, yhat, col = "#FF8080", lwd=4)

f.bnp <- rst.bnp$f_sim[[2]]
f.bnp <- matrix(f.bnp, 18)
eta.bnp <- rst.bnp$eta[2]
W.bnp <- rst.bnp$W_sim[[2]]
W.bnp <- mean(W.bnp)
yhat.bnp <- f.bnp + eta.bnp*log(W.bnp)
lines(Xm, yhat.bnp, col = "#80FF80", lwd=4)

## true line
Wtrue <- rlnorm(1, 0, sqrt(2))
ytrue <- 5 - 0.01*Xm - 0.05*(Xm)^2 + 1*log(Wtrue)
lines(Xm, ytrue, lwd = 4, col="#8080FF")


## sample size 100
## data set
Yall_100 <- get(load(paste(dir_odfl, "Simulated 200 Response (Y) n=100.RData", sep = '/')))
Uall_100 <- get(load(paste(dir_odfl, "Simulated 200 Response (Uall) n=100.RData", sep = '/')))
Zall_100 <- get(load(paste(dir_odfl, "Simulated 200 Response (Zall) n=100.RData", sep = '/')))
Xmall_100 <- get(load(paste(dir_odfl, "Simulated 200 Response (Xm) n=100.RData", sep = '/')))

## read data 
## use the 2nd 
Y <- Yall_100[[2]]
Z_all <- Zall_100[[2]]
Xm <- Xmall_100[[2]]
n <- ncol(Y)


## making plot
plot(c(1), c(1), type='n', ylim=c(min(Y), max(Y)), xlim = c(min(Xm), max(Xm)),
     xlab = "", ylab = "")

for(i in 1:n){
  lines(unlist(Z_all[[i]])[,2], Y[,i], type='l', lwd=0.5, col=rgb(0.3, 0.3, 0.3))
}

## read results
rst.100 <- get(load(paste(dir_rst, "Rst.100sz.NB.BP.Qua.200sim.RData", sep = '/')))
rst.bnp <- get(load(paste(dir_rst, "Rst.100sz.NB.Qua.BNP.200sim.RData", sep = '/')))


theta.bp <- rst.100$theta_sim[[2]]
theta.bp <- matrix(theta.bp, nrow = 3)
eta.bp <- rst.100$eta[2]
W.bp <- rst.100$W_sim[[2]]
W.bp <- mean(W.bp)
yhat <- theta.bp[1,] + theta.bp[2,]*Xm + theta.bp[3,]*(Xm)^2 + eta.bp*log(W.bp)

lines(Xm, yhat, col = '#FF8080', lwd=4)

f.bnp <- rst.bnp$f_sim[[2]]
f.bnp <- matrix(f.bnp, 18)
eta.bnp <- rst.bnp$eta[2]
W.bnp <- rst.bnp$W_sim[[2]]
W.bnp <- mean(W.bnp)
yhat.bnp <- f.bnp + eta.bnp*log(W.bnp)
lines(Xm, yhat.bnp, col = "#80FF80", lwd=4)

## true line
Wtrue <- rlnorm(1, 0, sqrt(2))
ytrue <- 5 - 0.01*Xm - 0.05*(Xm)^2 + 1*log(Wtrue)
lines(Xm, ytrue, lwd = 4, col="#8080FF")


