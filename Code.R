library(quantmod)
library(tseries)
library(sde)
library(fOptions)
library(stats4)
library(ggplot2)
library(VarianceGamma)
library(viridis)
library(moments)
library(PerformanceAnalytics)
library(fBasics)
library(Runuran)
library(zoo)
library(Quandl)
library(gsl)
library(copula)
library(qrng)
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(VineCopula)
library(TTR)
set.seed(123)

#User defined funtion made by me to have an overview of our return
report <- function(vec){
  
  count = length(vec)
  #central tendency
  mean = mean(vec)
  median = median(vec)
  
  #variability
  variance = var(vec)
  stdev = sd(vec)
  stdev_over_mean = stdev / mean
  skewness = e1071::skewness(vec)
  kurtosis = e1071::kurtosis(vec)
  
  #rank statistics
  min = min(vec)
  max = max(vec)
  range = max - min
  
  table = round(t(data.frame(count,
                             mean,
                             median,
                             variance,
                             stdev,
                             stdev_over_mean,
                             skewness,
                             kurtosis,
                             min,
                             max,
                             range)), 7)
  
  return(table)
}

getSymbols("KO", from="2010-01-01", to="2020-07-05")
lineChart(KO$KO.Adjusted,name="Coca-Cola",theme="black") 

COKE <- get.hist.quote("KO", start = "2010-01-01", end = "2020-07-05")
chartSeries(COKE, TA=c(addVo(), addBBands()), theme="black")
COKE <- COKE$Close
cpoint(COKE)
chartSeries(KO$KO.Adjusted, TA=c(addVo(), addBBands()), theme="black")
addVLine = function(dtlist) plot(addTA(xts(rep(TRUE,NROW(dtlist)),dtlist),on=1, col="red"))
#User defined function to plot the line where the cpoint begins
addVLine(cpoint(COKE)$tau0)

#all log returns
COKE <- as.numeric(COKE)
n <- length(COKE)
X <- log(COKE[-1]/COKE[-n])
plot(X, type = "l", main = "Coca-Cola stock log returns")
abline(v = cpoint(X)$tau0, col = "red")

#Now let's focus from the changing point
getSymbols("KO", from="2018-10-019", to="2020-07-05")
lineChart(KO$KO.Adjusted,name="Coca-Cola",theme="black") 
chartSeries(KO) #candlestick chart, can be optimized via plotly package
Coca <- KO$KO.Adjusted
X <- diff(log(Coca))
plot(X)
head(X)
sum(is.na(X)) #There should be only the first observation as NA
X <- na.omit(X)
head(X)
report(X)

#Beta of Cocacola during the period considered
#Beta 
getSymbols("^DJI", from="2018-01-01", to="2020-07-05")
Dow = DJI$DJI.Adjusted
D = na.omit(diff(log(Dow)))
CAPM.beta(X, D)

#VAR 99%
VaR(R = X$KO.Adjusted, p = 0.99, method = "historical")

#Sharpe Ratio
SharpeRatio(X$KO.Adjusted, Rf = 0.007, p = 0.99)

#Cumulative returns
#I didn't put this part in the project
{
  X$KO.Adjusted[1,1] = 0 #-0.006198556
  head(X)
  cum1 = sum(X$KO.Adjusted)
  cum = exp(cum1) - 1; cum

  XX = Delt(Coca)
  head(XX)
  XX[1,1] = 0
  XX$gross = 1+ XX$Delt.1.arithmetic
  XX$cump = cumprod(XX$gross)
  y.range = range(XX$cump)
  plot(XX$cump, type= "l", auto.grid = FALSE, xlab = "Date", ylab = "Value of investment ($)",
      ylim = y.range, main = "Coca-Cola performance based on total returns")
  ggplot(XX, aes(x = index(XX), y = XX$cump)) + geom_line(color = "red") + ggtitle("Coca-Cola performance based on total returns") + xlab("Date") + ylab("Value of investment ($)")

  X$KO.Adjusted[1,1] = -0.006198556
}

#Let's plot the density function and the qqplot in order to check if it looks 
#like a Gaussian distribution

#density plot:
plot(density(X), lwd=2, main="Coca-Cola stock log returns density plot")
f <- function(u) dnorm(u, mean=mean(X), sd=sd(X))
curve( f, -0.1, 0.1, add=TRUE, col="red",lwd=2)

#qq plot:
qqnorm(X, main="Coca-Cola stock Q-Q plot")
qqline(X, col="red",lwd=2)


####PARAMETER ESTIMATION####

#Historical volatility
Delta<- 1/252
alpha.hat <- mean(X,na.rm=TRUE)/Delta 
sigma.hat <- sqrt(var(X,na.rm=TRUE)/Delta) 
mu.hat<- alpha.hat + 0.5*sigma.hat^2 
mu.hat
sigma.hat

#Implied volatility
dim(X)
S0 <- as.numeric(Coca[441])
K <- 35 #deep in the money
T <- 14 * Delta
r <- 0.07
p <- 9.4 #KO200710C00035000	(yahoo)
sigma.imp1 <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.imp1 #result too much low

K <- 44.5 #moderately in the money
T <- 8 * Delta
r <- 0.07
p <- 0.76 #KO200710C00044500 (yahoo)
sigma.imp2 <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.imp2 #"more acceptable" result

K <- 42 #moderately in the money
T <- 8 * Delta
r <- 0.07
p <- 3.1 #KO200710C00042000 (yahoo)
sigma.imp3 <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.imp3 

K <- 49 #out of the money
T <- 8 * Delta
r <- 0.07
p <- 0.02 #KO200710C00042000 (yahoo)
sigma.impout <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.impout

K <- 42.5 
T <- 167 * Delta #let's try to expand the time horizon
r <- 0.07
p <- 4.5 #KO201218C00042500 (yahoo)
sigma.imp4 <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.imp4

K <- 20 #deep in the money
T <- 183 * Delta #let's try to expand the time horizon
r <- 0.07
p <- 26.8 #KO201218C00020000 (yahoo)
sigma.imp5 <- GBSVolatility(p, "c", S = S0, X = K, Time = T, r = r,  b = r) 
sigma.imp5


#MLE estimation
set.seed(123)
log.lik <- function(mu = mu.hat, sigma = sigma.hat) -sum(dnorm(X, mean = mu,
                                                               sd = sigma, log = TRUE)) 
fit <- mle(log.lik, lower = c(0.01,0.01), method = "L-BFGS-B") #"L-BFGS-B"
fit
fit2 <- mle(log.lik,  method = "Brent")  #Brent, gives error
fit2
fit3 <- mle(log.lik)
fit3
fit4 <- mle(log.lik,  method = "BFGS") #same as fit4, in fact it is the default algorithm for mle function in stats4
fit4
fit5 <- mle(log.lik,  method = "CG"); fit5
set.seed(123)
fit6 = mle (log.lik, method = "SANN"); fit6
#MLE didn't reach convergence
#Furthermore for the other algo it returns some warnings

### European Option pricing
#With Black & Scholes formula
S0 <- as.numeric(Coca[441])
sigma.hat <- as.numeric(sigma.hat)
K <- 42 #moderately in the money
T <- 8 * Delta
r <- 0.07
p <- 3.1 #KO200710C00042000 (yahoo)
p0 <- GBSOption("c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma.hat)@price
p0

#with MC
MCPrice <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, 
                    M = 1000, f) {
  h <- function(m) {
    u <- rnorm(m/2)
    tmp <- c(x * exp((r - 0.5 * sigma^2) * (T - t) + sigma * 
                       sqrt(T - t) * u), x * exp((r - 0.5 * sigma^2) * (T - 
                                                                          t) + sigma * sqrt(T - t) * (-u)))
    mean(sapply(tmp, function(xx) f(xx)))
  }
  p <- h(M)
  p * exp(-r * (T - t))
}

f <- function(x) max(0, x - K)

set.seed(123)
M <- 1000
MCPrice(x = S0, t = 0, T = T, r = r, sigma.hat, M = M, f = f)
M <- 50000
MCPrice(x = S0, t = 0, T = T, r = r, sigma.hat, M = M, f = f)
M <- 1e+06
MCPrice(x = S0, t = 0, T = T, r = r, sigma.hat, M = M, f = f)

#Speed of convergence
set.seed(123)
m <- c(10, 50, 100, 150, 200, 250, 500, 1000)
p1 <- NULL
err <- NULL
nM <- length(m)
repl <- 100
mat <- matrix(, repl, nM)
for (k in 1:nM) {
  tmp <- numeric(repl)
  for (i in 1:repl) tmp[i] <- MCPrice(x = S0, t = 0, T = T, 
                                      r = r, sigma.hat, M = m[k], f = f)
  mat[, k] <- tmp
  p1 <- c(p1, mean(tmp))
  err <- c(err, sd(tmp))
}
colnames(mat) <- m

p0 <- GBSOption(TypeFlag = "c", S = S0, X = K, Time = T, r = r, b = r, sigma = sigma.hat)@price
minP <- min(p1 - err)
maxP <- max(p1 + err)
plot(m, p1, type = "n", ylim = c(minP, maxP), axes = F, ylab = "MC price",  xlab = "MC replications")
lines(m, p1 + err, col = "blue")
lines(m, p1 - err, col = "blue")
axis(2, p0, "B&S price")
axis(1, m)
boxplot(mat, add = TRUE, at = m, boxwex = 15, col = "orange",  axes = F)
points(m, p1, col = "blue", lwd = 3, lty = 3)
abline(h = p0, lty = 2, col = "red", lwd = 3)

#FFT
FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
  m <- r - log(phi(-(0+1i)))
  phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
  psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 
                                                     1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * 
                                                                                                        alpha + 1) * v)
  lambda <- (2 * pi)/(N * eta)
  b <- 1/2 * N * lambda
  ku <- -b + lambda * (0:(N - 1))
  v <- eta * (0:(N - 1))
  tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - 
                                                 ((1:N) - 1 == 0))/3
  ft <- fft(tmp)
  res <- exp(-alpha * ku) * ft/pi
  inter <- spline(ku, Re(res), xout = log(K/S0))
  return(inter$y * S0)
}
mu=1
sigma <- 0.25
phiBS <- function(u) exp((0+1i) * u * (mu - 0.5 * sigma^2) - 0.5 * sigma^2 * u^2)
FFTcall.price(phiBS, S0 = S0, K = K, r = r, T = T)

#put-call parity
call.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, 
                       K = 1) {
  d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma * 
                                                      sqrt(T - t))
  d1 <- d2 + sigma * sqrt(T - t)
  x * pnorm(d1) - K * exp(-r * (T - t)) * pnorm(d2)
}
put.price <- function(x = 1, t = 0, T = 1, r = 1, sigma = 1, 
                      K = 1) {
  d2 <- (log(x/K) + (r - 0.5 * sigma^2) * (T - t))/(sigma * 
                                                      sqrt(T - t))
  d1 <- d2 + sigma * sqrt(T - t)
  K * exp(-r * (T - t)) * pnorm(-d2) - x * pnorm(-d1)
}
C <- call.price(x = S0, t = 0, T = T, r = r, K = K, sigma = sigma.hat)
P <- put.price(x = S0, t = 0, T = T, r = r, K = K, sigma = sigma.hat)
P
C = as.numeric(C)
C - S0 + K * exp(-r * T) 

#Greeks
GBSCharacteristics(TypeFlag = "c", S = S0, X = K, Time = T,
                   r = r, b = r, sigma = sigma.hat) 


###AMerican Option pricing
# Broadie and Glasserman Monte Carlo method

simTree <- function(b,d, S0, sigma, T, r){
  tot <- sum(b^(1:(d-1)))
  S <- numeric(tot+1) 
  S[1] <- S0
  dt <- T/d
  for(i in 0:(tot - b^(d-1))){
    for(j in 1:b){
      S[i*b+j +1] <- S[i+1]*exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1))
    }
  }
  S
}



upperBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  P <- S
  P[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  for(i in tot1:0){
    m <- mean(P[i*b+1:b+1])
    v <- f(S[i+1])
    P[i+1] <- max(v,m)
  }
  P
}

lowerBG <- function(S, b, d, f){
  tot <- sum(b^(1:(d-1)))
  start <- tot - b^(d-1) +1
  end <- tot +1
  p <- S 
  p[start:end] <- f(S[start:end])
  tot1 <- sum(b^(1:(d-2)))
  
  m <- numeric(b)
  for(i in tot1:0){
    v <- f(S[i+1])
    for(j in 1:b){
      m[j] <- mean(p[i*b+(1:b)[-j]+1])
      m[j] <- ifelse( v>m[j], v, p[i*b+(1:b)[j]+1])
    }
    p[i+1] <- mean(m)
  }
  p
}


b <- 6
d <- 4
M <- 10000
low <- 0
upp <- 0
f <- function(x) sapply(x, function(x) max(x-K,0))
sigma = sigma.hat
set.seed(123)
for(i in 1:M){
  S <- simTree(b,d, S0, sigma=sigma.hat, T, r)
  low <- low + lowerBG(S, b,d, f)[1]
  upp <- upp + upperBG(S, b,d, f)[1]
}
low/M
upp/M

#via regression
LSM <- function(n, d, S0, K, sigma, r, T){
  s0 <- S0/K
  dt <- T/d
  z <- rnorm(n)
  s.t <- s0*exp((r-1/2*sigma^2)*T+sigma*z*(T^0.5))
  s.t[(n+1):(2*n)] <- s0*exp((r-1/2*sigma^2)*T-sigma*z*(T^0.5))
  CC <- pmax(1-s.t, 0)
  payoffeu <- exp(-r*T)*(CC[1:n]+CC[(n+1):(2*n)])/2*K
  euprice <- mean(payoffeu)
  
  for(k in (d-1):1){
    z <- rnorm(n)
    mean <- (log(s0) + k*log(s.t[1:n]))/(k+1)
    vol <- (k*dt/(k+1))^0.5*z
    s.t.1 <- exp(mean+sigma*vol)
    mean <- (log(s0) + k*log( s.t[(n+1):(2*n)] )) / ( k + 1 )
    s.t.1[(n+1):(2*n)] <- exp(mean-sigma*vol)
    CE <- pmax(1-s.t.1,0)
    idx<-(1:(2*n))[CE>0]
    discountedCC<- CC[idx]*exp(-r*dt)
    basis1 <- exp(-s.t.1[idx]/2)
    basis2 <- basis1*(1-s.t.1[idx])
    basis3 <- basis1*(1-2*s.t.1[idx]+(s.t.1[idx]^2)/2)
    
    p <- lm(discountedCC ~ basis1+basis2+basis3)$coefficients
    estimatedCC <- p[1]+p[2]*basis1+p[3]*basis2+p[4]*basis3
    EF <- rep(0, 2*n)
    EF[idx] <- (CE[idx]>estimatedCC)
    CC <- (EF == 0)*CC*exp(-r*dt)+(EF == 1)*CE
    s.t <- s.t.1
  }
  
  payoff <- exp(-r*dt)*(CC[1:n]+CC[(n+1):(2*n)])/2
  usprice <- mean(payoff*K)
  error <- 1.96*sd(payoff*K)/sqrt(n)
  earlyex <- usprice-euprice
  data.frame(usprice, error, euprice)
}
sigma = sigma.hat
set.seed(123)
LSM(100000, 5, S0, K, sigma, r, T)

#explict finite difference method
AmericanPutExp <- function(Smin=38, Smax,  T=1, N=10, M=10, K, r=0.07, sigma=0.01){ #note that I used 0.07 as risk free rate, the same I always used
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)
  A <- function(j) (-0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt) 
  B <- function(j) (1-sigma^2*j^2*Dt)/(1+r*Dt) 
  C <- function(j) (0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt)
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,0:N+1] <- sapply(rev(S), function(x) max(K-x,0))
  optTime <- matrix(FALSE, M+1, N+1)
  optTime[M+1,] <- TRUE
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  
  for(i in (N-1):0){
    for(j in 1:(M-1)){
      J <- M+1-j
      I <- i+1
      P[J,I] <- A(j)*P[J+1,I+1] + B(j)*P[J,I+1] + C(j)*P[J-1,I+1]
      if(P[J,I] < P[J,N+1])
        optTime[J,I] <- TRUE
    }
  }
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}

plot.AmericanPut <- function( obj ){
  plot(range(obj$t),range(obj$S),type="n",axes=F,xlab="t", ylab="S")
  axis(1,obj$t,obj$t)
  axis(2,obj$S,obj$S)
  abline(v = obj$t, h = obj$S, col = "darkgray", lty = "dotted")
  for(i in 0:obj$N){
    for(j in 0:obj$M){
      J <- obj$M+1-j
      I <- i+1
      cl <- "green"; 
      if(obj$optTime[J,I])
        cl <- "black"
      text(obj$t[i+1],obj$S[j+1], round(obj$P[J,I],2),cex=0.75, col=cl)
    }
  }
  DS <- mean(obj$S[1:2])
  y <- as.numeric(apply(obj$optTime,2, function(x) which(x)[1]))
  lines(obj$t, obj$S[obj$M+2-y]+DS, lty=2)
}
put <- AmericanPutExp(Smax =50, sigma = sigma.hat, K = 42, T=T) 
round(put$P,2) 
par(mar=c(3,3,1,1))
plot(put)

myval <- round(put$P[which(rownames(put$P)==S0),1],2)


#Implicit finite difference method

AmericanPutImp <- function( Smin=38, Smax,  T=1, N=10, M=10, K, r=0.07, sigma=0.01){
  Dt = T/N 
  DS = (Smax-Smin)/M
  t <- seq(0, T, by =Dt) 
  S <- seq(Smin, Smax, by=DS)
  
  A <- function(j) 0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt 
  B <- function(j) 1+sigma^2*j^2*Dt+r*Dt
  C <- function(j) -0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt
  
  a <- sapply(0:M, A)
  b <- sapply(0:M, B)
  c <- sapply(0:M, C)
  
  P <- matrix(, M+1, N+1)
  colnames(P) <- round(t,2)
  rownames(P) <- round(rev(S),2)
  
  P[M+1, ] <- K   # C(,j=0) = K
  P[1,] <- 0   # C(,j=M) = 0
  P[,0:N+1] <- sapply(rev(S), function(x) max(K-x,0))
  
  AA <- matrix(0, M-1, M-1)
  for(j in 1:(M-1)){
    if(j>1) AA[j,j-1] <- A(j)
    if(j<M) AA[j,j] <- B(j)
    if(j<M-1) AA[j,j+1] <- C(j)
  }
  
  optTime <- matrix(FALSE, M+1, N+1)
  for(i in (N-1):0){
    I <- i+1
    bb <- P[M:2,I+1]
    bb[1] <- bb[1]-A(1)*P[M+1-0,I+1]
    bb[M-1] <- bb[M-1]-C(M-1)*P[M+1-M,I+1] 
    P[M:2,I] <- solve(AA,bb)
    idx <- which(P[,I] < P[,N+1])
    P[idx,I] <- P[idx,N+1] 
    optTime[idx, I] <- TRUE
  }
  optTime[M+1,] <- TRUE 
  optTime[which(P[,N+1]>0),N+1] <- TRUE
  colnames(optTime) <- colnames(P)
  rownames(optTime) <- rownames(P)
  ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
  class(ans) <- "AmericanPut"
  return(invisible(ans))
}



put <- AmericanPutImp(Smax = 50, sigma = sigma.hat, K= 42, T=T)
round(put$P,2)


par(mar=c(3,3,1,1))
plot(put)


###Lévy processes###

#Variance Gamma
vgFit(X) #esitmate VG parameters on the sample
str(vgFit(X))
vg_param <- vgFit(X)$param

c <- as.numeric(vg_param[1])
sigma <- as.numeric(vg_param[2])
theta <- as.numeric(vg_param[3])
nu <- as.numeric(vg_param[4])

N <- 100    
nsim <- 1000

#Variance Gamma function
VG=function(sigma, nu, theta, T, N, r) { 
  a=1/nu 
  b=1/nu 
  h=T/N
  t=(0:N)*T/N 
  X=rep(0, N+1) 
  I=rep(0,N) 
  X[1]=0 
  for(i in 1:N) { 
    I[i]=rgamma(1,a*h,b) 
    X[i+1]=X[i] + theta*I[i]+sigma*sqrt(I[i])*rnorm(1)
  }
  return((X)) }


set.seed(123)
VG_paths<-matrix(nrow = nsim, ncol=N+1)  
for(i in 1:nsim){                        
  VG_paths[i,]<-VG(sigma, nu, theta, T, N, r)
}

VG_paths


colori=viridis(nsim)
plot(VG_paths[1,], col=0, type="l", ylim = c(min(VG_paths),max(VG_paths)), 
     main = "Monte Carlo Simulation for VG returns", sub = "100 steps, 1000 paths", 
     xlab = "Time", ylab = "VG returns")
for(i in 2:nsim){
  lines(VG_paths[i,], col=colori[i], lwd = 2);
}

l_ret.s <- sort(as.numeric(X))  #sort the log returns

p <- ppoints(length(l_ret.s)) #plotting position

VG.q <- qvg(p, vgC=c, sigma=sigma, theta=theta, nu=nu) #compute the quantile

plot(VG.q, l_ret.s, main = "Variance-Gamma Q-Q Plot",
     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")

par(mfrow = c(2,1))

plot(density(X[-1,]), type = "l", lwd = 2, lty = 3, col = "coral2",
     xlim= c(-0.03,0.03), ylim=c(0,120), main ="", xlab ="", ylab = "")
legend ("topright", inset = .02, c("Kernel", "VG"),
        col=c("coral2","seagreen3"), lwd=c(2,1), lty=c(3,1), cex = 0.8, bty = "n")
points(seq(min(X[-1,]), max(X[-1,]), length.out=500), 
       dvg(seq(min(X[-1,]), max(X[-1,]), length.out=500),
           mean(X[-1,]), sd(X[-1,])), type="l", col="seagreen3")

#Log-density comparison
gridplot <- seq(min(X[-1,]), max(X[-1,]), length.out=length(X[-1,]))
plot(density(X[-1,])$x, log(density(X[-1,])$y), 
     type = "l", lwd = 2, lty = 3, col = "coral2",
     xlim= c(-0.025,0.025), ylim=c(-10,7.5),
     main ="", xlab ="", ylab = "")
legend ("topright", inset = .02, c("Kernel", "Normal"),
        col=c("coral2","seagreen3"), lwd=c(2,1), lty=c(3,1), cex = 0.8, bty = "n")
points(gridplot, log(dvg(gridplot, param = c(c, sigma, theta, nu))), 
       type="l", col="seagreen3")

par(mfrow=c(1,1))

#Hypotesis testing
#Chisquared test
chisq.test(l_ret.s, VG.q)
#K-S test
ks.test(as.numeric(X), rvg(length(as.numeric(X)), 
                           param = c(c, sigma, theta, nu)))
#summary statistics
final_retVG<-VG_paths[,N]
basicStats(final_retVG)
hist(final_retVG)


VGexp=function(sigma, nu, theta, T, N, r, S0) { 
  a=1/nu 
  b=1/nu 
  h=T/N
  t=(0:N)*T/N 
  X=rep(0, N+1) 
  I=rep(0,N) 
  X[1]=S0 
  for(i in 1:N) { 
    I[i]=rgamma(1,a*h,b) 
    X[i+1]=X[i]*exp(r*t+theta*I[i]+sigma*sqrt(I[i])*rnorm(1))
  }
  return(X)}


set.seed(123)
VGexp_paths<-matrix(nrow = nsim, ncol=N+1)
for(i in 1:nsim){
  VGexp_paths[i,]<-VGexp(sigma, nu, theta, T, N, r, S0)
}

VGexp_paths

plot(VGexp_paths[1,], col=0, type="l", ylim = c(min(VGexp_paths),max(VGexp_paths)), 
     main = "MC Simlation for VG stock prices", sub = "100 steps, 10 paths", 
     xlab = "Time", ylab = "Coke")
for(i in 2:nsim){
  lines(VGexp_paths[i,], col=colori[i], lwd = 2);
  
}
final_pricesVG<-VGexp_paths[,N]

#mean correcting martingale method on vg
rn_final_pricesVG<-S0*(final_pricesVG)*(exp(r*T)/(mean(final_pricesVG)))
rn_final_pricesVG

basicStats(rn_final_pricesVG)
hist(rn_final_pricesVG)


payoff_VG <- pmax(rn_final_pricesVG - K, 0)
optprice_VG <- mean(payoff_VG)*exp(-r*T)
optprice_VG

#FFT
alpha <- 1.65
phiVG <- function(u) {
  omega <- (1/nu) * (log(1 - theta * nu - sigma^2 * nu/2))
  tmp <- 1 - (0+1i) * theta * nu * u + 0.5 * sigma^2 * u^2 * nu
  tmp <- tmp^(-1/nu)
  exp((0+1i) * u * log(S0) + u * (r + omega) * (0+1i)) * tmp
}

FFTcall.price <- function(phi, S0, K, r, T, alpha = 1, N = 2^12, eta = 0.25) {
  m <- r - log(phi(-(0+1i)))
  phi.tilde <- function(u) (phi(u) * exp((0+1i) * u * m))^T
  psi <- function(v) exp(-r * T) * phi.tilde((v - (alpha + 
                                                     1) * (0+1i)))/(alpha^2 + alpha - v^2 + (0+1i) * (2 * 
                                                                                                        alpha + 1) * v)
  lambda <- (2 * pi)/(N * eta)
  b <- 1/2 * N * lambda
  ku <- -b + lambda * (0:(N - 1))
  v <- eta * (0:(N - 1))
  tmp <- exp((0+1i) * b * v) * psi(v) * eta * (3 + (-1)^(1:N) - 
                                                 ((1:N) - 1 == 0))/3
  ft <- fft(tmp)
  res <- exp(-alpha * ku) * ft/pi
  inter <- spline(ku, Re(res), xout = log(K/S0))
  return(inter$y * S0)
}


FFTcall.price(phiVG, S0 = S0, K = K, r = r, T = T)



#Meixner model 

#Unico metodo che funziona è la risk neutral trasformation
x   <-mean(X, na.rm = TRUE)
y   <-sd(X, na.rm = TRUE)
z <-as.numeric(skewness(X, na.rm = TRUE))
w <-as.numeric(kurtosis(X, na.rm = TRUE))

m <- x-((z*sqrt(y))/(w-(z^2)-3)) 
a <- sqrt(y*(2*w-3*(z^2)-6)) 
b <- 2*atan(-sqrt((z^2)/(2*w-3*(z^2)-6))) 
d <- 1/(w-(z^2)-3) 
nsim = 100
N <- 100


#qq plot
MX.q <- uq(pinvd.new(udmeixner(a, b, d, m)), p) #compute the quantile

plot(MX.q, l_ret.s, main = "Meixner Q-Q Plot",
     xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")

#esscher transform 
theta <- -1/a * (b + 2 * atan((-cos(a/2)+ exp((m-r)/2*d))/sin(a/2)))
b <- a*theta+b #it doesn't perform well

MX=function(a, b, d, m, N) {
  distr <- udmeixner(a, b, d, m) #meixner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  rdmMXgen <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  X=rep(0,N+1)
  for (i in 1:N){
    X[i+1]=X[i]+rdmMXgen[i]
  }
  return(X)
}

MX(a, b, d, m, N)


set.seed(123)
MX_paths<-matrix(nrow = nsim, ncol=N+1) #fill the matrix with random paths that follow
for(i in 1:nsim){                        #the function MX just created
  MX_paths[i,]<-MX(a,b,d,m,N)
}

MX_paths

plot(MX_paths[1,], col=0, type="l", ylim = c(min(MX_paths),max(MX_paths)), 
     main = "Monte Carlo Simulation for Meixner returns", sub = "100 steps, 100 paths", 
     xlab = "Time", ylab = "MXNR returns")
for(i in 2:nsim){
  lines(MX_paths[i,], col=colori[i], lwd = 2);
}

 

final_retMX<-MX_paths[,N]
basicStats(final_retMX)
hist(final_retMX) 

#function for stock price with Meixner returns
MXexp=function(a, b, d, m, N, T, r, S0) {
  distr <- udmeixner(a, b, d, m) #meixner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  generazioni <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  t=(0:N)*T/N
  X=rep(0,N+1)
  X[1]=S0
  for (i in 1:N){
    X[i+1]=X[i]*exp(r*t+generazioni[i])
  }
  return(X)
}


MXexp(a, b, d, m, N, T, r, S0)


set.seed(123)
MXexp_paths<-matrix(nrow = nsim, ncol=N+1)
for(i in 1:nsim){
  MXexp_paths[i,]<-MXexp(a,b,d,m,100,T,r,S0) #vengono tutte le linee uguali perch? MX non varia!
}

MXexp_paths

final_pricesMX<-MXexp_paths[,N]

payoff_MX <- pmax(final_pricesMX - K, 0)

optprice_MX <- mean(payoff_MX)*exp(-r*T)

optprice_MX



payoff_MXess <- pmax(rn_final_pricesMX - K, 0)
optprice_MXess <- mean(payoff_MXess)*exp(-r*T)
optprice_MXess


#mean correcting
#m <- r -2 *d*log(cos(b/2)/cos((a+b)/2)) is a mess
x   <-mean(X, na.rm = TRUE)
y   <-sd(X, na.rm = TRUE)
z <-as.numeric(skewness(X, na.rm = TRUE))
w <-as.numeric(kurtosis(X, na.rm = TRUE))

m <- x-((z*sqrt(y))/(w-(z^2)-3)) 
a <- sqrt(y*(2*w-3*(z^2)-6)) 
b <- 2*atan(-sqrt((z^2)/(2*w-3*(z^2)-6))) 
d <- 1/(w-(z^2)-3) 
nsim = 100
N <- 100

MX=function(a, b, d, m, N) {
  distr <- udmeixner(a, b, d, m) #meixner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  rdmMXgen <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  X=rep(0,N+1)
  for (i in 1:N){
    X[i+1]=X[i]+rdmMXgen[i]
  }
  return(X)
}

MX(a, b, d, m, N)


set.seed(123)
MX_paths<-matrix(nrow = nsim, ncol=N+1) #fill the matrix with random paths that follow
for(i in 1:nsim){                        #the function MX just created
  MX_paths[i,]<-MX(a,b,d,m,N)
}

MX_paths
plot(MX_paths[1,], col=0, type="l", ylim = c(min(MX_paths),max(MX_paths)), 
     main = "Monte Carlo Simulation for Meixner returns", sub = "100 steps, 100 paths", 
     xlab = "Time", ylab = "MXNR returns")
for(i in 2:nsim){
  lines(MX_paths[i,], col=colori[i], lwd = 2);
}



final_retMX<-MX_paths[,N]
basicStats(final_retMX)
hist(final_retMX) 

#function for stock price with Meixner returns
MXexp=function(a, b, d, m, N, T, r, S0) {
  distr <- udmeixner(a, b, d, m) #meixner distribution
  gen <- pinvd.new(distr) #Polynomial interpolation of INVerse CDF
  generazioni <- ur(gen,N) #randomly draws N objects from gen (from a Meixner distr)
  h=T/N
  t=(0:N)*T/N
  X=rep(0,N+1)
  X[1]=S0
  for (i in 1:N){
    X[i+1]=X[i]*exp(r*t+generazioni[i])
  }
  return(X)
}


MXexp(a, b, d, m, N, T, r, S0)


set.seed(123)
MXexp_paths<-matrix(nrow = nsim, ncol=N+1)
for(i in 1:nsim){
  MXexp_paths[i,]<-MXexp(a,b,d,m,100,T,r,S0) #vengono tutte le linee uguali perch? MX non varia!
}

MXexp_paths

final_pricesMX<-MXexp_paths[,N]

rn_final_pricesMX<-S0*(final_pricesMX)*(exp(r*T)/(mean(final_pricesMX)))

rn_final_pricesMX

payoff_MXmean <- pmax(rn_final_pricesMX - K, 0)

optprice_MXmean <- mean(payoff_MXmean)*exp(-r*T)

optprice_MXmean


### Multi asset options ###
#Rho correlation
GBM <- function(N, sigma, mu, S0, Wt = NULL) {
  if (is.null(Wt)) {
    Wt <- cumsum(rnorm(N, 0, 1))
  }
  t <- (1:N)/252
  p1 <- (mu - 0.5*(sigma*sigma)) * t
  p2 <- sigma * Wt
  St = S0 * exp(p1 + p2)
  return(St)
}

CorrelatedGBM <- function(N, S0, mu, sigma, cor.mat) {
  mu <- as.matrix(mu)
  sigma <- as.matrix(sigma)
  GBMs <- matrix(nrow = N, ncol = nrow(mu))
  Wt <- matrix(rnorm(N * nrow(mu), 0, 1), ncol = nrow(mu))
  Wt <- apply(Wt, 2, cumsum)
  chol.mat <- chol(cor.mat) # upper triangular cholesky decomposition
  Wt <- Wt %*% chol.mat   # key trick for creating correlated paths
  for (i in 1:nrow(mu)) {
    GBMs[,i] <- GBM(N, sigma[i], mu[i] , S0[i], Wt[, i])
  }
  return(GBMs)
}

GetPrices <- function(tickers, startDate='1992-01-02') {
  prices <- get.hist.quote(instrument = tickers[1], start = startDate, 
                           quote = 'AdjClose')
  for (tik in 2:length(tickers)) {
    tmp <- get.hist.quote(instrument = tickers[tik], 
                          start = start, quote = 'AdjClose')
    prices <- merge(prices, tmp)
  }
  return(prices)
}

set.seed(123)
N <- 2 * 252 
t <- (1:N)/252
start <- '2002-1-1'
tickers <- c('KO', 'AXP')
prices <- GetPrices(tickers, start)

returns.mat <- as.matrix(na.omit(diff(log(prices))))
mean.vec <- as.numeric(colMeans(returns.mat))
sigma.vec <- as.numeric(sqrt(apply(returns.mat, 2, var)))
prices.vec <- as.numeric(prices[nrow(prices)])
cor.mat <-as.matrix(cor(returns.mat))

paths <- CorrelatedGBM(N, prices.vec , mean.vec, sigma.vec, cor.mat)

colors <- c('red', 'darkgreen')
plot(t, paths[,1], type = 'l', ylim = c(0, max(paths)), xlab = "Year", 
     ylab = "Price", main = "Simulated Coca-Cola and American Express Prices", col = colors[1])
for (i in 2:ncol(paths)) {
  lines(t, paths[, i], col = colors[i])
}
legend(x = 0.5, y = 25, c('KO', 'AXP'), lty = c(1,1), col = colors, cex = 0.7)

cor.mat


#Copula

#loglikCopula(th.C,returns.mat, claytonCopula())
clay <- fitCopula(claytonCopula(dim = 2), returns.mat, method ="itau")
gumb <- fitCopula(gumbelCopula(dim = 2), returns.mat, method ="itau")
frank <- fitCopula(frankCopula(dim = 2), returns.mat, method ="itau")
summary(clay)
summary(gumb)
summary(frank)


getSymbols("KO", from="2018-10-019", to="2020-07-05")
Coca = KO$KO.Adjusted
B = na.omit(diff(log(Coca)))


getSymbols("AXP", from="2018-10-019", to="2020-07-05")
amex = AXP$AXP.Adjusted
C = na.omit(diff(log(amex)))


M = cbind(B,C)
corKendall(as.matrix(M))

risk.measures <- function(x, alpha)
{
  if(!is.matrix(x)) x <- rbind(x)
  n <- nrow(x)
  d <- ncol(x)
  
  aloss  <- rowSums(x) # n-vector of aggregated losses
  VaR <- quantile(aloss, probs=alpha, names=FALSE) # VaR estimate
  l <- sum(xcd <- aloss > VaR)
  ES <- mean(aloss[xcd]) / (1-alpha) # ES estimate
  
  Alloc.first <- x[aloss>=VaR, 1] %*% rep(1/n, l) / (l/n) # capital allocated to X_1
  Alloc.mid   <- x[aloss>=VaR, floor(d/2)] %*% rep(1/n, l) / (l/n) # capital allocated to X_{floor(d/2)}
  Alloc.last  <- x[aloss>=VaR,d] %*% rep(1/n, l) / (l/n) # Capital Allocated to X_d
  
  ## return estimated risk measures
  c(VaR = VaR, ES = ES, Alloc.first = Alloc.first, Alloc.mid = Alloc.mid,
    Alloc.last = Alloc.last)
}

#Gumbel KO AXP
family.G <- "Gumbel"
tau = 0.3577591
alpha = 0.01
d = 2
n <- 1e5
nu = 3 #degress of freedom of the Student's, but can't find its copula in the code
th.G <- iTau(getAcop(family.G), tau) #theta, the corresponding parameter
gumbel.cop <- onacopulaL(family.G, nacList=list(th.G, 1:d))
set.seed(123)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
U.C.CDM  <- cCopula(U.CDM,  cop=gumbel.cop, inverse=TRUE)

#erT <- exp(-r*T)
rm.C.CDM <- risk.measures(U.C.CDM, alpha)

res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Gumbel", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(rm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(rm.C.CDM[2]), ncol=1)
res



#Clayton
family.C <- "Clayton"
th.C <- iTau(getAcop(family.C), tau) # corresponding parameter
clayton.cop <- onacopulaL(family.C, nacList=list(th.C, 1:d))

set.seed(123)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
ClayU.C.CDM  <- cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE)
#erT <- exp(-r*T)
Clayrm.C.CDM <- risk.measures(ClayU.C.CDM, alpha)
res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Clayton", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(Clayrm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(Clayrm.C.CDM[2]), ncol=1)
res

#Frank
family.F <- "Frank"
th.F <- iTau(getAcop(family.F), tau) # corresponding parameter
frank.cop <- onacopulaL(family.F, nacList=list(th.F, 1:d)) # Frank copula
### 2.2 Sampling 
set.seed(123)
U.CDM  <- matrix(runif(n*d), ncol=d) # pseudo
FrankU.C.CDM  <- cCopula(U.CDM,  cop=frank.cop, inverse=TRUE) # pseudo
### 2.3 Functional Calculation
#erT <- exp(-r*T)
Frankrm.C.CDM <- risk.measures(FrankU.C.CDM, alpha)
### 2.4 Results
res <- array(dim=c(2,2,1), dimnames=list(type=c(paste0("VaR.", alpha), paste0("ES.", alpha)),
                                         copula=c("Frank", paste0("t", nu)),
                                         method=c("CDM")))
res[paste0("VaR.", alpha),,] <- matrix(c(Frankrm.C.CDM[1]), ncol=1)
res[paste0("ES.", alpha),,]  <- matrix(c(Frankrm.C.CDM[2]), ncol=1)
res

#to compare the time consuming
system.time(cCopula(U.CDM,  cop=frank.cop, inverse=TRUE))
system.time(cCopula(U.CDM,  cop=gumbel.cop, inverse=TRUE))
system.time(cCopula(U.CDM,  cop=clayton.cop, inverse=TRUE))


### Multi-Asset Options based on Lévy Processes ###

#Ballotta-Bonfiglioli model

VG=function(sigma, nu, mu, T, N) { 
  a=1/nu 
  b=1/nu 
  h=T/N 
  t=(0:T)/N 
  X=rep(0, N+1) 
  I=rep(0,N) 
  X[1]=0 
  for(i in 1:N) { 
    I[i]=rgamma(1,a*h,b) 
    X[i+1]=X[i] + mu*I[i]+sigma*sqrt(I[i])*rnorm(1) 
  } 
  return(X) 
}
time <- (0:10000) / 10000

sigma1 <- 0.75
nu1 <- 0.5
mu1 <- 0.1  

sigma2 <- 0.70
nu2 <- 0.3
mu2 <- 0.2

sigma3 <- 0.60
nu3 <- 0.2
mu3 <- 0.3

VG(sigma1, nu1, mu1, 1, 10000)
VG(sigma2, nu2, mu2, 1, 10000)
VG(sigma3, nu3, mu3, 1, 10000)
set.seed(1)
plot(x=time, y=VG(sigma1, nu1, mu1, 10, 10000)+VG(sigma3, nu3, mu3, 1, 10000), type="l",  ylab="Z1 (black) and Z2 (red)", main = "Simulations based on Ballotta Bonfiglioli method")
lines(x =time, y = VG(sigma2, nu2, mu2, 1, 10000)+VG(sigma3, nu3, mu3, 1, 10000), col="red")

# Simulation of a bivariate Levy process with finite variation and
# Clayton Levy copula.
# 
# Author: Alice Pignatelli di Cerchiara

# maximum number of jumps?
m <- 2000
n <- 10
d <- 1000
# parameters of Clayton Levy copula (chosen as in Figure 1)

# Marginal (1st variance gamma) parameters
c1 <- 10
lambda1_plus <- 1
lambda1_minus <- 1
c2 <- 10
lambda2_plus <- 1
lambda2_minus <- 1

eta <- 0.9
theta <- 0.3

## NOTE: In Tankov they use a different parameterisation?
## [cf. Example 2.1; especially eq. (2.4)]

# t in grid_pts
grid_pts <- (1:d) / d # has to be subset of [0,1]

# define the inverse of the conditional distribution function of U|xi
F_inv <- function(xi, u) {
  
  # define the two functions B and C internally
  B <- function(xi, u) {
    res <- sign(u - 1 + eta) * (xi >= 0) + sign(u - eta) * (xi < 0)
    return(res)
  }
  
  C <- function(xi, u) {
    res1 <- (((u - 1 + eta)/eta) * (u >= 1 - eta) + ((1 - eta - u)/(1-eta)) * (u < 1 - eta)) * (xi >= 0)
    res2 <- (((u - eta)/(1-eta)) * (u >= eta) + ((eta-u)/eta) * (u < eta)) * (xi < 0)
    return(res1 + res2)
  }
  
  res <- B(xi,u) * abs(xi) * (C(xi,u)^(-theta/(theta+1)) - 1)^(-1/theta)
  return(res)
}
#FROM https://cran.r-project.org/web/packages/copula/vignettes/NALC.html
# tail integral for variance gamma
nu_bar_vargamma <- function(x, th, kap, sig, c, lambda_plus, lambda_minus) {
  if (!hasArg(lambda_plus)) {
    lambda_plus <- (sqrt(th^2+2*sig^2/kap)-th)/sig^2
  }
  if (!hasArg(lambda_minus)) {
    lambda_minus <- (sqrt(th^2+2*sig^2/kap)+th)/sig^2
  }
  if (!hasArg(c)) {
    c <- 1/kap
  }
  lambda <- lambda_plus * (x > 0) + lambda_minus * (x < 0)
  -c*expint_Ei(-lambda*abs(x), give=FALSE)
}

## Inverse of the tail integral of a variance gamma Levy process
nu_bar_inv_vargamma <- function(Gamma, th, kap, sig, c, lambda_plus, lambda_minus, ...)
{
  if (!hasArg(lambda_plus)) {
    lambda_plus <- (sqrt(th^2+2*sig^2/kap)-th)/sig^2
  }
  if (!hasArg(lambda_minus)) {
    lambda_minus <- (sqrt(th^2+2*sig^2/kap)+th)/sig^2
  }
  if (!hasArg(c)) {
    c <- 1/kap
  }
  
  
  max.val <- nu_bar_vargamma(.Machine$double.xmin, th=th, kap=kap, sig=sig, c=c, lambda_plus=lambda_plus, lambda_minus=lambda_minus)
  res <- numeric(length(Gamma))
  large <- (Gamma >= max.val | Gamma <= -1*max.val)
  res[large] <- 0 # de facto indistinguishable from 0 anyways
  if(any(!large)) {
    #lambda <- (sqrt(th^2+2*sig^2/kap)-th)/sig^2
    nu_bar_vargamma_minus <- function(x, z) {
      lambda <- lambda_plus * (x > 0) + lambda_minus * (x < 0)
      -c*expint_Ei(-lambda*abs(x), give=FALSE) - abs(z)
    }
    
    if (any(Gamma > 0)) {
      res[!large & Gamma > 0] <- vapply(Gamma[!large & Gamma > 0], function(Gamma.)
        uniroot(nu_bar_vargamma_minus, z=Gamma.,
                interval=c(.Machine$double.xmin, 29))$root, NA_real_)    
    }
    
    if (any(Gamma < 0)) {
      res[!large & Gamma < 0] <- vapply(Gamma[!large & Gamma < 0], function(Gamma.)
        uniroot(nu_bar_vargamma_minus, z=Gamma.,
                interval=c(-29, -1*.Machine$double.xmin))$root, NA_real_)
    }
  }
  res
}

Z1 <- matrix(nrow = n, ncol = d)

Z2 <- matrix(nrow = n, ncol = d)


## [END] COPIED (and modified) FROM https://cran.r-project.org/web/packages/copula/vignettes/NALC.html
for (r in 1:n) {
  # now let's simulate the two levy processes?
  # along the lines of Theorem 4.3 in the Tankov-paper
  
  # define X: jump times of Poisson process with lambda = 2 
  inter_arrival_times <- rexp(m, 2)
  X <- cumsum(inter_arrival_times)
  
  # define sequence {Gamma_i^1 : i \in 1, ..., m} as in Remark 4.3
  Gamma1 <- (-1)^(1:m) * X
  
  # define the sequence {Gamma_i^2 : i \in 1, ..., m} as on page 14
  W <- runif(m) # W_i \sim U[0,1]
  
  # THIS IS BAD STYLE AND WE WILL IMPROVE IT LATER!
  Gamma2 <- c()
  for (i in 1:m) {
    Gamma2[i] <- F_inv(Gamma1[i], W[i])
  }
  
  U1_inv <- function(Gamma) {
    res <- nu_bar_inv_vargamma(Gamma, c = c1, lambda_plus = lambda1_plus, lambda_minus = lambda1_minus)
  }
  U1_inv <- Vectorize(U1_inv)
  
  
  U2_inv <- function(Gamma) {
    res <- nu_bar_inv_vargamma(Gamma, c = c2, lambda_plus = lambda2_plus, lambda_minus = lambda2_minus)
  }
  U2_inv <- Vectorize(U2_inv)
  
  # The next lines are according to eq. (4.9) in the Tankov-paper
  
  # use a common sequence {V_i : i \in 1, ..., m} as in Theorem 4.3
  V <- runif(m)
  
  # the function aux computes (4.9) for one t in [0,1] 
  U_inv_of_Gamma1 <- U1_inv(Gamma1) # otherwise it will be computed length(grid_pts) times
  aux1 <- function(time_t) {
    sum(U_inv_of_Gamma1 * (V <= time_t))
  }
  aux1 <- Vectorize(aux1)
  Z1[r,] <- sapply(grid_pts, aux1)
  
  U_inv_of_Gamma2 <- U2_inv(Gamma2) # otherwise it will be computed length(grid_pts) times
  aux2 <- function(time_t) {
    sum(U_inv_of_Gamma2 * (V <= time_t))
  }
  aux2 <- Vectorize(aux2)
  Z2[r,] <- sapply(grid_pts, aux2)
}
plot(x = grid_pts, y = Z1[1,], type="l", ylim=c(min(Z1, Z2), max(Z1, Z2)), ylab="Z1 (black) and Z2 (red))", main = "10 Simulation of two VG processes")
lines(x = grid_pts, y = Z2[1,], col="red")

for (r in 2:n) {
  lines(x = grid_pts, y = Z1[r,], col="black")
  lines(x = grid_pts, y = Z2[r,], col="red")
}



###My contribution: Best copula###
#rank scatter plot
scatco = -coredata(B) #I change the sign of the return, I find it more convenient in computing the VAR
scatcoca = pobs(scatco)
scatam = -coredata(C)
scatamex = pobs(scatam)
plot(scatcoca,scatamex, pch=19, col='red', cex=0.5, xlab='KO', ylab='AXP', main = "Scatter plot of Coca-Cola and American Express")

M = cbind(scatco,scatamex)
corKendall(as.matrix(M)) #same as before even if I changed the sign

#add sigma to the t distribution
dt_ls = function(x, df=1, mu=0, sigma=1) 1/sigma * dt((x - mu)/sigma, df)
pt_ls = function(q, df=1, mu=0, sigma=1)  pt((q - mu)/sigma, df)
qt_ls = function(p, df=1, mu=0, sigma=1)  qt(p, df)*sigma + mu
rt_ls = function(n, df=1, mu=0, sigma=1)  rt(n,df)*sigma + mu

par(mfrow=c(2,1))

#selection of best marginals: https://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
#for KO
x = unmatrix(scatco, byrow = TRUE)
x_fit <= fitDist(x, k = 2, type = "realline", trace = TRUE, try.gamlss = TRUE)
summary(x_fit)
x_fit$fits #look for the smallest AIC
barplot(x_fit$fits, main = "AIC score for KO") #TF the best: t family distribution
x_mod.t = fitdistrplus::fitdist(x, 't_ls', start = list(df=1,mu=mean(x),sigma=sd(x)))
summary(x_mod.t)

#for AXP
y = unmatrix(scatam, byrow = TRUE)
y_fit = fitDist(y, k = 2, type = "realline", trace = TRUE, try.gamlss = TRUE)
summary(y_fit)
y_fit$fits #look for the smallest AIC
barplot(y_fit$fits, main = "AIC score for AXP") #Johnson's SU-distribution the best, but the difference with students't is very small, so I chose the t for simplicity
y_mod.t = fitdistrplus::fitdist(y, 't_ls', start = list(df=1,mu=mean(x),sigma=sd(x)))
summary(y_mod.t)

#params
x_parameters =list('df'=x_mod.t$estimate[[1]],
                    'mu'=x_mod.t$estimate[[2]],
                    'sigma'=x_mod.t$estimate[[3]])
x_type <-'t_ls'
y_parameters =list('df'=y_mod.t$estimate[[1]],
                    'mu'=y_mod.t$estimate[[2]],
                    'sigma'=y_mod.t$estimate[[3]])
y_type <-'t_ls'

par(mfrow=c(1,1))

#plots
#x
hist(x, breaks = 100, probability = T,
     col='grey', border = 'white', main = "Density comparison of KO")
points(density(x), 
       type='l', col='red', lwd=2)
points(sort(x), dt_ls(sort(x), 
                      x_parameters[[1]], 
                      x_parameters[[2]], 
                      x_parameters[[3]]),
       type='l', col='purple', lwd=2)

#y
hist(y, breaks = 100, probability = T,
     col='grey', border = 'white', main = "Density comparison of AXP")
points(density(y), 
       type='l', col='red', lwd=2)
points(sort(y), dt_ls(sort(y), 
                      y_parameters[[1]], 
                      y_parameters[[2]], 
                      y_parameters[[3]]),
       type='l', col='purple', lwd=2)

#best copula selection
x_probs =pt_ls(x, 
                x_parameters[[1]], 
                x_parameters[[2]], 
                x_parameters[[3]])

y_probs =pt_ls(y, 
                y_parameters[[1]], 
                y_parameters[[2]], 
                y_parameters[[3]])

selectedCopula =BiCopSelect(x_probs, y_probs, familyset = c(0,1,2,3,4,5,6), se=T) #in theory the best copula was: BB7 (par = 1.72, par2 = 0.59, tau = 0.4)
                                                                      #But for the sake of simplicity and beacuse I don't know this copula I restricted the searching fot he best copula to the common ones
selectedCopula #best copula: t (par = 0.52, par2 = 2.03, tau = 0.34)

#best copula
copula_fit <-fitCopula(tCopula(dim=2, df=2, df.fixed = T), cbind(x_probs, y_probs), method='itau')
copula_fit 
summary(copula_fit)
t_copula <-copula_fit@copula

#normal copula (will be useful later)
copula_fit <-fitCopula(normalCopula(dim=2), cbind(x_probs, y_probs), method='irho')
copula_fit
summary(copula_fit)
normal_copula <-copula_fit@copula

#independent copula
independent_copula <-indepCopula(dim=2)

#User defined function to count how many times when the simulated KO var at 99%, AXP exceeds the value var of KO
  MonteCarlo <-function(copula, M, verbose){
    print(copula)
    counting <-matrix(0, M, 1)
    for (i in 1:M) {
      
      #simulation
      biv_probs = rCopula(10000, copula)
      x = qt_ls(biv_probs[,1], x_parameters[[1]], x_parameters[[2]], x_parameters[[3]])
      y = qt_ls(biv_probs[,2], y_parameters[[1]], y_parameters[[2]], y_parameters[[3]])
      data_sim = cbind(x, y)
      
      #find the index of the observations when Y=VaR_99(Y) (Empirical)
      data_sim_ord = data_sim[order(data_sim[,1]), ]
      idx = nrow(data_sim)*0.99
      
      #compare X and Y (if X>=VaR_99(Y) count 1, otherwise count 0)
      counting[i,1] = data_sim_ord[idx,2]>=data_sim_ord[idx,1]
      
      if (verbose) {
        cat(paste(i, ''))
      }
    }
    return(sum(counting)/nrow(counting))
  }

set.seed(123)
M = 1000
verbose = T  

#normal
prob_normal_copula = MonteCarlo(normal_copula, M, verbose)
cat('\n\n')

#tcopula
prob_t_copula = MonteCarlo(t_copula, M, verbose)
cat('\n\n')

#independent copula
prob_independent = MonteCarlo(independent_copula, M, verbose)
cat('\n\n')

Prob <- list(
  't-copula'=prob_t_copula,
  'Gaussian-copula'=prob_normal_copula,
  'Independent-copula'=prob_independent)
Prob

#RSI indicator

# Create an RSI with a 3-day lookback period
spy_rsi <- RSI(price = Cl(KO), n = 3)

# Plot the closing price of KO
plot(Cl(KO))

# Plot the RSI 2
plot(RSI(Cl(KO), n = 2))
