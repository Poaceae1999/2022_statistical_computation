# gogogo
# parameters
# hist ----
b_tau <- 1# i don't know what actually is,for cauchy plus(0,b_tau)
a <- 1 # for proposal distribution inverse gamma(a,b) 
b <- 1 # for proposal distribution inverse gamma(a,b)

Cmax <- Cauchy_plus(10^-8,b_tau)
target_dist <- function(x,a_1,b_1,b_cauchy){
  fx <- dinvgamma(x^2,shape =a_1,scale=b_1)*Cauchy_plus(x,b_cauchy=b_cauchy)*7.3
  return(fx)
}
dat <- collected_data_from(accept_reject_K_invG,times=10000,a_1=a,b_1=b,a = a,b = b,K = Cmax,fx=target_dist,b_cauchy=b_tau)

draw(target_dist,dat,b_cauchy=b_tau,a_1=a,b_1=b)

#QQplot----
percentile <- seq(from = 10^-8,to = 50,by = 0.1)
pr <- numeric(500L)
QR <- numeric(20L)
for (i in 1:500) {
  pr[i] <-  integrate(target_dist,lower=0,upper=percentile[i],b_1=b,a_1=a,b_cauchy=b_tau)#quantiles of Rayleigh
}
pr <- pr %>% unlist()
for (i in 1:20) {
  QR[i] <-  target_dist(percentile[6:25][i],b_1=b,a_1=a,b_cauchy=b_tau)#quantiles of Rayleigh
} 
Q <- quantile(dat, pr[6:25])

qqplot(QR, Q, main="",
           xlab="Target Quantiles", ylab="Sample Quantiles")

# bayes----
b <- .2 #actual value of beta
w <- .25 #width of the uniform support set
m <- 5000 #length of the chain
burn <- 1000 #burn-in time
days <- 250
x <- numeric(m) #the chain
# generate the observed frequencies of winners
   i <- sample(1:5, size=days, replace=TRUE,
                prob=c(1/3, (1-b)/3, (1-2*b)/3, (2*b)/3, b/3))
win <- tabulate(i)
print(win)


prob <- function(y, win) {
  # computes (without the constant) the target density
    if (y < 0 || y >= 0.5)
       return (0)
   return(((1/3)^win[1] *
              ((1-y)/3)^win[2] * ((1-2*y)/3)^win[3] *
             (2*y)/3)^win[4] * (y/3)^win[5])
}

u <- runif(m) #for accept/reject step
 v <- runif(m, -w, w) #proposal distribution
x[1] <- .25
for (i in 2:m) {
   y <- x[i-1] + v[i]
   if (u[i] <= prob(y, win) / prob(x[i-1], win))
     x[i] <- y else
       x[i] <- x[i-1]
}
print(win)
print(round(win/days, 3))
print(round(c(1, 1-b, 1-2*b, 2*b, b)/3, 3))
xb <- x[(burn+1):m]
print(mean(xb))

