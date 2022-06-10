# gogogo
# parameters

b_tau <- 1# i don't know what actually is,for cauchy plus(0,b_tau)
a <- 1 # for proposal distribution inverse gamma(a,b) 
b <- 1 # for proposal distribution inverse gamma(a,b)

Cmax <- Cauchy_plus(10^-8,b_tau)
target_dist <- function(x,a_1,b_1,b_cauchy){
  fx <- dinvgamma(x^2,shape =a_1,scale=b_1)*Cauchy_plus(x,b_cauchy=b_cauchy)*6.873381
  return(fx)
}
dat <- collected_data_from(accept_reject_K_invG,times=50000,a_1=a,b_1=b,a = a,b = b,K = Cmax,fx=target_dist,b_cauchy=b_tau)

draw(target_dist,dat,b_cauchy=b_tau,a_1=a,b_1=b)

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

