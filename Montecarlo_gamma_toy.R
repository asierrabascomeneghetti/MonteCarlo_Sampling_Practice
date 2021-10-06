##########################################
#Montecarlo simulation examples
##########################################

set.seed(32)
#Sampling from a gamma distribution to use as posterior for sampling
m <- 10000
a <- 2.0
b <- 1.0/3.0
#100 draws from a gamma distribution
theta <- rgamma(n=m,shape=a,rate=b)
#compare plots and draws
hist(theta, freq=FALSE)
curve(dgamma(x,shape=a,rate=b), col="blue",add=TRUE)
#we can use the simulated values for a MC through the mean estimator
theta_av <- mean(theta)      
theta_av_theory <- a/b
#error of 10% in approximation, we can make a function of different ms
m_numbers  <- c(100,300,600,900,1200,1500,1600,2000,5000,10000,50000)

theta_avg <- function(X,shape,rate){
  mean(rgamma(n=X, shape=shape, rate=rate))
}
thet_avg <- sapply(m_numbers, theta_avg, shape=a, rate=b) #sapply is faster than doing for loops within R
errors <- thet_avg - (a/b)
#Other characteristics can be approximated, for example probability that theta is less than 5:
#indicator variable (boolean):
ind <- theta < 5
#an approximation to the probability of theta being less than 5:
mean(ind)
#with pgamma we can get the CDF to compare (0.497 with m=10000 to true CDF value 0.4963)
pgamma(q=5.0, shape=a, rate=b)

#the quantile functions are
qgamma(shape=a, rate=b, p=0.9)
quantile(theta, probs=0.9)
#For other functions which might not be in the available so MC might be useful then