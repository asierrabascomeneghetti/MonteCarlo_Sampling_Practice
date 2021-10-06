##########################################################
#Montecarlo error computations
##########################################################
#through the use of the CLT we can approximate the accuracy of the montecarlo estimates
set.seed(32)
m <- 100
a <- 2
b <- 1/3

#the sample mean of the thetas will follow a normal distribution with the variance being the true theta variance and the mean being the true mean
#divided bz the montecarlo size
theta <- rgamma(n=m,shape=a,rate=b)
#to create the standard error:
se <- sd(theta)/sqrt(m)
# 95% ce: the true value of the estimnate is with 95% cinfidence within -2 se + mean / +2 se +mean
n95_CI <- 2*se
mean(theta) + n95_CI
mean(theta) - n95_CI
ind <- theta < 5
#esti9mate of the probability of being less than 5
mean(ind)
#the estimate still follows the CLT so, the error is still calculated the same way
se_mean_ind <- sd(ind)/sqrt(m)

2*se_mean_ind
#the true value of our estimate is then around 49-51%
#hierarchical model
#We had a binomial random variable, where y_1 represents success or failure with 10 trials, phi_1: success probability
#phi_1 has a beta distribution
#phi_1 <- beta(2,2)
#y_1 <- binomial(10,phi_1)


m <- 1e5
#function for draw simulator, we apply the rbinomial and rbeta

phi <-rbeta(m,shape1=2,shape2=2)
y <- rbinom(m, size = 10, prob = phi) #this is a beta-binomial, it is the marginal distribution of y
plot(table(y))
#simulated y mean
mean(y)


