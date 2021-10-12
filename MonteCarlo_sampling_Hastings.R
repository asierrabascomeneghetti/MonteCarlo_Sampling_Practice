#toy example for the Metropolis-Hastings Algorithm
#data: percentage change in total personel from one year to another for 10 companies
# yi| nu ~iid~ N(nu,1) : change of personel for company i given nu. The prior distribution on nu
#nu ~ t(0,1,1) t distribution (Cauchy distribution). The distribution is not conjugate so 
#the posterior does not have a closed form

#to get posteriors, we will need a markov chain whose stationary distribution is the distribution we want

#p(mu| y1,...yn) ~ exp(n (<y>mu - mu**2)/2)/(1+mu**2) which is a normal dist, the dist is symmetric
#therefore p(mu i+1 |mui)= p(mu i |mu i+1):
#log (g(nu)) = n((<y>mu - mu**2)/2)) - log (1+mu**2), this will be our g
library("coda")
lg <- function(ybar,n,mu){
  mu2 <- mu**2
  n*(ybar* mu - mu2/2) - log(1.0 + mu2)
}

mh <- function(n, ybar, n_iter,mu_init,cand_sd) {
  mu_out <- numeric(n_iter)
  accept <- 0
  mu_now <- mu_init
  lg_now <- lg(ybar=ybar,n=n,mu=mu_now)
  
  for (i in 1:n_iter) {
    mu_cand <- rnorm(1,mean=mu_now, sd=cand_sd)
    lg_cand <- lg(ybar=ybar, n=n,mu=mu_cand)
    alpha <- exp(lg_cand -lg_now)
    u <- runif(1)
    if (u< alpha) {
      mu_now <-mu_cand
      accept <- accept + 1
      lg_now <-lg_cand
    }
    mu_out[i] <- mu_now
  }
  list(mu=mu_out, accept=accept/n_iter)
}
#########################################
y <-c(1.2,1.4,-0.5,0.3,0.9,2.3,1.0,0.1,1.2,1.9)
ybar <- mean(y)
n <- length(y)
hist(y, freq = FALSE, xlim = c(-1.0,3.0))
points(y,rep(0.0,n))
points(ybar,0,pch=19)

#print the prior to our data
curve(dt(x, df=1), lty=2, add=TRUE) #discrepancy between prior mu and data mu, we expect the posterior
#to be a compromise between the sample data and the mu we already have

#posterior sampling
set.seed(43)
post <- mh(n=n,ybar=ybar, n_iter=1e3, mu_init=0, cand_sd = 3.0)
#acceptance rate of abt 10%
traceplot(as.mcmc(post$mu))
#shows the history of the markov chain, the acceptance range is low, so the acceptance has to be bigger
#therefore, we will use smaller sd
post <- mh(n=n,ybar=ybar, n_iter=1e3, mu_init=0, cand_sd = 0.05)
traceplot(as.mcmc(post$mu))
#the steps are too small;
post <- mh(n=n,ybar=ybar, n_iter=1e3, mu_init=0, cand_sd = 0.9)
traceplot(as.mcmc(post$mu))
#what happens if we use a different mu value?
post <- mh(n=n,ybar=ybar, n_iter=1e3, mu_init=20, cand_sd = 0.9)
str(post)
traceplot(as.mcmc(post$mu))
#we will remove exploratory values of the montecarlo algorithm, when the sampling is not yet stable
post$mu_keep <- post$mu[-c(1:100)]
#we need to compare the prior and posterior
plot(density(post$mu_keep), xlim= c(-1.0,3.0))
curve(dt(x, df=1), lty=2, add=TRUE)
points(ybar,0,pch=19)

####################################################################
#It still requires some fine-analysis to see if the sampled posterior is the converged MCMC 

