library(coda)
library(MASS)
library(BayesianTools)
library(latex2exp)
library(tidyverse)
x= c(0.070, 0.090, 0.120, 0.170, 0.179, 0.199, 0.200, 0.270, 0.280, 0.350, 0.352, 0.400, 0.440, 0.480, 0.593, 0.600, 0.680, 0.730, 0.781, 0.875, 0.880, 0.900, 1.037, 1.300, 1.430, 1.530, 1.750, 2.300)
y=c(69.0,  69.0,  68.6,  83.0,  75.0,  75.0,  72.9,  77.0,  88.8,  76.3,  83.0,  95.0,  82.6,  97.0, 104.0,  87.9,  92.0,  97.3, 105.0, 125.0,  90.0, 117.0,154.0, 168.0, 177.0, 140.0, 202.0, 224.0)
sd=c(19.6, 12.0, 26.2,  8.0,  4.0,  5.0, 29.6, 14.0, 36.6,  5.6, 14.0, 17.0,  7.8, 62.0, 13.0,  6.1,  8.0,  7.0, 12.0, 17.0, 40.0, 23.0, 20.0, 17.0, 18.0, 14.0, 40.0,8.0)
likelihood <- function(param)                                         
{
  a = param[1]
  b = param[2]
  pred=((sqrt((((1+x)^3)*a)+(1-a)))*b)
  chi=(sum((y-pred)*(y-pred)/(sd^2)))/26
  singlelikelihoods=-(((y-pred)*(y-pred)/(2*sd^2))+0.5*log(2*pi*(sd^2)))                    
  sumll = sum(singlelikelihoods)
  return(sumll)    
}
prior <- function(param)                                               
{                                            
  a = param[1]
  b = param[2]
  aprior = dunif(a, min=0, max=1, log = T)
  bprior = dunif(b, min=0, max=90, log = T)
  return(aprior+bprior)
}
proposalfunction <- function(param)                                  
{ 
  return(rnorm(2,mean = param, sd= c(0.05,0.5)))
}
run_metropolis_MCMC <- function(startvalue, iterations)                    
{            
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations)
  {
    proposal = proposalfunction(chain[i,])
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab)
    {
      chain[i+1,] = proposal
    }
    else
    {
      chain[i+1,] = chain[i,]
    }
  }
  return(mcmc(chain))
}
n = 100000
startvalue = c(0.2,70)
chain = run_metropolis_MCMC(startvalue, n)
burnin=0.2*n 
chain_burn=tail(chain,n-burnin) 
print(summary(chain))
samples_mat <- as.matrix(chain_burn)
plot(samples_mat,xlim = c(0.2,0.4),ylim = c(60,80))
y=samples_mat[,1]
x=samples_mat[,2]
den<-kde2d(x, y,n=200,lims=c(60,75,0.2,0.4))
den.z <-den$z
z <- array()
for (i in 1:n){
  z.x <- max(which(den$x < x[i]))
  z.y <- max(which(den$y < y[i]))
  z[i] <- den$z[z.x, z.y]
}
confidence.border <- quantile(z, probs=1-0.6827, na.rm = TRUE) # +-1sd
confidence.border1 <- quantile(z, probs=1-0.9516, na.rm = TRUE) # +-1sd
plot(x,y,xlab=TeX('$H_0$'),ylab = TeX('$\\Omega_{m0}$'))
par(new=TRUE)
contour(den, levels=confidence.border, col = "red", add = TRUE)
contour(den, levels=confidence.border1, col = "blue", add = TRUE)
ggs_pairs(ggs(chain_burn),mapping = aes(color = Chain))
