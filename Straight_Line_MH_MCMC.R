#===============================================================================
#=============================>>Load Requried Library<<=========================
#===============================================================================

library(coda)
library(ggmcmc)

#===============================================================================
#==========================>>Generate Some Test datapoints<<====================
#===============================================================================

x= c(1,2,3,4,5)
y=c(5,7,9,11,13)

#===============================================================================
#=================================>>Define Likelihood<<=========================
#===============================================================================

likelihood <- function(param) 
{
  m = param[1]
  c = param[2]
  y_model=m*x+c
  singlelikelihoods=dnorm(y_model,y,0.1,log = T)
  sumll =sum(singlelikelihoods)
  return(sumll) 
}

#===============================================================================
#===================================>>Define Priors<<===========================
#===============================================================================


prior <- function(param) 
{ 
  m = param[1]
  c = param[2]
  m_prior = dunif(m, min=0, max=4, log = T)
  c_prior = dunif(c, min=0, max=4, log = T)
  return(m_prior+c_prior)
}

#===============================================================================
#==============================>>Define Proposal Function<<=====================
#===============================================================================

proposalfunction <- function(param) 
{ 
  return(rnorm(2,mean = param, sd=c(0.01,0.01)))
}

#===============================================================================
#=================================>>Run Main MCMC Block<<=======================
#===============================================================================

run_metropolis_MCMC <- function(startvalue, iterations)
{
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations)
  {
    proposal = proposalfunction(chain[i,])
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- 
                   prior(chain[i,]))
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

n = 100000 #Total Number of Samples

startvalue =c(2,3)  #Initial Seeds
chain = run_metropolis_MCMC(startvalue, n)
burnin=0.1*n #Burnin steps=20% of Total Samples
chain_burn=tail(chain,n-burnin) #Remove burnin Phase

#===============================================================================
#==============================>>Plot and Print Chains<<========================
#===============================================================================
plot((chain_burn))
print(summary(chain_burn))
#===============================================================================
#==============================>>Draw Contour Plots<<=====================
#===============================================================================
ggs_pairs(ggs(chain_burn),mapping = aes(color = Chain))
