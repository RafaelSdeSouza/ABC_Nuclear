library(coda)

# assuming the data are 10 samples of a normal distribution
# with mean 5.3 and sd 2.7
data =  rnorm(10, mean =5.3, sd = 2.7)


# we want to use ABC to infer the parameters that were used. 
# we sample from the same model and use mean and variance
# as summary statstitics. We return true for ABC acceptance when
# the difference to the data is smaller than a certain threshold

meandata <- mean(data)
standarddeviationdata <- sd(data)

ABC_acceptance <- function(par){
  
  # prior to avoid negative standard deviation
  if (par[2] <= 0) return(F) 
  
  # stochastic model generates a sample for given par
  samples <- rnorm(10, mean = par[1], sd = par[2])
  
  # comparison with the observed summary statistics
  diffmean <- abs(mean(samples) - meandata)
  diffsd <- abs(sd(samples) - standarddeviationdata)
  if((diffmean < 0.1) & (diffsd < 0.2)) return(T) else return(F)
}


# we plug this in in a standard metropolis hastings MCMC, 
# with the metropolis acceptance exchanged for the ABC acceptance

run_MCMC_ABC <- function(startvalue, iterations){
  
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  
  for (i in 1:iterations){
    
    # proposalfunction
    proposal = rnorm(2,mean = chain[i,], sd= c(0.7,0.7))
    
    if(ABC_acceptance(proposal)){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(mcmc(chain))
}

posterior <- run_MCMC_ABC(c(4,2.3),300000)
plot(posterior)