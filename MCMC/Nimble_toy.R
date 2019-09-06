######################################################################
# preparation: remove all variables from the work space
rm(list=ls())
# import jags package
library(nimble)

toy_temp <- function(x1,x2,x3,x4){
  val1 <- 1e-1*x1^4 +
    2e-1*x2^3 +
    7e-1*x3^2 +
    4e-1*x4;
  val2 <- 3e-2*x1^4 +
    5e-2*x2^3 +
    9e-2*x3^2 + 
    1e-2*x4;
  val3 <- 7e-3*x1^4 +
    1e-3*x2^3 +
    4e-3*x3^2 +
    2e-3*x4;
  val4 <- 2e-4*x1^4 +
    5e-4*x2^3 +
    8e-4*x3^2 +
    8e-4*x4;
  val5 <- 9e-5*x1^4 +
    6e-5*x2^3 +
    3e-5*x3^2 +
    1e-5*x4
  return(c(val1,val2,val3,val4,val5))
}

# Vectorised sigma function
toy_model <- nimbleRcall(function(x1 = double(0),x2 = double(0),x3 = double(0),
                                  x4 = double(0)
){},
Rfun = "toy_temp", returnType = double(1))


######################################################################
## DATA SETS
######################################################################
# data input; the input is of the form: obsx, obsy, errobsy, where the
# latter is the individual statistical error of each datum [i];
# energy is in MeV, sqrt(Ecm)*sigma is in sqrt(MeV)b

# Parameters
obsy <- toy_model(34,12,456,78)
N <- 5
samplerCode <- nimbleCode({

  ###################
  # LIKELIHOODS
  ###################
 
     for (i in 1:N){
       obsy[i] ~ dlnorm(muL[i], sd = y.scat)
       muL[i] <- log(mu[i]/sqrt(1 + (y.scat^2)/mu[i]^2))
       mu[i] <- toy_model(x1,x2,x3,x4)[i] 
       }

  ###################
  # PRIORS
  ###################
  x1 ~ dunif(0,1e3)
  x2 ~ dunif(0,1e3)
  x3 ~ dunif(0,1e3)
  x4 ~ dunif(0,1e3)
  y.scat ~ dunif(0,10)
}
)
  
samplerData <- list(obsy = obsy)
samplerConst <- list(N = N)


samplerInits <- list(x1 = runif(1,0,500),x2 =runif(1,0,500),
                      x3=runif(1,0,500),x4=runif(1,0,500),
                      y.scat = 0.1)
ourmodel <- nimbleModel(code = samplerCode, constants = samplerConst,
                        data = samplerData, inits = samplerInits, check = FALSE
)
compileNimble(ourmodel)
# Always compile the model after you are done setting up with it

ourmodel$calculate()
# Calculate the log likelihood (logProb). If the value is not NA,
# we have successfully initialised the model (FULLY)
# One iteration: simulate -> new values -> calculate

conf <- configureMCMC(ourmodel,print=TRUE)
# print = TRUE tells you what kind of samplers are being used for each stochastic node

conf$addMonitors(c('x1','x2','x3','x4','y.scat'))
samplerMCMC <- buildMCMC(conf)
compiledMCMC <- compileNimble(samplerMCMC,project = ourmodel,showCompilerOutput = TRUE)

n.chains = 3
n.iter = 30000
n.burnin = 25000

system.time(
  mcmcChain <- runMCMC(compiledMCMC,niter = n.iter, nchains = n.chains, nburnin = n.burnin,
                       samplesAsCodaMCMC = TRUE)
)

summary(mcmcChain)


pdf("MCMC_toy.pdf")
plot(mcmcChain)
dev.off()
