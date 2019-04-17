require(EasyABC)
require(ggplot2)
# Toy programm


# ======================================================================
#   COMPUTE PREDICTED VALUES: BLACK BOX
# ======================================================================
#     the numbers in the following equations are "hidden"
toy_model <- function(par){
  val1 <- 1e-1*par[1]^4 +
          2e-1*par[2]^3.0 +
          7e-1*par[3]^2 +
          4e-1*par[4];
  val2 <- 3e-2*par[1]^4.0 +
          5e-2*par[2]^3.0 +
          9e-2*par[3]^2.0 + 
          1e-2*par[4];
  val3 <- 7e-3*par[1]^4.0 +
          1e-3*par[2]^3.0 +
          4e-3*par[3]^2 +
          2e-3*par[4];
  val4 <- 2e-4*par[1]^4 +
          5e-4*par[2]^3 +
          8e-4*par[3]^2 +
          8e-4*par[4]
#S1 <- mean(c(val1,val2,val3,val4)) 
#S2 <- sd(c(val1,val2,val3,val4))
  return(c(log(val1),log(val2),log(val3),log(val4)))
}

# Observed data
x1 = 34
x2 = 12
x3 = 456
x4 = 78

obs_data <- toy_model(c(x1,x2,x3,x4))
# Summary, just the data - no compression
sum_stat_obs <- obs_data

toy_prior <- list(c("unif",1e-1,1e3),c("unif",1e-1,1e3),
                  c("unif",1e-1,1e3),c("unif",1e-1,1e3))


# Fraction of acceptable data, measured via Euclidian distance between simulation and obs_data
tolerance=c(1,0.75,0.5,0.25,0.15)

# Run Sequencial Approximate Bayesian Computation
ABC_Beaumont <- ABC_sequential(method="Beaumont", model=toy_model,
                             prior=toy_prior, nb_simul=1e4, summary_stat_target=sum_stat_obs,
                             tolerance_tab=tolerance)

# Plot histogram
hist(ABC_Beaumont$param[,1])

# Plot density 2D
pdat <- as.data.frame(ABC_Beaumont$param)

ggplot(data = pdat,aes(x = V1,y = V2)) +
  geom_density_2d()
