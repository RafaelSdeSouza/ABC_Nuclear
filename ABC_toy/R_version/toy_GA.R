require(GA)
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
  val5 <- 9e-5*par[1]^4 +
          6e-5*par[2]^3 +
          3e-5*par[3]^2 +
          1e-5*par[4]
  

return(c(val1,val2,val3,val4,val5))

}



toy_model <- cmpfun(toy_model)
# Observed data
x1 = 34
x2 = 12
x3 = 456
x4 = 78

obs_data <- toy_model(c(x1,x2,x3,x4))
# Summary, just the data - no compression
sum_stat_obs <- obs_data

euc.dist <- function(x1, x2) sqrt(sum((log(x1) - log(x2))^2))
mah.dist <- function(x1,x2){
x <- cbind(x1,x2)
mean <- colMeans(x)
Sx <- cov(x)
sqrt(sum(mahalanobis(x,mean,Sx)^2))
}

loss <- function(xx,yy,zz,kk){mah.dist(toy_model(c(xx,yy,zz,kk)),obs_data)}

GA <- ga(type = "real-valued", fitness =  function(x) -loss(x[1],x[2],x[3],x[4]),
         lower = c(0.001, 0.001,0.001,0.001), upper = c(1e3, 1e3,1e3,1e3),
         popSize = 1000, maxiter = 5000)


summary(GA)
plot(GA)
