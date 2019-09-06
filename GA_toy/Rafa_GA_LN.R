require(GA)

# Toy programm

# ======================================================================
#   COMPUTE PREDICTED VALUES: BLACK BOX
# ======================================================================
#     the numbers in the following equations are "hidden"
toy_model <- function(par){
  val1 <- 1e-1*par[1]^4 +
          2e-1*par[2]^3 +
          7e-1*par[3]^2 +
          4e-1*par[4];
  val2 <- 3e-2*par[1]^4 +
          5e-2*par[2]^3 +
          9e-2*par[3]^2 + 
          1e-2*par[4];
  val3 <- 7e-3*par[1]^4 +
          1e-3*par[2]^3 +
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

# Observed data
#x1 = 34
#x2 = 12
#x3 = 456
#x4 = 78

obs_data <- toy_model(c(34,12,456,78))

muL <- function(mu,scat){log(mu/sqrt(1+(scat^2)/mu^2))}
fit <- function(x, mu,scat) sum((log(x) - muL(mu,scat))^2)


loss <- function(xx,yy,zz,kk,scat){fit(toy_model(c(xx,yy,zz,kk)),
  obs_data,scat)}

GA <- ga(type = "real-valued", 
         fitness =  function(x) -loss(x[1],x[2],x[3],x[4],x[5]),
         lower = c(0.00, 0.00, 0.00, 0.00, 0.00), 
         upper = c(1e3, 1e3, 1e3, 1e3,1),
         popSize = 100, 
         maxiter = 5000)

print(summary(GA))

pdf(file="Rafa_GA.pdf", width=8, height=5, onefile=F)
plot(GA)
dev.off( )
