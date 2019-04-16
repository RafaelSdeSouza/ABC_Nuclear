require(EasyABC)

# Toy programm


# ======================================================================
#   COMPUTE PREDICTED VALUES: BLACK BOX
# ======================================================================
#     the numbers in the following equations are "hidden"



model <- function(par){
  val1 <- 1e-1*par[1] +
          2e-1*par[2]^2.0 +
          7e-1*par[3] +
          4e-1*par[4];
  val2 <- 3e-2*par[1]^4.0 +
          5e-2*par[2]^3.0 +
          9e-2*par[3]^2.0 + 
          1e-2*par[4];
  val3 <- 7e-3*par[1]^4.0 +
          1e-3*par[2]^3.0 +
          4e-3*par[3]^2.0 +
          2e-3*par[4];
  val4 <- 2e-4*par[1]^4.0 +
          5e-4*par[2]^3.0 +
          8e-4*par[3]^2.0 +
          8e-4*par[4];
  val5 <- 9e-5*par[1]^4.0 +
          6e-5*par[2]^3.0 +
          3e-4*par[3]^2.0 +
          1e-4*par[4]
  
return(c(val1,val2,val3,val4,val5))  
}

data <- model(c(34,12,456,78))
summarydata <- data

toy_prior <- list(c("unif",1e-1,1000),c("unif",1e-1,1000),
                  c("unif",1e-1,1000),c("unif",1e-1,1000))

ABC_Marjoram_original <- ABC_mcmc(method="Marjoram", model=model, 
                                prior=list(c("unif",1e-1,1000),c("unif",1e-1,1000),
                                           c("unif",1e-1,1000),c("unif",1e-1,1000)), 
                                summary_stat_target=summarydata, n_rec = 1e5)


ABC_Drovandi <-ABC_sequential(method="Drovandi", model=model, prior=toy_prior,
                             nb_simul=5000, summary_stat_target=summarydata, tolerance_tab=0.5, c=0.7)

ABC_Beaumont <- ABC_sequential(method="Beaumont", model=model, prior=list(c("unif",1e-1,1000),c("unif",1e-1,1000),
                                                                        c("unif",1e-1,1000),c("unif",1e-1,1000)),
                             nb_simul=1e4, summary_stat_target=summarydata, tolerance_tab=c(1.5,0.5))
                             
pdat <- as.data.frame(ABC_Beaumont$param)

ggplot(data = pdat,aes(x = V1,y = V4)) +
  geom_density_2d()
