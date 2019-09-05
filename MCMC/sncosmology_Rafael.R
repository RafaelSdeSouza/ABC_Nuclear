# From Bayesian Models for Astrophysical Data 
# by Hilbe, de Souza & Ishida, 2016, Cambridge Univ. Press
#
# Code 10.26 Bayesian normal model for cosmological parameter 
#            inference from type Ia supernova data in R using Stan.
#
# Statistical Model: Gaussian regression in R using Stan
#                    example using ODE
#
# Astronomy case: Cosmological parameters inference from 
#                 type Ia supernovae data 
#
# Data: JLA sample, Betoule et al., 2014  
# http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html
#
# 1 response (obsy - observed magnitude)
# 5 explanatory variable (redshift - redshift,
#                         ObsMag   - apparent magnitude,
#                         x1       - stretch,
#                         color    - color,
#                         hmass    - host mass)


library(rstan)
library(mvtnorm)
library(SDMTools)
library(emdbook)
library(ggthemes)
library(coda)
library(tidyr)
library(dplyr)
# Preparation

# set initial conditions
z0 = 0                          # initial redshift
E0 = 0                          # integral(1/E) at z0

# physical constants
c = 3e5                         # speed of light
H0 = 73.24                         # Hubble constant

# Data
dataTrain <- read.csv("..//data/sncosmo_fit/sncosmo_fit_Ia_train.csv",
                  header = T, sep = "") # Emille


dataTest <- read.csv("..//data/sncosmo_fit/sncosmo_fit_Ia_test.csv",
                  header = T) # Emille

SnanaTest <- read.table("..//data/sncosmo_fit/plasticc_ia.FITRES.TEXT",
                     header = T,skip=4,sep="") # Emille


#data0 <- read.table("..//data/salt2_fit/sncosmo_fit_Ia_test_ 9000.csv",header=T) # Emille

#data2 <- read.table("/home/emille/Documents/git/BMAD/data/Section_10p11/jla_lcparams.txt",header=T)
data1 <- SnanaTest[!duplicated(SnanaTest$zCMB),] %>%  as_tibble() %>%
  filter(-0.3 < c & c < 0.3) %>%
  filter(-3 < x1 & x1 < 3) %>%
  filter(x1ERR < 1) %>%
  filter(FITPROB >= 0.995) 
#%>%
#  filter(zCMB > 0.051)



RedshiftAll <- data1$zCMB
ObsMagAll <- data1$mB
ColorAll <- data1$c
x1All <- data1$x1
x0All <- data1$x0



# prepare data for Stan
nobs         <- length(RedshiftAll)        # number of SNe
index        <- order(RedshiftAll)         # sort according to redshift
ObsMag       <- ObsMagAll[index]           # apparent magnitude
redshift     <- RedshiftAll[index]         # redshift
color        <- ColorAll[index]            # color
x1           <- x1All[index]               # stretch
x0           <- x0All[index]
obs_mag2 = -2.5*log(x0,10) + 10.635

plot(redshift,obs_mag2)

stan_data  <- list(nobs = nobs,
                   E0 = array(E0,dim=1),
                   z0 = z0,
                   c = c,
                   H0 = 70,
                   obs_mag = ObsMag,    
                   redshift = redshift, 
                   x1 = x1, 
                   color = color,
                   x0 =  x0,
                   obs_mag2 = obs_mag2)

# Fit
stan_model = "
functions {
     /** 
     * ODE for the inverse Hubble parameter. 
     * System State E is 1 dimensional.  
     * The system has 2 parameters theta = (om, w)
     * 
     * where 
     * 
     *   om:       dark matter energy density 
     *   w:        dark energy equation of state parameter
     *
     * The system redshift derivative is 
     * 
     * d.E[1] / d.z  =  
     *  1.0/sqrt(om * pow(1+z,3) + (1-om) * (1+z)^(3 * (1+w)))
     * 
     * @param z redshift at which derivatives are evaluated. 
     * @param E system state at which derivatives are evaluated. 
     * @param params parameters for system. 
     * @param x_r real constants for system (empty). 
     * @param x_i integer constants for system (empty). 
     */ 
     real[] Ez(real z,
               real[] H,
               real[] params,
               real[] x_r,
               int[] x_i) {
           real dEdz[1];

           dEdz[1] = 1.0/sqrt(params[1]*(1+z)^3
                     +(1-params[1])*(1+z)^(3*(1+params[2])));

           return dEdz;
    } 
}
data {
    int<lower=1> nobs;              // number of data points
    real E0[1];                     // integral(1/H) at z=0                           
    real z0;                        // initial redshift, 0
    real c;                         // speed of light
    vector[nobs] obs_mag;           // observed magnitude at B max
    real x1[nobs];                  // stretch
    real x0[nobs];                  // x0 has no name
    real color[nobs];               // color 
    real redshift[nobs];            // redshift
    real H0;                       // hubble parameter
    real obs_mag2[nobs]; 
    
}
transformed data {
      real x_r[0];                  // required by ODE (empty)
      int x_i[0]; 
}

parameters{
      real<lower=0, upper=1> om;    // dark matter energy density
      real alpha;                   // stretch coefficient   
      real beta;                    // color coefficient
      real Mint;                    // intrinsic magnitude
      real<lower=0> sigint;         // magnitude dispersion
      real<lower=-2, upper=0> w;    // dark matter equation of state parameter
}
transformed parameters{
      real DC[nobs,1];                        // co-moving distance 
      real pars[2];                           // ODE input = (om, w)
      vector[nobs] mag;                       // apparent magnitude
      real dl[nobs];                          // luminosity distance
      real DH;                                // Hubble distance = c/H0
 
      DH = (c/H0);

      pars[1] = om;
      pars[2] = w;

      // Integral of 1/E(z) 
      DC = integrate_ode_rk45(Ez, E0, z0, redshift, pars,  x_r, x_i);

      for (i in 1:nobs) {
            dl[i] = DH * (1 + redshift[i]) * DC[i, 1];
            mag[i] = 25 + 5 * log10(dl[i]) + Mint - alpha * x1[i] + beta * color[i];
      }
}
model {

      // priors and likelihood
      sigint ~ gamma(0.01, 0.01);
      Mint ~ normal(-20, 5.);
      beta ~ normal(0, 10);
      alpha ~ normal(0, 1);
      om ~  uniform(0,1);
      w ~ uniform(-2,-0.33);
      obs_mag2 ~ normal(mag, sigint);   
}
"

# run MCMC
fit <- stan(model_code = stan_model,
                data = stan_data,
                seed =  10,
                chains = 3,
                iter = 5000,
                cores = 3,
                warmup = 2000
)

# Output 
print(fit,pars=c("om", "Mint","w","alpha","beta","sigint"),
      intervals = c(0.025, 0.975), digits = 3)


# plot

gdata <- as.data.frame(as.matrix(fit, pars = c("om", "w")))

ggplot(gdata, aes(x =  om, y = w)) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom="polygon",n = 500) +
  scale_fill_viridis_c() +
  theme_pander() + ylab("w") + xlab(expression(Omega[m])) +
#  coord_cartesian(ylim=c(-1.55,-0.4),xlim=c(0,0.45)) +
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y = element_text(vjust = 0.75),
        axis.title.x = element_text(vjust = -0.25),
        text = element_text(size = 25))




require(rayshader)
plot_gg(ggSN, width = 5, height = 5, multicore = TRUE, scale = 250,
        zoom = 0.7, theta = 10, phi = 30, windowsize = c(800, 800))

#gdata <- data.frame(w = w, om= omega_m)

#write.csv(gdata,"posterior_LCDM_train.csv",row.names=F)

#gdata <- read.csv("posterior_LCDM_train.csv",header=T)


#c68 <- as.data.frame(HPDregionplot(mcmc(data.matrix(gdata)), prob = 0.50, n = 300))
#c95 <- as.data.frame(HPDregionplot(mcmc(data.matrix(gdata)), prob = 0.90,n = 300))








#library(RColorBrewer)
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#jet.colors <-
#  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
#                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))