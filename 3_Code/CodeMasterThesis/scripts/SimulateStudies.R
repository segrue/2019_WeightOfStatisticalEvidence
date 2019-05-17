### import packages ----
#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
require("data.table")

#custom functions
source("./functions/helper_functions.R")

### Define starting assumptions and conditions ----
# hypotheses:
# H0: p = mu0
# H1: p > mu0
n_sim <- 1e3 #number of times the setting below shall be simulated
seed <- 20190514
set.seed(seed)

n_studies <- c(5,10,20,30,40,50,100,200,500,1000) #number of subjects per study

mu0s <- 0
mu1s <- round(seq(-2,2,by=0.1),1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)

T_corr <- F #whether correction should be used or not when calculating vst

### simulate individual studies and calculate according summary and evidence statistics ----
df <- evidence_in_mean(mu0_vec=mu0s,mu1s,sgm0,alphas,T_corr=F,seed,n_studies,n_sim)
