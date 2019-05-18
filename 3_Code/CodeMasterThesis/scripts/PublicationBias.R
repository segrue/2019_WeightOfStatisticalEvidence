### import packages -----
#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
require("data.table")

#for plotting
require("ggplot2")
require("gridExtra")

#custom functions
source("./functions/helper_functions.R")

### Define starting assumptions and conditions ----
# hypotheses:
# H0: p = mu0
# H1: p > mu0

seed <- 20190514
set.seed(seed)

n_sim <- 1e2 #number of times the setting below shall be simulated
n_studies <- c(5,10,20,30,40,50,100,200,500,1000) #number of subjects per study

mu0s <- 0
mu1s <- seq(-2,2,by=0.1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)

T_corr <- F #whether correction should be used or not when calculating vst
correction <- ifelse(T_corr,"Corr","MLE")

name <- paste0("DF_Ev_Mean_",correction,"_",n_sim,"_",seed) #name to load and save data etc.
fig_name <- paste0("Funnel_Ev_Mean_",correction,"_",n_sim,"_",seed)

### load required data sets ----
load(paste0("data/",name,".RData"),verbose=TRUE) #load data set containing simulated values
mu1_select <- 0.2
n_study_min <- 10
n_study_max <- 200
rng <- range(c(n_study_min,n_study_max))
ctgs <- unique(dat$n_study)[unique(dat$n_study) %between% rng]

### simulate publication bias
## Situtation 1A:
## no publication bias - all studies where n_study between 10 and 200 are included (n=100 per study size group)
dat <- evidence_df[mu1==mu1_select & n_study %between% rng,]

clt <- dat[id == "clt",]
vst <- dat[id == "vst",]
stud <- dat[id == "stud",]

funnel_plotter(clt,vst,vst,xlim=2,figname=paste0(fig_name,"_1A.pdf"),ctgs)

## Situtation 1B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1==1,]
vst_sig <- vst[H1==1,]
stud_sig <- stud[H1==1,]

funnel_plotter(clt_sig,vst_sig,stud_sig,xlim=2,figname=paste0(fig_name,"_1B.pdf"),ctgs)

## Situtation 2A:
## only keep 50 studies in total, but weigh according to study size; 
## choose studies with lower sample size with higher probability
n_studies_selected <- n_studies[n_studies %between% rng]
probs <- c(0.29,0.24,0.2,0.14,0.09,0.03,0.01)
n_tot <- 50

clt <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,seed,id="clt")
vst <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,seed,id="vst")
stud <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,seed,id="stud")

stopifnot(sum(clt$mu1_hat-vst$mu1_hat)==0,sum(stud$mu1_hat-vst$mu1_hat)==0)

funnel_plotter(clt,vst,vst,xlim=2,figname=paste0(fig_name,"_2A.pdf"),ctgs)

## Situtation 2B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1==1,]
vst_sig <- vst[H1==1,]
stud_sig <- stud[H1==1,]

funnel_plotter(clt_sig,vst_sig,stud_sig,xlim=2,figname=paste0(fig_name,"_2B.pdf"),ctgs)

## Situtation 2C:
## only keep studies which turned out to be significant and small percentage of non-significant studies
n_nsig <- 5
probs_mix <- c(0,0,0,0,0.2,0.3,0.5)
clt_mix <- rbind(clt_sig,select_studies(clt[H1==0,],rev(probs),n_studies_selected,n_total=n_nsig,seed,id="clt"))
vst_mix <- rbind(vst_sig,select_studies(clt[H1==0,],rev(probs),n_studies_selected,n_total=n_nsig,seed,id="vst"))
stud_mix <- rbind(stud_sig,select_studies(clt[H1==0,],rev(probs),n_studies_selected,n_total=n_nsig,seed,id="stud"))

funnel_plotter(clt_mix,vst_mix,stud_mix,xlim=2,figname=paste0(fig_name,"_2C.pdf"),ctgs)

#function to calculate the aggregated mean of the studies
aggregate_mean <- function(dat){
  agg_mean <- sum(dat$mu1_hat*dat$n_study)/sum(dat$n_study)
  return(agg_mean)
}
aggregate_mean(stud)
