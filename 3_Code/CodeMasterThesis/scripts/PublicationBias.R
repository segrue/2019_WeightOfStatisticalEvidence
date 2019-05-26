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
x_lim <- 5
y_lim <- 8
rng <- range(c(n_study_min,n_study_max))
dat <- evidence_df[mu1==mu1_select & n_study %between% rng,]
ctgs <- unique(dat$n_study)[unique(dat$n_study) %between% rng]
sample_seed <- 10

### simulate publication bias
## Situtation 1A:
## no publication bias - all studies where n_study between 10 and 200 are included (n=100 per study size group)
label <- "1A"
n_tot <- length(ctgs)*n_sim

clt <- dat[id == "clt",]
vst <- dat[id == "vst",]
stud <- dat[id == "stud",]

funnel_plotter(clt,vst,vst,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

## Situtation 1B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1==1,]
vst_sig <- vst[H1==1,]
stud_sig <- stud[H1==1,]
label <- "1B"

funnel_plotter(clt_sig,vst_sig,stud_sig,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

## Situtation 2A:
## only keep 100 studies in total, but weigh according to study size; 
## choose studies with lower sample size with higher probability
n_studies_selected <- n_studies[n_studies %between% rng]
probs <- c(0.29,0.24,0.2,0.14,0.09,0.03,0.01)
n_tot <- 100
label <- "2A"

clt <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,sample_seed,T_id="clt")
vst <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,sample_seed,T_id="vst")
stud <- select_studies(dat,probs,n_studies_selected,n_total=n_tot,sample_seed,T_id="stud")

stopifnot(sum(clt$mu1_hat-vst$mu1_hat)==0,sum(stud$mu1_hat-vst$mu1_hat)==0)

funnel_plotter(clt,vst,vst,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

## Situtation 2B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1==1,]
vst_sig <- vst[H1==1,]
stud_sig <- stud[H1==1,]
label <- "2B"

funnel_plotter(clt_sig,vst_sig,stud_sig,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

## Situtation 2C:
## only keep studies which turned out to be significant and small percentage of non-significant studies
sel_prob <- 0.1 #percentage of significant studies to keep
probs_mix <- c(0,0,0.1,0.1,0.1,0.1,0.1)
clt_mix <- rbind(clt_sig,select_studies(clt[H1==0,],sel_prob,n_studies_selected,n_total=n_nsig,sample_seed,T_id="clt"))
vst_mix <- rbind(vst_sig,select_studies(vst[H1==0,],sel_prob,n_studies_selected,n_total=n_nsig,sample_seed,T_id="vst"))
stud_mix <- rbind(stud_sig,select_studies(stud[H1==0,],sel_prob,n_studies_selected,n_total=n_nsig,sample_seed,T_id="stud"))
label <- "2C"

funnel_plotter(clt_mix,vst_mix,stud_mix,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

#function to calculate the aggregated mean of the studies
aggregate_mean <- function(dat){
  agg_mean <- sum(dat$mu1_hat*dat$n_study)/sum(dat$n_study)
  return(agg_mean)
}
aggregate_mean(stud)

### Calculate number of papers stored in file drawer based Rosenthal 1984
calc_filedrawer <- function(Z_list,alph){
  min_filedrawer <- c()
  for (Zs in Z_list){
    k <- length(Zs)
    z <- qnorm(1-alph)
    X <- (k*mean(Zs)/z)^2-k
    min_filedrawer <- c(min_filedrawer,X)
  }
  return(min_filedrawer)
}

#problem: Filedrawer problem only works for checking whether there is no null effect; it doesn't work in situations in which there really is an effects
Z_list <- list(clt$Tn,vst$Tn,qnorm(stud$p,0,1,lower.tail=FALSE),clt_sig$Tn,vst_sig$Tn,qnorm(stud_sig$p,0,1,lower.tail=FALSE))
X_list <- calc_filedrawer(Z_list,alphas)


### Calculate publication probability based on empirical distribution function
#Hansen-Hurwitz estimator https://newonlinecourses.science.psu.edu/stat506/node/15/
sel_prob <- 0.1
quant_probs <- seq(0,1,by=sel_prob)
clt_mix_quants<- quantile(clt_mix$mu1_hat,quant_probs)
hans_hurwitz_estimator <- function(n,sel_prob,quants){
  #selection <- runif(n,min=0,max=1)
  samples <- runif(n,min(quants),max(head(quants,-2)))
  idx <- sapply(samples, function(s) min(which((s<tail(quants,-1))==TRUE)))
  diffs <- (tail(quants,-1)-head(quants,-1))
  probs <- diffs/sum(diffs)
  prob_inclusion <- 1-(1-probs)^n
  prob_inclusion <- sapply(idx,function(i) prob_inclusion[i])
  #med <- head(quants,-1)+(tail(quants,-1)-head(quants,-1))/2
  hans_hur <- mean(samples/prob_inclusion)
  
}

horvitz_thompson_estimator <- function(dat){
  y_i <- dat$mu1_hat*dat$n_study
  n <- length(y)
  probs <- 0.1+0.9*dat$H1
  pi_i <- 1-(1-probs)^n
  print(pi_i)
  horv_thomp <-sum(y_i/pi_i)/sum(1/pi_i*dat$n_study)
}
horv_thomp <- horvitz_thompson_estimator(clt_mix)
horv_thomp

hans_hurwitz_estimator <- function(n,sel_prob,quants){
  #selection <- runif(n,min=0,max=1)
  samples <- runif(n,min(quants),max(head(quants,-2)))
  idx <- sapply(samples, function(s) min(which((s<tail(quants,-1))==TRUE)))
  diffs <- (tail(quants,-1)-head(quants,-1))
  probs <- diffs/sum(diffs)
  prob_inclusion <- 1-(1-probs)^n
  prob_inclusion <- sapply(idx,function(i) prob_inclusion[i])
  #med <- head(quants,-1)+(tail(quants,-1)-head(quants,-1))/2
  hans_hur <- mean(samples/prob_inclusion)
  
}

clt <- select_studies(dat,0.2,n_studies_selected,n_total=n_tot,sample_seed,T_id="clt")
clt_sig <- clt[H1==1,]
clt_mix <- rbind(clt_sig,select_studies(clt[H1==0,],sel_prob,n_studies_selected,n_total=n_nsig,sample_seed,T_id="clt"))

hans_hurwitz_estimator <- function(dat){
  y_i <- dat$mu1_hat*dat$n_study
  n <- length(y_i)
  probs <- 0.1+0.9*dat$H1
  hans_hur <- sum(y_i/probs)/sum(1/probs*dat$n_study)
}

hans_hur <- hans_hurwitz_estimator(clt_mix)
hans_hur

simple_weighted_mean <- function(dat){
  y_i <- dat$mu1_hat*dat$n_study
  n <- length(y_i)
  probs <- 0.1+0.9*dat$H1
  weigh_mean <- sum(y_i/probs)/sum(1/probs*dat$n_study)
}

weigh_mean <- simple_weighted_mean(clt_mix)
weigh_mean

calc_pub_prob <- function(z,alph){
  pub_prob <- 0.1+ifelse(z>qnorm(alph,mean=0,sd=1,lower.tail=FALSE),0.9,0)
  return(pub_prob)
}

exp_pub_prob <- function(mu,alph){
  quant <- qnorm(alph,0,1,lower.tail=FALSE)
  expect <- 0.1*pnorm(quant,mu,1)+1*pnorm(quant,mu,1,lower.tail = FALSE)
  return(expect)
}
exp_pub_prob(0,0.05)

likeli <- function(z,mu){
  lik <- dnorm(z,mean=mu,1) #same as: dnorm(z-mu,0,1)
  return(lik)
}

trunc_sample_process <- function(z,mu,alph){
  lik <- calc_pub_prob(z,alph)/exp_pub_prob(mu,alph)*likeli(z,mu)
  return(lik)
}

mus <- seq(0,2,by=0.01)

z <- clt_mix$Tn

test <- sapply(mus,function(mu) prod(trunc_sample_process(z,mu,alph=0.05)))
plot(mus,test)
T_corr <- mus[which(max(test)==test)]
likelihood <- function(){
  
}

### Calculate publication probability of studies according to Andrew & Kasy: 


#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
require("fitdistrplus")
plotdist(clt_mix$mu1_hat,histo=TRUE,demp=TRUE)



