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

clt <- select_studies(dat,probs,n_studies_selected,n_select=n_tot,sample_seed,T_id="clt")
vst <- select_studies(dat,probs,n_studies_selected,n_select=n_tot,sample_seed,T_id="vst")
stud <- select_studies(dat,probs,n_studies_selected,n_select=n_tot,sample_seed,T_id="stud")

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
clt_mix <- rbind(clt_sig,select_studies(clt[H1==0,],sel_prob,n_studies_selected,n_select=NULL,sample_seed,T_id="clt"))
vst_mix <- rbind(vst_sig,select_studies(vst[H1==0,],sel_prob,n_studies_selected,n_select=NULL,sample_seed,T_id="vst"))
stud_mix <- rbind(stud_sig,select_studies(stud[H1==0,],sel_prob,n_studies_selected,n_select=NULL,sample_seed,T_id="stud"))
label <- "2C"

funnel_plotter(clt_mix,vst_mix,stud_mix,xlim=x_lim,ylim=y_lim,figname=paste0(fig_name,"_",label,"_",n_tot,".pdf"),ctgs)

#function to calculate the aggregated mean of the studies
aggregate_mean <- function(dat){
  agg_mean <- sum(dat$mu1_hat*dat$n_study)/sum(dat$n_study)
  return(agg_mean)
}
aggregate_mean(stud)
aggregate_mean(clt)
aggregate_mean(vst)

### Calculate number of papers stored in file drawer based Rosenthal 1984
calc_filedrawer <- function(Z_list,alph){
  min_filedrawer <- c()
  for (Zs in Z_list){
    k <- length(Zs)
    q <- qnorm(1-alph)
    X <- (k*mean(Zs)/q)^2-k
    #X <- (k/z^2)*(k*mean(Zs)^2-z^2)
    min_filedrawer <- c(min_filedrawer,X)
  }
  return(min_filedrawer)
}

#problem: Filedrawer problem only works for checking whether there is no null effect; it doesn't work in situations in which there really is an effects
Z_list <- list(clt$Tn,vst$Tn,qnorm(stud$p,0,1,lower.tail=FALSE),clt_sig$Tn,vst_sig$Tn,qnorm(stud_sig$p,0,1,lower.tail=FALSE))
X_list <- calc_filedrawer(Z_list,alphas)


### Calculate publication probability based on publication probability
#Hansen-Hurwitz estimator https://newonlinecourses.science.psu.edu/stat506/node/15/

calc_pub_prob <- function(dat,alph,p){
  z <- dat$Tn
  pub_prob <- p+ifelse(z>qnorm(alph,mean=0,sd=1,lower.tail=FALSE),1-p,0)
  return(pub_prob)
}

reweighted_mean <- function(dat,probs){
  y_i <- dat$mu1_hat*dat$n_study
  n <- length(y_i)
  if (missing(probs)){
    probs <- 0.1+0.9*dat$H1 
  }
  weigh_mean <- sum(y_i/probs)/sum(1/probs*dat$n_study)
  return(weigh_mean)
}

weigh_mean <- reweighted_mean(clt_mix)
weigh_mean

#Andrews-Kasy maximum likelihood estimator
exp_pub_prob <- function(mu,alph,x){
  quant <- qnorm(alph,0,1,lower.tail=FALSE)
  expect <- x*pnorm(quant,mu,1)+1*pnorm(quant,mu,1,lower.tail = FALSE)
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

# trunc_sample_process <- function(z,x,sgm,mu,alph){
#   lik <- calc_pub_prob(z,alph)/exp_pub_prob(mu,alph)*likeli(x,sgm,mu)
#   return(lik)
# }

mus <- seq(-4,4,by=0.01)
z <- clt_mix$Tn
x <- clt_mix$mu1_hat
sgm <- clt_mix$sgm_th/sqrt(clt_mix$n_study)

test <- sapply(mus,function(mu) prod(trunc_sample_process(z,mu,alph=0.05)))
plot(mus,test)
T_corr <- mus[which(max(test)==test)]
T_corr

Tn_to_mu_hat <- function(T_corr,dat){
  mu_hat <- T_corr*sqrt(sum(dat$n_study*dat$sgm_th^2))/sum(dat$n_study)
  return(mu_hat)
}

Tn_to_mu_hat(T_corr,clt_mix)

#question: how to translate corrected T intot correcte mu1_hats?
agg_sgm <- function(dat){
  var_mu <- (dat$sgm_th/sqrt(dat$n_study))^2
  agg_var <- sum(var_mu*(dat$n_study)^2)/(sum(dat$n_study))^2
  return(sqrt(agg_var))
}

sgm_agg <- agg_sgm(clt_mix)

#implement same as above but in order to correct mu1_hats:
likeli <- function(x,mu){
  lik <- dnorm(x,mean=mu,2) #same as: dnorm(z-mu,0,)
  return(lik)
}

#calculate total T

calc_T_tot <- function(dat){
  k <- length(dat$mu1_hat)
  mu_hat_tot <- sum(dat$mu1_hat*dat$n_study)/sum(dat$n_study)
  sgm_tot <- sum((dat$sgm_th*dat$n_study)^2)/sum(dat$n_study)^2
  T_tot <- mu_hat_tot/sgm_tot*sqrt(k)
  return(T_tot)
}

calc_T_tot(clt_mix)

aggregate_mean(clt_mix)

#calculate median bias
median_bias <- function(mu,alph,n){
  z <- rnorm(n,mean=mu,1)
  p <-calc_pub_prob(z,alph)
  idx <- which(1==ifelse(p==1,1,rbinom(n-sum(p==1),1,0.1)))
  med_bias <- median(z[idx])-mu
  return(med_bias)
}

thetas <- seq(0,5,by=0.01)
plot(thetas,sapply(thetas,function(theta) median_bias(theta,alph=0.05,n=1e6)),type="l")

### Calculate publication probability of studies according to Andrew & Kasy: 
#Assumption1: distribution of results X* in latent studies given the true effects
#Theta*, f(x*|Theta*) is known

#systemic replication studies

#https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
require("fitdistrplus")
plotdist(clt_mix$mu1_hat,histo=TRUE,demp=TRUE)
plotdist(clt_mix$Tn,histo=TRUE,demp=TRUE)
trunc_dens <- fitdist(clt_mix$mu1_hat,"norm")
trunc_dens_T <- fitdist(clt_mix$Tn,"norm")

#fit density: https://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf
trunc_dens_T <- density(clt_mix$Tn)
x <- trunc_dens_T$x
lik_trunc <- trunc_dens_T$y

lookup_id <- function(x,dens){
  id <- which(abs(x-dens$x)==min(abs(x-dens$x)))
  return(id)
}

calc_emp_density <- function(dat){
  z <- dat$Tn
  trunc_dens_emp <- density(z)
  ids <- sapply(z,function(z_i) lookup_id(z_i,trunc_dens_emp))
  z_dens_emp <- trunc_dens_emp$y[ids]
  return(z_dens_emp)
}

calc_th_density <- function(z,mu,probs,exp_p){
  return(likeli(z,mu)*probs/exp_p)
}

#install.packages("NMOF")
require("NMOF")

trunc_density_likeli <- function(x_mu,dat,alph){
  x <- x_mu[[1]]
  mu <- x_mu[[2]]
  z <- dat[alpha==alph,]$Tn
  probs <- calc_pub_prob(z,alph,x)
  exp_p <- exp_pub_prob(mu,alph,x)
  trunc_dens_th <- calc_th_density(z,mu,probs,exp_p)
  trunc_dens_emp <- calc_emp_density(dat)
  err <- sum(abs(trunc_dens_th-trunc_dens_emp))
  return(err)
}

z <- clt_mix$Tn
probs <- calc_pub_prob(z,0.05,0.1)
exp_p <- exp_pub_prob(1.67,0.05,0.1)
th_density <- calc_th_density(clt_mix$Tn,1.67,probs,exp_p)
emp_density <- calc_emp_density(clt_mix)

plot(z,calc_emp_density(clt_mix),ylim=c(0,1))
points(z,th_density,col="red")
plot(z,th_density)

trunc_density_likeli(list(0.1,1),clt_mix,0.05)

start.time <- Sys.time()
errs <- gridSearch(fun=trunc_density_likeli,levels=list(x,mus),dat=clt_mix,alph=0.05,method="multicore")
print(Sys.time()-start.time)

test_x <- seq(-2,2,by=0.01)
test <- likeli(test_x,0.2)*ifelse(test_x>1.64,1,0.1)/exp_pub_prob(0.2,0.05,0.1)
plot(test_x,test)

a <- optim(c(0.5,0),trunc_density_likeli,dat=clt_mix,alph=0.05)

x <- seq(0,1,by=0.01)
mus <- seq(-3,3,by=0.01)


plot(trunc_density_sim(clt_mix,2,0.1,0.05))


plot(x,calc_pub_prob(x,trunc_dens_T,mean(clt_mix$Tn)))
plot(z,calc_pub_prob(z,trunc_dens_T,0))

pub_prob <- calc_pub_prob(z,trunc_dens_T,0.2)

rew_mean <- reweighted_mean(clt_mix,pub_prob)
rew_mean

aggregate_mean(clt_mix)



pub_probs_lik <- sapply(mus, function(mu) prod(calc_pub_prob(x,trunc_dens_T,mu)))
T_corr <- mus[which(max(pub_probs_lik)==pub_probs_lik)]
T_corr

plot(x,calc_pub_prob(x,trunc_dens_T,2))

test <- sapply(mus,function(mu) prod(calc_pub_prob(z,lik_trunc,mu)))
plot(mus,test)
T_corr <- mus[which(max(test)==test)]
T_corr

equ <- function(z,mu){
  res <- dnorm(z-mu,0,1)*(z-mu)*((z-mu)-1)+1
  return(res)
}

results <- sapply(mus,function(mu) equ(z[1],mu))
plot(mus,results,type="l")
id <- which(abs(0-results)==min(abs(0-results)))
equ(z[1],2)

#draw from kernel density estimator: https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate?rq=1

#trim and fill method

trim_and_fill <- function(dat){
  n <- length(dat$Tn)
  T_sorted <- sort(dat$Tn)
  k0 <- 0
  k0_new <- 0
  find_ranks <- function(T_centered,avg_T){
    ix <- sort(abs(T_centered),index.return=T)$ix
    ranks <- sign(T_centered[ix])*sort(ix)
    return(ranks)
  }
  i <- 0
  while(TRUE){
    k0 <- k0_new
    if (k0==0){
      avg_T <- mean(T_sorted)
    } else {
      T_trimmed <- head(T_sorted,-k0)
      avg_T <- mean(T_trimmed)
    }
    T_centered <- T_sorted-avg_T
    ranks <- find_ranks(T_centered,avg_T)
    S_rank <- sum(ranks[ranks>0])
    k0_new <- round((4*S_rank-n*(n+1))/(2*n-1))
    if (k0_new <= k0){
      break
    }
  }
  T_filled <- c(T_centered,-tail(T_centered[T_centered>0],k0))
  to_add <- tail(T_sorted,k0)
  to_add <- dat[Tn %in% to_add,]
  to_add$mu1_hat <- -to_add$mu1_hat
  dat_filled <- rbind(dat,to_add)
  return(dat_filled)
}

clt_mix_filled <- trim_and_fill(clt_mix)
aggregate_mean(clt_mix_filled)
aggregate_mean(clt_mix)
#mean(clt_mix$mu1_hat)
#mean(clt_mix_filled$mu1_hat)
#trim and fill only works if assumption hold that most extreme values are omitted


#trim_and_fill with known selection probability

trim_and_fill_prob <- function(dat,p=0.1,alph){
  n <- length(dat$Tn)
  T_sorted <- sort(dat$Tn)
  k0_new <- min(sum(T_sorted>qnorm(alph,0,1,lower.tail=FALSE)),round(sum(T_sorted<qnorm(alph,0,1,lower.tail=FALSE))/p))
  find_ranks <- function(T_centered,avg_T){
    ix <- sort(abs(T_centered),index.return=T)$ix
    ranks <- sign(T_centered[ix])*sort(ix)
    return(ranks)
  }
  i <- 0
  while(TRUE){
    k0 <- k0_new
    if (k0==0){
      avg_T <- mean(T_sorted)
    } else {
      T_trimmed <- head(T_sorted,-k0)
      avg_T <- mean(T_trimmed)
    }
    T_centered <- T_sorted-avg_T
    ranks <- find_ranks(T_centered,avg_T)
    S_rank <- sum(ranks[ranks>0])
    k0_new <- round((4*S_rank-n*(n+1))/(2*n-1))
    if (k0_new == k0){
      break
    } else if (i > 100){
      break
    }
    i <- i+1
  }
  T_filled <- c(T_centered,-tail(T_centered[T_centered>0],k0))
  to_add <- tail(T_sorted,k0)
  to_add <- dat[Tn %in% to_add,]
  dat_filled <- rbind(dat,to_add)
  return(dat_filled)
}

clt_mix_filled <- trim_and_fill_prob(clt_mix,p=0.1,alph=0.05)
aggregate_mean(clt_mix_filled)
aggregate_mean(clt_mix)


##Gopas Selection
install.packages("metasens")
require("metasens")



###### Deprecated Functions ------

#Hansen-Hurwitz estimator & #Hansen-Hurwitz estimator 
#https://newonlinecourses.science.psu.edu/stat506/node/15/ https://newonlinecourses.science.psu.edu/stat506/node/15/
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


hans_hurwitz_estimator <- function(dat){
  y_i <- dat$mu1_hat*dat$n_study
  n <- length(y_i)
  probs <- 0.1+0.9*dat$H1
  hans_hur <- sum(y_i/probs)/sum(1/probs*dat$n_study)
}

hans_hur <- hans_hurwitz_estimator(clt_mix)
hans_hur


#####