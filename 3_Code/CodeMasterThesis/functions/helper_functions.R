## functions needed for simulation of binomial evidence


## functions needed for simulation of evidence in T statistic
vst_var_known <- function(mu0,mu1,sgm0,n_study){
  Tn <- (mu1-mu0)/sgm0*sqrt(n_study)
}

vst_var_est <- function(mu0,mu1,sgm_est,n_study,corr=F){
  if (corr==T){
    Tn <- sqrt(2*n_study)*asinh((mu1-mu0)/(sgm_est*sqrt(2)))*(1-0.7/(n_study-1)) 
  } else {
    Tn <- sqrt(2*n_study)*asinh((mu1-mu0)/(sgm_est*sqrt(2))) 
  }
}

calc_T <- function(mu0s,mu1s,sgm0s,n_study,func){
  if (missing(sgm0s)){
    Tn <- apply(mu1s, 1, function(p1) sapply(mu0s, function(p0) func(p0,p1,n_study)))
  } else {
    Tn <- sapply(1:dim(mu1s)[1], function(mu1) sapply(mu0s, function(mu0) func(mu0,mu1s[mu1,],sgm0s[mu1,],n_study)))
  }
  Tn[is.infinite(Tn)] <- NA
  return(Tn)
}



#brings data into the correct form for plottings
dat_transform <- function(T_avg,T_sd,th_emp,id,n_study,cols){
  dat <- merge(melt(T_avg),melt(T_sd),by=c("Var1","Var2"),sort=FALSE)
  dat <- cbind(dat,th_emp,id,n_study)
  colnames(dat) <- cols
  dat$mu1 <- rep(mu1s,times=1,each=length(mu0s))
  dat$mu0 <- rep(mu0s,times=length(mu1s))
  return(dat)
}

test_func <- function(val,crit_val){
  return((val>crit_val)*1)
}

ci_coverage <- function(Ts,T_avg,crit_val){
  T_avg <- as.vector(t(T_avg))
  vals <- t(abs(sapply(1:dim(Ts)[1], function(i) Ts[i,]-T_avg[i])))
  Ts_coverage <- 1-test_func(vals,crit_val)
  return(Ts_coverage)
}

T_averager <- function(Ts){
  Ts_avg <- matrix(apply(Ts,1,mean,na.rm=TRUE),length(mu0s),length(mu1s),byrow=TRUE)
  return(Ts_avg)
}

T_sd <- function(Ts){
  Ts_avg <- matrix(apply(Ts,1,sd,na.rm=TRUE),length(mu0s),length(mu1s),byrow=TRUE)
  return(Ts_avg)
}
