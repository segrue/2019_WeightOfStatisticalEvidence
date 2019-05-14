### functions needed for simulation of evidence in multiple scripts

#vst for binomial variable
vst <- function(p0,p1,n){
  Tn <- 2*sqrt(n)*(asin(sqrt(p1))-asin(sqrt(p0)))
  return(Tn)
}

#z-statistic based on CLT for binomial variable
clt <- function(p0,p1,n){
  Tn <- (p1-p0)/sqrt((p1*(1-p1)/n))
  return(Tn)
}


#vst for normally distributed variable - variance known (aka clt)
vst_var_known <- function(mu0,mu1,sgm0,n_study){
  Tn <- (mu1-mu0)/sgm0*sqrt(n_study)
}

#vast for normally distributed variable - variance unknown
vst_var_est <- function(mu0,mu1,sgm_est,n_study,corr=F){
  if (corr==T){
    Tn <- sqrt(2*n_study)*asinh((mu1-mu0)/(sgm_est*sqrt(2)))*(1-0.7/(n_study-1)) 
  } else {
    Tn <- sqrt(2*n_study)*asinh((mu1-mu0)/(sgm_est*sqrt(2))) 
  }
}

calc_T <- function(mu0s,mu1s,sgm1s,n_study,func){
  if (missing(sgm1s)){
    Tn <- apply(mu1s, 1, function(p1) sapply(mu0s, function(p0) func(p0,p1,n_study)))
  } else {
    Tn <- sapply(1:dim(mu1s)[1], function(mu1) sapply(mu0s, function(mu0) func(mu0,mu1s[mu1,],sgm1s[mu1,],n_study)))
  }
  Tn[is.infinite(Tn)] <- NA
  return(Tn)
}

#brings data into the correct form for plottings
dat_transform <- function(T_avg,T_sd,th_emp,id,n_study,cols,mu0s,mu1s){
  dat <- merge(melt(T_avg),melt(T_sd),by=c("Var1","Var2"),sort=FALSE)
  dat <- cbind(dat,th_emp,id,n_study)
  colnames(dat) <- cols
  dat[,1] <- rep(mu0s,times=length(mu1s))
  dat[,2] <- rep(mu1s,times=1,each=length(mu0s))
  return(dat)
}

test_func <- function(val,crit_val){
  return((val>crit_val)*1)
}

# ci_coverage <- function(Ts,T_avg,crit_val){
#   T_avg <- as.vector(t(T_avg))
#   vals <- t(abs(sapply(1:dim(Ts)[1], function(i) Ts[i,]-T_avg[i])))
#   Ts_coverage <- 1-test_func(vals,crit_val)
#   return(Ts_coverage)
# }

ci_coverage <- function(Ts,T_avg,crit_val){
  T_avg <- as.vector(t(T_avg))
  vals <- t(abs(sapply(1:dim(Ts)[1], function(i) Ts[i,]-T_avg[i])))
  Ts_coverage <- 1-test_func(vals,crit_val)
  return(Ts_coverage)
}

T_averager <- function(Ts,mu0s,mu1s){
  Ts_avg <- matrix(apply(Ts,1,mean,na.rm=TRUE),length(mu0s),length(mu1s),byrow=TRUE)
  return(Ts_avg)
}

T_sd <- function(Ts,mu0s,mu1s){
  Ts_avg <- matrix(apply(Ts,1,sd,na.rm=TRUE),length(mu0s),length(mu1s),byrow=TRUE)
  return(Ts_avg)
}

### functions needed for the script "SimulateStudies.R"

#brings data into the correct form to be handled in "SimulateStudies.R"
mat2df <- function(mats,id){
  if (length(mats) == 1){
    dat <- melt(mats[[1]]) 
  } else {
    dat <- melt(mats[[1]])
    for (i in 2:length(mats)){
      dat <- merge(dat,melt(mats[[i]]),by=c("Var1","Var2"),sort=FALSE)
    } 
  }
  dat <- cbind(dat,id)
  # colnames(dat) <- cols
  # dat[,1] <- rep(mu0s,times=length(mu1s))
  # dat[,2] <- rep(mu1s,times=1,each=length(mu0s))
  return(dat)
}

#function to simulate a lot of different data_sets to take as basis for additional calculations 
evidence_in_mean <- function(mu0_vec,mu1s,sgm0,alphas,T_corr=F,seed,n_studies,n_sim){
  cols_tot <- c("sim","mu1","mu1_hat","sgm_hat","sgm_th","Tn","id","T_avg_th","T_avg_emp",
                "T_sd_th","T_sd_emp","n_study","correction","mu0","mu1_hat_mean","sgm_hat_mean","CI","CI_crit_value","alpha","H1","H0_crit_val")
  
  evidence_df <- data.table(matrix(NA,nrow=0,ncol=length(cols_tot)))
  colnames(evidence_df) <- cols_tot
  
  correction <- ifelse(T_corr,"Corr","MLE")
  name <- paste0("DF_Ev_Mean_",correction,"_",n_sim,"_",seed)
  
  i <- 1
  for (mu0s in mu0_vec){
    for (n_study in n_studies) {
      start.time <- Sys.time()
      
      #calculate empirical and theoretical mus
      set.seed(seed)
      mu1_hats <- sapply(mu1s,function(mu1) sapply(1:n_sim, function(y) mean(rnorm(n_study,mu1,sgm0))))
      mu1_hats_mean <- T_averager(t(mu1_hats),mu0s,mu1s)
      
      #calculate empirical and theoretical sgm
      set.seed(seed)
      #note: one could also simply sample the variances from a chi-square distribution
      sgm_hats <- sapply(1:length(mu1s),
                         function(mu1) sapply(1:n_sim, 
                                              function(sim) sqrt(1/(n_study-1)*sum((rnorm(n_study,mu1s[mu1],sgm0)-mu1_hats[sim,mu1])^2))))
      sgm_hats_mean <- T_averager(t(sgm_hats),mu0s,mu1s)
      
      sgm_th <- matrix(rep(sgm0,length(mu1_hats)),n_sim,length(mu1s),byrow=TRUE)
      
      #calculate theoretical and empiral evidence based on clt-vst
      #this gives the distribution of the values if the clt holds, i.e. if we know the sd
      T_clt_emp <- calc_T(mu0s,mu1_hats,sgm_th,n_study,vst_var_known)
      
      T_clt_emp_avg <- T_averager(T_clt_emp,mu0s,mu1s)
      T_clt_emp_sd <- T_sd(T_clt_emp,mu0s,mu1s)
      
      T_clt_th_avg <- sapply(mu1s,function(mu1) vst_var_known(mu0s,mu1,sgm0,n_study))
      T_th_sd <- matrix(1,nrow=dim(matrix(T_clt_th_avg))[1],ncol=dim(matrix(T_clt_th_avg))[2])
      
      T_clt_averages <- cbind(rep(T_clt_th_avg,each=n_sim),rep(T_clt_emp_avg,each=n_sim),
                              rep(T_th_sd,each=n_sim),rep(T_clt_emp_sd,each=n_sim))
      
      T_clt_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_clt_emp)),id="clt"),T_clt_averages)
      
      #calculate theoretical and empirical evidence based on clt-vst, but using empirical sd instead of theoretical
      #this gives us values if we assume the clt holds, but we do not know the sd and just use the empirical sd instea
      T_clt_emp_stud <- calc_T(mu0s,mu1_hats,sgm_hats,n_study,vst_var_known)
      
      T_clt_emp_stud_avg <- T_averager(T_clt_emp_stud,mu0s,mu1s)
      T_clt_emp_stud_sd <- T_sd(T_clt_emp_stud,mu0s,mu1s)
      # 
      # #not that the "theoretical" distribution of T_clt_stud is the same as the theoretical
      # #distribution of T_clt, which admittedly is pretty nonsensical in this context
      # T_clt_stud <-  rbind(dat_transform(T_clt_emp_stud_avg,T_clt_emp_stud_sd,"emp","clt",n_study,cols_evidence,mu0s,mu1s),
      #                      dat_transform(T_clt_th_avg,T_th_sd,"th","clt",n_study,cols_evidence,mu0s,mu1s))
      
      
      T_clt_stud_averages <- cbind(rep(T_clt_th_avg,each=n_sim),rep(T_clt_emp_stud_avg,each=n_sim),
                                   rep(T_th_sd,each=n_sim),rep(T_clt_emp_stud_sd,each=n_sim))
      
      T_clt_stud_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_clt_emp_stud)),id="stud"),T_clt_stud_averages)
      
      #calculate theoretical and empiral evidence based on Student-t vst
      T_vst_emp <- calc_T(mu0s,mu1_hats,sgm_hats,n_study,vst_var_est)
      T_vst_df <- mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_vst_emp)),id="vst")
      
      T_vst_emp_avg <- T_averager(T_vst_emp,mu0s,mu1s)
      T_vst_emp_sd <- T_sd(T_vst_emp,mu0s,mu1s)
      
      T_vst_th_avg <- sapply(mu1s,function(mu1) vst_var_est(mu0s,mu1,sgm0,n_study))
      # 
      # T_vst <- rbind(dat_transform(T_vst_emp_avg,T_vst_emp_sd,"emp","vst",n_study,cols_evidence,mu0s,mu1s),
      #                dat_transform(T_vst_th_avg,T_th_sd,"th","vst",n_study,cols_evidence,mu0s,mu1s))
      
      T_vst_averages <- cbind(rep(T_vst_th_avg,each=n_sim),rep(T_vst_emp_avg,each=n_sim),
                              rep(T_th_sd,each=n_sim),rep(T_vst_emp_sd,each=n_sim))
      
      T_vst_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_vst_emp)),id="vst"),T_vst_averages)
      
      
      T_df <- cbind(rbind(T_clt_df,T_clt_stud_df,T_vst_df),n_study,correction,mu0s)
      T_df[,2] <- rep(mu1s,times=3,each=n_sim)
      T_df <- cbind(T_df,rep(mu1_hats_mean,times=3,each=n_sim),rep(sgm_hats_mean,times=3,each=n_sim))
      # colnames(T_df) <- c("sim","mu1","mu1_hat","sgm_hat","sgm_th","Tn","id","n_study","correction","mu0")
      
      #things to add:
      # column with critical value for CI -> done
      # column with critical value for power
      # mean over all mu1_hats 
      # mean over all sgm_hats > can be averaged, because mu and n are the same in each case
      
      #calculate power & confidence intervals
      for (alph in alphas){
        #confidence_intervals
        crit_val_ci <- qnorm((1-alph/2),0,1)
        crit_val_ci_stud <- qt((1-alph/2),n_study-1,0)
        
        #question: how to calculate empirical CI coverage
        emp_ci_clt <- t(ci_coverage(T_clt_emp,T_clt_th_avg,crit_val_ci))
        emp_ci_clt_stud <- t(ci_coverage(T_clt_emp_stud,T_clt_th_avg,crit_val_ci_stud))
        emp_ci_vst <- t(ci_coverage(T_vst_emp,T_vst_th_avg,crit_val_ci))
        
        ci_df <- rbind(cbind(melt(emp_ci_clt)$value,crit_val_ci),
                       cbind(melt(emp_ci_clt_stud)$value,crit_val_ci_stud),
                       cbind(melt(emp_ci_vst)$value,crit_val_ci))
        
        #th_ci <- matrix((1-alph),nrow=dim(emp_ci_clt)[1],ncol=dim(emp_ci_clt)[2])
        
        #ci_clt <- rbind(melt(emp_ci_clt),melt(th_ci))$value
        #ci_clt_stud <- rbind(melt(emp_ci_clt_stud),melt(th_ci))$value
        #ci_vst <- rbind(melt(emp_ci_vst),melt(th_ci))$value
        
        T_df <- cbind(T_df,ci_df,alph)
        #T_clt <- cbind(T_clt,alph,ci_clt)
        #T_clt_stud <- cbind(T_clt_stud,alph,ci_clt)
        #T_vst <- cbind(T_vst,alph,ci_vst)
        
        #colnames(T_vst) <- c(cols_evidence,"CI","alpha")
        
        #calculate theoretical and empirical power for test using the exact noncentral Student t test
        crit_val_stud <- qt(alph,n_study-1,0,lower.tail=FALSE)
        
        #pows_avg_stud_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pt(crit_val_stud,n_study-1,vst_var_known(mu0,mu1,sgm0,n_study))))
        pows_avg_stud_emp <- t(1*(T_clt_emp_stud>crit_val_stud))
        
        #pows_avg_stud <- rbind(melt(pows_avg_stud_emp),melt(pows_avg_stud_th))$value
        
        #calculate theoretical and empirical power for test using the clt and the normal distribution
        crit_val_norm <- qnorm(alph,mean=0,sd=1,lower.tail=FALSE)
        
        #pows_avg_clt_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,vst_var_known(mu0,mu1,sgm0,n_study),sd=1)))
        pows_avg_clt_emp <- t(1*(T_clt_emp>crit_val_norm))
        
        #pows_avg_clt <- rbind(melt(pows_avg_clt_emp),melt(pows_avg_clt_th))$value
        
        #calculate theoretical and empirical power for test using the vst (non-central Student t distribution)
        #pows_avg_vst_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,vst_var_est(mu0,mu1,sgm0,n_study),sd=1)))
        pows_avg_vst_emp <- t(1*(T_vst_emp>crit_val_norm))
        
        #pows_avg_vst <- rbind(melt(pows_avg_vst_emp),melt(pows_avg_vst_th))$value
        
        pows <- rbind(cbind(melt(pows_avg_clt_emp)$value,crit_val_norm),
                      cbind(melt(pows_avg_stud_emp)$value,crit_val_stud),
                      cbind(melt(pows_avg_vst_emp)$value,crit_val_norm))
        #pows <- rbind(pows_avg_clt_emp,pows_avg_stud_emp,pows_avg_vst_emp)
        
        T_df <- cbind(T_df,pows)
        colnames(T_df) <- cols_tot
        #T_clt <- cbind(T_clt,pows_avg_clt)
        #colnames(T_clt) <- cols_tot
        
        #T_clt_stud <- cbind(T_clt_stud,pows_avg_stud)
        #colnames(T_clt_stud) <- cols_tot
        
        #T_vst <- cbind(T_vst,pows_avg_vst)
        #colnames(T_vst) <- cols_tot
        
        evidence_df <- rbind(evidence_df,T_df)
      }
      
      diff.time <- Sys.time()-start.time
      print(paste0("round ",i," is done in ",round(diff.time,2), " ",attributes(diff.time)$units))
      i <- i+1
    }  
  }
  save(evidence_df, file = paste0("data/",name,".RData")) #save data so that it doesn't have to be computed again and again
  return(evidence_df)
}
