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
n_studies <- c(5,10,20,30,40,50,100)
n_sim <- 1e5
mu0s <- seq(-2,2,by=2)
mu1s <- seq(-2,2,by=0.1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)
evd_corr <- T #whether finite sample correction should be used or not when calculating vst
seed <- 20190504
set.seed(seed)
correction <- ifelse(evd_corr,"Corr","MLE")
name <- paste0("Ev_Student_",correction,"_",n_sim,"_",seed)

### simulate values ----
#https://stackoverflow.com/questions/34999019/apply-a-function-to-all-pairwise-combinations-of-list-elements-in-r
cols_evidence <- c("mu0","mu1","evd_mean","evd_sd","th_emp","id","n_study")
cols_tot <- c("mu0","mu1","evd_mean","evd_sd","th_emp","id","n_study","alpha","CI","power")
evidence_student <- data.table(matrix(NA,nrow=0,ncol=10))
colnames(evidence_student) <- cols_tot

i <- 1
for (n_study in n_studies) {
  start.time <- Sys.time()
  #calculate empirical and theoretical mus
  set.seed(seed)
  mu1_hats <- sapply(mu1s,function(mu1) sapply(1:n_sim, function(y) mean(rnorm(n_study,mu1,sgm0))))

  #calculate empirical and theoretical sgm
  set.seed(seed)
  #note: one could also simply sample the variances from a chi-square distribution
  sgm_hats <- sapply(1:length(mu1s),
                     function(mu1) sapply(1:n_sim, 
                                          function(sim) sqrt(1/(n_study-1)*sum((rnorm(n_study,mu1s[mu1],sgm0)-mu1_hats[sim,mu1])^2))))
  sgm_th <- matrix(rep(sgm0,length(mu1_hats)),n_sim,length(mu1s),byrow=TRUE)
  
  #calculate theoretical and empiral evidence based on Zn-vst
  #this gives the distribution of the values if the Zn holds, i.e. if we know the sd
  #question: can I use empirical SE to calculate Zn_th > how is distribution of this statistic?
  Zn_emp <- calc_evidence(mu0s,mu1_hats,sgm_hats,n_study,z_stat_norm,evd_corr)
  Zn_emp_avg <- avg_evidence(Zn_emp,mu0s,mu1s)
  Zn_emp_sd <- sd_evidence(Zn_emp,mu0s,mu1s)
  
  Zn_th_avg <- sapply(mu1s,function(mu1) z_stat_norm(mu0s,mu1,sgm0,n_study))
  Zn_th_sd <- matrix(1,nrow=dim(Zn_th_avg)[1],ncol=dim(Zn_th_avg)[2])
  
  Zn <- rbind(dat_transform(Zn_emp_avg,Zn_emp_sd,"emp","Zn",n_study,cols_evidence,mu0s,mu1s),
                 dat_transform(Zn_th_avg,Zn_th_sd,"th","Zn",n_study,cols_evidence,mu0s,mu1s))
  
  #calculate theoretical and empirical evidence based on clt-vst, but using empirical sd instead of theoretical
  #this gives us values if we assume the clt holds, but we do not know the sd and just use the empirical sd instead
  
  Tn_emp <- calc_evidence(mu0s,mu1_hats,sgm_hats,n_study,z_stat_norm,evd_corr=F)
  
  #Tn_th_avg <- avg_evidence(Tn_th)
  Tn_emp_avg <- avg_evidence(Tn_emp,mu0s,mu1s)
  Tn_emp_sd <- sd_evidence(Tn_emp,mu0s,mu1s)
  
  Tn_th_avg <- Zn_th_avg
  Tn_th_sd <- Zn_th_sd
  #note that the "theoretical" distribution of Tn is the same as the theoretical
  #distribution of Zn, which admittedly is pretty nonsensical in this context
  Tn <-  rbind(dat_transform(Tn_emp_avg,Tn_emp_sd,"emp","Tn",n_study,cols_evidence,mu0s,mu1s),
                       dat_transform(Tn_th_avg,Tn_th_sd,"th","Tn",n_study,cols_evidence,mu0s,mu1s))

  #calculate theoretical and empiral evidence based on Student-t vst
  #question: does it make sense to use theoretical mu1 but empirical sigma_est? 
  Vn_emp <- calc_evidence(mu0s,mu1_hats,sgm_hats,n_study,vst_student,evd_corr=evd_corr)
  Vn_emp_avg <- avg_evidence(Vn_emp,mu0s,mu1s)
  Vn_emp_sd <- sd_evidence(Vn_emp,mu0s,mu1s)
  
  Vn_th_avg <- sapply(mu1s,function(mu1) vst_student(mu0s,mu1,sgm0,n_study))
  Vn_th_sd <- Zn_th_sd
  
  Vn <- rbind(dat_transform(Vn_emp_avg,Vn_emp_sd,"emp","Vn",n_study,cols_evidence,mu0s,mu1s),
                 dat_transform(Vn_th_avg,Vn_th_sd,"th","Vn",n_study,cols_evidence,mu0s,mu1s))
  
  #calculate power & confidence intervals
  for (alph in alphas){
    #confidence_intervals
    crit_val_ci_norm <- qnorm((1-alph/2),0,1)
    crit_val_ci_stud <- qt((1-alph/2),n_study-1,0)
    
    emp_ci_Zn <- avg_evidence(ci_coverage(Zn_emp,Zn_th_avg,crit_val_ci_norm),mu0s,mu1s)
    emp_ci_Tn <- avg_evidence(ci_coverage(Tn_emp,Zn_th_avg,crit_val_ci_stud),mu0s,mu1s)
    emp_ci_Vn <- avg_evidence(ci_coverage(Vn_emp,Vn_th_avg,crit_val_ci_norm),mu0s,mu1s)
    
    th_ci <- matrix((1-alph),nrow=dim(emp_ci_Zn)[1],ncol=dim(emp_ci_Zn)[2])
    
    ci_Zn <- rbind(melt(emp_ci_Zn),melt(th_ci))$value
    ci_Tn <- rbind(melt(emp_ci_Tn),melt(th_ci))$value
    ci_Vn <- rbind(melt(emp_ci_Vn),melt(th_ci))$value
    
    Zn <- cbind(Zn,alph,ci_Zn)
    Tn <- cbind(Tn,alph,ci_Tn)
    Vn <- cbind(Vn,alph,ci_Vn)
    
    colnames(Vn) <- c(cols_evidence,"alpha","CI")
    
    #calculate theoretical and empirical power for test using the exact noncentral Student t test
    crit_val_stud <- qt(alph,n_study-1,0,lower.tail=FALSE)
    
    pows_avg_Tn_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pt(crit_val_stud,n_study-1,z_stat_norm(mu0,mu1,sgm0,n_study))))
    pows_avg_Tn_emp <- avg_evidence(1*(Tn_emp>crit_val_stud),mu0s,mu1s)
    
    pows_avg_Tn <- rbind(melt(pows_avg_Tn_emp),melt(pows_avg_Tn_th))$value
    
    #calculate theoretical and empirical power for test using the Tn and the normal distribution
    #"abused" Z-test, because Z was calculated using the empirical variance, not the theoretical
    #note: I always use sd=1 to calculate power > but: this assumption only holds in certain areas (question)
    crit_val_norm <- qnorm(alph,mean=0,sd=1,lower.tail=FALSE)
    
    pows_avg_Zn_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,z_stat_norm(mu0,mu1,sgm0,n_study),sd=1)))
    pows_avg_Zn_emp <- avg_evidence(1*(Zn_emp>crit_val_norm),mu0s,mu1s)
    
    pows_avg_Zn <- rbind(melt(pows_avg_Zn_emp),melt(pows_avg_Zn_th))$value
    
    #calculate theoretical and empirical power for test using the vst (non-central Student t distribution)
    #question: should I be using the theoretical or empirical variance here?
    pows_avg_Vn_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,vst_student(mu0,mu1,sgm0,n_study),sd=1)))
    pows_avg_Vn_emp <- avg_evidence(1*(Vn_emp>crit_val_norm),mu0s,mu1s)
    
    pows_avg_Vn <- rbind(melt(pows_avg_Vn_emp),melt(pows_avg_Vn_th))$value
    
    #maybe replace sd_th with sd_emp?
    #question: in this case, isn't the maximum evidence 1?
    #question: what standard deviation do we use? the empirical or the theoretical?
    #question: how do we know that approximate standard normality holds for vst?
    #question: is it prudent to use the theoretical mu1s to calculate the theoretical power for the student vst?
    
    #combine coverages & power to dataframe
    Zn <- cbind(Zn,pows_avg_Zn)
    colnames(Zn) <- cols_tot
    
    Tn <- cbind(Tn,pows_avg_Tn)
    colnames(Tn) <- cols_tot
    
    Vn <- cbind(Vn,pows_avg_Vn)
    colnames(Vn) <- cols_tot
    
    evidence_student <- rbind(evidence_student,Zn,Tn,Vn)
  }
  
  diff.time <- Sys.time()-start.time
  print(paste0("round ",i," is done in ",round(diff.time,2), " ",attributes(diff.time)$units))
  i <- i+1
}

save(evidence_student, file = paste0("data/",name,".RData")) #save data so that it doesn't have to be computed again and again