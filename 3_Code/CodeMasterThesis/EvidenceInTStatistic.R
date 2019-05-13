  #import packages -----
  
  #for plotting
  require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
  
  #for data handling
  require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
  
  #custom functions
  source("./functions/helper_functions.R")
  
  # Define starting assumptions and conditions ----
  # hypotheses:
  # H0: p = mu0
  # H1: p > mu0
  n_studies <- c(5,10,20,30,40,50,100,200,500,1000)
  n_sim <- 1e4
  mu0s <- seq(-2,2,by=1)
  mu1s <- seq(-2,2,by=0.1)
  sgm0 <- sgm1 <- 2
  alphas <- c(0.05)
  T_corr <- F #whether correction should be used or not when calculating vst
  seed <- 20190504
  set.seed(seed)
  correction <- ifelse(T_corr,"Corr","MLE")
  name <- paste0("Ev_Mean_",correction,"_",n_sim,"_",seed)
  
  #simulate values ----
  #https://stackoverflow.com/questions/34999019/apply-a-function-to-all-pairwise-combinations-of-list-elements-in-r
  cols_evidence <- c("mu0","mu1","Tn","Tn_sd","th_emp","id","n_study")
  cols_tot <- c("mu0","mu1","Tn","Tn_sd","th_emp","id","n_study","alpha","CI","power")
  evidence_mean <- data.frame(matrix(NA,nrow=0,ncol=12))
  colnames(evidence_mean) <- cols_tot
  
  i <- 1
  for (n_study in n_studies) {
    start.time <- Sys.time()
    #calculate empirical and theoretical mus
    set.seed(seed)
    mu1_hats <- sapply(mu1s,function(mu1) sapply(1:n_sim, function(y) mean(rnorm(n_study,mu1,sgm0)))) #regular estimator:
    #p1_hats <- (sapply(p1s,function(x) rbinom(n_sim,n_study,x))+3/8)/(n_study+3/4) #corrected estimater - see page 160ff in book
    
    #keep for the moment, but probably not needed anymore
    #mu1s_th <- matrix(rep(mu1s,n_sim),n_sim,length(mu1s),byrow=TRUE)
    
    #calculate empirical and theoretical sgm
    set.seed(seed)
    #note: one could also simply sample the variances from a chi-square distribution
    sgm_hats <- sapply(1:length(mu1s),
                       function(mu1) sapply(1:n_sim, 
                                            function(sim) sqrt(1/(n_study-1)*sum((rnorm(n_study,mu1s[mu1],sgm0)-mu1_hats[sim,mu1])^2))))
    sgm_th <- matrix(rep(sgm0,length(mu1_hats)),n_sim,length(mu1s),byrow=TRUE)
    
    #calculate theoretical and empiral evidence based on clt-vst
    #this gives the distribution of the values if the clt holds, i.e. if we know the sd
    #question: can I use empirical SE to calculate T_clt_th > how is distribution of this statistic?
    T_clt_emp <- calc_T(mu0s,mu1_hats,sgm_th,n_study,vst_var_known)
    T_clt_emp_avg <- T_averager(T_clt_emp)
    T_clt_emp_sd <- T_sd(T_clt_emp)
    
    T_clt_th_avg <- sapply(mu1s,function(mu1) vst_var_known(mu0s,mu1,sgm0,n_study))
    T_th_sd <- matrix(1,nrow=dim(T_clt_th_avg)[1],ncol=dim(T_clt_th_avg)[2])
    
    T_clt <- rbind(dat_transform(T_clt_emp_avg,T_clt_emp_sd,"emp","clt",n_study,cols_evidence),
                   dat_transform(T_clt_th_avg,T_th_sd,"th","clt",n_study,cols_evidence))
    
    #T_avg_clt <- dat_transform_th_emp(list(T_clt_th_avg,T_clt_emp_avg))
    
    #calculate theoretical and empirical evidence based on clt-vst, but using empirical sd instead of theoretical
    #this gives us values if we assume the clt holds, but we do not know the sd and just use the empirical sd instead
    #T_clt_th_stud <- calc_T(mu0s,mu1s_th,sgm_hats,n_study,vst_var_known)
    
    T_clt_emp_stud <- calc_T(mu0s,mu1_hats,sgm_hats,n_study,vst_var_known)
    
    #T_clt_th_stud_avg <- T_averager(T_clt_th_stud)
    T_clt_emp_stud_avg <- T_averager(T_clt_emp_stud)
    T_clt_emp_stud_sd <- T_sd(T_clt_emp_stud)
    
    #not that the "theoretical" distribution of T_clt_stud is the same as the theoretical
    #distribution of T_clt, which admittedly is pretty nonsensical in this context
    T_clt_stud <-  rbind(dat_transform(T_clt_emp_stud_avg,T_clt_emp_stud_sd,"emp","clt",n_study,cols_evidence),
                         dat_transform(T_clt_th_avg,T_th_sd,"th","clt",n_study,cols_evidence))
    
    #T_avg_clt_stud <- dat_transform_th_emp(list(T_clt_emp_stud_avg),c("emp"))
    
    #calculate theoretical and empiral evidence based on Student-t vst
    #question: does it make sense to use theoretical mu1 but empirical sigma_est? 
    T_vst_emp <- calc_T(mu0s,mu1_hats,sgm_hats,n_study,vst_var_est)
    T_vst_emp_avg <- T_averager(T_vst_emp)
    T_vst_emp_sd <- T_sd(T_vst_emp)
    
    T_vst_th_avg <- sapply(mu1s,function(mu1) vst_var_est(mu0s,mu1,sgm0,n_study))
    
    T_vst <- rbind(dat_transform(T_vst_emp_avg,T_vst_emp_sd,"emp","vst",n_study,cols_evidence),
                   dat_transform(T_vst_th_avg,T_th_sd,"th","vst",n_study,cols_evidence))
    
    #T_avg_vst <- dat_transform_th_emp(list(T_vst_th_avg,T_vst_emp_avg))
    
  
    #calculate power & confidence intervals
    j <- 1
    for (alph in alphas){
      #confidence_intervals
      crit_val_ci <- qnorm((1-alph/2),0,1)
      crit_val_ci_stud <- qt((1-alph/2),n_study-1,0)
      
      emp_ci_clt <- T_averager(ci_coverage(T_clt_emp,T_clt_emp_avg,crit_val_ci))
      emp_ci_clt_stud <- T_averager(ci_coverage(T_clt_emp_stud,T_clt_emp_stud_avg,crit_val_ci_stud))
      emp_ci_vst <- T_averager(ci_coverage(T_vst_emp,T_vst_emp_avg,crit_val_ci))
      
      th_ci <- matrix((1-alph),nrow=dim(emp_ci_clt)[1],ncol=dim(emp_ci_clt)[2])
      
      ci_clt <- rbind(melt(emp_ci_clt),melt(th_ci))$value
      ci_clt_stud <- rbind(melt(emp_ci_clt_stud),melt(th_ci))$value
      ci_vst <- rbind(melt(emp_ci_vst),melt(th_ci))$value
      
      T_clt <- cbind(T_clt,alph,ci_clt)
      T_clt_stud <- cbind(T_clt_stud,alph,ci_clt)
      T_vst <- cbind(T_vst,alph,ci_vst)
      
      colnames(T_vst) <- c(cols_evidence,"alpha","CI")
      
      #calculate theoretical and empirical power for test using the exact noncentral Student t test
      crit_val_stud <- qt(alph,n_study-1,0,lower.tail=FALSE)
      
      pows_avg_stud_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pt(crit_val_stud,n_study-1,vst_var_known(mu0,mu1,sgm0,n_study))))
      pows_avg_stud_emp <- T_averager(1*(T_clt_emp_stud>crit_val_stud))
      
      pows_avg_stud <- rbind(melt(pows_avg_stud_emp),melt(pows_avg_stud_th))$value
      
      #calculate theoretical and empirical power for test using the clt and the normal distribution
      #note: I always use sd=1 to calculate power > but: this assumption only holds in certain areas (question)
      crit_val_norm <- qnorm(alph,mean=0,sd=1,lower.tail=FALSE)
      
      pows_avg_clt_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,vst_var_known(mu0,mu1,sgm0,n_study),sd=1)))
      pows_avg_clt_emp <- T_averager(1*(T_clt_emp>crit_val_norm))
      
      pows_avg_clt <- rbind(melt(pows_avg_clt_emp),melt(pows_avg_clt_th))$value
      
      #calculate theoretical and empirical power for test using the vst (non-central Student t distribution)
      #question: should I be using the theoretical or empirical variance here?
      pows_avg_vst_th <- sapply(mu1s, function(mu1) sapply(mu0s, function(mu0) 1-pnorm(crit_val_norm,vst_var_est(mu0,mu1,sgm0,n_study),sd=1)))
      pows_avg_vst_emp <- T_averager(1*(T_vst_emp>crit_val_norm))
      
      pows_avg_vst <- rbind(melt(pows_avg_vst_emp),melt(pows_avg_vst_th))$value
  
      #maybe replace sd_th with sd_emp?
      #question: in this case, isn't the maximum evidence 1?
      #question: what standard deviation do we use? the empirical or the theoretical?
      #question: how do we know that approximate standard normality holds for vst?
      #question: is it prudent to use the theoretical mu1s to calculate the theoretical power for the student vst?
      
      #combine coverages & power to dataframe
      T_clt <- cbind(T_clt,pows_avg_clt)
      colnames(T_clt) <- cols_tot
      
      T_clt_stud <- cbind(T_clt_stud,pows_avg_stud)
      colnames(T_clt_stud) <- cols_tot
  
      T_vst <- cbind(T_vst,pows_avg_vst)
      colnames(T_vst) <- cols_tot
  
      evidence_mean <- rbind(evidence_mean,T_clt,T_clt_stud,T_vst)
    }
    
    diff.time <- Sys.time()-start.time
    print(paste0("round ",i," is done in ",round(diff.time,2), " ",attributes(diff.time)$units))
    i <- i+1
  }
  
  #save data so that it doesn't have to be computed again and again:
  save(evidence_mean, file = paste0("data/",name,".RData"))



pdf(paste0("./figs/",name,".pdf"),onefile=TRUE)
i <- 1
for (sim in sims){
  T_clt <- sim[[1]]
  T_clt_sd_emp <- sim[[2]]
  T_vst <- sim[[3]]
  
  Ts <- merger(dats =list(T_clt,T_clt_sd_emp,T_vst),ids=c("T_clt","T_clt_sd_emp","T_vst"))
  #Ts <- Ts[Ts$th_emp=="th",]
  
  T_plot <- ggplot(data=Ts,aes(x=mu1,y=Tn/sqrt(n_studies[i]),col=id,group=interaction(mu0,th_emp,id),linetype=th_emp)) + 
    geom_line() + scale_color_manual(values=c("T_clt"=2,"T_vst"=3,"T_clt_sd_emp"=4))+
    geom_point(data=Ts[Ts$th_emp=="emp",], aes(shape=factor(mu0)),alpha=0.5)+
    ylim(-3,3) + xlim(min(mu1s),max(mu1s)) + ggtitle(bquote("Study size"==.(n_studies[i]))) + labs(x= expression(mu[1]), y = "Evidence")
  print(T_plot)
  
  
  #+ scale_linetype_manual(values=c("T_clt_sd_emp"=2,"T_clt"=1,"T_vst"=4)) +
  #scale_color_manual(values=c("T_clt"=2,"T_vst"=3,"T_clt_sd_emp"=4))+
  
  #plot the power
  pows <- sim[[4]]
  j <-1
  for (pow in pows){
    pow_clt <- pow[[1]] 
    pow_clt_sd_emp <- pow[[2]] #Student distribution!
    pow_vst <- pow[[3]]
    
    pows_merged <- merger(dats=list(pow_clt,pow_clt_sd_emp,pow_vst),ids=c("pow_clt","pow_sd_emp","pow_vst"))
    #pows_merged <- pows_merged[pows_merged$th_emp=="th",]
    
    cutoff <- data.frame(x = c(min(mu1s), max(mu1s)), y = alphas[j], cutoff = factor(alphas[j]))
    mu0_cutoff <- data.frame(x = rep(mu0s,each=2),y=rep(c(0,1),length(mu0s)))
    
    pow_plot <- ggplot(data=pows_merged,aes(x=mu1,y=Tn,col=id,group=interaction(mu0,th_emp,id),linetype=th_emp)) + 
      geom_line() + scale_color_manual(values=c("pow_clt"=2,"pow_vst"=3,"pow_sd_emp"=4))+ 
      geom_point(data=pows_merged[pows_merged$th_emp=="emp",], aes(shape=factor(mu0)),alpha=.5) +
      geom_line(aes(x, y,group=factor(x)),mu0_cutoff,inherit.aes=FALSE,alpha=.2) + geom_line(aes(x, y),cutoff,inherit.aes=FALSE,alpha=.2) +
      ggtitle(bquote("Study size"==.(n_studies[i])~", alpha"==.(alphas[j]))) + labs(x= expression(mu[1]), y = "power")
    print(pow_plot)
    
    #+ scale_linetype_manual(values=c("pow_sd_emp"=2,"pow_clt"=1,"pow_vst"=4)) + scale_alpha_manual(values=c("emp"=1,"th"=0.5))+ 
    j <- j+1
  }
  
  
  i <- i+1
}
dev.off()

#questions: 
# 1) why is fit between the two vst better for larger n (is it central limit theorem?)
# 2) why is fit wors a the extremes (-2,2) as opposed to the middle section (0) and why does it get worse
#    the more extreme the values get?
# 3) why is corrected vst *worse* than uncorrected?
# 4) when using the clt-vst with the empirical sd, the values converge quickly to the clt-vst with the theoretical sd > why is studen-t-vst necessary?
# 5) why does clt_sd_emp (student t distribution) have a higher evidence but a lower power?
# 6) why does correction deteriorate fit? > it is supposed to increas it because it is supposed to match the p-values with the regular p-values
# A6) it does not! By correcting the the value, it makes the normal approximation better; this way, when calculating the power, the power calculated
#     based on the vst is closer to the power calculated based on the clt_sd_emp (student t); however, this only seems to hold for the empirical, not
#     the theoretical values
# 7) why is correction better for empirical than for theoretical values?
# 8) why does clt show higher power than both clt_sd_emp and vst? Probably because clt gets this addition of power at the expense of Type I error?

#to add:
# 1) error between pow_stud and pow_vst? 


#discontinued functions ----
# dat_transform_th_emp <- function(dats,th_emps=c("th", "emp")){
#   for (i in 1:length(th_emps)){
#     dat <- dat_transform(dats[[i]])
#     th_emp <- rep(th_emps[i],dim(dat)[1])
#     dat$th_emp <- th_emp
#     dats[[i]] <- dat
#   }
#   dats <- do.call("rbind",dats)
#   return(dats)
# }

# dat_transform <- function(T_avg){
#   dat <- melt(T_avg)
#   colnames(dat) <- c("mu0","mu1","Tn")
#   dat$mu1 <- rep(mu1s,times=1,each=length(mu0s))
#   dat$mu0 <- rep(mu0s,times=length(mu1s))
#   return(dat)
# }

# merger <- function(dats,ids){
#   for (i in 1:length(ids)){
#     id <- rep(ids[i],dim(dats[[i]])[1])
#     dat <- dats[[i]]
#     dat$id <- id
#     dats[[i]] <- dat
#   }
#   dats <- do.call("rbind",dats)
#   return(dats)
# }