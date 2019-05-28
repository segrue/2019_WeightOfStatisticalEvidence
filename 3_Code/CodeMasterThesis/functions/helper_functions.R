### packages needed for functions below ----
#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
require("data.table")

#for plotting
require("ggplot2")
require("gridExtra")
#require("cowplot") #needed for blankplots, currently not used
require("RColorBrewer")

### functions needed for simulation of evidence in multiple scripts ----
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

calc_T <- function(mu0s,mu1s,sgm1s,n_study,func,T_corr){
  if (missing(sgm1s)){
    Tn <- apply(mu1s, 1, function(p1) sapply(mu0s, function(p0) func(p0,p1,n_study)))
  } else {
    if (missing(T_corr)){
      Tn <- sapply(1:dim(mu1s)[1], function(mu1) sapply(mu0s, function(mu0) func(mu0,mu1s[mu1,],sgm1s[mu1,],n_study)))  
    } else {
      Tn <- sapply(1:dim(mu1s)[1], function(mu1) sapply(mu0s, function(mu0) func(mu0,mu1s[mu1,],sgm1s[mu1,],n_study,corr=T_corr)))
    }
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

### functions needed for plotting in several scripts
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  leg <- tmp$grobs[[leg]]
  return(leg)
}

### functions needed for the script "SimulateStudies.R" ----

# brings data into the correct form to be handled in "SimulateStudies.R"
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
evidence_in_mean <- function(mu0_vec,mu1s,sgm0,alphas,T_corr=F,seed,n_studies,n_sim,fig_name){
  cols_tot <- c("sim","mu1","mu1_hat","sgm_hat","sgm_th","Tn","p","id","T_avg_th","T_avg_emp",
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
      
      p_clt <- pnorm(T_clt_emp,0,1,lower.tail=FALSE)
      T_clt_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_clt_emp),t(p_clt)),id="clt"),T_clt_averages)
      
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
      
      p_clt_stud <- pt(T_clt_emp_stud,n_study-1,0,lower.tail=FALSE)
      
      T_clt_stud_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_clt_emp_stud),t(p_clt_stud)),id="stud"),T_clt_stud_averages)

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
      
      p_vst <- pnorm(T_vst_emp,0,1,lower.tail=FALSE)
      T_vst_df <- cbind(mat2df(mats=list(mu1_hats,sgm_hats,sgm_th,t(T_vst_emp),t(p_vst)),id="vst"),T_vst_averages)
      
      
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

### functions needed for script "PublicationBias.R"
funnel_plotter <- function(clt,vst,stud,xlim=2,ylim=8,figname,ctgs){
  #define colour code & legend that is consistent over all plots with "ctgs" (regardless of whether they appear in plot or not)
  n_ctgs <- length(ctgs)
  clrs <- brewer.pal(n_ctgs,"Paired")
  names(clrs) <- ctgs
  sc_col <- scale_color_manual(values=clrs, drop=FALSE)
  
  std_errs_inv <- seq(0.1,ylim,by=0.1)
  x_cutoff <- qnorm(0.05,mean=0,lower.tail=F)/std_errs_inv
  cutoff <- data.table(x_cutoff,std_errs_inv)
  cutoff_norm <- data.table(x_cutoff=rep(qnorm(0.05,mean=0,lower.tail=F),length(std_errs_inv)),std_errs_inv)
  

  #create plots
  p1 <- ggplot(data=clt,aes(x=Tn,y=n_study,col=factor(n_study))) + geom_point() + xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(paste(italic("z"),"-statistic")),y= "n") + sc_col + ylim(0,max(ctgs)+50)
  p2 <- ggplot(data=clt,aes(x=Tn,y=1/(sgm_hat/sqrt(n_study)),col=factor(n_study))) + geom_point()+ xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(paste(italic("z"),"-statistic")),y= expression(1/SE(bar(X)))) + sc_col + 
        ylim(0,ylim) + geom_line(data=cutoff_norm,aes(x=x_cutoff,y=std_errs_inv),inherit.aes=F,linetype=2)
  p3 <- ggplot(data=stud,aes(x=Tn,y=n_study,col=factor(n_study))) + geom_point() + xlim(-xlim,xlim)+
        theme(legend.position = "none") + labs(x = expression(paste(italic("t"),"-statistic")),y= "n") + sc_col + ylim(0,max(ctgs)+50)
  p4 <- ggplot(data=stud,aes(x=Tn,y=1/(sgm_hat/sqrt(n_study)),col=factor(n_study))) + geom_point()+ xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(paste(italic("t"),"-statistic")),y= expression(1/SE(bar(X)))) + sc_col +
        ylim(0,ylim)
  p5 <- ggplot(data=vst,aes(x=Tn,y=n_study,col=factor(n_study))) + geom_point() + xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(italic(T[vst])),y= "n") + sc_col + ylim(0,max(ctgs)+50)
  p6 <- ggplot(data=vst,aes(x=Tn,y=1/(sgm_hat/sqrt(n_study)),col=factor(n_study))) + geom_point()+ xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(italic(T[vst])),y= expression(1/SE(bar(X)))) + sc_col +
        ylim(0,ylim) + geom_line(data=cutoff_norm,aes(x=x_cutoff,y=std_errs_inv),inherit.aes=F,linetype=2)
  p7 <- ggplot(data=stud,aes(x=mu1_hat,y=n_study,col=factor(n_study))) + geom_point()+ xlim(-xlim,xlim) +
        theme(legend.position = "none") + labs(x = expression(bar(X)),y= "n") + sc_col + ylim(0,max(ctgs)+50)
  p8 <- ggplot(data=stud,aes(x=mu1_hat,y=1/(sgm_hat/sqrt(n_study)),col=factor(n_study))) + geom_point()+ xlim(-xlim,xlim) + 
        guides(color = guide_legend(reverse = TRUE)) + labs(x = expression(bar(X)),y= expression(1/SE(bar(X))),color="n") +
        theme(legend.title.align=0.5) + sc_col + ylim(0,ylim) + geom_line(data=cutoff,aes(x=x_cutoff,y=std_errs_inv),inherit.aes=F,linetype=2)+
        geom_line(data=vst,aes(x=mu1,y=seq(0,ylim,length.out=length(mu1))),inherit.aes=F,linetype=1)
   
  #see here for more details on gridarrange: http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
  leg <- get_legend(p8)
  #blankPlot <- ggplot()+geom_blank(aes(1,1))+cowplot::theme_nothing()
  pdf(paste0("./figs/",figname),onefile=TRUE)
  grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8+theme(legend.position="none"),leg,
               nrow=4,ncol=3,widths=c(3,3,1),layout_matrix=cbind(c(1,3,5,7),c(2,4,6,8),c(9,9,9,9)))
  dev.off()
}

select_studies <- function(dat,probs,n_studies,n_select,seed,T_id) {
  set.seed(seed)
  if (sum(probs)>1){
    stop("Sum of probabilties must not exceed 1.")
  } 
  if (length(probs)!= length(n_studies) & (length(probs)>1)){
    stop("Probs must be skalra or a vector of the same as length of n_studies vector.")
  }
  # if (n_select  > n_sim*length(probs)){
  #   stop("The total number of entries to be selected exceeds total number of studies available.")
  # }
  
  dat_selected <- data.table(matrix(NA,nrow=0,ncol=ncol(dat)))
  colnames(dat_selected) <- colnames(dat)
  
  if (length(probs)>1){
    n_select <- ceiling(n_select*probs)
    for (i in 1:length(n_studies)){
      temp_dat <- dat[n_study==n_studies[i] & id==T_id,]
      n <- dim(temp_dat)[1]
      
      if (n == 0){
        next
      } else if (n < n_select[i]){
        selected_idx <- seq(1,n,1)
      } else {
        selected_idx <- sample(n,n_select[i])
      }
      dat_selected <- rbind(dat_selected,temp_dat[selected_idx,])
    }
  } else {
    temp_dat <- dat[id==T_id & n_study %in% n_studies,]
    n <- dim(temp_dat)[1]
    n_select <- ceiling(n*probs)
    selected_idx <- sample(n,n_select)
    dat_selected <- rbind(dat_selected,temp_dat[selected_idx,])
  }
  
  return(dat_selected)
}


### Functions needed diagnosing and correcting publication bias --------------

## Calculate number of papers stored in file drawer based Rosenthal 1984 ------
## problem: Filedrawer problem only works for checking whether there is no null 
## effect; it doesn't work in situations in which there really is an effects

calc_filedrawer <- function(dat){
  if (length(unique(dat$alpha))>1 | length(unique(dat$id))>1){
    stop("Only test results evaluated at the same alpha threshold 
         and using one evidence metric are permitted.")
  } else {
    alph <- dat$alpha[1]
  }
  Z <- qnorm(dat$p,0,1,lower.tail=FALSE)
  k <- length(Z)
  q <- qnorm(1-alph)
  filedrawer <- (k*mean(Z)/q)^2-k
  return(filedrawer)
}


## Calculate reweighted mean based on publication probability -----------------
## Insired by the Hansen-Hurwitz estimator: 
## https://newonlinecourses.science.psu.edu/stat506/node/15/

# Calculate publication probability of individual studies ---------------------

# calculates publication probability of studies below and above significance 
# threshold, respectively. Studies with significant results are published with 
# probability 1, studies with non-significant results are published with 
# probability p. 
# If max_N is given, studies with higher study size have probability of 
# publication of p+min(1-p,n/N), where we assume that a study of size N has a 
# probability of 1 of being published regardless of whether the results
# are significant. If max_N is set to 0, all publication probabilities 
# are set to 1, hence you can use these publications to calculate the 
# unweighted mean of the data.

calculate_pub_prob <- function(dat,p,max_N){
  if (length(unique(dat$alpha))>1 | length(unique(dat$id))>1){
    stop("Only test results evaluated at the same alpha threshold 
         and using one evidence metric are permitted.")
  } else {
    alph <- dat$alpha[1]
  }
  
  z <- dat$Tn
  if (missing(max_N)){
    pub_prob <- p+ifelse(z>qnorm(alph,mean=0,sd=1,lower.tail=FALSE),1-p,0)
  } else {
    n <- dat$n_study
    pub_prob <- sapply(1:length(z), 
                       function(i) p+ifelse(z[i]>qnorm(alph,mean=0,sd=1,lower.tail=FALSE),1-p,min(1-p,n[i]/max_N)))
  }
  return(pub_prob)
}

# Reweight mean of all studies based on publication probability ---------------

# Calculates a reweighted of the means based on study size and publication 
# probability given Evidence value; if probability values are not submitted, 
# the step-function is assumed with p = 0.1 if Z non-significant
# and p = 1 if Z is significant

reweight_mean <- function(dat, probs){
  mu_hat <- dat$mu1_hat*dat$n_study
  n <- length(mu_hat)
  
  if (missing(probs)){
    probs <- 0.1+0.9*dat$H1 
  }
  
  if (n!=length(probs)){
    stop("Vector of probabilities must be the same length as number of studiy results to be reweighted.")
  }
  
  weigh_mean <- sum(mu_hat/probs)/sum(dat$n_study/probs)
  return(weigh_mean)
}

# Check calculate_pub_prob and reweight_mean ----------------------------------
check_reweighting <- function(dat) {
  # reweight with default probability 
  weigh_mean_1 <- reweight_mean(dat)
  
  # reweight with default probability, but use "calculate_pub_prob"
  pub_prob <- calculate_pub_prob(dat, p = 0.1) 
  weigh_mean_2 <- reweight_mean(dat,pub_prob)
  
  # reweight probability based on study size
  pub_prob <- calculate_pub_prob(dat, p = 0.1,max_N = 1000) 
  weigh_mean_3 <- reweight_mean(dat,pub_prob)
  
  stopifnot(weigh_mean_1==weigh_mean_2,weigh_mean_1!=weigh_mean_3)
  print("Functions 'calculate_pub_prob' and 'reweight_mean' work as expected.")
  return(list(weigh_mean_1,weigh_mean_3))
}

## Trim and fill methods and variations thereof -------------------------------

# Traditional trim and fill method according to Duval and Tweedie 2000 --------
trim_and_fill <- function(dat,pub_prob) {
  if (length(unique(dat$alpha))>1 | length(unique(dat$id))>1){
    stop("Only test results evaluated at the same alpha threshold 
         and using one evidence metric are permitted.")
  } else {
    alph <- dat$alpha[1]
  }
  
  # helper function to find signed ranks
  find_ranks <- function(T_centered,avg_T){
    ix <- sort(abs(T_centered),index.return=T)$ix
    ranks <- sign(T_centered[ix])*sort(ix)
    return(ranks)
  }
  
  n <- length(dat$Tn)
  T_sorted <- sort(dat$Tn,index.return=T)
  k0 <- 0
  
  if (missing(pub_prob)){
    k0_new <- 0
  } else if (length(pub_prob) == 1) {
    # if publication probability is known, trim off number of actually omitted
    # studies in the beginning; if this number is larger than the number of 
    # significant results, trim all significant results off.
    k0_new <- min(sum(T_sorted$x>qnorm(alph,0,1,lower.tail=FALSE)),
                  round(sum(T_sorted$x<qnorm(alph,0,1,lower.tail=FALSE))/pub_prob))
  } else {
    pub_prob_sorted <- pub_prob[T_sorted$ix]
    significant <- T_sorted$x>qnorm(alph,0,1,lower.tail=FALSE)
    k0_new <- min(sum(significant),
                  round(sum(rep(1,n)[!significant]/pub_prob_sorted[!significant])))
  }
  
  i <- 0
  while(TRUE){
    diff <- k0-k0_new
    k0 <- k0_new #number of supposedly omitted studies
    
    if (k0==0){
      avg_T <- mean(T_sorted$x)
    } else {
      T_trimmed <- head(T_sorted$x,-k0)
      avg_T <- mean(T_trimmed)
    }
    
    T_centered <- T_sorted$x-avg_T
    ranks <- find_ranks(T_centered,avg_T)
    S_rank <- sum(ranks[ranks>0])
    
    k0_new <- round((4*S_rank-n*(n+1))/(2*n-1))
    
    i <- i+1
    
    if (k0_new == k0 | i > 100 | diff == -(k0-k0_new)) {
      break
    }
  }
  
  T_filled <- c(T_centered,-tail(T_centered[T_centered>0],k0))
  to_add <- dat[Tn %in% tail(T_sorted$x,k0),] # add k0 studies with highest positive rank
  to_add$mu1_hat <- -to_add$mu1_hat # invert sign of mu_hat
  dat_filled <- rbind(dat,to_add)
  
  return(dat_filled)
}

trim_and_fill(clt_mix,pub_prob=0.1)

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