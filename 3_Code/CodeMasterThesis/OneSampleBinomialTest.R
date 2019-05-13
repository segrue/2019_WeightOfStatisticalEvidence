#import packages -----

#for plotting
require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# Define starting assumptions and conditions ----
# hypotheses:
# H0: p = p0
# H1: p > p0
n_studies <- c(5,10,20,30,40,50,100,200,500,1000)
n_sim <- 1e4
p0s <- seq(0.1,0.9,by=0.2)
p1s <- seq(0.01,0.99,by=0.01)
alphas <- c(0.05)
T_corr <- T
seed <- 20190504
set.seed(seed)
correction <- ifelse(T_corr,"Ans","MLE")
name <- paste0("Ev_Binom_",correction,"_",n_sim,"_",seed)


# Define vst and helper functions  ----
vst <- function(p0,p1,n){
  Tn <- 2*sqrt(n)*(asin(sqrt(p1))-asin(sqrt(p0)))
  return(Tn)
}

#z-statistic based on CLT
clt <- function(p0,p1,n){
  Tn <- (p1-p0)/sqrt((p1*(1-p1)/n))
  return(Tn)
}
#question: what changes when we use the theoretical variance instead of the empirical? fit should get better, shouldn't it?
#also: where do peaks come from?

calc_T <- function(p0s,p1s,n_study,func){
  Tn <- apply(p1s, 1, function(p1) sapply(p0s, function(p0) func(p0,p1,n_study)))
  Tn[is.infinite(Tn)] <- NA
  return(Tn)
}

#brings data into the correct form for plottings
dat_transform <- function(T_avg,T_sd,th_emp,id,n_study,cols){
  dat <- merge(melt(T_avg),melt(T_sd),by=c("Var1","Var2"),sort=FALSE)
  dat <- cbind(dat,th_emp,id,n_study)
  colnames(dat) <- cols
  dat$p1 <- rep(p1s,times=1,each=length(p0s))
  dat$p0 <- rep(p0s,times=length(p1s))
  return(dat)
}

dat_transform_power <- function(pow_avg,th_emp,id,alph,n_study,cols){
  dat <- cbind(melt(pow_avg),th_emp,id,alph,n_study)
  colnames(dat) <- cols
  dat$p1 <- rep(p1s,times=1,each=length(p0s))
  dat$p0 <- rep(p0s,times=length(p1s))
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
  Ts_avg <- matrix(apply(Ts,1,mean,na.rm=TRUE),length(p0s),length(p1s),byrow=TRUE)
  return(Ts_avg)
}

T_sd <- function(Ts){
  Ts_avg <- matrix(apply(Ts,1,sd,na.rm=TRUE),length(p0s),length(p1s),byrow=TRUE)
  return(Ts_avg)
}

#simulate values
#https://stackoverflow.com/questions/34999019/apply-a-function-to-all-pairwise-combinations-of-list-elements-in-r
sims <- vector("list",length(n_studies))
cols_evidence <- c("p0","p1","Tn","Tn_sd","th_emp","id","n_study")
cols_tot <- c("p0","p1","Tn","Tn_sd","th_emp","id","n_study","alpha","CI","power")
evidence_binom <- data.frame(matrix(NA,nrow=0,ncol=12))
colnames(evidence_binom) <- cols_tot
#power_binom <- data.frame(matrix(NA,nrow=0,ncol=6))
#cols_power <- c("p0","p1","power","th_emp","id","alpha","n_study")
#colnames(power_binom) <- cols_power

i <- 1
for (n_study in n_studies) {
  start.time <- Sys.time()
  #question: does transformation make sense in the theoretical case?
  if (T_corr==T){
    p1_hats <- (sapply(p1s,function(x) rbinom(n_sim,n_study,x))+3/8)/(n_study+3/4) #Anscombe estimator (included continuity correction) p_anscombe <- (x+3/8)/(n+3/4)  
    p1s_th <- matrix(rep((p1s*n_study+3/8)/(n_study+3/4),n_sim),n_sim,length(p1s),byrow=TRUE)
  } else {
    p1_hats <- sapply(p1s,function(x) rbinom(n_sim,n_study,x))/n_study #regular estimator: p_hat <- x/n_study
    p1s_th <- matrix(rep(p1s,n_sim),n_sim,length(p1s),byrow=TRUE)
    }
  
  #calculate empirical and theoretical SE
  #se_hats <- SE_calculator(p1_hats,n_study)
  #se_th <- SE_calculator(p1s_th, n_study)
  
  #calculate theoretical and empiral evidence based on clt-vst
  #question: can I use empirical SE to calculate T_clt_th > how is distribution of this statistic?
  T_clt_emp <- calc_T(p0s,p1_hats,n_study,clt)
  T_clt_emp_avg <- T_averager(T_clt_emp)
  T_clt_emp_sd <- T_sd(T_clt_emp)
  
  # if (T_corr==T){
  #   T_clt_th_avg <- sapply((p1s*n_study+3/8)/(n_study+3/4),function(p1) clt(p0s,p1,n_study))
  # } else{
  #   T_clt_th_avg <- sapply(p1s,function(p1) clt(p0s,p1,n_study)) 
  # }
  #T_avg_clt <- dat_transform_th_emp(list(T_clt_th_avg,T_clt_emp_avg))
  #T_clt <- rbind(dat_transform(T_clt_th_avg,"th","clt",n_study,cols_evidence),
  #                   dat_transform(T_clt_emp_avg,"emp","clt",n_study,cols_evidence))
  
  T_clt_th_avg <- sapply(p1s,function(p1) clt(p0s,p1,n_study))
  T_th_sd <- matrix(1,nrow=dim(T_clt_th_avg)[1],ncol=dim(T_clt_th_avg)[2])
  
  T_clt <- rbind(dat_transform(T_clt_emp_avg,T_clt_emp_sd,"emp","clt",n_study,cols_evidence),
                 dat_transform(T_clt_th_avg,T_th_sd,"th","clt",n_study,cols_evidence))
  
  #calculate theoretical and empiral evidence based on binomial vst
  T_vst_emp <- calc_T(p0s,p1_hats,n_study,vst)
  T_vst_emp_avg <- T_averager(T_vst_emp)
  T_vst_emp_sd <- T_sd(T_vst_emp)
  
  # if (T_corr==T){
  #   T_vst_th_avg <- sapply((p1s*n_study+3/8)/(n_study+3/4),function(p1) vst(p0s,p1,n_study))
  # } else{
  #   T_vst_th_avg <- sapply(p1s,function(p1) vst(p0s,p1,n_study))
  # }
  # 
  #T_avg_vst <- dat_transform_th_emp(list(T_vst_th_avg,T_vst_emp_avg))
  #T_avg_vst <- rbind(dat_transform(T_vst_th_avg,"th","vst",n_study,cols_evidence),
  #                   dat_transform(T_vst_emp_avg,"emp","vst",n_study,cols_evidence))
  
  T_vst_th_avg <- sapply(p1s,function(p1) vst(p0s,p1,n_study))
  
  T_vst <- rbind(dat_transform(T_vst_emp_avg,T_vst_emp_sd,"emp","vst",n_study,cols_evidence),
                 dat_transform(T_vst_th_avg,T_th_sd,"th","vst",n_study,cols_evidence))

  #calculate power & confidence intervals
  j <- 1
  sims_pow <- vector("list",length(alphas))
  for (alph in alphas){
    #confidence_intervals
    crit_val_ci <- qnorm((1-alph/2),0,1)
    
    emp_ci_clt <- T_averager(ci_coverage(T_clt_emp,T_clt_emp_avg,crit_val_ci))
    emp_ci_vst <- T_averager(ci_coverage(T_vst_emp,T_vst_emp_avg,crit_val_ci))
    th_ci <- matrix((1-alph),nrow=dim(emp_ci_clt)[1],ncol=dim(emp_ci_clt)[2])
    ci_clt <- rbind(melt(emp_ci_clt),melt(th_ci))$value
    ci_vst <- rbind(melt(emp_ci_vst),melt(th_ci))$value
    T_clt <- cbind(T_clt,alph,ci_clt)
    T_vst <- cbind(T_vst,alph,ci_vst)
    colnames(T_vst) <- c(cols_evidence,"alpha","CI")
    
    #calculate theoretical and empirical power for test using the exact binomial test
    crit_vals_binom <- sapply(p0s,function(p0) qbinom(alph,n_study,p0,lower.tail=FALSE))
    pows_avg_binom_th <- sapply(1:dim(p1s_th)[2], function(i) apply(sapply(p1s_th[,i], function(p1) 1-pbinom(crit_vals_binom,n_study,p1,lower.tail = TRUE)),1,mean))
    
    #the following is wrong! > reason: it mixes theoretical distribution with indivdual simulations
    #pows_avg_binom_emp <- sapply(1:dim(p1s_th)[2], function(i) apply(sapply(p1_hats[,i], function(p1) 1-pbinom(crit_vals,n_study,p1,lower.tail = TRUE)),1,mean))
    #UPPER TAIL TRUE, DANN: PBINOM(...) 

    pows_avg_binom_emp <- t(sapply(crit_vals_binom, 
                                 function(crit_val) apply(test_func(p1_hats*n_study,crit_val),2,mean)))
    
    #pows_avg_binom <- dat_transform_th_emp(list(pows_avg_binom_th,pows_avg_binom_emp))
    #pows_avg_binom <- rbind(dat_transform_power(pows_avg_binom_emp,"emp","binom",alph,n_study,cols_power),
                            #dat_transform_power(pows_avg_binom_th,"th","binom",alph,n_study,cols_power))$power
    pows_avg_binom <- rbind(melt(pows_avg_binom_emp),melt(pows_avg_binom_th))$value
    
    #calculate theoretical and empirical power for the test using the clt (z-score)
    crit_val_norm <- qnorm(alph,mean=0,sd=1,lower.tail = FALSE)
    pows_avg_clt_th <- sapply(1:dim(p1s_th)[2], 
                              function(i) sapply(p0s, 
                                                 function(p0) mean(sapply(p1s_th[,i],
                                                                          function(p1) 1-pnorm(crit_val_norm,clt(p0,p1,n_study),sd=1)))))
    pows_avg_clt_emp <- T_averager(1*(T_clt_emp>crit_val_norm))
    
    #pows_avg_clt <- dat_transform_th_emp(list(pows_avg_clt_th,pows_avg_clt_emp))
    
    #pows_avg_clt <- rbind(dat_transform_power(pows_avg_clt_emp,"emp","clt",alph,n_study,cols_power),
    #                      dat_transform_power(pows_avg_clt_th,"th","clt",alph,n_study,cols_power))$power
    
    pows_avg_clt <- rbind(melt(pows_avg_clt_emp),melt(pows_avg_clt_th))$value
    #calculate theoretical and empirical power for test using the vst
    pows_avg_vst_th <- sapply(1:dim(p1s_th)[2], 
                           function(i) sapply(p0s, 
                                              function(p0) mean(sapply(p1s_th[,i],
                                                                       function(p1) 1-pnorm(crit_val_norm,vst(p0,p1,n_study),sd=1)))))
    pows_avg_vst_emp <- T_averager(1*(T_vst_emp>crit_val_norm))
  
    #pows_avg_vst <- dat_transform_th_emp(list(pows_avg_vst_th,pows_avg_vst_emp))
    #pows_avg_vst <- rbind(dat_transform_power(pows_avg_vst_emp,"emp","vst",alph,n_study,cols_power),
    #                      dat_transform_power(pows_avg_vst_th,"th","vst",alph,n_study,cols_power))$power
    
    pows_avg_vst <- rbind(melt(pows_avg_vst_emp),melt(pows_avg_vst_th))$value
    #maybe replace sd_th with sd_emp?
    #question: in this case, isn't the maximum evidence 1?
    #question: what standard deviation do we use? the empirical or the theoretical?
    #question: how do we know that approximate standard normality holds for vst?
    
    #sims_pow[[j]] <- list(pows_avg_binom,pows_avg_clt,pows_avg_vst)
    j <- j+1
    
    T_clt <- cbind(T_clt,pows_avg_clt)
    pow_binom <- cbind(T_clt$p0,T_clt$p1,NA,NA,T_clt$th_emp,"binom",n_study,alph,NA,pows_avg_binom)
    colnames(pow_binom) <- cols_tot
    colnames(T_clt) <- cols_tot
    T_vst <- cbind(T_vst,pows_avg_vst)
    colnames(T_vst) <- cols_tot
    evidence_binom <- rbind(evidence_binom,T_clt,T_vst,pow_binom)
  }
  
  #sims[[i]] <- list(T_avg_clt,T_avg_vst,T_avg_std_norm,sims_pow)
  #power_binom <- rbind(power_binom,pows_avg_binom,pows_avg_clt,pows_avg_vst)
  
  diff.time <- Sys.time()-start.time
  print(paste0("round ",i," is done in ",round(diff.time,2), " ",attributes(diff.time)$units))
  i <- i+1
}

#save data so that it doesn't have to be computed again and again:
#save(sims, file = paste0("data/",name,".RData"))
save(evidence_binom, file = paste0("data/",name,".RData"))


#question/todo: a lot of 0 values are created empirical SE and low n_study > leads to infinity values for the z-scores based on the empirical variance.
#can be amended by the continuity correction.

# plot values
#plot
#for math notation, see here: https://www.r-bloggers.com/math-notation-for-r-plot-titles-expression-and-bquote/

#NamingConvention: Ev_Binom_corr_1e4.pdf; EV_DIST_CORR_NUMSIM.pdf

pdf(paste0("./figs/",name,".pdf"),onefile=TRUE)
i <- 1
for (sim in sims){
  T_clt <- sim[[1]]
  T_vst <- sim[[2]]
  T_std_norm <- sim[[3]]
  T_std_norm$Tn <- T_std_norm$Tn
  
  Ts <- merger(dats =list(T_clt,T_vst,T_std_norm),ids=c("T_clt","T_vst","T_std_norm"))
  Ts <- Ts[Ts$th_emp=="th",]
  
  T_plot <- ggplot(data=Ts,aes(x=p1,y=Tn/sqrt(n_studies[i]),col=id,group=interaction(p0,th_emp,id),linetype=th_emp)) + 
    geom_line() + scale_color_manual(values=c("T_clt"=2,"T_vst"=3,"T_std_norm"=4))+
    geom_point(data=Ts[Ts$th_emp=="emp",], aes(shape=factor(p0)),alpha=0.5)+
    ylim(-3,3) + xlim(0,1) + ggtitle(bquote("Study size"==.(n_studies[i]))) + labs(x= expression(p[1]), y = "Evidence")
  print(T_plot)
  
  # T_plot <- ggplot(data=Ts,aes(x=p1,y=Tn/sqrt(n_studies[i]),col=factor(p0),group=interaction(p0,th_emp,id),linetype=id,alpha=th_emp)) + 
  #   geom_line() + scale_linetype_manual(values=c("T_clt"=2,"T_vst"=1)) + scale_alpha_manual(values=c("emp"=1,"th"=0.5))+
  #   ylim(-3,3) + xlim(0,1) + ggtitle(bquote("Study size"==.(n_studies[i]))) + labs(x= expression(p[1]), y = "vst (solid), clt (dashed)")
  # print(T_plot)

  #plot the power
  pows <- sim[[4]]
  j <-1
  for (pow in pows){
    pow_bin <- pow[[1]]
    pow_clt <- pow[[2]]
    pow_vst <- pow[[3]]

    pows_merged <- merger(dats=list(pow_bin,pow_clt,pow_vst),ids=c("pow_bin","pow_clt","pow_vst"))
    #pows_merged <- pows_merged[pows_merged$th_emp=="th",]
    
    #pows_merged <- pows_merged[pows_merged$th_emp=="emp" & pows_merged$id =="pow_bin" & pows_merged$p0 == "0.7",]
    
    cutoff <- data.frame(x = c(0, 1), y = alphas[j], cutoff = factor(alphas[j]))
    p0_cutoff <- data.frame(x = rep(p0s,each=2),y=rep(c(0,1),length(p0s)))
    
    pow_plot <- ggplot(data=pows_merged,aes(x=p1,y=Tn,col=id,group=interaction(p0,th_emp,id),linetype=th_emp)) + 
      geom_line() + scale_color_manual(values=c("pow_clt"=2,"pow_vst"=3,"pow_bin"=4))+ 
      geom_point(data=pows_merged[pows_merged$th_emp=="emp",], aes(shape=factor(p0)),alpha=.5) +
      geom_line(aes(x, y,group=factor(x)),p0_cutoff,inherit.aes=FALSE,alpha=.2) + geom_line(aes(x, y),cutoff,inherit.aes=FALSE,alpha=.2) +
      ggtitle(bquote("Study size"==.(n_studies[i])~", alpha"==.(alphas[j]))) + labs(x= expression(mu[1]), y = "power")
    print(pow_plot)
    
    # pow_plot <- ggplot(data=pows_merged,aes(x=p1,y=Tn,col=factor(p0),group=interaction(p0,th_emp,id),linetype=id,alpha=th_emp)) + 
    #   geom_line() + scale_linetype_manual(values=c("pow_bin"=2,"pow_vst"=1)) +  scale_alpha_manual(values=c("emp"=1,"th"=0.5))+
    #   geom_line(aes(x, y,group=factor(x)),p0_cutoff,inherit.aes=FALSE,alpha=.2) + geom_line(aes(x, y),cutoff,inherit.aes=FALSE,alpha=.2) +
    #   ggtitle(bquote("Study size"==.(n_studies[i])~", alpha"==.(alphas[j]))) + labs(x= expression(p[1]), y = "power")
    # print(pow_plot)
    j <- j+1
  }
  
  i <- i+1
}
dev.off()
#questions: 
# 1) where do the gaps between 0.1 and 0 as well as between 0.9 and 1 come from?
# 2) why do vst and clt diverge so much at these values?
# 3) why does exact binomial test have lower power for p0 and p1 large, but higher power for p0 low and p1 large?

# Calculate p-values and Type II errors for the different statistics
# https://www.cyclismo.org/tutorial/R/power.html

#discontinued functions ----
dat_transform_th_emp <- function(dats,th_emps=c("th", "emp")){
  for (i in 1:length(th_emps)){
    dat <- dat_transform(dats[[i]])
    th_emp <- rep(th_emps[i],dim(dat)[1])
    dat$th_emp <- th_emp
    dats[[i]] <- dat
  }
  dats <- do.call("rbind",dats)
  return(dats)
}

merger <- function(dats,ids){
  for (i in 1:length(ids)){
    id <- rep(ids[i],dim(dats[[i]])[1])
    dat <- dats[[i]]
    dat$id <- id
    dats[[i]] <- dat
  }
  dats <- do.call("rbind",dats)
  return(dats)
}

#calcualtes SE of binomial variable X/n
SE_calculator <- function(p,n){
  SE <- sqrt(p*(1-p)/n)
  return(SE)
}


calc_T_clt <- function(p0s,p1s,n_study,func){
  Tn <- sapply(1:dim(p1s)[1], function(p1) sapply(p0s, function(p0) func(p0,p1s[p1,],n_study,SEs[p1,])))
  Tn[is.infinite(Tn)] <- NA
  return(Tn)
}
