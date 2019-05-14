#simulation of binomials to asses normality
#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# Define starting assumptions and conditions ----
# hypotheses:
# H0: p = p0
# H1: p > p0
# Note: all calculations be
n_studies <- seq(5,1000)
n_sim <- 5e3
n_sim_rss <- 1000
alpha <- 0.001
p1s <- seq(0.01,0.99,by=0.01)
p0 <- 0
probs <- seq(0.01,0.99,length.out=100) #probs for which to calculate quantiles
seed <- 20190505
set.seed(20190505)

# Define vst and helper functions  ----
vst <- function(p0,p1,n){
  Tn <- 2*sqrt(n)*(asin(sqrt(p1))-asin(sqrt(p0)))
  return(Tn)
}

#z-statistic based on CLT
clt <- function(p0,p1,n){
  Tn <- (p1-p0)/sqrt(p1*(1-p1)/n)
  return(Tn)
}

#brings data into the correct form for plottings
dat_transform <- function(T_avg,id,n_study){
  dat <- cbind(melt(T_avg),id,n_study)
  colnames(dat) <- cols
  dat$prob <- rep(probs,length(p1s))
  dat$p1 <- rep(p1s,each=length(probs))
  return(dat)
}


#calculates rss
rss_calculator <- function(x,y,avg=T){
  if (avg==T){
    rss <- sum((x-y)^2)/length(x)
  } else {
    rss <- sum((x-y)^2) 
  }
  return(rss)
}

#calculate and save theoretical distributions T_clt and T_vst
p1s <- c(0.1,0.5,0.9)
p0 <- 0
n_studies <- c(5,10,20,30)
probs <- seq(0.01,0.99,length.out=100)
n_sim <- 5e3
cols <- c("prob","p1","quantile","id","n_study")
seed <- 20190505
corr <- T
correction <- ifelse(corr,"Ans","MLE")
name <- paste0("quantiles_binom_",correction,"_",n_sim,"_",seed)

quantiles <- data.frame(matrix(NA,nrow=0,ncol=5))
colnames(quantiles) <- cols
#need to adjust for the expectation of T_Vst!! E(Tn_vst) = sqrt(n)*K(theta)

set.seed(seed)
for (i in 1:length(n_studies)){
  idx <- 1:length(p1s)
  n_study <- n_studies[i]
  ses <- sqrt(p1s*(1-p1s)/n_study)
  mus_clt <- (p1s-p0)/ses
  mus_vst <- vst(p0,p1s,n_study)
  diffs <- mus_clt-mus_vst
    
  Qnorms_th <- dat_transform((sapply(mus_clt,function(mu) qnorm(probs,mu,1))),"th",n_study) #theoretical quantiles of normal distribution with mean mu
  
  #calculate quantiles of T_vst and T_clt
  if (corr==T){
    p1_hats <- (sapply(p1s,function(x) rbinom(n_sim,n_study,x))+3/8)/(n_study+3/4) #Anscombe estimator (included continuity correction) p_anscombe <- (x+3/8)/(n+3/4)  
  } else {
    p1_hats <- sapply(p1s, function(p1) rbinom(n_sim,n_study,p1))/n_study 
  }
  T_vst <- sapply(idx, function(j) vst(p0,p1_hats[,j],n_study)+diffs[j]) #difference needs to be added so that T_vst and T_clt have the same expectation
  T_clt <- sapply(idx, function(j) clt(p0,p1_hats[,j],n_study))
  T_clt[is.infinite(T_clt)] <- NA
  
  T_clt_q <- dat_transform((sapply(idx, function(j) quantile(T_clt[,j],probs,na.rm=T))),"clt",n_study)
  T_vst_q <- dat_transform((sapply(idx, function(j) quantile(T_vst[,j],probs,na.rm=T))),"vst",n_study)
  quantiles <- rbind(quantiles,Qnorms_th,T_clt_q,T_vst_q)
}

save(quantiles,file=paste0("data/",name,".RData"))

#simulate values and calculate RSS
#Qnorm_th <- qnorm(probs,0,1) #theoretical quantiles of standard normal
Qnorms_th <- sapply(p1s,function(p1) qnorm(probs,p1,1))#theoretical quantiles of normal distribution with mean p1
set.seed(seed)
rss_norm_means <- sapply(1:length(p1s), function(i) mean(sapply(1:n_sim_rss, function(j) rss_calculator(quantile(rnorm(n_sim,p1s[i],1),probs),Qnorms_th[,i],avg=T))))
set.seed(seed)
rss_norm_sds <- sapply(1:length(p1s), function(i) sd(sapply(1:n_sim_rss, function(j) rss_calculator(quantile(rnorm(n_sim,p1s[i],1),probs),Qnorms_th[,i],avg=T))))

rss_vals <- data.frame(matrix(NA,nrow=length(p1s)*2*length(n_studies),ncol=4))
cols <- c("rss","Tn","p1","n_study")
colnames(rss_vals) <- cols

for (i in 1:length(n_studies)){
  start.time <- Sys.time()
  n_study <- n_studies[i]
  set.seed(seed)
  binoms <- sapply(p1s, function(p1) rbinom(n_sim,n_study,p1))/n_study
  T_vst <- sapply(1:length(p1s), function(i) vst(p0,binoms[,i],n_study))
  T_clt <- sapply(1:length(p1s), function(i) clt(p0,binoms[,i],n_study))
  T_clt[is.infinite(T_clt)] <- NA
  
  rss_clt <- data.frame(sapply(1:length(p1s), function(i) rss_calculator(quantile(T_clt[,i],probs,na.rm=T),Qnorms_th[,i],avg=T)))
  rss_vst <- data.frame(sapply(1:length(p1s), function(i) rss_calculator(quantile(T_vst[,i],probs,na.rm=T),Qnorms_th[,i],avg=T)))
  rss_clt <- cbind(rss_clt,"clt",p1s,n_study)
  rss_vst <- cbind(rss_vst,"vst",p1s,n_study)
  colnames(rss_clt) <- colnames(rss_vst) <- cols
  idx <- (dim(rss_vst)[1]*2*(i-1))+1:(dim(rss_vst)[1]*2*i)
  rss_vals[idx,] <- rbind(rss_clt,rss_vst)
  
  diff.time <- Sys.time()-start.time
  print(paste0("round ",i," is done in ",round(diff.time,2), " ",attributes(diff.time)$units))
}

save(rss_vals, file = paste0("data/rss_values.RData"))
#next steps: compare rss_vals with the rss_vals from the theoretical distribution and reject all values which are x standard deviations beyond the mean
#problem: rss based on T_Clt and T_vst is always much larger than rss based on rnorm > maybe use T_clt(rnorm) as basis?



shapiro.test(rbinom(200,100,0.5))

test <- rbinom(n_sim,5,0.5)
test <- test/5

shapiro.test(vst(0,rbinom(n_sim,100,0.5)/100,100))

n_p_mat <- matrix(data=c(n_studies,p1s),ncol=length(p1s),nrow=length(n_studies))

T_vst <- sapply(p1s,function(p1) sapply(n_studies, function(n_study) vst(0,p1,n_study)))
T_clt <- sapply(p1s,function(p1) sapply(n_studies, function(n_study) clt(0,p1,n_study)))

n_sim_rss <- 100
p <- 0.5
mu <- p*n
sd <- sqrt(n*p*(1-p))
se <- sqrt(p*(1-p)/n)

rss_df <- data.frame(matrix(NA,nrow=n_sim_rss,ncol=2))
colnames(rss_df) <- c("RSS_bino", "RSS_norm")


test <- t(sapply(1:n_sim_rss, function(i) q_calculator(100,0.5,n_sim,probs)))

rss_df[,] <- test
mean(rss_df$RSS_bino)
mean(rss_df$RSS_norm)
sd(rss_df$RSS_norm)

z_bino <- (rbinom(n_sim,n,p)-mu)/n/se
Fz_bino <- ecdf(z_bino)
Qz_bino <- Fz_bino(probs) 

normdist <- rnorm(n_sim,0,1)
Fnorm_emp <- ecdf(normdist)
Qnorm_emp <- Fnorm_emp(probs)
  
rss_calculator(Qz_bino,Qnorm_th,avg=T)
rss_calculator(Qnorm_emp,Qnorm_th,avg=T)


