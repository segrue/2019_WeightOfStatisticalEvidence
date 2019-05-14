# Simulate p-values - Binomial case ----

# set up initialisation values ----
require("RColorBrewer")
n <- 100
p0 <- 0.5
ps <- seq(0.0,1,by=0.1)
cols <- brewer.pal(length(ps), "Set3")
sims <- 10000
n_bins <- 10
hist_breaks=seq(0,1,l=n_bins+1)
mid_points <- seq(1/n_bins/2,1-(1/n_bins/2),by=1/n_bins)

# the easy-peasy way using the "TeachingDemos" package
require("TeachingDemos")
require("likelihood")

helper_pval <- function(p){
  pval <- Pvalue.binom.sim(n=n,p=p,p0=p0,test="exact",alternative="greater",
                         alpha=0.05,B=sims)
  #p_dens <- density(pval,bw=0.01,kernel="gaussian",weights=rep(1/sims,sims),from=0,to=1,n=hist_bins+1,
                    #cut=0,adjust=1)
  #plot(p_dens,ylim=c(0,1))
  pval_density <- hist(pval,breaks=hist_breaks,plot=FALSE)$density/n_bins
  return(pval_density)
}

pvals <- data.frame(sapply(ps,helper_pval))

helper_Bayes <- function(p){
  X <- rbinom(sims,n,p)
  likeli0 <- dbinom(X,n,p0)
  likeli1 <- dbinom(X,n,p)
  bf <- likeli0/likeli1
  bf_density <- hist(bf,breaks=hist_breaks,plot=FALSE)$density/n_bins
  return(bf_density)
}

bfs <- data.frame(sapply(ps,helper_Bayes))

# plot values
require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

plotter <- function(vals,abscissa){
  vals <- cbind(abscissa,vals)
  colnames(vals) <- c("x",as.character(ps))
  dat <- melt(vals,id.vars = "x")
  ggplot(dat,aes(x=x,y=value,color=variable))+geom_line()+scale_color_manual(values=cols)
}

plotter(pvals,mid_points)
plotter(bfs,mid_points)

# the more tedious, but also more insightful way


# Simulate p-values - approximate normal case ----
mu0 = p0*n
mus <- ps*n
sigma0 <- sigma <- n*p0*(1-p0)

helper_pval_norm <- function(mu){
  pval <- Pvalue.norm.sim(n=n,mu=mu,mu0=mu0, sigma=sigma,sigma0=sigma0,
                          test="z",alternative="greater",alpha=0.05,B=sims)
  #p_dens <- density(pval,bw=0.01,kernel="gaussian",weights=rep(1/sims,sims),from=0,to=1,n=hist_bins+1,
  #cut=0,adjust=1)
  #plot(p_dens,ylim=c(0,1))
  pval_density <- hist(pval,breaks=hist_breaks,plot=FALSE)$density/n_bins
  return(pval_density)
}

pvals <- data.frame(sapply(mus,helper_pval_norm))
plotter(pvals,mid_points)


# Simulate random p-value - binomial case ----

random_pv <- function(mu){
  F_pv <- 1 - pnorm(qnorm(1-p,n,p0),n,p)
}

helper_random_pv <- function(mu){
  pval <- Pvalue.norm.sim(n=n,mu=mu,mu0=mu0, sigma=sigma,sigma0=sigma0,
                          test="z",alternative="greater",alpha=0.05,B=sims)
  #p_dens <- density(pval,bw=0.01,kernel="gaussian",weights=rep(1/sims,sims),from=0,to=1,n=hist_bins+1,
  #cut=0,adjust=1)
  #plot(p_dens,ylim=c(0,1))
  q_ <- qnorm(1-pval,mean=mu0,sd=sigma0)
  random_pvs <- 1-pnorm(q_,mean=mu,sd=sigma)
  pval_density <- hist(random_pvs,breaks=hist_breaks,plot=FALSE)$density/n_bins
  #pval_density <- dnorm(q_-mu)/dnorm(q_)
  return(pval_density)
}

random_pvals <- data.frame(sapply(mus,helper_random_pv))
plotter(random_pvals,mid_points)
