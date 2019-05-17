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

### load required data sets ----
load(paste0("data/",name,".RData"),verbose=TRUE) #load data set containing simulated values

### simulate publication bias
## Situtation 1:
## no publication bias - all studies (n=100) are included
clt <- evidence_df[mu1=="0.2" & id == "clt",]
vst <- evidence_df[mu1=="0.2" & id == "vst",]
stud <- evidence_df[mu1=="0.2" & id == "stud",]
p1 <- ggplot(data=clt,aes(x=Tn/sqrt(dat$n_study),y=n_study,col=factor(n_study))) + geom_point() + xlim(-3,3)
p2 <- ggplot(data=clt,aes(x=Tn/sqrt(dat$n_study),y=1/(sgm_hat/sqrt(dat$n_study)),col=factor(n_study))) + geom_point()+ xlim(-3,3)
p3 <- ggplot(data=stud,aes(x=Tn/sqrt(dat$n_study),y=n_study,col=factor(n_study))) + geom_point() + xlim(-3,3)
p4 <- ggplot(data=stud,aes(x=Tn/sqrt(dat$n_study),y=1/(sgm_hat/sqrt(dat$n_study)),col=factor(n_study))) + geom_point()+ xlim(-3,3)
p5 <- ggplot(data=vst,aes(x=Tn/sqrt(dat$n_study),y=n_study,col=factor(n_study))) + geom_point() + xlim(-3,3)
p6 <- ggplot(data=vst,aes(x=Tn/sqrt(dat$n_study),y=1/(sgm_hat/sqrt(dat$n_study)),col=factor(n_study))) + geom_point()+ xlim(-3,3)
p7 <- ggplot(data=clt,aes(x=mu1_hat,y=n_study,col=factor(n_study))) + geom_point()+ xlim(-3,3)
p8 <- ggplot(data=clt,aes(x=mu1_hat,y=1/(sgm_hat/sqrt(dat$n_study)),col=factor(n_study))) + geom_point()+ xlim(-3,3)

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,nrow=4,ncol=2)

## Situtation 1:
## only keep studies which turned out to be significant
clt_sig <- clt[]
  
geom_line() + scale_color_manual(values=c("clt"=3,"vst"=4,"std_norm"=2))+
  ylim(-3,3) + xlim(0,1) + ggtitle(bquote("Study size"==.(dat$n_study[1]))) + labs(x= expression(p[1]), y = "Evidence/sqrt(n)")
return(T_plot)



