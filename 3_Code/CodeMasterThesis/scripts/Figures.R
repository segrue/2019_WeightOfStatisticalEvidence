### import packages ------------------------------------------------------------

#for plotting
require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
require("gridExtra")

#for data handling
require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

### Figure fig:quantiles_binom, describing approximation of z-score based on the binomial ----
corr <- "MLE"
load(paste0("data/quantiles_binom_",corr,"_5000_20190505.RData"),verbose=TRUE) #load quantiles

#plotting multiple plots in one window: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
#n=5

plotter <- function(dat){
  T_plot <- ggplot(data=dat,aes(x=quantile,y=prob, col=id, group=interaction(p1,id),linetype=factor(p1))) + 
    geom_line() + scale_color_manual(values=c("th"=2,"clt"=3,"vst"=4))+
    ylim(0,1) + ggtitle(bquote("Study size"==.(dat$n_study[1]))) + labs(x= "Quantiles", y = "Cumulative Probability")
}

#task: plot all four figurs (one for each n_study) on the same page
pdf(paste0("./figs/quantiles_binom_",corr,".pdf"),onefile=TRUE)
p1 <- plotter(quantiles[quantiles$n_study==5,])
p2 <- plotter(quantiles[quantiles$n_study==10,])
p3 <- plotter(quantiles[quantiles$n_study==20,])
p4 <- plotter(quantiles[quantiles$n_study==30,])
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()


### Figure fig:evidence_binom - figure of theoretical evidence using the z-score and the clt based on the Binomial ----
corr <- "Ans"
load(paste0("data/Ev_Binom_",corr,"_100_20190504.RData"),verbose=TRUE)
selection <- ((evidence_binom$p0 %in% c(0.1,0.5,0.9)) & 
                ((evidence_binom$th_emp == "emp" & evidence_binom$id %in% c("clt","vst")) 
              | evidence_binom$id == "std_norm"))
evidence_binom <- evidence_binom[selection,]

plotter <- function(dat){
  T_plot <- ggplot(data=dat,aes(x=p1,y=Tn/sqrt(dat$n_study[1]),col=id,group=interaction(p0,id),linetype=factor(p0))) + 
    geom_line() + scale_color_manual(values=c("clt"=3,"vst"=4,"std_norm"=2))+
    ylim(-3,3) + xlim(0,1) + ggtitle(bquote("Study size"==.(dat$n_study[1]))) + labs(x= expression(p[1]), y = "Evidence/sqrt(n)")
  return(T_plot)
  }

pdf(paste0("./figs/evidence_binom_",corr,".pdf"),onefile=TRUE)
dat <- evidence_binom[(evidence_binom$p0 %in% c(0.1,0.5,0.9)) & evidence_binom$n_study==5 & evidence_binom$th_emp=="emp",]

#plot to plot CI
p1 <- ggplot(data=dat,aes(x=p1,y=Tn/sqrt(dat$n_study[1]),
                          ymin=(Tn-Tn_sd)/sqrt(dat$n_study[1]),ymax=(Tn+Tn_sd)/sqrt(dat$n_study[1]),
                          col=id,group=interaction(p0,id),linetype=factor(p0))) + geom_line() +
      geom_ribbon(aes(fill=id),alpha=0.2,linetype=0) + 
      scale_color_manual(values=c("clt"=3,"vst"=4,"std_norm"=2))+scale_fill_manual(values=c("clt"=3,"vst"=4,"std_norm"=2))+
      ylim(-4,4) + xlim(0,1) + ggtitle(bquote("Study size"==.(dat$n_study[1]))) + labs(x= expression(p[1]), y = "Evidence/sqrt(n)")

print(p1)

p1 <- plotter(evidence_binom[evidence_binom$n_study==5,])
p2 <- plotter(evidence_binom[evidence_binom$n_study==10,])
p3 <- plotter(evidence_binom[evidence_binom$n_study==20,])
p4 <- plotter(evidence_binom[evidence_binom$n_study==30,])
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()
#tochange: need to change how clt is calculated. as of now, infinity values are excluded > this artificially
#distorts the final result
