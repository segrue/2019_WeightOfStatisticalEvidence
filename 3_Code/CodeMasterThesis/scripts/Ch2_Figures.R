### import packages ------------------------------------------------------------

#for plotting
require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-essentials
#http://r-statistics.co/ggplot2-Tutorial-With-R.html
require("gridExtra")
#require("grid")
require("ggpubr")

#for data handling
require("data.table")
#require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# custom functions
source("./functions/helper_functions.R")

### define input paths (for data) and outpout paths (for figures)---------------
in_path <- "./data/"
out_path <- ".figs/chapter2/"

### Figure 1: Type I and Type II error explained (fig:error_types)
x_max = 5
x <- seq(-x_max,x_max,by=0.001)
mu0 <- 0
mu1 <- 2
p_alpha <- 0.05
q_alpha <- qnorm(p_alpha,mu0,1,lower.tail = FALSE)
p_beta <- pnorm(q_alpha,mean=mu1,1)

p1 <- ggplot(data=data.table(x),aes(x=x)) + 
      stat_function(fun=dnorm,args = list(mean = mu0),aes(linetype="1")) + 
      stat_function(fun=dnorm, xlim=c(q_alpha,x_max), geom="area", 
                    args = list(mean = mu0),alpha=0.5,aes(fill="1")) +
      stat_function(fun=dnorm,args = list(mean = mu1),aes(linetype="2")) + 
      stat_function(fun=dnorm, xlim=c(-x_max,q_alpha), geom="area", 
                args = list(mean = mu1),alpha=0.5,aes(fill="2")) +
      scale_linetype(name = "Mean", 
                     labels = c(bquote(mu[0]==~.(mu0)),bquote(mu[1]==~.(mu1)))) +
      scale_fill_manual(name = "Error", 
                        labels = c(bquote(Type~I==~.(p_alpha)),
                                   bquote(Type~II==~.(round(p_beta,2)))),
                        values = c("1"=alpha("red",0.1),"2"=alpha("blue",1))) +
      theme(axis.title.x = element_blank(), axis.text.x = element_blank())
leg <- get_legend(p1)

p2 <- ggplot(data=data.table(x),aes(x=x)) + 
  stat_function(fun=pnorm,args = list(mean = mu0),aes(linetype="1")) + 
  stat_function(fun=pnorm,args = list(mean = mu1),aes(linetype="2")) +
  #geom_vline(xintercept=q_alpha, linetype=3, size=1) +
  geom_segment(aes(x = q_alpha, y = 0, xend = q_alpha, yend = 1),size=0.1) +
  geom_segment(aes(x = -x_max, y = 1-p_alpha, xend = q_alpha, yend=1-p_alpha),
               col="red",size=0.05) +
  geom_segment(aes(x = -x_max, y = p_beta, xend = q_alpha, yend=p_beta),
               col="blue",size=0.05)

# make sure that all axis are nicely aligned: 
# http://www.sthda.com/english/wiki/print.php?id=177
ggarrange(p1,p2+theme(legend.position="none"),ncol=1,nrow=2,
          align="hv",legend="right",common.legend = T)