### import packages ------------------------------------------------------------

#for plotting
require("ggplot2") #http://www.sthda.com/english/wiki/ggplot2-essentials
#http://r-statistics.co/ggplot2-Tutorial-With-R.html
#require("gridExtra")
#require("grid")
require("ggpubr")
require("extrafont") # https://cran.r-project.org/web/packages/extrafont/README.html
# I used the Palatino font for the graphic - you probably need to install it if you're using Linux
#font_import("~/.local/share/fonts") # only run the first time / when new fonts were installed
#loadfonts() # only run the first time / when new fonts were installed

#for data handling
require("data.table")
#require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# custom functions
source("./functions/helper_functions.R")

### define input paths (for data) and outpout paths (for figures)---------------
in_path <- "data/"
out_path <- "figs/chapter2/"

### Figure 1: Type I and Type II error explained (fig:error_types)
figname <- "ch2_fig1_error_types"
x_max = 5
mu0 <- 0
mu1 <- 2
p_alpha <- 0.05

# define x and y ranges
x <- seq(-x_max,x_max,by=0.001)
x_major_breaks <- c(seq(-x_max+1,x_max-1,by=2))
x_minor_breaks <- c(x_major_breaks-1)[-1]
y_major_breaks <- seq(0,1,by=0.25)

q_alpha <- qnorm(p_alpha,mu0,1,lower.tail = FALSE)
p_beta <- pnorm(q_alpha,mean=mu1,1)

# plot density function

p1 <- ggplot(data=data.table(x),aes(x=x)) + 
      stat_function(fun=dnorm,args = list(mean = mu0),aes(linetype="1")) + 
      stat_function(fun=dnorm, xlim=c(q_alpha,x_max), geom="area", 
                    args = list(mean = mu0),alpha=0.5,aes(fill="1")) +
      stat_function(fun=dnorm,args = list(mean = mu1),aes(linetype="2")) + 
      stat_function(fun=dnorm, xlim=c(-x_max,q_alpha), geom="area", 
                args = list(mean = mu1),alpha=0.5,aes(fill="2")) +
      scale_linetype(name = "Mean:", 
                     labels = c(bquote(mu==.(mu0)),bquote(mu==.(mu1)))) +
      scale_fill_manual(name = "Error:", 
                        labels = c(bquote(alpha==.(p_alpha)),
                                   bquote(beta==.(round(p_beta,2)))),
                        values = c("1"=alpha("red",0.5),"2"=alpha("blue",0.5))) +
      scale_x_continuous(breaks = c(x_major_breaks, q_alpha), 
                         labels = c(x_major_breaks, bquote(z[(1-alpha)])),
                         minor_breaks = x_minor_breaks) +
      theme(panel.grid.major.x = element_line(color = c(rep("white",length(x_major_breaks)),NA)),
            panel.grid.minor.y = element_blank(),
            text=element_text(size=16, family="Palatino Linotype")) +
      ylab(label = bquote(phi(x-mu)))#+
      #theme(axis.title.x = element_blank(), axis.text.x = element_blank())
#leg <- get_legend(p1)

# plot cdf
p2 <- ggplot(data=data.table(x),aes(x=x)) + 
  stat_function(fun=pnorm,args = list(mean = mu0),aes(linetype="1")) + 
  stat_function(fun=pnorm,args = list(mean = mu1),aes(linetype="2")) +
  geom_segment(aes(x = -x_max, y = 1-p_alpha, xend = q_alpha, yend=1-p_alpha),
               col="red",size=0.3) +
  geom_segment(aes(x = q_alpha, y = 1-p_alpha, xend = q_alpha, yend = 1),
               col="red",size=0.3) +
  geom_segment(aes(x = -x_max, y = p_beta, xend = q_alpha, yend=p_beta),
               col="blue",size=0.3) +
  geom_segment(aes(x = q_alpha, y = 0, xend = q_alpha, yend = p_beta),
               col="blue",size=0.3) +
  scale_x_continuous(breaks = c(x_major_breaks, q_alpha), 
                     labels = c(x_major_breaks, bquote(z[(1-alpha)])),
                     minor_breaks = x_minor_breaks) +
  theme(panel.grid.major.x = element_line(color = c(rep("white",length(x_major_breaks)),NA)),
        panel.grid.minor.y = element_blank(),
        text=element_text(size=16, family="Palatino Linotype")) +
  ylab(label = bquote(Phi(x-mu)))

#width and height of A4 page in inches: width="8.27", height="11.69"
#width and height of an A5 page: width = "5.8", height = "8.3"
# make sure that all axis are nicely aligned: 
# http://www.sthda.com/english/wiki/print.php?id=177
fig1 <- ggarrange(p1,p2+theme(legend.position="none"),ncol=1,nrow=2,
          align="hv",legend="top",common.legend = T,labels = c("A", "B"))
ggsave(filename= paste0(out_path,figname,".pdf"),plot=fig1,
       width = 8.3, height = 5.8, device = cairo_pdf)

