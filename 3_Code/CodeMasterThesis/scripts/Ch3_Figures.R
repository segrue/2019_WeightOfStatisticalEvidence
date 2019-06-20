### import packages ------------------------------------------------------------

# for plotting
require("ggplot2") # http://www.sthda.com/english/wiki/ggplot2-essentials
# http://r-statistics.co/ggplot2-Tutorial-With-R.html
# require("gridExtra")
# require("grid")
require("ggpubr")
require("extrafont") # https://cran.r-project.org/web/packages/extrafont/README.html
# I used the Palatino font for the graphic - you probably need to install it if you're using Linux
# font_import("~/.local/share/fonts") # only run the first time / when new fonts were installed
# loadfonts() # only run the first time / when new fonts were installed

# for data handling
require("data.table")
# require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# custom functions
source("./functions/helper_functions.R")

### define input paths (for data) and outpout paths (for figures)---------------
in_path <- "data/"
out_path <- "figs/chapter3/"

### properties that should be the same across all figures ---------------
font_size <- 16
font_family <- "Palatino Linotype"
A4 <- c(8.27, 11.69) # width and height of A4 page in inches
A5 <- c(5.8, 8.3) # width and height of an A5 page
corr_student <- "MLE" # either "MLE" for no correction or "Corr" for finite sample correction
n_sim <- 1000

n_studies <- c(5, 10, 20, 30, 40, 50, 100, 200, 500, 1000) # number of subjects per study

mu0s <- 0
mu1s <- seq(-2, 2, by = 0.1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)

data_name <- paste0("DF_Ev_Mean_", corr_student, "_", n_sim, "_20190514.RData") # name to load and save data etc.
data_path <- paste0(in_path, data_name)

load(data_path, verbose = TRUE) # load data set containing simulated values
dat_orig <- evidence_df


n_study_min <- 5
n_study_max <- 100
x_lim <- 5
y_lim <- 8
rng <- range(c(n_study_min, n_study_max))
ctgs <- unique(dat_orig$n_study)[unique(dat_orig$n_study) %between% rng]
sampling_seed <- 20190612

### Figure 3.1: Funnel plots explained -----------------------------------------
### (fig:funnel_plot)
figname <- "ch3_fig1_funnel_plot"

funnel_plotter <- function(dat, xlim = 2, ylim = 8, ctgs) {
  # define colour code & legend that is consistent over all plots with "ctgs" (regardless of whether they appear in plot or not)
  n_ctgs <- length(ctgs)
  clrs <- brewer.pal(n_ctgs, "Paired")
  names(clrs) <- ctgs
  sc_col <- scale_color_manual(values = clrs, drop = FALSE)
  
  std_errs_inv <- seq(0.1, ylim, by = 0.1)
  x_cutoff <- qnorm(0.05, mean = 0, lower.tail = F) / std_errs_inv
  cutoff <- data.table(x_cutoff, std_errs_inv)
  cutoff_norm <- data.table(x_cutoff = rep(qnorm(0.05, mean = 0, lower.tail = F), length(std_errs_inv)), std_errs_inv)
  funnel_plot <- ggplot(data = dat, 
                        aes(x = mu1_hat, y = 1 / (sgm_hat / sqrt(n_study)), 
                            col = factor(n_study))) + 
    geom_point() + 
    coord_cartesian(xlim = c(-xlim, xlim), ylim = c(0,ylim)) +
    guides(color = guide_legend(reverse = TRUE)) + 
    labs(x = expression(bar(X)), y = expression(1 / SE(bar(X))), color = "n:") +
    theme(
      text = element_text(size = font_size, family = font_family),
      plot.title = element_text(size = font_size)
    ) + 
    sc_col + ylim(0, ylim) + 
    geom_line(data = cutoff, aes(x = x_cutoff, y = std_errs_inv), inherit.aes = F, linetype = "dotted") +
    geom_line(aes(x = mu1, y = seq(0, ylim, length.out = length(mu1))), inherit.aes = F, linetype = "solid")
}

### simulate publication bias
## Situtation 1A:
n_select <- 200
dat_0 <-  evidence_df[mu1 == 0 & n_study %between% rng, ]
dat_0 <- select_studies(dat_0,probs=FALSE,n_select,seed=sampling_seed,evd_id="Vn")
dat_0.3 <- evidence_df[mu1 == 0.3 & n_study %between% rng, ]
dat_0.3 <- select_studies(dat_0.3,probs=FALSE,n_select,seed=sampling_seed,evd_id="Vn")

## Situtation 1B:
## only keep studies which turned out to be significant
dat_0_sig <- dat_0[H1 == 1,]
dat_0.3_sig <- dat_0.3[H1 == 1]

## Situtation 1C:
## only keep studies which turned out to be significant and small percentage of 
## non-significant studies
sel_prob <- 0.1 # percentage of non-significant studies to keep
#probs_mix <- c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1)
dat_0_mix <- rbind(dat_0_sig, 
                 select_studies(dat_0[H1 == 0, ], sel_prob, 
                                n_select = NULL, sampling_seed, evd_id = "Vn"))
dat_0.3_mix <- rbind(dat_0.3_sig, 
                   select_studies(dat_0.3[H1 == 0, ], sel_prob, 
                                  n_select = NULL, sampling_seed, evd_id = "Vn"))


p1 <- funnel_plotter(dat_0, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0))
p2 <- funnel_plotter(dat_0.3, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0.3))
p3 <- funnel_plotter(dat_0_sig, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0))
p4 <- funnel_plotter(dat_0.3_sig, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0.3))
p5 <- funnel_plotter(dat_0_mix, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0))
p6 <- funnel_plotter(dat_0.3_mix, xlim =2, ylim = 8, ctgs = ctgs) + ggtitle(bquote(mu[1]==0.3))

fig1 <- ggarrange(p1, p2, p3, p4, p5, p6,
                  ncol = 2, nrow = 3,
                  align = "hv", legend = "top", common.legend = T,
                  labels = c("A", "", "B", "", "C")
)

ggsave(
  filename = paste0(out_path, figname, ".pdf"), plot = fig1,
  width = A4[1], height = 0.85 * A4[2], device = cairo_pdf
)

### Table 3.1: Filedrawer calculation for data shown above -----------------------------------------
### (tab:file_drawer)
tablename <- "ch3_tab1_filedrawer"

dat_nams <- c("dat_0","dat_0.3","dat_0_sig","dat_0.3_sig","dat_0_mix","dat_0.3_mix")
dats <- list(dat_0,dat_0.3,dat_0_sig,dat_0.3_sig,dat_0_mix,dat_0.3_mix)

for (i in 1:length(dats)){
  print(dat_nams[i])
  calc_filedrawer(dats[[i]])
  cat("\n")
}

for (i in 1:length(dats)){
  print(dat_nams[i])
  calc_filedrawer(dats[[i]], exact = T)
  cat("\n")
}

### Table 3.2: Expected significance -----------------------------------------
### (tab:expected_significance)
expected_significance <- function(dat){
  k <- dim(dat)[1]
  print(paste0("Number of studies: ",k))
  glob_mean <- weighted_mean_var(dat)
  print(paste0("Global mean: ", glob_mean))
  crit_val_norm <- qnorm(0.95,0,1,lower.tail = T)
  E <- sum(pnorm(glob_mean/(dat$sgm_th)*sqrt(dat$n_study)-crit_val_norm,0,1))
  O <- sum(dat$H1)
  A <- ((O-E)^2/E + (O-E)^2/(k-E))
  print(paste0("A: ", A))
}


for (i in 1:length(dats)){
  print(dat_nams[i])
  expected_significance(dats[[i]])
  cat("\n")
}


### Table 3.3: Rank correlation for funnel plot described above -----------------------------------------
### (tab:rank_correlation)

weighted_mean_var <- function(dat) {
  weights <- dat$n_study / dat$sgm_hat^2
  agg_mean <- sum(dat$mu1_hat * weights) / sum(weights)
  return(agg_mean)
}

rank_corr_stat <- function(dat) {
  k <- dim(dat)[1]
  print(paste0("Number of studies: ",k))
  weights <- dat$n_study / dat$sgm_hat^2
  glob_mean <- weighted_mean_var(dat)
  print(paste0("Global mean: ", glob_mean))
  Var_js <- weights - 1/sum(weights)
  Z_js <- (dat$mu1_hat-glob_mean)/sqrt(Var_js)
  diff_Z_js <- outer(Z_js,Z_js,"-")
  diff_Var_js <- outer(Var_js,Var_js,"-")
  extract_idx <- upper.tri(diff_Var_js, diag = FALSE)
  S_k <- sum(sign(diff_Var_js[extract_idx])*sign(diff_Z_js[extract_idx])) / sqrt((2*k+5)*k*(k-1)/18)
  print(paste0("S_k = ", S_k))
  return(S_k)
}

for (i in 1:length(dats)){
  print(dat_nams[i])
  rank_corr_stat(dats[[i]])
  cat("\n")
}

### Table 3.4: Egger regression for funnel plot described above -----------------------------------------
### (tab:Egger_regression)
Egger_regression <- function(dat) {
  k <- dim(dat)[1]
  print(paste0("Number of studies: ",k))
  weights <- dat$n_study / dat$sgm_hat^2
  glob_mean <- weighted_mean_var(dat)
  print(paste0("Global mean: ", glob_mean))
  Var_js <- weights - 1/sum(weights)
  Z_js <- (dat$mu1_hat-glob_mean)/sqrt(Var_js)
  prec <- 1/sqrt(Var_js)
  m1 <- lm(Z_js ~ 1+prec)
  print(summary(m1)$coefficients)
}

for (i in 1:length(dats)){
  print(dat_nams[i])
  Egger_regression(dats[[i]])
  cat("\n")
}

### Table 3.5: Caliper test for funnel plot described above -----------------------------------------
### (tab:caliper_test)
caliper <- function(dat,e= 0.2) {
  k <- dim(dat)[1]
  print(paste0("Number of studies: ",k))
  weights <- dat$n_study / dat$sgm_hat^2
  glob_mean <- weighted_mean_var(dat)
  print(paste0("Global mean: ", glob_mean))
  crit_val_norm <- qnorm(0.95,0,1,lower.tail = T)
  upper <- crit_val_norm + e
  lower <- crit_val_norm - e
  dat_selected <- dat[Tn < upper & Tn > lower,]
  k_prime <- dim(dat_selected)[1]
  print(paste0("Number of studies around trheshold: ",k_prime))
  dat_sig <- dat_selected[H1==1,]
  k_2prime <-  dim(dat_sig)[1]
  print(paste0("Number of significant studies around threshold: ",k_2prime))
  crit_val_bin <- qbinom(0.975,k_prime,0.5)
  print(paste0("Caliper threshold value ",crit_val_bin))
}

for (i in 1:length(dats)){
  print(dat_nams[i])
  caliper(dats[[i]])
  cat("\n")
}



### Table 3.6: The p-curve -----------------------------------------
### (tab:p_curve)
p_curve <- function(dat,under_null =T) {
  k <- dim(dat)[1]
  print(paste0("Number of studies: ",k))
  weights <- dat$n_study / dat$sgm_hat^2
  glob_mean <- weighted_mean_var(dat)
  print(paste0("Global mean: ", glob_mean))
  dat_sig <- dat[H1==1,]
  l <- dim(dat_sig)[1]
  skew <- 1/l*sum((dat_sig$p-mean(dat_sig$p))^3)/(1/(l-1)*sum((dat_sig$p-mean(dat_sig$p))^2))^(3/2)
  print(paste0("skew: ",skew))
  if (under_null == T) {
    pp_val <- dat_sig$p/0.05
  }
  else {
    pp_val <- dat_sig$p/0.05
  }
  print(paste0("l: ", l))
  pp_comb <- -2*sum(log(pp_val))
  print(paste0("pp_comb: ", pp_comb))
  crit_val <- qchisq(0.95,2*l)
  print(paste0("Critical value: ", crit_val))
}

for (i in 1:length(dats)){
  print(dat_nams[i])
  p_curve(dats[[i]])
  cat("\n")
}

