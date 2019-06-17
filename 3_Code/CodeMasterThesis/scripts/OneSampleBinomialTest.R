### import packages ------------------------------------------------------------

# for plotting
require("ggplot2") # http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

# for data handling
require("reshape2") # https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
require("data.table")

# custom functions
source("functions/helper_functions.R")

### Define starting assumptions and conditions ---------------------------------
# hypotheses:
# H0: p = p0
# H1: p > p0
n_studies <- c(5, 10, 20, 30, 40, 50, 100, 200, 500, 1000)
n_sim <- 1e5
p0s <- c(0.1,0.5,0.9)
p1s <- seq(0.01, 0.99, by = 0.01)
alphas <- c(0.05)
corr <- T
seed <- 20190504
set.seed(seed)
correction <- ifelse(corr, "Ans", "MLE")
name <- paste0("Ev_Binom_", correction, "_", n_sim, "_", seed)

### simulate values ------------------------------------------------------------
cols_evidence <- c("p0", "p1", "evd_mean", "evd_sd", "th_emp", "id", "n_study")
cols_tot <- c("p0", "p1", "evd_mean", "evd_sd", "th_emp", "id", "n_study", "alpha", "CI", "power")
evidence_binom <- data.table(matrix(NA, nrow = 0, ncol = 10))
colnames(evidence_binom) <- cols_tot

i <- 1
for (n_study in n_studies) {
  start.time <- Sys.time()
  # question: does transformation make sense in the theoretical case?
  x_binom <- sapply(p1s, function(x) rbinom(n_sim, n_study, x))
  if (corr == T) {
    # Anscombe estimator (includes continuity correction) p_anscombe <- (x+3/8)/(n+3/4)
    p1_hats <- (x_binom + 3 / 8) / 
               (n_study + 3 / 4)
    #p1s_th <- matrix(rep((p1s * n_study + 3 / 8) / (n_study + 3 / 4), n_sim), 
                     #n_sim, length(p1s), byrow = TRUE)
    p1s_th <- (p1s*n_study + 3 / 8) / (n_study + 3 / 4)
  } else {
    # MLE estimator without continuity correction p_hat <- x/n
    p1_hats <- x_binom / n_study
    p1_hats[p1_hats==0] <- 1e-6 #add small value to prevent infinity values
    p1_hats[p1_hats==1] <- 1-1e-6 #remove small value to prevent infinity values
    #p1s_th <- matrix(rep(p1s, n_sim), n_sim, length(p1s), byrow = TRUE)
    p1s_th <- p1s
  }
  # calculate theoretical and empiral evidence based on z-statistic
  # question: can I use empirical SE to calculate Zn_th > how is distribution of this statistic?
  Zn_emp <- calc_evidence(mu0s = p0s, mu1s = p1_hats, n_study = n_study, func = z_stat_binom)
  Zn_emp_avg <- avg_evidence(Zn_emp, p0s, p1s_th)
  Zn_emp_sd <- sd_evidence(Zn_emp, p0s, p1s_th)

  Zn_th_avg <- sapply(p1s_th, function(p1) z_stat_binom(p0s, p1, n_study))
  Zn_th_sd <- matrix(1, nrow = dim(Zn_th_avg)[1], ncol = dim(Zn_th_avg)[2])

  Zn <- rbind(
    dat_transform(Zn_emp_avg, Zn_emp_sd, "emp", "Zn", 
                  n_study, cols_evidence, p0s, p1s),
    dat_transform(Zn_th_avg, Zn_th_sd, "th", "Zn", 
                  n_study, cols_evidence, p0s, p1s)
  )

  # calculate theoretical and empiral evidence based on binomial vst
  Vn_emp <- calc_evidence(mu0s = p0s, mu1s = p1_hats, n_study = n_study, func = vst_binom)
  Vn_emp_avg <- avg_evidence(Vn_emp, p0s, p1s_th)
  Vn_emp_sd <- sd_evidence(Vn_emp, p0s, p1s_th)

  Vn_th_avg <- sapply(p1s_th, function(p1) vst_binom(p0s, p1, n_study))
  Vn_th_sd <- Zn_th_sd

  Vn <- rbind(
    dat_transform(Vn_emp_avg, Vn_emp_sd, "emp", "Vn", 
                  n_study, cols_evidence, p0s, p1s),
    dat_transform(Vn_th_avg, Vn_th_sd, "th", "Vn", 
                  n_study, cols_evidence, p0s, p1s)
  )

  # calculate power & confidence intervals
  for (alph in alphas) {
    # confidence_intervals
    crit_val_ci <- qnorm((1 - alph / 2), 0, 1)

    emp_ci_Zn <- avg_evidence(ci_coverage(Zn_emp, Zn_th_avg, crit_val_ci), 
                               p0s, p1s_th)
    emp_ci_Vn <- avg_evidence(ci_coverage(Vn_emp, Vn_th_avg, crit_val_ci), 
                               p0s, p1s_th)
    th_ci <- matrix((1 - alph), 
                    nrow = dim(emp_ci_Zn)[1], ncol = dim(emp_ci_Zn)[2])
    ci_Zn <- rbind(melt(emp_ci_Zn), melt(th_ci))$value
    ci_Vn <- rbind(melt(emp_ci_Vn), melt(th_ci))$value
    Zn <- cbind(Zn, alph, ci_Zn)
    Vn <- cbind(Vn, alph, ci_Vn)
    colnames(Vn) <- c(cols_evidence, "alpha", "CI")

    # calculate theoretical and empirical power for test using the exact binomial test
    crit_vals_binom <- sapply(p0s, function(p0) qbinom(alph, n_study, 
                                                       p0, lower.tail = FALSE))
    pows_avg_binom_th <- sapply(
      p1s, 
      function(p1) 1-pbinom(crit_vals_binom, 
                              n_study, p1, lower.tail = TRUE))

    pows_avg_binom_emp <- t(sapply(
      crit_vals_binom,
      function(crit_val) apply(test_func(x_binom, crit_val), 2, mean)
    ))

    pows_avg_binom <- rbind(melt(pows_avg_binom_emp), melt(pows_avg_binom_th))$value

    # calculate theoretical and empirical power for the test using the z-score
    crit_val_norm <- qnorm(alph, mean = 0, sd = 1, lower.tail = FALSE)
    
    pows_avg_Zn_th <- sapply(
      p1s_th, 
      function(p1) sapply(
        p0s, 
        function(p0) 1-pnorm(crit_val_norm, z_stat_binom(p0, p1, n_study),sd=1)))

    pows_avg_Zn_emp <- avg_evidence(1 * (Zn_emp > crit_val_norm), p0s, p1s_th)
    pows_avg_Zn <- rbind(melt(pows_avg_Zn_emp), melt(pows_avg_Zn_th))$value
    
    # calculate theoretical and empirical power for test using the vst
    pows_avg_Vn_th <- sapply(
      p1s_th, 
      function(p1) sapply(
        p0s, 
        function(p0) 1-pnorm(crit_val_norm, vst_binom(p0, p1, n_study),sd=1))
      )
    
    pows_avg_Vn_emp <- avg_evidence(1 * (Vn_emp > crit_val_norm), p0s, p1s_th)
    pows_avg_Vn <- rbind(melt(pows_avg_Vn_emp), melt(pows_avg_Vn_th))$value
    # maybe replace sd_th with sd_emp?
    # question: in this case, isn't the maximum evidence 1?
    # question: what standard deviation do we use? the empirical or the theoretical?
    # question: how do we know that approximate standard normality holds for vst?

    Zn <- cbind(Zn, pows_avg_Zn)
    pow_binom <- cbind.data.frame(Zn$p0, Zn$p1, NA, NA, Zn$th_emp, "binom", n_study, alph, NA, pows_avg_binom)
    colnames(pow_binom) <- cols_tot
    colnames(Zn) <- cols_tot
    Vn <- cbind(Vn, pows_avg_Vn)
    colnames(Vn) <- cols_tot
    evidence_binom <- rbind(evidence_binom, Zn, Vn, pow_binom)
  }

  diff.time <- Sys.time() - start.time
  print(paste0("round ", i, " is done in ", round(diff.time, 2), " ", attributes(diff.time)$units))
  i <- i + 1
}

save(evidence_binom, file = paste0("data/", name, ".RData")) # save data so that it doesn't have to be computed again and again


# question/todo: a lot of 0 values are created empirical SE and low n_study > leads to infinity values for the z-scores based on the empirical variance.
# can be amended by the continuity correction.

### plot values (to be transferred to a different script!) ----
# plot
# for math notation, see here: https://www.r-bloggers.com/math-notation-for-r-plot-titles-expression-and-bquote/

# NamingConvention: Ev_Binom_corr_1e4.pdf; EV_DIST_CORR_NUMSIM.pdf

pdf(paste0("./figs/", name, ".pdf"), onefile = TRUE)
i <- 1
for (sim in sims) {
  Zn <- sim[[1]]
  Vn <- sim[[2]]
  T_std_norm <- sim[[3]]
  T_std_norm$evd_mean <- T_std_norm$evd_mean

  Ts <- merger(dats = list(Zn, Vn, T_std_norm), ids = c("Zn", "Vn", "T_std_norm"))
  Ts <- Ts[Ts$th_emp == "th", ]

  T_plot <- ggplot(data = Ts, aes(x = p1, y = evd_mean / sqrt(n_studies[i]), col = id, group = interaction(p0, th_emp, id), linetype = th_emp)) +
    geom_line() + scale_color_manual(values = c("Zn" = 2, "Vn" = 3, "T_std_norm" = 4)) +
    geom_point(data = Ts[Ts$th_emp == "emp", ], aes(shape = factor(p0)), alpha = 0.5) +
    ylim(-3, 3) + xlim(0, 1) + ggtitle(bquote("Study size" == .(n_studies[i]))) + labs(x = expression(p[1]), y = "Evidence")
  print(T_plot)

  # T_plot <- ggplot(data=Ts,aes(x=p1,y=evd_mean/sqrt(n_studies[i]),col=factor(p0),group=interaction(p0,th_emp,id),linetype=id,alpha=th_emp)) +
  #   geom_line() + scale_linetype_manual(values=c("Zn"=2,"Vn"=1)) + scale_alpha_manual(values=c("emp"=1,"th"=0.5))+
  #   ylim(-3,3) + xlim(0,1) + ggtitle(bquote("Study size"==.(n_studies[i]))) + labs(x= expression(p[1]), y = "vst (solid), clt (dashed)")
  # print(T_plot)

  # plot the power
  pows <- sim[[4]]
  j <- 1
  for (pow in pows) {
    pow_bin <- pow[[1]]
    pow_Zn <- pow[[2]]
    pow_Vn <- pow[[3]]

    pows_merged <- merger(dats = list(pow_bin, pow_Zn, pow_Vn), ids = c("pow_bin", "pow_Zn", "pow_Vn"))
    # pows_merged <- pows_merged[pows_merged$th_emp=="th",]

    # pows_merged <- pows_merged[pows_merged$th_emp=="emp" & pows_merged$id =="pow_bin" & pows_merged$p0 == "0.7",]

    cutoff <- data.frame(x = c(0, 1), y = alphas[j], cutoff = factor(alphas[j]))
    p0_cutoff <- data.frame(x = rep(p0s, each = 2), y = rep(c(0, 1), length(p0s)))

    pow_plot <- ggplot(data = pows_merged, aes(x = p1, y = evd_mean, col = id, group = interaction(p0, th_emp, id), linetype = th_emp)) +
      geom_line() + scale_color_manual(values = c("pow_Zn" = 2, "pow_Vn" = 3, "pow_bin" = 4)) +
      geom_point(data = pows_merged[pows_merged$th_emp == "emp", ], aes(shape = factor(p0)), alpha = .5) +
      geom_line(aes(x, y, group = factor(x)), p0_cutoff, inherit.aes = FALSE, alpha = .2) + geom_line(aes(x, y), cutoff, inherit.aes = FALSE, alpha = .2) +
      ggtitle(bquote("Study size" == .(n_studies[i]) ~ ", alpha" == .(alphas[j]))) + labs(x = expression(mu[1]), y = "power")
    print(pow_plot)

    # pow_plot <- ggplot(data=pows_merged,aes(x=p1,y=evd_mean,col=factor(p0),group=interaction(p0,th_emp,id),linetype=id,alpha=th_emp)) +
    #   geom_line() + scale_linetype_manual(values=c("pow_bin"=2,"pow_Vn"=1)) +  scale_alpha_manual(values=c("emp"=1,"th"=0.5))+
    #   geom_line(aes(x, y,group=factor(x)),p0_cutoff,inherit.aes=FALSE,alpha=.2) + geom_line(aes(x, y),cutoff,inherit.aes=FALSE,alpha=.2) +
    #   ggtitle(bquote("Study size"==.(n_studies[i])~", alpha"==.(alphas[j]))) + labs(x= expression(p[1]), y = "power")
    # print(pow_plot)
    j <- j + 1
  }

  i <- i + 1
}
dev.off()
# questions:
# 1) where do the gaps between 0.1 and 0 as well as between 0.9 and 1 come from?
# 2) why do vst and clt diverge so much at these values?
# 3) why does exact binomial test have lower power for p0 and p1 large, but higher power for p0 low and p1 large?

# Calculate p-values and Type II errors for the different statistics
# https://www.cyclismo.org/tutorial/R/power.html

### discontinued functions ----
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
#
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

# calcualtes SE of binomial variable X/n
# SE_calculator <- function(p,n){
#   SE <- sqrt(p*(1-p)/n)
#   return(SE)
# }

# calc_Zn <- function(p0s,p1s,n_study,func){
#   evd_mean <- sapply(1:dim(p1s)[1], function(p1) sapply(p0s, function(p0) func(p0,p1s[p1,],n_study,SEs[p1,])))
#   evd_mean[is.infinite(evd_mean)] <- NA
#   return(evd_mean)
# }

# vst_binom <- function(p0,p1,n){
#   evd_mean <- 2*sqrt(n)*(asin(sqrt(p1))-asin(sqrt(p0)))
#   return(evd_mean)
# }
#
# #z-statistic based on CLT
# z_stat_binom <- function(p0,p1,n){
#   evd_mean <- (p1-p0)/sqrt((p1*(1-p1)/n))
#   return(evd_mean)
# }
# question: what changes when we use the theoretical variance instead of the empirical? fit should get better, shouldn't it?
# also: where do peaks come from?

# calc_evidenc <- function(p0s,p1s,n_study,func){
#   evd_mean <- apply(p1s, 1, function(p1) sapply(p0s, function(p0) func(p0,p1,n_study)))
#   evd_mean[is.infinite(evd_mean)] <- NA
#   return(evd_mean)
# }

# brings data into the correct form for plottings
# dat_transform <- function(T_avg,sd_evidence,th_emp,id,n_study,cols){
#   dat <- merge(melt(T_avg),melt(sd_evidence),by=c("Var1","Var2"),sort=FALSE)
#   dat <- cbind(dat,th_emp,id,n_study)
#   colnames(dat) <- cols
#   dat$p1 <- rep(p1s,times=1,each=length(p0s))
#   dat$p0 <- rep(p0s,times=length(p1s))
#   return(dat)
# }

# test_func <- function(val,crit_val){
#   return((val>crit_val)*1)
# }
#
# ci_coverage <- function(Ts,T_avg,crit_val){
#   T_avg <- as.vector(t(T_avg))
#   vals <- t(abs(sapply(1:dim(Ts)[1], function(i) Ts[i,]-T_avg[i])))
#   Ts_coverage <- 1-test_func(vals,crit_val)
#   return(Ts_coverage)
# }

# avg_evidence <- function(Ts){
#   Ts_avg <- matrix(apply(Ts,1,mean,na.rm=TRUE),length(p0s),length(p1s),byrow=TRUE)
#   return(Ts_avg)
# }
#
# sd_evidence <- function(Ts){
#   Ts_avg <- matrix(apply(Ts,1,sd,na.rm=TRUE),length(p0s),length(p1s),byrow=TRUE)
#   return(Ts_avg)
# }
# dat_transform_power <- function(pow_avg,th_emp,id,alph,n_study,cols){
#   dat <- cbind(melt(pow_avg),th_emp,id,alph,n_study)
#   colnames(dat) <- cols
#   dat$p1 <- rep(p1s,times=1,each=length(p0s))
#   dat$p0 <- rep(p0s,times=length(p1s))
#   return(dat)
# }