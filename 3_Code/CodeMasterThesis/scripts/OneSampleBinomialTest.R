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
n_studies <- c(5, 10, 20, 30, 40, 50, 100)
n_sim <- 1e5
p0s <- c(0.1,0.5,0.9)
p1s <- seq(0.01, 0.99, by = 0.01)
alphas <- c(0.05)
corr <- F
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

# NamingConvention: Ev_Binom_corr_1e4.pdf; EV_DIST_CORR_NUMSIM.pdf
