### import packages -----
# for data handling
require("reshape2") # https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r
require("data.table")

# for plotting
require("ggplot2")
require("gridExtra")

# custom functions
source("./functions/helper_functions.R")

### Define starting assumptions and conditions ----
# hypotheses:
# H0: p = mu0
# H1: p > mu0

seed <- 20190514
set.seed(seed)

n_sim <- 1e2 # number of times the setting below shall be simulated
n_studies <- c(5, 10, 20, 30, 40, 50, 100, 200, 500, 1000) # number of subjects per study

mu0s <- 0
mu1s <- seq(-2, 2, by = 0.1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)

T_corr <- F # whether correction should be used or not when calculating vst
correction <- ifelse(T_corr, "Corr", "MLE")

name <- paste0("DF_Ev_Mean_", correction, "_", n_sim, "_", seed) # name to load and save data etc.
fig_name <- paste0("Funnel_Ev_Mean_", correction, "_", n_sim, "_", seed)

### load required data sets ----
load(paste0("data/", name, ".RData"), verbose = TRUE) # load data set containing simulated values
mu1_select <- 0.2
n_study_min <- 10
n_study_max <- 200
x_lim <- 5
y_lim <- 8
rng <- range(c(n_study_min, n_study_max))
dat <- evidence_df[mu1 == mu1_select & n_study %between% rng, ]
ctgs <- unique(dat$n_study)[unique(dat$n_study) %between% rng]
sample_seed <- 10

### simulate publication bias
## Situtation 1A:
## no publication bias - all studies where n_study between 10 and 200 are included (n=100 per study size group)
label <- "1A"
n_tot <- length(ctgs) * n_sim

clt <- dat[id == "clt", ]
vst <- dat[id == "vst", ]
stud <- dat[id == "stud", ]

funnel_plotter(clt, vst, vst, xlim = x_lim, ylim = y_lim, figname = paste0(fig_name, "_", label, "_", n_tot, ".pdf"), ctgs)

## Situtation 1B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1 == 1, ]
vst_sig <- vst[H1 == 1, ]
stud_sig <- stud[H1 == 1, ]
label <- "1B"

funnel_plotter(clt_sig, vst_sig, stud_sig, xlim = x_lim, ylim = y_lim, figname = paste0(fig_name, "_", label, "_", n_tot, ".pdf"), ctgs)

## Situtation 2A:
## only keep 100 studies in total, but weigh according to study size;
## choose studies with lower sample size with higher probability
n_studies_selected <- n_studies[n_studies %between% rng]
probs <- c(0.29, 0.24, 0.2, 0.14, 0.09, 0.03, 0.01)
n_tot <- 100
label <- "2A"

clt <- select_studies(dat, probs, n_studies_selected, n_select = n_tot, sample_seed, T_id = "clt")
vst <- select_studies(dat, probs, n_studies_selected, n_select = n_tot, sample_seed, T_id = "vst")
stud <- select_studies(dat, probs, n_studies_selected, n_select = n_tot, sample_seed, T_id = "stud")

stopifnot(sum(clt$mu1_hat - vst$mu1_hat) == 0, sum(stud$mu1_hat - vst$mu1_hat) == 0)

funnel_plotter(clt, vst, vst, xlim = x_lim, ylim = y_lim, figname = paste0(fig_name, "_", label, "_", n_tot, ".pdf"), ctgs)

## Situtation 2B:
## only keep studies which turned out to be significant
clt_sig <- clt[H1 == 1, ]
vst_sig <- vst[H1 == 1, ]
stud_sig <- stud[H1 == 1, ]
label <- "2B"

funnel_plotter(clt_sig, vst_sig, stud_sig, xlim = x_lim, ylim = y_lim, figname = paste0(fig_name, "_", label, "_", n_tot, ".pdf"), ctgs)

## Situtation 2C:
## only keep studies which turned out to be significant and small percentage of non-significant studies
sel_prob <- 0.1 # percentage of significant studies to keep
probs_mix <- c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1)
clt_mix <- rbind(clt_sig, select_studies(clt[H1 == 0, ], sel_prob, n_studies_selected, n_select = NULL, sample_seed, T_id = "clt"))
vst_mix <- rbind(vst_sig, select_studies(vst[H1 == 0, ], sel_prob, n_studies_selected, n_select = NULL, sample_seed, T_id = "vst"))
stud_mix <- rbind(stud_sig, select_studies(stud[H1 == 0, ], sel_prob, n_studies_selected, n_select = NULL, sample_seed, T_id = "stud"))
label <- "2C"

funnel_plotter(clt_mix, vst_mix, stud_mix, xlim = x_lim, ylim = y_lim, figname = paste0(fig_name, "_", label, "_", n_tot, ".pdf"), ctgs)

# function to calculate the aggregated mean of the studies
aggregate_mean <- function(dat) {
  agg_mean <- sum(dat$mu1_hat * dat$n_study) / sum(dat$n_study)
  return(agg_mean)
}
aggregate_mean(stud)
aggregate_mean(clt)
aggregate_mean(vst)

### Calculate number of papers stored in file drawer based Rosenthal 1984 -----
# problem: Filedrawer problem only works for checking whether there is no null effect;
# it doesn't work in situations in which there really is an effect
dat_list <- list(clt, vst, stud, clt_sig, vst_sig, stud_sig)

filedrawer <- c()
i <- 0
for (dat in dat_list) {
  filedrawer <- c(filedrawer, calc_filedrawer(dat))
}

filedrawer

### Reweight means by publication probability ---------------------------------
check_reweighting(clt_mix)

pub_prob <- calculate_pub_prob(dat, p = 0.1)
weigh_mean <- reweight_mean(dat, pub_prob)

### Trim and fill based methods -----------------------------------------------
# Note: trim and fill only works if assumption that most extrem values are
# omitted holds
clt_mix_filled <- trim_and_fill(clt_mix)
clt_mix_filled_2 <- trim_and_fill(clt_mix, pub_prob = 0.1)

aggregate_mean(clt_mix)
aggregate_mean(clt_mix_filled)
aggregate_mean(clt_mix_filled_2)

### Calculate reweighted mean based MLE estimator -----------------------------

# Andrews-Kasy maximum likelihood estimator, assuming pub_prob is known

# expected publication probability
exp_pub_prob <- function(theta, alph, p) {
  quant <- qnorm(alph, 0, 1, lower.tail = FALSE)
  expect <- p * pnorm(quant, theta, 1) + 1 * pnorm(quant, theta, 1, lower.tail = FALSE)
  return(expect)
}
exp_pub_prob(0, 0.05, 0.2)

# likelihood function
likeli <- function(z, theta) {
  lik <- dnorm(z, mean = theta, 1) # same as: dnorm(z-theta,0,1)
  return(lik)
}

# calculate truncated likelihood of theta given the data
trunc_likeli <- function(dat, theta, p) {
  alph <- dat$alpha[1]
  z <- dat$Tn / sqrt(dat$n_study)
  # checking whether all alpha levels and all transformations are the same
  # is done by "calculate_pub_prob"
  lik <- calculate_pub_prob(dat, p) / exp_pub_prob(theta, alph, p) * likeli(z, theta)
  #lik <- calculate_pub_prob(dat, p) * likeli(z, theta)
  #idea: EM algorithm for this value
  return(lik)
}

thetas <- seq(-4, 4, by = 0.01)

likelihoods <- sapply(thetas, function(theta) prod(trunc_likeli(clt_mix, theta, p = 0.1)))
plot(thetas, likelihoods)
T_corr <- thetas[which(max(likelihoods) == likelihoods)]
T_corr

prod(calculate_pub_prob(clt_mix, 0.1) / exp_pub_prob(0.54, 0.05, 0.1) * likeli(clt_mix$Tn, 0.54))

theta <- 2

tes <- calculate_pub_prob(clt_mix,0.1)/exp_pub_prob(theta,0.05,0.1)*likeli(clt_mix$Tn,theta)
plot(clt_mix$Tn,tes)

# transform Z scores into estimates of mu
z_to_mu <- function(z, n, sgm_X) {
  mu_hat <- z * sgm_X / sqrt(n)
  return(mu_hat)
}

z_to_mu(T_corr,sum(clt_mix$n_study),clt_mix$sgm_th[1])



Tn_to_mu_hat <- function(T_corr, dat) {
  mu_hat <- T_corr * sqrt(sum(dat$n_study * dat$sgm_th^2)) / sum(dat$n_study)
  return(mu_hat)
}

Tn_to_mu_hat(T_corr, clt_mix)

# question: how to translate corrected T in tot correcte mu1_hats?
agg_sgm <- function(dat) {
  var_mu <- (dat$sgm_th / sqrt(dat$n_study))^2
  agg_var <- sum(var_mu * (dat$n_study)^2) / (sum(dat$n_study))^2
  return(sqrt(agg_var))
}

sgm_agg <- agg_sgm(clt_mix)

# implement same as above but in order to correct mu1_hats:
likeli <- function(x, mu) {
  lik <- dnorm(x, mean = mu, 2) # same as: dnorm(z-mu,0,)
  return(lik)
}

# calculate total T

calc_T_tot <- function(dat) {
  k <- length(dat$mu1_hat)
  mu_hat_tot <- sum(dat$mu1_hat * dat$n_study) / sum(dat$n_study)
  sgm_tot <- sum((dat$sgm_th * dat$n_study)^2) / sum(dat$n_study)^2
  T_tot <- mu_hat_tot / sgm_tot * sqrt(k)
  return(T_tot)
}

calc_T_tot(clt_mix)

aggregate_mean(clt_mix)

# calculate median bias
median_bias <- function(mu, alph, n) {
  z <- rnorm(n, mean = mu, 1)
  p <- calc_pub_prob(z, alph)
  idx <- which(1 == ifelse(p == 1, 1, rbinom(n - sum(p == 1), 1, 0.1)))
  med_bias <- median(z[idx]) - mu
  return(med_bias)
}

thetas <- seq(0, 5, by = 0.01)
plot(thetas, sapply(thetas, function(theta) median_bias(theta, alph = 0.05, n = 1e6)), type = "l")

### Calculate publication probability of studies according to Andrew & Kasy:
# Assumption1: distribution of results X* in latent studies given the true effects
# Theta*, f(x*|Theta*) is known

# systemic replication studies

# https://cran.r-project.org/web/packages/fitdistrplus/vignettes/paper2JSS.pdf
require("fitdistrplus")
plotdist(clt_mix$mu1_hat, histo = TRUE, demp = TRUE)
plotdist(clt_mix$Tn, histo = TRUE, demp = TRUE)
trunc_dens <- fitdist(clt_mix$mu1_hat, "norm")
trunc_dens_T <- fitdist(clt_mix$Tn, "norm")

# fit density: https://cran.r-project.org/doc/contrib/Ricci-distributions-en.pdf
trunc_dens_T <- density(clt_mix$Tn)
x <- trunc_dens_T$x
lik_trunc <- trunc_dens_T$y

lookup_id <- function(x, dens) {
  id <- which(abs(x - dens$x) == min(abs(x - dens$x)))
  return(id)
}

calc_emp_density <- function(dat) {
  z <- dat$Tn
  trunc_dens_emp <- density(z)
  ids <- sapply(z, function(z_i) lookup_id(z_i, trunc_dens_emp))
  z_dens_emp <- trunc_dens_emp$y[ids]
  return(z_dens_emp)
}

calc_th_density <- function(z, mu, probs, exp_p) {
  return(likeli(z, mu) * probs / exp_p)
}

# install.packages("NMOF")
require("NMOF")

trunc_density_likeli <- function(x_mu, dat, alph) {
  x <- x_mu[[1]]
  mu <- x_mu[[2]]
  z <- dat[alpha == alph, ]$Tn
  probs <- calc_pub_prob(z, alph, x)
  exp_p <- exp_pub_prob(mu, alph, x)
  trunc_dens_th <- calc_th_density(z, mu, probs, exp_p)
  trunc_dens_emp <- calc_emp_density(dat)
  err <- sum(abs(trunc_dens_th - trunc_dens_emp))
  return(err)
}

z <- clt_mix$Tn
probs <- calc_pub_prob(z, 0.05, 0.1)
exp_p <- exp_pub_prob(1.67, 0.05, 0.1)
th_density <- calc_th_density(clt_mix$Tn, 1.67, probs, exp_p)
emp_density <- calc_emp_density(clt_mix)

plot(z, calc_emp_density(clt_mix), ylim = c(0, 1))
points(z, th_density, col = "red")
plot(z, th_density)

trunc_density_likeli(list(0.1, 1), clt_mix, 0.05)

start.time <- Sys.time()
errs <- gridSearch(fun = trunc_density_likeli, levels = list(x, mus), dat = clt_mix, alph = 0.05, method = "multicore")
print(Sys.time() - start.time)

test_x <- seq(-2, 2, by = 0.01)
test <- likeli(test_x, 0.2) * ifelse(test_x > 1.64, 1, 0.1) / exp_pub_prob(0.2, 0.05, 0.1)
plot(test_x, test)

a <- optim(c(0.5, 0), trunc_density_likeli, dat = clt_mix, alph = 0.05)

x <- seq(0, 1, by = 0.01)
mus <- seq(-3, 3, by = 0.01)


plot(trunc_density_sim(clt_mix, 2, 0.1, 0.05))


plot(x, calc_pub_prob(x, trunc_dens_T, mean(clt_mix$Tn)))
plot(z, calc_pub_prob(z, trunc_dens_T, 0))

pub_prob <- calc_pub_prob(z, trunc_dens_T, 0.2)

rew_mean <- reweighted_mean(clt_mix, pub_prob)
rew_mean

aggregate_mean(clt_mix)



pub_probs_lik <- sapply(mus, function(mu) prod(calc_pub_prob(x, trunc_dens_T, mu)))
T_corr <- mus[which(max(pub_probs_lik) == pub_probs_lik)]
T_corr

plot(x, calc_pub_prob(x, trunc_dens_T, 2))

test <- sapply(mus, function(mu) prod(calc_pub_prob(z, lik_trunc, mu)))
plot(mus, test)
T_corr <- mus[which(max(test) == test)]
T_corr

equ <- function(z, mu) {
  res <- dnorm(z - mu, 0, 1) * (z - mu) * ((z - mu) - 1) + 1
  return(res)
}

results <- sapply(mus, function(mu) equ(z[1], mu))
plot(mus, results, type = "l")
id <- which(abs(0 - results) == min(abs(0 - results)))
equ(z[1], 2)

# draw from kernel density estimator: https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate?rq=1



## Gopas Selection
install.packages("metasens")
require("metasens")



###### Deprecated Functions ------

# Hansen-Hurwitz estimator & #Hansen-Hurwitz estimator
# https://newonlinecourses.science.psu.edu/stat506/node/15/ https://newonlinecourses.science.psu.edu/stat506/node/15/
sel_prob <- 0.1
quant_probs <- seq(0, 1, by = sel_prob)
clt_mix_quants <- quantile(clt_mix$mu1_hat, quant_probs)
hans_hurwitz_estimator <- function(n, sel_prob, quants) {
  # selection <- runif(n,min=0,max=1)
  samples <- runif(n, min(quants), max(head(quants, -2)))
  idx <- sapply(samples, function(s) min(which((s < tail(quants, -1)) == TRUE)))
  diffs <- (tail(quants, -1) - head(quants, -1))
  probs <- diffs / sum(diffs)
  prob_inclusion <- 1 - (1 - probs)^n
  prob_inclusion <- sapply(idx, function(i) prob_inclusion[i])
  # med <- head(quants,-1)+(tail(quants,-1)-head(quants,-1))/2
  hans_hur <- mean(samples / prob_inclusion)
}

horvitz_thompson_estimator <- function(dat) {
  y_i <- dat$mu1_hat * dat$n_study
  n <- length(y)
  probs <- 0.1 + 0.9 * dat$H1
  pi_i <- 1 - (1 - probs)^n
  print(pi_i)
  horv_thomp <- sum(y_i / pi_i) / sum(1 / pi_i * dat$n_study)
}
horv_thomp <- horvitz_thompson_estimator(clt_mix)
horv_thomp

hans_hurwitz_estimator <- function(n, sel_prob, quants) {
  # selection <- runif(n,min=0,max=1)
  samples <- runif(n, min(quants), max(head(quants, -2)))
  idx <- sapply(samples, function(s) min(which((s < tail(quants, -1)) == TRUE)))
  diffs <- (tail(quants, -1) - head(quants, -1))
  probs <- diffs / sum(diffs)
  prob_inclusion <- 1 - (1 - probs)^n
  prob_inclusion <- sapply(idx, function(i) prob_inclusion[i])
  # med <- head(quants,-1)+(tail(quants,-1)-head(quants,-1))/2
  hans_hur <- mean(samples / prob_inclusion)
}


hans_hurwitz_estimator <- function(dat) {
  y_i <- dat$mu1_hat * dat$n_study
  n <- length(y_i)
  probs <- 0.1 + 0.9 * dat$H1
  hans_hur <- sum(y_i / probs) / sum(1 / probs * dat$n_study)
}

hans_hur <- hans_hurwitz_estimator(clt_mix)
hans_hur


#####