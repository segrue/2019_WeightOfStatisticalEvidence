### import packages ------------------------------------------------------------
# for plotting
require("ggplot2") # http://www.sthda.com/english/wiki/ggplot2-essentials
require("ggpubr")
require("extrafont") # https://cran.r-project.org/web/packages/extrafont/README.html
# I used the Fira Sans font for the graphics - you probably need to install it if you're using Linux
# font_import(c("~/.local/share/fonts","/home/drosoneuro/.fonts")) # only run the first time / when new fonts were installed
# loadfonts() # only run the first time / when new fonts were installed

# for data handling
require("data.table")
# require("reshape2") #https://stackoverflow.com/questions/21563864/ggplot2-overlay-density-plots-r

# custom functions
source("./functions/helper_functions.R")

### define input paths (for data) and outpout paths (for figures)---------------
in_path <- "data/"
out_path <- "figs/presentation/"

### properties that should be the same across all figures ---------------
font_size <- 16
font_family <- "Fira Sans"
A4 <- c(8.27, 11.69) # width and height of A4 page in inches
A5 <- c(5.8, 8.3) # width and height of an A5 page
corr_student <- "MLE" # either "MLE" for no correction or "Corr" for finite sample correction
n_sim <- 1e5

n_studies <- c(5, 10, 20, 30, 40, 50, 100, 200, 500, 1000) # number of subjects per study

mu0s <- 0
mu1s <- seq(-2, 2, by = 0.1)
sgm0 <- sgm1 <- 2
alphas <- c(0.05)


### Figure 1A: Plot CI of Zn, Vn & Student T ----------------------------------
### (fig:CI_student)

# define figure name
figname <- paste0("presentation_fig1A_CI_student")

# load data
load(paste0("data/Ev_Student_MLE_", n_sim, "_20190504.RData"), verbose = TRUE)
ev_stud_mle <- evidence_student

load(paste0("data/Ev_Student_Corr_", n_sim, "_20190504.RData"), verbose = TRUE)
ev_stud_corr <- evidence_student

# plot üpwer curves
CI_plotter <- function(dat) {
  dat <- dat[th_emp == "emp" & (id %in% c("Vn", "Zn")), ]
  CI_plot <- ggplot(data = dat, aes(
    x = mu1, y = CI, col = id, group = interaction(mu0, id, th_emp, n_study)
  )) +
    geom_line(alpha = 1) +
    scale_color_manual(
      name = "CI:",
      labels = c(
        "Zn" = bquote(T[n] ~ "\U00B1" * z[(1 - alpha / 2)]),
        "Vn" = bquote(V[n] ~ "\U00B1" * z[(1 - alpha / 2)]),
        "Tn" = T[n] ~ "\U00B1" * t[(n - 1 * "," * 1 - alpha / 2)]
      ),
      values = c("Zn" = "blue", "Vn" = "red", "Tn" = "black")
    ) +
    scale_shape_manual(
      name = bquote(H[0] * ":"),
      labels = c("-2" = bquote(mu == -2), "0" = bquote(mu == 0), "2" = bquote(mu == 2)),
      values = c("-2" = 3, "0" = 1, "2" = 4)
    ) +
    # scale_linetype_manual(
    #   name = bquote(H[0] * ":"),
    #   labels = c(bquote(mu == -2), bquote(p == 0.5), bquote(p == 0.9)),
    #   values = c("0.1" = "solid", "0.5" = "dashed", "0.9" = "dotdash")
    # ) +
    geom_point(
      data = dat[mu1 %in% seq(-2, 2, 0.4), ],
      aes(shape = factor(mu0)), alpha = 0.5, size = 1
    ) +
    labs(x = bquote(mu[1]), y = "Coverage") +
    theme(
      text = element_text(size = font_size, family = font_family),
      plot.title = element_text(size = font_size)
    ) +
    geom_segment(aes(x = -2, y = .95, xend = 2, yend = 0.95),
                 col = "black", size = 0.3, linetype = "dotted"
    ) + coord_cartesian(ylim = c(0.75, 1))
  return(CI_plot)
}

p1 <- CI_plotter(ev_stud_mle[n_study == 5 & mu0 == 0, ]) +
  ggtitle(bquote("n" == 5))
p2 <- CI_plotter(ev_stud_corr[n_study == 5 & mu0 == 0, ])
# p3 <- CI_plotter(ev_stud_mle[n_study == 10 & mu0 == 0, ]) +
#   ggtitle(bquote("n" == 10))
# p4 <- CI_plotter(ev_stud_corr[n_study == 10 & mu0 == 0, ])
# p5 <- CI_plotter(ev_stud_mle[n_study == 30 & mu0 == 0, ]) +
#   ggtitle(bquote("n" == 30))
# p6 <- CI_plotter(ev_stud_corr[n_study == 30 & mu0 == 0, ])
# p7 <- CI_plotter(ev_stud_mle[n_study == 50 & mu0 == 0, ]) +
#   ggtitle(bquote("n" == 50))
# p8 <- CI_plotter(ev_stud_corr[n_study == 100 & mu0 == 0, ])

fig1A <- ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  align = "hv", legend = "top", common.legend = T
)

ggsave(
  filename = paste0(out_path, figname, ".pdf"), plot = fig1A,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)



### Figure 1B: Normal fit of Tn & Vn based on student distribution ------------
### (fig:normal_fit_student_[Corr])

# define figure name
figname <- paste0("presentation_fig1B_normal_fit_student")

# variables are
# cdf: cumulative probability
# p0: p under H0
# p1: p under H1
# quantile: x-value corresponding to cumulative probability prob
# id: th = cdf of the normal distribution with mu = ; Zn = cdf of z-score based on binomial
#     Vn = cdf of variance stabilised evidence measure Vn
# n_study: number of samples per simulated study

# plotting multiple plots in one window: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
# n=5

cdf_plotter <- function(dat) {
  # q_min <- min(dat[id == "th", ]$quantile)
  # q_max <- max(dat[id == "th", ]$quantile)
  q_min <- -6
  q_max <- 6
  cdf_plot <- ggplot(data = dat, aes(
    x = quantile, y = cdf, col = id,
    group = interaction(mu1, id), linetype = factor(mu1)
  )) +
    geom_line() + scale_color_manual(
      name = "F(x):",
      labels = c(
        bquote(Phi(x - E * "[" * T[n] * "]")), bquote(T[n]),
        bquote(V[n] - E * "[" * V[n] * "]" + E * "[" * T[n] * "]")
      ),
      values = c("th" = "black", "Tn" = "blue", "Vn" = "red")
    ) +
    scale_linetype_manual(
      name = bquote(H[1] ~ ":"),
      labels = c(
        "-2" =
          bquote(mu == -2), "0" = bquote(mu == 0),
        "2" = bquote(mu == 2)
      ),
      values = c("-2" = "solid", "0" = "dashed", "2" = "dotdash")
    ) +
    ylim(0, 1) + ggtitle(bquote("n" == .(dat$n_study[1]))) +
    labs(x = "x", y = "F(x)") + coord_cartesian(xlim = c(q_min, q_max)) +
    theme(
      text = element_text(size = font_size, family = font_family),
      plot.title = element_text(size = font_size)
    )
}

qq_plotter <- function(dat) {
  q_min <- min(dat[id == "th", ]$quantile)
  q_max <- max(dat[id == "th", ]$quantile)
  th_quantiles <- qnorm(seq(0.01, 0.99,
                            along.with = dat[id == "th" & mu1 == -2, quantile]
  ), 0, 1)
  dat$q_theoretical <- rep(th_quantiles, dim(unique(dat[, .(id, mu1)]))[1])
  q_plot <- ggplot(data = dat, aes(
    x = q_theoretical, y = quantile, col = id,
    group = interaction(mu1, id), linetype = factor(mu1)
  )) +
    geom_line() +
    scale_color_manual(
      name = "F(x):",
      labels = c(
        bquote(Phi(x - E * "[" * T[n] * "]")), bquote(T[n]),
        bquote(V[n] - E * "[" * V[n] * "]" + E * "[" * T[n] * "]")
      ),
      values = c("th" = "black", "Tn" = "blue", "Vn" = "red")
    ) +
    scale_linetype_manual(
      name = bquote(H[1] ~ ":"),
      labels = c(
        "-2" =
          bquote(mu == -2), "0" = bquote(mu == 0),
        "2" = bquote(mu == 2)
      ),
      values = c("-2" = "solid", "0" = "dashed", "2" = "dotdash")
    ) +
    labs(x = bquote(Phi^{
      -1
    } * (x)), y = bquote(F^{
      -1
    } * (x))) +
    coord_cartesian(ylim = c(q_min, q_max)) +
    theme(
      text = element_text(size = font_size, family = font_family),
      legend.position = "none"
    )
}

# plot all figures (two for each n_study) on the same page
# load data
load(paste0(in_path, "student_quantiles_", "MLE", "_5000_20190505.RData"),
     verbose = TRUE
)
p1 <- qq_plotter(student_quantiles[n_study == 5, ]) +
  ggtitle(bquote("n" == 5))
load(paste0(in_path, "student_quantiles_", "Corr", "_5000_20190505.RData"),
     verbose = TRUE
)
p2 <- qq_plotter(student_quantiles[n_study == 5, ])
# p3 <- cdf_plotter(student_quantiles[n_study == 10, ])
# p4 <- qq_plotter(student_quantiles[n_study == 10, ])
# p5 <- cdf_plotter(student_quantiles[n_study == 20, ])
# p6 <- qq_plotter(student_quantiles[n_study == 20, ])
# p7 <- cdf_plotter(quantiles[quantiles$n_study==30,])
# p8 <- qq_plotter(quantiles[quantiles$n_study==30,])
fig1B <- ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  align = "hv", legend = "top", common.legend = T
)
ggsave(
  filename = paste0(out_path, figname, ".pdf"), plot = fig1B,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)

### Figure 1C: Plot Tn & Vn against mu1 & compare with normal distribution ----
### (fig:evidence_student_[Corr])

# define figure name
figname <- paste0("presentation_fig1C_evidence_student")

# plot fit of empirical distribution of Zn and Vn with theoretical dist
fit_plotter <- function(dat) {
  dat <- dat[(id %in% c("Zn", "Vn")), ]
  y_min <- min(dat[th_emp == "th", ]$evd_mean / sqrt(dat$n_study[1]), na.rm = T)
  y_max <- max(dat[th_emp == "th", ]$evd_mean / sqrt(dat$n_study[1]), na.rm = T)
  fit_plot <- ggplot(data = dat, aes(
    x = mu1, y = evd_mean / sqrt(n_study),
    col = id, group = interaction(factor(mu0), th_emp, id),
    linetype = th_emp
  )) +
    geom_line() +
    geom_point(
      data = dat[mu1 %in% seq(-2, 2, 0.4) & th_emp == "emp", ],
      aes(shape = factor(mu0)), alpha = 0.5, size = 1
    ) +
    scale_color_manual(
      name = "Statistic:",
      labels = c("Vn" = bquote(V[n]), "Zn" = bquote(T[n])),
      values = c("Vn" = "red", "Zn" = "blue")
    ) +
    scale_linetype_manual(
      name = "emp/th:",
      values = c("th" = "dashed", "emp" = "solid")
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    scale_shape_manual(
      name = bquote(H[0] * ":"),
      labels = c("-2" = bquote(mu == -2), "0" = bquote(mu == 0), "2" = bquote(mu == 2)),
      values = c("-2" = 3, "0" = 1, "2" = 4)
    ) +
    labs(x = bquote(mu[1]), y = bquote(tau * "/" * sqrt(n))) +
    ggtitle(bquote("n" == .(dat$n_study[1]))) +
    theme(
      text = element_text(size = font_size, family = font_family),
      plot.title = element_text(size = font_size)
    )
  return(fit_plot)
}

# plot empirical variance around Zn and Vn
sd_plotter <- function(dat) {
  dat <- dat[(id %in% c("Zn", "Vn")) & mu0 == 0, ]
  y_min <- min(dat[th_emp == "th", ]$evd_mean / sqrt(dat$n_study[1]), na.rm = T)
  y_max <- max(dat[th_emp == "th", ]$evd_mean / sqrt(dat$n_study[1]), na.rm = T)
  # dat_th <- dat[th_emp == "th" & id=="Zn" & mu0 == 0, ]
  # dat_emp <- dat[(id %in% c("Tn", "Vn")) & mu0==0 & th_emp == "emp", ]
  sd_plot <- ggplot(
    data = dat[th_emp == "emp", ],
    aes(
      x = mu1, y = evd_mean / sqrt(n_study),
      col = id, group = interaction(factor(mu0), th_emp, id),
      linetype = th_emp
    )
  ) + geom_line() + geom_ribbon(aes(
    ymin = (evd_mean - evd_sd) / sqrt(n_study),
    ymax = (evd_mean + evd_sd) / sqrt(n_study), fill = id
  ),
  alpha = 0.5, linetype = "blank"
  ) +
    geom_ribbon(
      data = dat[th_emp == "th", ],
      aes(
        ymin = (evd_mean - evd_sd) / sqrt(n_study),
        ymax = (evd_mean + evd_sd) / sqrt(n_study)
      ),
      alpha = 1, fill = NA
    ) +
    geom_point(
      data = dat[mu1 %in% seq(-2, 2, 0.4) & th_emp == "emp", ],
      aes(shape = factor(mu0)), alpha = 0.5, size = 1
    ) +
    scale_color_manual(
      name = "Statistic:",
      labels = c("Zn" = bquote(T[n]), "Vn" = bquote(V[n])),
      values = c("Zn" = "blue", "Vn" = "red"),
      guide = F
    ) +
    scale_fill_manual(
      name = "Statistic:",
      labels = c("Zn" = bquote(T[n]), "Vn" = bquote(V[n])),
      values = c("Zn" = "blue", "Vn" = "red")
    ) +
    scale_shape_manual(
      name = bquote(H[0] * ":"),
      labels = c("-2" = bquote(mu == -2), "0" = bquote(mu == 0), "2" = bquote(mu == 2)),
      values = c("-2" = 3, "0" = 1, "2" = 4)
    ) +
    scale_linetype_manual(
      name = "emp/th:",
      values = c("th" = "dashed", "emp" = "solid")
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = bquote(mu[1]), y = bquote(tau * "/" * sqrt(n))) +
    theme(text = element_text(size = font_size, family = font_family))
  return(sd_plot)
}

# load data
load(paste0("data/Ev_Student_", "MLE", "_", n_sim, "_20190504.RData"), verbose = TRUE)
dat <- evidence_student
# p1 <- fit_plotter(dat[dat$n_study == 5, ]) +
#   ggtitle(bquote("n" == 5))
p1 <- sd_plotter(dat[dat$n_study == 5, ]) +
  ggtitle(bquote("n" == 5))
load(paste0("data/Ev_Student_", "Corr", "_", n_sim, "_20190504.RData"), verbose = TRUE)
dat <- evidence_student
p2 <- sd_plotter(dat[dat$n_study == 5, ])
# p3 <- fit_plotter(dat[dat$n_study == 10, ])
# p4 <- sd_plotter(dat[dat$n_study == 10, ])
# p5 <- fit_plotter(dat[dat$n_study == 20, ])
# p6 <- sd_plotter(dat[dat$n_study == 20, ])

fig1C <- ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  align = "hv", legend = "top", common.legend = T
)
ggsave(
  filename = paste0(out_path, figname, ".pdf"), plot = fig1C,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)


### Figure 2: Funnel plots explained -----------------------------------------
n_sim <- 1000
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

figname <- "presentation_funnel_plot"

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

fig2A <- ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  align = "hv", legend = "top", common.legend = T
)

fig2B <- ggarrange(p3, p4,
                   ncol = 2, nrow = 1,
                   align = "hv", legend = "top", common.legend = T
)

fig2C <- ggarrange(p5, p6,
                   ncol = 2, nrow = 1,
                   align = "hv", legend = "top", common.legend = T
)

ggsave(
  filename = paste0(out_path, figname, "_2A.pdf"), plot = fig2A,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)

ggsave(
  filename = paste0(out_path, figname, "_2B.pdf"), plot = fig2B,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)

ggsave(
  filename = paste0(out_path, figname, "_2C.pdf"), plot = fig2C,
  width = A4[1], height = 0.4 * A4[2], device = cairo_pdf
)

