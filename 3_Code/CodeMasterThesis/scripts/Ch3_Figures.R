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
out_path <- "figs/chapter2/"

### properties that should be the same across all figures ---------------
font_size <- 16
font_family <- "Palatino Linotype"
A4 <- c(8.27, 11.69) # width and height of A4 page in inches
A5 <- c(5.8, 8.3) # width and height of an A5 page
corr_student <- "Corr" # either "MLE" for no correction or "Corr" for finite sample correction 
corr_bin <- "Ans" # either "MLE" for no correction or "Ans" for Anscombe correction
n_sim <- "1e+05"