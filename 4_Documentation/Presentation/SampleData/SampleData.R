#set environment to environemnt where R-script is stored at
isRStudio <- Sys.getenv("RSTUDIO") == "1"
old_wd <- getwd()
if (isRStudio){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else{
  setwd(getSrcDirectory(function(x) {x}))
}

#number of tweets in each data set
ntweets_tot <- 9481484/4*1200 + 2594472 #each bigmemory matrix combines four binary files of tweet ratings; there are a total of 300 full quartets, with 
#1200&1201 representing an imcomplete quartet
ntweets_100 <- 17*2370371 + 2314697
ntweets_sick <- 2370371 + 1761279
ntweets <- c(ntweets_tot,ntweets_100,ntweets_sick)

#time window
time_window <- list(start=as.Date("2011-03-05"),end=as.Date("2015-07-11"))

#full set aggregated
load("df_agg_merged.RData")
df_agg_sample_full <- df_agg_merged[sample(.N, 6)]
df_agg_sample_raw <- df_agg_sample_full[,-c(7,8)]
df_agg_sample <- list(df_agg_sample_full,df_agg_sample_raw)

save(head_raw,ntweets,time_window,df_agg_sample,file="SampleData.RData")

#set environment back to env where script was called from
setwd(old_wd)
# getSrcDirectory(function(x) {x})
# sys.frame(1)3
# parent.frame(2)$ofile
# 
# this.dir <- dirname(parent.frame(2)$ofile) # frame(3) also works.
# setwd(this.dir)
# 
# getSrcDirectory(function(x) {x})
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
