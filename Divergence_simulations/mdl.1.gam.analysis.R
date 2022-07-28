# rm(list = ls())
library(mgcv)
scripts<-"sim_scripts"
source(paste0(scripts,"/gam.plotting.funcs.R"))
# load("Variable.Ne.m.mut.splitTime.300sim_TN93.Rdata")
load ("Variable.Ne.m.mut.splitTime.sim_16400bp_TN93_20yr.Rdata")
# load("Variable.Ne.m.mut.splitTime.300sim_JC69_dAPD.Rdata")

plottable.df <- cbind(df[-which(is.na(df$diagnosability)),],selected.sets[-which(is.na(df$diagnosability)),])
plottable.df$l.Ne <- log10(plottable.df$deme0)
plottable.df$l.mig.rate <- log10(plottable.df$m)
plottable.df$l.mut.rate <- log10(plottable.df$mut)
plottable.df$l.div.time <- log10(plottable.df$split.time)

label <- "300.sim"


#next two control tensor splint fit
bs <- "tp"
k <- 4

dA.zlim <- range(plottable.df$log10.dA)

gt0.1.form.diag <- cbind(right, wrong) ~  te(l.Ne, bs = bs, k = k) +  te(l.mig.rate, bs = bs, k = k) +  
  te(l.mut.rate, bs = bs, k = k) +  te(l.div.time, bs = bs, k = k)
eq0.1.form.diag <- cbind(right, wrong) ~  te(l.Ne, bs = bs, k = k) +  te(l.mut.rate, bs = bs, k = k) + 
  te(l.div.time, bs = bs, k = k)

mdl1.all.diag <- sim.gam(eq0.1.form.diag, plottable.df, paste(label, ".all", sep = ""))
mdl1.all.dA <- gam(form = log10.dA ~  te(l.Ne, bs = bs, k = k) +  te(l.mut.rate, bs = bs, k = k) +  
                     te(l.div.time, bs = bs, k = k), data = plottable.df, family = gaussian)
mdl1.mig.eq0.diag <- sim.gam(eq0.1.form.diag, subset(plottable.df, m == 0), paste(label, ".mig.eq0", sep = ""))
mdl1.mig.eq0.dA <- gam(form = log10.dA ~  te(l.Ne, bs = bs, k = k) +  te(l.mut.rate, bs = bs, k = k) +  
                     te(l.div.time, bs = bs, k = k), data = subset(plottable.df, m==0), family = gaussian)
mdl1.mig.gt0.diag <- sim.gam(gt0.1.form.diag, subset(plottable.df, m > 0), paste(label, ".mig.gt0", sep = ""))
mdl1.mig.gt0.dA <- gam(form = log10.dA ~  te(l.Ne, bs = bs, k = k) +  te(l.mig.rate, bs = bs, k = k) +  
                         te(l.mut.rate, bs = bs, k = k) +  te(l.div.time, bs = bs, k = k), 
                       data = subset(plottable.df, m > 0), family = gaussian)

fname <- paste(label, ".gam.results.rdata", sep = "")
save(df, dA.zlim, mdl1.all.diag, mdl1.all.dA, mdl1.mig.eq0.diag, mdl1.mig.eq0.dA, mdl1.mig.gt0.diag, mdl1.mig.gt0.dA, file = fname)

# run script "Figure 3 - pairwise_PM.R"