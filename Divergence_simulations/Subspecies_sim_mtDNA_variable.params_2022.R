rm(list = ls())

# install from Github:
# strataG
# devtools::install_github("ericarcher/stratag", ref = "redevel2018")
##  https://github.com/EricArcher/strataG/tree/redevel2018

# geneticRF
# devtools::install_github("EricArcher/geneticRF")

# requires fastsimcoal2 (fsc26: http://cmpg.unibe.ch/software/fastsimcoal26/):
# fastsimcoal2 is not included with 'strataG' and must be downloaded separately. 
# Additionally, it must be installed such that it can be run from the command 
# line in the current working directory. The function fscTutorial() will open 
# a detailed tutorial on the interface in your web browser.

library(strataG)
library(tidyverse)
library(rfPermute) # from cran
library(geneticRF) # from github (see above).
library(adegenet)
library(swfscMisc)
library(latticeExtra)
library(ape)
# setwd("/Users/philmorin/Documents/Mol_Ecol_Lab/_Projects/\ Subspecies_mitogenomes/Simulations/Subsp_1ksims_16.4kbp")
wd<-getwd()
scripts <- paste0(wd,"/sim_scripts") # assumes all scripts in the same directory

dA_PD <- file.path(scripts,"calc.dA.and.PD_KM280621.R")
simco <- file.path(scripts, "generate.simcoal.subspecies.paper_KM280621.R")
source(dA_PD)
source(simco)

dA.lower.bound <- 0.0006  # lower boundary on dA for qualifying as a subspecies
dA.upper.bound <- 0.008   # upper boundary on dA for qualifying as a subspecies

NumSims<-1000 # 100000 after testing
sims.per.category <- 100 # 1000 after testing
####################################################
# Generate parameter sets

split.time <- floor(10^runif(NumSims, 1.398, 5.398)) #~25 to ~250,000 generations! (25*15=375; 250k*15=3.75M)
Ne <- floor(10^runif(NumSims, 2, 5)) # 100 to 100,000
m <- 10^runif(NumSims, -9, -2) # 10e-9 to 10e-2
mut <- 10^runif(NumSims, -8, -5.5) # 10e-8 to 10e-5.5 subst/site/gen (publication average = 1.0878E-2 subst/site/Myr (=1.08E-8/yr; 1.63E-7/gen))
# for CR (from Duchene et al 2011 (Table 3), CR 95% HPD=6.08E-3 Sub/site/MY = 6.08E-9 sub/site/yr = 1.216E-7 sub/site/gen, 20yr gen)
# for mitogenome(from Duchene et al, CR 95% HPD=4.24E-3 Sub/site/MY = 4.24E-9 sub/site/yr = 8.48E-8 sub/site/gen, 20yr gen)
pop.change <- 2
length <- 16400 #simulated sequence length
gen.time <- 20
pop.num <- 2
sample.size <- 20
num.cores <- 1
model <- "TN93"

all.params <- data.frame(cbind(split.time, deme0=Ne, deme1=Ne, pop.change, m, 
                               mut, Nem=Ne*m, theta=Ne*mut, length, gen.time, 
                               pop.num, sample.size, num.cores))

title = paste0("Variable.Ne.m.mut.splitTime.sim_",length,"bp_",model,"_",gen.time,"yr")

####################################################
# Run simulations and calculate diagnosability and dA (generate.simcoal.subspecies.paper now also calculates PD and dA)

# Nem = 0 simulations
Nem.0.indices <- sample(1:nrow(all.params), (sims.per.category*2))
Nem.0.sets <- all.params[Nem.0.indices,]
Nem.0.sets$m <- Nem.0.sets$Nem <- 0
Nem.0.sims <- generate.simcoal.subspecies.paper(Nem.0.sets,dA.lower.bound,dA.upper.bound,sims.per.category)
Nem.0.sets <- Nem.0.sets[1:length(Nem.0.sims[[1]]),]

all.params <- all.params[-Nem.0.indices,]

# Nem < 1 simulations
Nem.lt1.indices <- sample(which(all.params$Nem < 1), (sims.per.category*2))
Nem.lt1.sets <- all.params[Nem.lt1.indices,]
Nem.lt1.sims <- generate.simcoal.subspecies.paper(Nem.lt1.sets,dA.lower.bound,dA.upper.bound,sims.per.category)
Nem.lt1.sets <- Nem.lt1.sets[1:length(Nem.lt1.sims[[1]]),]

all.params <- all.params[-Nem.lt1.indices,]

# Nem >= 1 simulations
Nem.gt1.sets <- all.params[sample(which(all.params$Nem >= 1), (sims.per.category*2)),]
Nem.gt1.sims <- generate.simcoal.subspecies.paper(Nem.gt1.sets,dA.lower.bound,dA.upper.bound,sims.per.category)
Nem.gt1.sets <- Nem.gt1.sets[1:length(Nem.gt1.sims[[1]]),]

param.p <- c(Nem.0.sims$param.p, Nem.lt1.sims$param.p, Nem.gt1.sims$param.p)
mito.g <- c(Nem.0.sims$mito.g, Nem.lt1.sims$mito.g, Nem.gt1.sims$mito.g)
df <- rbind(Nem.0.sims$df, Nem.lt1.sims$df, Nem.gt1.sims$df)
selected.sets <- rbind(Nem.0.sets, Nem.lt1.sets, Nem.gt1.sets)

var.sites <- do.call('rbind', lapply(mito.g, function(g){
  x <- variableSites(g)
  v <- dim(x$site.freqs)[2]
  c(v,base.freq(getSequences(g)))
}))

var.sites.and.params <- data.frame(selected.sets, var.sites)
save(title, scripts, all.params, selected.sets, param.p, mito.g, df, var.sites.and.params, file=paste(title,"Rdata",sep="."))

# save.image(file = paste0(title,"_sim.Rdata"))
date()

# Run mdl.1.gam.analysis.R
