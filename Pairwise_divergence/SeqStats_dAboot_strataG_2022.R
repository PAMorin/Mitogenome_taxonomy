#Creating and anlyzing a gtypes object for mtDNA data, using strataG v2.5.01

## install strataG from GitHub (v2.5.01)
# if (!require('devtools')) install.packages('devtools')
# devtools::install_github('ericarcher/strataG', build_vignettes = TRUE)
## install from CRAN
# rfPermute v2.1.7
# geneticRF v0.9 (from R3.4; DO NOT INSTALL from github, as Eric has updated it without changing the version number)

## data files required:
# fasta file of unique haplotypes (fasta_file)
# csv file of samples with associated haplotypes, and 1 column of strata designations
#   for each pairwise comparison (stratifications)

rm(list = ls())
options(stringsAsFactors = F, mc.cores = 2)
require(strataG)
require(swfscMisc)
require(geneticRF)
require(rfPermute)
require(phangorn)
require(adegenet)
require(parallel)

##################### PARAMETERS AND FILE NAMES ####
description <- "Bphy_mtGen"
stratifications <- c("NASH", "NANP", "NPSH")
stats <- c("PHIst", "Fst", "CHIsq")  # "Fst", "Fst_prime", "CHIsq", "D", "Fis", "Gst", "Gst_prime", "Gst_dbl_prime", "PHIst"
fasta_file <- "4_Bphy_mtGenHaps.fasta"
strata_file <- "Bphy.csv"
hap_header <- "Mitogenome_Hap" #column header for haplotype data in strata_file
PD_trees <- 500 # number of trees for random forest percent diagnosibility (suggest 5000 after testing)
boot_rep <- 100 # number of replicates for bootstrapping (used in dA_boot function; suggest 1000 after testing)
mdl <- c("TN93") #run model test first (below), then select best model for the data
gamma <- 1.5
pldy <- 1 #1 for mtDNA sequence, 2 for genotypes
num.cores <- 2
mydate <- Sys.Date()

# use jmodeltest (https://doi.org/10.1093/molbev/msn083) to determine, or select from list:
# can specify model as "raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81",
# "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock".
# If the model isn't listed as one of these, then we need to pick the closest model
# When comparing models fitted by maximum likelihood to the same data, the smaller
#  the AIC or BIC, the better the fit.
# TN93 = TrN in jmodeltest (Tamura-Nei 1993) is most similar to HKY and GTR models from jModelTest. 
# more info: http://evomics.org/resources/substitution-models/nucleotide-substitution-models/)
#####################

######## dA_boot function (caculates and bootstraps dA)
boot.dA <- function(g, nrep = boot_rep, num.cores = 1, ...) {
  .bootFunc <- function(rep.i, g) {
    ran.i <- unlist(tapply(
      getIndNames(g),
      getStrata(g),
      function(i) sample(i, length(i), replace = TRUE)
    ))
    boot.g <- g[ran.i, , ]
    nucleotideDivergence(boot.g, ...)$between$dA
  }

  dA.boot <- mclapply(1:nrep, .bootFunc, g = g, mc.cores = num.cores)
  dA.boot <- do.call(cbind, dA.boot)

  dA.95 <- t(apply(dA.boot, 1, swfscMisc::central.quantile))
  colnames(dA.95) <- paste("dA", colnames(dA.95), sep = ".")
  nd <- nucleotideDivergence(g, ...)
  nd$between <- cbind(nd$between[, 1:4], dA.95, nd$between[, 5:ncol(nd$between)])
  nd
}

#load all haplotype sequences as fasta file
Haps <- read.fasta(fasta_file) # need to decompress fasta.gz files first

#load strata file
strata.df <- readGenData(strata_file)
strata.df<-as.data.frame(strata.df)

###########################################################
# pairwise comparisons for all strata pairs
for(i in stratifications) {
stratum <- i

#Create gtypes object for hap sequences. Need to specify the sample id and haplotype columns from the
#strata file, so that haplotypes match the haplotype names in the fasta file.
row.names(strata.df) <- strata.df$id
g <- df2gtypes(strata.df[, c("id", hap_header)],
               ploidy = pldy,
               strata.col = NULL,
               loc.col = 2,
               schemes = strata.df,
               sequences = Haps,
               description = description)
#schemes = strata.df allows any column in the strata file to be used for re-stratification
#prior to analysis.

#set the stratum to be used in downstream analyses (see "stratum" setting above)
g <- stratify(g,stratum)

# Calculate Percent Diagnosability (Random Forest)
 rf <- gtypesRF(g, ntree = PD_trees)
 write.csv(confusionMatrix(rf$rp$rf, threshold=c(0.8, 0.9, 0.95)), file = paste(stratum,"_RF_PD.csv",sep = ""))

##############
#run modelTest {phangorn} to select the optimal model (based on AIC) for the data
gphyDat <- gtypes2phyDat(g, locus = 1)
modeltest <- modelTest(gphyDat)
write.csv(modeltest, paste(stratum, "_modeltest.csv",sep=""))


####################
#haplotypic diversity (overall for all samples in stratum, not divided by strata)
HapDiv <- heterozygosity(g, by.strata = FALSE, type = "expected") # exptdHet(g)
write.csv(HapDiv, paste(stratum,"_HapDiv.csv",sep = ""))

####################
#summarize data for unique haplotypes (No. samples, No. haplotypes, Hap diversity,
#% unique haps)
summaryfile<- summary(g)

####################
#fixed and variable sites:
varsites <- variableSites(g, bases = c("a", "c", "g", "t", "-"))
write.csv(varsites$site.freqs, paste(stratum,"_varsites.csv", sep=""))

fixedDif <- fixedDifferences(g, count.indels = TRUE, consec.indels.as.one = TRUE,
                             bases = c("a", "c", "g", "t", "-"))
write.csv(fixedDif$num.fixed, paste(stratum,"_fixedDif.csv", sep=""))

####################
#nucleotide divergence (should be calculated on the whole population set of sequences, 
  #not only unique haplotypes.
# nucleotideDivergence requires number of samples in each stratum to be >1 (or error "n < m").

nd <- nucleotideDivergence(g, model = mdl, gamma = gamma, variance = FALSE,
                           pairwise.deletion = TRUE) 
# write table of within-stratum nucleotide diversity (Nei's pi = nucleotide diversity)
write.csv(nd$within, paste(stratum,"_",mdl,"_NucDiv_w-in.csv", 
                                          sep=""))
# write table of dA (net divergence) between strata (Nei's dA between strata)
write.csv(nd$between, paste(stratum,"_",mdl,"_NucDiv_dA.csv", 
                                          sep=""))
#################
#population differentiation (divergence metrics based on "stats" selected above):
#Conduct overall and pairwise tests of population differentiation.
#set nrep to 1000 after testing, and select stats = "all" or selected stats from list

pop_overall <- overallTest(g, nrep = 1000, model = mdl)
write.csv(pop_overall$result, paste(stratum, "_", mdl, "_overall.summary.csv", sep=""))

pop_pairwise <- pairwiseTest(g, nrep = 1000, model = mdl) 
pw_sum <- pairwiseSummary(pop_pairwise)
pw_PHIst <- pairwiseMatrix(pop_pairwise, stat = "PHIst")

write.csv(pw_sum, paste(stratum, "_", mdl, "_pairwise.summary.csv", sep=""))

write.csv(pw_PHIst, paste(stratum, "_phist.pairwise.matrix.csv", sep = ""))
}

#################
# bootstrapped dA for each strata pair
row.names(strata.df) <- strata.df$id
g <- df2gtypes(strata.df[, c("id", hap_header)], #this creates a df with 2 columns for id and hap
               ploidy = pldy,
               strata.col = NULL,
               loc.col = 2,
               schemes = strata.df,
               sequences = Haps,
               description = description)
# bootstrap dA (in addition to calculating dA above for each strata. This does the bootstrap
# for each and combines all into a single output file "_dA_boot.csv"
# instead of overwriting 'g' with each call to 'stratify', add results to list:
da.ci1 <- sapply(stratifications, function(st) {
  cat(st, "\n") #prints to new lines for each stratification
  boot.dA(stratify(g, st), model = mdl, pairwise.deletion = TRUE)$between
})
da.ci1t <- t(da.ci1)
write.csv(da.ci1t, paste0(description,"_",mdl,"_dA_boot_",mydate,".csv"))

#################

# save the rdata
save.image(file=paste0(description, ".Rdata"))
