calc.dA.and.PD <- function(g, dA.lower.bound, dA.upper.bound){
  
  rf <- geneticRF::gtypesRF(g, gene = 1, pairwise = FALSE, conf.level = 0.95)
  
  nuc.div <- nucleotideDivergence(g, model = "TN93")$between
  
  if(is.na(nuc.div$dA)) dA.ci <- c(NA, NA) else dA.ci <- central.quantile(nuc.div[6:10])
  if(is.na(nuc.div$dA) || is.null(rf)) return(c(diagnosability=NA,diag.lci=NA,diag.uci=NA, dA=NA, lower_0.95=NA, upper_0.95=NA)) else {
    return(data.frame(c(rf$smry[3:5]), dA=nuc.div$dA, lower_0.95=dA.ci[1], upper_0.95=dA.ci[2]))}
}

