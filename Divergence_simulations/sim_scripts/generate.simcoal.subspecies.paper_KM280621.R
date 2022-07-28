generate.simcoal.subspecies.paper <- function(params, dA.lower.bound, dA.upper.bound, n.sims){
  
  param.p <- mito.g <- df <- list()
  successful.sims <- 0
  i <- 1
  
  while(successful.sims < n.sims) {
    
    p <- params[i,]
    demes = fscSettingsDemes(
      pop0 = fscDeme(deme.size = p$deme0, sample.size = p$sample.size),  
      pop1 = fscDeme(deme.size = p$deme1, sample.size = p$sample.size),  
      ploidy = 1
    )
    mtGen = fscBlock_dna(
      sequence.length = p$length,
      mut.rate = p$mut,
      recomb.rate = 0,
      transition.rate = 1/3,
      chromosome = 1
    )
    genetics = fscSettingsGenetics(
      mtGen, 
      num.chrom = 1
    )
    migration <- matrix(c(0, p$m, p$m, 0), nrow = 2)
    events = fscSettingsEvents(
      fscEvent(
        event.time = p$split.time,
        source = 0,
        sink = 1,
        prop.migrants = 1,
        new.size = p$pop.change,
        new.growth = 0,
        migr.mat = 0
      )
    )
    
    ex1.params <- fscWrite(demes = demes, 
                           genetics = genetics, 
                           migration = fscSettingsMigration(migration),
                           events = events,
                           label = paste("ex1.",i,sep=""))
    
    param.p[[i]] <- fscRun(ex1.params)
    
    arp.file <- fscReadArp(param.p[[i]])
    strata<-arp.file[,c(1,2,1)]  
    colnames(strata)[3] <- "mtGenome"  
    seqs <- strsplit(arp.file$C1B1_DNA, "")
    names(seqs) <- strata$id   
    rownames(strata) <- strata$id
    g <- df2gtypes(strata, 
                   ploidy = 1,
                   strata.col = NULL,
                   loc.col = 3,   
                   schemes = strata,
                   sequences = seqs,
                   description = title)
    mito.g[[i]] <- stratify(g,"deme")
    
    print("mito.g")
    
    df[[i]] <- calc.dA.and.PD(mito.g[[i]], dA.lower.bound, dA.upper.bound)
    if(!is.na(df[[i]][1])) successful.sims <- successful.sims + 1
    i <- i+1
    print(paste("i = ", i, "; successful = ", successful.sims, sep=""))
    
  } #end while
  
  df <- do.call('rbind', df)
  df$dA[which(df$dA <= 0)] <- range(na.omit(df$dA[which(df$dA >0)]))[1] 
  df <- cbind(df, log10(df[,c("dA","lower_0.95","upper_0.95")]))
  names(df)[5:9] <- c("dA.lci","dA.uci","log10.dA","log10.dA.lci","log10.dA.uci")
  
  # Label each simulation as to whether it would qualify as a population (1), subspecies (2), or species (3) based
  # on the dA and PD criteria
  
  df$dA.label <- 1
  df$dA.label[which(df$dA > dA.lower.bound)] <- 2
  df$dA.label[which(df$dA > dA.upper.bound)] <- 3
  df$dA.label[which(is.na(df$dA))] <- 0
  df$PD.label <- 1
  df$PD.label[which(df$diagnosability > 80)] <- 2
  df$PD.label[which(df$diagnosability > 95)] <- 3
  df$combo.label <- apply(df[,c("dA.label","PD.label")], 1, min)
  
  return(list(mito.g=mito.g,param.p=param.p, df=df))
}

