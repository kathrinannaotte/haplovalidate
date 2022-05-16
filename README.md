# haplovalidate

## Installation

Before installing haplovalidate you need to make sure that all the dependencies are available. Please get the latest haploReconstruct from github https://github.com/popgenvienna/haploReconstruct . Haplovalidate will NOT WORK with the current CRAN haploreconstruct version 0.1.2 . 

     R (>= 3.6.0)
     psych (>= 1.8.12)
     stringr (>= 1.4.0)
     haploReconstruct (>= 0.1.3)
     data.table (>= 1.12.0)     


For now you need to install these manually. Once this is done you can proceed by downloading the latest release of haplovalidate. After the download you can install haplovalidate with the following R command:


        install.packages("/Path/To/haplovalidate_x.x.x.tar.gz", repos=NULL, type="source"
     
## Input Data
You need an object containing your allele frequencies in haploReconstruct format, which is created from a sync file.
   
    ## make sure to load version v0.1.3_3 from https://github.com/popgenvienna/haploReconstruct
    library(haploReconstruct)  
    
    ## you need to specify the order of your replicates/generations in the original sync file, below an example for an experiment with 5 replicates and 4 generations ordered by generations (i.e. F0_rep1, F0_rep2, ..., F59_rep4, F59_rep5)
    
    repl <- 1:5
    gens <- c(0,15,37,59)
    
    ### define which columns of the sync file contain the base population (e.g. the first 5 columns)
    base.pops <- c(rep(TRUE, length(repl)),rep(FALSE,length(repl)*(length(gens)-1)))
    
    ### define which columns should be used to polarize for the rising allele
    list(c(F0Rep1,F59Rep1),c(F0Rep2,F59Rep2),...))
    
    polaRise = list(c(1,26),c(2,27),c(3,28),c(4,29),c(5,30)) 
    
    ###define which columns should be included in the analysis
    compare <- rep(TRUE,length(repl)*length(gens))
    
    ### define which column contain what replicate 
    pop.ident <- rep(repl,length(gens))
    
    ### define which column contain which generation
    pop.generation  <- rep(gens, each=length(repl))

    ### load frequency data
    cands.all <- sync_to_frequencies(syncfile,base.pops=base.pops,header=FALSE,mincov=15,polaRise = polaRise)
   
You also need an object with the results of a CMH-Test (for a method incorporating genetic drift and poolSeq noise see https://github.com/popgenvienna/ACER).
The object should be a data.frame or data.table object containing chromosome and position (matching the candidates) and the corresponding CMH score (-log10(p-value)). Please remove missing values. 
 
    ## column names should be chr, pos and score 
    cmh <- readRDS("cmh.rds") 
    
You need to filter for SNPs with a significant allele frequency change

    cmh05 <- cmh[score> 1.3,.(chr,pos)] ## p-value < 0.05 
    cands <- merge(cands.all,cmh05,by=.(chr,pos))
    saveRDS(cands,"cands.rds")
    
## Usage

     ## install.packages("../haplovalidate_x.x.x.tar.gz",type="source",repos=NULL)

     library(haplovalidate)

     repl <- 1:5
     gens <- c(0,15,37,59)

     base.pops <- c(rep(TRUE, length(repl)),rep(FALSE,length(repl)*(length(gens)-1)))
     polaRise <- list(c(1,26),c(2,27),c(3,28),c(4,29),c(5,30)) 
     compare <- rep(TRUE,length(repl)*length(gens))
     pop.ident <- rep(repl,length(gens))
     pop.generation  <- rep(gens, each=length(repl))

    
     ## cands <- sync_to_frequencies(syncfile,base.pops=base.pops,header=FALSE,mincov=15,polaRise = polaRise)

     cands <- readRDS("candidates.rds")

     cmh <- readRDS("cmh.rds")
     ## columns chr(=chromosome), pos(=position) and (cmh) score needed

     ## get haplovalidate parameters
     parameters <- get.mncs.win(cands,cmh,wins=seq(0.1,10,0.05),mncs=0.01)
     print(parameters)

    ## run haplovalidate
    ## parameter takerandom specifies the number of SNPs that are sampled for cluster comparisons, decrease to reduce computational time at but with reduced accuracy 
    ## parameter filterrange defines the range for identifying dominant clusters in bp
    
    happy <- haplovalidate(cands=cands,cmh=cmh,parameters=parameters,pop.ident=pop.ident,pop.generation=pop.generation,base.pops=base.pops,compare=compare,takerandom=2000,filterrange=5000)
  
    ## generate summary plot
     plot.haplovalidate(blocks=happy$dominant_haplotypes,cmh,title="My beautiful haplotype blocks",label=TRUE)
     
## Citations

Otte, K. A., & C. Schlötterer, 2020. Detecting selected haplotype blocks in evolve and resequence experiments. Molecular Ecology Resources 93–109. https://doi.org/10.1111/1755-0998.13244

Susanne U. Franssen, Nicholas H. Barton, Christian Schlötterer, Reconstruction of Haplotype-Blocks Selected during Experimental Evolution, Molecular Biology and Evolution, Volume 34, Issue 1, January 2017, Pages 174–184, https://doi.org/10.1093/molbev/msw210

