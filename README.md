# haplovalidate

## Installation

Before installing haplovalidate you need to make sure that all the dependencies are available. Please make sure to get the latest haploReconstruct from github https://github.com/popgenvienna/haploReconstruct

     R (>= 3.6.0)
     psych (>= 1.8.12)
     stringr (>= 1.4.0)
     haploReconstruct (>= 0.1.3)
     data.table (>= 1.12.0)     


For now you need to install these manually. Once this is done you can proceed by downloading the latest release of haplovalidate. After the download you can install haplovalidate with the following R command:

     install.packages("/Path/To/haplovalidate_0.1.0.tar.gz", repos=NULL, type="source")

## Usage

     ## install.packages("../haplovalidate_0.1.1.tar.gz",type="source",repos=NULL)

     library(haplovalidate)

     repl <- 1:5
     gens <- c(0,15,37,59)

     ## base.pops <- c(rep(TRUE, length(repl)),rep(FALSE,length(repl)*(length(gens)-1)))
     ## polaRise = list(c(1,26),c(2,27),c(3,28),c(4,29),c(5,30)) 
     ## cands <- sync_to_frequencies(syncfile,base.pops=base.pops,header=FALSE,mincov=15,polaRise = polaRise)

     cands <- readRDS("candidates.rds")

     cmh <- readRDS("cmh.rds")
     ## columns chr(=chromosome), pos(=position) and (cmh) score needed

     parameters <- get.mncs.win(cands,cmh,wins=seq(0.1,10,0.05),wincut=0.01)
     print(parameters)

     happy <- haplovalidate(cands,cmh,parameters,repl,gens,takerandom=2000,filterrang=5000)

     plot.haplovalidate(blocks=happy$dominant_haplotypes,cmh,title="port",label=F)
     
## Citations

Kathrin A. Otte, Christian Schlötterer A generalised approach to detect selected haplotype blocks in Evolve and Resequence experiments bioRxiv 691659; doi: https://doi.org/10.1101/691659

Susanne U. Franssen, Nicholas H. Barton, Christian Schlötterer, Reconstruction of Haplotype-Blocks Selected during Experimental Evolution, Molecular Biology and Evolution, Volume 34, Issue 1, January 2017, Pages 174–184, https://doi.org/10.1093/molbev/msw210

