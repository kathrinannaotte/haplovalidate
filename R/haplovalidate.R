#' @import haploReconstruct
#' @export



haplovalidate <- function(cands,cmh,parameters,repl,gens,takerandom=2000,filterrange=5000){

    base.pops <- c(rep(TRUE, length(repl)), rep(FALSE, length(repl) * (length(gens) -  1)))
    compare <- c(rep(rep(TRUE, length(gens)), length(repl)))
### haploReconstruct parameters
    min.minor.freq <- 0
    max.minor.freq <- 1
    minfreqchange <- 0
    minrepl <- 1
    min.lib.frac <- 0.75
    thres.ttest <- 0.025
    min.cl.size <- 20
    min.inter <- ceiling(min.cl.size/5)
    ####add this parameters to the function for more user freedom?
    findthreshold=10
    cor.low=0.3
    cor.high=0.9
    rawsteps=0.1
    finesteps=0.01
###functions
    ##asin sqrt transform allele frequency for normality 
    transform.af <- function(af) {
        af.sqrt <- asin(sqrt(af))
        af.transf <- t(af.sqrt)
        af.scale <- scale(af.transf, center = TRUE, scale = TRUE)
        return(af.scale)
    }


### calculate cluster correlations
    comclust <- function(x,cand.clust,takerandom,afcolis) {
        sub <- cand.clust[tag == x, ]
        if (nrow(sub) > takerandom) {
            red.pos <- sample(sub$pos, takerandom)
            pos.indi <- sub$pos %in% red.pos
            sub <- sub[pos.indi, ]
        }
        cluster.af <- sub[, afcolis, with = FALSE]
        cluster.scale <- transform.af(cluster.af)
        clustcor.sub <- median(cor(cluster.scale.i, cluster.scale), 
                               na.rm = TRUE)
        return(clustcor.sub)
    }


### cluster validation
    chromis <- unique(cands$chr)
    final <- c()
    hapval.result <- list()
    afcolis <- colnames(cands)[7:ncol(cands)]
    if (length(chromis) > 1) {
        ts.all <- list()
        for (d in chromis){
            winsize <- round(as.numeric(parameters[chr == d, win]) * 1000000,0)
            ts.all[[d]] <- try(initialize_SNP_time_series(chr = cands[chr == d, chr],
                                                          pos = cands[chr == d, pos],
                                                          base.freq = cands[chr == d, basePops],
                                                          lib.freqs = cands[chr == d, afcolis, with = F],
                                                          pop.ident = c(rep(repl, each = length(gens))),
                                                          pop.generation = rep(gens, length(repl)), 
                                                          use.libs = compare, min.minor.freq = min.minor.freq, 
                                                          max.minor.freq = max.minor.freq, winsize = winsize, 
                                                          win.scale = "bp", min.lib.frac = min.lib.frac, minfreqchange = minfreqchange, 
                                                          minrepl = minrepl),silent=T)
        }

        ts.filter <- ts.all[unlist(lapply(ts.all,isS4))]
        chromis <- names(ts.filter)

        cluster.snps <- c()
        for (chromo in chromis) {
            ts <- ts.filter[[chromo]]
            for (c in rev(seq(cor.low, cor.high, 0.1))) {
                print(paste("chr ", chromo, " cluster corr ", c, sep = ""))
                hbs <- reconstruct_hb(ts, chrom = chromo, min.cl.size = min.cl.size, 
                                      min.cl.cor = c, min.inter = min.inter, single.win = T, 
                                      transf = TRUE, arcsine = TRUE, scaleSNP = TRUE)
                print(paste("found ", number_hbr(hbs), " clusters", sep = ""))
                if (number_hbr(hbs) > 0) {
                    summi <- nrow(summary(hbs))
                    for (k in 1:summi) {
                        snpis <- data.table(cbind(chromo, c, k, markers(hbs, k)))
                        colnames(snpis) <- c("chr", "corr", "clust", "pos")
                        ids <- snpis[, .(chr, corr, clust)]
                        taggis <- apply(ids, 1, paste, collapse = "_")
                        snpis[, `:=`(tag, taggis)]
                        cluster.snps <- rbind(cluster.snps, snpis)
                    }
                }
            }
        }

        cluster.snps <- data.table(cluster.snps)
        check.cor.raw <- unique(cluster.snps[, .(chr, clust, corr)])
        check.cor.raw <- data.table(check.cor.raw)
        check.cor <- check.cor.raw[, .N, by = .(chr, corr)]

        for (a in chromis) {
            ## find corr to start
            check.cor.sub <- check.cor[chr != a & N != 0, corr]
            min.cl.cor <- as.numeric(check.cor[chr == a & corr %in% check.cor.sub, 
                                               corr[which.max(N)]])
            ## step parameters for raw and fine comparison
            ## rawsteps <- 0.05
            ## finesteps <- 0.01
            ## ## objects needed for comparison
            clusteron <- TRUE
            noT <- TRUE
            search <- FALSE
            raw <- TRUE
            succesful <- FALSE
            old.cl.cor <- 0
            countup <- 0

            
            while (clusteron) {
                alreadythere <- cluster.snps[corr == min.cl.cor]
                if (nrow(alreadythere) == 0) {
### get clusters per chromosome for focal correlation
                    for (b in chromis) {
                        ts <- ts.filter[[b]]
                        hbs <- reconstruct_hb(ts, chrom = b, min.cl.size = min.cl.size, 
                                              min.cl.cor = min.cl.cor, min.inter = min.inter, single.win = T, 
                                              transf = TRUE, arcsine = TRUE, scaleSNP = TRUE)
                                        # print(paste('corr',min.cl.cor))
                        print(paste("found ", number_hbr(hbs), " clusters", sep = ""))
                        if (number_hbr(hbs) > 0) {
                            summi <- nrow(summary(hbs))
                            for (k in 1:summi) {
                                snpis <- data.table(cbind(b, min.cl.cor, k, markers(hbs, k)))
                                colnames(snpis) <- c("chr", "corr", "clust", "pos")
                                ids <- snpis[, .(chr, corr, clust)]
                                taggis <- apply(ids, 1, paste, collapse = "_")
                                snpis[, `:=`(tag, taggis)]
                                cluster.snps <- rbind(cluster.snps, snpis)
                            }
                        }
                    }
                }            else{
                    print(paste("chr",a,"corr",unique(alreadythere$corr)))
                    print(paste("found ", length(unique(alreadythere$clust)), " clusters", sep = ""))
                }           
                ## ### clusters found on both chromosomes?
                ids <- cluster.snps[, .(chr, corr, clust)]
                taggis <- apply(ids, 1, paste, collapse = "_")
                cluster.snps[, `:=`(tag, taggis)]
                cluster.snps.sub <- cluster.snps[corr == min.cl.cor]
                check.cluster.snps.sub <- cluster.snps.sub[, .N, by = chr]
                if (nrow(check.cluster.snps.sub) > 1) {
                    cluster.snps.sub[, `:=`(pos, as.numeric(pos))]
                    cluster.ord <- cluster.snps.sub[order(as.numeric(pos)), 
                                                    .SD, by = .(chr, corr)]
                    ids <- cluster.ord[, .(chr, corr, clust)]
                    taggis <- apply(ids, 1, paste, collapse = "_")
                    cand.clust <- merge(cluster.ord, cands, by = c("chr", 
                                                                   "pos"))
                    clust.tags.foc <- unique(cand.clust[chr == a, tag])
                    clust.tags <- unique(cand.clust[, tag])
                    combi.tab <- c()
                    corfoc <- c()
                    corcomp <- c()
### compare focal tag to neighbouring tags in winsize + all tags from
### other chr
                    print(paste("corr ", min.cl.cor))
                    for (i in clust.tags.foc) {
                        cat(".")
                        cluster.sub.i.raw <- cand.clust[tag == i, ]
### reduce number of SNPs if too many per cluster
                        if (nrow(cluster.sub.i.raw) > takerandom) {
                            red.pos <- sample(cluster.sub.i.raw$pos, takerandom)
                            pos.indi <- cluster.sub.i.raw$pos %in% red.pos
                            cluster.sub.i <- cluster.sub.i.raw[pos.indi, ]
                        } else cluster.sub.i <- cluster.sub.i.raw
                        cluster.af.i <- cluster.sub.i[, afcolis, with = FALSE]
                        cluster.scale.i <- transform.af(cluster.af.i)
### find the neighbouring clusters in winsize
                        minmax <- cand.clust[tag == i, .(min(pos), max(pos))]
                        left <- na.exclude(unique(c(rev(cand.clust[chr == a & 
                                                                   pos >= minmax[, V1] - winsize & pos < minmax[, V2], 
                                                                   tag]))))
                        right <- na.exclude(unique(cand.clust[chr == a & pos > 
                                                              minmax[, V2] & pos <= minmax[, V2] + winsize, tag]))
                        neighbours <- unique(c(left[left != i], right))
### check if comparison is already, only do it once!
                        if (length(neighbours) > 0) {
                            combi.tab.sub <- cbind(i, neighbours)
                            if (length(combi.tab) > 0) {
                                check1 <- combi.tab.sub[, 1] %in% combi.tab[, 2]
                                check2 <- combi.tab.sub[, 2] %in% combi.tab[, 1]
                                red.indi <- apply(cbind(check1, check2), 1, all)
                                combis <- combi.tab.sub[red.indi == F, "neighbours"]
                                combi.tab <- rbind(combi.tab, combi.tab.sub[red.indi == 
                                                                            F, ])
                            } else {
                                combis <- combi.tab.sub[, "neighbours"]
                                combi.tab <- rbind(combi.tab, combi.tab.sub)
                            }
                        } else {
                            combis <- c()
                        }
### add clusters from other chr for comparison
                        compare.raw <- clust.tags[grep(paste("^", a, sep = ""), 
                                                       clust.tags, invert = T)]
                        compare.clust <- c(combis, compare.raw)
### 
                        clustcor <- sapply(compare.clust, comclust,cand.clust,takerandom,afcolis)
### separate intra and inter cluster comparison
                        corfoc <- c(corfoc,na.exclude(clustcor[compare.clust %in% combis]))
                        corcomp <- c(corcomp,na.exclude(clustcor[compare.clust %in% combis == F]))
### only perform when you have some values, fisher transform for
### normality
                    }
                    if (length(corfoc) > 3 & length(corcomp) > 3) {
                        t <- try(t.test(fisherz(corfoc), fisherz(corcomp), "greater"),silent=TRUE)
                        if (length(grep("Error",t))>0)
                            p <- NA
                        else
                            p <- t$p.value
                        print(paste(" pval ", round(p, 4), " for chr ", a, sep = ""))
### not if pvalue threshold was at least crossed once (because I want to
### see the threshold)
                        if (is.na(p))
                            p <- 1
                        old.cl.cor <- min.cl.cor
                        if (p <= thres.ttest){
                            if (raw) {
                                min.cl.cor <- min.cl.cor - rawsteps
                                noT <- F
                            } else {
                                min.cl.cor <- min.cl.cor - finesteps
                                raw <- FALSE
                                noT <- FALSE
                            }
                        }
                        else {
                            if (raw & noT == F) {
                                raw <- F
                                min.cl.cor <- min.cl.cor + rawsteps - finesteps
                            } else {
                                if (raw == F & noT == F) {
                                    print("Success!")
                                    clusteron <- FALSE
                                    successful <- TRUE
                                    cluster.final.sub<-  cluster.snps[corr == old.cl.cor &chr == a]
                                    final <- rbind(final, cluster.final.sub)
                                    print(paste("chr ", a, " done", sep = ""))
                                    print(final)
                                }
                            }
                        }
                        if (raw & noT) {
                            min.cl.cor <- min.cl.cor + rawsteps
                            countup <- countup + 1
                            if (countup > findthreshold )
                                search <- T
                        }
                    }
                    else{
                        if (raw) 
                            min.cl.cor <- min.cl.cor - rawsteps
                        else 
                            min.cl.cor <- min.cl.cor - finesteps
                    }
                }
                else{
                    if (raw) 
                        min.cl.cor <- min.cl.cor - rawsteps
                    else 
                        min.cl.cor <- min.cl.cor - finesteps
                }
                if (search) {
                    print(paste("chr ", a, " done; no  threshold!", sep = ""))
                    clusteron <- FALSE
                    cluster.final.sub<-  cluster.snps[corr == old.cl.cor &chr == a]
                    cluster.final.sub[,tag:=paste(tag,"_noThreshold",sep="")]
                    final <- rbind(final, cluster.final.sub)
                }
                ## else {
                ##     if (noT==F)
                ##         min.cl.cor <- min.cl.cor - finesteps
                ## }
                if (min.cl.cor <= 0.1 | min.cl.cor >= 1) {
                clusteron <- FALSE
                successful <- TRUE
                print(paste("chr ", a, " done; no  threshold!", sep = ""))
                cluster.final.sub<-  cluster.snps[corr == old.cl.cor &chr == a]
                    cluster.final.sub[,tag:=paste(tag,"_noThreshold",sep="")]
                    final <- rbind(final, cluster.final.sub)
                }
            }
        }
    }
    if (length(final)>0){
        final[,pos:=as.numeric(pos)]
        hapval.result[["all_haplotypes"]] <- final[,.(chr,pos,tag)]
        red <- list()
        cands.cmh <- merge(cands, cmh, by = c("chr", "pos"))
        datcmh <- merge(final, cands.cmh, by = c("chr", "pos"), all = T)
        datcmh.ord <- datcmh[order(score, decreasing = T)]
### check overlapping clusters ###
        for (i in na.exclude(unique(datcmh.ord$tag))) {
            x.sub <- final[tag == i]
            minmax <- x.sub[, .(min(pos), max(pos))]
            chr.sub <- unique(x.sub$chr)
            maxtag <- datcmh.ord[pos >= minmax[, V1] - filterrange & pos <= minmax[, 
                                                                                   V2] + filterrange & chr == chr.sub & is.na(tag) == F, tag[which.max(score)]]
            if (i != maxtag) 
                datcmh.ord[tag == i, `:=`(tag, NA)]
        }
        mergomat <- datcmh.ord[, .(chr, pos, score, tag)]
        mergomat.final <- mergomat[is.na(tag) == F]
        hapval.result[["dominant_haplotypes"]] <- mergomat.final[,.(chr,pos,tag)]
    }
    return(hapval.result)
}

