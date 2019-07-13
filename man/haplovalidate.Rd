\name{haplovalidate}
\alias{haplovalidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validate reconstructed haplotypes}
%%  ~~function to do ... ~~

\description{under consrtruction
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
haplovalidate(cands, cmh, parameters, repl, gens, takerandom, filterrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cands}{
%%     ~~Describe \code{cands} here~~
}
  \item{cmh}{
%%     ~~Describe \code{cmh} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{repl}{
%%     ~~Describe \code{repl} here~~
}
  \item{gens}{
%%     ~~Describe \code{gens} here~~
}
  \item{takerandom}{
%%     ~~Describe \code{takerandom} here~~
}
  \item{filterrange}{
%%     ~~Describe \code{filterrange} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cands, cmh, parameters, repl, gens, takerandom, filterrange)
{
    base.pops <- c(rep(TRUE, length(repl)), rep(FALSE, length(repl) *
        (length(gens) - 1)))
    compare <- c(rep(rep(TRUE, length(gens)), length(repl)))
    min.minor.freq <- 0
    max.minor.freq <- 1
    minfreqchange <- 0
    minrepl <- 1
    min.lib.frac <- 0.75
    thres.ttest <- 0.025
    min.cl.size <- 20
    min.inter <- ceiling(min.cl.size/5)
    transform.af <- function(af) {
        af.sqrt <- asin(sqrt(af))
        af.transf <- t(af.sqrt)
        af.scale <- scale(af.transf, center = TRUE, scale = TRUE)
        return(af.scale)
    }
    chromis <- unique(cands$chr)
    final <- c()
    hapval.result <- list()
    if (length(chromis) > 1) {
        cluster.snps <- c()
        for (chromo in chromis) {
            winsize <- round(as.numeric(parameters[chr == chromo,
                win]) * 1e+06, 0)
            ts <- initialize_SNP_time_series(chr = cands[chr ==
                chromo, chr], pos = cands[chr == chromo, pos],
                base.freq = cands[chr == chromo, basePops], lib.freqs = cands[chr ==
                  chromo, 7:(ncol(cands)), with = F], pop.ident = c(rep(repl,
                  each = length(gens))), pop.generation = rep(gens,
                  length(repl)), use.libs = compare, min.minor.freq = min.minor.freq,
                max.minor.freq = max.minor.freq, winsize = winsize,
                win.scale = "bp", min.lib.frac = min.lib.frac,
                minfreqchange = minfreqchange, minrepl = minrepl)
            for (c in rev(seq(0.3, 0.9, 0.1))) {
                print(paste("chr ", chromo, " cluster corr ",
                  c, sep = ""))
                hbs <- reconstruct_hb(ts, chrom = chromo, min.cl.size = min.cl.size,
                  min.cl.cor = c, min.inter = min.inter, single.win = T,
                  transf = TRUE, arcsine = TRUE, scaleSNP = TRUE)
                print(paste("found ", number_hbr(hbs), " clusters",
                  sep = ""))
                if (number_hbr(hbs) > 0) {
                  summi <- nrow(summary(hbs))
                  for (k in 1:summi) {
                    snpis <- data.table(cbind(chromo, c, k, markers(hbs,
                      k)))
                    colnames(snpis) <- c("chr", "corr", "clust",
                      "pos")
                    ids <- snpis[, .(chr, corr, clust)]
                    taggis <- apply(ids, 1, paste, collapse = "_")
                    snpis[, `:=`(tag, taggis)]
                    cluster.snps <- rbind(cluster.snps, snpis)
                  }
                }
            }
        }
        cluster.snps <- data.table(cluster.snps)
        check.cor.raw <- unique(cluster.snps[, .(chr, clust,
            corr)])
        check.cor.raw <- data.table(check.cor.raw)
        check.cor <- check.cor.raw[, .N, by = .(chr, corr)]
        for (a in chromis) {
            check.cor.sub <- check.cor[chr != a & N != 0, corr]
            min.cl.cor <- as.numeric(check.cor[chr == a & corr \%in\%
                check.cor.sub, corr[which.max(N)]])
            rawsteps <- 0.05
            finesteps <- 0.01
            clusteron <- TRUE
            noT <- TRUE
            search <- FALSE
            raw <- TRUE
            succesful <- FALSE
            old.cl.cor <- 0
            while (clusteron) {
                alreadythere <- cluster.snps[corr == min.cl.cor]
                if (nrow(alreadythere) == 0) {
                  for (b in chromis) {
                    winsize <- round(as.numeric(parameters[chr ==
                      b, win]) * 1e+06, 0)
                    min.inter <- ceiling(min.cl.size/5)
                    ts <- initialize_SNP_time_series(chr = cands[chr ==
                      b, chr], pos = cands[chr == b, pos], base.freq = cands[chr ==
                      b, basePops], lib.freqs = cands[chr ==
                      b, 7:(ncol(cands)), with = F], pop.ident = c(rep(repl,
                      each = length(gens))), pop.generation = rep(gens,
                      length(repl)), use.libs = compare, min.minor.freq = min.minor.freq,
                      max.minor.freq = max.minor.freq, winsize = winsize,
                      win.scale = "bp", min.lib.frac = min.lib.frac,
                      minfreqchange = minfreqchange, minrepl = minrepl)
                    hbs <- reconstruct_hb(ts, chrom = b, min.cl.size = min.cl.size,
                      min.cl.cor = min.cl.cor, min.inter = min.inter,
                      single.win = T, transf = TRUE, arcsine = TRUE,
                      scaleSNP = TRUE)
                    print(paste("found ", number_hbr(hbs), " clusters",
                      sep = ""))
                    if (number_hbr(hbs) > 0) {
                      summi <- nrow(summary(hbs))
                      for (k in 1:summi) {
                        snpis <- data.table(cbind(b, min.cl.cor,
                          k, markers(hbs, k)))
                        colnames(snpis) <- c("chr", "corr", "clust",
                          "pos")
                        ids <- snpis[, .(chr, corr, clust)]
                        taggis <- apply(ids, 1, paste, collapse = "_")
                        snpis[, `:=`(tag, taggis)]
                        cluster.snps <- rbind(cluster.snps, snpis)
                      }
                    }
                  }
                }
                ids <- cluster.snps[, .(chr, corr, clust)]
                taggis <- apply(ids, 1, paste, collapse = "_")
                cluster.snps[, `:=`(tag, taggis)]
                cluster.snps.sub <- cluster.snps[corr == min.cl.cor]
                check.cluster.snps.sub <- cluster.snps.sub[,
                  .N, by = chr]
                if (nrow(check.cluster.snps.sub) > 1) {
                  cluster.snps.sub[, `:=`(pos, as.numeric(pos))]
                  cluster.ord <- cluster.snps.sub[order(as.numeric(pos)),
                    .SD, by = .(chr, corr)]
                  ids <- cluster.ord[, .(chr, corr, clust)]
                  taggis <- apply(ids, 1, paste, collapse = "_")
                  cand.clust <- merge(cluster.ord, cands, by = c("chr",
                    "pos"))
                  afcolis <- grep("L[0-9]", colnames(cand.clust),
                    value = TRUE)
                  clust.tags.foc <- unique(cand.clust[chr ==
                    a, tag])
                  clust.tags <- unique(cand.clust[, tag])
                  combi.tab <- c()
                  corfoc <- c()
                  corcomp <- c()
                  print(paste("corr ", min.cl.cor))
                  for (i in clust.tags.foc) {
                    cat(".")
                    cluster.sub.i.raw <- cand.clust[tag == i,
                      ]
                    if (nrow(cluster.sub.i.raw) > takerandom) {
                      red.pos <- sample(cluster.sub.i.raw$pos,
                        takerandom)
                      pos.indi <- cluster.sub.i.raw$pos \%in\%
                        red.pos
                      cluster.sub.i <- cluster.sub.i.raw[pos.indi,
                        ]
                    }
                    else cluster.sub.i <- cluster.sub.i.raw
                    cluster.af.i <- cluster.sub.i[, afcolis,
                      with = FALSE]
                    cluster.scale.i <- transform.af(cluster.af.i)
                    minmax <- cand.clust[tag == i, .(min(pos),
                      max(pos))]
                    left <- na.exclude(unique(c(rev(cand.clust[chr ==
                      a & pos >= minmax[, V1] - winsize & pos <
                      minmax[, V2], tag]))))
                    right <- na.exclude(unique(cand.clust[chr ==
                      a & pos > minmax[, V2] & pos <= minmax[,
                      V2] + winsize, tag]))
                    neighbours <- unique(c(left[left != i], right))
                    if (length(neighbours) > 0) {
                      combi.tab.sub <- cbind(i, neighbours)
                      if (length(combi.tab) > 0) {
                        check1 <- combi.tab.sub[, 1] \%in\% combi.tab[,
                          2]
                        check2 <- combi.tab.sub[, 2] \%in\% combi.tab[,
                          1]
                        red.indi <- apply(cbind(check1, check2),
                          1, all)
                        combis <- combi.tab.sub[red.indi == F,
                          "neighbours"]
                        combi.tab <- rbind(combi.tab, combi.tab.sub[red.indi ==
                          F, ])
                      }
                      else {
                        combis <- combi.tab.sub[, "neighbours"]
                        combi.tab <- rbind(combi.tab, combi.tab.sub)
                      }
                    }
                    else {
                      combis <- c()
                    }
                    compare.raw <- clust.tags[grep(paste("^",
                      a, sep = ""), clust.tags, invert = T)]
                    compare.clust <- c(combis, compare.raw)
                    comclust <- function(x) {
                      sub <- cand.clust[tag == x, ]
                      if (nrow(sub) > takerandom) {
                        red.pos <- sample(sub$pos, takerandom)
                        pos.indi <- sub$pos \%in\% red.pos
                        sub <- sub[pos.indi, ]
                      }
                      cluster.af <- sub[, afcolis, with = FALSE]
                      cluster.scale <- transform.af(cluster.af)
                      clustcor.sub <- median(cor(cluster.scale.i,
                        cluster.scale), na.rm = TRUE)
                      return(clustcor.sub)
                    }
                    clustcor <- sapply(compare.clust, comclust)
                    corfoc <- c(corfoc, clustcor[compare.clust \%in\%
                      combis])
                    corcomp <- c(corcomp, clustcor[compare.clust \%in\%
                      combis == F])
                  }
                }
                if (length(corfoc) > 3 & length(corcomp) > 3) {
                  t <- t.test(fisherz(corfoc), fisherz(corcomp),
                    "greater")
                  p <- t$p.value
                  print(paste(" pval ", round(p, 4), " for chr ",
                    a, sep = ""))
                  if (p <= thres.ttest) {
                    old.cl.cor <- min.cl.cor
                    if (raw & search == F) {
                      min.cl.cor <- min.cl.cor - rawsteps
                      noT <- F
                    }
                    else {
                      min.cl.cor <- min.cl.cor - finesteps
                      raw <- FALSE
                      noT <- F
                    }
                  }
                  else {
                    if (raw & noT == F) {
                      raw <- F
                      min.cl.cor <- min.cl.cor + rawsteps - finesteps
                    }
                    else {
                      if (raw == F & noT == F) {
                        print("Normal")
                        clusteron <- FALSE
                        successful <- TRUE
                        cluster.final.sub <- cluster.snps[corr ==
                          old.cl.cor & chr == a]
                        final <- rbind(final, cluster.final.sub)
                        print(paste("chr ", a, " done", sep = ""))
                      }
                    }
                    if (raw & noT) {
                      min.cl.cor <- min.cl.cor + finesteps
                      search <- T
                    }
                  }
                }
                else {
                  if (search) {
                    print("search")
                    clusteron <- FALSE
                    successful <- TRUE
                    final.raw <- cluster.snps[corr == old.cl.cor &
                      chr == a, `:=`(tag, paste(tag, "_search",
                      sep = ""))]
                    final <- rbind(final, final.raw)
                    print(paste("chr ", a, " done", sep = ""))
                  }
                  else {
                    if (noT)
                      min.cl.cor <- min.cl.cor + finesteps
                    else min.cl.cor <- min.cl.cor - finesteps
                  }
                }
                if (min.cl.cor < 0.1 | min.cl.cor > 1) {
                  clusteron <- FALSE
                  successful <- TRUE
                  final.raw <- cluster.snps[corr == old.cl.cor &
                    chr == a, `:=`(tag, paste(tag, "_no_conv",
                    sep = ""))]
                  final <- rbind(final, final.raw[chr == a])
                  print(paste("chr ", a, " done; no  convergence!",
                    sep = ""))
                }
            }
        }
    }
    final[, `:=`(pos, as.numeric(pos))]
    hapval.result[["all_haplotypes"]] <- final
    red <- list()
    if (final != "clustering impossible") {
        corris <- unique(final[, .(chr, corr)])[, corr]
        finish <- unique(str_extract(final[, tag], "no_.*"))
        if (length(finish) == 1)
            finish <- rep(finish, 2)
        cands.cmh <- merge(cands, cmh, by = c("chr", "pos"))
        datcmh <- merge(final, cands.cmh, by = c("chr", "pos"),
            all = T)
        datcmh.ord <- datcmh[order(score, decreasing = T)]
        for (i in na.exclude(unique(datcmh.ord$tag))) {
            x.sub <- final[tag == i]
            minmax <- x.sub[, .(min(pos), max(pos))]
            chr.sub <- unique(x.sub$chr)
            maxtag <- datcmh.ord[pos >= minmax[, V1] - filterrange &
                pos <= minmax[, V2] + filterrange & chr == chr.sub &
                is.na(tag) == F, tag[which.max(score)]]
            if (i != maxtag)
                datcmh.ord[tag == i, `:=`(tag, NA)]
        }
        mergomat <- datcmh.ord[, .(chr, pos, score, tag)]
        mergomat.final <- mergomat[is.na(tag) == F]
        hapval.result[["dominant_haplotypes"]] <- mergomat.final
    }
    return(hapval.result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
