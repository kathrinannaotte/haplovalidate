\name{plot.haplovalidate}
\alias{plot.haplovalidate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make manhattanplot highlighting reconstructed haplotype blocks}
%%  ~~function to do ... ~~

\description{under construction
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.haplovalidate(blocks, cmh, label = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blocks}{
%%     ~~Describe \code{blocks} here~~
}
  \item{cmh}{
%%     ~~Describe \code{cmh} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
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
function (blocks, cmh, label = TRUE)
{
    neuefarbe.raw <- unique(c(brewer.pal(8, "Set1"), brewer.pal(7,
        "Dark2"), brewer.pal(7, "Set2"), brewer.pal(12, "Set3"),
        brewer.pal(12, "Paired"), brewer.pal(7, "Accent"), brewer.pal(7,
            "Spectral")))
    neuefarbe <- colorRampPalette(neuefarbe.raw)(max(as.numeric(factor(blocks$tag)),
        na.rm = TRUE))
    col.sub <- neuefarbe[as.numeric(factor(blocks$tag))]
    blocks[, `:=`(groupcol, col.sub)]
    for (j in unique(blocks$chr)) {
        cluster.sub <- blocks[chr == j, ]
        png(paste(j, "_haplovalidate.png", sep = ""), 2400, 1200)
        maxy <- round(max(as.numeric(cmh[chr == j, score]), na.rm = TRUE),
            -1) + round(max(as.numeric(cmh[chr == j, score]),
            na.rm = TRUE), -1) * 0.05
        cmh[chr == j, plot(pos/1e+06, as.numeric(score), pch = 19,
            col = "#414547", ylim = c(0, maxy), xlab = "position Mb",
            ylab = "cmh score", lwd = 3, cex.axis = 1.5, cex.lab = 2,
            main = paste("sim ", no, " merge chr", j, sep = ""),
            xaxp = c(0, 50, 50))]
        cluster.sub[, points(as.numeric(pos)/1e+06, as.numeric(score),
            cex = 2, lwd = 3, pch = 19, col = groupcol)]
        count <- 0
        cluster.sub.ord <- cluster.sub[order(score, decreasing = TRUE),
            ]
        for (k in na.exclude(unique(cluster.sub.ord$groupcol))) {
            group.sub <- cluster.sub.ord[groupcol == k, ]
            distance <- group.sub$pos/1e+06
            levi <- maxy - count
            points(distance, rep(levi, length(distance)), col = k,
                pch = 19, cex = 3)
            if (label)
                text(distance[1], levi, unique(group.sub[, tag]),
                  cex = 2)
            count <- count + round(max(as.numeric(cmh[chr ==
                j, score]), na.rm = TRUE), -1) * 0.01
        }
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
