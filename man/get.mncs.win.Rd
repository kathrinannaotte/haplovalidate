\name{get.mncs.win}
\alias{get.mncs.win}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MNCS for window size
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{under construction
get.mncs.win(cands, cmh, wins, wincut)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cands}{
%%     ~~Describe \code{cands} here~~
}
  \item{cmh}{
%%     ~~Describe \code{cmh} here~~
}
  \item{wins}{
%%     ~~Describe \code{wins} here~~
}
  \item{wincut}{
%%     ~~Describe \code{wincut} here~~
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
function (cands, cmh, wins, wincut)
{
    snp.cmh <- merge(cands, cmh, by = c("chr", "pos"))
    scorewin <- c()
    for (k in wins) {
        snp.cmh[, `:=`(pos.bin, cut(pos, seq(0, max(pos), k *
            1e+06), include.lowest = T))]
        scorewin <- rbind(scorewin, (cbind(k, snp.cmh[, sum(as.numeric(score)),
            by = c("chr", "pos.bin")])))
    }
    scorewin <- data.table(scorewin)
    parameters <- c()
    for (c in unique(cands$chr)) {
        scorewin[chr == c, `:=`(score.sum, V1/sum(as.numeric(snp.cmh[chr ==
            c, score])))]
        cutoff <- scorewin[chr == c, median(score.sum), by = k]
        win <- cutoff[V1 >= wincut, k][1]
        parameters <- rbind(parameters, c(c, win))
    }
    colnames(parameters) <- c("chr", "win")
    parameters <- data.table(parameters)
    return(parameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
