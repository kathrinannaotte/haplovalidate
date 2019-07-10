

get.mncs.win <- function(cands,cmh,wins,wincut){
    snp.cmh <- merge(cands,cmh,by=c("chr","pos"))
    scorewin <- c()
    for (k in wins){
        snp.cmh[,pos.bin := cut(pos,seq(0,max(pos),k * 1000000),include.lowest=T)]
        scorewin <-rbind(scorewin,(cbind(k,snp.cmh[,sum(as.numeric(score)),by=c("chr","pos.bin")])))
    }
    scorewin <- data.table(scorewin)
    parameters <- c()
    for (c in unique(cands$chr)){
        scorewin[chr==c,score.sum:=V1/sum(as.numeric(snp.cmh[chr==c,score]))]
        cutoff<- scorewin[chr==c,median(score.sum),by=k]
        win <- cutoff[V1>=wincut,k][1]
        parameters <- rbind(parameters,c(c,win))
    }
    colnames(parameters) <- c("chr","win")
    parameters <- data.table(parameters)
    return(parameters)
}


