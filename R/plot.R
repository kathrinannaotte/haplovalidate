#' @ import RColorBrewer brewer.pal colorRampPalette

#' @export
 plot.haplovalidate<- function(blocks,cmh,title="",label=TRUE){
     neuefarbe.raw <- unique(c(brewer.pal(8,"Set1"),brewer.pal(7,"Dark2"),brewer.pal(7,"Set2"),brewer.pal(12,"Set3"),brewer.pal(12,"Paired"),brewer.pal(7,"Accent"),brewer.pal(7,"Spectral")))
     neuefarbe  <- colorRampPalette(neuefarbe.raw)(max(as.numeric(factor(blocks$tag)),na.rm=TRUE))
    ## as.numeric removes the non-grouped tags!
    col.sub <- neuefarbe[as.numeric(factor(blocks$tag))]
    blocks[,groupcol := col.sub]
        for (j in unique(blocks$chr)){
            cluster.sub <- merge(blocks[chr==j,],cmh,by=c("chr","pos"))
            png(paste(title,"_",j,"_haplovalidate.png",sep=""),2400,1200)
            maxy <- round(max(as.numeric(cmh[chr==j,score]),na.rm=TRUE),0) + round(max(as.numeric(cmh[chr==j,score]),na.rm=TRUE),0) * 0.02
            cmh[chr==j,plot(pos/1000000,as.numeric(score),pch=19,col="#414547",ylim=c(0,maxy),xlab="position Mb",ylab="cmh score",lwd=3,cex.axis=1.5,cex.lab=2,main=title,xaxp=c(0,50,50))]
            cluster.sub[,points(as.numeric(pos)/1000000,as.numeric(score),cex=2,lwd=3,pch=19,col=groupcol)]
            count <- 0
            cluster.sub.ord <- cluster.sub[order(score,decreasing=TRUE),]
            for (k in na.exclude(unique(cluster.sub.ord$groupcol))){
                group.sub <- cluster.sub.ord[groupcol==k,]
                distance <- group.sub$pos/1000000
                levi <- maxy  - count 
                points(distance,rep(levi,length(distance)),col=k,pch=19,cex=3)
                if (label)
                    text(distance[1],levi,unique(group.sub[,tag]),cex =2)
                count <- count +  round(max(as.numeric(cmh[chr==j,score]),na.rm=TRUE),0) * 0.01
            }
            dev.off()
        }
     return(blocks)
 }

