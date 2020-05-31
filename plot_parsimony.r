### (c) rucorbet@ucsc.edu
### this script is distributed freely but with no warranty
### 

### USAGE: Rscript --vanilla plot_extremal.r [extremal_file] [output pdf] 

args <- commandArgs(TRUE)
srcFile <- args[1]
outFile <- args[2] 

require(plyr)

pdf(outFile,height=12,width=8)
par(mfrow=c(3,1))

d <- as.data.frame(read.table(srcFile))

## generate counts of sites by parsimony
counts_all <- ddply(d, .(d$V2, d$V3), nrow)
names(counts_all) <- c("alt", "parsimony", "Freq")

## plot for all sites 
plot(counts_all$alt,counts_all$parsimony,log="x",pch=19,cex=0,xlab="Alternate allele count",ylab="Parsimony score")
for ( i in 1:length(counts_all$alt)+1 ) {
         points(counts_all$alt,counts_all$parsimony,pch=19,cex=c(2*log10(counts_all$Freq+1)),col=rgb(0,0,0,0.005))
}

## grab extremal
extremal <- d[d$V7==1,]

## add the linear fit to the graph
fit <- lm(extremal$V3 ~ log10( extremal$V2 ) )
abline( a = fit$coefficients[1], b = fit$coefficients[2], lty= 3 )

## plot extremal points
points( extremal$V2, extremal$V3, col="red", pch=18)

## label the slope
text(labels=c("slope (all)=\n\n", fit$coefficients[2]/3.32116635809), 10,max(d$V3)-2)

## again but no C>U
x <- d[d$V4 == 0,]
counts <- ddply(x, .(x$V2, x$V3), nrow)
names(counts) <- c("alt", "parsimony", "Freq")

## plot 
plot(counts_all$alt,counts_all$parsimony,log="x",pch=19,cex=0,xlab="Alternate allele count",ylab="Parsimony score")
for ( i in 1:length(counts$alt)+1 ) {
         points(counts$alt,counts$parsimony,pch=19,cex=c(2*log10(counts$Freq+1)),col=rgb(0,0,0,0.005))
}

## grab extremal
extremal <- x[x$V9==1,]

## add the linear fit to the graph
fit <- lm(extremal$V3 ~ log10( extremal$V2 ) )
abline( a = fit$coefficients[1], b = fit$coefficients[2], lty= 3 )

## plot extremal points
points( extremal$V2, extremal$V3, col="red", pch=18)

## label the slope
text(labels=c("slope (no C>U)=\n\n", fit$coefficients[2]/3.32116635809), 10,max(x$V3)-2)



## again but no C>U or G>U
x <- d[d$V4 == 0 & d$V5 == 0,]
counts <- ddply(x, .(x$V2, x$V3), nrow)
names(counts) <- c("alt", "parsimony", "Freq")

## plot
plot(counts_all$alt,counts_all$parsimony,log="x",pch=19,cex=0,xlab="Alternate allele count",ylab="Parsimony score")
for ( i in 1:length(counts$alt)+1 ) {
         points(counts$alt,counts$parsimony,pch=19,cex=c(2*log10(counts$Freq+1)),col=rgb(0,0,0,0.005))
}

## grab extremal
extremal <- x[x$V9==1,]

## add the linear fit to the graph
fit <- lm(extremal$V3 ~ log10( extremal$V2 ) )
abline( a = fit$coefficients[1], b = fit$coefficients[2], lty= 3 )

## plot extremal points
points( extremal$V2, extremal$V3, col="red", pch=18)

## label the slope
text(labels=c("slope (no C>U, no G>U)=\n\n", fit$coefficients[2]/3.32116635809), 10,max(x$V3)-2)
dev.off()
