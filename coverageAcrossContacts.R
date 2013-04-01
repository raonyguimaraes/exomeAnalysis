contacts = c("nBMB","sjo","carlia","gillies")
allcov <- data.frame()

for (i in 1:length(contacts)) {

	file <- paste("/Users/singhal/thesisWork/introgression/coverage/",contacts[i],".depth.summarized",sep="")
	d <- read.table(file,header=F)

	coverage <- d$V2 + d$V3 + d$V4 + d$V5 + d$V6 + d$V7 + d$V8 + d$V9 + d$V10
	coverage <- coverage/9
	
	loci <- d$V1
	loci <- sub("_\\d+","",loci,perl=T)	
	loci <- sub("_\\d+","",loci,perl=T)	
		
	tmp <- data.frame(loci,coverage)
	names(tmp) <- c("locus","cov")
	tmp <- aggregate(tmp$cov,by=list(tmp$locus),mean)
	names(tmp) <- c("locus",paste("cov",contacts[i],sep="_"))
	
	if (i == 1) {
		allcov = tmp
		}
	else {
		allcov = merge(allcov,tmp,by.x="locus",by.y="locus",all=T)
		}
	}

par(mfrow=c(2,2))
for (i in 2:5) {
	j = i + 1
	if (i == 5) { j = 2 }
	plot(allcov[,i],allcov[,j],xlim=c(0,2500),ylim=c(0,2500),pch=16,cex=0.5,xlab=paste("coverage of",contacts[i -1],sep=" "), ylab=paste("coverage of",contacts[j-1],sep=" "))
	abline(a=0,b=1,lty=2,col=2)
	corr <- round(cor(allcov[,i],allcov[,j],use="complete.obs"),3)
	text(500,2300,paste("r",corr,sep="="))	
}	
	
	