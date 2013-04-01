contacts = c("nBMB","sjo","carlia","gillies")

par(mfrow=c(2,2))
for (i in 1:length(contacts)) {

	file1 <- paste("/Users/singhal/thesisWork/introgression/coverage/",contacts[i],".depth.summarized",sep="")
	d <- read.table(file1,header=F)

	coverage <- d$V2 + d$V3 + d$V4 + d$V5 + d$V6 + d$V7 + d$V8 + d$V9 + d$V10
	coverage <- coverage/9
	tmp <- data.frame(d$V1,coverage)
	names(tmp) <- c("locus","cov")
	
	file2 = paste("/Users/singhal/Desktop/",contacts[i],"_divpoly.out",sep="")
	d <- read.table(file2,header=T)
	
	new = merge(tmp,d,by.x="locus",by.y="contig",all=T)

	plot(new$net,new$cov,pch=16,cex=0.5,xlab="divergence", ylab="coverage",ylim=c(0,2500))
	corr <- round(cor(new$net,new$cov,use="complete.obs"),3)
	text(range(new$net,na.rm=T)[2]*0.9,2300,paste("r",corr,sep="="))	
	title(paste(contacts[i]))
	}

