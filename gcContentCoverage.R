library(seqinr)
contacts = c("nBMB","sjo","carlia","gillies")

par(mfrow=c(2,2))
for (i in 1:length(contacts)) {

	file <- paste("/Users/singhal/Desktop/activeWork/prelimClines/",contacts[i],".depth.summarized",sep="")
	d <- read.table(file,header=F)

	coverage <- d$V2 + d$V3 + d$V4 + d$V5 + d$V6 + d$V7 + d$V8 + d$V9 + d$V10
	coverage <- coverage/9
			
	a <- data.frame(d$V1,coverage)
	names(a) <- c("locus","cov")
	
	gc <- c()
	seqfile <- paste("/Users/singhal/thesisWork/introgression/targetSequences/final/",contacts[i], "_targets.fa.final",sep="")
	seq <- read.fasta(file=seqfile)
	for (j in 1:length(seq)) {
		l <- length(seq[[j]])
		c <- length(grep("c",seq[[j]]))
		g <- length(grep("g",seq[[j]]))
		cg <- (c+g)/l
		gc[j] <- cg
		}
	b <- data.frame(names(seq),gc)
	names(b) <- c("locus","gc")
	
	new = merge(a,b,by.x="locus",by.y="locus",all=T)
	names(new) <- c("locus","cov","gc")
	plot(new$gc,new$cov,xlim=c(0,1),ylim=c(0,2500),pch=16,cex=0.5,xlab="GC-content", ylab="coverage")
	title(paste(contacts[i]))
	}
