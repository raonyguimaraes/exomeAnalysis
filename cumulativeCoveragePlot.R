contacts = c("nBMB","sjo","carlia","gillies")
par(mfrow=c(1,4),mar=c(5,4,3,5))
for (i in 1:length(contacts)) {

	file <- paste("/Users/singhal/Desktop/activeWork/prelimClines/",contacts[i],".depth.summarized",sep="")
	d <- read.table(file,header=F)

	coverage <- d$V2 + d$V3 + d$V4 + d$V5 + d$V6 + d$V7 + d$V8 + d$V9 + d$V10
	coverage <- coverage/9
	coverage <- coverage[coverage < 900]

	# create vector of number of coverages greater than coverage 0X to (max.coverage+1)X
	cumulative <- c(NA)

	for (i in 1:(ceiling(max(coverage))+1)) {
		cumulative[i] <- length(which(coverage >= i-1))
		i <- i+1
	}
	
	cumulative.freq <- (cumulative[1:length(cumulative)])/length(coverage) # Make vector elements proportions

	scale = 1/( max(density(coverage)$y)*1.1 )
	cumulative.freq.scale <- cumulative.freq/scale #scale for plotting

	cover <- c(0:(ceiling(max(coverage)))) # vector of coverages ranging from 0 to avg. max coverage rounded up 

	#plot, must plot density first
	plot(NULL,xlim=range(coverage),ylim=range	(cumulative.freq.scale),xlab="Coverage",ylab="Percent contigs",main=NA,axes=T) 
	polygon(density(coverage),col="gray",border="gray")
	lines(cover,cumulative.freq.scale,type="l",lty=2,col=1)
	mtext("cumulative % coverage",side=4,line=2.9,cex=0.8)
	ymax = max(cumulative.freq.scale)
	axis(side=4,at=c(0 * ymax ,0.2 * ymax, 0.4 * ymax, 0.6 * ymax,0.8 * ymax, 1 * 	ymax),labels=c(0,20,40,60,80,100))
	box(which="plot")
	}