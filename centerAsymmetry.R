library(moments)

fit <- read.table("/Users/singhal/thesisWork/introgression/clineData",header=T)
fit <- fit[complete.cases(fit$type),]
fit <- fit[fit$type=="cline",]

contact = c("nBMB","sjo","carlia","gillies")

par(mfrow=c(2,2),mai=c(0.8,0.8,0.4,0.4))
for (i in 1:length(contact)) {
	tmp <- fit[fit$contact==contact[i],]
	tmp <- tmp[complete.cases(tmp$center),]
	
	dens <- density(tmp$center)$y / sum( density(tmp$center)$y )
	
	tot = 0
	for (j in 1:length(dens)) {
		tot = tot + dens[j]
		if (tot >= 0.05) {
			xmin <- density(tmp$center)$x[j]
			break
			}				
		}
		
	tot = 0
	for (j in 1:length(dens)) {
		tot = tot + dens[j]
		if (tot >= 0.95) {
			xmax <- density(tmp$center)$x[j]
			break
			}				
		}			
		
	center <- tmp$center
	center <- center[center <= xmax]
	center <- center[center >= xmin]
	center <- center - median(center)
	
	numBreaks = length(center)/50

	a <- hist(center,breaks=numBreaks,col="gray",border=F,main=paste(contact[i]),freq=F)	
	a

	med <- median(center)
	mean <- mean(center)
	agost <- agostino.test(center,alternative="two.sided")
	
	abline(v=mean,lty=3,col="black")
	text(x=max(center)*0.7, y=max(a$density)*0.8, labels=paste("skew: ",round(agost$statistic[1],2),sep=""))
	text(x=max(center)*0.7, y=max(a$density)*0.74, labels=paste("p-val < 1e-20"))
}
