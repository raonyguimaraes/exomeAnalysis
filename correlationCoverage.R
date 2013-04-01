d <- read.table("/Users/singhal/Desktop/activeWork/prelimClines/carlia.depth.summarized",header=F)

par(mfrow=c(1,4),oma=c(0,0,0,0),mai=c(0.6,0.6,0.2,0.2))

plot(d$V2,d$V3,xlab="coverage, library 1", ylab="coverage, library 2",xlim=c(0,4500),ylim=c(0,4500),pch=16,cex=0.5)
abline(a=0,b=1,lty=2,col=2)
corr <- round(cor(d$V2,d$V3),3)
text(500,4200,paste("r",corr,sep="="))

plot(d$V4,d$V5,xlab="coverage, library 1", ylab="coverage, library 2",xlim=c(0,4500),ylim=c(0,4500),pch=16,cex=0.5)
abline(a=0,b=1,lty=2,col=2)
corr <- round(cor(d$V4,d$V5),3)
text(500,4200,paste("r",corr,sep="="))

plot(d$V6,d$V7,xlab="coverage, library 1", ylab="coverage, library 2",xlim=c(0,4500),ylim=c(0,4500),pch=16,cex=0.5)
abline(a=0,b=1,lty=2,col=2)
corr <- round(cor(d$V6,d$V7),3)
text(500,4200,paste("r",corr,sep="="))

plot(d$V8,d$V9,xlab="coverage, library 1", ylab="coverage, library 2",xlim=c(0,4500),ylim=c(0,4500),pch=16,cex=0.5)
abline(a=0,b=1,lty=2,col=2)
corr <- round(cor(d$V8,d$V9),3)
text(500,4200,paste("r",corr,sep="="))
