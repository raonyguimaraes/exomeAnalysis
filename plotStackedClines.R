clines <- 100
contacts <- c("nBMB","sjo","carlia","gillies")

fit <- read.table("/Users/singhal/thesisWork/introgression/clineData",header=T)
fit <- fit[complete.cases(fit$type),]
fit <- fit[fit$type=="cline",]

widthfit <- fit
centerfit <- fit

centerfit <- centerfit[complete.cases(centerfit$center),]
centerfit <- centerfit[centerfit$center >= 0,]
centerfit <- centerfit[centerfit$center <= 25000,]
centerfit <- data.frame(centerfit, rep(NA, dim(centerfit)[1]))
names(centerfit)[length(centerfit)] <- c("center")
for (i in 1:length(contacts)) {
	centermedian <- median(centerfit[centerfit$contact==contacts[i],]$center)
	centerfit[centerfit$contact==contacts[i],]$center = centerfit[centerfit$contact==contacts[i],]$center	- centermedian	
}

widthfit <- widthfit[complete.cases(widthfit$width),]
widthfit <- widthfit[widthfit$width >= 0,]
widthfit <- widthfit[widthfit$width <= 100000,]

quartz("A", 3,7)
par(mfrow=c(4,1),mai=c(0.4,0.4,0.1,0.1),oma=c(2,2,0,0))
for (i in 1:length(contacts)) {
	plot(NULL,xlim=c(0,25e3),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n")
	axis(1,at=c(2.5e3,7.5e3,12.5e3,17.5e3,22.5e3),labels=c(-10e3,-5e3,0,5e3,10e3))
	axis(2,at=c(0,0.25,0.5,0.75,1),labels=c(0,0.25,0.5,0.75,1))
	tmp <- widthfit[widthfit$contact == contacts[i],]
	tmp <- tmp[sample(nrow(tmp),clines*1.5),]
	count = 0
	
	for (n in 1:dim(tmp)[1]) {
		if (count < clines) {
			a <- tmp[n,]
			b <- centerfit[centerfit$name==a$name & centerfit$position==a$position,]
			if (dim(b)[1]  == 1) {
				center <- b$center + 12.5e3
				if (center > 0) {
					x <- seq(0,25e3,by=50)
					width <- a$width
					fit <- ((1+tanh(2*(x-center)/width))/2)
					lines(x,fit,col=rgb(0,0,0,.15,maxColorValue=1))
					count = count + 1
					}						 
				}	
			}		
		}	
	}	
mtext("distance (m)",side=1,outer=T)	
mtext("allele frequency",side=2,outer=T)	
quartz.save("/Users/singhal/Desktop/cline.pdf",type="pdf")