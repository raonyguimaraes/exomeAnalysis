d <- read.table("/Users/singhal/thesisWork/introgression/LD/moransI.out",header=F)
names(d) <- c("contact","dist","moransI")

contacts <- c("sjo","nBMB","carlia", "gillies")
colors <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")

d[d$dist == 0,]$dist = 1.1
plot(NULL,xlim=range(d$dist),ylim=range(d$moransI),xlab="distance (bps)", ylab="Moran's I")

for (i in 1:length(contacts)) {
	tmp <- d[d$contact==contacts[i],]
	lines(tmp$dist,tmp$moransI,col=colors[i],lwd=2)	
}
legend(x=700,y=0.55,legend=contacts,fill=colors,col=colors,bty="n")