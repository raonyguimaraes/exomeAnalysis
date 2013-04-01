d <- read.table("/Users/singhal/thesisWork/introgression/clineAndSummaryStats.out",header=T)

d<-d[complete.cases(d$type),]
d<-d[d$type=="widerange",]

dist <- read.table("/Users/singhal/thesisWork/introgression/distances/distances",header=T)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')

c <- split(d,d$contact)
contacts <- c("sjo","carlia","gillies","nBMB")

for (i in 3:3) {
	af_file <- paste("/Users/singhal/thesisWork/introgression/clineAF/",contacts[i],".cline.out",sep="")
	af <- read.table(af_file,header=F,stringsAsFactors=F)
	names(af) <- c("locus","pos","pop","af")

	#get rid of control loci
	af<-af[grep('ND4',af$locus,invert=T),]
	af<-af[grep('16s',af$locus,invert=T),]
	af<-af[grep('utr',af$locus,invert=T),]
	
	#add distance
	dist_contact = dist[dist$contact==contacts[i],]
	af <- data.frame(af,rep(NA,dim(af)[1]))
	names(af)[5] <- c("dist")
	for (p in 1:length(pops)) {
		af[af$pop == pops[p],]$dist = dist_contact[dist_contact$library==pops[p],]$distance
		}
	
	plot(NULL,xlim=c(-40000,40000),ylim=c(-0.1,1.1))
	loci <- c[[contacts[i]]]
	for (j in 1:dim(loci)[1]) {
		contig <- loci[j,]$contig
		pos <- loci[j,]$pos
		tmp <- af[af$locus==contig & af$pos == pos,]
		
		line <- lm(tmp[2:10,]$af~tmp[2:10,]$dist)
		b <- line$coefficients[1]
		m <- line$coefficients[2]
		
		xstart <- min(tmp[2:10,]$dist)
		xend <- max(tmp[2:10,]$dist)
		ystart <- m * xstart + b
		yend <- m * xend + b
		
		#if (ystart < 0.3 | ystart > 0.7) {
		#	if (yend > 0.7 | yend < 0.3) {			
				segments(tmp[1,]$dist,tmp[1,]$af,x1=xstart,y1=ystart,col=rgb(0,0,0,.15,maxColorValue=1))
				segments(xstart,ystart,x1=xend,y1=yend,col=rgb(0,0,0,.15,maxColorValue=1))
				segments(xend,yend,x1=tmp[11,]$dist,y1=tmp[11,]$af,col=rgb(0,0,0,.15,maxColorValue=1))
				}
		#	}	
		#}		
	}
	
	
	
}