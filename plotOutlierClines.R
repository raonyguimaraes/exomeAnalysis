clines <- 5

fit <- read.table("/Users/singhal/thesisWork/introgression/clineAndSummaryStats.out",header=T)
fit$outtype <- factor(fit$outtype, levels = c(levels(fit$outtype), "sweep"))
indices <- which(fit$type == "widerange")
fit$outtype = replace(fit$outtype,indices,"sweep")

contact <- "gillies"

fit <- fit[complete.cases(fit$outtype),]
fit <- fit[fit$contact==contact,]

types <- c("sweep","narrow","normal","wide")
colors <- c("#00BFC4","#F8766D","gray","#C77CFF")

mean_center = mean(fit$center,na.rm=T)

x <- seq(0,10e3,by=50)
plot(NULL,xlim=range(x),ylim=c(-0.1,1.1),xaxt='n',xlab="distance (m)",ylab="allele freq.")
axis(1,at=c(1000,3000,5000,7000,9000),labels=c(-4000,-2000,0,2000,4000))

for (i in 2:length(types)) {
	tmp <- fit[fit$outtype==types[i],]	
	count = 0	
	while (count < clines) {
		n = round(runif(1,min=0,max=dim(tmp)[1]))
		a <- tmp[n,]
		center <- 5000 #this is so all clines have the same center
		width <- a$width
		xfit <- ((1+tanh(2*(x-center)/width))/2)
		lines(x,xfit,col=colors[i])
		count = count + 1		
		}	
	}

dist <- read.table("/Users/singhal/thesisWork/introgression/distances/distances",header=T)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')

for (i in 1:1) {
	tmp <- fit[fit$outtype==types[i],]	
	count = 0	
	
	af_file <- paste("/Users/singhal/thesisWork/introgression/clineAF/", contact, ".cline.out",sep="")	
	af <- read.table(af_file,header=F,stringsAsFactors=F)
	names(af) <- c("locus","pos","pop","af")
	
	#add distance
	dist_contact = dist[dist$contact==contact,]
	af <- data.frame(af,rep(NA,dim(af)[1]))
	names(af)[5] <- c("dist")
	for (p in 1:length(pops)) {
		af[af$pop == pops[p],]$dist = dist_contact[dist_contact$library==pops[p],]$distance
		}
	
	while (count < clines) {
		n = round(runif(1,min=0,max=dim(tmp)[1]))
		a <- tmp[n,]
		contig <- a$contig
		pos <- a$pos
		
		
		
		tmp_af <- af[af$locus==contig & af$pos == pos,]	
		line <- lm(tmp_af[2:10,]$af~tmp_af[2:10,]$dist)
		b <- line$coefficients[1]
		m <- line$coefficients[2]
		
		xstart <- min(x)
		xend <- max(x)
		ystart <- m * xstart + b
		yend <- m * xend + b
		
		if (ystart < 0.3 | ystart > 0.7) {
			if (yend > 0.7 | yend < 0.3) {			
				segments(xstart,ystart,x1=xend,y1=yend,col=colors[i])
				count = count + 1
				}
			}	
		}		
	}	












