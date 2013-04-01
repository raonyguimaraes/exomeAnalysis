d <- read.table("/Users/singhal/thesisWork/introgression/clineData",header=T)
d <- data.frame(d,rep(NA,dim(d)[1]))
names(d)[length(d)] = "outlier"

contact <- d$contact
c <- split(d,d$contact)

par(mfrow=c(2,2))
for (i in 1:length(c)) {
	tmp <- c[[i]]
	tmp <- tmp[complete.cases(tmp$type),]
	cline <- tmp[tmp$type=="cline",]
	wide <- tmp[tmp$type=="widerange",]
	cat(names(c)[i],dim(cline)[1],dim(wide)[1],"\n",sep="\t")	
	
	width <- cline$width
	width <- width[width >= 0]
	width <- sort(width)		
					
	xmax <- width[round(length(width)*0.95,digits=0)]		
	width <- width[width <= xmax]
	xmin <- width[round(length(width)*0.05,digits=0)]
	xmax <- width[round(length(width)*0.95,digits=0)]	
	
	numBreaks = length(center)/50
	a <- hist(width,breaks=numBreaks,col="gray",border=F,main=paste(names(c)[i]),xlab="width (m)")	
	a
	
	small <- width[width <= xmin]
	hist(small,breaks=a$breaks,col="red",add=T,border=F)
	
	big <- width[width >= xmax]
	hist(big,breaks=a$breaks,col="red",add=T,border=F)
					
	}