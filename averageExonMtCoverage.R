contacts = c("nBMB","sjo","carlia","gillies")
allcov <- data.frame()

for (n in 1:length(contacts)) {

	file <- paste("/Users/singhal/thesisWork/introgression/coverage/",contacts[n],".depth.summarized",sep="")
	d <- read.table(file,header=F)

	a <- d[grep("ND4",d$V1),]
	print(a)

	names(d) <- c("locus","10kN","2kN","1kN","nTail","center","sTail","1kS","2kS","10kS")

	pops <- c("10kN","2kN","1kN","nTail","center","sTail","1kS","2kS","10kS")
	cov <- data.frame()
	for (i in 1:length(pops)) {
		tmp <- data.frame(d[,i+1],rep(pops[i],length(d[,i+1])))
		cov <- rbind(cov,tmp)		
	}
	names(cov) <- c("coverage","pop")

	b <- aggregate(cov$coverage,by=list(cov$pop),mean)
	print(b)
	
	}	
	