 library(hzar)
 
 d <- read.table("~/thesisWork/lamcogg/gilliesHZ/formalClines/alleleFrequencies2.txt",header=T)
 data <- read.table("/Users/singhal/thesisWork/lamcogg/gilliesHZ/formalClines/clineFitting.out",header=T)

par(mfrow=c(2,6))
loci <- 11
for (i in 7:17) {

	plot(d$dist,d[,i],pch=16, xlim=c(-100, 2500), col="gray", ylim=c(0,1),main=paste(names(d)[i]), xlab="location (m)", ylab="p")

	tmp <- data.frame(d$dist,d[,i])
	names(tmp) <- c("dist","p")
	cline1 <- nls( p ~ (pmin + ((1+tanh(2*(dist-c)/w))/2) * (pmax - pmin)), data=tmp,start=list(c=1100,w=150,pmin=0,pmax=1) )
	c<-summary(cline1)$coefficients[1]
	w<-summary(cline1)$coefficients[2]
	pmin<-summary(cline1)$coefficients[3]
	pmax<-summary(cline1)$coefficients[4]
	dist <- seq(-100,7000,10)
	fit <- ((1+tanh(2*(dist-c)/w))/2)
	pscaled <- pmin + fit * (pmax - pmin)
	lines(dist,pscaled)
	
	c<- abs(data[data$locus==names(d)[i],]$c - max(d$northing))
	w<-data[data$locus==names(d)[i],]$w
	pmin<-1 - data[data$locus==names(d)[i],]$pmax
	pmax<-1- data[data$locus==names(d)[i],]$pmin
	fit <- ((1+tanh(2*(dist-c)/w))/2)
	pscaled <- pmin + fit * (pmax - pmin)
	lines(dist,pscaled,col="red")
	
	a <- hzar.doMolecularData1DPops(d$dist, d[,i], d$N)
	a_model <- hzar.makeCline1DFreq(a, scaling="free",tails="none"); 
	a_model <- hzar.model.addBoxReq(a_model, min(d$dist),max(d$dist)); 
	a_modelFitR <- hzar.first.fitRequest.old.ML(model= a_model , a, verbose=FALSE); 
	a_modelFitR$mcmcParam$chainLength <- 2e3;
	a_modelFitR$mcmcParam$burnin <- 5e2; 
	a_modelFit <- hzar.doFit(a_modelFitR)
	a_modelData <- hzar.dataGroup.add(a_modelFit);		
	c <- a_modelData$ML.cline$param.all$center
	w <- a_modelData$ML.cline$param.all$width
	pmin <- a_modelData$ML.cline$param.all$pMin
	pmax <- a_modelData$ML.cline$param.all$pMax
	fit <- ((1+tanh(2*(dist-c)/w))/2)
	pscaled <- pmin + fit * (pmax - pmin)
	lines(dist,pscaled,col="blue")
	
	} 
#dev.off()		
