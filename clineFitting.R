library(nls2)
library(reshape)

contacts <- c("nBMB","gillies","carlia","sjo")
dist <- read.table("/Users/singhal/thesisWork/introgression/distances/distances",header=T)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')

#only try to fit clines if difference in allele frequency is this much
diff = 0.7
diffwiderange = 0.3

adjustments = c(0,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)
fitting <- list()

paramfile <- paste("/Users/singhal/thesisWork/introgression/clineAF/starting",sep="")
param <- read.table(paramfile,header=T,stringsAsFactors=F)

for (i in 1:length(contacts)) {
	#read in file
	file <- paste("/Users/singhal/thesisWork/introgression/clineAF/",contacts[i],".cline.out",sep="")
	d <- read.table(file,header=F,stringsAsFactors=F)
	names(d) <- c("locus","pos","pop","af")

	#get rid of control loci
	d<-d[grep('ND4',d$locus,invert=T),]
	d<-d[grep('16s',d$locus,invert=T),]
	d<-d[grep('utr',d$locus,invert=T),]
	
	#add distance
	dist_contact = dist[dist$contact==contacts[i],]
	d <- data.frame(d,rep(NA,dim(d)[1]))
	names(d)[5] <- c("dist")
	for (p in 1:length(pops)) {
		d[d$pop == pops[p],]$dist = dist_contact[dist_contact$library==pops[p],]$distance
		}
	
	#notekeeping
	center <- c()
	width <- c()
	name <- c()
	position <- c()
	ancDiff <- c()
	transDiff <- c()
	type <- c()
	k <- 1
		
	c_param = param[param$contact==contacts[i] & param$measure == "c",]$value
	w_param = param[param$contact==contacts[i] & param$measure == "w",]$value

	loci <- split(d,d$locus)
	for (a in 1:length(loci)) {
		locus <- loci[[a]]
		pos <- split(locus,locus$pos)
		
		for (b in 1:length(pos)) {
			tmp <- pos[[b]]			
			c_new = c_param			
			w_new = w_param
			fit <- NA
			try <- NA
			
			if (is.na(tmp[tmp$pop=='ancN',]$af) || is.na(tmp[tmp$pop=='ancN',]$af) ) {
				diff1 = NA;
				}
			else {
				diff1 = abs(tmp[tmp$pop=="ancN",]$af - tmp[tmp$pop=="ancS",]$af)	
				}	
					
			tmp <- tmp[tmp$pop!="ancN",]
			tmp <- tmp[tmp$pop!="ancS",]		
			min <- min(tmp$af,na.rm=T)
			max <- max(tmp$af,na.rm=T)
			diff2 = max - min
						
			if (is.na(diff1)) { 
				if (!is.na(tmp[tmp$pop=='10kN',]$af)) {
					if (tmp[tmp$pop=='10kN',]$af >= 0.5) {
						tmp$af = 1 - tmp$af
						}
					}						
				}
			else if (diff1 <= diff) {
				if (!is.na(tmp[tmp$pop=='10kN',]$af)) {
					if (tmp[tmp$pop=='10kN',]$af >= 0.5) {
						tmp$af = 1 - tmp$af
						}
					}						
				}		
															
			tmp <- tmp[tmp$pop!="10kN",]
			tmp <- tmp[tmp$pop!="10kS",]			
			tmp <- tmp[complete.cases(tmp$af),]
			
			if (dim(tmp)[1] == 7) {
				if (diff2 >= diff) {															
					for (m in 1:length(adjustments)) {
						for (n in 1:length(adjustments)) {
							c_new = c_param + c_param * adjustments[m]
							w_new = w_param + w_param * adjustments[n]									
							if(is.na(fit)) {
								fit = tryCatch(nls( af ~ (((1+tanh(2*(dist-c)/w))/2)), data=tmp,start=list(c=c_new,w=w_new)), error=function(e) NA);
								}
							}								
						}
							
					if (is.na(fit)) {
						fit = tryCatch(nls2( af ~ (((1+tanh(2*(dist-c)/w))/2)), data=tmp,start=list(c=c_new,w=w_new)), error=function(e) NA);	
						}
								
					if (is.na(fit)) {
						type[k] <- NA
						center[k] <- NA								
						width[k] <- NA
						ancDiff[k] <- diff1
						transDiff[k] <- diff2
						name[k] <- names(loci)[a]
						position[k] <- names(pos)[b]
						k = k+1
						}
					else {
						type[k] <- "cline"
						center[k] <- summary(fit)$coefficients[1]
						width[k] <- summary(fit)$coefficients[2]
						ancDiff[k] <- diff1
						transDiff[k] <- diff2
						name[k] <- names(loci)[a]
						position[k] <- names(pos)[b]
						k = k+1
						}	
					}
				else {
					if(!is.na(diff1)) {
						if (diff1 >= diff & diff2 <= diffwiderange) {
							type[k] <- "widerange"
 							center[k] <- NA								
							width[k] <- NA
							ancDiff[k] <- diff1
							transDiff[k] <- diff2
							name[k] <- names(loci)[a]
							position[k] <- names(pos)[b]
							k = k+1	
							}
						else {
							type[k] <- NA
 							center[k] <- NA								
							width[k] <- NA
							ancDiff[k] <- diff1
							transDiff[k] <- diff2
							name[k] <- names(loci)[a]
							position[k] <- names(pos)[b]
							k = k+1
							}							
						}
					}		
				}
			}
		}		
	fitting[[i]] <- data.frame(name,position,center,width,ancDiff,transDiff,type)
	names(fitting)[i] <- contacts[i]
	}	

fit <- data.frame()
for (i in 1:length(fitting)) {
	tmp <- data.frame(fitting[[i]],rep(names(fitting)[i],dim(fitting[[i]])[1]))
	names(tmp)[dim(tmp)[2]] <- "contact"
	fit <- rbind(fit,tmp)
}

fit <- read.table("/Users/singhal/thesisWork/introgression/clineData",header=T)
fit <- fit[complete.cases(fit$type),]
fit <- fit[fit$type=="cline",]

widthfit <- fit
centerfit <- fit

centerfit <- centerfit[complete.cases(centerfit$center),]
centerfit <- centerfit[centerfit$center > 0,]
centerfit <- centerfit[centerfit$center < 30000,]
centerfit <- data.frame(centerfit, rep(NA, dim(centerfit)[1]))
names(centerfit)[length(centerfit)] <- c("center_standarized")
for (i in 1:length(contacts)) {
	centermedian <- median(centerfit[centerfit$contact==contacts[i],]$center)
	centerfit[centerfit$contact==contacts[i],]$center_standarized = centerfit[centerfit$contact==contacts[i],]$center	- centermedian	
}

widthfit$contact <- factor(widthfit$contact,levels=c("sjo","nBMB","carlia","gillies"))
a<-ggplot(widthfit, aes(width, colour = contact,fill= contact)) + geom_density(aes(y = ..scaled..),size=0.75)+scale_x_continuous(limits = c(0, 100000))+xlab("width (m)")+ylab("scaled density")
centerfit$contact <- factor(centerfit$contact,levels=c("sjo","nBMB","carlia","gillies"))
b<-ggplot(centerfit, aes(center_standarized, colour = contact,fill=contact)) + geom_density(aes(y = ..scaled..),size=0.75)+scale_x_continuous(limits = c(-15e3, 15e3))+xlab("center (m)")+ylab("scaled density")
multiplot(a,b,cols=1)				
				
				