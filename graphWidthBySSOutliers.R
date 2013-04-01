library(ggplot2)

d <- read.table("/Users/singhal/thesisWork/introgression/clineAndSummaryStats.out",header=T)

d$outtype <- factor(d$outtype, levels = c(levels(d$outtype), "sweep"))
indices <- which(d$type == "widerange")
d$outtype = replace(d$outtype,indices,"sweep")

d<-d[complete.cases(d$outtype),]

a <- read.table("/Users/singhal/thesisWork/introgression/summaryStatistics/outlierBounds.out",header=F,stringsAsFactors=F)

contacts <- c("nBMB","sjo","carlia","gillies")
types <- c("dnds","fst","divpoly")
name <- c("full_dnds","full_Fst","full_net")

pdf("/Users/singhal/Desktop/summaryStats_boxplots.pdf")
for (j in 1:length(types)) {
	ss <- a[a$V2==types[j],]
	bounds <- unique(sub("^\\d","",ss$V3,perl=T))
	for (k in 1:length(bounds)) {	
		par(mfrow=c(2,2))	
		for (i in 1:length(contacts)) {
			tmp <- d
			tmp <- tmp[tmp$contact==contacts[i],]
			tmp <- tmp[complete.cases(tmp[name[1]][,1]),]	
			tmp <- data.frame(tmp,rep(NA,dim(tmp)[1]))
			names(tmp)[length(tmp)] <- "group"
			
			sstmp <- ss[ss$V1==contacts[i],]
			sstmp <- sstmp[grep(bounds[k],sstmp$V3),]

			for (n in 1:dim(sstmp)[1]) {
				upper <- sstmp[n,]$V5
				lower <- sstmp[n,]$V4
				id <- sstmp[n,]$V3
					
				indices <- which(tmp[name[j]][,1] >= lower & tmp[name[j]][,1] < upper)
				tmp$group = replace(tmp$group,indices,id)				
				}
		
			b <- boxplot(tmp$width ~ tmp$group,ylab="width",xlab=paste(types[j]),outline=F)			
			title(paste(contacts[i]))				
			
			aov <- aov(tmp$width ~ as.factor(tmp$group))
			tukey <- TukeyHSD(aov)
			groups <- multcompLetters(extract_p(tukey[[1]]))
			namesT <- b$names
			
			pval<-round(summary(aov)[[1]][["Pr(>F)"]][1],3)
			
			if (pval < 0.05) {
				letters <- rep(NA,length(namesT))
				for (n in 1:length(namesT)) {
					index <- which( names ( groups[[1]] ) == namesT[n]) 
					letters[n] <- groups[[1]][index]  
				}
				text(x=1:length(namesT), y = max(b$stats) * 0.8, letters )
			}	
		
		
		
			}
		}		
	}
dev.off()	