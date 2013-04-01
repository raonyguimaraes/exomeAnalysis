d <- read.table("/Users/singhal/thesisWork/introgression/clineAndSummaryStats.out",header=T)
dv <- c("transDiff","center","width","trans_net","trans_raw","outlier","outtype")
iv <- c("Fst_snp","contig_Fst","contig_dnds","contig_net","contig_raw","full_Fst","full_dnds","full_net","full_raw","ancDiff")

c <- split(d,d$contact)

pdf("/Users/singhal/Desktop/summaryStats_scatterplots.pdf")
for (i in 1:5) {
	for (j in 1:length(iv)) {
		par(mfrow=c(2,2))			
		for (k in 1:length(c)) {
			tmp <- c[[k]]
			tmp <- tmp[complete.cases(tmp$outtype),]
			plot(tmp[iv[j]][,1], tmp[dv[i]][,1], pch=16, cex=0.5,xlab=paste(iv[j]),ylab=paste(dv[i]))
			
			corr <- round(cor(tmp[iv[j]][,1], tmp[dv[i]][,1],use="complete.obs"),3)
			text(x = max(tmp[iv[j]][,1],na.rm=T) * 0.8 , y = max(tmp[dv[i]][,1],na.rm=T) * 0.8, corr)
			
			title(paste(names(c)[k]))										
		}		
	}
}
dev.off()

d$outtype <- factor(d$outtype, levels = c(levels(d$outtype), "sweep"))
indices <- which(d$type == "widerange")
d$outtype = replace(d$outtype,indices,"sweep")
d <- d[d$ancDiff >= 0.7,] #so no ascertainment bias
c <- split(d,d$contact)

library(multcompView)

pdf("/Users/singhal/Desktop/summaryStats_boxplots.pdf")
for (i in 6:7) {
	for (j in 1:length(iv)) {
		par(mfrow=c(2,2))			
		for (k in 1:length(c)) {
			tmp <- c[[k]]
			tmp <- tmp[complete.cases(tmp$outtype),]
			b <- boxplot(tmp[iv[j]][,1] ~ tmp[dv[i]][,1],ylab=paste(iv[j]),xlab=paste(dv[i]),outline=F)					
			title(paste(names(c)[k]))				
			
			a <- aov(tmp[iv[j]][,1] ~ tmp[dv[i]][,1])
			tukey <- TukeyHSD(a)
			groups <- multcompLetters(extract_p(tukey[[1]]))
			namesT <- b$names
			
			pval<-round(summary(a)[[1]][["Pr(>F)"]][1],3)
			
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
