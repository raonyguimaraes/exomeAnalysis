d <- read.table("/Users/singhal/thesisWork/introgression/clineAndSummaryStats.out",header=T)
dv <- c("transDiff","center","width","trans_net","trans_raw","outlier","outtype")

c <- split(d,d$contact)

contacts <- c("nBMB","carlia","sjo", "gillies")
colors <- c("#7CAE00","#00BFC4","#F8766D","#C77CFF")

graph <- c("Fst_snp","contig_Fst","contig_net","contig_dnds")
par(mfrow=c(2,2),mar=c(4,4,2,1),oma=c(1,1,1,1))	
for (i in 3:3) {
	for (j in 1:length(graph)) {
		tmp <- d		
		tmp <- tmp[complete.cases(tmp$outtype),]	
		plot(NULL,xlim=range(tmp[graph[j]][,1],na.rm=T),ylim=c(0,10000),xlab=paste(graph[j]),ylab=paste(dv[i]))
		for (k in 1:length(c)) {
			tmp <- c[[contacts[k]]]
			points(tmp[graph[j]][,1], tmp[dv[i]][,1], pch=16, cex=0.5, col=colors[k]) 		}							
		}
	}

d$outtype <- factor(d$outtype, levels = c(levels(d$outtype), "sweep"))
indices <- which(d$type == "widerange")
d$outtype = replace(d$outtype,indices,"sweep")
d <- d[complete.cases(d$ancDiff),] #so no ascertainment bias
d <- d[d$ancDiff >= 0.7,] #so no ascertainment bias
c <- split(d,d$contact)

library(multcompView)

pdf("/Users/singhal/Desktop/formalSummaryStatsBoxplots.pdf")
graph <- c("Fst_snp","contig_Fst","contig_net","contig_dnds")
par(mfrow=c(2,2),mar=c(3,4,3,3),oma=c(1,1,1,1))	
for (i in 7:7) {
	for (j in 1:length(graph)) {		
	#	for (k in 1:length(c)) {
			tmp <- d
			tmp <- tmp[complete.cases(tmp$outtype),]
			b <- boxplot(tmp[graph[j]][,1] ~ tmp[dv[i]][,1],ylab=paste(graph[j]),xlab=paste(dv[i]),outline=F)					
			#title(paste(names(c)[k]))				
			
			a <- aov(tmp[graph[j]][,1] ~ tmp[dv[i]][,1])
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
	#}
}
dev.off()
