d <- read.table("/Users/singhal/thesisWork/introgression/benchtruth/benchtruthing",header=T)

x <- split(d,d$contact)

contacts <- c("sapro","nbmb","carlia", "gillies")
colors <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")

plot(NULL,xlim=c(0,1),ylim=c(0,1),xlab="actual allele freq.",ylab="estimated allele freq.")

for (i in 1:length(x)) {
	tmp <- x[[contacts[i]]]
	points(tmp$actual,tmp$estimated,pch=16,col=colors[i])
	text(0.8,0.05+0.05*i,paste(round(cor(tmp$actual,tmp$estimated),2),contacts[i],sep=": "),col=colors[i])	
}

text(0.8,0.05,paste(round(cor(d$actual,d$estimated),2),"all",sep=": "))	