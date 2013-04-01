library(seqinr)
contacts = c("carlia")
type = c("divpoly","dnds","fst")
column = c("net","dnds","Fst")

for (i in 1:length(contacts)) {
	par(mfrow=c(1,3))
	
	for (j in 1:length(type)) {
 
		file <- paste("/Users/singhal/thesisWork/introgression/summaryStatistics/",contacts[i],"_full_", type[j],".out2",sep="")
		d <- read.table(file,header=T)
		
		if (j == 2) {
			d<- d[d$dnds < 2,] 		
		}
		if (j == 1) {
			d<- d[d$net < 0.04,] 		
		}
	
		seqfile <- paste("/Users/singhal/thesisWork/introgression/targetSequences/final/",contacts[i], "_targets.fa.final",sep="")
		seq <- read.fasta(file=seqfile)
		loci <- names(seq)
		loci <- sub("_\\d+","",loci,perl=T)	
		loci <- sub("_\\d+","",loci,perl=T)	
		loci <- sub("_exon\\d+","",loci,perl=T)	
		loci <- unique(loci)
		
		subset <- d[which( d$loci %in% loci),]
				
		a<- density(d[column[j]][,1])
		b<- density(subset[column[j]][,1])
		
		plot(a,col="white",xlab=paste(type[j]),title=NULL)
		polygon(a,col=rgb(0,0,1,1/4),border=F)
		polygon(b,col=rgb(1,0,0,1/4),border=F) 

	
		}
	}
