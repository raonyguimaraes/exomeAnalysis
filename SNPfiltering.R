contacts <- c("nBMB","gillies","carlia","sjo")
dist <- read.table("/Users/singhal/thesisWork/introgression/distances/distances",header=T)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')

#only try to fit clines if difference in allele frequency is this much
diff = 0.8

for (i in 1:length(contacts)) {
	#read in file
	file <- paste("/Users/singhal/thesisWork/introgression/clineAF/",contacts[i],".cline.out",sep="")
	d <- read.table(file,header=F,stringsAsFactors=F)
	names(d) <- c("locus","pos","pop","af")
		
	k <- 1

	loci <- split(d,d$locus)
	for (a in 1:length(loci)) {
		locus <- loci[[a]]
		pos <- split(locus,locus$pos)
		
		for (b in 1:length(pos)) {
			tmp <- pos[[b]]			

			min <- min(tmp$af,na.rm=T)
			max <- max(tmp$af,na.rm=T)
			diff2 = max - min
					
			tmp <- tmp[tmp$pop!="ancN",]
			tmp <- tmp[tmp$pop!="ancS",]
			tmp <- tmp[tmp$pop!="10kN",]
			tmp <- tmp[tmp$pop!="10kS",]			
			tmp <- tmp[complete.cases(tmp$af),]	
									
			if (diff2 >= diff) {
				try = 1	
				}			
											
			if (dim(tmp)[1] == 7) {
				if (!is.na(try)) {	
					k <- k + 1												
					}	
				}
			}
		}		
	cat(contacts[i],k,"\n",sep="\t")
	}	
			
				
				