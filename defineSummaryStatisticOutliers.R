contacts <- c("carlia","gillies","nBMB","sjo")
type <- c("dnds","fst","divpoly")
names <- c("dnds","Fst","net")

for (i in 1:length(contacts)) {
	for (j in 1:length(type)) {
		file <- paste("/Users/singhal/thesisWork/introgression/summaryStatistics/",contacts[i],"_full_",type[j],".out",sep="")
		d <- read.table(file,header=T)	
		
		vec <- d[names[j]][,1]
		vec <- vec[order(vec)]
		
		out05 = round(length(vec) * 0.95)
		out05 = vec[out05]		
		cat(contacts[i],type[j],"outlier1",out05,max(vec),"\n",sep="\t")
		cat(contacts[i],type[j],"notoutlier1",min(vec),out05,"\n",sep="\t")
		
		if (j == 1) {
			cat(contacts[i],type[j],"outlier2",0,1,"\n",sep="\t")
			cat(contacts[i],type[j],"notoutlier2",1,max(vec),"\n",sep="\t")		
		}
		
		split <- c(0.0,0.2,0.4,0.6,0.8,1.0)
		
		for (n in 1:(length(split) - 1)) {
			limit1 <- max(vec) * split[n]
			limit2 <- max(vec) * split[n+1]
			
			name <- paste(n,"value",sep="")
			cat(contacts[i],type[j],name,limit1,limit2,"\n",sep="\t")
												
			}
			
		oldvec <- vec	
		vec <- vec[vec > 0]	
			
		for (n in 1:(length(split) - 1)) {
			val1 <- round(length(vec) * split[n])
			if (n == 1) { limit1 <- oldvec[1] }
			else { limit1 <- vec[val1] }
			
			limit2 <- vec[round(length(vec) * split[n+1])]
			
			name <- paste(n,"loci",sep="")
			cat(contacts[i],type[j],name,limit1,limit2,"\n",sep="\t")
												
			}	
					
		}
	}