files <- list.files(path="/Users/singhal/Desktop/activeWork/prelimClines/",pattern="*.out",full.names=T)

contacts <- c("nBMB","gillies","carlia","sjo")
dist <- read.table("/Users/singhal/thesisWork/introgression/distances/distances",header=T)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')

#only try to fit clines if difference in allele frequency is this much
diff = 0.8

allelefreq <- list()

for (i in 1:length(contacts)) {
	#read in file
	file <- paste("/Users/singhal/Desktop/activeWork/prelimClines/",contacts[i],".cline.out",sep="")
	d <- read.table(file,header=F,stringsAsFactors=F)
	names(d) <- c("locus","pos","pop","af")

	#get rid of control loci
	d<-d[grep('ND4',d$locus),]
	
	if(i==2) {
		d<-d[grep('ND4_Lampro1_1',d$locus),]
	}
	
	#add distance
	dist_contact = dist[dist$contact==contacts[i],]
	d <- data.frame(d,rep(NA,dim(d)[1]))
	names(d)[5] <- c("dist")
	for (p in 1:length(pops)) {
		d[d$pop == pops[p],]$dist = dist_contact[dist_contact$library==pops[p],]$distance
		}
	
	#notekeeping
	tmp_af = data.frame()

	loci <- split(d,d$locus)
	for (a in 1:length(loci)) {
		locus <- loci[[a]]
		pos <- split(locus,locus$pos)
		
		for (b in 1:length(pos)) {
			tmp <- pos[[b]]

			if (is.na(tmp[tmp$pop=='ancN',]$af) || is.na(tmp[tmp$pop=='ancN',]$af) ) {
				if (is.na(tmp[tmp$pop=='10kN',]$af) || is.na(tmp[tmp$pop=='10kS',]$af) ) {
					}
				else {
					diff2 = abs(tmp[tmp$pop=="10kS",]$af - tmp[tmp$pop=="10kN",]$af)
					tmp <- tmp[complete.cases(tmp$af),]
					if (dim(tmp)[1] == 11) {
						if (diff2 >= diff & tmp[tmp$pop=="10kN",]$af < 0.05 & tmp[tmp$pop=="10kS",]$af > 0.95 & tmp[tmp$pop=="ancN",]$af < 0.05 & tmp[tmp$pop=="ancS",]$af > 0.95) {
							if (i == 4) {
								if (tmp[tmp$pop=="2kN",]$af < 0.5) {
									tmp_af = rbind(tmp_af,tmp)	
									}
								}
							else {
								tmp_af = rbind(tmp_af,tmp)
								}
							}	
						}	
					}
				}	
			else {
				#this is where i want to be!
				diff1 = abs(tmp[tmp$pop=="ancN",]$af - tmp[tmp$pop=="ancS",]$af)
				if (is.na(tmp[tmp$pop=='10kN',]$af) || is.na(tmp[tmp$pop=='10kS',]$af) ) {
					diff2 = 0
					}
				else {	
					diff2 = abs(tmp[tmp$pop=="10kS",]$af - tmp[tmp$pop=="10kN",]$af)
					}
				tmp <- tmp[complete.cases(tmp$af),]
				if (dim(tmp)[1] == 11) {
					if (diff2 >= diff || diff1 >= diff) {						
						if (tmp[tmp$pop=="10kN",]$af < 0.05 & tmp[tmp$pop=="10kS",]$af > 0.95 & tmp[tmp$pop=="ancN",]$af < 0.05 & tmp[tmp$pop=="ancS",]$af > 0.95) {
							if (i == 4) {
								if (tmp[tmp$pop=="2kN",]$af < 0.5) {
									tmp_af = rbind(tmp_af,tmp)	
									}
								}
							else {
								tmp_af = rbind(tmp_af,tmp)
								}
							}	
						}		
					}	
				}	
			}	
		}
	allelefreq[[i]] = tmp_af	
	}	

par(mfrow=c(2,2),mai=c(0.7,0.7,0.3,0.3),oma=c(0,0,0,0))
for (i in 1:length(allelefreq)) {
	dist_contact = dist[dist$contact==contacts[i],]
	dist_contact = dist_contact[order(dist_contact$distance),]
	boxplot(allelefreq[[i]]$af~allelefreq[[i]]$dist,names=seq(1,11,by=1),outline=F,col="gray",xlab="population along transect",ylab="allele frequency")
	text(1,0.9,paste(contacts[i]))	
}