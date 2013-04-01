d <- read.table("/Users/singhal/Desktop/activeWork/prelimClines/gillies.cline.out",header=F,stringsAsFactors=F)
names(d) <- c("locus","pos","pop","af")

#sapro
#dist <- c(4494,7781,0,8893,679,9414,12707)
#carlia
#dist <- c(764,1278,2334,2676,3175,3689,6732)
#gillies
#dist <- c(0,3730,4764,4803,4829,5319,11295)
#nbmb
dist <- c(3497,4988,4075,5546,6029,6419,6740)
pops <- c('10kN','10kS','1kN','1kS','2kN','2kS','ancN','ancS','center','nTail','sTail')
order <- c(2,10,5,8,3,9,1,11,6,4,7)

num <- c(16,16,16,16,16,16,16)

ND4 <- d[grep('ND4',d$locus),]
ND4 <- data.frame(ND4,rep(NA,dim(ND4)[1]))
names(ND4)[5] <- c("order")

for (i in 1:length(pops)) {
	 ND4[ND4$pop == pops[i],]$order = order[i]
}

boxplot(ND4$af ~ ND4$order)



#######################################
# number of loci for fitting analysis #
#######################################

d<-d[grep('ND4',d$locus,invert=T),]
d<-d[grep('16s',d$locus,invert=T),]
d<-d[grep('utr',d$locus,invert=T),]

cline = 1
bad = 1
longrange = 1
incomplete = 1
toosmall = 1
diff = 0.8

locus <- split(d,d$locus)
for (i in 1:length(locus)) {
	loci <- locus[[i]]
	pos <- split(loci,loci$pos)
	
	for (j in 1:length(pos)) {
		tmp <- pos[[j]]
		
		#have no ancestral info
		if (is.na(tmp[tmp$pop=='ancN',]$af) || is.na(tmp[tmp$pop=='ancN',]$af) ) {
			if (is.na(tmp[tmp$pop=='10kN',]$af) || is.na(tmp[tmp$pop=='10kS',]$af) ) {
				#hmm, this shouldn't happen
				bad = bad + 1
				}
			else {
				diff2 = abs(tmp[tmp$pop=="10kS",]$af - tmp[tmp$pop=="10kN",]$af)
				tmp <- tmp[tmp$pop!="10kN",]
				tmp <- tmp[tmp$pop!="10kS",]
				tmp <- tmp[tmp$pop!="ancN",]
				tmp <- tmp[tmp$pop!="ancS",]
				tmp <- tmp[complete.cases(tmp$af),]
				if (dim(tmp)[1] == 7) {
					if (diff2 >= diff) {
						cline = cline + 1
					}			
					else {
						toosmall = toosmall + 1
					}
				}	
				else {
					incomplete = incomplete + 1
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
			tmp <- tmp[tmp$pop!="10kN",]
			tmp <- tmp[tmp$pop!="10kS",]
			tmp <- tmp[tmp$pop!="ancN",]
			tmp <- tmp[tmp$pop!="ancS",]
			tmp <- tmp[complete.cases(tmp$af),]
			longrange = longrange + 1
			if (dim(tmp)[1] == 7) {
				if (diff2 >= diff || diff1 >= diff) {
					cline = cline + 1
					}			
				else {
					toosmall = toosmall + 1
					}
				}	
			else {
				incomplete = incomplete + 1
			}
		}	
	}
}		

cat("cline",cline,"\n",sep=" ")
cat("bad",bad,"\n",sep=" ")
cat("longrange",longrange,"\n",sep=" ")
cat("incomplete",incomplete,"\n",sep=" ")
cat("toosmall",toosmall,"\n",sep=" ")
