d <- read.table("/Users/singhal/thesisWork/introgression/metrics/annotatedOrphans.out",header=F)
pie(d$V2,labels=d$V1)