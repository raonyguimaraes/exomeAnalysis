d <- read.table("/Users/singhal/Desktop/alignmentSuccess.out",header=F)
boxplot(d$V3~d$V1,ylim=c(0,1),ylab="specificity",xlab="contact")