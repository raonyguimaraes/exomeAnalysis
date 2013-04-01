library(ggplot2)
d<-read.table("/Users/singhal/thesisWork/introgression/clineTypes.txt",header=T)
d$contact <- factor(d$contact,levels=c("nBMB","sjo","carlia","gillies"))
qplot(factor(contact),number,data=d,geom="bar",fill=factor(type),stat="identity")