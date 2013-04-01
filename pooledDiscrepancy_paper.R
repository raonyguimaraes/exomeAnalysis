nsim <- 100000
af_c <- c(0.01,0.1,0.25,0.5)
sample_c <- c(40) # number of chromosomes
coverage_c <- c(50,100,200,500)
pdf("/Users/singhal/Desktop/af.pdf")

par(mfrow=c(length(coverage_c),length(af_c)),oma=c(2,2,0,0),cex=0.5,mar=c(4,4,2,2))

for (c in 1:length(coverage_c)) {
	for (a in 1:length(af_c)) {
		af <- af_c[a]
		sample <- sample_c[1]
		coverage <- coverage_c[c]
	
		af_pop_true <- rbinom(nsim,sample,af) / sample

		af_pop_pool <- rep(NA,nsim)
		for (i in 1:length(af_pop_true)) {
			af_pop_pool[i] = rbinom(1,coverage,af_pop_true[i])/coverage
			}
		
		hist(af_pop_pool,freq=T,col="gray",border="gray",xlab="allele freq.",main=NULL,xlim=c(0,1))
		abline(v=af,col="red",lty=2)
		title(main=paste("coverage:",coverage,sep=""))
		}	
	}	
dev.off()