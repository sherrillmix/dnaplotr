set.seed(12345)

png('dnaPlotExample.png',width=1600,height=800)
	fakeSeqs<-createFakeDNA(1000)
	refSeq<-fakeSeqs[1]
	fakeSeqs<-fakeSeqs[-1]
	species<-sprintf('Species %s',sub(' [0-9]+$','',names(fakeSeqs)))
	par(mar=c(3.7,4.6,.5,9.5),cex.axis=1.8,cex.lab=1.8,mgp=c(2.5,.9,0))
	plotDNA(fakeSeqs,groups=species,groupCexScale=TRUE)
dev.off()
