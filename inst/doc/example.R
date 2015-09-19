### R code from vignette source 'example.Rnw'

###################################################
### code chunk number 1: package
###################################################
options(keep.source = TRUE, width = 60)
packageInfo <- packageDescription("DNAPlotR")
library(DNAPlotR)
packageKeywords<-"visualization, alignment, display, genome, DNA, sequence, multiple sequence alignment, reference sequence"


###################################################
### code chunk number 2: plotDNA (eval = FALSE)
###################################################
##   set.seed(1234)
##   fakeSeqs<-createFakeDNA(5000,500)
##   refSeq<-fakeSeqs[1]
##   fakeSeqs<-fakeSeqs[-1]
##   species<-sprintf('Species %s',sub(' [0-9]+$','',names(fakeSeqs)))
##   par(mar=c(3.5,4.4,.5,7),mgp=c(2.5,1,0))
##   dummy<-plotDNA(fakeSeqs,groups=species)


###################################################
### code chunk number 3: showPlotDNA
###################################################
  set.seed(1234)
  fakeSeqs<-createFakeDNA(5000,500)
  refSeq<-fakeSeqs[1]
  fakeSeqs<-fakeSeqs[-1]
  species<-sprintf('Species %s',sub(' [0-9]+$','',names(fakeSeqs)))
  par(mar=c(3.5,4.4,.5,7),mgp=c(2.5,1,0))
  dummy<-plotDNA(fakeSeqs,groups=species)


