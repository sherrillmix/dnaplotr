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


###################################################
### code chunk number 4: plotShortReads (eval = FALSE)
###################################################
##   set.seed(1234)
##   seqLength<-1000
##   fakeSeqs<-createFakeDNA(5000,seqLength,pGap=0)
##   refSeq<-fakeSeqs[1]
##   fakeSeqs<-fakeSeqs[-1]
##   potentialStarts<-1:(seqLength-99)
##   #leave a gap in sequences
##   potentialStarts<-potentialStarts[potentialStarts>550|potentialStarts<400]
##   startCoords<-sort(sample(potentialStarts,5000,TRUE))
##   endCoords<-startCoords+99
##   dummy<-paste(rep('-',seqLength),collapse='')
##   substring(fakeSeqs,1,startCoords)<-substring(dummy,1,startCoords)
##   substring(fakeSeqs,endCoords+1,seqLength)<-substring(dummy,endCoords+1,seqLength)
##   fakeSeqs<-replaceOuterGaps(fakeSeqs)
##   par(mar=c(3.5,4.4,.5,7),mgp=c(2.5,1,0))
##   dummy<-plotDNA(fakeSeqs)


###################################################
### code chunk number 5: showShortReads
###################################################
  set.seed(1234)
  seqLength<-1000
  fakeSeqs<-createFakeDNA(5000,seqLength,pGap=0)
  refSeq<-fakeSeqs[1]
  fakeSeqs<-fakeSeqs[-1]
  potentialStarts<-1:(seqLength-99)
  #leave a gap in sequences
  potentialStarts<-potentialStarts[potentialStarts>550|potentialStarts<400]
  startCoords<-sort(sample(potentialStarts,5000,TRUE))
  endCoords<-startCoords+99
  dummy<-paste(rep('-',seqLength),collapse='')
  substring(fakeSeqs,1,startCoords)<-substring(dummy,1,startCoords)
  substring(fakeSeqs,endCoords+1,seqLength)<-substring(dummy,endCoords+1,seqLength)
  fakeSeqs<-replaceOuterGaps(fakeSeqs)
  par(mar=c(3.5,4.4,.5,7),mgp=c(2.5,1,0))
  dummy<-plotDNA(fakeSeqs)


