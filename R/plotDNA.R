#' Functions to plot out a collection of DNA sequences
#'
#' Produce a visual representation of a vector of DNA sequences
#'
#' The main function is:
#'      \describe{
#'        \item{\code{\link{plotDNA}}:}{to produce a plot from a vector of DNA sequences}
#'      }
#'
#' And main helper functions are:
#'      \describe{
#'        \item{\code{\link{replaceOuterGaps}}:}{to mark gaps at the ends of sequences differently than internal indels}
#'      }
#'
#' A vignette is a available using \code{vignette('example',package='dnaplotr')}.
#'
#'
#'
#' @docType package
#' @name dnaplotr
#' @author Scott Sherrill-Mix, \email{shescott@@upenn.edu}
#' @examples
#' seqs<-c('ACTGGA','ACTGCA','ACTGGC')
#' plotDNA(seqs)
NULL

#' Find contiguous ranges within a vector of indices
#' 
#' Take a vector of indices and returns ranges of contiguous regions.
#'
#' @param index Indices to be condensed e.g. from \code{\link{which}} or \code{\link{grep}}
#'
#' @return A data frame with columns: 
#'      \describe{
#'        \item{start:}{Start of a contiguous region}
#'        \item{end:}{End of a contiguous region}
#'      }
#'
#' @export
#' 
#' @examples
#' indexToRange(c(1:10,11,14,16,17:20))
indexToRange<-function(index){
  index<-sort(unique(index))
  diffs<-c(diff(index),Inf)
  ends<-c(which(diffs>1))
  starts<-c(1,ends[-length(ends)]+1)
  return(data.frame('start'=index[starts],'end'=index[ends]))
}

#' Plot a bunch of DNA sequences
#' 
#' Take a vector of strings representing DNA sequences and plot them to the
#' current device. A, C, T and G are colored, - are colored gray and all other
#' characters are white.
#'
#' @param seqs A character vector containing DNA sequences
#' @param seqCounts A integer vector with the number of counts for each sequence. This can be used to improve run time and file size if some sequences are duplicated multiple times (default: 1 for each entry in seqs)
#' @param cols A named vector with names corresponding to the DNA bases and values showing the appropriate color (default: A: green, T: red, C: blue, G: yellow, -: grey)
#' @param xlab A string specifying the x-axis label (default: Position)
#' @param ylab A string specifying the y-axis label (default: Sequence read)
#' @param display A logical vector with element names in 'legend', 'xAxis', 'yAxis', 'groups' where a FALSE suppresses outputting that plot element (default: TRUE for all or any missing elements)
#' @param xStart First base in plot should be labelled as this (default: 1)
#' @param groups Group sequences by group and show label on right side of plot. Note that any prior ordering of sequences will be disrupted. Use a factor and reorder the levels to set a particular order of groups.
#' @param groupCexScale A logical whether to scale group label size by the number of sequences. Useful to highlight more abundant groups and help squeeze in labels on smaller groups.
#' @param refSeq Reference sequence used for numbering the x-axis without counting gaps present in this sequence (note that for further annotations outside this function, e.g. abline(v=3), the axis will be from xStart:xStart+max(nchar(seqs)) without any adjustments to ignore reference gaps
#' @param ... Additional arguments to plot
#'
#' @return NULL
#'
#' @export
#' 
#' @examples
#' plotDNA(c('ACACA','ACACA','ACACT'))
#' refSeq<-'AC---A'
#' seqs<-c('ACTGGA','ACTGCA','ACTGGC','GCTGGG','GGGG',refSeq)
#' par(mar=c(4,5,.5,6),cex.axis=2,cex.lab=2)
#' groupOrder<-c('Group2','Group1','Group3','Reference')
#' groups<-factor(c('Group1','Group2','Group3','Group1','Group3','Reference'),levels=groupOrder)
#' seqCounts<-c(30,10,10,15,5,5)
#' plotDNA(seqs,seqCounts=seqCounts,groups=groups,xStart=10,groupCexScale=TRUE,refSeq=refSeq)
#' fakeSeqs<-createFakeDNA(1000)
#' refSeq<-fakeSeqs[1]
#' fakeSeqs<-fakeSeqs[-1]
#' species<-sprintf('Species %s',sub(' [0-9]+$','',names(fakeSeqs)))
#' par(mar=c(3.5,4.4,.5,7))
#' plotDNA(fakeSeqs,groups=species,groupCexScale=TRUE)
#' fakeAA<-c('MALWTRLRPLLALLALWPPPPARAFVNQHLCGSHLVEALY',
#' 'MALWTRLRPLLALLALWPLPPARAFVNQHLCGSHLVEALY',
#' 'MALWTRLRPLLALLALWPPPPARAFVNX')
#' plotAA(fakeAA,groups=c('Ref','Sub','Stop'))
#things to add back:
# gapTrim
# orderBy?
# refSeq matched - into white
# refseq display
# distOrder
plotDNA<-function(seqs,seqCounts=rep(1,length(seqs)),cols=c('A'='green','T'='red','C'='blue','G'='yellow','-'='grey','default'='white'),xlab='Position',ylab='Sequence Read',display=c('groups'=!is.null(groups)),xStart=1,groups=NULL,groupCexScale=FALSE, refSeq=NULL,...){
  if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
  if(length(seqs)!=length(seqCounts))stop(simpleError('Lengths of seqs and seqCounts not equal'))
  if(!is.null(groups)&&length(seqs)!=length(groups))stop(simpleError('Lengths of seqs and groups not equal'))
  if(!is.null(groups)&&any(is.na(groups)))stop(simpleError('NAs in sequence groups'))
  seqs<-toupper(seqs)
  displayOptions<-c('legend','xAxis','yAxis','groups')
  missingOptions<-!displayOptions %in% names(display)
  if(any(missingOptions))display[displayOptions[missingOptions]]<-TRUE

  if(!is.null(groups)){
    groups<-as.factor(groups)
    groupOrder<-order(as.numeric(groups),1:length(seqs)) #preserve original order within groups
    seqs<-seqs[groupOrder]
    seqCounts<-seqCounts[groupOrder]
    groups<-groups[groupOrder]
  }

  #fill any trailing gaps
  seqMat<-seqSplit(seqs,fill='-')
  seqNum<-seqMat
  if(!'default' %in% names(cols))cols['default']<-'white'
  seqNum[,]<-cols['default']
  for(ii in names(cols))seqNum[seqMat==ii]<-cols[ii]

  #Converting to first base as 0 for ease of use
  xStart<-xStart-1

  graphics::plot(1,1,xlim=xStart+c(0.5,ncol(seqNum)+.5),ylim=c(0.5,sum(seqCounts)+.5),ylab="",xlab=xlab,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',...)
  
  #y axis
  prettyY<-pretty(1:min(sum(seqCounts)))
  prettyY<-prettyY[round(prettyY)==prettyY]
  if(display['yAxis'])graphics::axis(2,prettyY,format(prettyY,scientific=FALSE,big.mark=','),mgp=c(3,.6,0),las=1)
  graphics::title(ylab=ylab,line=3.25,las=3)

  if(!is.null(refSeq)){
    maxNoGap<-gapToNoGap(refSeq,ncol(seqNum))
    prettyX<-pretty(xStart+c(1,maxNoGap))
    prettyX<-prettyX[prettyX<=xStart+maxNoGap]
    prettyX[prettyX==0]<-1
    prettyX<-unique(prettyX[prettyX>xStart])
    prettyXPos<-noGapToGap(refSeq,prettyX-xStart)+xStart
  }else{
    prettyX<-pretty(xStart+c(1,ncol(seqNum)))
    prettyXPos<-prettyX
  }
  if(display['xAxis'])graphics::axis(1,prettyXPos,prettyX)
  #needs to be slight overlap to avoid stupid white line problem
  spacer<-.001
  for(ii in 1:ncol(seqNum)){
    #1 rectangle per read
    #rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
    bottoms<-c(0,cumsum(seqCounts)[-length(seqCounts)])
    tops<-cumsum(seqCounts)
    thisCols<-seqNum[,ii]
    #1 rectangle per repped read 
    #rect(i-.5,bottoms+.5,i+.5,tops+.5+spacer,col=cols,border=NA)
    uniqCols<-unique(thisCols)
    colRanges<-do.call(rbind,lapply(unique(thisCols),function(x){
      out<-indexToRange(which(thisCols==x))
      out$bottom<-bottoms[out$start]
      out$top<-tops[out$end]
      out$col<-x
      return(out)
    }))
    #1 rectangle per string of identical bases
    graphics::rect(xStart+ii-.5,colRanges$bottom+.5,xStart+ii+.5,colRanges$top+.5+spacer,col=colRanges$col,border=NA)
  }
  if(!is.null(groups)){
    groupOrder<-rep(groups,seqCounts)
    maxGroupCount<-max(table(groupOrder))
    uniqueGroups<-unique(groups)
    for(ii in uniqueGroups){
      thisMin<-min(which(groupOrder==ii))
      thisMax<-max(which(groupOrder==ii))
      if(groupCexScale)cexScale<-((diff(c(thisMin,thisMax))+1)/maxGroupCount)^.5
      else cexScale<-1
      if(display['groups'])graphics::mtext(sub('^[$^]','',ii),4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*graphics::par('cex.axis')),line=.5,las=2)
      graphics::abline(h=c(thisMin-.5,thisMax+.5))
    }
  }
  graphics::box()
  if(display['legend']){
    insetPos<-c(graphics::grconvertX(1,'nfc','user'),graphics::grconvertY(0,'nfc','user')) #-.01 could cause trouble here
    legendCols<-cols[!names(cols) %in% c('default','-')]
    graphics::legend(insetPos[1],insetPos[2], names(legendCols),col=legendCols, pt.bg=legendCols,pch = 22,ncol=max(4,length(legendCols)/2),bty='n',xjust=1,yjust=0,xpd=NA,cex=graphics::par('cex.axis'),x.intersp=0.75)
  }
  invisible(NULL)
}

#' @describeIn plotDNA Plot a bunch of AA sequences
#' @param mar margin sizes as in \code{\link{par}}. If left at default then use the current margin settings increasing margin[1] (bottom) to 6.5 if less than 6.5 (needed to give the amino acid legend more space by default). Set explicitly if this behavior is undesired.
#' @export
plotAA<-function(...,mar=NULL,cols=c(dnaplotr::aminoCols,'-'='grey')){
  if(is.null(mar)){
    mar<-graphics::par('mar')
    if(mar[1]<6.5)mar[1]<-6.5
    graphics::par(mar=mar)
  }
  plotDNA(...,cols=cols)
}

#' Convenience function for binding a bunch of sequences together
#' 
#' Take a vector of strings and return a matrix with each row corresponding to
#' a string and each column a position in those strings. If fill is set then
#' strings are padded to an equal length. Otherwise errors out if strings are
#' unequal length.
#'
#' @param ... Character vectors to turn into a matrix
#' @param fill A character to fill the ends of short strings with. If NULL then all strings must be the same length.
#'
#' @return A matrix with a single character in each cell with rows corresponding to the number of strings and columns corresponding to the position in the string
#'
#' @export
#'
#' @examples
#' seqSplit('AAAT','AA','ATT',fill='-')
seqSplit<-function(...,fill=NULL){
  seqs<-c(...)
  if(any(is.na(seqs)))stop(simpleError('NA sequence found in seqSplit'))
  seqN<-nchar(seqs)
  maxN<-max(seqN)
  dummy<-paste(rep(fill,maxN-min(seqN)+100),collapse='')
  if(is.null(fill)&&any(seqN!=maxN))stop(simpleError('All sequences not same length'))
  else seqs<-sprintf('%s%s',seqs,substring(dummy,1,maxN-seqN))
  return(do.call(rbind,strsplit(seqs,'')))
}

#' Convenience function to convert gapped coordinates to what the coordinates would be without gaps
#'
#' Take a reference sequence with gaps and coordinates in that gapped reference
#' sequence and translate the coordinates to the corresponding positions in the
#' gap free sequence. Note that positions landing on a gap are converted to the
#' first 5' nongap position e.g. the fifth position of 'G-----G' is converted
#' to 1 and the second position of '--G--GG' is converted to 0
#'
#' @param refSeq the sequence containing gaps
#' @param coords coordinates on the gapped refSeq to be converted into equivalent nongap coordinates
#' @param gapChars characters interpreted as gaps
#'
#' @return A vector of coordinates in the gap free reference sequence
#'
#' @export
#'
#' @examples
#' gapToNoGap('AA--AA-A',c(1:8))
gapToNoGap<-function(refSeq,coords,gapChars=c('*','.','-')){
  gapSeqSplit<-strsplit(refSeq,'')[[1]]
  nonDash<-!gapSeqSplit %in% gapChars
  newCoords<-cumsum(nonDash)
  coords[coords<1|coords>length(newCoords)]<-NA
  return(newCoords[c(0,coords)]) #using 0 to prevent x[NA] returning everything
}

#' Convenience function to convert ungapped coordinates in a reference function to what the coordinates would be in the gapped reference sequence
#'
#' Take a reference sequence with gaps and coordinates in the ungapped
#' reference sequence and translate the coordinates to the corresponding
#' positions in the gapped sequence
#'
#' @param refSeq the sequence containing gaps
#' @param coords coordinates on the ungapped refSeq to be converted into equivalent gapped coordinates
#' @param gapChars characters interpreted as gaps
#'
#' @return A vector of coordinates in the gapped reference sequence
#'
#' @export
#'
#' @examples
#' noGapToGap('AA--AA-A',c(1:5))
noGapToGap<-function(refSeq,coords,gapChars=c('*','.','-')){
  gapSeqSplit<-strsplit(refSeq,'')[[1]]
  nonDash<-which(!gapSeqSplit %in% gapChars)
  coords[coords<1|coords>length(nonDash)]<-NA
  return(nonDash[c(0,coords)])
}

#' Convenience function to replace gaps at the start and end of a sequence with a different character
#'
#' Replace gaps at the start and/or end of a sequence with another character
#' (e.g. ---A-A--A--- to ...A-A--A...). This can be used to indicate the
#' difference between indels within known sequences and unknown sequence
#' surrounding a sequence.
#'
#' @param seqs a character vector of sequences
#' @param leftEnd logical indicating whether gaps should be replaced at the start of the sequence
#' @param rightEnd logical indicating whether gaps should be replaced at the end of the sequence
#' @param gapChars a vector of single characters that count as gaps
#' @param replaceChar character to replace start/end gaps with
#'
#' @return A character vector of sequences with start and/or end gaps replaced with replaceChar
#'
#' @export
#'
#' @examples
#' replaceOuterGaps(c('--A-A--','AAA--AAA','--A-A','A-A--'))
replaceOuterGaps<-function(seqs,leftEnd=TRUE,rightEnd=TRUE,gapChars=c('*','-'),replaceChar='.'){
  if(any(is.na(seqs)))stop(simpleError('NA sequence found in replaceOuterGaps'))
  gapRegex<-sprintf('[%s]',paste(gapChars,collapse=''))
  startGapLength<-attr(regexpr(sprintf('^%s+',gapRegex),seqs),'match.length')
  endGapLength<-attr(regexpr(sprintf('%s+$',gapRegex),seqs),'match.length')
  nChars<-nchar(seqs)
  dummy<-paste(rep(replaceChar,max(c(startGapLength,endGapLength))+100),collapse='')
  if(leftEnd)substring(seqs,1,startGapLength)<-substring(dummy,1,startGapLength)
  if(rightEnd)substring(seqs,nChars-endGapLength+1)<-substring(dummy,1,endGapLength)
  return(seqs)
}

#' Convenience function to remove columns composed mostly of gaps
#'
#' Remove columns in an alignment that contain more than a give proportion of gaps
#' (e.g. A-A, G-G to AA,GG). 
#'
#' @param seqs a character vector of sequences
#' @param gapChars a vector of single characters that count as gaps
#' @param maxGapProp remove columns with greater than this proportion of gaps
#' @param ignoreChars a vector of single characters that are ignored i.e. count neither as gaps or towards the total sequence count
#'
#' @return A character vector of sequences with gap columns removed
#'
#' @export
#'
#' @examples
#' removeGapCols(c('A-A-','A-AA','A-AT','A-AG'))
removeGapCols<-function(seqs,gapChars=c('*','-','.'),maxGapProp=.9,ignoreChars=c()){
  if(any(is.na(seqs)))stop(simpleError('NA sequence found in replaceOuterGaps'))
  mat<-seqSplit(seqs)
  gapProp<-apply(mat,2,function(x)mean(x[!x %in% ignoreChars] %in% gapChars))
  gapCols<-gapProp>maxGapProp
  out<-apply(mat[,!gapCols],1,paste,collapse='')
  return(out)
}

#' Convenience function to remove amino acids following a stop codon 
#'
#' Replace any non X characters following a stop codon (marked as X)
#' with Xs (e.g. LSYXAAA to LSYXXXX). This can be used to frame shift
#' mutations and early terminations
#'
#' @param seqs a character vector of sequences
#' @param stopChars a vector of single characters that count as stop codons
#' @param replaceChar character to replace with (default: X)
#'
#' @return A character vector of sequences with amino acids following a stop replaced
#'
#' @export
#'
#' @examples
#' replaceAfterStop(c('LYSXAAA','LYSRAAA','AXAAA','AAAAX'))
replaceAfterStop<-function(seqs,stopChars='X',replaceChar='X'){
  if(any(is.na(seqs)))stop(simpleError('NA sequence found in replaceAfterStop'))
  stopRegex<-sprintf('[%s].*$',paste(stopChars,collapse=''))
  stopLength<-attr(regexpr(stopRegex,seqs),'match.length')-1
  dummy<-paste(rep(replaceChar,max(stopLength)+100),collapse='')
  substring(seqs,nchar(seqs)-stopLength+1)<-substring(dummy,1,stopLength)
  return(seqs)
}



#' Create fake DNA or amino acid sequences
#'
#' Creates a random reference sequence then adds mutations and indels to
#' hierarchical sets of the sequences and random noise
#'
#' @param n number of fake sequences to generate
#' @param nChar character length of the output fake sequences
#' @param nSplit The number of hierarchical splits to make in the data e.g. 3 splits produces 2^3=8 "species"
#' @param pGap probability of a large insertion or deletion in each grouping
#' @param pNoise probability of a random substitution at each base
#' @param pMutation probability of a substitution at each base in each hierarchical grouping
#' @param bases the bases used in generating the sequence (bases listed in excludeBases are excluded from the initial reference sequence generation
#' @param excludeBases bases excluded from the initial reference sequence generation (default: '-' and 'X')
#'
#' @return A n+1 length character vector of fake sequences. The first sequence is the reference. Names of the remaining sequences indicate their hierarchical groupings follow by an arbitrary id
#'
#' @export
#'
#' @examples
#' createFakeDNA(10,10)
createFakeDNA<-function(n=500,nChar=400,nSplit=3,pGap=.3,pNoise=.01,pMutation=.005,bases=c('A','C','T','G','-'),excludeBases=c('-','X')){
  if(nSplit==0){
    nSplit<-1
    pGap=0
    pMutation<-0
  }
  refSeq<-sample(bases[!bases %in% excludeBases],nChar,TRUE)
  seqMat<-matrix(refSeq,nrow=n,ncol=nChar,byrow=TRUE)
  groupAssign<-list(rep('0',n))
  for(ii in 1:nSplit){
    #groupSplits<-1:n%%2^ii
    groupSplits<-stats::ave(groupAssign[[ii]],groupAssign[[ii]],FUN=function(x){
      pGroup<-exp(stats::rnorm(2))
      pGroup<-pGroup/sum(pGroup)
      paste(x,sample(0:1,length(x),TRUE,pGroup),sep='')
    })
    groupAssign[[ii+1]]<-groupSplits
    for(jj in unique(groupSplits)){
      nSubs<-stats::rbinom(1,nChar,pMutation)  
      selector<-groupSplits==jj
      seqMat[selector,sample(1:nChar,nSubs)]<-matrix(sample(bases,nSubs,TRUE),nrow=sum(selector),ncol=nSubs,byrow=TRUE)
    }
  }
  nNoise<-stats::rbinom(1,n*nChar,pNoise)
  randomNoiseLoc<-sample(1:(n*nChar),nNoise)
  seqMat[randomNoiseLoc]<-sample(bases,nNoise,TRUE)
  if(pGap>0){
    for(ii in 1:nSplit){
      groupSplits<-groupAssign[[ii+1]]
      for(jj in unique(groupSplits)){
        selector<-groupSplits==jj
        if(stats::runif(1)<pGap){
          gapStart<-sample(1:nChar,1)
          gapEnd<-min(nChar,gapStart+sample(floor(nChar/30):ceiling(nChar/10),1))
          if(stats::runif(1)<.5){
            #deletion
            seqMat[selector,gapStart:gapEnd]<-'-'
          } else {
            #insertion
            seqMat[!selector,gapStart:gapEnd]<-'-' 
            refSeq[gapStart:gapEnd]<-'-' 
          }
        }
      }
    }
  }
  rownames(seqMat)<-sprintf('%s %d',groupAssign[[nSplit+1]],1:n)
  seqMat<-seqMat[do.call(order,c(groupAssign,list(apply(seqMat,1,paste,collapse='')))),]
  out<-apply(rbind(refSeq,seqMat),1,paste,collapse='')
  return(out)
}

#' @rdname createFakeDNA 
#' @param ... additional arguments to createFakeDNA
#' @export
#' @examples
#' createFakeAA(10,10)
createFakeAA<-function(n=100,nChar=100,...,pGap=.2,pNoise=0,pMutation=.01,bases=c(names(dnaplotr::aminoCols),'-')){
  createFakeDNA(n,nChar,bases=bases,pGap=pGap,pNoise=pNoise)
}

#' Some colors for amino acids
#'
#' A vector indexed by single letter amino acid code giving color codes based on colors provided by Jmol tweaked slightly so that no two colors are identical
#'
#' @docType data
#' @format A vector with each element giving a color for an amino acid. Vector names correspond to single letter amino acid codes.
#' @references \url{http://jmol.sourceforge.net/jscolors/}
#' @source system.file("data-raw", "makeAminoColors.R", package = "dnaplotr")
"aminoCols"
