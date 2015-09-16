#' Functions to plot out a collection of DNA sequences
#'
#' Produce a visual representation of a vector of DNA sequences
#'
#' The main function is:
#'      \describe{
#'        \item{\code{\link{plotDNA}}:}{to produce a plot from a vector of DNA sequences}
#'      }
##'
#' @docType package
#' @name DNAPlotR
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
#' 
#' @examples
#' DNAPlotR:::indexToRange(c(1:10,11,14,16,17:20))
indexToRange<-function(index){
	index<-sort(unique(index))
	diffs<-c(diff(index),Inf)
	ends<-c(which(diffs>1))
	starts<-c(1,ends[-length(ends)]+1)
	return(data.frame('start'=index[starts],'end'=index[ends]))
}

#Function: plotSeq(c('ACTAAATTT','GGTAAGTTT'), 'test.png') #Produce a red, green, blue, yellow plot of sequences
#Arguments:
	#seqs = a vector of sequences (strings)
	#gapTrim = Delete any columns with fewer than or equal  gapTrim non[-*.] chars
	#groups = Group sequences by group and show lable on right side of plot
	#groupOrdering = preferred order for groups
	#xstart = First base should be numbered as xstart
	#distOrderDecreasing = Sort by distance in decreasing order?
	#refSeq = Reference sequence used for determining gap free base position and distance
	#groupCex = cex for group labels
	#lineStagger = offset group labels
	#groupCexScale = scale group label cex by number of sequences?
	#convertGap2NoGap = Use refSeq to determine nonGap base positions for x-axis? requires dna.R
	#seqCounts = vector of counts 
	#noText = Suppress margin text (e.g. for embedding somewhere else)
	#xlab = label for x axis
	#ylab = label for y axis
	#seqCountDisplay = display left sequence count axis?
	#... = arguments passed to plot()

#Returns: invisible logical vector indicating whether a columns was plotted
#Side effect: Produces plot in outFile

#' Plot a bunch of DNA sequences
#' 
#' Take a vector of strings representing DNA sequences and plot them to the current device. A, C, T and G are colored, - are colored gray and all other characters are white.
#'
#' @param seqs A character vector containing DNA sequences
#' @param seqCounts A integer vector with the number of counts for each sequence. This can be used to improve run time and file size if some sequences are duplicated multiple times (default: 1 for each entry in seqs)
#' @param cols A named vector with names corresponding to the DNA bases and values showing the appropriate color (default: A: green, T: red, C: blue, G: yellow, -: grey)
#' @param xlab A string specifying the x-axis label (default: Position)
#' @param ylab A string specifying the y-axis label (default: Sequence read)
#' @param display A logical vector with element names in 'legend', 'xAxis', 'yAxis', 'groups' where a FALSE suppresses outputting that plot element (default: TRUE for all or any missing elements)
#' @param xStart First base in plot should be labelled as this (default: 1)
#' @param groups Group sequences by group and show label on right side of plot. Note that any prior ordering of sequences will be corrupted
#' @param groupCexScale A logical wheter to scale group label size by the number of sequences. Useful to highlight more abundant groups and help squeeze in labels on smaller groups.
#' @param refSeq Reference sequence used for numbering the x-axis without counting gaps present in this sequence
#' @param ... Additional arguments to plot
#'
#' @return NULL
#'
#' @export
#' 
#' @examples
#' plotDNA(c('ACACA','ACACA','ACACT'))

#things to add back:
# gapTrim
# endGapTrim
# orderBy?
# groups
# refSeq matched - into white
# refseq display
# distOrder

plotDNA<-function(seqs,seqCounts=rep(1,length(seqs)),cols=c('A'='green','T'='red','C'='blue','G'='yellow','-'='grey','default'='white'),xlab='Position',ylab='Sequence Read',display=c('groups'=!is.null(groups)),xStart=1,groups=NULL,groupCexScale=FALSE,
	refSeq=NULL,...){
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	if(length(seqs)!=length(seqCounts))stop(simpleError('Lengths of seqs and seqCounts not equal'))
	if(!is.null(groups)&&length(seqs)!=length(groups))stop(simpleError('Lengths of seqs and groups not equal'))
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

	plot(1,1,xlim=c(0.5,ncol(seqNum)+.5),ylim=c(0.5,sum(seqCounts)+.5),ylab="",xlab=xlab,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',mgp=c(2.6,.75,0),...)
	
	#y axis
	prettyY<-pretty(1:min(sum(seqCounts)))
	prettyY<-prettyY[round(prettyY)==prettyY]
	if(display['yAxis'])axis(2,prettyY,format(prettyY,scientific=FALSE,big.mark=','),mgp=c(3,.75,0))
	title(ylab=ylab,line=3.5,las=3)

	#Converting to first base as 0 for ease of use
	xStart<-xStart-1
	if(!is.null(refSeq)){
		maxNoGap<-gapToNoGap(refSeq,ncol(seqNum))
		prettyX<-pretty(xStart+c(1,maxNoGap))
		prettyX<-prettyX[prettyX<=xStart+maxNoGap]
		prettyX[prettyX==0]<-1
		prettyX<-prettyX[prettyX-xStart>0]
		if(length(prettyX<4)){prettyX<-pretty(prettyX);}
		prettyX<-prettyX[prettyX<=xStart+maxNoGap]; prettyX[prettyX==0]<-1 #axis(1,noGapToGap(refSeq,prettyX-xStart),noGapToGap(refSeq,prettyX-xStart),cex.axis=3)
		if(display['xAxis'])axis(1,noGapToGap(refSeq,prettyX-xStart),prettyX,mgp=c(3,1,0))
	}else{
		prettyX<-pretty(xStart+c(1,ncol(seqNum)))
		if(display['xAxis'])axis(1,prettyX-xStart,prettyX,mgp=c(3,1,0))
	}
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
		rect(ii-.5,colRanges$bottom+.5,ii+.5,colRanges$top+.5+spacer,col=colRanges$col,border=NA)
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
			if(display['groups'])mtext(sub('^[$^]','',ii),4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*par('cex.axis')),line=.5)
			segments(-.5,thisMin-.5,ncol(seqNum)+.5,thisMin-.5)
			segments(-.5,thisMax+.5,ncol(seqNum)+.5,thisMax+.5)
		}
	}
	box()
	if(display['legend']){
		insetPos<-c(grconvertX(1,'nfc','user'),grconvertY(0,'nfc','user')) #-.01 could cause trouble here
		print(insetPos)
		legendCols<-cols[!names(cols) %in% c('default','-')]
		legend(insetPos[1],insetPos[2], names(legendCols),col=legendCols, pt.bg=legendCols,pch = 22,ncol=4,bty='n',xjust=1,yjust=0,xpd=NA,cex=par('cex.axis'),y.intersp=0)
	}
	return(NULL)
}
#seqs<-c('ACTGGA','ACTGCA','ACTGGC','GCTGGG','GGGG')
#par(mar=c(4,5,.5,5),cex.axis=2,cex.lab=2)
#plotDNA(seqs,seqCounts=c(30,rep(10,2),15,5),groups=factor(c('Aaaa','Bbbb','Cccc','Aaaa','Cccc'),levels=c('ZZZ','Bbbb','Aaaa','Cccc')),xStart=10,groupCexScale=TRUE)

#convenience function for binding a bunch of sequences together
#...: various sequences to split into a matrix
#fill: fill sequence to pad ends of sequences

#' Convenience function for binding a bunch of sequences together
#' 
#' Take a vector of strings and return a matrix with each row corresponding to a string and each column a position in those strings. If fill is set then strings are padded to an equal length. Otherwise errors out if strings are unequal length.
#'
#' @param ... Character vectors to turn into a matrix
#' @param fill A character to fill the ends of short strings with. If NULL then all strings must be the same length.
#'
#' @return A matrix with a single character in each cell with rows corresponding to the number of strings and columns corresponding to the position in the string
#'
#' @export
seqSplit<-function(...,fill=NULL){
	seqs<-c(...)
	if(any(is.na(seqs)))stop(simpleError('NA sequences found in seqSplit'))
	seqN<-nchar(seqs)
	maxN<-max(seqN)
	dummy<-paste(rep(fill,maxN-min(seqN)+100),collapse='')
	if(is.null(fill)&&any(seqN!=maxN))stop(simpleError('All sequences not same length'))
	else seqs<-sprintf('%s%s',seqs,substring(dummy,1,maxN-seqN))
	return(do.call(rbind,strsplit(seqs,'')))
}

#' Convenience function to convert gapped coordinates to what the coordinates would be without gaps
#'
#' Take a reference sequence with gaps and coordinates in that gapped reference sequence and translate the coordinates to the corresponding positions in the gap free sequence. Note that positions landing on a gap are converted to the first 5' nongap position e.g. the fifth position of 'G-----G' is converted to 1 and the second position of '--G--GG' is converted to 0
#'
#' @param refSeq the sequence containing gaps
#' @param coords coordinates on the gapped refSeq to be converted into equivalent nongap coordinates
#' @param gapChars characters interpreted as gaps
#'
#' @return A vector of coordinates in the gap free reference sequence
#'
#' @export
gapToNoGap<-function(refSeq,coords,gapChars=c('*','.','-')){
	gapSeqSplit<-strsplit(refSeq,'')[[1]]
	nonDash<-!gapSeqSplit %in% gapChars
	newCoords<-cumsum(nonDash)
	coords[coords<1|coords>length(newCoords)]<-NA
	return(newCoords[c(0,coords)]) #using 0 to prevent x[NA] returning everything
}

#' Convenience function to convert ungapped coordinates in a reference function to what the coordinates would be in the gapped reference sequence
#'
#' Take a reference sequence with gaps and coordinates in the ungapped reference sequence and translate the coordinates to the corresponding positions in the gapped sequence
#' @param refSeq the sequence containing gaps
#' @param coords coordinates on the ungapped refSeq to be converted into equivalent gapped coordinates
#' @param gapChars characters interpreted as gaps
#' @return A vector of coordinates in the gapped reference sequence
noGapToGap<-function(refSeq,coords,gapChars=c('*','.','-')){
	gapSeqSplit<-strsplit(refSeq,'')[[1]]
	nonDash<-which(!gapSeqSplit %in% gapChars)
	coords[coords<1|coords>length(nonDash)]<-NA
	return(nonDash[c(0,coords)])
}

