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
#' DNAPlotR:::index2range(c(1:10,11,14,16,17:20))
index2range<-function(index){
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
#'
#' @return An invisible
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

plotDNA<-function(seqs,seqCounts=rep(1,length(seqs)),cols=c('A'='green','T'='red','C'='blue','G'='yellow','-'='grey','default'='white'),xlab='Position',ylab='Sequence Read',display=c('groups'=!is.null(groups)),xStart=1,groups=NULL,
	refSeq=NULL,lineStagger=FALSE,groupCexScale=FALSE,convertGap2NoGap=FALSE,...){
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	if(length(seqs)!=length(seqCounts))stop(simpleError('Lengths of seqs and seqCounts not equal'))
	seqs<-toupper(seqs)
	displayOptions<-c('legend','xAxis','yAxis','groups')
	missingOptions<-!displayOptions %in% names(display)
	if(any(missingOptions))display[displayOptions[missingOptions]]<-TRUE

	#fill any trailing gaps
	seqList<-strsplit(seqs,'')
	lengths<-sapply(seqList,length)
	if(any(lengths[1]!=lengths)){
		maxLength<-max(lengths)
		seqList<-lapply(seqList,function(x,num){return(c(x,rep('-',maxLength-length(x))))},maxLength)
	}
	seqMat<-do.call(rbind,seqList)
	gapSelector<-rep(TRUE,ncol(seqMat))
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
	if(convertGap2NoGap&!is.null(refSeq)){
		if(!exists('gap2NoGap'))source('~/scripts/R/dna.R')
		maxNoGap<-gap2NoGap(refSeq,ncol(seqNum))
		prettyX<-pretty(xStart+c(1,maxNoGap))
		prettyX<-prettyX[prettyX<=xStart+maxNoGap]
		prettyX[prettyX==0]<-1
		prettyX<-prettyX[prettyX-xStart>0]
		if(length(prettyX<4)){prettyX<-pretty(prettyX);}
		prettyX<-prettyX[prettyX<=xStart+maxNoGap]; prettyX[prettyX==0]<-1
		#axis(1,noGap2Gap(refSeq,prettyX-xStart),noGap2Gap(refSeq,prettyX-xStart),cex.axis=3)
		if(display['xAxis'])axis(1,noGap2Gap(refSeq,prettyX-xStart),prettyX,mgp=c(3,1,0))
	}else{
		prettyX<-pretty(xStart+c(1,ncol(seqNum)))
		if(display['xAxis'])axis(1,prettyX-xStart,prettyX,mgp=c(3,1,0))
	}
	#needs to be slight overlap to avoid stupid white line problem
	spacer<-.001
	for(i in 1:ncol(seqNum)){
		#1 rectangle per read
		#rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
		bottoms<-c(0,cumsum(seqCounts)[-length(seqCounts)])
		tops<-cumsum(seqCounts)
		thisCols<-seqNum[,i]
		#1 rectangle per repped read 
		#rect(i-.5,bottoms+.5,i+.5,tops+.5+spacer,col=cols,border=NA)
		uniqCols<-unique(thisCols)
		colRanges<-do.call(rbind,lapply(unique(thisCols),function(x){
			out<-index2range(which(thisCols==x))
			out$bottom<-bottoms[out$start]
			out$top<-tops[out$end]
			out$col<-x
			return(out)
		}))
		#1 rectangle per string of identical bases
		rect(i-.5,colRanges$bottom+.5,i+.5,colRanges$top+.5+spacer,col=colRanges$col,border=NA)
	}
	if(!is.null(groups)){
		groupOrder<-rep(groups,seqCounts)
		counter<-0
		maxGroupCount<-max(table(groupOrder))
		uniqueGroups<-unique(groups)
		for(ii in uniqueGroups){
			counter<-counter+1
			thisMin<-min(which(groupOrder==ii))
			thisMax<-max(which(groupOrder==ii))
			if(lineStagger)line=(counter-1)*.3
			else line=.5
			if(groupCexScale)cexScale<-((diff(c(thisMin,thisMax))+1)/maxGroupCount)^.5
			else cexScale<-1
			if(display['groups'])mtext(sub('^[$^]','',ii),4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*par('cex.axis')),line=line)
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
	invisible(gapSelector)
}
par(mar=c(4,5,.5,5),cex.axis=2,cex.lab=2)
plotDNA(seqs,seqCounts=rep(10,3),groups=c('A','B','C'),xStart=10)


