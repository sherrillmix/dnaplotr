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
#' @name plotDNA
#' @author Scott Sherrill-Mix, \email{shescott@@upenn.edu}
#' @examples
#' seqs<-c('ACTGGA','ACTGCA','ACTGGC')
#' plotDNA(seqs)
NULL


#' Find contiguous ranges within a vector of indices
#' 
#' Takes a vector of indices and returns ranges of contiguous regions
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
#' index2range(c(1:10,11,15,16,17:20))
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
	#emptyTrim = Delete any columns with all -, *, .'s 
	#gapTrim = Delete any columns with fewer than or equal  gapTrim non[-*.] chars
	#groups = Group sequences by group and show lable on right side of plot
	#groupOrdering = preferred order for groups
	#legend = Show A C T G legend?
	#endGapRemove = Change -'s on ends to .'s  (e.g --AAA-CCC-- becomes ..AAA-CC..) (dots are displayed white and dashes grey)
	#orderBy = Column used to order sequences
	#pause = Call browser() near end for manual adjustments/debugging
	#extraCmds = Extra commands added after plotting sequences (e.g. lines arrows etc)
	#xstart = First base should be numbered as xstart
	#distOrderDecreasing = Sort by distance in decreasing order?
	#refSeq = Reference sequence used for determining gap free base position and distance
	#groupCex = cex for group labels
	#lineStagger = offset group labels
	#groupCexScale = scale group label cex by number of sequences?
	#convertGap2NoGap = Use refSeq to determine nonGap base positions for x-axis? requires dna.R
	#seqCounts = vector of counts 
	#fixedAxis = vector of x-axis label positions (for fine-tuning axis labelling)
	#refGapWhite = Make '-' characters white where refSeq also gap
	#noText = Suppress margin text (e.g. for embedding somewhere else)
	#xlab = label for x axis
	#ylab = label for y axis
	#noTick = suppress ticks? can be length 1 for both or 2 for x then y
	#seqCountDisplay = display left sequence count axis?
	#maxAxis = maximum lab to display on y axis (e.g. not to scale stuff on top of this)
	#... = arguments passed to plot()

#Returns: invisible logical vector indicating whether a columns was plotted
#Side effect: Produces plot in outFile

plotSeq<-function(seqs,emptyTrim=TRUE,gapTrim=0,groups=NULL,groupOrdering=c(),legend=!noText,endGapRemove=FALSE,orderBy=NULL,pause=FALSE,extraCmds=NULL,xstart=1,distOrderDecreasing=FALSE,refSeq=NULL,groupCex=NULL,lineStagger=FALSE,groupCexScale=FALSE,convertGap2NoGap=FALSE,seqCounts=rep(1,length(seqs)),fixedAxis=NULL,refGapWhite=FALSE,noText=FALSE,xlab='Position',ylab='Sequence Read',noTick=FALSE,seqCountDisplay=TRUE,maxAxis=Inf,...){
	if(length(noTick)==1)noTick<-rep(noTick,2)
	gapChars<-c('-','*','.') #need to standardize throughout
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	if(length(seqs)!=length(seqCounts))stop(simpleError('Lengths of seqs and seqCounts not equal'))
	seqs<-toupper(seqs)
	if(!is.null(refSeq))refGaps<-strsplit(refSeq,'')[[1]] %in% gapChars

	#ordering
	if(!is.null(groups)){
		if (any(grep('^[^$]',groups)))dummy<-paste(ifelse(substring(groups,1,1)=='^','0','1'),ifelse(substring(groups,1,1)=='$','1','0'),sub('^\\^','0',sub('^\\$','Z',groups)),sep='')
		else dummy<-groups
		if(length(groupOrdering)==0)groupRank<-rank(dummy)
		else groupRank<-orderIn(groups,groupOrdering,orderFunc=rank)
	}else{
		groupRank<-rep(0,length(seqs))
	}
	if(!is.null(orderBy)){
		if(!is.list(orderBy))orderBy<-list(orderBy)
		orderByRank<-do.call(rank,orderBy)
	} else orderByRank<-rep(0,length(seqs))
	seqRank<-rank(gsub('[.*-]','Z',seqs))
	if(any(c(orderByRank,groupRank)!=0))thisOrder<-order(groupRank,orderByRank,seqRank)
	else thisOrder<-1:length(seqs)

	seqList<-strsplit(seqs,'')
	lengths<-sapply(seqList,length)
	if(any(lengths[1]!=lengths)){
		maxLength<-max(lengths)
		seqList<-lapply(seqList,function(x,num){return(c(x,rep('-',maxLength-length(x))))},maxLength)
	}
	seqMat<-do.call(rbind,seqList)
	gapSelector<-rep(TRUE,ncol(seqMat))
	if(emptyTrim){
		selector<-!apply(seqMat,2,function(x){all(x %in% gapChars)})
		if(convertGap2NoGap & !is.null(refSeq)){
			edges<-range(which(selector))
			selector[edges[1]:edges[2]]<-selector[edges[1]:edges[2]]|!refGaps[edges[1]:edges[2]]
			if(edges[1]>1)xstart<-xstart+sum(!selector[1:(edges[1]-1)]&!refGaps[1:(edges[1]-1)])
		}
		seqMat<-seqMat[,selector,drop=FALSE]
		if(!is.null(refSeq))refSeq<-paste(strsplit(refSeq,'')[[1]][selector],collapse='')
		gapSelector<-selector
	}
	if(gapTrim>0){
		selector<-apply(seqMat,2,function(x,y){sum(!rep(x,y) %in% gapChars)},seqCounts)>gapTrim
		if(!is.null(refSeq)){
			selector<-selector|!strsplit(refSeq,'')[[1]] %in% gapChars
			refSeq<-paste(strsplit(refSeq,'')[[1]][selector],collapse='')
		}
		seqMat<-seqMat[,selector,drop=FALSE]
		gapSelector[gapSelector]<-selector
	}
	if(endGapRemove){
		seqMat<-t(apply(seqMat,1,function(x){lims<-range(which(x!='-'));x[-(lims[1]:lims[2])]<-'.';return(x)}))
	}
	if(refGapWhite&!is.null(refSeq)){
		seqMat<-t(apply(seqMat,1,function(x,y){x[x=='-'&y]<-'*';return(x)},strsplit(refSeq,'')[[1]] %in% gapChars))
	}
	seqNum<-seqMat
	seqNum[,]<- 0
	seqNum[seqMat=='A']<-'green'
	seqNum[seqMat=='T']<-'red'
	seqNum[seqMat=='C']<-'blue'
	seqNum[seqMat=='G']<-'yellow'
	seqNum[seqMat=='-']<-'grey'
	digits<-ceiling(log10(sum(seqCounts)+1))
	digits<-digits+(ceiling(digits/3)-1) #account for , in 1,000
	axisCex<-1.6
	groupCex<-ifelse(is.null(groupCex),axisCex,groupCex)
	#if(!is.null(groups))seqNum<-seqNum[order(groups),]
		#add some space to the right margin if annotating groups or distance
		marRightPad<-ifelse(is.null(groups),0,max(c(nchar(groups),0))*2*groupCex/3)
		mars<-c(4.1,2+digits*1.06,1,0+marRightPad/1.2)
		if(!seqCountDisplay)mars[2]<-.2
		par(mar=mars,las=1)
		plot(1,1,xlim=c(0.5,ncol(seqNum)+.5),ylim=c(0.5,sum(seqCounts)+.5),ylab="",xlab=ifelse(noText,'',xlab),type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',cex.axis=axisCex,cex.lab=axisCex,mgp=c(mars[1]-1.5,.75,0),...)
		prettyY<-pretty(1:min(maxAxis,sum(seqCounts)))
		prettyY<-prettyY[prettyY<=maxAxis&round(prettyY)==prettyY]
		if(!noTick[2]&seqCountDisplay)axis(2,prettyY,ifelse(rep(noText,length(prettyY)),rep('',length(prettyY)),format(prettyY,scientific=FALSE,big.mark=',')),cex.axis=axisCex,mgp=c(3,.75,0))
		if(!noText)mtext(ylab,2,line=digits^1.03*1.05,las=3,cex=axisCex)

		#Converting to first base as 0 for ease of use
		xstart<-xstart-1
		if(convertGap2NoGap&!is.null(refSeq)){
			if(!exists('gap2NoGap'))source('~/scripts/R/dna.R')
			maxNoGap<-gap2NoGap(refSeq,ncol(seqNum))
			prettyX<-pretty(xstart+c(1,maxNoGap))
			prettyX<-prettyX[prettyX<=xstart+maxNoGap]
			prettyX[prettyX==0]<-1
			prettyX<-prettyX[prettyX-xstart>0]
			if(length(prettyX<4)){prettyX<-pretty(prettyX);}
			if(!is.null(fixedAxis))prettyX<-fixedAxis	
			prettyX<-prettyX[prettyX<=xstart+maxNoGap]; prettyX[prettyX==0]<-1
			#axis(1,noGap2Gap(refSeq,prettyX-xstart),noGap2Gap(refSeq,prettyX-xstart),cex.axis=3)
			if(!noTick[1])axis(1,noGap2Gap(refSeq,prettyX-xstart),ifelse(rep(noText,length(prettyX)),rep('',length(prettyX)),prettyX),cex.axis=axisCex,mgp=c(3,1,0))
		}else{
			prettyX<-pretty(xstart+c(1,ncol(seqNum)))
			if(!noTick[1])axis(1,prettyX-xstart,ifelse(rep(noText,length(prettyX)),rep('',length(prettyX)),prettyX),cex.axis=axisCex,mgp=c(3,1,0))
		}
		#needs to be slight overlap to avoid stupid white line problem
		spacer<-.001
		for(i in 1:ncol(seqNum)){
			#1 rectangle per read
			#rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
			bottoms<-c(0,cumsum(seqCounts[thisOrder])[-length(seqCounts)])
			tops<-cumsum(seqCounts[thisOrder])
			cols<-seqNum[,i]
			#1 rectangle per repped read 
			#rect(i-.5,bottoms+.5,i+.5,tops+.5+spacer,col=cols,border=NA)
			uniqCols<-unique(cols)
			colRanges<-do.call(rbind,lapply(unique(cols),function(x){
				out<-index2range(which(cols==x))
				out$bottom<-bottoms[out$start]
				out$top<-tops[out$end]
				out$col<-x
				return(out)
			}))
			#1 rectangle per string of identical bases
			rect(i-.5,colRanges$bottom+.5,i+.5,colRanges$top+.5+spacer,col=colRanges$col,border=NA)
		}
		if(!is.null(groups)){
			groupOrder<-rep(groups[thisOrder],seqCounts[thisOrder])
			counter<-0
			maxGroupCount<-max(table(groupOrder))
			uniqueGroups<-unique(groups)
			for(i in uniqueGroups){
				counter<-counter+1
				thisMin<-min(which(groupOrder==i))
				thisMax<-max(which(groupOrder==i))
				if(lineStagger)line=(counter-1)*.3
				else line=.5
				if(groupCexScale)cexScale<-((diff(c(thisMin,thisMax))+1)/maxGroupCount)^.5
				else cexScale<-1
				if(!noText)mtext(sub('^[$^]','',i),4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*groupCex),line=line)
				segments(-.5,thisMin-.5,ncol(seqNum)+.5,thisMin-.5)
				segments(-.5,thisMax+.5,ncol(seqNum)+.5,thisMax+.5)
			}
		}
		if(!is.null(extraCmds))eval(parse(text=extraCmds))
		box()
		if(legend){
			ypos<- -sum(seqCounts)*.04
			xpos<-ncol(seqNum)*.98	
			#adj=c(1,0)
			bottomMarHeight<-(par('din')[2]*diff(par('fig')[3:4])-par('pin')[2])/sum(par('mar')[c(1,3)])*par('mar')[1]/par('pin')[2]
			rightMarWidth<-(par('din')[1]*diff(par('fig')[1:2])-par('pin')[1])/sum(par('mar')[c(2,4)])*par('mar')[4]/par('pin')[1]
			insetPos<-c(-rightMarWidth,-bottomMarHeight)
			legend('bottomright', c("A", "T", "C","G"),col=c('green','red','blue','yellow'), pt.bg= c('green','red','blue','yellow'),pch = c(22,22,22,22),ncol=4,bty='n',cex=axisCex*.75,xjust=1,yjust=1,xpd=NA,inset=insetPos)
		}
		if(pause)browser()
	invisible(gapSelector)
}

plotSeqMuscle<-function(seqs,...,extra.args=''){
	library(bio3d)
	if(!exists('seqSplit'))source('~/scripts/R/dna.R')
	aligned<-apply(seqaln(seqSplit(seqs,fill='-'),protein=FALSE,extra.args=extra.args)$ali,1,paste,collapse='')
	plotSeq(aligned,...)
}

plotSeqPair<-function(seqs,...,type='global'){
	library(Biostrings)
	aligned<-pairwiseAlignment(seqs[1],seqs[2],type=type)
	plotSeq(c(as.character(aligned@subject),as.character(aligned@pattern)),...)
}
