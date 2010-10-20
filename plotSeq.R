#Function: plotSeq(c('ACTAAATTT','GGTAAGTTT'), 'test.png')
	#Produce a red, green, blue, yellow plot of sequences

#Arguments:
	#seqs = a vector of sequences (strings)
	#outFile = a file name (.eps and .png used to determine file type, NULL for no file generation and dev.off)
	#distOrder = order sequences by Levenshtein distance from refSeq if defined otherwise most abundant sequence? (requires levenshtein.R)
	#homoLimit = If calculating distance, ignore difference in homopolymers longer than homoLimt
	#emptyTrim = Delete any columns with all -, *, .'s 
	#gapTrim = Delete any columns with fewer than or equal  gapTrim non[-*.] chars
	#groups = Group sequences by group and show lable on right side of plot
	#groupTrim = Delete any groups with sequence counts <= groupTrim
	#groupOrdering = preferred order for groups
	#distShow = Show distances on right side of plot (can be messy when few sequences for given distances)
	#vocal = Show message every vocal sequences when calculating levenshtein (for monitoring progress on long calculations)
	#legend = Show A C T G legend?
	#endGapRemove = Change -'s on ends to .'s  (e.g --AAA-CCC-- becomes ..AAA-CC..) (dots are displayed white and dashes grey)
	#orderBy = Column used to order sequences
	#pause = Call browser() near end for manual adjustments/debugging
	#plotPng = Output file as png or eps (set automatically if filenames end in .eps or .png)
	#extraCmds = Extra commands added after plotting sequences (e.g. lines arrows etc)
	#xstart = First base should be numbered as xstart
	#distOrderDecreasing = Sort by distance in decreasing order?
	#refSeq = Reference sequence used for determining gap free base position and distance
	#res = resolution for png's 1600x900 * res
	#groupCex = cex for group labels
	#lineStagger = offset group labels
	#groupCexScale = scale group label cex by number of sequences?
	#convertGap2NoGap = Use refSeq to determine nonGap base positions for x-axis? requires dna.R
	#seqCounts = vector of counts 
	#fixedAxis = vector of x-axis label positions (for fine-tuning axis labelling)
	#refGapWhite = Make '-' characters white where refSeq also gap
	#noText = Suppress margin text (e.g. for embedding somewhere else)
	#verticalLines = Nongapped refSeq (if refSeq non null) or plain coordinates to draw vertical dotted lines _after_
	#verticalLty = Line type for vertical lines (NULL to calculate but not plot)
	#xlab = label for x axis
	#ylab = label for y axis
	#noTick = suppress ticks?
	#cache = do not create plot if file already exists
	#seqCountDisplay = display left sequence count axis?
	#... = arguments passed to plot()

#Returns: invisible logical vector indicating whether a columns was plotted
#Side effect: Produces plot in outFile

plotSeq<-function(seqs,outFile="test.eps",distOrder=FALSE,homoLimit=0,emptyTrim=TRUE,gapTrim=0,groups=NULL,groupTrim=0,groupOrdering=c(),distShow=TRUE,vocal=2,legend=!noText,endGapRemove=FALSE,orderBy=NULL,pause=FALSE,plotPng=FALSE,extraCmds=NULL,xstart=1,distOrderDecreasing=FALSE,refSeq=NULL,res=1,groupCex=NULL,lineStagger=FALSE,groupCexScale=FALSE,convertGap2NoGap=FALSE,seqCounts=rep(1,length(seqs)),fixedAxis=NULL,refGapWhite=FALSE,noText=FALSE,verticalLines=NULL,verticalLty=2,xlab='Position',ylab='Sequence Read',noTick=FALSE,cache=FALSE,seqCountDisplay=TRUE,...){
	if(cache&&file.exists(outFile))return('CACHED')
	gapChars<-c('-','*','.') #need to standardize throughout
	if(any(grep('.png$',outFile)))plotPng=TRUE
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	if(length(seqs)!=length(seqCounts))stop(simpleError('Lengths of seqs and seqCounts not equal'))
	seqs<-toupper(seqs)
	if(!is.null(refSeq))refGaps<-strsplit(refSeq,'')[[1]] %in% gapChars

	if(groupTrim>0){
			groupNums<-tapply(seqCounts,groups,sum)
			groupNums<-groupNums[groupNums>groupTrim]
			seqs<-seqs[groups %in% names(groupNums)]
			groups<-groups[groups %in% names(groupNums)]
	}

	#ordering
	if(distOrder){
		source('~/scripts/R/levenshtein.R')
		#pick out most common sequence
		if(is.null(refSeq)){
			seqCount<-tapply(seqCounts,seqs,sum)
			maxSeq<-names(seqCount)[seqCount==max(seqCount)][1]
		}else maxSeq<-refSeq
		dists<-levenStringsToStrings(gsub('[*.-]+','',maxSeq),gsub('[-.*]','',seqs,perl=TRUE),substring1=TRUE,homoLimit=homoLimit,vocal=vocal,substring2=TRUE)
		multiplier<-ifelse(distOrderDecreasing,-1,1)
		distRank<-rank(dists*multiplier)
	}else distRank<-rep(0,length(seqs))
	if(!is.null(groups)){
		if(length(groupOrdering)==0)groupRank<-rank(sub('^\\^','0',sub('^\\$','Z',groups)))
		else groupRank<-orderIn(groups,groupOrdering,orderFunc=rank)
	}else{
		groupRank<-rep(0,length(seqs))
	}
	if(!is.null(orderBy)){
		if(!is.list(orderBy))orderBy<-list(orderBy)
		orderByRank<-do.call(rank,orderBy)
	} else orderByRank<-rep(0,length(seqs))
	seqRank<-rank(gsub('[.*-]','Z',seqs))
	if(any(c(distRank,orderByRank,groupRank)!=0))thisOrder<-order(groupRank,orderByRank,distRank,seqRank)
	else thisOrder<-1:length(seqs)

	seqList<-strsplit(seqs,'')
	lengths<-sapply(seqList,length)
	if(any(lengths[1]!=lengths)){
		maxLength<-max(lengths)
		seqList<-lapply(seqList,function(x,num){return(c(x,rep('-',maxLength-length(x))))},maxLength)
	}
	seqMat<-do.call(rbind,seqList)
	if(distOrder|!is.null(groups)|!is.null(orderBy))seqMat<-seqMat[thisOrder,,drop=FALSE]
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
	axisCex<-ifelse(plotPng,3,1.6)
	groupCex<-ifelse(is.null(groupCex),axisCex,groupCex)
	#if(!is.null(groups))seqNum<-seqNum[order(groups),]
	if(!is.null(outFile)){
		if(plotPng) png(outFile,width=round(1600*res/2),height=round(900*res/2),res=80*res/2,type='cairo',antialias='subpixel')
		else postscript(outFile,horizontal=FALSE,width=10,height=6,paper='special')
	}
		#add some space to the right margin if annotating groups or distance
		marRightPad<-ifelse(is.null(groups),ifelse(distShow,3,0),max(nchar(groups))*1.05)
		if(plotPng){
			mars<-c(6,5.1+digits*1.06,1,4+marRightPad)
		} else {
			mars<-c(4.1,2+digits*1.06,1,0+marRightPad/1.2)
		}
		if(!seqCountDisplay)mars[2]<-.2
		par(mar=mars,las=1)
		plot(1,1,xlim=c(0.5,ncol(seqNum)+.5),ylim=c(0.5,sum(seqCounts)+.5),ylab="",xlab=ifelse(noText,'',xlab),type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',cex.axis=axisCex,cex.lab=axisCex,mgp=c(mars[1]-1.5,ifelse(plotPng,1,.75),0),...)
		prettyY<-pretty(1:sum(seqCounts))
		if(!noTick&seqCountDisplay)axis(2,prettyY,ifelse(rep(noText,length(prettyY)),rep('',length(prettyY)),format(prettyY,scientific=FALSE)),cex.axis=axisCex,mgp=c(3,ifelse(plotPng,1,.75),0))
		if(!noText)mtext(ylab,2,line=ifelse(plotPng,2.8,0)+digits^1.03*1.05,las=3,cex=axisCex)

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
			if(!noTick)axis(1,noGap2Gap(refSeq,prettyX-xstart),ifelse(rep(noText,length(prettyX)),rep('',length(prettyX)),prettyX),cex.axis=axisCex,mgp=c(3,ifelse(plotPng,ifelse(plotPng,1.6,1),1),0))
		}else{
			prettyX<-pretty(xstart+c(1,ncol(seqNum)))
			if(!noTick)axis(1,prettyX-xstart,ifelse(rep(noText,length(prettyX)),rep('',length(prettyX)),prettyX),cex.axis=axisCex,mgp=c(3,ifelse(plotPng,1.6,1),0))
		}
		#needs to be slight overlap to avoid stupid white line problem
		spacer<-.1
		for(i in 1:ncol(seqNum)){
			#rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
			rect(i-.5,c(0,cumsum(seqCounts[thisOrder])[-length(seqCounts)])+.5,i+.5,cumsum(seqCounts[thisOrder])+.5+spacer,col=seqNum[,i],border=NA)
		}
		if(distOrder&is.null(groups)&distShow){
			dists<-dists[thisOrder]
			cumSums<-cumsum(seqCounts[thisOrder])
			#first bin goes to 1
			cumSums[1]<-1
			for(i in unique(dists)){
				mtext(i,4,at=cumSums[min(which(dists==i))],cex=2)	
			}
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
		if(!is.null(verticalLines)){
			if(!is.null(refSeq))verticalLines<-noGap2Gap(refSeq,verticalLines-xstart)+.5
			if(!is.null(verticalLty))segments(verticalLines,.5,verticalLines,sum(seqCounts)+.5,lty=verticalLty,lwd=1)
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
	if(!is.null(outFile))dev.off()
	invisible(gapSelector)
}
