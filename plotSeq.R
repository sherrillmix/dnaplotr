plotSeq<-function(seqs,outFile="test.eps",segs=FALSE,distOrder=FALSE,homoLimit=0,emptyTrim=TRUE,gapTrim=0,groups=NULL,groupTrim=0,distShow=TRUE,vocal=0,legend=TRUE,endGapRemove=FALSE,orderBy=NULL,pause=FALSE,plotPng=FALSE,extraCmds=NULL,xstart=1,distOrderDecreasing=FALSE,refSeq=NULL,res=1,groupCex=3,lineStagger=FALSE,groupCexScale=FALSE,convertGap2NoGap=FALSE,seqCounts=rep(1,length(seqs)),...){
	if(any(grep('.png$',outFile)))plotPng=TRUE
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	seqs<-toupper(seqs)

	if(groupTrim>0){
			groupNums<-tapply(seqCounts,groups,sum)
			groupNums<-groupNums[groupNums>groupTrim]
			seqs<-seqs[groups %in% names(groupNums)]
			groups<-groups[groups %in% names(groupNums)]
	}
	if(distOrder){
		source('~/scripts/R/levenshtein.R')
		if(is.null(refSeq)){
			seqCount<-tapply(seqCounts,seqs,sum)
			maxSeq<-names(seqCount)[seqCount==max(seqCount)][1]
		}else maxSeq<-refSeq
		dists<-levenOnetoMany(gsub('[*.-]+','',maxSeq),gsub('[-.*]','',seqs,perl=TRUE),subString=TRUE,homoLimit=homoLimit,vocal=vocal,subBoth=TRUE)
		thisOrder<-order(dists,seqs)
	}
	if(!is.null(groups)){
		if(distOrderDecreasing) multiplier<- -1
		else multiplier<-1
		if(distOrder)thisOrder<-order(groups,dists*multiplier,seqs)
		else{
			if (is.null(orderBy))thisOrder<-order(groups,gsub('[.*-]+','Z',seqs))
			else thisOrder<-order(groups,orderBy,gsub('[.*-]','Z',seqs))
		}
	}
	if(!exists('thisOrder'))thisOrder<-1:length(seqs)
	seqList<-strsplit(seqs,'')
	lengths<-unlist(lapply(seqList,length))
	if(any(lengths[1]!=lengths)){
		maxLength<-max(lengths)
		seqList<-lapply(seqList,function(x,num){return(c(x,rep('-',maxLength-length(x))))},maxLength)
	}
	seqMat<-do.call(rbind,seqList)
	if(distOrder|!is.null(groups))seqMat<-seqMat[thisOrder,]
	if(emptyTrim){
		selector<-!apply(seqMat,2,function(x){all(x=='-'|x=='.'|x=='*')})
		seqMat<-seqMat[,selector]
		if(!is.null(refSeq))refSeq<-paste(strsplit(refSeq,'')[[1]][selector],collapse='')
	}
	if(gapTrim>0){
		selector<-apply(seqMat,2,function(x){sum(x!='-'&x!='.')})>gapTrim
		seqMat<-seqMat[,selector]
		if(!is.null(refSeq))refSeq<-paste(strsplit(refSeq,'')[[1]][selector],collapse='')
	}
	if(endGapRemove){
		seqMat<-t(apply(seqMat,1,function(x){lims<-range(which(x!='-'));x[-(lims[1]:lims[2])]<-'.';return(x)}))
	}
	seqNum<-seqMat
	seqNum[,]<- 0
	seqNum[seqMat=='A']<-'green'
	seqNum[seqMat=='T']<-'red'
	seqNum[seqMat=='C']<-'blue'
	seqNum[seqMat=='G']<-'yellow'
	seqNum[seqMat=='-']<-'grey'
	digits<-ceiling(log10(sum(seqCounts)+1))
	#if(!is.null(groups))seqNum<-seqNum[order(groups),]
	if(plotPng) png(outFile,width=round(1600*res),height=round(900*res),res=80*res,type='cairo',antialias='subpixel')
	else postscript(outFile,horizontal=FALSE,width=75,height=25,paper='special')
		par(mar=c(7,5+digits*1,1,7),las=1)
		plot(1,1,xlim=c(0.5,ncol(seqNum)+.5),ylim=c(0.5,sum(seqCounts)+.5),ylab="",xlab="Position",type='n',xaxs='i',yaxs='i',xaxt='n',cex.axis=3,cex.lab=3,mgp=c(6,1.2,0),...)
		mtext('Sequence Read',2,line=3+digits,las=3,cex=3)
		xstart<-xstart-1
		if(convertGap2NoGap&!is.null(refSeq)){
			source('~/scripts/R/dna.R')
			maxNoGap<-gap2NoGap(refSeq,ncol(seqNum))
			prettyX<-pretty(xstart+c(1,maxNoGap))
			prettyX<-prettyX[prettyX<=maxNoGap]
			prettyX[prettyX==0]<-1
			#axis(1,noGap2Gap(refSeq,prettyX-xstart),noGap2Gap(refSeq,prettyX-xstart),cex.axis=3)
			axis(1,noGap2Gap(refSeq,prettyX-xstart),prettyX,cex.axis=3)
		}else{
			prettyX<-pretty(xstart+c(1,ncol(seqNum)))
			axis(1,prettyX-xstart,prettyX,cex.axis=3)
		}
		#needs to be slight overlap to avoid stupid white line problem
		spacer<-.1
		for(i in 1:ncol(seqNum)){
			#rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
			rect(i-.5,c(0,cumsum(seqCounts[thisOrder])[-length(seqCounts)])+.5,i+.5,cumsum(seqCounts[thisOrder])+.5+spacer,col=seqNum[,i],border=NA)
		}
		if(distOrder&is.null(groups)&distShow){
			dists<-dists[thisOrder]
			for(i in unique(dists)){
				mtext(i,4,at=cumsum(seqCounts)[min(which(dists==i))],cex=2)	
			}
		}
		if(!is.null(groups)){
			groupOrder<-rep(groups[thisOrder],seqCounts[thisOrder])
			counter<-0
			maxGroupCount<-max(table(groupOrder))
			for(i in sort(unique(groups))){
				counter<-counter+1
				thisMin<-min(which(groupOrder==i))
				thisMax<-max(which(groupOrder==i))
				if(lineStagger)line=(counter-1)*.3
				else line=.5
				if(groupCexScale)cexScale<-((diff(c(thisMin,thisMax))+1)/maxGroupCount)^.5
				else cexScale<-1
				mtext(i,4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*groupCex),line=line)
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
			legend(xpos, ypos, c("A", "T", "C","G"),col=c('green','red','blue','yellow'), pt.bg= c('green','red','blue','yellow'),pch = c(22,22,22,22),ncol=4,bty='n',cex=2,xjust=1,yjust=1,xpd=NA)
		}
		if(pause)browser()
	dev.off()
}

multiPlotSeq<-function(seqs,amps,outputPath='',segs=FALSE,groups=NULL){
	for (i in unique(amps)){
		if(!is.null(groups))thisGroups<-groups[amps==i]
		else thisGroups<-NULL
		plotSeq(seqs[amps==i],outFile=paste(outputPath,i,'.eps',sep=""),segs=segs,groups=thisGroups)
	}
}

