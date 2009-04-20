plotSeq<-function(seqs,outFile="test.eps",segs=FALSE,groups=NULL,distOrder=FALSE,homoLimit=0,emptyTrim=TRUE,gapTrim=0,groupOrder=NULL,groupTrim=0,distShow=TRUE,vocal=0,legend=TRUE,endGapRemove=FALSE,orderBy=NULL,pause=FALSE,plotPng=FALSE,extraCmds=NULL,xstart=0,distOrderDecreasing=FALSE,refSeq=NULL,res=1,groupCex=3,lineStagger=FALSE,groupCexScale=FALSE,...){
	if(any(grep('.png$',outFile)))plotPng=TRUE
	if(length(seqs)<1|is.null(seqs))stop(simpleError("Seqs missing"))
	seqs<-toupper(seqs)
	if(groupTrim>0){
			groupNums<-tapply(1:length(groupOrder),groupOrder,length)
			groupNums<-groupNums[groupNums>groupTrim]
			seqs<-seqs[groupOrder %in% names(groupNums)]
			groupOrder<-groupOrder[groupOrder %in% names(groupNums)]
	}
	if(distOrder){
		source('~/scripts/R/levenshtein.R')
		if(is.null(refSeq)){
			seqCount<-tapply(seqs,seqs,length)
			maxSeq<-names(seqCount)[seqCount==max(seqCount)][1]
		}else maxSeq<-refSeq
		dists<-levenOnetoMany(gsub('[*.-]+','',maxSeq),gsub('[-.*]','',seqs,perl=TRUE),subString=TRUE,homoLimit=homoLimit,vocal=vocal,subBoth=TRUE)
		thisOrder<-order(dists,seqs)
	}
	if(!is.null(groupOrder)){
		if(distOrderDecreasing) multiplier<- -1
		else multiplier<-1
		if(distOrder)thisOrder<-order(groupOrder,dists*multiplier,seqs)
		else{
			if (is.null(orderBy))thisOrder<-order(groupOrder,gsub('[.*-]','Z',seqs))
			else thisOrder<-order(groupOrder,orderBy,gsub('[.*-]','Z',seqs))
		}
	}
	seqList<-strsplit(seqs,'')
	lengths<-unlist(lapply(seqList,length))
	if(any(lengths[1]!=lengths)){
		maxLength<-max(lengths)
		seqList<-lapply(seqList,function(x,num){return(c(x,rep('-',maxLength-length(x))))},maxLength)
	}
	seqMat<-do.call(rbind,seqList)
	if(distOrder|!is.null(groupOrder))seqMat<-seqMat[thisOrder,]
	if(emptyTrim)seqMat<-seqMat[,!apply(seqMat,2,function(x){all(x=='-'|x=='.'|x=='*')})]
	if(gapTrim>0)seqMat<-seqMat[,apply(seqMat,2,function(x){sum(x!='-'&x!='.')})>gapTrim]
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
	digits<-ceiling(log10(length(seqs)+1))
	if(!is.null(groups))seqNum<-seqNum[order(groups),]
	if(plotPng) png(outFile,width=round(1600*res),height=round(900*res),res=80*res,type='cairo',antialias='subpixel')
	else postscript(outFile,horizontal=FALSE,width=75,height=25,paper='special')
		par(mar=c(7,5+digits*1,1,7),las=1)
		plot(1,1,xlim=c(0.5,ncol(seqNum)+.5),ylim=c(0.5,nrow(seqNum)+.5),ylab="",xlab="Position",type='n',xaxs='i',yaxs='i',xaxt='n',cex.axis=3,cex.lab=3,mgp=c(6,1.2,0),...)
		mtext('Sequence Read',2,line=3+digits,las=3,cex=3)
		prettyX<-pretty(xstart+c(1,ncol(seqNum)))
		axis(1,prettyX-xstart,prettyX,cex.axis=3)
		#needs to be slight overlap to avoid stupid white line problem
		spacer<-.1
		for(i in 1:nrow(seqNum)){
			rect(1:ncol(seqNum)-.5,i-.5,1:ncol(seqNum)+.5,i+.5+spacer,col=seqNum[i,],border=NA)
		}
		#segments(0:ncol(seqNum)+.5,0.5,0:ncol(seqNum)+.5,nrow(seqNum)+.5,lwd=.1)
		#if(segs)segments(0.5,0:nrow(seqNum)+.5,ncol(seqNum)+.5,0:nrow(seqNum)+.5,lwd=.1)
		if(distOrder&is.null(groupOrder)&distShow){
			dists<-dists[thisOrder]
			for(i in unique(dists)){
				mtext(i,4,at=min(which(dists==i)),cex=2)	
			}
		}
		if(!is.null(groupOrder)){
			groups<-groupOrder[thisOrder]
			counter<-0
			maxGroupCount<-max(table(groupOrder))
			print(maxGroupCount)
			for(i in sort(unique(groupOrder))){
				counter<-counter+1
				thisMin<-min(which(groups==i))
				thisMax<-max(which(groups==i))
				if(lineStagger)line=(counter-1)*.3
				else line=.5
				if(groupCexScale)cexScale<-((diff(c(thisMin,thisMax))+1)/maxGroupCount)^.5
				else cexScale<-1
				print(cexScale)
				mtext(i,4,at=mean(c(thisMin,thisMax)),cex=max(.3,cexScale*groupCex),line=line)
				segments(-.5,thisMin-.5,ncol(seqNum)+.5,thisMin-.5)
				segments(-.5,thisMax+.5,ncol(seqNum)+.5,thisMax+.5)
			}
		}
		if(!is.null(extraCmds))eval(parse(text=extraCmds))
		box()
		if(legend){
			ypos<- -nrow(seqNum)*.04
			xpos<-ncol(seqNum)*.98	
			#adj=c(1,0)
			legend(xpos, ypos, c("A", "T", "C","G"),col=c('green','red','blue','yellow'), pt.bg= c('green','red','blue','yellow'),pch = c(22,22,22,22),ncol=4,bty='n',cex=2,xjust=1,yjust=1,xpd=NA)
		}
		if(pause)browser()
	dev.off()
	#num<-nrow(seqMat)
	#As<-apply(seqMat,2,function(x){return(sum(x=='A'))})/num
	#Ts<-apply(seqMat,2,function(x){return(sum(x=='T'))})/num
	#Cs<-apply(seqMat,2,function(x){return(sum(x=='C'))})/num
	#Gs<-apply(seqMat,2,function(x){return(sum(x=='G'))})/num
	#dashes<-apply(seqMat,2,function(x){return(sum(x=='-'))})/num
	#percs<-rbind(As,Ts,Cs,Gs,dashes)
	#browser()
}

multiPlotSeq<-function(seqs,amps,outputPath='',segs=FALSE,groups=NULL){
	for (i in unique(amps)){
		if(!is.null(groups))thisGroups<-groups[amps==i]
		else thisGroups<-NULL
		plotSeq(seqs[amps==i],outFile=paste(outputPath,i,'.eps',sep=""),segs=segs,groups=thisGroups)
	}
}

