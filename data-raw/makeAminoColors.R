
#http://jmol.sourceforge.net/jscolors/
aminoColors<-data.frame('code'=c("D","E","C","M","K","R","S","T","F","Y","N","Q","G","L","V","I","A","W","H","P"),'abbr'=c("ASP","GLU","CYS","MET","LYS","ARG","SER","THR","PHE","TYR","ASN","GLN","GLY","LEU","VAL","ILE","ALA","TRP","HIS","PRO"),'col'=c("#E60A0A","#E60A0A","#E6E600","#E6E600","#145AFF","#145AFF","#FA9600","#FA9600","#3232AA","#3232AA","#00DCDC","#00DCDC","#EBEBEB","#0F820F","#0F820F","#0F820F","#C8C8C8","#B45AB4","#8282D2","#DC9682"),stringsAsFactors=FALSE)
tmpAngles1<-cos(1+1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles2<-sin(1:nrow(aminoColors)/nrow(aminoColors)*pi)
tmpAngles1<-tapply(tmpAngles1,aminoColors$col,c)
tmpAngles2<-tapply(tmpAngles2,aminoColors$col,c)
aminoColors$spreadCol<-ave(aminoColors$col,aminoColors$col,FUN=function(x){
	if(length(x)==1)return(x)
	spacer<-20
	angles1<-tmpAngles1[[x[1]]]
	angles2<-tmpAngles2[[x[1]]]
	spacing<-seq((length(x)-1)*-spacer,(length(x)-1)*spacer,length.out=length(x))
	lab<-convertColor(t(col2rgb(x[1])),from='sRGB',to='Lab',scale.in=255)[rep(1,length(x)),]
	lab[,'b']<-lab[,'b']-spacing*angles1
	lab[,'a.x']<-lab[,'a.x']+spacing*angles2
	rgbs<-convertColor(lab,from='Lab',to='sRGB')
	return(rgb(rgbs))
})
aminoCols<-aminoColors$spreadCol
names(aminoCols)<-aminoColors$code
aminoCols['-']<-'white'
aminoCols['X']<-'black'

save(aminoCols,file='data/aminoColors.RData')
tools::resaveRdaFiles('data/aminoColors.RData')
