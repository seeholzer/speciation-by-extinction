##################################################
#	Description: Figure S2 - plot of mclust and BFD* results
#	1.	Figure S2: Number of Inferred Species vs. Number of Extinction Events	
#	2.	Figure S3: Number of Inferred Species vs. Latitudinal Midpoint of Extinction Event	
#	3.	Sequential, Nested Extinction Scenarios (Not in manuscript)
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github')

#load packages
library(RColorBrewer)
library(colorspace)
library(plyr)
library(stringi)
source('supporting.functions/FUN.add.alpha.R', chdir = TRUE)

#load extinction scenarios used in mclust
#Extinction Scenarios
popdata = read.delim('cranioleuca/cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]

#load summary tables for bfd
bfd = read.delim('cranioleuca/results.bfd.txt',stringsAsFactors=F)

#load summary tables for mclust
mclust = read.delim('cranioleuca/results.mclust.txt',stringsAsFactors=F)
colnames(mclust)[1] = 'scenario'
#only use combined dataset
mclust = mclust[mclust$trait.set %in% 'combined', ]


#get the number of extinction events for each scenario 
scenarios = c(paste0('FULL_v',1:8),paste0('EXT_v',1:28))
tmp = popdata[,scenarios]
nExtinctions = apply(tmp,2,function(x){
	foo = data.frame(group=popdata$group,scenario=x)
	foo = unique(foo)
	nExtinctions = length(which(foo$scenario %in% 'extinct'))
	return(nExtinctions)
})
nExtinctions = data.frame(scenario=names(nExtinctions),nExtinctions=nExtinctions)


#get the midpoint of the extinction event 
# or
#point of the transition from antisiensis to baroni for the full dataset (nExtinctions=0)

x = tmp[,1]
extMidpoint = apply(tmp,2,function(x){
	foo = data.frame(group=popdata$group,scenario=x)
	foo = unique(foo)
	extMidpoint = mean(foo[foo$scenario %in% 'extinct','group'])
	if(is.na(extMidpoint)){	#transition point from ant to bar
		ant = max(foo[foo$scenario %in% 'ant','group'])
		bar = min(foo[foo$scenario %in% 'bar','group'])
		extMidpoint = mean(c(ant,bar))
	}
	return(extMidpoint)
})
extMidpoint = data.frame(scenario=names(extMidpoint), extMidpoint =extMidpoint)
extMidpoint$extMidpoint[is.nan(extMidpoint$extMidpoint)] = NA


#merge nExtinctionEvents and extMidpoint with bfd results
bfd = join_all(list(bfd,nExtinctions,extMidpoint),by='scenario',match='all')

#merge nExtinctionEvents and extMidpoint with mclust results
mclust = join_all(list(mclust,nExtinctions,extMidpoint),by='scenario',match='all')
mclust[mclust$scenario == 'FULL','nExtinctions'] = 0 

#make a key for the latitude of each extMidpoint
foo = aggregate(lat ~ group, data=popdata, mean)
colnames(foo)[1] = 'extMidpoint'
coo = data.frame(extMidpoint=foo$extMidpoint[-length(foo$extMidpoint)] + .5,lat= (diff(foo$lat)/2) + foo$lat[-length(foo$lat)])
lat = rbind(foo,coo)
key  = lat[order(lat[,1]),]

#shades points based on where in the transect the midpoint of the extinction event occurs or the species transition
#colfunc = colorRampPalette(c("#425BD3", "#E8B618"))
colfunc = colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
key$col = colfunc(length(seq(1,9,.5)))

######################################################################
######################################################################
######################################################################
#	1.	Figure S2: Number of Inferred Species vs. Number of Extinction Events	
######################################################################
######################################################################
######################################################################
png('cranioleuca/Fig.S2.png',width=7,height=6,units='in',res=300,bg='transparent')
#dev.new(width=7,height=6)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),
widths=c(7,7), heights=c(3,3))

xlim = c(-0.5,7.5); ylim = c(-50,110)

#bfd
	par(mar=c(3.5,5,0,1))
	plot(bfd[ ,'BF'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	axis(1,at=c(0:7),labels=c(0:7),tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)
	# mtext('Number of Extinctions',1,line=2.5,cex=1.25,col='grey25')

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext("Bayes Factor Support",2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-0.5,105,'A. Genetic Species Delimitation',adj=c(0,0.5),col='grey25')
	
	#add dashed line connecting scenarios in Figure 2
	xy = bfd[bfd$scenario  %in% c('FULL_v5','EXT_v4','EXT_v16','EXT_v24'),c('nExtinctions','BF')]
	xy = xy[order(xy[,1]), ]
	lines(xy[,1],xy[,2],col='grey50',cex=1.5,lty=3)
	
	#shade points based on where in the transect the midpoint of the extinction event occurs
	bfd$col = mapvalues(bfd$extMidpoint,key[,'extMidpoint'], key[,'col'])
	col = rep(add.alpha('grey25',0.75),nrow(bfd)); col[grep('FULL',bfd$scenario)] = bfd$col[grep('FULL',bfd$scenario)]  
	bg = add.alpha(bfd$col,1); bg[grep('FULL',bfd$scenario)] = 'transparent'
	pch = 21 	
	lty = 3
	lwd = rep(0.5,nrow(bfd)); lwd[grep('FULL',bfd$scenario)] = 3  
	x = jitter(bfd$nExtinctions,0.5)
	y = jitter(bfd$BF,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=lwd)


	#Add legend
	tmp = key[key$extMidpoint %in% c(1,2,3,4,5,6,7,8,9),]
	tmp$x = seq(5.5,7,length.out=nrow(tmp))
	tmp$y = 85
	x = tmp$x
	y = tmp$y
	col = 'grey25'
	points(x,y,pch=15,col=tmp$col,cex=1.5)
	text(mean(x),y+22.5,'Species Boundary / Extinction Midpoint',cex=.6,col=col)
	latLabels = round(tmp[tmp$extMidpoint %in% c(1,2,3,4,5,6,7,8,9),'lat'],1)
	latLabels = paste0(format(latLabels, nsmall = 1),'°')
	latX = tmp[tmp$extMidpoint %in% c(1,2,3,4,5,6,7,8,9),'x']
	text(latX,y+10,latLabels,cex=.5,srt=45,adj=c(0.5,0.5),col=col)
	text(min(x),y-7.5,'North',cex=.5,col=col)
	text(tmp[tmp$extMidpoint %in% max(tmp$extMidpoint,na.rm=T),'x'],y-7.5,'South',cex=.5,col=col)
	meanX = mean(c(min(x),tmp[tmp$extMidpoint %in% max(tmp$extMidpoint,na.rm=T),'x']))
	text(meanX,y-7.5,'<————————>',cex=.5,col=col)

	
#mclust
	par(mar=c(4,5,0,1))
	plot(mclust[ ,'deltaBIC'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	axis(1,at=c(0:7),labels=c(0:7),tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)
	mtext('Extinction Magnitude',1,line=2,cex=1.25,col='black')
	mtext('(number of extinct population groups)',1,line=3,cex=.75,col='black')

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	#axis(4,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext(expression(paste(Delta,"BIC")),2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Cluster',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Clusters',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-0.5,105,'B. Phenotypic Species Delimitation',adj=c(0,0.5),col='grey25')
	
	#add dashed line connecting scenarios in Figure 2
	xy = mclust[mclust$scenario  %in% c('FULL','EXT_v4','EXT_v16','EXT_v24'),c('nExtinctions','deltaBIC')]
	xy = xy[order(xy[,1]), ]
	lines(xy[,1],xy[,2],col='grey50',cex=1.5,lty=3)

	mclust$col = mapvalues(mclust$extMidpoint,key[,'extMidpoint'], key[,'col'])
	col = add.alpha('grey25',.75) 
	#col = bfd$col
	#col = rep('green',nrow(bfd))
	bg = mclust$col
	bg[is.na(bg)] = 'grey75'
	pch = 24 	
	lty = 3
	x = jitter(mclust$nExtinctions,0.5)
	y = jitter(mclust$deltaBIC,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=.5)

dev.off()




######################################################################
######################################################################
######################################################################
#	2.	Figure S3: Number of Inferred Species vs. Latitudinal Midpoint of Extinction Event	
######################################################################
######################################################################
######################################################################

png('cranioleuca/Fig.S3.png',width=7,height=6,units='in',res=300,bg='transparent')
#dev.new(width=7,height=6)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),
widths=c(7,7), heights=c(3,3))


#bfd
	bfd$lat = mapvalues(bfd$extMidpoint,key[,'extMidpoint'],key[,'lat'])	
	#range(bfd$lat,na.rm=T)
	xlim = c(-10,-5.75); ylim = c(-50,110)

	par(mar=c(3.5,5,0,1))
	plot(bfd[ ,'BF'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	at = seq(xlim[1],xlim[2],0.5)
	axis(1,at=at,labels=at,tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext("Bayes Factor Support",2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-10.1,105,'A. Genetic Species Delimitation',adj=c(0,0.5),col='grey25')

	#shade points based on number of extinction events
	colfunc = colorRampPalette(c("grey95", "grey5"))
	keyExt = data.frame(nExt=0:7,colExt=colfunc(length(0:7)))
	bfd$colExt = mapvalues(bfd$nExtinctions,keyExt$nExt,keyExt$colExt) 
	col = add.alpha('grey25',.75) 
	#col = bfd$col
	#col = rep('green',nrow(bfd))
	bg = bfd$colExt
	pch = 21 	
	lty = 3
	x = jitter(bfd$lat,0.5)
	y = jitter(bfd$BF,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=.5)
	#text(x,y,bfd$nExtinctions,cex=.5)
	
	#Add legend for shading
	keyExt$x = seq(-9.9,-9.2,length.out=nrow(keyExt))
	keyExt$y = 60
	x = keyExt$x
	y = keyExt$y
	col = 'grey25'
	points(x,y,pch=15,col=keyExt$col,cex=1.5)
	text(mean(x),y+17.5,'N Extinct Pop. Groups ',cex=.75,col=col)
	latLabels = keyExt$nExt
	latX = keyExt$x
	text(latX,y+8,latLabels,cex=.65,adj=c(0.5,0.5),col=col)
	
#mclust
	mclust$lat = mapvalues(mclust$extMidpoint,key[,'extMidpoint'],key[,'lat'])	
	#range(mclust$lat,na.rm=T)
	xlim = c(-10,-5.75); ylim = c(-50,110)

	par(mar=c(4,5,0,1))
	plot(mclust[ ,'deltaBIC'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	at = seq(xlim[1],xlim[2],0.5)
	axis(1,at=at,labels=at,tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)
	mtext('Latitudinal Midpoint of Extinction Event ',1,line=2.5,cex=1,col='black')
	mtext('South ',1,line=2.5,cex=1,col='grey25',at=-10)
	mtext('North ',1,line=2.5,cex=1,col='grey25',at=-6)


	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	#axis(4,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext(expression(paste(Delta,"BIC")),2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Cluster',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Clusters',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-10.1,105,'B. Phenotypic Species Delimitation',adj=c(0,0.5),col='grey25')

	#shade points based on number of extinction events
	colfunc = colorRampPalette(c("grey95", "grey5"))
	mclust$colExt = mapvalues(mclust$nExtinctions,0:7,colfunc(length(0:7))) 
	col = add.alpha('grey25',.75) 
	#col = bfd$col
	#col = rep('green',nrow(bfd))
	bg = mclust$colExt
	pch = 24 	
	lty = 3
	x = jitter(mclust$lat,0.5)
	y = jitter(mclust$deltaBIC,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=.5)
	#text(x,y,mclust$nExtinctions,cex=.5)

	#Add legend for shading
	keyExt$x = seq(-9.9,-9.2,length.out=nrow(keyExt))
	keyExt$y = 60
	x = keyExt$x
	y = keyExt$y
	col = 'grey25'
	points(x,y,pch=15,col=keyExt$col,cex=1.5)
	text(mean(x),y+17.5,'N Extinct Pop. Groups ',cex=.75,col=col)
	latLabels = keyExt$nExt
	latX = keyExt$x
	text(latX,y+8,latLabels,cex=.65,adj=c(0.5,0.5),col=col)
	

dev.off()







######################################################################
######################################################################
######################################################################
#	3.	Sequential, Nested Extinction Scenarios (Not in manuscript)
######################################################################
######################################################################
######################################################################

#
#	GOAL: create a list of sequential scenario names 
#		the represent all possible sequences of 
#		expanding nested extinction scenarios
#	For example - c('FULL','EXT_v4','EXT_v16','EXT_v24') is
#		what is in Figure 2
#	QUESTION: There is a general trend of increasing support for two species
#			with more extinction events
#			Does this hold for extinction trajectories?
#			Which trajectories have the strongest trend? 

#list of groups for each extinction scenario
tmp = apply(popdata[,!colnames(popdata) %in% c('popNumber','pop','long','lat','group')],2,function(x){
	data.frame(unique(data.frame(popdata$group,x)))
})
tmp = join_all(tmp,by='popdata.group')
colnames(tmp) = c('group',colnames(popdata)[6:ncol(popdata)])
tmp = apply(tmp[,-1],2,function(x)tmp$group[x%in%'extinct'])
#convert these numeric sequences into character vectors
x = tmp[[9]]
foo = sapply(tmp,function(x){
			if(length(x) == 0) return(0)	
			if(length(x) == 1) return(x)	
			if(length(x) > 1) return(paste0(min(x),':',max(x))	)
		})
#create a key
foo = data.frame(scenario = names(tmp), set = unlist(foo))

#Extinction Trajectories
#north to south
NS1 = c('2','2:3','2:4','2:5','2:6','2:7','2:8')
NS2 = c('3','3:4','3:5','3:6','3:7','3:8')
NS3 = c('4','4:5','4:6','4:7','4:8')
NS4 = c('5','5:6','5:7','5:8')
NS5 = c('6','6:7','6:8')
NS6 = c('7','7:8')
#south to north
SN1 = c('8','8:7','8:6','8:5','8:4','8:3','8:2')
SN2 = c('7','7:6','7:5','7:4','7:3','7:2')
SN3 = c('6','6:5','6:4','6:3','6:2')
SN4 = c('5','5:4','5:3','5:2')
SN5 = c('4','4:3','4:2')
SN6 = c('3','3:2')
#symetrical
Sym1 = c('3','2:4')
Sym2 = c('4','3:5','2:6')
Sym3 = c('5','4:6','3:7','2:8')
Sym4 = c('6','5:7','4:8')
Sym5 = c('7','6:8')

trajectories = list(NS1=NS1,NS2=NS2,NS3=NS3,NS4=NS4,NS5=NS5,NS6=NS6,SN1=SN1,SN2=SN2,SN3=SN3,SN4=SN4,SN5=SN5,SN6=SN6,Sym1=Sym1,Sym2=Sym2,Sym3=Sym3,Sym4=Sym4,Sym5=Sym5)

#reverse character order in strings for SN trajectories
#need to correspond to character vectors of the foo key
trajectories[grep('SN',names(trajectories))] = sapply(trajectories[grep('SN',names(trajectories))],function(x) stri_reverse(x))

#List of Trajectories of sequential, expanding extinction scenarios
trajectories = sapply(trajectories,function(x) mapvalues(x,foo[,2],foo[,1]))

#add in the FULL scenarios
#	two FULL scenarios can begin any given trajectory
#	because the first population to go extinct in a trajectory
#	can be assigned as either ant or bar in the FULL scenario
#	the first FULL scenario corresponds to bar being assigned to this population
#	the second FULL scenario corresponds to ant being assigned to this population
NS1 = c('FULL_v1','FULL_v2')
NS2 = c('FULL_v2','FULL_v3')
NS3 = c('FULL_v3','FULL_v4')
NS4 = c('FULL_v5','FULL_v6')
NS5 = c('FULL_v6','FULL_v7')
NS6 = c('FULL_v7','FULL_v8')
#south to north
SN1 = c('FULL_v7','FULL_v8')
SN2 = c('FULL_v6','FULL_v7')
SN3 = c('FULL_v5','FULL_v6')
SN4 = c('FULL_v4','FULL_v5')
SN5 = c('FULL_v3','FULL_v4')
SN6 = c('FULL_v2','FULL_v3')
#symetrical
Sym1 = c('FULL_v2','FULL_v3')
Sym2 = c('FULL_v3','FULL_v4')
Sym3 = c('FULL_v5','FULL_v6')
Sym4 = c('FULL_v6','FULL_v7')
Sym5 = c('FULL_v7','FULL_v8')

FULLkey = data.frame(rbind(NS1=NS1,NS2=NS2,NS3=NS3,NS4=NS4,NS5=NS5,NS6=NS6,SN1=SN1,SN2=SN2,SN3=SN3,SN4=SN4,SN5=SN5,SN6=SN6,Sym1=Sym1,Sym2=Sym2,Sym3=Sym3,Sym4=Sym4,Sym5=Sym5))
colnames(FULLkey) = c('bar','ant')
FULLkey$trajectory = row.names(FULLkey)

trajectoriesBar = list() #list of trajectories where bar assigned to the first population to go extinct in the FULL dataset
trajectoriesAnt = list() #list of trajectories where ant assigned to the first population to go extinct in the FULL dataset

for(i in 1:nrow(FULLkey)){
	trj = FULLkey[i,'trajectory']
	bar = FULLkey[i,'bar']
	ant = FULLkey[i,'ant']

	trajectoriesBar[[trj]] = c(bar,trajectories[[trj]])
	trajectoriesAnt[[trj]] = c(ant,trajectories[[trj]])
	
}


###################################
#	Number of Inferred Species vs. Number of Extinction Events
#	TRAJECTORIES
png('cranioleuca/Fig.S4.png',width=8,height=6,units='in',res=300,bg='transparent')
#dev.new(width=8,height=6)
layout(matrix(c(1,2), 2, 1, byrow = TRUE),
widths=c(7,7), heights=c(3,3))

xlim = c(-0.5,7.5); ylim = c(-50,110)
	if(min(mclust$deltaBIC) < ylim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	if(max(mclust$deltaBIC) > ylim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	if(min(bfd$BF,na.rm=T) < ylim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	if(max(bfd$BF,na.rm=T) > ylim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')


#bfd
	par(mar=c(3.5,5,0,1))
	plot(bfd[ ,'BF'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	axis(1,at=c(0:7),labels=c(0:7),tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)
	# mtext('Number of Extinctions',1,line=2.5,cex=1.25,col='grey25')

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext("Bayes Factor Support",2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Species',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-0.5,105,'A. Genetic Species Delimitation',adj=c(0,0.5),col='grey25')
	
	
	lty = c(1,3,5)
	cols = brewer.pal(6, 'Set1')
	lineKey = expand.grid(lty,cols)

	
	trajectories = trajectoriesBar
	#trajectories = trajectories[!grepl('Sym',names(trajectories))]
	trajectories = trajectories[grepl('Sym3',names(trajectories))]
	#trajectories = trajectoriesAnt
	for(i in 1:length(trajectories)){
		trj = trajectories[[i]]
		tmp = bfd[bfd$scenario %in% trj, ]
		col = lineKey[i,2]	
		lty = lineKey[i,1]
		lines(tmp$nExtinctions,tmp$BF,col=col,lty=lty,lwd=2)
	}
	
		
	
	#shade points based on where in the transect the midpoint of the extinction event occurs
	bfd$col = mapvalues(bfd$extMidpoint,key[,'extMidpoint'], key[,'col'])
	col = add.alpha('grey25',.75) 
	#col = bfd$col
	#col = rep('green',nrow(bfd))
	bg = add.alpha(bfd$col,1); #bg[grep('FULL',bfd$scenario)] = 'transparent'
	pch = 21 	
	lty = 3
	x = jitter(bfd$nExtinctions,0)
	y = jitter(bfd$BF,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=.5)


#Add legend for shading
	key$x = seq(5,7,length.out=nrow(key))
	key$y = 85
	x = key$x
	y = key$y
	col = 'grey25'
	points(x,y,pch=15,col=key$col,cex=1.5)
	text(mean(x),y+22.5,'Extinction Midpoint',cex=.75,col=col)
	latLabels = round(key[key$extMidpoint %in% c(2,3,4,5,6,7,8),'lat'],1)
	latLabels = c(paste0(format(latLabels, nsmall = 1),'°'),'no ext.')
	latX = key[key$extMidpoint %in% c(2,3,4,5,6,7,8,NA),'x']
	text(latX,y+6,latLabels,cex=.5,srt=45,adj=c(0,0.5),col=col)
	text(min(x),y-7.5,'North',cex=.5,col=col)
	text(key[key$extMidpoint %in% max(key$extMidpoint,na.rm=T),'x'],y-7.5,'South',cex=.5,col=col)
	meanX = mean(c(min(x),key[key$extMidpoint %in% max(key$extMidpoint,na.rm=T),'x']))
	text(meanX,y-7.5,'<—————————————>',cex=.5,col=col)
	
#mclust
	par(mar=c(4,5,0,1))
	plot(mclust[ ,'deltaBIC'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	axis(1,at=c(0:7),labels=c(0:7),tick=T,lty=NULL,col.axis = 'black', col.ticks = 1,cex=.75)
	mtext('Number of Extinctions',1,line=2.5,cex=1.25,col='black')

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	#axis(4,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext(expression(paste(Delta,"BIC")),2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Cluster',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Clusters',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
	text(-0.5,105,'B. Phenotypic Species Delimitation',adj=c(0,0.5),col='grey25')
	
	#add dashed line connecting scenarios in Figure 2
	xy = mclust[mclust$scenario  %in% c('FULL','EXT_v4','EXT_v16','EXT_v24'),c('nExtinctions','deltaBIC')]
	xy = xy[order(xy[,1]), ]
	lines(xy[,1],xy[,2],col='grey50')

	mclust$col = mapvalues(mclust$extMidpoint,key[,'extMidpoint'], key[,'col'])
	col = add.alpha('grey25',.75) 
	#col = bfd$col
	#col = rep('green',nrow(bfd))
	bg = mclust$col
	pch = 24 	
	lty = 3
	x = jitter(mclust$nExtinctions,0)
	y = jitter(mclust$deltaBIC,0)
	points(x,y,col=col,bg=bg,pch=pch,cex=2,lwd=.5)

dev.off()







