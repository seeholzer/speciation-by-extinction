##################################################
#	Description: 
#	1. Figure S1a - map of population groups
#	2. Figure S1b - schematic of extinction scenarios  
#	3. create maps of individual scenarios (not in paper)
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')


#load packages
library(maptools)
library(raster)
library(plyr)
library(RColorBrewer)
source('FUN.add.alpha.R', chdir = TRUE)


#load layers
ex = c(-81.5, -75, -12.5, -2) #raster extent for Cranioleuca antisiensis/baroni
PE0 = readShapeSpatial('map_layers/PER_adm_shp/PER_adm0.shp')
PE1 = readShapeSpatial('map_layers/PER_adm_shp/PER_adm1.shp')
EC0 = readShapeSpatial('map_layers/ECU_adm_shp/ECU_adm0.shp')
EC1 = readShapeSpatial('map_layers/ECU_adm_shp/ECU_adm1.shp')
alt = raster('map_layers/alt_CranCrop.grd')
alt = crop(alt,extent(ex))

range = readShapeSpatial('map_layers/C.ant.scenarios.shp')
C.ant.FULL  = readShapeSpatial('map_layers/C.ant.FULL.shp')

colfunc = colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
#dev.new(width=8,height=3); plot(rep(1,9),col=colfunc(9),pch=19,cex=5)


#population data
popdata = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]

scenarios = colnames(popdata)[!colnames(popdata) %in% c('popNumber','pop','long','lat','group')]


#####################################################################################
#####################################################################################
#####################################################################################
#	1. Figure S1a - map displaying population groups
#####################################################################################
#####################################################################################
#####################################################################################

w <- ncol(alt)/max(dim(alt))
h <- nrow(alt)/max(dim(alt))

filename = 'Fig.S1a.png'
png(filename,width=5*w,height=5*h,units='in',res=300,bg='transparent')

	#set up blank plot space, necessary to avoid white blank space around raster plot
	#dev.new(width= 5*w, height= 5*h)
	plot.new()
	par(mar=c(2,2,0,0), oma=c(0,0,0,0),new=T)
	plot.window(xlim=extent(alt)[1:2], ylim=extent(alt)[3:4], xaxs="i",yaxs="i")
	
	breakpoints <- c(minValue(alt),1000,3500,maxValue(alt))
	#colors <- add.alpha(c("grey100","grey50","grey0"),.5)
	colors <- add.alpha(c("grey100","grey60","grey20"),.5)
	plot(alt,add=T,breaks=breakpoints,col=colors,useRaster=T,
		axes=F,box=F,legend=F,interpolate=F,ylim=c(-12.5,0),xlim=c(-81.5,-72))

	axis(1,at=c(-81,-80,-79,-78,-77,-76),labels=NA,line=.1,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-81,-80,-79,-78,-77,-76),1,at=c(-81,-80,-79,-78,-77,-76),line=.35,cex=.75)
	mtext(expression(paste(degree,"Longitude")),1,line=1,cex=.75)

	axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
	mtext(expression(paste(degree,"Latitude")),2,line=1,cex=.75)

	#Figure S1a label
	text(-81.25,-2.35,'A.',cex=1.25)
	
	admin.lines = 'grey50' 
	plot(PE0,add=T,lwd=0.75,border=admin.lines)
	plot(EC0,add=T,lwd=0.75,border=admin.lines)
	plot(PE1,add=T,lwd=0.1,border=admin.lines)
	plot(EC1,add=T,lwd=0.1,border=admin.lines)
	
	plot(C.ant.FULL,col=add.alpha('#d01c8b',.1),border='grey25',add=T,lwd=0.75)

	tmp = popdata[,c('popNumber','pop','long','lat','group')]
	tmp$col = mapvalues(tmp$group,1:9,colfunc(9))
	
	#Group 1
	text(-80.6,-5,'G1',cex=.75)
	pop=1;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -80.6; 	y1 = -4.8
		segments(x0,y0,x1,y1,col='black')
	pop=2;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -80.3; y1 = -5
		segments(x0,y0,x1,y1,col='black')
	#Group 2
	text(-80.5,-5.8,'G2',cex=.75)
	pop=3;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -80.2; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 3
	text(-80.7,-6.5,'G3',cex=.75)
	pop=4;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -80.4; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 4
	text(-80.1,-7.4,'G4',cex=.75)
	group=4;foo = tmp[tmp$group == group, ]
		foo = foo[match(foo$popNumber,c(6,5,7)),]
		x=foo[,'long'];y=foo[,'lat']
		lines(x,y)
	pop=6;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -79.8; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 5
	text(-79.8,-7.77,'G5',cex=.75)
	pop=8;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -79.5; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 6
	text(-79.3,-8.75,'G6',cex=.75)
	group=6;foo = tmp[tmp$group == group, ]
		x=foo[,'long'];y=foo[,'lat']
		lines(x,y)
	pop=10;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -79; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 7
	text(-78.7,-10.075,'G7',cex=.75)
	group=7;foo = tmp[tmp$group == group, ]
		foo = foo[match(foo$popNumber,c(13,12,14)),]
		x=foo[,'long'];y=foo[,'lat']
		lines(x,y)
	pop=14;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -78.4; 	y1 = y0		
		segments(x0,y0,x1,y1,col='black')
	#Group 8
	text(-78.4,-10.581657,'G8',cex=.75)
	pop=15;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -77.49522; 	y1 = -10.581657	
		segments(x0,y0,x1,y1,col='black')
		x0 = -77.49522; 	y0 = -10.581657	
		x1 = -78.1; 	y1 = y0	
		segments(x0,y0,x1,y1,col='black')
	#Group 9
	text(-77.9,-11.51,'G9',cex=.75)
	group=9;foo = tmp[tmp$group == group, ]
		foo = foo[match(foo$popNumber,c(16,17,19,18)),]
		foo = rbind(foo,foo[foo$popNumber == 16, ])
		x=foo[,'long'];y=foo[,'lat']
		lines(x,y)
	pop=19;x0=tmp[tmp$popNumber == pop,'long'];y0=tmp[tmp$popNumber == pop,'lat']
		x1 = -77.6; 	y1 = y0	
		segments(x0,y0,x1,y1,col='black')
	
	#add population points
	points(tmp[ ,'long'], tmp[ ,'lat'],pch=21,cex=2,bg=tmp$col,col='grey50',lwd=1)
	text(tmp[ ,'long'], tmp[ ,'lat'],1:nrow(tmp),cex=.5,col='black')		
		

#elevation legend
	scalebar(d = 200, xy = c(extent(alt)[1]+.5,extent(alt)[3]+0.1),type = "line", lonlat = TRUE,lwd = 2,label='')
	text(extent(alt)[1]+.5,extent(alt)[3]+0.1,'200 km',adj=c(-1,-0.75),cex=.5)

	xo = -81.25
	yo = -12
	x=xo; 		rect(x,yo,x+0.75,yo+0.15,col=colors[1])
	x=xo+0.75; 	rect(x,yo,x+0.75,yo+0.15,col=colors[2])
	x=xo+1.5; 	rect(x,yo,x+0.75,yo+0.15,col=colors[3])
	text(xo,yo+0.2,'0',adj=c(.5,0),cex=.5)
	text(xo+0.75,yo+0.2,'1000',adj=c(.5,0),cex=.5)
	text(xo+1.5,yo+0.2,'3500',adj=c(.5,0),cex=.5)
	text(xo+1.5+1,yo+0.2,'3500+ m',adj=c(.5,0),cex=.5)	

dev.off()




#####################################################################################
#####################################################################################
#####################################################################################
#	2. Figure S1b - schematic of extinction scenarios
#####################################################################################
#####################################################################################
#####################################################################################
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#population data
popdata = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]
scenarios = colnames(popdata)[!colnames(popdata) %in% c('popNumber','pop','long','lat','group')]
colfunc = colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))


filename = 'Fig.S1b.png'
png(filename,width=3,height=5,units='in',res=300,bg='transparent')

#dev.new(width=3,height=5)
par(mar=c(0,0,0,0))
xlim = c(0,46)
ylim = c(9.51,-0.4)
plot(1, xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)

#Figure S1b label
text(1,-.6,'B.',cex=1.25)

#Group
text(41,-.25,'Population',cex=.75)
text(41,0.025,'Group',cex=.75)

bg = colfunc(9)
col = 'grey25'
x = rep(40,9)
y = 1:9
points(x,y,pch=21,bg=bg,col=col,cex=2.5)
text(x+4,y,paste0('G',1:9))

abline(h=seq(1.5,8.5,1),lty=3,col='grey25',lwd=.75)

#Pre-Extinction
text(4.5,-.25,'Pre-',cex=.75)
text(4.5,0,'Extinction',cex=.75)

i = 1
for(i in grep('FULL',scenarios)){
	scn = scenarios[i]
	
	foo = unique(popdata[,c('group',scn)])
	ant = max(foo[foo[,2] == 'ant','group'])
	bar = min(foo[foo[,2] == 'bar','group'])
	transition = mean(c(ant,bar))
	lwd = 1
	#ant
	col = '#d7191c'
	x0 = i; y0 = 0.5
	x1 = i; y1 = transition	
	segments(x0,y0,x1,y1,lty=3,lwd=lwd,col=col)
	points(x0,y0,pch=16,col=col,cex=.5)
	#bar
	col = '#2c7bb6'
	x0 = i; y0 = transition
	x1 = i; 	y1 = 9.5	
	segments(x0,y0,x1,y1,lty=3,lwd=lwd,col=col)	
	points(x1,y1,pch=16,col=col,cex=.5)
	#transition
	x = i; y = transition
	points(x,y,pch=15,cex=.75)	
	#label
	lab = gsub('FULL_v','',scn)
	text(x,y-.2,lab,cex=.5)	
}


#Post-Extinction
text(mean(c(9,length(scenarios)))+1,-.15,'Post-Extinction',cex=.75)
i = 16
for(i in grep('EXT',scenarios)){
	scn = scenarios[i]
	
	#Post-Extinction
	foo = unique(popdata[,c('group',scn)])
	ant = max(foo[foo[,2] == 'ant','group']) + 0.5
	bar = min(foo[foo[,2] == 'bar','group']) - 0.5
	
	lwd = 1
	#ant
	col = '#d7191c'
	x0 = i+1; y0 = 0.5
	x1 = i+1; y1 = ant	
	segments(x0,y0,x1,y1,lty=3,lwd=lwd,col=col)
	points(x0,y0,pch=16,col=col,cex=.5)

	#bar
	col = '#2c7bb6'
	x0 = i+1; y0 = bar
	x1 = i+1; 	y1 = 9.5	
	segments(x0,y0,x1,y1,lty=3,lwd=lwd,col=col)	
	points(x1,y1,pch=16,col=col,cex=.5)

	#Extinction gap
	gap = foo[foo[,2] == 'extinct','group']
	min = min(gap) - 0.5
	max = max(gap) + 0.5
	x0 = i+1; y0 = min
	x1 = i+1; y1 = max	
	lwd = 2.5
	points(x0,y0,pch=16,cex=.5)
	points(x1,y1,pch=16,cex=.5)
	segments(x0,y0,x1,y1,lty=1,lwd=lwd)
	
	#label
	lab = gsub('EXT_v','',scn)
	text(x0,y0-.2,lab,cex=.5)
}

#add black border to separate it from Figure S1a
abline(v=-1.6,lwd=1,lty=2)


dev.off()





#####################################################################################
#####################################################################################
#####################################################################################
#	3. Map for each scenario (not in paper)
#####################################################################################
#####################################################################################
#####################################################################################

#scaling factors of height and width 
w <- ncol(alt)/max(dim(alt))
h <- nrow(alt)/max(dim(alt))

## Set up appropriately sized device with no borders and required coordinate system    
## png("eg.png", width=480*w, height=480*h)
i = 'FULL_v8'
for(i in scenarios){
	cat(i,'\n')
	filename = paste0('scenario.maps/scenario.map.',i,'.png') 
	png(filename,width=5*w,height=5*h,units='in',res=300,bg='transparent')
	#set up blank plot space, necessary to avoid white blank space around raster plot
	#dev.new(width= 5*w, height= 5*h)
	plot.new()
	par(mar=c(2,2,0,0), oma=c(0,0,0,0),new=T)
	plot.window(xlim=extent(alt)[1:2], ylim=extent(alt)[3:4], xaxs="i",yaxs="i")
	
	breakpoints <- c(minValue(alt),1000,3500,maxValue(alt))
	#colors <- add.alpha(c("grey100","grey50","grey0"),.5)
	colors <- add.alpha(c("grey100","grey60","grey20"),.5)
	plot(alt,add=T,breaks=breakpoints,col=colors,useRaster=T,
		axes=F,box=F,legend=F,interpolate=F,ylim=c(-12.5,0),xlim=c(-81.5,-72))

	axis(1,at=c(-81,-80,-79,-78,-77,-76),labels=NA,line=.1,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-81,-80,-79,-78,-77,-76),1,at=c(-81,-80,-79,-78,-77,-76),line=.35,cex=.75)
	mtext(expression(paste(degree,"Latitude")),1,line=1,cex=.75)

	axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
	mtext(expression(paste(degree,"Longitude")),2,line=1,cex=.75)

	# if(i == 'FULL'){
		# axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black',col='transparent',cex.axis=.75)
		# mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
		# mtext(expression(paste(degree,"Longitude")),2,line=1,cex=.75)
	# }
	
	admin.lines = 'grey50' 
	plot(PE0,add=T,lwd=0.75,border=admin.lines)
	plot(EC0,add=T,lwd=0.75,border=admin.lines)
	plot(PE1,add=T,lwd=0.1,border=admin.lines)
	plot(EC1,add=T,lwd=0.1,border=admin.lines)

	plot(C.ant.FULL,col=add.alpha('#d01c8b',.1),border='grey25',add=T,lwd=0.75)
	
	ant.color = '#d7191c'
	ext.color = 'grey90'
	bar.color = '#2c7bb6'
	
	bg = rep(NA,nrow(popdata))
	bg[popdata[ ,i] %in% 'ant'] = ant.color
	bg[popdata[ ,i] %in% 'extinct'] = ext.color
	bg[popdata[ ,i] %in% 'bar'] = bar.color
	text.col = 'black'
		
	points(popdata[ ,'long'],popdata[ ,'lat'],pch=21,cex=2,bg=bg,col='grey50',lwd=1)
	text(popdata[ ,'long'],popdata[ ,'lat'],1:nrow(popdata),cex=.5,col=text.col)		

		
#ant/extinct/bar legend + scenario label
	x = -80.2
	y = -10.6
	offset = .5

	points(x - offset,y+0.5,pch=21,cex=2,bg= ant.color,col='grey50',lwd=1)
	text(x,y+0.5,'C. antisiensis',cex=.75,col='black',adj=c(0,.5))

	points(x - offset,y,pch=21,cex=2,bg= ext.color,col='grey50',lwd=1)
	text(x,y,'extinct',cex=.75,col='black',adj=c(0,.5))

	points(x - offset,y-0.5,pch=21,cex=2,bg= bar.color,col='grey50',lwd=1)
	text(x,y-0.5,'C. baroni',cex=.75,col='black',adj=c(0,.5))
	

#scenario label
	text(x,y+1,i,cex=1)
	# a = strsplit(foo,'-')[[1]][1]
	# b = strsplit(foo,'-')[[1]][2]
	# text(x,y+0.5,bquote(paste(bold(.(a)),'-',.(b))),cex=1)

#elevation legend

	scalebar(d = 200, xy = c(extent(alt)[1]+.5,extent(alt)[3]+0.1),type = "line", lonlat = TRUE,lwd = 2,label='')
	text(extent(alt)[1]+.5,extent(alt)[3]+0.1,'200 km',adj=c(-1,-0.75),cex=.5)

	xo = -81.25
	yo = -12
	x=xo; 		rect(x,yo,x+0.75,yo+0.15,col=colors[1])
	x=xo+0.75; 	rect(x,yo,x+0.75,yo+0.15,col=colors[2])
	x=xo+1.5; 	rect(x,yo,x+0.75,yo+0.15,col=colors[3])
	text(xo,yo+0.2,'0',adj=c(.5,0),cex=.5)
	text(xo+0.75,yo+0.2,'1000',adj=c(.5,0),cex=.5)
	text(xo+1.5,yo+0.2,'3500',adj=c(.5,0),cex=.5)
	text(xo+1.5+1,yo+0.2,'3500+ m',adj=c(.5,0),cex=.5)

			
	dev.off()

}
