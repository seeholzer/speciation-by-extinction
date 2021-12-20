##################################################
#	SCRIPT.Fig.2.maps
#	Description: Create four maps for panel A of Figure 2 
#
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

#colfunc = colorRampPalette(c("#425BD3", "#E8B618"))
colfunc = colorRampPalette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))
# colfunc(7)
# plot(rep(1,7),col=colfunc(7),pch=19,cex=3)

#create a color key for the scenarios
color.key = data.frame(scenario=c("C.ant.SBE3","C.ant.SBE2","C.ant.SBE1","FULL","C.bar.SBE1","C.bar.SBE2","C.bar.SBE3"),col=colfunc(7))
color.key$col = add.alpha(color.key$col,0.75)

#population data
popdata = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
#subset to a sequential, nested series of extinction scenarios of increasing magnitude 
popdata = popdata[!popdata$pop %in% 'curtata',c(1:5,grep(paste(c("FULL_v5","EXT_v4","EXT_v16","EXT_v24"),collapse='|'),colnames(popdata)))]

#modify popdata so that aligns with the labelling of shapefile C.ant.FULL 
popdata$scenario = 'FULL'
popdata$scenario[popdata$EXT_v24 %in% 'extinct'] = 'SBE3'
popdata$scenario[popdata$EXT_v16 %in% 'extinct'] = 'SBE2'
popdata$scenario[popdata$EXT_v4 %in% 'extinct'] = 'SBE1'

#change names of scenarios
scenario.key = cbind(c('FULL','SBE1','SBE2','SBE3'),c('Pre-Extinction 5','Post-Extinction 4','Post-Extinction 16','Post-Extinction 24'))


#scaling factors of height and width 
w <- ncol(alt)/max(dim(alt))
h <- nrow(alt)/max(dim(alt))

## Set up appropriately sized device with no borders and required coordinate system    
## png("eg.png", width=480*w, height=480*h)
i = 'FULL'
for(i in c('FULL','SBE1','SBE2','SBE3')){

	scenario = scenario.key[scenario.key[,1] %in% i,2]

	filename = paste0('Fig.2a.',scenario,'.png') 
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

	if(i == 'FULL'){
		axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black',col='transparent',cex.axis=.75)
		mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
		mtext(expression(paste(degree,"Longitude")),2,line=1,cex=.75)
	}
	
	admin.lines = 'grey50' 
	plot(PE0,add=T,lwd=0.75,border=admin.lines)
	plot(EC0,add=T,lwd=0.75,border=admin.lines)
	plot(PE1,add=T,lwd=0.1,border=admin.lines)
	plot(EC1,add=T,lwd=0.1,border=admin.lines)

	if(i == 'FULL'){
	
		col = mapvalues(as.character(range@data$scenario),color.key[,1],color.key[,2])
	
		plot(range,add=T,col=col,border='transparent',lwd=0.001)
		extinct.pops.removed = range
		border = unionSpatialPolygons(extinct.pops.removed,rep(1,nrow(extinct.pops.removed@data)))
		plot(border,add=T,col='transparent',border='black',lwd=.5)
				
		bg = '#FFEDA0'
		text.col = 'black'

		points(popdata[ ,'long'],popdata[ ,'lat'],pch=21,cex=2,bg=bg,col='grey50',lwd=1)
		text(popdata[ ,'long'],popdata[ ,'lat'],1:nrow(popdata),cex=.5,col=text.col)
		
	}else if(i == 'SBE1'){
		
		tmp.color.key = color.key
		tmp.color.key[grepl(paste(c('FULL'),collapse='|'),tmp.color.key$scenario),'col'] = 'transparent'
		col = mapvalues(as.character(range@data$scenario), tmp.color.key[,1], tmp.color.key[,2])

		plot(range,add=T,col=col,border='transparent',lwd=0.001)
		
		extinct.pops.removed = range[!grepl(paste(c('FULL'),collapse='|'),range@data$scenario),]
		border = unionSpatialPolygons(extinct.pops.removed,rep(1,nrow(extinct.pops.removed@data)))
		plot(border,add=T,col='transparent',border='black',lwd=.5)
		
		bg = rep('#FFEDA0',nrow(popdata))
		bg[popdata$EXT_v4 %in% 'extinct'] = 'grey90'
		text.col = 'black'
			
		points(popdata[ ,'long'],popdata[ ,'lat'],pch=21,cex=2,bg=bg,col='grey50',lwd=1)
		text(popdata[ ,'long'],popdata[ ,'lat'],1:nrow(popdata),cex=.5,col=text.col)		
	
	}else if(i == 'SBE2'){
		
		tmp.color.key = color.key
		tmp.color.key[grepl(paste(c('FULL','SBE1'),collapse='|'),tmp.color.key$scenario),'col'] = 'transparent'
		col = mapvalues(as.character(range@data$scenario), tmp.color.key[,1], tmp.color.key[,2])

		plot(range,add=T,col=col,border='transparent',lwd=0.001)
		
		extinct.pops.removed = range[!grepl(paste(c('FULL','SBE1'),collapse='|'),range@data$scenario),]
		border = unionSpatialPolygons(extinct.pops.removed,rep(1,nrow(extinct.pops.removed@data)))
		plot(border,add=T,col='transparent',border='black',lwd=.5)

		bg = rep('#FFEDA0',nrow(popdata))
		bg[popdata$EXT_v16 %in% 'extinct'] = 'grey90'
		text.col = 'black'
			
		points(popdata[ ,'long'],popdata[ ,'lat'],pch=21,cex=2,bg=bg,col='grey50',lwd=1)
		text(popdata[ ,'long'],popdata[ ,'lat'],1:nrow(popdata),cex=.5,col=text.col)		
		
	}else if(i == 'SBE3'){
		
		tmp.color.key = color.key
		tmp.color.key[grepl(paste(c('FULL','SBE1','SBE2'),collapse='|'),tmp.color.key$scenario),'col'] = 'transparent'
		col = mapvalues(as.character(range@data$scenario), tmp.color.key[,1], tmp.color.key[,2])

		plot(range,add=T,col=col,border='transparent',lwd=0.001)
		
		extinct.pops.removed = range[!grepl(paste(c('FULL','SBE1','SBE2'),collapse='|'),range@data$scenario),]
		border = unionSpatialPolygons(extinct.pops.removed,rep(1,nrow(extinct.pops.removed@data)))
		plot(border,add=T,col='transparent',border='black',lwd=.5)

		bg = rep('#FFEDA0',nrow(popdata))
		bg[popdata$EXT_v24 %in% 'extinct'] = 'grey90'
		text.col = 'black'
			
		points(popdata[ ,'long'],popdata[ ,'lat'],pch=21,cex=2,bg=bg,col='grey50',lwd=1)
		text(popdata[ ,'long'],popdata[ ,'lat'],1:nrow(popdata),cex=.5,col=text.col)		
		
	}
	
	
	if(i == 'FULL'){
		
		#dividing line between antisiensis and baroni
		x0 = -76.5; y0 = -7.7 
		x1 = -79.4; 	y1 = -8.9		
		segments(x0,y0,x1,y1,col='grey25',lty=2)
		text(x0,y0,'< antisiensis  ',cex=.75,col='#d7191c',srt=295,adj=c(1,0))
		text(x0,y0,'  baroni >',cex=.75,col='#2c7bb6',srt=295,adj=c(0,0))
		
		#scale bar
		scalebar(d = 200, xy = c(extent(alt)[1]+.5,extent(alt)[3]+0.1),type = "line", lonlat = TRUE,lwd = 2,label='')
		text(extent(alt)[1]+.5,extent(alt)[3]+0.1,'200 km',adj=c(-1,-0.75),cex=.5)

		#elevation legend
		xo = -81.25
		yo = -12
		x=xo; 		rect(x,yo,x+0.75,yo+0.15,col=colors[1])
		x=xo+0.75; 	rect(x,yo,x+0.75,yo+0.15,col=colors[2])
		x=xo+1.5; 	rect(x,yo,x+0.75,yo+0.15,col=colors[3])
		text(xo,yo+0.2,'0',adj=c(.5,0),cex=.5)
		text(xo+0.75,yo+0.2,'1000',adj=c(.5,0),cex=.5)
		text(xo+1.5,yo+0.2,'3500',adj=c(.5,0),cex=.5)
		text(xo+1.5+1,yo+0.2,'3500+ m',adj=c(.5,0),cex=.5)
	
		}
	
	#extinct/extant legend + scenario label
	x = -79.9
	y = -10.6
	offset = 0.5

	points(x - offset,y,pch=21,cex=2,bg='#FFEDA0',col='grey50',lwd=1)
	text(x,y,'extant',cex=.75,col='black',adj=c(.25,.5))

	points(x - offset,y-0.5,pch=21,cex=2,bg='grey90',col='grey50',lwd=1)
	text(x,y-0.5,'extinct',cex=.75,col='black',adj=c(.25,.5))
	
	#scenario label
	text(x,y+0.5,scenario,cex=1)
	# a = strsplit(foo,'-')[[1]][1]
	# b = strsplit(foo,'-')[[1]][2]
	# text(x,y+0.5,bquote(paste(bold(.(a)),'-',.(b))),cex=1)
	
	if(i == 'SBE3'){
	
	par(fig=c(.75,1,.75,1),mar = rep(0, 4),new=T)
		xlim=c(-82,-68)
		ylim=c(-16,2)
		data(wrld_simpl)
	#	rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col='white')
		plot(wrld_simpl[wrld_simpl$NAME %in% c('Peru','Ecuador','Colombia'), ],col='transparent',border='transparent',xlim=xlim,ylim=ylim)
		rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
		plot(wrld_simpl[wrld_simpl$NAME %in% c('Peru','Ecuador','Colombia','Bolivia','Chile'), ],col='white',border='black',xlim=xlim,ylim=ylim,add=T)
		rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4])
		rect(extent(alt)[1],extent(alt)[3],extent(alt)[2],extent(alt)[4],lty=2,col=add.alpha('grey50',.5))
		text(-73,3,labels='Colombia',cex=.5)
		text(-78,-1,labels='Ecuador',cex=.4)
		text(-73,-13,labels='Peru',cex=.5)
	}
			
	dev.off()

}
