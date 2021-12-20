##################################################
#	Description: Create maps for Figure 2
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

#colfunc = colorRampPalette(c("#E19DC9", "#C9E19D"))
#colfunc = colorRampPalette(c("#ACAFD3", "#EAD898"))
colfunc = colorRampPalette(c("#425BD3", "#E8B618"))
# colfunc(7)
# plot(rep(1,7),col=colfunc(7),pch=19,cex=3)

color.key = data.frame(scenario=c('C.ant.SBE3','C.ant.SBE2','C.ant.SBE1','FULL','C.bar.SBE1','C.bar.SBE2','C.bar.SBE3'),col=colfunc(7))


#population data
popdata = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
popdata = popdata[!popdata$pop %in% 'curtata', ]

scenarios = colnames(popdata)[!colnames(popdata) %in% c('popNumber','pop','long','lat','group')]


#####################################################################################
#	Map displaying extinction groups
#####################################################################################


#####################################################################################
#	Map for each scenario
#####################################################################################

#scaling factors of height and width 
w <- ncol(alt)/max(dim(alt))
h <- nrow(alt)/max(dim(alt))

## Set up appropriately sized device with no borders and required coordinate system    
## png("eg.png", width=480*w, height=480*h)
i = 'FULL_v8'
for(i in scenarios){

	filename = paste0('Fig.S1.maps/',i,'.png') 
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

	
	ant.color = colfunc(7)[2]
	ext.color = 'grey90'
	bar.color = colfunc(7)[6]
	
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
