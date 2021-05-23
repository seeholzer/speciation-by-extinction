##################################################
#	Description: Create Figure S1
#		sampling maps for the three variants of the 	simulated extinction scenarios 
#		for Cranioleuca antisiensisSpecies delimitation using Normal Mixture Models 
#		for scenarios of simulated extinction in Cranioleuca antisiensis.
#
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')


#load packages
library(maptools)
library(raster)
library(plyr)
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

colfunc = colorRampPalette(c("#425BD3", "#E8B618"))

color.key = data.frame(scenario=c('C.ant.SBE3','C.ant.SBE2','C.ant.SBE1','FULL','C.bar.SBE1','C.bar.SBE2','C.bar.SBE3'),col=colfunc(7))

scenario.key = cbind(c('FULL','SBE1','SBE2','SBE3'),c('Pre-Extinction','Post-Extinction 1','Post-Extinction 2','Post-Extinction 3'))

#plot sampling map for each geographic variant of the simulated extinction scenarios SBE1-SBE3
#Figure S1 = v1, v2, v3
v = 'v1'
for(v in c('v1','v2','v3')){

	popdata = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
	popdata = popdata[!popdata$pop %in% 'curtata',c(1:3,grep(v,colnames(popdata)))]

	popdata$scenario = 'FULL'
	popdata$scenario[popdata$SBE3 %in% 'extinct'] = 'SBE3'
	popdata$scenario[popdata$SBE2 %in% 'extinct'] = 'SBE2'
	popdata$scenario[popdata$SBE1 %in% 'extinct'] = 'SBE1'
	
	w <- ncol(alt)/max(dim(alt))
	h <- nrow(alt)/max(dim(alt))
	## Set up appropriately sized device with no borders and required coordinate system    
	## png("eg.png", width=480*w, height=480*h)
#	if(v == 'v2'){filename = 'Fig.2.png'}else{filename = paste0('Fig.S1.',v,'.png')}
	filename = paste0('Fig.S1.',v,'.png')
	png(filename,width=5*w,height=5*h,units='in',res=300,bg='transparent')
	#set up blank plot space, necessary to avoid white blank space around raster plot
	#dev.new(width= 5*w, height= 5*h)
	plot.new()
	par(mar=c(2,2,0,0), oma=c(0,0,0,0),new=T)
	plot.window(xlim=extent(alt)[1:2], ylim=extent(alt)[3:4], xaxs="i",yaxs="i")
	
	breakpoints <- c(minValue(alt),1000,3500,maxValue(alt))
	colors <- add.alpha(c("grey100","grey60","grey20"),.5)
	plot(alt,add=T,breaks=breakpoints,col=colors,useRaster=T,
		axes=F,box=F,legend=F,interpolate=F,ylim=c(-12.5,0),xlim=c(-81.5,-72))
	axis(1,at=c(-81,-80,-79,-78,-77,-76),labels=NA,line=.1,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-81,-80,-79,-78,-77,-76),1,at=c(-81,-80,-79,-78,-77,-76),line=.35,cex=.75)
	mtext(expression(paste(degree,"Latitude")),1,line=1,cex=.75)
	
	axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
	mtext(expression(paste(degree,"Longitude")),2,line=1,cex=.75)
	
	admin.lines = 'grey50' 
	plot(PE0,add=T,lwd=0.75,border=admin.lines)
	plot(EC0,add=T,lwd=0.75,border=admin.lines)
	plot(PE1,add=T,lwd=0.1,border=admin.lines)
	plot(EC1,add=T,lwd=0.1,border=admin.lines)
	plot(C.ant.FULL,add=T,col=add.alpha('#9ECAE1',.5),border=add.alpha('#9ECAE1',.5))
	
	key = data.frame(scenario=c('FULL','SBE3','SBE2','SBE1'),
						new.name=c('Pre-Extinction','Post-Extinction 3','Post-Extinction 2','Post-Extinction 1'),
						col=c("#FFEDA0","#E6AB02","#A6761D","#543005"),
						text.col = c('black','black','white','white'),
						x=-81.25,
						y=c(-9,-11.25,-10.5,-9.75)+0.5)
	
	bg = mapvalues(popdata$scenario,key[,'scenario'],key[,'col'])
	text.col = mapvalues(popdata$scenario,key[,'scenario'],key[,'text.col'])
	for(i in 1:nrow(popdata)){
		points(popdata[i,'long'],popdata[i,'lat'],pch=21,cex=2,bg=bg[i],col='grey50',lwd=1)
		text(popdata[i,'long'],popdata[i,'lat'],c(1:nrow(popdata))[i],cex=.5,col=text.col[i])
	}
	
	if(v == 'v1'){
		x = -77.5
		y = -6.25
		text(x,y,'Origin of',cex=.75,adj=c(0,.5),col='#543005')
		text(x,y-.25,'Local Extinction',cex=.75,adj=c(0,.5),col='#543005')
	}

	if(v == 'v2'){
		x = -77.4
		y = -7.5
		text(x,y,'Origin of',cex=.75,adj=c(0,.5),col='#543005')
		text(x,y-.25,'Local Extinction',cex=.75,adj=c(0,.5),col='#543005')
	}
	
	if(v == 'v3'){
		x = -77.5
		y = -8.5
		text(x,y,'Origin of',cex=.75,adj=c(0,.5),col='#543005')
		text(x,y-.25,'Local Extinction',cex=.75,adj=c(0,.5),col='#543005')
	}	
	
	#populations removed in each scenario
	text(key$x,key$y,key[,'new.name'],cex=.75,adj=c(0,0))
	offset = .25
	
	x = key[key$scenario %in% c('FULL'),]
	text(x$x+offset,x$y-0.25,expression('\U2012'),c(.5,1))
	
	x = key[key$scenario %in% c('SBE1'),]
	points(unique(x$x)+seq(0,0,0.5)+offset,rep(x[x$scenario %in% 'SBE1','y'],1)-0.25,pch=21,bg=rev(x$col),col='grey50',cex=1.75)
	
	x = key[key$scenario %in% c('SBE1','SBE2'),]
	points(unique(x$x)+seq(0,0.5,0.5)+offset,rep(x[x$scenario %in% 'SBE2','y'],2)-0.25,pch=21,bg=rev(x$col),col='grey50',cex=1.75)
	
	x = key[key$scenario %in% c('SBE1','SBE2','SBE3'),]
	points(unique(x$x)+seq(0,1,0.5)+offset,rep(x[x$scenario %in% 'SBE3','y'],3)-0.25,pch=21,bg=rev(x$col),col='grey50',cex=1.75)
	
	x = key[key$scenario %in% c('FULL'),]
	text(x$x+offset+0.5,x$y+.6+.2,expression(italic('in silico')),adj=c(.25,.5),cex=.75)
	text(x$x+offset+0.5,x$y+.35+.2,bquote(underline('extinction')),adj=c(.25,.5),cex=.75)
	
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
			
	dev.off()

}