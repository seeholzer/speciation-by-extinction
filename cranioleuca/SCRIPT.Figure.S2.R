##################################################
#	Description: Figure S2 - sampling maps for three variants of the 
#	simulated extinction scenarios for Cranioleuca antisiensis
#		- produces maps of mclust cluster assignments to each population for all extinction scenarios and versions
#		- these are combined into a single png file in adobe illustrator to make Figure S2
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#load packages
library(maptools)
library(raster)
library(mapplots)
library(RColorBrewer)
source('~/Dropbox/FUN.add.alpha.R', chdir = TRUE)
source('~/Dropbox/FUN.add.pie.new.R', chdir = TRUE)


# library(colorspace)
# source('~/Dropbox/FUN.add.alpha.R', chdir = TRUE)

results.list = get(load('mclust/results.list.rda'))


ex = c(-81.5, -75, -12.5, -2) #raster extent for Cranioleuca antisiensis/baroni
PE0 = readShapeSpatial('map_layers/PER_adm_shp/PER_adm0.shp')
PE1 = readShapeSpatial('map_layers/PER_adm_shp/PER_adm1.shp')
EC0 = readShapeSpatial('map_layers/ECU_adm_shp/ECU_adm0.shp')
EC1 = readShapeSpatial('map_layers/ECU_adm_shp/ECU_adm1.shp')
alt = raster('map_layers/alt_CranCrop.grd')
alt = crop(alt,extent(ex))
range = readShapeSpatial('map_layers/range.map.shp')

pops.no.plumage.data = c('Chinchan','Huamatanga','Pacar')

colkey = data.frame(class=1:3,col=brewer.pal(3,'Set1'),stringsAsFactors=F)

dir = 'mclust.sampling.maps'
if(!dir.exists(dir)){
		dir.create(dir)
	}else if(dir.exists(dir)){
		unlink(dir,recursive=T);dir.create(dir)
	}

i = 13
for(i in 1:length(results.list)){
	tmp = results.list[[i]]
	scenario.name = tmp$scenario.name
	trait.set = tmp$trait.set	
	p = tmp$p
	
	w <- ncol(alt)/max(dim(alt))
	h <- nrow(alt)/max(dim(alt))

	png(paste0('mclust.sampling.maps/',trait.set,'.', scenario.name,'_map.png'),		width=5*w,height=5*h,units='in',res=300,bg='transparent')
	#dev.new(width=5*w, height=5*h)

	plot.new()
	par(mar=c(2,2,0,0), oma=c(0,0,0,0),new=T)
	plot.window(xlim=extent(alt)[1:2], ylim=extent(alt)[3:4], xaxs="i",yaxs="i")
	
	breakpoints <- c(minValue(alt),1000,3500,maxValue(alt))
	colors <- add.alpha(c("grey100","grey50","grey0"),.5)
	plot(alt,add=T,breaks=breakpoints,col=colors,useRaster=T,
		axes=F,box=F,legend=F,interpolate=F,ylim=c(-12.5,0),xlim=c(-81.5,-72))
	axis(1,at=c(-81,-80,-79,-78,-77,-76),labels=NA,line=.1,col.ticks='black',col='transparent',cex.axis=.75)
	mtext(c(-81,-80,-79,-78,-77,-76),1,at=c(-81,-80,-79,-78,-77,-76),line=.35,cex=.75)
	mtext(expression(paste(degree,"Latitude")),1,line=1,cex=.75)
	
	axis(2,at=c(-12,-10,-8,-6,-4),labels=NA,line=0,col.ticks='black', col='transparent',cex.axis=.75)
	mtext(c(-12,-10,-8,-6,-4),2,at=c(-12,-10,-8,-6,-4),line=.4,cex=.75)
	mtext(expression(paste(degree,"Longitude")),2,line=1,cex=.75)
	
	admin.lines = 'grey50' 
	plot(PE0,add=T,lwd=0.75,border=admin.lines)
	plot(EC0,add=T,lwd=0.75,border=admin.lines)
	plot(PE1,add=T,lwd=0.1,border=admin.lines)
	plot(EC1,add=T,lwd=0.1,border=admin.lines)
	plot(range,add=T,col=add.alpha(brewer.pal(9,'Blues')[4],.5),border=add.alpha(brewer.pal(9,'Blues')[9],.5))

	f = 16	
	for(f in 1:nrow(p)){
		scenario.col = colnames(p)[apply(p,2,function(x)!all(!(grepl('extinct|ant|bar',x))))]
		cat(f,'\n')
		if(p[f,'pop'] %in% pops.no.plumage.data & !grepl('morphology', trait.set)){
			pie.columns = grep('N\\d',colnames(p))
			z = as.matrix(p[f,pie.columns])
			z[1,1] = 1; z[is.na(z)] = 0 
			x=p[f,'long']; y=p[f,'lat']
			add.pie.new(z=z,x=x,y=y,col='transparent', radius=0.15,labels=NA,edge=600,lwd=.25)			
		}else if(p[f,scenario.col] %in% 'extinct'){
			pie.columns = grep('N\\d',colnames(p))
			z = matrix(nrow=1,ncol=length(pie.columns))
			z[1,1] = 1; z[is.na(z)] = 0
			x=p[f,'long']
			y=p[f,'lat']
			add.pie.new(z=z,x=x,y=y,col=add.alpha('grey75',1),radius=0.15,labels=NA,edge=600,lwd=.5)
		
		}else if(p[f,scenario.col] %in% c('ant','bar')){
			pie.columns = grep('N\\d',colnames(p))
			z = as.matrix(p[f,pie.columns])
			x=p[f,'long']
			y=p[f,'lat']
			add.pie.new(z=z,x=x,y=y,col=colkey[1:2,'col'], radius=0.15,labels=NA,edge=600,lwd=.25)
			}
		}

	text(-76.5,-3.25, trait.set,cex=1.5)
	text(-76.5,-4, scenario.name,cex=1.5)
	
	dev.off()

}

