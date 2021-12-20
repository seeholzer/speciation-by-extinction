##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	Script 3 - Figure 3
#
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza')

#load packages
library(plyr)
library(mixtools)
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza')
source('FUN.add.alpha.R', chdir = TRUE)

##################################################
# 3.1) Read data
data = read.table("Geospiza.data.csv", header=T, sep=",")
data$ID = paste0(data$Institution,data$Museum.Number)
colnames(data)[1] = 'Taxa.Lack'
colnames(data)[6] = 'Taxa.New'
colnames.to.keep = c('ID','Taxa.Lack','Taxa.New','Wing','Tail','Blength','Bdepth','Bwidth','Tarsus')
data = data[,colnames.to.keep]

#load mclust results
mclust.data = get(load('mclust.Bdepth.SML.groundfinches.rda'))
plots.dir = 'plots.Bdepth.SML.groundfinches'


##################################################
# 3.2) Merge all mclust classifications and other relevant data

#create list of dataframes with the mclust classification of each extant specimen
l = lapply(mclust.data,function(x){
	mclust = x$mclust
	extant = x$extant
	tmp = data.frame(ID = extant, class = mclust$classification)
	return(tmp)	
	})

#rename class column as the extinction scenario
for(i in 1:length(l)) colnames(l[[i]])[2] = names(l)[i]

#choose scenario
scenario = "7.00-8.25"

#extract pre.extinction and post.extinction
pre.extinction = l[['no.extinction']]
post.extinction = l[[scenario]]

#join all the dataframes where NA means the specimens fell into the bill depth extinction range
foo = join_all(list(pre.extinction = pre.extinction, post.extinction = post.extinction),by='ID',type='full')
colnames(foo) = c('ID','pre.extinction','post.extinction')

d = merge(data[,c('ID','Taxa.Lack','Taxa.New')],foo,by='ID',all=T)

#change class columns to factors with levels equal to classes
for(i in 4:ncol(d)){
	x = d[,i]
	levels = sort(as.numeric(unique(x[!is.na(x)])))
	d[,i] = factor(x,levels=levels)	
}


##################################################
# 2.3) Fig. 3 plotting

#color key for morphspecies
color.key = data.frame(col=c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854'),PreExt=c(1,2,3,4,NA),PostExt=c(1,2,3,4,5))

#shape key for all Geospiza species
taxa.key = data.frame(Taxa.New = c("G. acutirostris","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. propinqua","G. scandens","G. septentrionalis"),Taxa.Lack = c("G. difficilis","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. conirostris","G. scandens","G. difficilis"),pch=1:9)

#merge PCA data for pre.extinction dataset with mclust clustering data
mclust = mclust.data[['no.extinction']]$mclust
extant = mclust.data[['no.extinction']]$extant
pre = data.frame(ID=extant,mclust$data)
data = merge(d,pre[,c('ID','PC1','PC2','PC3','PC4')],by='ID',all=T)

#key for pre-extinction and post-extinction names
pre.ext.key = data.frame(mclust=1:4,name=c('G. fuliginosa','G. fortis 2','G. fortis 1','G. magnirostris'),pch=c(21,3,4,8),col=c('black','grey75','grey75','grey75'))
post.ext.key = data.frame(mclust=1:5,name=c('G. fuliginosa 2','G. fortis 2','G. fuliginosa 1','G. fortis 1','G. magnirostris'),pch=c(22,3,24,4,8),col=c('black','grey75','black','grey75','grey75'))

png(paste0('Fig.3.png'),width=6,height=3,units='in',res=600,bg='transparent')
	#dev.new(width=6,height=3.25)
	m = matrix(nrow=2,ncol=2)
	m[1,] = c(1,2)
	m[2,] = c(3,3)
	layout(m,widths=c(2,2),heights=c(2,0.25))
	
	par(mar=c(1,3.0,1.5,0.5))
	
	trait.x = 'PC1'
	trait.y = 'PC2'
	
	#Pre-extinction
	pre.ext = data[,c('ID','Taxa.New','PC1','PC2','PC3','PC4','pre.extinction')]
	col = mapvalues(as.numeric(pre.ext[,'pre.extinction']), pre.ext.key[,'mclust'],pre.ext.key[,'col'])
	pch = mapvalues(as.numeric(pre.ext[,'pre.extinction']), pre.ext.key[,'mclust'],pre.ext.key[,'pch'])
	
	xlim = c(range(pre.ext[,trait.x],na.rm=T)[2]*1.2,range(pre.ext[,trait.x],na.rm=T)[1]*1.1)
	ylim = c(range(pre.ext[,trait.y],na.rm=T)[1]*1.5,range(pre.ext[,trait.y],na.rm=T)[2]*1.3)
	
	plot(pre.ext[,trait.x], pre.ext[,trait.y],xlim=xlim,ylim=ylim,type="n", xlab="", ylab="",axes=F)
	axis(1,padj=-1,cex.axis=.75);		mtext(trait.x,1,line=1.75,cex=1)
	axis(2,padj= 1,cex.axis=.75);		mtext(trait.y,2,line=1.5,cex=1)
	mtext('A',3,cex=1.25,adj=0)
	mtext('Pre-Extinction',3)
	
	points(pre.ext[ ,trait.x], pre.ext[ ,trait.y],cex=1.5, pch=pch, col=col, bg='transparent',lwd=1)
	
	e = 2
	for(e in 1:length(levels(pre.ext[,'pre.extinction']))){
		x = pre.ext[pre.ext[,'pre.extinction'] %in% e, ]
		mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
		sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
		ellipse = ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=0.5,lty=1,col='black')
		
		name = pre.ext.key[pre.ext.key$mclust %in% e, 'name']
		col = pre.ext.key[pre.ext.key$mclust %in% e, 'col']
		if(name == 'G. magnirostris'){y.offset = .21}else{y.offset = .2}
		text(mu[1],mu[2]-y.offset,name,col= col,font=3,cex=.9)
		#segments(ellipse[which.min(ellipse[,2]),1],min(ellipse[,2]),mu[1],mu[2]-y.offset+.015,col=col)
	}
	
	#Post-extinction
	post.ext = data[,c('ID','Taxa.New','PC1','PC2','PC3','PC4','post.extinction')]
	col = mapvalues(as.numeric(post.ext[,'post.extinction']), post.ext.key[,'mclust'],post.ext.key[,'col'])
	pch = mapvalues(as.numeric(post.ext[,'post.extinction']), post.ext.key[,'mclust'],post.ext.key[,'pch'])
	
	xlim = c(range(post.ext[,trait.x],na.rm=T)[2]*1.2,range(post.ext[,trait.x],na.rm=T)[1]*1.4)
	ylim = c(range(post.ext[,trait.y],na.rm=T)[1]*1.5,range(post.ext[,trait.y],na.rm=T)[2]*1.3)
	
	plot(post.ext[,trait.x], post.ext[,trait.y],xlim=xlim,ylim=ylim,type="n", xlab="", ylab="",axes=F)
	axis(1,padj=-1,cex.axis=.75);		mtext(trait.x,1,line=1.75,cex=1)
	axis(2,padj= 1,cex.axis=.75);		mtext(trait.y,2,line=1.5,cex=1)
	mtext('B',3,cex=1.25,adj=0)
	mtext('Post-Extinction',3)

	extinct = is.na(data[,'post.extinction'])
	
	geo = read.table("Geospiza.data.csv", header=T, sep=",")
	geo[paste0(geo$Institution,geo$Museum.Number) %in% post.ext[extinct & post.ext$Taxa.New %in% 'G. fuliginosa','ID'], ]
	head(post.ext)
	
	
	#label extinct Geospiza fuliginosa
	x = post.ext[post.ext[extinct,'PC2'] %in% max(post.ext[extinct,'PC2'],na.rm = TRUE),'PC1']
	y = max(post.ext[extinct,'PC2'],na.rm = TRUE)
	segments(x+.1,y,x+0.3,y+.015,col='grey50')
	text(x+0.6,y+.055,'extinct',cex=.65,adj=0,col='grey50')
	text(x+0.6,y+.03,'G. fuliginosa',cex=.65,adj=0,font=3,col='grey50')
	
	#plot extinct Geospiza fuliginosa
	points(post.ext[extinct,trait.x], post.ext[extinct,trait.y],cex=1.5 , pch=21,col='grey75',bg='transparent',lwd=.5)
	#increase line width of G. fuliginosa 1 and 2
	lwd = rep(1,nrow(post.ext))
	lwd[post.ext$post.extinction %in% c(1,3)] = 1.5
	#plot extant Geospiza	
	points(post.ext[ ,trait.x], post.ext[ ,trait.y],cex=1.5, pch=pch, col=col, bg='transparent', lwd=lwd)

	e = 2
	for(e in 1:length(levels(post.ext[,'post.extinction']))){
		x = post.ext[post.ext[,'post.extinction'] %in% e, ]
		mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
		sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
		ellipse = ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=0.5,lty=1,col='black')
		
		name = post.ext.key[post.ext.key$mclust %in% e, 'name']
		col = post.ext.key[post.ext.key$mclust %in% e, 'col']
		if(name == 'G. magnirostris'){y.offset = .21}else{y.offset = .2}
		text(mu[1],mu[2]-y.offset,name,col= col,font=3,cex=.9)
		#segments(ellipse[which.min(ellipse[,2]),1],min(ellipse[,2]),mu[1],mu[2]-y.offset+.015,col=col)
	}

	dev.off()
