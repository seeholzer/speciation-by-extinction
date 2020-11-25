##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	Script 3 - Figure 4
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza')

#load packages
library(plyr)
library(mixtools)
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza')
source('~/Dropbox/SBE/speciation-by-extinction_github/supporting.functions/FUN.add.alpha.R', chdir = TRUE)

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
#join all the dataframes where NA means the specimens fell into the bill depth extinction range
foo = join_all(l,by='ID',type='full')

d = merge(data[,c('ID','Taxa.Lack','Taxa.New')],foo,by='ID',all=T)

#change class columns to factors with levels equal to classes
for(i in 4:ncol(d)){
	x = d[,i]
	levels = sort(as.numeric(unique(x[!is.na(x)])))
	d[,i] = factor(x,levels=levels)	
}




##################################################
# 2.3) Figure 4 plotting

#color key for morphspecies
color.key = data.frame(col=c('#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854'),PreExt=c(1,2,3,4,NA),PostExt=c(1,2,3,4,5))

#shape key for all Geospiza species
taxa.key = data.frame(Taxa.New = c("G. acutirostris","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. propinqua","G. scandens","G. septentrionalis"),Taxa.Lack = c("G. difficilis","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. conirostris","G. scandens","G. difficilis"),pch=1:9)

#merge PCA data for pre.extinction dataset with mclust clustering data
mclust = mclust.data[['no.extinction']]$mclust
extant = mclust.data[['no.extinction']]$extant
pre = data.frame(ID=extant,mclust$data)
foo = merge(d,pre[,c('ID','PC1','PC2','PC3','PC4')],by='ID',all=T)

#specify scenario for plotting
scenario = "7.00-8.25"

png(paste0('Figure.4.png'),width=6,height=3,units='in',res=600,bg='transparent')
	#dev.new(width=6,height=3.25)
	m = matrix(nrow=2,ncol=2)
	m[1,] = c(1,2)
	m[2,] = c(3,3)
	layout(m,widths=c(2,2),heights=c(2,0.25))
	
	par(mar=c(2,3.0,1.5,0.5))
	
	trait.x = 'PC1'
	trait.y = 'PC2'
	
	#Pre-extinction
	noext = foo[,c('ID','Taxa.New','PC1','PC2','PC3','PC4','no.extinction')]
	cols = as.character(mapvalues(noext[,'no.extinction'], color.key$PreExt, color.key$col))
	pch = as.numeric(mapvalues(noext[,'Taxa.New'],taxa.key[,1],taxa.key[,3]))
	
	xlim = c(range(noext[,trait.x],na.rm=T)[1]*1.1,range(noext[,trait.x],na.rm=T)[2]*1.1)
	ylim = c(range(noext[,trait.y],na.rm=T)[1]*1.3,range(noext[,trait.y],na.rm=T)[2]*1.1)
	
	plot(noext[,trait.x], noext[,trait.y],xlim=xlim,ylim=ylim,type="n", xlab="", ylab="",axes=F)
	axis(1,padj=-1,cex.axis=.75);		mtext(trait.x,1,line=1.5,cex=1)
	axis(2,padj= 1,cex.axis=.75);		mtext(trait.y,2,line=1.5,cex=1)
	mtext('Pre-Extinction',3)

	points(noext[ ,trait.x], noext[ ,trait.y],cex=1.5, pch=pch, col=cols ,lwd=1.5)
	
	for(e in 1:length(levels(noext[,'no.extinction']))){
		x = noext[noext[,'no.extinction'] %in% e, ]
		mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
		sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
		ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=0.5,lty=1,col='black')
	}
	
	x = -0.5	;	y = -0.21
	text(x,y+.045,'Pre-Extinction',adj=c(0.5,0.5),cex=0.7)
	text(x,y+.025,'Morphospecies',adj=c(0.5,0.5),cex=0.7)
	col = color.key[color.key$PreExt %in% levels(noext[,'no.extinction']),'col']
	legend(x,y,legend=levels(noext[,'no.extinction']),pch=16,col= col,bty='n',ncol=4,cex=0.75,pt.cex=1.25,xjust=0.5,yjust=0.5,adj=0.65,text.width=0.01)
	

	#Post-extinction	
	tmp = foo[,c('ID','Taxa.New','PC1','PC2','PC3','PC4',scenario)]
	#change names of Post-extinction morphospecies so they line up with corresponding Pre-Extinction morphospecies
	key = data.frame(old = c(1,2,3,4,5), new = c(1,2,5,3,4))
	tmp[,scenario] = mapvalues(tmp[,scenario],key$old,key$new)

	cols = as.character(mapvalues(tmp[,scenario], color.key$PostExt, color.key$col))
	pch = as.numeric(mapvalues(tmp[,'Taxa.New'],taxa.key[,1],taxa.key[,3]))
	
	plot(tmp[,trait.x], tmp[,trait.y],xlim=xlim,ylim=ylim,type="n", xlab="", ylab="",axes=F)
	axis(1,padj=-1,cex.axis=.75);		mtext(trait.x,1,line=1.5,cex=1)
	axis(2,padj= 1,cex.axis=.75);		mtext(trait.y,2,line=1.5,cex=1)
	mtext('Post-Extinction',3)
	
	extinct = is.na(tmp[,scenario])
	points(tmp[ extinct,trait.x], tmp[ extinct,trait.y],cex=1.25, pch=21,col='black',bg='grey95',lwd=.15)
	points(tmp[!extinct,trait.x], tmp[!extinct,trait.y],cex=1.5 , pch=pch[!extinct],col=cols[!extinct],lwd=1.5)
		
	for(e in 1:length(levels(tmp[,scenario]))){
		x = tmp[tmp[,scenario] %in% e, ]
		mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
		sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
		ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=.5,lty=1,col='black')
	}

	x = -0.5	;	y = -0.21
	text(x,y+.045,'Post-Extinction',adj=c(0.5,0.5),cex=0.7)
	text(x,y+.025,'Morphospecies',adj=c(0.5,0.5),cex=0.7)
	col = color.key[color.key$PostExt %in% levels(tmp[,scenario]),'col']
	legend(x,y,legend=sort(levels(tmp[,scenario])),pch=16,col= col,bty='n',ncol=5,cex=0.75,pt.cex=1.25,xjust=0.5,yjust=0.5,adj=0.65,text.width=0.01)

	
 	#Taxonomic Legend
	par(mar=c(0, 0, 0, 0))
	plot(c(0,1), c(0,1),type='n',axes=F)

	x = 0.5;	y = 0.75
	text(x,y,'Current Taxonomy',adj=c(0.5,0.5),cex=0.75)
	spp = c('G. fuliginosa','G. fortis','G. magnirostris')
	legend = c(spp,'extinct')
	pch = c(as.numeric(mapvalues(spp,taxa.key[,'Taxa.New'],taxa.key[,'pch'])),21)
	legend(x,y - 0.4,legend=legend,text.font=c(3,3,3,1),pch=pch,pt.bg=add.alpha('grey75',.25),bty='n',cex=0.7,pt.cex=1.25,pt.lwd=1,xjust=0.5,yjust=0.5,ncol=4)


	dev.off()




