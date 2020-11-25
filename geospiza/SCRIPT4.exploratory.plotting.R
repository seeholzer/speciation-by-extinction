##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	Script 4 - exploratory plotting
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
# 4.1) Read data
data = read.table("Geospiza.data.csv", header=T, sep=",")
data$ID = paste0(data$Institution,data$Museum.Number)
colnames(data)[1] = 'Taxa.Lack'
colnames(data)[6] = 'Taxa.New'
colnames.to.keep = c('ID','Taxa.Lack','Taxa.New','Wing','Tail','Blength','Bdepth','Bwidth','Tarsus')
data = data[,colnames.to.keep]

#load mclust results
mclust.data = get(load('mclust.Bdepth.SML.groundfinches.rda'))
plots.dir = 'plots.Bdepth.SML.groundfinches'

#Data for all Geospiza species (LEGACY)
# mclust.data = get(load('mclust.Bdepth.rda')) 
# plots.dir = 'plots.Bdepth.all'

##################################################
# 4.2) Check correlation between pre and post PCA space
#pre and post extinction PCA morphospace are virtually identical
#strongest correlation for PC1 reduced for PC2-PC4 but still strong
#NO NEED TO PLOT IN POST-EXTINCTION MORPHOSPACE

mclust = mclust.data[['no.extinction']]$mclust
extant = mclust.data[['no.extinction']]$extant
pre = data.frame(ID=extant,mclust$data)
colnames(pre)[-1] = paste0('pre.',colnames(pre)[-1])

colnames = c('scenario','PC1','PC2','PC3','PC4')
sum = data.frame(matrix(nrow=length(names(mclust.data)[-1]),ncol=length(colnames)))
colnames(sum) = colnames
sum$scenario = names(mclust.data)[-1]

for(i in 1:nrow(sum)){

	mclust = mclust.data[[i]]$mclust
	extant = mclust.data[[i]]$extant
	post = data.frame(ID=extant,mclust$data)
	colnames(post)[-1] = paste0('post.',colnames(post)[-1])
	
	x = merge(pre,post,by='ID')
	for(A in 1:4){
		sum[i,paste0('PC',A)] = cor(x[,paste0('pre.PC',A)],x[,paste0('post.PC',A)],method='pearson')
	}
	
}

apply(sum[,-1],2,mean)
apply(sum[,-1],2,var)



##################################################
# 4.3) Merge all mclust classifications and other relevant data

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
# 4.4) Figure 4 plotting
#	plot the extinction events and old and new morphogroups in the morphospace for pre-extinction dataset

#color key for morphspecies
color.key = data.frame(group = c(1:10),col=c(brewer.pal(8, 'Set2'),'darkorchid4','cadetblue1'))

#shape key for all Geospiza species
taxa.key = data.frame(Taxa.New = c("G. acutirostris","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. propinqua","G. scandens","G. septentrionalis"),Taxa.Lack = c("G. difficilis","G. conirostris","G. difficilis","G. fortis","G. fuliginosa","G. magnirostris","G. conirostris","G. scandens","G. difficilis"),pch=1:9)

#merge PCA data for pre.extinction dataset with mclust clustering data
mclust = mclust.data[['no.extinction']]$mclust
extant = mclust.data[['no.extinction']]$extant
pre = data.frame(ID=extant,mclust$data)
foo = merge(d,pre[,c('ID','PC1','PC2','PC3','PC4')],by='ID',all=T)

scenarios = names(mclust.data)

for(i in scenarios[-1]){
	cat(i,'\n')
	gap = as.numeric(strsplit(i,'-')[[1]][2]) - as.numeric(strsplit(i,'-')[[1]][1])
	
	if(!dir.exists(plots.dir)) dir.create(paste0(plots.dir,'/',gap))
	png(paste0(plots.dir,'/',i,'.png'),width=6,height=6,units='in',res=500,bg='transparent')
	#dev.new(width=6,height=6)

	m = matrix(nrow=3,ncol=3)
	m[1,] = c(1,7,9)
	m[2,] = c(2,4,8)
	m[3,] = c(3,5,6)
	layout(m,widths=c(2,2,2),heights=c(2,2,2))

	par(mar=c(3,3,0,0))
	noext = foo[,c('ID','Taxa.New','PC1','PC2','PC3','PC4','no.extinction')]
	tmp = foo[,c('ID','Taxa.New','PC1','PC2','PC3','PC4',i)]
	cols = as.character(mapvalues(tmp[,i],color.key[,1],color.key[,2]))
	pch = as.numeric(mapvalues(tmp[,'Taxa.New'],taxa.key[,1],taxa.key[,3]))

	#all pairwise combos of PC axes
	traits.x = c('PC1','PC1','PC1','PC2','PC2','PC3')
	traits.y = c('PC2','PC3','PC4','PC3','PC4','PC4')
	
	t = 1
	for(t in 1:6){
		trait.x = traits.x[t]
		trait.y = traits.y[t]
		
		plot(tmp[,trait.x], tmp[,trait.y],type="n", xlab="", ylab="",axes=F)
		axis(1,padj=-1);mtext(trait.x,1,line=1.5,cex=.75)
		axis(2,padj= 1);mtext(trait.y,2,line=1.5,cex=.75)
		
		extinct = is.na(tmp[,i])
		points(tmp[ extinct,trait.x], tmp[ extinct,trait.y],cex=1.25, pch=21,bg=add.alpha('grey75',.5),lwd=.75)
		points(tmp[!extinct,trait.x], tmp[!extinct,trait.y],cex=1.5 , pch=pch[!extinct],col=cols[!extinct],lwd=.5)
			
		for(e in 1:length(levels(noext[,'no.extinction']))){
			x = noext[noext[,'no.extinction'] %in% e, ]
			mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
			sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
			ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=1,lty=3,col='black')
		}

		e = 1
		for(e in 1:length(levels(tmp[,i]))){
			x = tmp[tmp[,i] %in% e, ]
			mu = c(mean(x[,trait.x]), mean(x[,trait.y]))
			sigma = var(x[,c(trait.x,trait.y)])  # returns a variance-covariance matrix.
			if(!all(is.na(sigma))){
				ellipse(mu, sigma, alpha=0.05, npoints = 200, newplot = FALSE,lwd=.5,col='black')
				}
		}
	}

	par(mar=c(0, 0, 0, 0))
	plot(c(0,1), c(0,1),type='n',axes=F)
	legend(.5,.5,legend=c(taxa.key[,'Taxa.New'],'extinct'),pch=c(taxa.key[,'pch'],21),pt.bg=add.alpha('grey75',1),bty='n',cex=.9,pt.cex=1.25,pt.lwd=.75,xjust=.5,yjust=.5)

	par(mar=c(0, 0, 0, 0))
	plot(c(0,1), c(0,1),type='n',axes=F)
	text(.5,1,'morpho-species',adj=c(.5,.5),cex=1.25)
	text(.5,.9,'diagnosed after',adj=c(.5,.5),cex=1.25)
	text(.5,.8,'simulated extinction',adj=c(.5,.5),cex=1.25)
	
	legend(.5,.6,legend=levels(tmp[,i]),pch=16,col=color.key[color.key$group %in% levels(tmp[,i]),2],bty='n',ncol=4,cex=1.5,pt.cex=2,xjust=.5,yjust=.5,adj=.65,text.width=.01)

	par(mar=c(0, 0, 0, 0))
	plot(c(0,1), c(0,1),type='n',axes=F)
	text(.5,.9,paste0(i,' mm'),cex=2)

	x = 0.25; y = 0.4 
	draw.ellipse(x,y,a=.225,b=.15,lty=2)
	text(x,y,paste0('N = ',length(levels(noext[,'no.extinction']))),cex=1.5)
	text(x,y+0.3,'morpho-species',cex=.75)
	text(x,y+0.25,'before',cex=.75)
	text(x,y+0.2,'sim. extinction',cex=.75)

	x = 0.75; y = 0.4 
	draw.ellipse(x,y,a=.225,b=.15,lty=1)
	text(x,y,paste0('N = ',length(levels(tmp[,i]))),cex=1.5)
	text(x,y+0.3,'morpho-species',cex=.75)
	text(x,y+0.25,'after',cex=.75)
	text(x,y+0.2,'sim. extinction',cex=.75)

	dev.off()

}

