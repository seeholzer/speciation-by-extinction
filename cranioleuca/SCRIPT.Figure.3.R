##################################################
#	Description: Figure 3 - plot of mclust and BFD* results
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#load packages
library(RColorBrewer)
library(colorspace)
source('~/Dropbox/SBE/speciation-by-extinction_github/supporting.functions/FUN.add.alpha.R', chdir = TRUE)

#load summary tables for mclust and bfd
mclust = read.delim('~/Dropbox/SBE/SBE.in.silico/mclust/results.mclust.txt',stringsAsFactors=F)
bfd = read.delim('~/Dropbox/SBE/SBE.in.silico/results.bfd.txt')

#modify bfd table 
bfd$scenario.name = paste0(bfd$Dataset,'_',bfd$Version)
bfd$trait.set = 'SNPs'
bfd = bfd[,c('scenario.name','trait.set','MLE.1sp','MLE.2sp','BF')]

#key for plotting parameters 
combos = expand.grid(c('v1','v2','v3'),c('combined','SNPs'))
key = data.frame(expand.grid(c('v1','v2','v3'),c('combined','SNPs')),col=NA,pch=NA,lty=NA,x=NA,y=NA,stringsAsFactors=F)
colnames(key) = c('version','trait.set','col','pch','lty','x','y')

key[key$version %in% 'v1','col'] = brewer.pal(8,'Set2')[1]
key[key$version %in% 'v2','col'] = brewer.pal(8,'Set2')[6]
key[key$version %in% 'v3','col'] = brewer.pal(8,'Set2')[2]

key[key$trait.set %in% 'combined','pch'] = 24
key[key$trait.set %in% 'SNPs','pch'] = 21

key[key$trait.set %in% 'combined','lty'] = 3
key[key$trait.set %in% 'SNPs','lty'] = 1

key[key$version %in% 'v1','x'] = 2
key[key$version %in% 'v2','x'] = 2.5
key[key$version %in% 'v3','x'] = 3

key[key$trait.set %in% 'combined','y'] = 90
key[key$trait.set %in% 'SNPs','y'] = 80


##########################
#Figure 3
##########################

png(paste0('~/Dropbox/SBE/figures/Fig.3.png'),
	width=5,height=5,units='in',res=300,bg='transparent')
#dev.new(width=5,height=5)
	par(mar=c(3.5,5,0,1))
	#xlim = c(0.5,4.5); ylim = c(-50,150)
	xlim = c(0.5,4.5); ylim = c(-50,110)
	if(min(mclust$deltaBIC) < ylim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	if(max(mclust$deltaBIC) > ylim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	
	plot(mclust[ ,'deltaBIC'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
	abline(h=0,lty=2)
	axis(1,at=c(1,2,3,4),labels=c('FULL','SBE1','SBE2','SBE3'),tick=T,lty=NULL,col.axis = 'grey25', col.ticks = 1)
	mtext('Exctintion Scenarios',1,line=2.5,cex=1.25,col='black')

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	#axis(4,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	mtext(expression(paste(Delta,"BIC / Bayes Factor Support")),2,line=3.5,cex=1.25,col='black')
	mtext('One',2,at=-25,line=2.5,col='grey25',cex=1.25)
		mtext('Cluster/Species',2,at=-25,line=2,col='grey25',cex=.75)
		#col='#E41A1C'
	mtext('Two',2,at= 50,line=2.5,col='grey25',cex=1.25)
		mtext('Clusters/Species',2,at=50,line=2,col='grey25',cex=.75)
		#col='#377EB8'
		
v = 'v2'
for(v in c('v3','v2','v1')){
	
	col = key[key$version %in% v & key$trait.set %in% 'combined','col']	
	pch = key[key$version %in% v & key$trait.set %in% 'combined','pch']	
	lty = key[key$version %in% v & key$trait.set %in% 'combined','lty']
	deltaBIC = mclust[mclust$trait.set %in% 'combined' & grepl(v,mclust$scenario.name),'deltaBIC']
	deltaBIC = jitter(deltaBIC,1.5)
	lines(1:4, deltaBIC,col=col,lty=lty)
	points(1:4, deltaBIC,col='grey50',bg=col,pch=pch,cex=2,lwd=.5)

	col = key[key$version %in% v & key$trait.set %in% 'SNPs', "col"]
	pch = key[key$version %in% v & key$trait.set %in% 'SNPs', "pch"]
	lty = key[key$version %in% v & key$trait.set %in% 'SNPs', "lty"]
	BF = bfd[bfd$trait.set %in% 'SNPs' & grepl(v, bfd$scenario.name), "BF"]
	lines(1:4, BF, col = col, lty = lty)
	points(1:4, BF, col = 'grey50', bg = col, pch = pch, cex = 2, lwd=.5)

}

#legend
	tmp = key
	tmp[,'x'] = tmp[,'x'] - 0.25 
	tmp[,'y'] = tmp[,'y'] + 10 
	xlines = sort(unique(tmp[,'x']))
	ylines = sort(unique(tmp[,'y']))
	text(0.45,ylines,c('SNPs','Phenotype'),adj=c(0,.5))
	text(xlines,max(ylines) + 10 ,c('V1','V2','V3'),adj=c(.5,.5))
	for(i in 1:nrow(tmp)){
		x = tmp[i,'x']
		y = tmp[i,'y']
		lty = tmp[i,'lty']
		pch = tmp[i,'pch']
		col = tmp[i,'col']
		bg	= tmp[i,'bg']
		segments(x-0.2,y,x+0.2,y,lty=lty,col=col)
		points(x,y,pch=pch,col='grey50',bg=col,cex=1.25,lwd=.5)
		}


dev.off()


