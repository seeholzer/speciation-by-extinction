##################################################
#	Description: Figure 3 - plot of mclust and BFD* results
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Change working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/')

#load packages
library(RColorBrewer)
library(colorspace)
source('supporting.functions/FUN.add.alpha.R', chdir = TRUE)

#load summary tables for mclust and bfd
mclust = read.delim('cranioleuca/results.mclust.txt',stringsAsFactors=F)
bfd = read.delim('cranioleuca/results.bfd.txt')

#modify bfd table 
bfd$scenario.name = paste0(bfd$Dataset,'_',bfd$Version)
bfd$trait.set = 'SNPs'
bfd = bfd[,c('scenario.name','trait.set','MLE.1sp','MLE.2sp','BF')]

scenario.key = cbind(c('FULL','SBE1','SBE2','SBE3'),c('Pre-Extinction','Post-Extinction 1','Post-Extinction 2','Post-Extinction 3'))

##########################
#Figure 2 Horizontal
##########################

png(paste0('cranioleuca/Fig.2.empirical.support.horizontal.png'),
	width=5.5,height=1.5,units='in',res=300,bg='transparent')
	#dev.new(width=5.5,height=1.5)
	par(mar=c(1.5,3,0.1,1))
	#xlim = c(0.5,4.5); ylim = c(-50,150)
	xlim = c(0.9,4.3); ylim = c(-50,110)
	if(min(mclust$deltaBIC) < ylim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	if(max(mclust$deltaBIC) > ylim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	
	plot(1:4, xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F, frame=F)
	segments(xlim[1]-0.5,0,5,0,,lty=2)
	labels = scenario.key[,2]
	a = sapply(strsplit(labels,'-'),'[',1)
	b = sapply(strsplit(labels,'-'),'[',2)
	axis(1,at=c(1,2,3,4),labels=paste0(a,'-'),tick=F,lty=NULL,col.axis = 'black', col.ticks = 1, las=1, line=-1.5,cex.axis=.65)
	axis(1,at=c(1,2,3,4),labels=b,tick=F,lty=NULL,col.axis = 'black', col.ticks = 1, las=1, line=-1.0,cex.axis=.65)
	mtext('Extinction Scenarios',1,line=.6,cex=.75,col='grey25',adj=.475)

	axis(2,at=c(-50,0,50,100),labels=c(-50,0,50,100),col.axis = 'black',cex.axis=.75,padj=1.25)
	mtext("Bayes Factor or",2,line=2,cex=.75,col='grey25',adj=c(0))
	mtext(expression(paste(Delta,"BIC Support")),2,line=1.25,cex=.75,col='grey25',adj=c(0))

	#One species, two species 		
	x = .8
	text(x, 50,expression(bold('Two')),col='#E41A1C',cex=.5,srt=90)
	text(x+.05, 50,expression(bold('Species')),col='#E41A1C',cex=.5,srt=90)

	text(x, -30,expression(bold('One')),col='#377EB8',cex=.5,srt=90)
	text(x+.05, -30,expression(bold('Species')),col='#377EB8',cex=.5,srt=90)

	
	v = 'v2'

	col = 'grey50'
	pch = 24	
	lty = 3
	deltaBIC = mclust[mclust$trait.set %in% 'combined' & grepl(v,mclust$scenario.name),'deltaBIC']
	deltaBIC = jitter(deltaBIC,1.5)

	text(4.25,deltaBIC[4]+10,,labels='phenotype',cex=.5,col='grey50')
	segments(5,deltaBIC[4],4,deltaBIC[4],col='grey50')
	
	lines(1:4,deltaBIC,col=col,lty=lty)
	points(1:4,deltaBIC,col=col,bg='grey50',pch=pch,cex=1.5,lwd=.5)
	
	col = 'grey50'
	pch = 21
	lty = 1
	BF = bfd[bfd$trait.set %in% 'SNPs' & grepl(v, bfd$scenario.name), "BF"]
	
	text(4.25,BF[4]+10,labels='SNPs',cex=.5,col='grey50')
	segments(5,BF[4],4,BF[4],,col='grey50')
	
	lines(1:4,BF, col = col, lty = lty)
	points(1:4,BF, col = col, bg = 'grey50', pch = pch, cex = 1.5, lwd=.5)


dev.off()



# ##########################
# #Figure 2 vertical
# ##########################

# png(paste0('~/Dropbox/SBE/figures/Fig.empirical.support.png'),
	# width=2.5,height=5,units='in',res=300,bg='transparent')
	# #dev.new(width=2.5,height=5)
	# par(mar=c(5,4,0,1))
	# #xlim = c(0.5,4.5); ylim = c(-50,150)
	# xlim = c(-50,110); ylim = c(1,4.25)
	# if(min(mclust$deltaBIC) < xlim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	# if(max(mclust$deltaBIC) > xlim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
	
	# plot(1:4, xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F, frame=F)
	# segments(0,0,ylim[1],4,lty=2)
	# axis(2,at=c(4,3,2,1),labels=c('FULL','SBE1','SBE2','SBE3'),tick=T,lty=NULL,col.axis = 'grey25', col.ticks = 1, las=1, hadj=.8)
	# mtext('Extinction Scenarios',2,line=2.9,cex=1.25,col='black')

	# axis(1,at=c(-50,0,50,100),labels=c(-50,0,50,100))
	# mtext(expression(paste(Delta,"BIC / Bayes Factor Support")),1,line=4,cex=1,col='black',adj=c(1))
	# mtext(expression(bold('One')),1,at=-40,line=2,col='grey25',cex=1)
		# mtext(expression(bold('Species')),1,at=-40,line=2.75,col='grey25',cex=.75)
		# # #col='#E41A1C'
	# mtext(expression(bold('Two')),1,at=60,line=2,col='grey25',cex=1)
		# mtext(expression(bold('Species')),1,at=60,line=2.75,col='grey25',cex=.75)
		# # #col='#377EB8'
		
# v = 'v2'
# #for(v in c('v3','v2','v1')){
	
	# #col = key[key$version %in% v & key$trait.set %in% 'combined','col']	
	# col = 'grey50'
	# pch = key[key$version %in% v & key$trait.set %in% 'combined','pch']	
	# lty = key[key$version %in% v & key$trait.set %in% 'combined','lty']
	# deltaBIC = mclust[mclust$trait.set %in% 'combined' & grepl(v,mclust$scenario.name),'deltaBIC']
	# deltaBIC = jitter(deltaBIC,1.5)

	# text(-30,4.2,labels='phenotype',cex=.5)
	# segments(-30,4.15,deltaBIC[1],4)
	
	# lines(deltaBIC,4:1,col=col,lty=lty)
	# points(deltaBIC,4:1,col=col,bg='grey50',pch=pch,cex=2,lwd=.5)

	# #col = key[key$version %in% v & key$trait.set %in% 'SNPs', "col"]
	# col = 'grey50'
	# pch = key[key$version %in% v & key$trait.set %in% 'SNPs', "pch"]
	# lty = key[key$version %in% v & key$trait.set %in% 'SNPs', "lty"]
	# BF = bfd[bfd$trait.set %in% 'SNPs' & grepl(v, bfd$scenario.name), "BF"]
	
	# text(15,4.2,labels='SNPs',cex=.5)
	# segments(15,4.15,BF[1],4)
	
	# lines(BF,4:1, col = col, lty = lty)
	# points(BF,4:1, col = col, bg = 'grey50', pch = pch, cex = 2, lwd=.5)



# dev.off()

