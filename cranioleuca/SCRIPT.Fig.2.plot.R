##################################################
#	SCRIPT.Fig.2.plots
#	Description: plot of mclust and BFD* results for scenarios FULL_v5, EXT_v4, EXT_v16, EXT_v24
#
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
bfd$trait.set = 'SNPs'
bfd = bfd[,c('scenario','trait.set','MLE_1sp','MLE_2sp','BF')]



scenario.key = cbind(c('FULL_v5','EXT_v4','EXT_v16','EXT_v24'),c('Pre-Extinction','Post-Extinction 1','Post-Extinction 2','Post-Extinction 3'))

##########################
#Figure 2 Horizontal
##########################

png(paste0('cranioleuca/Fig.2b.png'),
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

	scenarios = c('FULL','EXT_v4','EXT_v16','EXT_v24') 
	col = 'grey50'
	pch = 24	
	lty = 3
	foo = mclust[mclust$trait.set %in% 'combined' & mclust$scenario.name %in% scenarios, ]
	deltaBIC = foo[match(scenarios,foo$scenario.name),'deltaBIC']
	deltaBIC = jitter(deltaBIC,1.5)

	text(4.25,deltaBIC[4]+10,,labels='phenotype',cex=.5,col='grey50')
	segments(5,deltaBIC[4],4,deltaBIC[4],col='grey50')
	
	lines(1:4,deltaBIC,col=col,lty=lty)
	points(1:4,deltaBIC,col=col,bg='grey50',pch=pch,cex=1.5,lwd=.5)
	
	scenarios = c('FULL_v5','EXT_v4','EXT_v16','EXT_v24') 
	col = 'grey50'
	pch = 21
	lty = 1
	BF = bfd[bfd$scenario %in% scenarios, "BF"]
	
	text(4.25,BF[4]+10,labels='SNPs',cex=.5,col='grey50')
	segments(5,BF[4],4,BF[4],,col='grey50')
	
	lines(1:4,BF, col = col, lty = lty)
	points(1:4,BF, col = col, bg = 'grey50', pch = pch, cex = 1.5, lwd=.5)


dev.off()