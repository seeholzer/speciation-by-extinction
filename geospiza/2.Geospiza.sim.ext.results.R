##################################################
#	Description: Plot the results of simulated extinction in Darwin's Finches from script 1.Geospiza.sim.ext.R
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/')

#load results
results = get(load('results.rda'))

#Version 1
# png('Figure.Geospiza.extinction.gaps.1.2.3.png',width=3,height=6.25,units='in',res=300,bg='transparent')
# gaps = c(1,2,3)

png('Figure.Geospiza.extinction.gaps.1.3.5.png',width=3.25,height=6.25,units='in',res=300,bg='transparent')
gaps = c(1,3,5)

#Version 2
#dev.new(width=3.25,height=6.5)
m = t(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3))
m = rbind(m,c(7,7)) # create space for x axis label

layout(mat = m, heights = c(2, 2, 2, 0.5), widths = c(.75,2.5))
#layout.show(7)

i = 1
for(i in 1:length(gaps)){
	gap = gaps[i]	
	tmp = results[results$gap %in% gap, ]

	par(mar=c(3,4,1,0))
	xlim = c(0,1)
	ylim = c(5,11)
	xat = seq(0 +.5,xlim[2]-.5,1)
	plot(1,1,xlim=xlim,ylim=ylim,xlab='',ylab='',axes=F)
	abline(h=8,lty=1,col='grey')
	segments(xat,0,xat,8,col='grey')
	points(xat,8,pch=21,cex=1.5,col='black',bg='grey85',lwd=1)
	axis(1,labels='',at=xat,cex.axis=1,lwd.ticks=2,las=1,padj=-0.5)
	axis(2,at=c(5:10),label=c(5:10))
	text(xat-.1,par("usr")[3], labels = 'No', 			xpd = TRUE,srt=45, adj = c(2.25, 1  ),cex=.8) 
	text(xat-.1,par("usr")[3], labels = 'Extinction', 	xpd = TRUE,srt=45, adj = c(1  , 2.5),cex=.8) 
	mtext('N Inferred Species', 2,2.5,cex=.75)

	text(1,10.85,paste0(LETTERS[i],'.'),pos=2,cex=1.25)

	par(mar=c(3,0,1,2))
	xlim = c(0,nrow(tmp)-1)
	ylim = c(5,11)
	xat = seq(xlim[1],xlim[2],1)
	
	plot(1,1,xlim=xlim,ylim=ylim,xlab='',ylab='',axes=F)
	abline(h=8,lty=1,col='grey')
	segments(xat,0,xat,tmp$max.BIC.G,col='grey')
	points(xat,tmp$max.BIC.G,pch=16,cex=1.5)
	axis(1,at=xat,label=tmp[,'extinction.range'],cex.axis=1,lwd.ticks=2,las=2)

	text(-.9,10.85,paste0(gap,' mm ext. range'),pos=4,cex=1.25)

}

	par(mar=c(0,0,0,0))
	plot(0:1,0:1,xlab='',ylab='',axes=F,type='n')
	text(.6,.4,'Extinction Events by Bill Depth (mm)',cex=1.3)
	
dev.off()

