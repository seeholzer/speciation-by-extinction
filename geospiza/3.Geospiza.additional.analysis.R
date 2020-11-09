mclust.data = get(load('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/mclust.data.rda'))
data = read.table("Geospiza.data.csv", header=T, sep=",")
head(data)

names(mclust.data)
lapply(mclust.data,names)

gapsize = 1
for(gapsize in 1:5){
	
	tmp = mclust.data[[paste0(gapsize,'mm.gap')]]
	gap = '6-7'
	for(gap in names(tmp)){
		
		mclust = tmp[[gap]]$mclust
		extant = tmp[[gap]]$extant
		head(extant)
		cbind(extant$Museum.No.,mclust$classification)
		
		
		
		
		head(extant)
		unique(extant$New_Taxonomy)
		
		tmp[[gap]]$mclust$data
		tmp[[gap]]$mclust$classification
				
		trait.x <- 1
		trait.y <- 2
		plot(Mcluster.morpho.data.ln.pca.subset$data[,trait.x], Mcluster.morpho.data.ln.pca.subset$data[,trait.y],
		xlab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.y],
		cex.axis=1.5, cex.lab=1.5, bty="n", col=col.BestModel, pch=21,
		cex=1.5, lwd=2)
		for(i in 1:Mcluster.morpho.data.ln.pca.subset$G)
		{
			points(ellipse(Z.all[c(trait.x,trait.y),c(trait.x,trait.y),i], centre = m.all[c(trait.x,trait.y),i], level = 0.95, npoints = 10000), type="l")
		}

		
		
		
		
		
		
		
	}
	
}







