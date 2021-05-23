##################################################
#	Description: Create summary tables of reuslts of quantitative species delimitation for 
#	phenotypic and genetic data. Combined manually into Table S1.
#
##################################################

library(plyr)
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

key = cbind(c('FULL','SBE1','SBE2','SBE3'),c('Pre-Extinction','Post-Extinction 1','Post-Extinction 2','Post-Extinction 3'))

#load summary tables for mclust and bfd
mclust = read.delim('results.mclust.txt',stringsAsFactors=F)
mclust = mclust[mclust$trait.set %in% 'combined', ]

v1 = mclust[grep('v1',mclust$scenario.name),c('scenario.name','BIC.G1','BIC.G2','deltaBIC')]
colnames(v1)[-1] = paste('v1',colnames(v1)[-1])

v2 = mclust[grep('v2',mclust$scenario.name),c('scenario.name','BIC.G1','BIC.G2','deltaBIC')]
colnames(v2)[-1] = paste('v2',colnames(v2)[-1])

v3 = mclust[grep('v3',mclust$scenario.name),c('scenario.name','BIC.G1','BIC.G2','deltaBIC')]
colnames(v3)[-1] = paste('v3',colnames(v3)[-1])

mor = cbind(v1,v2[,-1],v3[,-1])
colnames(mor)[1] = 'Scenario'
colnames(mor) = gsub('BIC.G1','BIC.1sp',colnames(mor))
colnames(mor) = gsub('BIC.G2','BIC.2sp',colnames(mor))
mor$Scenario = gsub('_v1','',mor$Scenario)

mor[,-1] = apply(mor[,-1],2,round)

mor[,1] = mapvalues(mor[,1],key[,1],key[,2])

####################################################################
bfd = read.delim('results.bfd.txt')

v1 = bfd[bfd$Version %in% 'v1',c('Dataset','MLE.1sp','MLE.2sp','BF')]
colnames(v1)[-1] = paste('v1',colnames(v1)[-1])

v2 = bfd[bfd$Version %in% 'v2',c('Dataset','MLE.1sp','MLE.2sp','BF')]
colnames(v2)[-1] = paste('v2',colnames(v2)[-1])

v3 = bfd[bfd$Version %in% 'v3',c('Dataset','MLE.1sp','MLE.2sp','BF')]
colnames(v3)[-1] = paste('v3',colnames(v3)[-1])

gen = cbind(v1,v2[,-1],v3[,-1])
colnames(gen)[1] = 'Scenario'

gen[,-1] = apply(gen[,-1],2,round)

gen[,1] = mapvalues(gen[,1],key[,1],key[,2])


#These tables are combined manually to create Table S1
write.table(mor,'results.mclust.formatted.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(gen,'results.bfd.formatted.txt',row.names=F,col.names=T,sep='\t',quote=F)




