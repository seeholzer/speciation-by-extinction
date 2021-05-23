##################################################
#	Description: Script to calculate the number of bird species with three or more subspecies
#
#	Author: N/A
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################
setwd('~/Dropbox/SBE/speciation-by-extinction_github/bird.subspecies')

#load Clements checklist
clem = read.csv('Clements-Checklist-v2019-August-2019.csv',stringsAsFactors=F)

#get total list of species
sp  = clem[clem$category %in% 'species','scientific.name']

#get counts of number of subspecies for all species that have subspecies
ssp = clem[clem$category %in% c('group (monotypic)','subspecies'),'scientific.name'] 
x = t(sapply(strsplit(ssp,' '),'['))
ssp.summary = data.frame(table(apply(x[,1:2],1,paste,collapse=' ')))
colnames(ssp.summary) = c('sp','Nssp')


#number of species with 3 or more subspecies
nrow(ssp.summary[ssp.summary$Nssp >= 3, ])
#total number of species
length(sp)