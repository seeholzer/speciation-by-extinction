clem = read.csv('~/Dropbox/Clements/Clements-Checklist-v2019-August-2019.csv',stringsAsFactors=F)
table(clem$category)

#clem = clem[grepl('Cranioleuca',clem$scientific.name), ]

#get total list of species
sp  = clem[clem$category %in% 'species','scientific.name']

#get counts of number of species for all species that have subspecies
ssp = clem[clem$category %in% c('group (monotypic)','subspecies'),'scientific.name'] 
x = t(sapply(strsplit(ssp,' '),'['))
ssp.summary = data.frame(table(apply(x[,1:2],1,paste,collapse=' ')))
colnames(ssp.summary) = c('sp','Nssp')

#number of species with 3 or more subspecies
nrow(ssp.summary[ssp.summary$Nssp >= 3, ])
#total number of species
length(sp)



