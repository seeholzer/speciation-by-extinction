##################################################
#	Description: BFD species delimitation with simulated extinction in Cranioleuca antisiensis
#	Script 1 - Reformat Genetic Data 
#	Goal: 	1. Create nexus file of 2 individuals per population and 100 loci for input into Beauti 
#				to create template .xml for use in SCRIPT2.xml.editing.R 
#
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#load packages
library(dartR)
library(adegenet)
library(plyr)


#DO NOT RUN
#old code used to reformat genetic data
####	import phenotypic and meta-data
# meta = read.delim('~/Dropbox/SBE/cran/bfd/cran.meta.data.txt',stringsAsFactors=F)
####	import genetic data
# gen = get(load('~/Dropbox/SBE/cran/bfd/cran.snp.data.genind.Rdata'))
# gl = gi2gl(gen, parallel = TRUE)
# gl@pop = as.factor(mapvalues(gl@ind.names,meta[,'ID'],meta[,'Population']))
# save(gl,file='~/Dropbox/SBE/cran/bfd/cran.snp.data.genlight.Rdata')
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/bfd')

#	F4
dataset = 'F4'
nperpop = 2
missingness = 0
nloci = 100
file.name = paste0(dataset,'.nex')

#load genetic data
g = get(load('cran.snp.data.genlight.Rdata'))

#Include two curtata samples (excluding cra.cur.B43876)
g = g[!grepl('cra.cur.B43876',indNames(g)), ]


#subsample down to n individuals per population
#must do before sampling loci so that you get 100% coverage
x = data.frame(ID=indNames(g),pop=g@pop)
foo = split(x,as.factor(x$pop))
foo = lapply(foo,function(f) if(nrow(f) >= nperpop){	set.seed(12);f[sample(1:nrow(f),nperpop), ]} else f)
x = do.call(rbind,foo)
g = g[indNames(g) %in% x[,1],]


percent.missing = round(glNA(g) / (nInd(g)*2),2)
x = length(which(percent.missing <= missingness ))
cat('\n\nN loci with % missing <= 0    -    ',x,'\n\n\n',sep='')

#	% missing data <= X
g = g[, (glNA(g) / (nInd(g)*2)) <= missingness]

#subsample loci
set.seed(123)
g = g[,sample(locNames(g),nloci)]

#create nexus file to create base xml file in Beauti following instructions in 
#		http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2016/01/BFDstar-tutorial1.pdf
#pop names are the species names
#
gl2snapp(g, outfile = file.name, outpath = getwd(), v = 2)

