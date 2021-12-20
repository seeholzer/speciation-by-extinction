##################################################
#	Description: Creates extinction scenarios used in the genetic (BFD) and phenotypic (NMM) species delimitation analyses  
#
##################################################
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#Load Population Data
popdata = read.delim('cran.popdata.txt',stringsAsFactors=F)
cur = popdata[popdata$pop %in% 'curtata', ]
popdata = popdata[!popdata$pop %in% 'curtata',c('popNumber','pop','long','lat','group')]

#create all possible unique combinations of extinction groups 2:8
#extinction groups 1 and 9 the populations of the subspecies antisiensis (ant) and baroni (bar), respectively, 
#	that are always extent 
v1 = 2:8
sets = do.call("c", lapply(seq_along(v1), function(i) combn(v1, i, FUN = list)))
#subset this list to only those combinations with consecutive numbers
#this ensures that each extinction scenario creates on a single distributional gap
sets = sets[unlist(lapply(sets,function(x) all(unique(diff(x)) == 1)))]


#Extinction scenarios with no extinction 
#	these 9 variants are only used in the bfd analysis
#	NMM only requires a single scenario with all populations extant
i = 2
for(i in 2:9){
	setName = paste0('FULL_v',i-1)
	popdata[ ,setName] = NA
	popdata[popdata$group < i,setName] = 'ant'
	popdata[popdata$group >= i,setName] = 'bar'
}

#Extinction scenarios with extinction
for(i in 1:length(sets)){
	set = sets[[i]]
	setName = paste0('EXT_v',i)
	popdata[ ,setName] = NA
	popdata[popdata$group %in% set,paste0('EXT_v',i)] = 'extinct'
	popdata[popdata$group < min(set),paste0('EXT_v',i)] = 'ant'
	popdata[popdata$group > max(set),paste0('EXT_v',i)] = 'bar'
}

cur[(ncol(cur)+1):ncol(popdata)] = 'cur'
colnames(cur) = colnames(popdata)
popdata = rbind(popdata,cur)

write.table(popdata,'cran.SBE.scenarios.txt',col.names=T,row.names=F,sep='\t',quote=F)

