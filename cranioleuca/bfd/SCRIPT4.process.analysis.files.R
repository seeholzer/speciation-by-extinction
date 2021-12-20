##################################################
#	Description: BFD species delimitation with simulated extinction in Cranioleuca antisiensis
#	Script 4 - BFD processing
#	Goal: 	1. Check run times and convergence.
#			2. Compute Bayes Factors for each 2sp vs. 1sp comparison
#
#
##################################################
#load packages
library(dplyr)

#setwd to analysis files for a given dataset/run
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#path to analysis files
path = 'bfd/runs/F5/analysis.files'


##################################################
# 1) Get paths of analysis files
outs = list.files(path,pattern='.out',full.names=T)
logs = list.files(path,pattern='.log',full.names=T)
trees = list.files(path,pattern='.trees',full.names=T)

##################################################
# 2) Check to make sure all the models ran

#full list of model names
	models = c('FULL_1sp','FULL_v1_2sp','FULL_v2_2sp','FULL_v3_2sp','FULL_v4_2sp','FULL_v5_2sp','FULL_v6_2sp','FULL_v7_2sp','FULL_v8_2sp','EXT_v1_1sp','EXT_v1_2sp','EXT_v2_1sp','EXT_v2_2sp','EXT_v3_1sp','EXT_v3_2sp','EXT_v4_1sp','EXT_v4_2sp','EXT_v5_1sp','EXT_v5_2sp','EXT_v6_1sp','EXT_v6_2sp','EXT_v7_1sp','EXT_v7_2sp','EXT_v8_1sp','EXT_v8_2sp','EXT_v9_1sp','EXT_v9_2sp','EXT_v10_1sp','EXT_v10_2sp','EXT_v11_1sp','EXT_v11_2sp','EXT_v12_1sp','EXT_v12_2sp','EXT_v13_1sp','EXT_v13_2sp','EXT_v14_1sp','EXT_v14_2sp','EXT_v15_1sp','EXT_v15_2sp','EXT_v16_1sp','EXT_v16_2sp','EXT_v17_1sp','EXT_v17_2sp','EXT_v18_1sp','EXT_v18_2sp','EXT_v19_1sp','EXT_v19_2sp','EXT_v20_1sp','EXT_v20_2sp','EXT_v21_1sp','EXT_v21_2sp','EXT_v22_1sp','EXT_v22_2sp','EXT_v23_1sp','EXT_v23_2sp','EXT_v24_1sp','EXT_v24_2sp','EXT_v25_1sp','EXT_v25_2sp','EXT_v26_1sp','EXT_v26_2sp','EXT_v27_1sp','EXT_v27_2sp','EXT_v28_1sp','EXT_v28_2sp')

#check to see if all the models have an out file
models[!models %in% gsub('F5.|.r1.out','',basename(outs))]



##################################################
# 3) Correct memory issues
#		Almost all runs will work on the first try but a few 
#			will have a fatal error during the computation of marginal likelihoods.
#		This error is related to memory issues.
#		This code will determine which runs have the memory issue.
#		Solution: Just need to drop the number of cores (ncpus) to 8 in the .job file and rerun. Things should work.
#			(make sure to change -threads flag in beast command in addition to ncpus) 
#
#		There is second error that always occurs right when trying to finish a 16th step (although this maybe any of the labelled steps).
#		Dropping the number of cores to 8 seems to work as well.
#
#
worked = c()
for(i in outs){
	o = readLines(i)
	dataset.hyp.run  = gsub(paste0(path,'/','|.out'),'',i)
	if(length(grep('marginal L estimate =',o)) > 0){
		cat(dataset.hyp.run,'\n')
		worked = c(worked,dataset.hyp.run)
	 }
}

incomplete = models[!models %in% gsub('F5.|.r1','',worked)] #both failed and runs that haven't finished
failed = outs[!grepl(paste0(incomplete,'.r1.out',collapse='|'), outs)] #failed
failed = gsub('.out','',failed)


##################################################
# 4) Load MLE, run times, and convergence
colnames = c('dataset.hyp.run','dataset','hyp','run','MLE','start','duration.hrs')
sum = data.frame(matrix(nrow=length(models),ncol=length(colnames)))
colnames(sum) = colnames
sum$hyp = models
sum$dataset.hyp.run = basename(outs)[match(models,gsub('F5.|.r1.out','',basename(outs)))]
sum$dataset = sapply(strsplit(sum$dataset.hyp.run,'\\.'),'[',1)
sum$run = gsub('r','',sapply(strsplit(sum$dataset.hyp.run,'\\.'),'[',3))


for(i in 1:nrow(sum)){
	if(!any(grepl(sum$hyp[i],outs))) next
	o = readLines(outs[grepl(sum$hyp[i],outs)])
	MLE = as.numeric(gsub('marginal L estimate = ','', o[grep('marginal L estimate',o)]))
	
	if(length(MLE) == 0 ){
		sum[i,'MLE'] = NA
	}else{
		sum[i,'MLE'] = MLE
		}
	 
	sum[i,'start'] = o[4]
	foo = gsub('real\t','',o[grep('real',o)[2]])
	sum[i,'duration.hrs'] = as.numeric(gsub('m(.*)','',foo))/60
}

#restart jobs that failed because "beast.xml.state (No such file or directory)"
# x = sum$dataset.hyp.run[is.na(sum$MLE)][!sum$dataset.hyp.run[is.na(sum$MLE)] %in% failed]
# for(i in x) cat(paste0('rm -r ',i),sep='\n')
# for(i in paste0(x,'.out')) cat(paste0('rm ',i),sep='\n')
# for(i in paste0(x,'.job')) cat(paste0('qsub ',i),sep='\n')

##################################################
# 6) Check Convergence
#check for convergence by getting range of MLE values, should all be very small
colnames = c('hyp','MLE.range')
foo = data.frame(matrix(ncol=2,nrow=length(unique(sum$hyp))))
colnames(foo) = colnames
foo$hyp = unique(sum$hyp)
for(i in 1:nrow(foo)){
	tmp = sum[sum$hyp %in% foo[i,'hyp'],'MLE']
	foo[i,'MLE.range'] = max(tmp) - min(tmp)
}



##################################################
# 6) Compute Bayes Factors for each scenario (i.e. 2sp vs. 1sp comparison) 

#create list of scenario names
fullModels = models[grep('FULL_v',models)] #only 2 species models because will compare each of these to the FULL_1sp
extModels = models[grep('EXT',models)]
scenarioList = unique(c(gsub('_\\dsp','',fullModels),gsub('_\\dsp','',extModels)))

#compile MLEs for 1sp and 2sp models
colnames = c('scenario','MLE_1sp','MLE_2sp','BF')
foo = data.frame(matrix(nrow=length(scenarioList),ncol=length(colnames)))
colnames(foo) = colnames
foo$scenario = scenarioList

for(i in scenarioList){
	rowIndex = which(foo$scenario == i) 

	#one species hypothesis (MLE_1sp)
	#all FULL 2sp models are compared to FULL_1sp
	if(grepl('FULL',i)) 	hyp = 'FULL_1sp' else hyp = paste0(i,'_1sp')
	foo[rowIndex,'MLE_1sp'] = sum[sum$hyp %in% hyp,'MLE']

	#two species hypothesis (MLE_1sp)
	hyp = paste0(i,'_2sp')
	if(any(grepl(hyp,sum$hyp))){
		foo[rowIndex,'MLE_2sp'] = sum[sum$hyp %in% hyp,'MLE']	
	}
}


##################################################
# 7) Calculate Bayes Factors

# From BFD* tutorial
# BF = 2 x (MLE1 - MLE0)
#	MLE0 is the reference model, in all cases, 1 species
# 	A positive BF value indicates support in favor of model 1 (MLE1).
#	A negative BF value indicates support in favor of model 0 (MLE0).

# The BF scale is as follows:
# 0 < BF < 2 is not worth more than a bare mention
# 2 < BF < 6 is positive evidence
# 6 < BF < 10 is strong support
# BF > 10 is decisive.

foo$BF = 2*(foo$MLE_2sp - foo$MLE_1sp)

write.table(foo,'results.bfd.txt',col.names=T,row.names=F,quote=F,sep='\t')