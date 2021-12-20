##################################################
#	Description: BFD species delimitation with simulated extinction in Cranioleuca antisiensis
#	Script 4b - batch restart of jobs that ran into memory issues
#
#
##################################################
#load packages
library(dplyr)

#setwd to analysis files for a given dataset/run
root = '~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/'

#path to analysis files
path = 'bfd/runs/F1r1/analysis.files'
path = 'bfd/runs/F2/analysis.files'
path = 'bfd/runs/F3/analysis.files'
path = 'bfd/runs/F4/analysis.files'
path = 'bfd/runs/F5/analysis.files'
path = 'bfd/runs/F5/analysis.files2'
path = 'bfd/runs/F5/analysis.files3'

setwd(paste0(root,path))

##################################################
# 1) Get paths of analysis files
outs = list.files(pattern='.out')
logs = list.files(pattern='.log')
trees = list.files(pattern='.trees')


##################################################
# 2) Check to make sure all the models ran

#full list of model names
	models = c('FULL_1sp','FULL_v1_2sp','FULL_v2_2sp','FULL_v3_2sp','FULL_v4_2sp','FULL_v5_2sp','FULL_v6_2sp','FULL_v7_2sp','FULL_v8_2sp','EXT_v1_1sp','EXT_v1_2sp','EXT_v2_1sp','EXT_v2_2sp','EXT_v3_1sp','EXT_v3_2sp','EXT_v4_1sp','EXT_v4_2sp','EXT_v5_1sp','EXT_v5_2sp','EXT_v6_1sp','EXT_v6_2sp','EXT_v7_1sp','EXT_v7_2sp','EXT_v8_1sp','EXT_v8_2sp','EXT_v9_1sp','EXT_v9_2sp','EXT_v10_1sp','EXT_v10_2sp','EXT_v11_1sp','EXT_v11_2sp','EXT_v12_1sp','EXT_v12_2sp','EXT_v13_1sp','EXT_v13_2sp','EXT_v14_1sp','EXT_v14_2sp','EXT_v15_1sp','EXT_v15_2sp','EXT_v16_1sp','EXT_v16_2sp','EXT_v17_1sp','EXT_v17_2sp','EXT_v18_1sp','EXT_v18_2sp','EXT_v19_1sp','EXT_v19_2sp','EXT_v20_1sp','EXT_v20_2sp','EXT_v21_1sp','EXT_v21_2sp','EXT_v22_1sp','EXT_v22_2sp','EXT_v23_1sp','EXT_v23_2sp','EXT_v24_1sp','EXT_v24_2sp','EXT_v25_1sp','EXT_v25_2sp','EXT_v26_1sp','EXT_v26_2sp','EXT_v27_1sp','EXT_v27_2sp','EXT_v28_1sp','EXT_v28_2sp')

#check to see if all the models have an out file
models[!models %in% gsub('F5.|.r1.out','',basename(outs))]


##################################################
# 3) Determine which runs failed due to memory issues
#		Almost all runs will work on the first try but a few (often in SBE3) 
#			will have a fatal error during the computation of marginal liklihoods.
#		This error is related to memory issues.
#		This code will determine which runs have the memory issue.
#		Solution: Just need to drop the number of cores (ncpus) from 16 to 8 in the .job file and rerun. Things should work.
#			(make sure to change -threads flag in beast command in addition to ncpus) 
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
failed = outs[grep(paste0(incomplete,'.r1.out',collapse='|'), outs)] #failed
failed = gsub('.out','',failed)

setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/bfd/runs/F5/')
files = list.files('analysis.files3',full.names=T)
from = files[grep(paste0(failed,'.r1.out',collapse='|'),files)]
to = gsub('analysis.files3/','failed.jobs/',from)
file.copy(from,to)


#commands to remove files associated with the failed runs
for(i in failed) cat(paste0('rm -r ',i),sep='\n')
for(i in paste0(failed,'.out')) cat(paste0('rm ',i),sep='\n')



newpath = gsub('analysis.files3','',getwd())
setwd(newpath)

dataset.hyp.run = failed[1]
qsubCommands = c()
for(dataset.hyp.run in failed){
	
	lines = readLines(paste0(newpath,'/',dataset.hyp.run,'.job'))		
	#edit qsub command
	lines[grep('#PBS -l select=1:ncpus=16',lines)] = '#PBS -l select=1:ncpus=8'
	#edit -threads command
	beastIndex = grep('time beast -threads',lines)
	lines[beastIndex] = gsub('-threads 16','-threads 8',lines[beastIndex])
	writeLines(lines,paste0(newpath,dataset.hyp.run,'.job'))

	#write command to delete existing files for failed dataset.hyp.run
	#combine with qsub command for batch submission
	qsub = paste0('qsub ',dataset.hyp.run,'.job')
	qsubCommands = c(qsubCommands,qsub)
}


#copy the new jobs to another folder so easy to transfer
files = list.files()
newJobs = paste0(failed,'.job')
newXML = paste0(failed,'.xml')


newDir = 'resubmitFailedJobs'
dir.create(newDir)
from = files[grep(paste(c(newJobs,newXML),collapse='|'),files)]
to   = file.path(newDir,files[grep(paste(c(newJobs,newXML),collapse='|'),files)])
file.copy(from,to)

writeLines(qsubCommands,'resubmitFailedJobs/resubmitFailedJobs.sh')



