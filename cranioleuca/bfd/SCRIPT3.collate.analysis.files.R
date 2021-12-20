##################################################
#	Description: BFD species delimitation with simulated extinction in Cranioleuca antisiensis
#	Script 3 - collate analysis files 
#	Goal: 	1. Move all files necessary for downstream BFD* analysis to single folder in cluster. 
#
#	Usage: in any directory on AMNH cluster (but preferably in target.directory)
#			target.directory = root directory with all runs for a given dataset
#
#	>module load R-3.4.1
#	>Rscript SCRIPT3.collate.analysis.files.R [target.directory]
#
##################################################
options(scipen=999) #suppress scientific notation

##################################################
# 1) Get target directory from commandArgs()

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

##################################################
# 2) Set up empty folder for analysis files (or clear contents of said folder if present)

dir = 'analysis.files'
if(!dir.exists(dir)){
		dir.create(dir)
	}else if(dir.exists(dir)){
		unlink(dir,recursive=T);dir.create(dir)
	}

##################################################
# 3) get file paths of target files  

f = list.files(recursive=T)
#remove files already in previously made analysis.files folders
f = f[!grepl('analysis.files',f)]


outs	= f[grep('out',f)]
logs 	= f[grep('step0.*log',f)]
logs 	= logs[!grepl('likelihood',logs)]
trees 	= f[grep('step0.*trees',f)]

##################################################
# 4) copy target files to analysis.files directory  

cat('\n\n\n\ncopying outs\n\n')
file.copy(outs,paste0(dir,'/',outs))

warnings()

cat('\n\n\n\ncopying logs\n\n')
file.copy(logs,paste0(dir,'/',sapply(strsplit(logs,'/'),'[',3)))
