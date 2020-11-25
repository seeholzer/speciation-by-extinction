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
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
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

outs	= f[grep('out',f)]
logs 	= f[grep('step0.*log',f)]
logs 	= logs[!grepl('likelihood',logs)]
trees 	= f[grep('step0.*trees',f)]

##################################################
# 4) copy target files to analysis.files directory  

file.copy(outs,paste0(dir,'/',outs))
file.copy(logs,paste0(dir,'/',sapply(strsplit(logs,'/'),'[',3)))
file.copy(trees,paste0(dir,'/',sapply(strsplit(trees,'/'),'[',3)))

##################################################
# 5) run treeannotator on tree files in analysis.files directory  

setwd(dir)

f = list.files(pattern='trees')

for(t in f){
	I = t
	O = gsub('trees','consensus.tree',I)
	system(paste('treeannotator -burnin 20 -heights median',I,O))
}


