##################################################
#	Description: BFD species delimitation with simulated extinction in Cranioleuca antisiensis
#	Script 2 - XML editing
#	Goal: 	1. Create Beast .xml file for each replicate of each scenario for a given input genetic dataset.
#			2. Create PBS .job files for each Beast .xml file to run on AMNH cluster
#			3. Create a bash script to batch submit all .job files to AMNH cluster
#
#	- xml editing follows instructions in 
#		http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2016/01/BFDstar-tutorial1.pdf
#	- a template path sampling xml file must be first generate in Beauti following tutorial instructions
#	- a template PBS job must also be present
#	- both are named [dataset].xml or [dataset].job
#
##################################################
#load packages
library(dplyr)

#setwd to analysis files for a given dataset/run
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/')

#path to analysis files
path = 'bfd.old/runs/F1r1/analysis.files'
path = 'bfd.old/runs/F2/analysis.files'
path = 'bfd.old/runs/F3/analysis.files'
path = 'bfd.old/runs/F4/analysis.files'

##################################################
# 1) Get paths of analysis files
outs = list.files(path,pattern='.out',full.names=T)
logs = list.files(path,pattern='.log',full.names=T)
trees = list.files(path,pattern='.trees',full.names=T)

##################################################
# 2) Correct memory issues
#		Almost all runs will work on the first try but a few (often in SBE3) 
#			will have a fatal error during the computation of marginal liklihoods.
#		This error is related to memory issues.
#		This code will determine which runs have the memory issue.
#		Solution: Just need to drop the number of cores (ncpus) from 16 to 8 in the .job file and rerun. Things should work.
#			(make sure to change -threads flag in beast command in addition to ncpus) 
i = outs[1]
for(i in outs){
	o = readLines(i)
	dataset.hyp.run  = gsub(paste0(path,'/','|.out'),'',i)
	if(length(grep('Fatal exception',o)) > 0) cat(dataset.hyp.run,'\n')
}

##################################################
# 2) Run times and convergence
colnames = c('dataset.hyp.run','dataset','hyp','run','MLE','start','duration.hrs')
sum = data.frame(matrix(nrow=length(outs),ncol=length(colnames)))
colnames(sum) = colnames
sum$dataset.hyp.run = gsub(paste0(path,'/','|.out'),'',outs)
sum$dataset = sapply(strsplit(sum$dataset.hyp.run,'\\.'),'[',1)
sum$hyp = sapply(strsplit(sum$dataset.hyp.run,'\\.'),'[',2)
sum$run = gsub('r','',sapply(strsplit(sum$dataset.hyp.run,'\\.'),'[',3))

i = 1
for(i in 1:nrow(sum)){
	dataset.hyp.run = sum$dataset.hyp.run[i] 
	o = readLines(outs[grep(dataset.hyp.run,outs)])
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
# 3) Compute Bayes Factors for each 2sp vs. 1sp comparison

#empty list for each geographical variant (i.e. V1, V2, V3) of the nested extinction scenarios 
versions = vector(mode = "list", length = length(c('v1','v2','v3')))	
names(versions) = c('v1','v2','v3')

#loop through each geogeographic variant (version)
v = 1
for(v in 1:length(versions)){
	version = names(versions)[v]

	#empty list for each extinction scenario
	model.set.names = c('FULL','SBE1','SBE2','SBE3')
	model.sets = vector("list", length(model.set.names))				
	names(model.sets) = model.set.names
	
	s = 'FULL'
	for(s in model.set.names){
		cat(s,'_',version,'\n',sep='')
		
		#populate dataframe with maximum-likelihood estimates (MLE) for 
		#comparions of 1 vs. 2 species models for a given extinction scenario
		Nsp = 2		#1 vs. 2 species models
		#colnames = c('Model','MLE','Rank','BF')
		colnames = c('Model','MLE')
		d = data.frame(matrix(nrow=Nsp,ncol=length(colnames)))
		colnames(d) = colnames
		d$Model = paste0(s,'_',version,'_',paste0(1:Nsp,'sp'))
		
		i = 1
		for(i in 1:nrow(d)){
			m = d[i,'Model']
			x = sum[sum$hyp %in% m, ]
			#r = x[which.max(x$MLE),'run']		#select max MLE estimate across all runs DON'T USE
			r = 1								#only select MLE estimate from run 1			
			o = readLines(outs[grep(paste0(m,'.r',r),outs)])
			marg.L.line = o[grep('marginal L estimate',o)]
			marg.L = as.numeric(gsub('marginal L estimate = ','', marg.L.line))
			d[i,'MLE'] = marg.L
		}
		
		model.sets[[s]] = d		
				
	}	
	versions[[v]] = model.sets
}	
versions



##################################################
# 4) Reformat list of model MLEs so that each extinction scenario has it's own row

#collate list of dataframes into single dataframe
tmp = bind_rows(bind_rows(versions[[1]]),bind_rows(versions[[2]]),bind_rows(versions[[3]]))

colnames = c('Dataset','Version','MLE.1sp','MLE.2sp','BF')
sum = data.frame(matrix(nrow = length(grep('1sp',tmp$Model)), ncol=length(colnames)))
colnames(sum) = colnames

sum[,c('Dataset','Version')] = unique(t(sapply(strsplit(tmp$Model,'_'),'[',1:2)))

for(i in 1:nrow(sum)){
	one = paste0(sum[i,'Dataset'],'_',sum[i,'Version'],'_1sp')
	two = paste0(sum[i,'Dataset'],'_',sum[i,'Version'],'_2sp')
	
	sum[i,'MLE.1sp'] = tmp[grep(one,tmp$Model),'MLE']
	sum[i,'MLE.2sp'] = tmp[grep(two,tmp$Model),'MLE']
	
}

##################################################
# 5) Calculate Bayes Factors

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

sum$BF = 2*(sum$MLE.2sp - sum$MLE.1sp)

##################################################
# 6) Write results

write.table(sum,'results.bfd.txt',col.names=T,row.names=F,quote=F,sep='\t')









##################################################
# 7) Quickly plot results to check that things look right
library(RColorBrewer)

key = data.frame(trait.set=c('v1','v2','v3'),col=brewer.pal(3,'Set2'),pch=c(15,16,17),stringsAsFactors=F)

dev.new(width=5,height=5)
par(mar=c(3,5,0,1))
xlim = c(0.5,4.5); ylim = c(-50,150)
if(min(sum$BF) < ylim[1]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')
if(max(sum$BF) > ylim[2]) cat('\n\n\n\n\n\n ylim out of bounds \n\n\n\n\n\n ')

plot(sum[ ,'BF'], xlim=xlim, ylim=ylim, xlab='', ylab='', type='n', axes=F)
abline(h=0,lty=2)
axis(1,at=c(1,2,3,4),labels=c('FULL','SBE1','SBE2','SBE3'),tick=T,lty=NULL,col = NA, col.ticks = 1)
axis(2,at=c(-50,-25,0,25,50,70),labels=c(-50,-25,0,25,50,70))
mtext('Bayes Factor',2,line=3.5,cex=1.25)
mtext('One Species',2,at=-25,line=2,col='#E41A1C')
mtext('Two Species',2,at= 75/2,line=2,col='#377EB8')
i = 'v1'
for(i in unique(sum$Version)){
	col = key[key[,1] %in% i,'col']	
	pch = key[key[,1] %in% i,'pch']	
	lines(1:4,sum[sum$Version %in% i,'BF'],col=col)
	points(1:4,sum[sum$Version %in% i,'BF'],col=col,pch=pch,cex=2)
}

legend(0.5, 150, legend=key[,1], col=key[,2], pch=key[,3], lty=1, cex=1, box.lty=0)









