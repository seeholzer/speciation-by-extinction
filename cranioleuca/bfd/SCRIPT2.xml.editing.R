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
options(scipen=999) #suppress scientific notation

##################################################
# 1) Set up hypotheses 
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/bfd')

#Reformat hypotheses for code used for XML editing
	h = read.delim('cran.SBE.scenarios.txt',stringsAsFactors=F)
	tmp = h[ ,-grep('pop|long|lat',colnames(h))]
	#create 1 species hypotheses
		tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
		tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
	#create 2 species hypotheses
		tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
	#interleave 1 and 2 species hypotheses
		s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp) #creates interleave index
		h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s]) #interleave hypotheses and add population
	#minor reformating
		h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x))) #change 'extinct' to NA
		h[,-1] = lapply(h[,-1], factor)  #Make all hypotheses factors # as.factor() could also be used

#Set model names and order
#for(i in colnames(h[,-1])) cat("'",i,"',",sep='') #formatted model names
	models = c('FULL_v1_1sp','FULL_v1_2sp','SBE1_v1_1sp','SBE1_v1_2sp','SBE2_v1_1sp','SBE2_v1_2sp','SBE3_v1_1sp','SBE3_v1_2sp','FULL_v2_1sp','FULL_v2_2sp','SBE1_v2_1sp','SBE1_v2_2sp','SBE2_v2_1sp','SBE2_v2_2sp','SBE3_v2_1sp','SBE3_v2_2sp','FULL_v3_1sp','FULL_v3_2sp','SBE1_v3_1sp','SBE1_v3_2sp','SBE2_v3_1sp','SBE2_v3_2sp','SBE3_v3_1sp','SBE3_v3_2sp')
#Get column indices in h for each model
	indices = grep(paste(models,collapse="|"),colnames(h))


##################################################
# 2) Set up run parameters
#	-full chainlength, 3 replicates, 2perpop, 100 loci 

#name of dataset
	dataset = 'F4'
	setwd(paste0('~/Dropbox/SBE/SBE.in.silico/bfd/runs/',dataset))
	l  = readLines(paste0(dataset,'.xml'))
#setup path sampling variables
	chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24
#number of replicates for each model to assess convergence
	nreps = 3


##################################################
# 3) Add path sampler block to XML template

	block = c(
	'<run spec="beast.inference.PathSampler"',
	'	chainLength="1000"',
	'	alpha="0.3"',
	'	rootdir="/home/gseeholzer/nas3/bfd/runs/dataset/dataset.hyp.run"',
	'	burnInPercentage="0"',
	'	preBurnin="0"',
	'	deleteOldLogs="true"',
	'	nrOfSteps="24">',
	'	cd $(dir)',
	'	java -cp $(java.class.path) beast.app.beastapp.BeastMain $(resume/overwrite) -java -seed $(seed) beast.xml'
	)
	
	#replace path sampling variables in path sampling block
	block = gsub('chainLength="[0-9]+"',paste0('chainLength="',chainLength,'"'),block)
	block = gsub('burnInPercentage="[0-9]+"',paste0('burnInPercentage="',burnInPercentage,'"'),block)
	block = gsub('preBurnin="[0-9]+"',paste0('preBurnin="',preBurnin,'"'),block)
	block = gsub('nrOfSteps="[0-9]+"',paste0('nrOfSteps="',nrOfSteps,'"'),block)
	block = gsub('nrOfSteps="[0-9]+"',paste0('nrOfSteps="',nrOfSteps,'"'),block)
	
	#insert path sampling block in to template
	l = c(l[1:(grep('<run',l)-2)],block,l[(grep('<run',l)-1):length(l)])
	
	#replace <run with <mcmc
	l[grep('<run id="mcmc"',l)] = gsub('<run','<mcmc',l[grep('<run id="mcmc"',l)])
	
	#add terminal </mcmc>
	l[(grep('</run>',l)-1)] = "    </mcmc>"


#		DO NOT NEED TO SET UP MCMC CHAIN LENGTH
#	MCMC variables over written with the path sampling MCMC variables
#	when the path sampling program generates each path step


##################################################
# 3) Edit XML template for each specific scenario and write 
#		[dataset].[scenario]_[replicate]_[Nsp].xml
#		[dataset].[scenario]_[replicate]_[Nsp].job

 
start.taxonset	= l[grep('<taxonset',l)[1]] 	#start line index of taxon set
end.taxonset 	= l[grep('</taxonset>',l)[1]] 	#end line index of taxon set
taxa = l[grep('<taxon id=',l)]					#taxon sets

commands = c()									#create empty commands vector

for(i in indices){
	cat(i,'\n')
	hyp = colnames(h)[i]
	tmp = h[,c('pop',hyp)]

	dataset.hyp = paste0(dataset,'.',hyp)

	#create taxon block for one species + outgroup
	if(length(levels(tmp[,2])) == 2){ 
		
		taxon = levels(tmp[,2])[1]
		taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		taxonblock1 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
		
		taxon = levels(tmp[,2])[2]
		taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		taxonblock2 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
		
		l.hyp = c(
		l[1:grep('<rawdata',l)],				#pre-taxonblock
		taxonblock1,							#taxonblock1
		taxonblock2,							#taxonblock2 - outgroup
		l[grep('</taxa>',l):length(l)]) 		#pos-taxonblock	

	}

	#create taxon block for two species + outgroup
	if(length(levels(tmp[,2])) == 3){ 
		
		taxon = levels(tmp[,2])[1]
		taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		taxonblock1 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
		
		taxon = levels(tmp[,2])[2]
		taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		taxonblock2 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)

		taxon = levels(tmp[,2])[3]
		taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		taxonblock3 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
	
		l.hyp = c(
		l[1:grep('<rawdata',l)],				#pre-taxonblock
		taxonblock1,							#taxonblock1
		taxonblock2,							#taxonblock2
		taxonblock3,							#taxonblock3 - outgroup
		l[grep('</taxa>',l):length(l)]) 		#pos-taxonblock	

	}
	
	#DO NOT RUN - for possible future use
	#create taxon block for three species + outgroup
	# if(length(levels(tmp[,2])) == 4){ 
				
		# taxon = levels(tmp[,2])[1]
		# taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		# taxonblock1 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
		
		# taxon = levels(tmp[,2])[2]
		# taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		# taxonblock2 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)

		# taxon = levels(tmp[,2])[3]
		# taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		# taxonblock3 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
	
		# taxon = levels(tmp[,2])[4]
		# taxonset = taxa[grep(paste(tmp[!is.na(tmp[,2]) & tmp[,2] %in% taxon,'pop'],collapse='|'),taxa)]
		# taxonblock4 = c(gsub('id=\"[a-zA-Z]+\"',paste0('id=\"',taxon,'\"'),start.taxonset), taxonset, end.taxonset)
	
		# l.hyp = c(
		# l[1:grep('<rawdata',l)],		#pre-taxonblock
		# taxonblock1,						#taxonblock1
		# taxonblock2,						#taxonblock2
		# taxonblock3,						#taxonblock3
		# taxonblock4,						#taxonblock3
		# l[grep('</taxa>',l):length(l)]) #pos-taxonblock	

	# }

	#create three replicates per scenario to assess convergence
	rep = 1
	for(rep in 1:nreps){
		dataset.hyp.run = paste0(dataset.hyp,'.r',rep)
		
		job = readLines('~/Dropbox/SBE/SBE.in.silico/bfd/runs/bfd.template.job')
	
		job = gsub('dataset.hyp.run',dataset.hyp.run,job)
		job = gsub('dataset.hyp',dataset.hyp,job)
		job = gsub('dataset',dataset,job)
		
		if(grepl('SBE3',dataset.hyp.run)){
			job = gsub('ncpus=32','ncpus=16',job)
			job = gsub('-threads 32','-threads 16',job)
		}
			
		writeLines(job,paste0(dataset.hyp.run,'.job'))
		
		path = l.hyp[grep('rootdir',l.hyp)]
		new.path = path
		new.path = gsub('dataset.hyp.run',dataset.hyp.run, new.path)
		new.path = gsub('dataset',dataset, new.path)
	
		l.hyp.run = l.hyp
		l.hyp.run = gsub(path,new.path,l.hyp.run)
		l.hyp.run = gsub('snap.log',paste0(dataset.hyp.run,'.log'),l.hyp.run)	
		l.hyp.run = gsub('snap.trees',paste0(dataset.hyp.run,'.trees'),l.hyp.run)	
		
		writeLines(l.hyp.run,paste0(dataset.hyp.run,'.xml'))
			
		qsub 	= paste0("qsub ", dataset.hyp.run,'.job')
		commands = c(commands,	qsub)	
		}
		
}

##################################################
# 4) reorder commands so longest runs start first
# FULL_2sp, SBE1_2sp, FULL_1sp, SBE1_1sp, SBE2_2sp, SBE2_1sp, SBE3_2sp, SBE3_1sp, 

order = c(grep('FULL.+2sp',commands),
			grep('SBE1.+2sp',commands),
			grep('FULL.+1sp',commands),
			grep('SBE1.+1sp',commands),
			grep('SBE2.+2sp',commands),
			grep('SBE2.+1sp',commands),
			grep('SBE3.+2sp',commands),
			grep('SBE3.+1sp',commands))

writeLines(commands[order],'qsub.commands.sh')
#To batch submit
#	chmod u+x qsub.commands.sh
#	./qsub.commands.sh

