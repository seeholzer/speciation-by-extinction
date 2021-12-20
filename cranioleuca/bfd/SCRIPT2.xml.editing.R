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
setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/bfd/runs')

#Reformat hypotheses for code used for XML editing
	h = read.delim('../../cran.SBE.scenarios.txt',stringsAsFactors=F)
	tmp = h[ ,-grep('popNumber|pop|long|lat|group',colnames(h))]
	#create 1 species hypotheses
		tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
		tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
		#remove redundant FULL scenarios for single species
		tmp1sp = data.frame(t(unique(t(tmp1sp))),stringsAsFactors=F)
		tmp1sp.FULL = tmp1sp[grep('FULL',colnames(tmp1sp))]; colnames(tmp1sp.FULL) = 'FULL_1sp'
		tmp1sp.EXT = tmp1sp[grep('EXT',colnames(tmp1sp))]		
	#create 2 species hypotheses
		tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
		tmp2sp.FULL = tmp2sp[grep('FULL',colnames(tmp2sp))]
		tmp2sp.EXT = tmp2sp[grep('EXT',colnames(tmp2sp))]		
	#interleave 1 and 2 species EXT hypotheses
		s = rep(1:ncol(tmp1sp.EXT), each = 2) + (0:1) * ncol(tmp1sp.EXT) #creates interleave index
		EXT = cbind(tmp1sp.EXT,tmp2sp.EXT)[,s] #interleave hypotheses
		FULL = cbind(tmp1sp.FULL,tmp2sp.FULL) #combines 1 and 2 species hypotheses with no extinction
		h = cbind(pop=h[ ,'pop'],cbind(EXT,FULL))  #combine and add population data

	#minor reformating
		h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x))) #change 'extinct' to NA
		h[,-1] = lapply(h[,-1], factor)  #Make all hypotheses factors # as.factor() could also be used
		

#Set model names and order
#for(i in colnames(h[,-1])) cat("'",i,"',",sep='') #formatted model names
	models = c('FULL_1sp','FULL_v1_2sp','FULL_v2_2sp','FULL_v3_2sp','FULL_v4_2sp','FULL_v5_2sp','FULL_v6_2sp','FULL_v7_2sp','FULL_v8_2sp','EXT_v1_1sp','EXT_v1_2sp','EXT_v2_1sp','EXT_v2_2sp','EXT_v3_1sp','EXT_v3_2sp','EXT_v4_1sp','EXT_v4_2sp','EXT_v5_1sp','EXT_v5_2sp','EXT_v6_1sp','EXT_v6_2sp','EXT_v7_1sp','EXT_v7_2sp','EXT_v8_1sp','EXT_v8_2sp','EXT_v9_1sp','EXT_v9_2sp','EXT_v10_1sp','EXT_v10_2sp','EXT_v11_1sp','EXT_v11_2sp','EXT_v12_1sp','EXT_v12_2sp','EXT_v13_1sp','EXT_v13_2sp','EXT_v14_1sp','EXT_v14_2sp','EXT_v15_1sp','EXT_v15_2sp','EXT_v16_1sp','EXT_v16_2sp','EXT_v17_1sp','EXT_v17_2sp','EXT_v18_1sp','EXT_v18_2sp','EXT_v19_1sp','EXT_v19_2sp','EXT_v20_1sp','EXT_v20_2sp','EXT_v21_1sp','EXT_v21_2sp','EXT_v22_1sp','EXT_v22_2sp','EXT_v23_1sp','EXT_v23_2sp','EXT_v24_1sp','EXT_v24_2sp','EXT_v25_1sp','EXT_v25_2sp','EXT_v26_1sp','EXT_v26_2sp','EXT_v27_1sp','EXT_v27_2sp','EXT_v28_1sp','EXT_v28_2sp')
#Get column indices in h for each model
	indices = sapply(models,function(x)grep(x,colnames(h)))


##################################################
# 2) Set up run parameters
#	-full chainlength, 3 replicates, 2perpop, 100 loci 

#name of dataset
	dataset = 'F5'
	l  = readLines(paste0(dataset,'/',dataset,'.xml'))
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
#	'	rootdir="/home/gseeholzer/nas3/bfd/runs/dataset/dataset.hyp.run"',
	'	rootdir="/home/lmoreira/nas3/bfd/runs/dataset/dataset.hyp.run"',
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

i = indices['FULL_v8_2sp']
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
		
		job = readLines('bfd.template.job')
	
		job = gsub('dataset.hyp.run',dataset.hyp.run,job)
		job = gsub('dataset.hyp',dataset.hyp,job)
		job = gsub('dataset',dataset,job)
		
		if(grepl('SBE3',dataset.hyp.run)){
			job = gsub('ncpus=32','ncpus=16',job)
			job = gsub('-threads 32','-threads 16',job)
		}
			
		#writeLines(job,paste0(dataset,'/',dataset.hyp.run,'.job'))
		
		path = l.hyp[grep('rootdir',l.hyp)]
		new.path = path
		new.path = gsub('dataset.hyp.run',dataset.hyp.run, new.path)
		new.path = gsub('dataset',dataset, new.path)
	
		l.hyp.run = l.hyp
		l.hyp.run = gsub(path,new.path,l.hyp.run)
		l.hyp.run = gsub('snap.log',paste0(dataset.hyp.run,'.log'),l.hyp.run)	
		l.hyp.run = gsub('snap.trees',paste0(dataset.hyp.run,'.trees'),l.hyp.run)	
		
		#writeLines(l.hyp.run,paste0(dataset,'/',dataset.hyp.run,'.xml'))
			
		qsub 	= paste0("qsub ", dataset.hyp.run,'.job')
		commands = c(commands,	qsub)	
		}
		
}

##################################################
# 4) reorder commands so longest runs start first
# FULL_2sp, SBE1_2sp, FULL_1sp, SBE1_1sp, SBE2_2sp, SBE2_1sp, SBE3_2sp, SBE3_1sp, 

scenarios = c('FULL',paste0('EXT_v',1:28))
popdata = read.delim('../../cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = popdata[,grep(paste(scenarios,collapse='|'),colnames(popdata))]
x = tmp[1]
nExtinctions = apply(tmp,2,function(x){
	foo = data.frame(group=popdata$group,scenario=x)
	foo = unique(foo)
	nExtinctions = length(which(foo$scenario %in% 'extinct'))
	return(nExtinctions)
})
nExtinctions = data.frame(scenario.name=names(nExtinctions),nExtinctions=nExtinctions)
#add a line for Full_1sp
nExtinctions = rbind(data.frame(scenario.name='FULL',nExtinctions = 0),nExtinctions)


#reorder to get the ones with the most data started first
order = c()
for(i in 1:nrow(nExtinctions)){
	order = c(order,grep(nExtinctions[i,'scenario.name'],commands))
}

writeLines(commands[order],paste0(dataset,'/qsub.commands.sh'))


#To batch submit
#	chmod u+x qsub.commands.sh
#	./qsub.commands.sh




#subset to only run1 to start 
commands = readLines(paste0(dataset,'/qsub.commands.sh'))
command.subset = commands[grep('r1',commands)]
writeLines(command.subset,paste0(dataset,'/qsub.commands.r1.sh'))

#	chmod u+x qsub.commands.r1.sh
#	./qsub.commands.r1.sh
