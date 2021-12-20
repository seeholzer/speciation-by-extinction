#Log of test runs to determine right balance of number of loci, chainlength, individuals per pop, etc.



options(scipen=999)
#BFD hypotheses
h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
hraw = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])


##################################################################
######	First Trial Run with full chainlength and 3 replicates
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/nexus.files.remote/bfd.2perpop_100loci/')
models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 3
#setup path sample variables
chainLength = 100000; burnInPercentage = 50; preBurnin = 10000; nrOfSteps = 24

#conclusions from first trial run
#	- FULL_2sp model took about 2 days to run, everything else can only be shorter than this
#	- replicates have very high consistency in ML estimates (+/- 0.1)
#		- don't need to run replicates for future trail runs


##################################################################
######	Second Trial Run with full chainlength, 1 replicate, 
#	addition of SBE2 and noZar scenarios
#	changed burnin percentage from 50 to 20 and Pre-burnin from 10000 to 50000
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/nexus.files.remote/bfd.2perpop_100loci/')
models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp","SBE2_1sp","SBE2_2sp","FULLnoZar_1sp","FULLnoZar_2sp","SBE1noZar_1sp","SBE1noZar_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions from first trial run
#	- in the FULL dataset, 2sp is consistently ranked higher than 1 sp with a decent BF of 3.47 ('positive evidence')
#		- Why are they being split?
#	- Removal of Zaratensis decreases support for two species by about 1 BF to 2.35. Pretty marginal support for splitting. SEB1noZar_2sp is ranked higher than 1sp, but not as dramatic an increase in support for 2sp as for the FULL dataset.
#		-probably should focus on FULL dataset
#	- Don't need to really need test SBE models again, pretty much know that if an intermediate population is removed BF support will dramatically increase. Goal is find a justifiable FULL dataset that supports 1 sp species.
#	-Tests
#		- Remove Amotape (or make outgroup) 
#		- Increase loci
#		- Increase N per pop


##################################################################
######	Third Trial Run with full chainlength, 1 replicate, 2perpop, 200 loci 
#	increased loci to 200
#	only include SBE1 scenario
#	allows direct comparison with Second Trial to see how increasing loci affects computation time and results
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T3/')
dataset = 'T3'
l  = readLines(paste0(dataset,'.xml'))

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions from third trial run
#	-doubling N loci doubles run times but results in higher BF support
#	- pretty much same pattern as T2 just with more BF support
#	- in the full dataset, 2 sp is ranked higher than 1 sp.


##################################################################
######	Fourth Trial Run with full chainlength, 1 replicate, 3perpop, 100 loci
#	increased N per pop to 3
#	need to include loci with some missing data (<= 0.025) to get to 100 loci
#	only include SBE1 scenario
#	allows direct comparison with Second Trial to see how increasing population sample size affects computation time and results
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T4/')
dataset = 'T4'
l  = readLines(paste0(dataset,'.xml'))

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions from four trial run
#	- didn't work at all, probably something wrong with priors
#	- also, prohibitively long runtimes so just stick with two per population, which is still around 27-29 individual per sp, so plenty of individuals


##################################################################
######	Fifth Trial Run with full chainlength, 1 replicate, 2perpop, 100 loci 
#	increased path sampling steps to 48
#	only include SBE1 scenario
#	allows direct comparison with Second Trial to see how increasing nrOfSteps affects computation time and results
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T5/')
dataset = 'T5'
l  = readLines(paste0(dataset,'.xml'))

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 48

#conclusions from same results and run times as T3 trial but with higher BF support overall with the pattern as T3
#	- likelihood curve for path sampling graphs looks pretty much the same as 28 (flattening out at the end) so probably not worth bumping up to 48 steps


##################################################################
######	Sixth Trial Run with full chainlength, 1 replicate, 2perpop, 100 loci 
#	Amotape might be adding more divergent loci to antisiensis that are not bridged by intermediates due to it's geographic isolation, disruption of isolation-by-distance patterns have been shown increase likelihood of recovering more species (Mason et al. 2020)
#	Removed Amotape
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T6/')
dataset = 'T6'
l  = readLines(paste0(dataset,'.xml'))

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- this one was interesting. 
#	- when you once again, 2sp model ranked highest in both FULL and SBE1 however, the BF support overall for lumping is much lower
#	- when amotape was removed, it decreased the support for but did not eliminate the pattern seen in T3. 
#	- The relative increase in BF support for 2sp from FULL to SBE1 was also much lower compared to T2 and especially T3


##################################################################
######	Seventh Trial Run with full chainlength, 1 replicate, 2perpop, 100 loci 
#	Removed curtata as outgroup and used Amotape instead
options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T7/')
dataset = 'T7'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
hraw = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = hraw
h = h[!(h$pop %in% 'curtata'), ]
h[h$pop %in% 'Amotape', -1] = 'amo' 
h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp","SBE2_1sp","SBE2_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- relative to T2, increased support for splitting in the Full dataset
#	- I expected it to decrease support for splitting because it made the rest of antisiensis more similar to baroni, but I guess not. 
#	- probably not worth persueing this further

##################################################################
######	Eigth Trial Run with full chainlength, 1 replicate, 2perpop, 100 loci 
#	Moved origin of local extinction one population cluster to the north
#	curtata is outgroup
options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T8/')
dataset = 'T8'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.T8.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp","SBE2_1sp","SBE2_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- by shifting the origin of simulated extinction N, this had basically the same pattern as T2
#	- moderate support for splitting in the full with progressively stronger support with extinction 


##################################################################
###### Trial Run 9 with full chainlength, 1 replicate, 2perpop, 100 loci 
#	Moved origin of local extinction one population cluster to the south
options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T9/')
dataset = 'T9'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.T9.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

models = c("FULL_1sp","FULL_2sp","SBE1_1sp","SBE1_2sp","SBE2_1sp","SBE2_2sp")
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- by shifting the origin of simulated extinction N, this had basically the same pattern as T2
#	- moderate support for splitting in the full with progressively stronger support with extinction 

#conclusions from T8 & T9
#	- shifting the arbitrary geographic boundary between antisiensis N or S yields moderate support for splitting, thus demonstrating that the boundary is arbitrary and the BFD* results are not taxonomically actionable for the full dataset, particularly because the combined phenotypic NMM analysis supports a single cluster
#	- However, as extinction increases, we recover very strong BF support for two species, this, combined with a clear distributional gap and a trend towards phenotypic clustering of samples on either side of the boundary, strongly suggests that extinction would result the recognition of two species and thus a speciation event. Proof-of-concept for speciation-by-extinction!

#to do
#	- set up runs for three different origins of local extinction with 100 loci for all scenarios
#	- rework the NMM to examine the three different origins of local extinction as well     

##################################################################
###### Trial Run 10 with full chainlength, 1 replicate, 2perpop, 100 loci 
#	Combined all three versions of the local extinction progression into a single run
#	Each versions varies in the geographic location of the divide between antisiensis and baroni and thus the origin of local extinction. 
#	v1 - taxonomic border and extinction origin N of lat midpoint
#	v2 - taxonomic border and extinction origin at lat midpoint (original version) 
#	v3 - taxonomic border and extinction origin S of lat midpoint
 
options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/T10/')
dataset = 'T10'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

#for(i in colnames(h[,-1])) cat("'",i,"',",sep='')
models = c('FULL_v1_1sp','FULL_v1_2sp','SBE1_v1_1sp','SBE1_v1_2sp','SBE2_v1_1sp','SBE2_v1_2sp','SBE3_v1_1sp','SBE3_v1_2sp','FULL_v2_1sp','FULL_v2_2sp','SBE1_v2_1sp','SBE1_v2_2sp','SBE2_v2_1sp','SBE2_v2_2sp','SBE3_v2_1sp','SBE3_v2_2sp','FULL_v3_1sp','FULL_v3_2sp','SBE1_v3_1sp','SBE1_v3_2sp','SBE2_v3_1sp','SBE2_v3_2sp','SBE3_v3_1sp','SBE3_v3_2sp')
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- worked largely as expected
#	- version 1 recovered must strong BF support for 2 species in the FULL dataset than version 2 or 3
#	- all versions showed increasing BF support for 2 species with increased extinction
#	- T10.SBE3_v2_2sp.r1 and T10.SBE3_v3_2sp.r1 ran out of memory for some reason, decreasing the number of nodes to 16 fixed the issue without greatly increasing run time
#		- decrease nodes to 16 for all SBE3 runs


##################################################################
###### Final Run 1 with full chainlength, 1 replicate, 2perpop, 200 loci 
#	Combined all three versions of the local extinction progression into a single run
#	Each versions varies in the geographic location of the divide between antisiensis and baroni and thus the origin of local extinction. 
#	v1 - taxonomic border and extinction origin N of lat midpoint
#	v2 - taxonomic border and extinction origin at lat midpoint (original version) 
#	v3 - taxonomic border and extinction origin S of lat midpoint

#	- decrease nodes to 16 for all SBE3 runs to avoid memory issues

options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/F1r1/')
dataset = 'F1r1'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

#for(i in colnames(h[,-1])) cat("'",i,"',",sep='')
models = c('FULL_v1_1sp','FULL_v1_2sp','SBE1_v1_1sp','SBE1_v1_2sp','SBE2_v1_1sp','SBE2_v1_2sp','SBE3_v1_1sp','SBE3_v1_2sp','FULL_v2_1sp','FULL_v2_2sp','SBE1_v2_1sp','SBE1_v2_2sp','SBE2_v2_1sp','SBE2_v2_2sp','SBE3_v2_1sp','SBE3_v2_2sp','FULL_v3_1sp','FULL_v3_2sp','SBE1_v3_1sp','SBE1_v3_2sp','SBE2_v3_1sp','SBE2_v3_2sp','SBE3_v3_1sp','SBE3_v3_2sp')
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 1
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24


#conclusions
#	- everything worked as expected and ran to completion
#	- took much longer than for 100 loci
#	- pattern for all three versions looks almost exactly the same as for 100 loci
#	- only difference is the BF values are about 3x greater than for 100 loci
#	- publish with the runs for 100 loci because more computational tractable but results do not differ for 200 loci


##################################################################
###### Final Runs 2 with full chainlength, 3 replicates, 2perpop, 100 loci 
#	Combined all three versions of the local extinction progression into a single run
#	Each versions varies in the geographic location of the divide between antisiensis and baroni and thus the origin of local extinction. 
#	v1 - taxonomic border and extinction origin N of lat midpoint
#	v2 - taxonomic border and extinction origin at lat midpoint (original version) 
#	v3 - taxonomic border and extinction origin S of lat midpoint

#	- decrease nodes to 16 for all SBE3 runs to avoid memory issues

options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/F2/')
dataset = 'F2'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

#for(i in colnames(h[,-1])) cat("'",i,"',",sep='')
models = c('FULL_v1_1sp','FULL_v1_2sp','SBE1_v1_1sp','SBE1_v1_2sp','SBE2_v1_1sp','SBE2_v1_2sp','SBE3_v1_1sp','SBE3_v1_2sp','FULL_v2_1sp','FULL_v2_2sp','SBE1_v2_1sp','SBE1_v2_2sp','SBE2_v2_1sp','SBE2_v2_2sp','SBE3_v2_1sp','SBE3_v2_2sp','FULL_v3_1sp','FULL_v3_2sp','SBE1_v3_1sp','SBE1_v3_2sp','SBE2_v3_1sp','SBE2_v3_2sp','SBE3_v3_1sp','SBE3_v3_2sp')
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 3
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- almost all runs worked on the first try, but 3 in SBE3 had a fatal error during the computation of marginal liklihoods
#	- just had to drop the number of cores to 8 and they worked after starting again
#	- results differed slightly from T10, same general pattern but weaker pattern of increasing support 
#	- go back to individuals and loci selected for T10, do three replicates and call it a day


##################################################################
###### Final Runs 3 with full chainlength, 3 replicates, 2perpop, 100 loci 
#	Used T10 template (individual + loci) rather than randomly selecting again
#	everything else like F2

options(scipen=999)
setwd('~/Dropbox/SBE/SBE.in.silico/bfd/runs/F3/')
dataset = 'F3'
l  = readLines(paste0(dataset,'.xml'))

h = read.delim('~/Dropbox/SBE/SBE.in.silico/cran.SBE.scenarios.txt',stringsAsFactors=F)
tmp = h[ ,-grep('pop|long|lat',colnames(h))]
tmp1sp = tmp; colnames(tmp1sp) = paste0(colnames(tmp1sp),'_1sp')
tmp1sp = data.frame(apply(tmp1sp,2,function(x) gsub('bar','ant',x)),stringsAsFactors=F)
tmp2sp = tmp; colnames(tmp2sp) = paste0(colnames(tmp2sp),'_2sp')
s = rep(1:ncol(tmp), each = 2) + (0:1) * ncol(tmp)
h = cbind(pop=h[ ,'pop'],cbind(tmp1sp,tmp2sp)[,s])

h = data.frame(apply(h,2,function(x) gsub('extinct',NA,x)))
h[,-1] = lapply(h[,-1], factor)  ## as.factor() could also be used

#for(i in colnames(h[,-1])) cat("'",i,"',",sep='')
models = c('FULL_v1_1sp','FULL_v1_2sp','SBE1_v1_1sp','SBE1_v1_2sp','SBE2_v1_1sp','SBE2_v1_2sp','SBE3_v1_1sp','SBE3_v1_2sp','FULL_v2_1sp','FULL_v2_2sp','SBE1_v2_1sp','SBE1_v2_2sp','SBE2_v2_1sp','SBE2_v2_2sp','SBE3_v2_1sp','SBE3_v2_2sp','FULL_v3_1sp','FULL_v3_2sp','SBE1_v3_1sp','SBE1_v3_2sp','SBE2_v3_1sp','SBE2_v3_2sp','SBE3_v3_1sp','SBE3_v3_2sp')
indices = grep(paste(models,collapse="|"),colnames(h))

nreps = 3
#setup path sample variables
chainLength = 500000; burnInPercentage = 20; preBurnin = 50000; nrOfSteps = 24

#conclusions
#	- almost all runs worked on the first try, but a few in SBE3 had a fatal error during the computation of marginal liklihoods
#	- just had to drop the number of cores to 8 and they worked after starting again
#	- but, didn't go back to fix everything while it was fresh in my head and 5 months later need to restart everthing








#Records of old genind editing parameters from SCRIPT1.genind2snapp.R

#	Trial 3
trial = 3
nperpop = 2
missingness = 0.01
nloci = 200
file.name = paste0('T',trial,'.nex')

#	Trial 4
trial = 4
nperpop = 3
missingness = 0.025
nloci = 100
file.name = paste0('T',trial,'.nex')

#	Trial 5
trial = 5
nperpop = 2
missingness = 0.01
nloci = 100
file.name = paste0('T',trial,'.nex')

#	Trial 6
trial = 6
nperpop = 2
missingness = 0.01
nloci = 100
file.name = paste0('T',trial,'.nex')

#	T7
trial = 7
nperpop = 2
missingness = 0.01
nloci = 100
file.name = paste0('T',trial,'.nex')

#	T8
trial = 8
nperpop = 2
missingness = 0.01
nloci = 100
file.name = paste0('T',trial,'.nex')

#	T9
trial = 9
nperpop = 2
missingness = 0.01
nloci = 100
file.name = paste0('T',trial,'.nex')


#	F1r1
trial = 'F1r1'
nperpop = 2
missingness = 0
nloci = 200
file.name = paste0(trial,'.nex')


#	F2
trial = 'F2'
nperpop = 2
missingness = 0
nloci = 100
file.name = paste0(trial,'.nex')
