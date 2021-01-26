##################################################
#	Description: Species delimitation of Darwin's Finch morphological dataset with simulated extinction using normal mixture models
#	Script 2 - mclust model processing, identify cases of SBE, create Figure S3
#	Author: Glenn F. Seeholzer
#	Update history at github.com/seeholzer/speciation-by-extinction
##################################################

#	***Set working directory to local path***
setwd('~/Dropbox/SBE/speciation-by-extinction_github/geospiza')

#load packages
library(plyr)

##################################################
# 2.1) Read data
data = read.table("Geospiza.data.csv", header=T, sep=",")
data$ID = paste0(data$Institution,data$Museum.Number)
colnames(data)[1] = 'Taxa.Lack'
colnames(data)[6] = 'Taxa.New'
colnames.to.keep = c('ID','Taxa.Lack','Taxa.New','Wing','Tail','Blength','Bdepth','Bwidth','Tarsus')
data = data[,colnames.to.keep]

#load mclust results
mclust.data = get(load('mclust.Bdepth.SML.groundfinches.rda'))
plots.dir = 'plots.Bdepth.SML.groundfinches'

#Data for all Geospiza species (LEGACY)
# mclust.data = get(load('mclust.Bdepth.rda')) 
# plots.dir = 'plots.Bdepth.all'


##################################################
# 2.2) Merge all mclust classifications and other relevant data

#create list of dataframes with the mclust classification of each extant specimen
l = lapply(mclust.data,function(x){
	mclust = x$mclust
	extant = x$extant
	tmp = data.frame(ID = extant, class = mclust$classification)
	return(tmp)	
	})

#rename class column as the extinction scenario
for(i in 1:length(l)) colnames(l[[i]])[2] = names(l)[i]
#join all the dataframes where NA means the specimens was removed in the simulated extinction event 
foo = join_all(l,by='ID',type='full')

d = merge(data[,c('ID','Taxa.Lack','Taxa.New')],foo,by='ID',all=T)

#change class columns to factors with levels equal to classes
for(i in 4:ncol(d)){
	x = d[,i]
	levels = sort(as.numeric(unique(x[!is.na(x)])))
	d[,i] = factor(x,levels=levels)	
}

table(d[,c('no.extinction','7.00-8.25')])

##################################################
# 2.3) Compare extinction scenario classifications to the classifications with no extinction
#			- determine which extinction scenarios resulted in speciation-by-extinction
scenarios = names(mclust.data)

list = rep(list(NA),length(scenarios))
names(list) = scenarios

i = grep('6.25-9.25',scenarios)
for(i in 1:length(scenarios)){
	
	tmp = d[,c('ID','no.extinction',scenarios[i])]
	#table of number of specimens assigned to all pairwise combinations of pre- and post-extinction morphospecies
	#rows = pre-extinction morphspecies PreExt.MS
	#columns = post-extinction morphospecies PostExt.MS
	t = table(tmp[,c(2,3)]) 
	
	#Splits and Lumps
	#split 
	#	more than one PostExt.MS for a given PreExt.MS
	#lump 
	#	more than one PreExt.MS for a given PostExt.MS
	#
	#create two empty matrixs of same dimensions as t for splits and lumps
	#determine if each cell was part of a splitting or lumping event
	#if part of a splitting event, cell gets an 's'
	#if part of a lumping event, cell gets an 'l'
	tsplit = matrix(NA,nrow(t),ncol(t))
	tlump = matrix(NA,nrow(t),ncol(t))
	r = 1; c = 1
	for(r in 1:nrow(t)){
		for(c in 1:ncol(t)){
			if(length(which(t[r,] > 0)) > 1) tsplit[r,][which(t[r,] > 0)] = 's'
			if(length(which(t[,c] > 0)) > 1) tlump[,c][which(t[,c] > 0)] = 'l'
		}
	}
	
	#paste tsplit and tlump matrices together
	#four possible outcomes
	#NA = either no individuals for given combination of PreExt.MS and PostExt.MS or all individuals from PreExt.MS assigned to PostExt.MS
	#'s l' = individuals involved in both a splitting and lumping event (e.g. a shuffle)
	#'s' = individuals involved in only a splitting event (must be paired rowwise with one or more of either 's' or 's l')
	#'l' = individuals involved in a lumping event (must be paired columnwise with one or more of either 'l' or 's l')
	m = matrix(paste(tsplit,tlump),nrow=nrow(t),ncol=ncol(t))
	m = gsub('NA | NA','',m)
	m = gsub('NA',NA,m)
	t
	
	#Extinction
	#	greater than or equal to 5% of individuals in pre-ext morphospecies removed from analysis 
	pre = table(tmp[,'no.extinction'])#number of individuals in each pre-extinction morphospecies pre-extinction
	post = table(tmp[is.na(tmp[,scenarios[i]]),'no.extinction']) #number of individuals removed from each pre-extinction morphospecies
	extinction = rep(0,length(pre))
	extinction[post/pre >= 0.05] = 1 #at least 5% of pre-ext morphospecies individuals extinct
		
	#SBE
	#	for a given PreExt.MS in which extinction occured, one or more instances of splitting (2 or more 's' in row) without involvement in a lumping event
	SBE = apply(m,1,function(x){
		if(all(is.na(x))){	#PreExt.MS with no splitting or lumping
				0
			}else if(all(table(x) >= 2 & names(table(x)) == 's')){ #PreExt.MS with SBE
				1
			}else{	#PreExt.MS with splitting and/or lumping
				0
				}
	})
	
	#remove any SBE for which no extinction occured
	SBE[paste(extinction,SBE) == "0 1"] = 0
	
	sum = cbind(t,ext=extinction,SBE=SBE)	
	
	list[[scenarios[i]]] = sum

}


#Scenarios with SBE 
SBE.names = sort(names(which(unlist(lapply(list,function(x) any(x[,'SBE'] == 1))))))
SBE.names



##################################################
# 2.4) percent of individuals removed by extinction scenario

l = list()

for(i in SBE.names){
	tmp = d[,c('ID','no.extinction',i)]
	pre = table(tmp[,'no.extinction'])#number of individuals in each pre-extinction morphospecies pre-extinction
	post = table(tmp[is.na(tmp[,i]),'no.extinction']) #number of individuals removed from each pre-extinction morphospecies
	l[[i]]$pre = pre
	l[[i]]$post = post
}
x = l[[1]]
sapply(l,function(x) x$post['1']/x$pre['1'])


##################################################
# 2.5) How many times did SBE occur?
N.SBE = length(sort(names(which(unlist(lapply(list,function(x) any(x[,'SBE'] == 1)))))))
N.Total = length(scenarios) - 1 #remove no.extinction
N.SBE/N.Total

##################################################
# 2.6) SBE occured within which pre-extinction morphospecies?

x = lapply(list,function(x){
	names(x[ ,'SBE'])[x[ ,'SBE'] == 1]	
})

x = unlist(x)
tmp = data.frame(scenario=names(x),PreE.MS=x)
rownames(tmp) = NULL
tmp = tmp[order(tmp[,2]),]
t = table(tmp[,2])


##################################################
# 2.7) At what extinction ranges did SBE occur?
#		Figure S3

png(paste0('Fig.S3.png'),width=5,height=5,units='in',res=600,bg='transparent')
#dev.new(width=5,height=5)
par(mar=c(4,4,1,1))

names = names(list)[-1]

min = as.numeric(sapply(strsplit(names,'-'),'[',1))
max = as.numeric(sapply(strsplit(names,'-'),'[',2))

plot(min, max,type="n", xlab="", ylab="",axes=F)

at1 = seq(range(min)[1],range(max)[2],by=1)
at2 = seq(range(min)[1],range(max)[2],by=0.25)
at2 = at2[!at2 %in% at1]
label = rep('',length(at2))

#X axis
axis(1,at=at1,label=at1,cex.axis=1,las=2,tck=-0.025)
axis(1, at2, label=label, tck=-0.01)
mtext('min Bill Depth (mm)',1,line=2.5,cex=1)

#Y axis
axis(2,at=at1,label=at1,cex.axis=1,las=2,tck=-0.025)
axis(2, at2, label=label, tck=-0.01)
mtext('max Bill Depth (mm)',2,line=2.5,cex=1)


clip(par('usr')[1], max(max), par('usr')[3], max(max)) #keep abline within max xlim and ylim
abline(h=at1,lwd=0.25,lty=3,col='black')
abline(v=at1,lwd=0.25,lty=3,col='black')

SBE.names = sort(names(which(unlist(lapply(list,function(x) any(x[,'SBE'] == 1))))))
SBE = rep(F,length(list))
SBE[names %in% SBE.names] = T


points(min[!SBE], max[!SBE],cex=.9, pch=22, col='transparent',bg='grey90',lwd=.5)
points(min[SBE], max[SBE],cex=.9, pch=22, col='transparent',bg='grey25',lwd=.5)


dev.off()











##################################################
# 2.8) Extra code DO NOT RUN


#location of each Geospiza species in PC morphospace
PreE.MS.key = data.frame(PreE.MS=c(1:8),Taxa.new=c('fortis (large-billed) + most conirostris + some magnirostris and propinqua','magnirostris + 1/2 propinqua + some conirostris','propinqua','fortis (small-billed)','septentrionalis + some difficilis','difficilis','fuliginosa + acutirostris','scandens'),location=c('BottomCenter, L of 4','LowerLeft','LeftCenter','BottomCenter, R of 1','btw 6,8','above 7','LowerRight','TopCenter'))

#transfer plots showing SBE from plots producted by 4.Geospiza.plotting.R

file.remove(list.files('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/plots.Bdepth.SML.groundfinches/SBE',full.names=T))
files = list.files('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/plots.Bdepth.SML.groundfinches',full.names=T)
files.to.move = files[grep(paste(SBE.names,collapse='|'),files)]
file.copy(files.to.move,gsub('plots.Bdepth.SML.groundfinches/','plots.Bdepth.SML.groundfinches/SBE/', files.to.move))

file.remove(list.files('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/plots.Bdepth.all/SBE1',full.names=T))
files = list.files('~/Dropbox/SBE/speciation-by-extinction_github/geospiza/plots.Bdepth.all',full.names=T)
files.to.move = files[grep(paste(SBE.names,collapse='|'),files)]
file.copy(files.to.move,gsub('plots.Bdepth.all/','plots.Bdepth.all/SBE1/', files.to.move))


