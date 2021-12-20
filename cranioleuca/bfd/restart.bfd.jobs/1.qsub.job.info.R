#	use  
#	qstat -f | cat > ~/qstat.info
#	to get full details on all qstat jobs. some munging and you can info on what jobs are running and make deletion commands

setwd('~/Dropbox/SBE/speciation-by-extinction_github/cranioleuca/bfd')
lines = readLines('restart.bfd.jobs/qstat.info')

start = grep('Job Id:',lines)
end = c(start[-1] - 1,length(lines))

coords = cbind(start,end)

l = list()
i = 1
for(i in 1:nrow(coords)){
	l[[i]] = lines[coords[i,1]:coords[i,2]]
}

gseeholzer = unlist(lapply(l,function(x) length(grep('gseeholzer',x)) > 1 ))


tmp = l[gseeholzer]

ID = unlist(lapply(tmp,function(x) x[grep('Job Id:',x)] ))
name = unlist(lapply(tmp,function(x) x[grep('Job_Name',x)] ))
stat = unlist(lapply(tmp,function(x) x[grep('job_state',x)] ))

tmp = data.frame(cbind(ID=ID,name=name,stat=stat),stringsAsFactors=F)

tmp$ID   = trimws(gsub('Job Id: ','',tmp$ID))
tmp$stat = trimws(gsub('job_state = ','',tmp$stat))
tmp$name = trimws(gsub('Job_Name = ','',tmp$name))

head(tmp)
tmp$name

running = tmp[tmp$stat == 'R','name']
inqueue = tmp[tmp$stat == 'Q','name']


#create list of delete commands for jobs in queue
#delete = tmp[tmp$stat == 'Q','ID']
delete = tmp[tmp$stat == 'R','ID']; for(i in delete) cat('qdel ',i,'\n',sep='')
delete = c(tmp[tmp$stat == 'Q','ID'],tmp[grep('F5.EXT_v1_1sp.r1',name),'ID'])
writeLines(paste0('qdel ',delete),'restart.bfd.jobs/delete.jobs.sh')
#	chmod u+x delete.jobs.sh
#	./delete.jobs.sh


#resubmit deleted jobs
foo = readLines('runs/F5/qsub.commands.sh') #this is the original full qsub command submission
foo = foo[grep('r1',foo)]

#restart all the ones in queue plus F5.EXT_v1_1sp.r1 to free up a node for Lucas' run
restart = c(foo[grep(paste(inqueue,collapse="|"),foo)],foo[grep('F5.EXT_v1_1sp.r1',foo)])

writeLines(restart,'restart.bfd.jobs/resubmit.qsub.commands.sh')
#	chmod u+x resubmit.qsub.commands.sh
#	./resubmit.qsub.commands.sh



