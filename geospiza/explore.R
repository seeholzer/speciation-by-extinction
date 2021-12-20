extinct = tmp[!is.na(tmp[,2]) & is.na(tmp[,3]),'ID']


data = read.table("Geospiza.data.csv", header=T, sep=",")
data$ID = paste0(data$Institution,data$Museum.Number)
colnames(data)[1] = 'Taxa.Lack'
colnames(data)[6] = 'Taxa.New'
colnames.to.keep = c('ID','Taxa.Lack','Taxa.New','Wing','Tail','Blength','Bdepth','Bwidth','Tarsus')
data = data[,colnames.to.keep]


x = data[grep('magni|fortis|fulig',data$Taxa.Lack), ]
table(x$Taxa.Lack)

hist(data[data$ID %in% extinct,'Bdepth'])
data[data$ID %in% extinct, ]


x = data[grep('fuliginosa',data$Taxa.Lack), ]
x = x[order(x$Bdepth), ]
boxplot(Bdepth ~ Island, data=x)
x[x$Bdepth > 8.0, ]

x = data[grepl('fuliginosa',data$Taxa.Lack) & !data$ID %in% extinct, ]
boxplot(Bdepth ~ Island, data=x)
head(x)
hist(x$Bdepth)