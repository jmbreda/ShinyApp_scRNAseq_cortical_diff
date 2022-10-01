library(tidyverse)
Datasets=c('ES','iPS','Nas2')

for (d in Datasets){
	print(d)
	
	mydata <- read_delim(paste0('~/SingleCell_Analysis/TrueVariance/output/Gene/Zahra/',d,'/log_transcription_quotients.txt'),'\t',col_types=cols("GeneID"=col_character(),.default=col_double()))
	#mydata = read.table(paste0('~/SingleCell_Analysis/TrueVariance/output/Gene/Zahra/',d,'/log_transcription_quotients.txt'),header = TRUE, row.names=1, sep = "\t")	
	save(mydata, file = paste0('ltq/',d,'/log_transcription_quotients.RData'))
}
