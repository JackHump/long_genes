##Script that takes a biomart list of transcripts, calculates the intron length and then extracts the largest transcript size for each gene.

data <- read.table('Downloads/mart_export.txt',header=T,sep = ",")
#order data by gene (ascending) and by transcript length (descending)
data <- data[with(data, order(Ensembl.Gene.ID,-Transcript.length)),]
#get logical vector where first occurence of gene is TRUE and those following are FALSE:
duplicates <- !(duplicated(data$Ensembl.Gene.ID))
data <- data[duplicates,]
#calculate intron length for each gene
data$Intron.length <- (data$Gene.End..bp.-data$Gene.Start..bp.) - data$Transcript.length
