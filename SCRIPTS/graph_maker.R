## Script to make graphs for all my lovely RNA seq datasets
setwd('project/Humphrey_C9orf72_longGenes/long_genes/DESEQ/')
file_list <- list.files(getwd())

for(i in file_list){
  sample <- gsub('deseq_','',x = i)
  sample <- gsub('_differential_expression.tab','',sample)
  print("processing for:")
  print(sample)
  data <- read.table(i,header=T)
  #create new columns of length and Z-score
  data$zscore <- qnorm(p = 1 - data$pvalue/2)
  data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
  data$length <- data$end_position - data$start_position
  
  
  
  #remove small genes (below 1000bp)
  data <- subset(data, length > 1000)
  
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  average.z.score.by.group <- tapply(data$zscore, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  pdf(paste('../',sample,"/figs/", sample, "_ranked_Zscore_bylength_1kb_plus.pdf",sep = ''))
  plot(x = mean.length.by.group, xaxt = "n", xlab ="", ylab = "Mean gene length by group / base pairs", main = paste(sample, "grouped genes ranked by length (>1kb)", sep = ' '))
  labels <- c("Most downregulated","Not regulated","Most upregulated")
  labelpos <- c(25,100,175) / 192 # so that the labels are constant despite how many groups there are
  mtext(labels, side = 1,line = 0.5,at = (labelpos * length(mean.length.by.group)))
  dev.off()
  
  #run again but set cutoff at 5kb
  data <- subset(data, length > 5000)
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  average.z.score.by.group <- tapply(data$zscore, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  pdf(paste('../',sample,"/figs/",sample,"_ranked_Zscore_bylength_5kb_plus.pdf",sep = ''))
  plot(x = mean.length.by.group, xaxt = "n", xlab ="", ylab = "Mean gene length by group / base pairs", main = paste(sample, "grouped genes ranked by length (>5kb)", sep = ' '))
  labels <- c("Most downregulated","Not regulated","Most upregulated")
  labelpos <- c(25,100,175) / 192 # so that the labels are constant despite how many groups there are
  mtext(labels, side = 1,line = 0.5,at = (labelpos * length(mean.length.by.group)))
  dev.off()
  
}
#compile lists of long genes

for(i in file_list){
  sample <- gsub('deseq_','',x = i)
  sample <- gsub('_differential_expression.tab','',sample)
  print("processing for:")
  print(sample)
  data <- read.table(i,header=T)
  #create new columns of length and Z-score
  data$zscore <- qnorm(p = 1 - data$pvalue/2)
  data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
  data$length <- data$end_position - data$start_position
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  average.z.score.by.group <- tapply(data$zscore, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  long_genes <- data[order(subset(data, ranked.group == 1 | ranked.group == 2)$length,decreasing = T),c(8)]
  long_genes
  write.table(x = long_genes, file = paste('../',"HITS/",sample,"_top_hits.txt",sep = ''), row.names = F, quote = FALSE, sep = '\t',col.names = "")
}

#compare Z-scores


# setwd('project/Humphrey_C9orf72_longGenes/long_genes/')
# file_list <- list.files(getwd())