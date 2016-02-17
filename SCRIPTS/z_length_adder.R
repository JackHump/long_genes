#script that takes whatever is defined as data and adds the zscore and length and other jazz 
data <- read.table(paste("deseq_",sample,"_differential_expression.tab",sep=''),header=T)
data$zscore <- qnorm(p = (data$pvalue)/2)
data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
names(data)[names(data) == "zscore"] <- paste(sample,"_z_score",sep='')
data$length <- data$end_position - data$start_position
