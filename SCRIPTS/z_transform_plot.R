#
data <- "../../scratch/IoN_RNASeq//Fratta_RNASeq/Fratta_mouse_F210I_adult/output/deseq2/WT_HE/deseq_Fratta_deseq_WTvsHE_differential_expression.tab"
#data <- "../../scratch/IoN_RNASeq//Cleveland/TDP43/processed//deseq2/control_KD//deseq_TDP43_differential_expression.tab"
data <- read.table(data,header=T)
data$Z <- ifelse(data$log2FoldChange < 0, qnorm(data$pvalue/2,mean=0,sd=1),qnorm(1-data$pvalue/2,mean=0,sd=1))
data$Z <- ifelse(data$log2FoldChange < 0, qnorm(data$basic.pval,mean=0,sd=1),qnorm(1-data$basic.pval,mean=0,sd=1))

plot(data$end_position - data$start_position, data$Z,xlab="gene length/bp",ylab="Degree of regulation",main ="Cleveland TDP",type = "p")
data$length <- data$end_position - data$start_position
lmdata <- data[order(data$Z),]
data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
average.z.score.by.group <- tapply(data$Z, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
par(mfrow=c(2,1))
plot(x = average.z.score.by.group, y = mean.length.by.group)
plot(mean.length.by.group)
lm