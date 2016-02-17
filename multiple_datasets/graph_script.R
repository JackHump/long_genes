# script for creating multipanel plot of all datasets
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
# to convert list in the terminal to an R list: sed 's/^/"/g' | sed 's/$/",/g' | tr "\n" " "
datasets <-c("project/Humphrey_C9orf72_longGenes/long_genes/cleveland_fus/deseq_cleveland_fus_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/cleveland_tdp/deseq_cleveland_tdp_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/deseq_fratta_embryo_f210i_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/deseq_fratta_adult_f210i_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/isaacs_c9/deseq_isaacs_c9_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/isaacs_tau/deseq_isaacs_tau_differential_expression.tab")


codes <- c("cleveland_fus","cleveland_tdp","fratta_f210i_1_embryo","fratta_f210i_2_adult","isaacs_c9","isaacs_tau")
output.pdf <- "project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/rank_all_datasets.pdf"
directory <- "project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/ranks/"
title_list <- c("Mouse FUS knockdown", "Mouse TDP43 knockdown", "Mouse Embryo TDP43 F210I", "Mouse Adult TDP43 F210I","Human C9orf72 FTD", "Human Tau FTD")


##Looking at removing each sample sequentially
title_list <- c("HET1_removed", "HET2_removed", "HET3_removed", "HET4_removed", "All samples, 1 Surrogate Variable", "All samples, 2 Surrogate Variables", "All samples", "WT1_removed", "WT2_removed", "WT3_removed, 1 Surrogate Variable", "WT3_removed, 2 Surrogate Variables", "WT3_removed", "WT4_removed")
codes <- c("HET1_removed", "HET2_removed", "HET3_removed", "HET4_removed", "1_sva", "2_sva", "untouched", "WT1_removed", "WT2_removed", "WT3_removed_1_sva", "WT3_removed_2_sva", "WT3_removed", "WT4_removed")
datasets <- c("fratta_f210i_adult_HET1_removed_results.tab", "fratta_f210i_adult_HET2_removed_results.tab", "fratta_f210i_adult_HET3_removed_results.tab", "fratta_f210i_adult_HET4_removed_results.tab", "fratta_f210i_adult1_sva_results.tab", "fratta_f210i_adult2_sva_results.tab", "fratta_f210i_adult_results.tab", "fratta_f210i_adult_WT1_removed_results.tab", "fratta_f210i_adult_WT2_removed_results.tab", "fratta_f210i_adult_WT3_removed1_sva_results.tab", "fratta_f210i_adult_WT3_removed2_sva_results.tab", "fratta_f210i_adult_WT3_removed_results.tab", "fratta_f210i_adult_WT4_removed_results.tab")
rankings <-c("1_sva_mean_length_by_group.tab", "2_sva_mean_length_by_group.tab", "HET1_removed_mean_length_by_group.tab", "HET2_removed_mean_length_by_group.tab", "HET3_removed_mean_length_by_group.tab", "HET4_removed_mean_length_by_group.tab", "untouched_mean_length_by_group.tab", "WT1_removed_mean_length_by_group.tab", "WT2_removed_mean_length_by_group.tab", "WT3_removed_1_sva_mean_length_by_group.tab", "WT3_removed_2_sva_mean_length_by_group.tab", "WT3_removed_mean_length_by_group.tab", "WT4_removed_mean_length_by_group.tab")
setwd("project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/wt3_removed/")
output.pdf <- "length_ranking_removal.pdf"
directory <- getwd()


for(i in 1:length(codes)){
  data <- read.table(datasets[i],header=T)
  data$zscore <- qnorm(p = 1 - data$pvalue/2)
  data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
  data$length <- data$end_position - data$start_position
  #remove genes shorter than 1kb
  data <- subset(data,length > 1000)
  data <- subset(data, baseMean > 0.1)
  data <- data[ order(data$zscore), ]
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  #for writing out the mean lengths per group
  output.rank <- paste0(directory, "/", codes[i], "_mean_length_by_group.tab")
  write.table(mean.length.by.group, output.rank)
}
##START HERE
output.pdf <- "length_ranking_WT3_removal_sva.pdf"
pdf(output.pdf, onefile = T)

par(mfrow = c(3,2))
#need to readjust margins for this setup

#files_for_graphing <- list.files(directory)[grep("mean_length", list.files(directory))]
for(i in c(12,10,11)){
  data <- read.table(paste0(directory, "/", rankings[i]))
  y_limit = 250
  #if (i ==2){y_limit = 500}
  plot(x = data$x/1000, xaxt = "n", xlab ="", ylab = "Mean gene length by group / kb", cex.lab = 0.8, ylim = c(0,y_limit), pch = 20, col = ggplotColours(n=length(codes))[i], )
  labels <- c("Most downregulated","Not regulated","Most upregulated")
  position <- length(mean.length.by.group)
  mtext(labels, side = 1,line = 0.1,at = c(0.125 * position,position/2,0.875 * position),cex = 0.5)
  title(main = title_list[i],cex.main = 0.95)  
}
dev.off()

#datasets 2 and 3
#load each dataset, tack on Z-scores, exclude any genes with z-scores inside -1.96 - +1.96
#create table of genes that significantly up in Fratta embryo and down in Cleveland TDP
codes <- c("blank","tdp_kd","tdp_")

for(i in 2:3){
  data <- read.table(datasets[i],header=T)
  data$zscore <- qnorm(p = 1 - data$pvalue/2)
  data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
  data$length <- data$end_position - data$start_position
  #remove genes shorter than 1kb
  #data <- subset(data,length > 1000)
  data <- subset(data, zscore > 1.96 | zscore < -1.96)
  data <- subset(data, baseMean > 0.1)
  data <- data[ order(data$zscore), ]
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  #for writing out the mean lengths per group
  output.rank <- paste0(directory, "/", codes[i], "_mean_length_by_group.tab")
  write.table(mean.length.by.group, output.rank)
}

#Comparing TDP knockdown and TDP F210I datasets
#comparison table made by loading both datasets, calculating zscores and then merging the relevent columns together
data <- read.table(datasets[i],header=T)
data$zscore <- qnorm(p = 1 - data$pvalue/2)
data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
data$length <- data$end_position - data$start_position
data <- subset(data, zscore > 1.96 OR zscore < -1.96)
data <- subset(data, zscore > 1.96 | zscore < -1.96)
comp_table <- data[,c(1,8,14,6,7,13)]
names(comp_table)[4:6] <- c("TDP.kd.pvalue","TDP.kd.padj","TDP.kd.zscore")
comp_table$TDP.f210i.pvalue <- data$pvalue[match(x = comp_table$EnsemblID,table = data$EnsemblID)]
comp_table$TDP.f210i.paj <- data$padj[match(x = comp_table$EnsemblID,table = data$EnsemblID)]
comp_table$TDP.f210i.zscore <- data$zscore[match(x = comp_table$EnsemblID,table = data$EnsemblID)]
comp_clean <- comp_table[which(!is.na(comp_table$TDP.f210i.zscore)),]
comp_clean_anno <- comp_clean
comp_clean_anno$Gene.Type <- genes$Gene.type[match(x = comp_clean$EnsemblID, table = genes$Ensembl.Gene.ID)]

#comparison plots

ggplot(data = comp_clean_anno, aes(y=log10(length))) + 
  geom_point(aes(x=TDP.kd.zscore),colour = "red") + 
  geom_point(aes(x=TDP.f210i.zscore),colour="blue") + 
  ggtitle(label = "TDP43 KD in adult mouse brain(red) vs. TDP43 F210I embryonic mouse brain (blue)") + 
  xlab(label = "Z score") + 
  ylab("log10(Gene length/bp)")

ggplot(comp_clean_anno, aes(x = Direction, y = log10(length)),colour="blue") + 
  geom_jitter(na.rm = T,aes(colour=Direction)) +
  scale_colour_manual(values=c("blue","darkgreen"))+
  #scale_fill_brewer(palette="Spectral") +
  xlab("Direction of regulation in both datasets") + 
  ylab("log10(Gene length/bp)") + 
  ggtitle("Comparing TDP43 knockdown and TDP43 F210I datasets")
