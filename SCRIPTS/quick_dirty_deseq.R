library(DESeq2)
#Quick and Dirty DESEQ 
#a shorter deseq script adapted from Kitty's script
#requires a counts file and a support file, annotation is optional (if you want gene names and lengths etc)
#counts must be raw numbers. Colnames = sample IDs, rownames = ensemblIDS
#support: 1 column is sample IDs (must match counts), 2nd column is "condition" - put your control/wildtype first alphabetically as Deseq turns the values into factors with a hierarchy based on alphabetical order
#support must contain the same samples as in counts. If you don't want to test a particular sample then make its condition NA.
#if you have covariates then put "covariates" column, as many as you like.

# this hasn't been cleaned up.

#rm(list=ls())

getArgs <- function() {
    myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
    myargs <- lapply(myargs.list,function(x) x[2] )
    names(myargs) <- lapply(myargs.list,function(x) x[1])
    return (myargs)
  }

myArgs <- getArgs()
  if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
  if ('wt3' %in% names(myArgs)) wt3 <- myArgs[['wt3']]
print(wt3)

# covariates = 2
# WT_removed = 0
# HET_removed = 0
# WT3_removed = T

#support.frame <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/design_wt3_removed_2sva.tab"
# Fratta Embryonic Homozygous F210I
support.frame <- 
deseq.counts <- "counts_embryo_f210i.RData"
load(deseq.counts)

##Running DESeq2 on fratta F210I and creating various graphs
library(DESeq2)
#library(qqman)
#library(ggplot2)
#library(gridExtra)

# Pietro Adult F210I
#support.frame <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/design_jack.tab"
annotation.file <- "project/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"
#loc.code <- "fratta_f210i_adult"
loc.code <- "fratta_f210i_embryo_hom"
#deseq2.figs <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/wt3_removed/"
plots <- F
extra <- F
# if(WT_removed > 0){loc.code <- paste0(loc.code,"_WT",WT_removed,"_removed")}
# if(HET_removed > 0){loc.code <- paste0(loc.code,"_HET",HET_removed,"_removed")}
# if(WT3_removed == T){loc.code <-paste0(loc.code,"wt3_removed")}
# if(covariates > 0){loc.code <- paste0(loc.code,covariates,"_sva")}
#regular counts
deseq.counts <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/counts_Fratta_deseq_WTvsHE.RData"
deseq.counts <- "../../counts_Fratta_deseq_WTvsHE.RData"
if(wt3 == "remove"){deseq.counts <- "../../genes_counts_wt3_removed.RData"}
##support containing 3 SVAs as covariates
#support.frame <- "long_genes/fratta_adult_f210i/SVA/fratta_f210i_adult_support_3_sva.tab"
#clean counts from SVA
#counts <- "/Users/Jack/long_genes/fratta_adult_f210i/SVA/counts_fratta_f210i_adult_2_sva_clean.RData"
#counts <- "long_genes/fratta_adult_f210i/SVA/counts_fratta_f210i_adult_1_sva_clean.RData"

#set all adjusted count values as integers
#load(counts)
#genes.counts <- apply(X=clean.genes.counts, MARGIN = c(1,2),FUN = function(x) as.integer(x))

### HOMOZYGOUS F210I vs. CONTROL

support.frame <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/ctl_hom/design.tab"
deseq.counts <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/counts_embryo_f210i.RData"
annotation.file <- "project/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"
loc.code <- "fratta_f210i_embryo_hom"
deseq2.figs <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/ctl_hom/"

## read in support table, create lists of conditions and covariates
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
# if(WT_removed > 0){support[WT_removed,2] <- NA}
# if(HET_removed > 0){support[(HET_removed+4),2] <- NA}
# if(covariates == 0){support <- support[,1:2]
# }else {support <- support[,seq(1:(covariates+2))]}
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE) #looks for what conditions are referred to, either "condition" or "conditions"
list.covars <- grep(names(support), pattern = '^covar.*', value  = TRUE)

## read in annotation file

annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('', 'NA'), quote = "" )
names(annotation) <- ifelse (names(annotation) == "external_gene_name", "external_gene_id", names(annotation)) # trying to agree on the column names
# 
print(head(annotation))

#deseq.counts <- counts
load(deseq.counts)

### Loop over all proposed conditions
for (condition in list.conditions) {
  num.cond <- FALSE 
  message("Processing for ", condition) 
  samples.to.use   <- !is.na(support[,condition]) #creates a vector of TRUE/FALSE values to use to filter out samples you're not interested in comparing
  
  genes.counts.loc <- genes.counts[, samples.to.use ] #from gene counts RData, take only the columns of the samples you're interested in
  
  support.loc <-  support[  samples.to.use, ] #extracts only the relevent samples from the support file
  
  ## determines if the condition be treated as a factor or a numeric variable 
  if (substr(condition, 10, 12) == "Num") { 
    num.cond <- TRUE 
    support.loc$condition <- as.numeric(support.loc[, condition])
    loc.code <-  condition
  } else { 
    support.loc$condition <- factor(support.loc[, condition]) #ensures that each condition is treated as a factor
    #loc.code <-  paste(unique(support.loc$condition), collapse = '_') #creates a code "condition1_condition2". This will be in the alphabetical order of the conditions
  } 

  use.covars <- FALSE
  if (length(list.covars) > 0) {
    use.covars <- TRUE
  }
  
  if (use.covars) {
    message("Using ", length(list.covars), " covariates") 
    formula1 <- paste(" ~ ", paste(list.covars, collapse = "+"), " + ", condition, sep = "") 
    formula1 <- as.formula(formula1) 
    formula0 <- paste(" ~ ", paste(list.covars, collapse = "+"), sep = "") 
    formula0 <- as.formula(formula0) 
    
    design.deseq <- support.loc[, which(names(support.loc) %in% c(condition, list.covars))]
  } else {
    formula1 <- as.formula(paste0("~ ", condition))
    formula0 <- as.formula("~ 1")
    design.deseq <- support.loc[, c(condition), drop = FALSE] #take just the condition column of the relevent samples from the support file
  }
  head(genes.counts.loc)
  print(design.deseq)
  
    CDS <- DESeqDataSetFromMatrix(countData = genes.counts.loc, colData = design.deseq, design = formula1) #create the DESeqDataSet from the counts tables according to the conditions in the design table and according to the formula.
  
  #################### Do the actual model fitting 
  CDS <- DESeq(CDS, test = "LRT", reduced = formula0, 
               minReplicatesForReplace = 5 ) 
  deseq.res <- results(CDS)  ## extract the results from the DESEq table
  
  ############# Make the results table into a sensible format 
  deseq.res.df <- data.frame(deseq.res) ## convert the results table into a data frame 
  print(head(deseq.res.df)) 
  deseq.res.df$EnsemblID <- row.names( deseq.res.df) #add EnsemblID column by transferring across the row names
  deseq.res.df <- merge(deseq.res.df, annotation, by = 'EnsemblID', all.x = TRUE) #merge with annotation table to give gene name and length info
  deseq.res.df <- deseq.res.df[order(deseq.res.df$pvalue),] #sort by p value
  deseq.res.df$length <- deseq.res.df$end_position - deseq.res.df$start_position
  print(head(deseq.res.df)) 
}

output.table <- paste0(deseq2.figs, "/", loc.code, "_results.tab")
write.table(deseq.res.df,file = output.table)
sig <- table(deseq.res.df$padj < 0.1)
print(sig)
outfile <- paste0("sig_genes/", support.frame,".sig_genes.tab")
write.table(sig,outfile)
#output.sig <- paste0(deseq2.figs, "/", loc.code, "_num_sig.tab")
#write.table(sig, output.sig)

#output.qq <- paste(deseq2.figs,"/",loc.code, "_qq.pdf",sep= '')

# for creating qq plots

# pdf(output.qq)
qq(deseq.res.df$pvalue, main = loc.code)
# dev.off()

# for creating a grand plot
# load results as 'data' and write name of dataset as `loc.code`
# if (extra = T){
# 
datasets <-c("project/Humphrey_C9orf72_longGenes/long_genes/cleveland_fus/deseq_cleveland_fus_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/cleveland_tdp/deseq_cleveland_tdp_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/deseq_fratta_embryo_f210i_differential_expression.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/fratta_embryo_f210i/ctl_hom/fratta_f210i_embryo_hom_results.tab",
             "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/deseq_fratta_adult_f210i_differential_expression.tab",
            "project/Humphrey_C9orf72_longGenes/long_genes/isaacs_c9/deseq_isaacs_c9_differential_expression.tab",
            "project/Humphrey_C9orf72_longGenes/long_genes/isaacs_tau/deseq_isaacs_tau_differential_expression.tab")
codes <- c("cleveland_fus","cleveland_tdp","fratta_f210i_1_embryo_het","fratta_f210i_embryo_hom","fratta_f210i_3_adult_het","isaacs_c9","isaacs_tau")
# 
 for(i in 1:7){
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
  output.rank <- paste0("project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/f210i_hom_included/", codes[i], "_mean_length_by_group.tab")
  write.table(mean.length.by.group, output.rank)
}
#Just the TDP datasets 
dataset_name <- "TDP_datasets"
dataset_nums <- c(3,4,5,2)
# All datasets
dataset_name <- "All_datasets"
dataset_nums <-c(3,4,5,2,1,6,7)

pdf(paste0("project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/f210i_hom_included/",dataset_name,"_long_genes.pdf"))
title_list <- c("Mouse FUS knockdown", "Mouse TDP43 knockdown","Mouse Embryo TDP43 F210I Heterozygous","Mouse Embryo TDP43 F210I Homozygous", "Mouse Adult TDP43 F210I Heterozygous","Human C9orf72 FTD", "Human Tau FTD")
par(mfrow = c(3,2))
for(i in dataset_nums){
  data <- read.table(paste0("project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/f210i_hom_included/", list.files("project/Humphrey_C9orf72_longGenes/long_genes/rank_graph/f210i_hom_included/")[i]))
  y_limit = 250
  #if (i ==2){y_limit = 500}
  plot(x = data$x/1000, xaxt = "n", xlab ="", ylab = "Mean gene length by group / kb", cex.lab = 0.8, ylim = c(0,y_limit), pch = 20, col = ggplotColours(n=6)[i], )
  labels <- c("Most downregulated","Not regulated","Most upregulated")
  position <- length(mean.length.by.group)
  mtext(labels, side = 1,line = 0.1,at = c(0.125 * position,position/2,0.875 * position),cex = 0.5)
  title(main = title_list[i],cex.main = 0.95)  
}
dev.off()


# 
# 
# 
# 
# }

# if (plots == T){
#   
#   output.pca <- paste(deseq2.figs, '/', loc.code, '_pca.pdf', sep = '')
#   
#   vsd <- varianceStabilizingTransformation(CDS)
#   #source('/cluster/project8/vyp/Humphrey_C9orf72_longGenes/c9_deseq//tools.R')
#   source('project/Humphrey_C9orf72_longGenes/c9_deseq/tools.R')
#   myPCA(vsd,output.pca)
#   disp.plot <- paste0(deseq2.figs,"/",loc.code,"_dispersion.pdf")  
#   pdf(disp.plot) 
#   plotDispEsts(CDS)  
#   dev.off() 
#   
#   ###GGplot for labelled counts
#   
#   library(gridExtra)
#   library(ggplot2)
#   #based on this code
#   ggplot(data, aes(x=condition, y=count, color=sample)) +
#     #scale_y_log10() + 
#     geom_point(position=position_jitter(width=.1,height=0.5),size =5) +
#     theme(title = deseq.res.df$external_gene_id[i])
#   
 # sig_pdf <- paste0(deseq2.figs,"/",loc.code,"_sig_genes.pdf")
  
#Significant Genes charts
# 
# sig.genes <- head(which(deseq.res.df$padj < 0.1),5)
#   #pdf(sig_pdf,onefile = T)
#   for (i in sig.genes) {
#     data <- plotCounts(CDS, gene= deseq.res.df$Ensembl[i], intgroup=condition, returnData=TRUE)
#     #data$sample <- c("WT1","WT2","WT3","WT4","HET1","HET2","HET3","HET4")
#     title = paste0(toupper(deseq.res.df$external_gene_id[i]),"\n P: ", round(deseq.res.df$padj[i],10),"\n Length: ", round(deseq.res.df$length[i] / 1e3), "kb")
#     plot <- ggplot(data, aes(x=condition, y=count, color = condition)) + geom_point(position=position_jitter(width=.1,height=0.5),size=5) + ggtitle(title) + theme(plot.title = element_text(lineheight=.8, face="bold"))
#     grid.arrange(plot)
#   } 
#   #dev.off()
#   
  ## create Z score vs length plots
  data <- deseq.res.df
  data$zscore <- qnorm(p = 1 - data$pvalue/2)
  data$zscore<- ifelse ( data$log2FoldChange > 0, abs(data$zscore), - abs(data$zscore))
  data$length <- data$end_position - data$start_position
  #remove genes shorter than 1kb
  data <- subset(data,length > 1000)
  
  data <- data[ order(data$zscore), ]
  data$ranked.group <- rep(1:1000, each = 100)[ 1:nrow(data) ]
  mean.length.by.group <- tapply(data$length, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  average.z.score.by.group <- tapply(data$zscore, FUN = mean, INDEX = data$ranked.group, na.rm = TRUE)
  
  z.plot <- paste(deseq2.figs, '/',loc.code, "_average_z_score_v_length.pdf", sep ='')
  z.ranked.plot <- paste(deseq2.figs, '/', loc.code, "_ranked_z_score_v_length.pdf", sep ='')
  z.length.plot <- paste(deseq2.figs, '/', loc.code, "_z_score_v_length.pdf", sep = '')
  
  pdf(z.plot)
  plot(x = average.z.score.by.group, y = mean.length.by.group, ylab = "mean length by group/bp", main = paste(loc.code, "Z score against length, grouped genes", sep = ' '))
  dev.off()
  
  #this is the one Vincent wants!
  pdf(z.ranked.plot)
  plot(x = mean.length.by.group/1000, xaxt = "n", xlab ="", ylab = "Mean gene length by group / kb", pch = 20, col = ggplotColours(n=3)[1], main = paste(loc.code, "grouped genes ranked by length", sep = ' '))
  labels <- c("Most downregulated","Not regulated","Most upregulated")
  position <- length(mean.length.by.group)
  mtext(labels, side = 1,line = 0.5,at = c(0.125 * position,position/2,0.875 * position))
  dev.off()
  
  df <- data.frame(x = seq(1,length(mean.length.by.group)),length = mean.length.by.group)
  
  plot <- ggplot(data = mean.length.by.group, aes(y = mean.length.by.group) + geom_point())
  
  pdf(z.length.plot)
  plot(y = data$zscore, x = data$length,xlab = "gene length (bp)", ylab="Signed Z score", main = paste(loc.code, "Z score against length, all genes", sep = ' '))
  dev.off()
  
  #function to get those nice ggplot2 colours
ggplotColours <- function(n=6, h=c(0, 360) +15){
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)}
  

# End of extra plots  
#else {stop}




# pdf("pairwiseMAs.pdf")
# MA.idx = t(combn(1:8, 2))
# for( i in  seq_along( MA.idx[,1])){ 
#   MDPlot(counts(CDS, normalized = T)[idx.nz ,], 
#          c(MA.idx[i,1],MA.idx[i,2]), 
#          main = paste( colnames(CDS)[MA.idx[i,1]], " vs ",
#                        colnames(CDS)[MA.idx[i,2]] ), ylim = c(-3,3))
# }
# dev.off()
# 
# dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
#                     PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
#                     sampleNO = colData(rld)$sampleNO,
#                     condition = colData(rld)$condition)
# 
# (qplot(PC1, PC2, data = dataGG, color =  condition, 
#        main = "PC1 vs PC2, top variable genes", size = I(6))
#  + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
#         y = paste0("PC2, VarExp:", round(percentVar[2],4)))
#  + scale_colour_brewer(type="qual", palette=2)
# )
# 
# Pvars <- rowVars(assay(rld))
# select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, 
#                                                       length(Pvars)))]
# 
# PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
# percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
# 
# dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], 
#                     PC3 = PCA$x[,3], PC4 = PCA$x[,4], 
#                     sampleNO = colData(rld)$sampleNO,
#                     condition = colData(rld)$condition)
# 
# (qplot(PC1, PC2, data = dataGG, color =  condition, 
#        main = "Fratta Adult F210I PC1 vs PC2, top variable genes, no outliers", size = I(6))
#  + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)),
#         y = paste0("PC2, VarExp:", round(percentVar[2],4)))
#  + scale_colour_brewer(type="qual", palette=2)
# )