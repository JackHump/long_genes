#library(ggplot2) 
# 
# getArgs <- function() {
#   myargssvseq.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
#   myargs <- lapply(myargs.list,function(x) x[2] )
#   names(myargs) <- lapply(myargs.list,function(x) x[1])
#   return (myargs)
# # }
# 
# ### Just for debugging 
# support.frame <- "/cluster/project8/vyp/Tabrizi_Huntington_RNASeq/support/htt_combined.tab"
# code <- "htt_combined"
# keep.dups <- TRUE  
# annotation.file <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human/biomart/biomart_annotations_human.tab"
# iFolder <- "/scratch2/vyp-scratch2/Tabrizi_Huntington_RNASeq/processed/Combined/"

#Jack debugging
support.frame <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/design_jack.tab"
deseq.counts <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/counts_Fratta_deseq_WTvsHE.RData"
#deseq.counts <- "scratch/IoN_RNASeq/Fratta_RNASeq/Fratta_mouse_F210I_adult/output/deseq2/rpkm_values.csv"
#annotation.file <- "project/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"

remove_WT3 <- T
sv_num <- 1

# First SV with WT3 removed: 0.3264174  0.5272082 -0.1186758 -0.2566292 -0.6870269  0.2488792 -0.0401729
# First SV with all samples: 0.340689204  0.522768253 -0.353942358  0.007053147 -0.160883346 -0.620246847  0.273229841 -0.008667895
# Two SVs with WT3 removed:

########################## read arguments
# 
# myArgs <- getArgs()
# if ('support.frame' %in% names(myArgs)) support.frame <- myArgs[['support.frame']]
# if ('code' %in% names(myArgs)) code <- myArgs[['code']]
# if ('iFolder' %in% names(myArgs)) iFolder <- myArgs[['iFolder']]
# if ('annotation.file' %in% names(myArgs)) annotation.file <- myArgs[['annotation.file']]
# if ('keep.dups' %in% names(myArgs)) keep.dups <- as.logical(myArgs[['keep.dups']])

###check input files and data frame

support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
load(deseq.counts)
# JACK - to remove WT3 from SVA calculation
if(remove_WT3 == T){support <- support[which(support$sample != "3-SC-F210-WT3"),]}
if(remove_WT3 == T){genes.counts <- genes.counts[,which(colnames(genes.counts) != "3-SC-F210-WT3") ]}

# annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('', 'NA'), quote = "" )
# # trying to agree on the column names
# names(annotation) <- ifelse (names(annotation) == "external_gene_name", "external_gene_id", names(annotation)) 
# 
# ### deseq output folders and files
# deseq2.folder <- paste(iFolder, '/deseq2', sep = '')
# for (folder in c(deseq2.folder)) {
#   if (! file.exists(folder)) dir.create(folder)
# }
# deseq2.folder <- "long_genes/fratta_adult_f210i/SVA/"
# ########## load the count data
# 
# if (keep.dups) deseq.counts <- paste(deseq2.folder, '/deseq_counts_', code, '_keep_dups.RData', sep = '')  
# if (!keep.dups) deseq.counts <- paste(deseq2.folder, '/deseq_counts_', code, '.RData', sep = '')  
#JACK -to remove WT3 from gene counts:


## Check that the support file samples are in the same order as the genes.counts 
print(table(support[,1] == colnames(genes.counts))) 
print(head(support)) 
print(head(genes.counts)) 
support <- support[match(colnames(genes.counts), support[,1]), ]  

print(table(support[,1] == colnames(genes.counts))) 

total.per.sample <- colSums(genes.counts) 

genes.ratios <- 1e6*genes.counts/total.per.sample

mean.ratios <- rowMeans(genes.ratios) 
low.counts <- which(mean.ratios < 0.5) 

mean.ratios <- rowMeans(genes.ratios[-low.counts,]) 
sd.ratios <- apply(genes.ratios[-low.counts,], 1, FUN = sd)

genes.ratios <- genes.ratios[-low.counts,]
# genes.zscore <- genes.ratios 
# for(i in 1:nrow(genes.ratios) ) { 
#    genes.zscore[i,] <- (genes.ratios[i,]-mean.ratios[i])/sd.ratios[i]
# } 
# 
# fig.folder <- paste0(iFolder, "/rna-seqc/figs/") 
# if (! file.exists(fig.folder)) dir.create(fig.folder)
# 
# plot.fname <- paste0(fig.folder, "zscores_distribution.pdf") 
# 
# pdf(plot.fname) 
# breaks <- seq(-5, 15, 0.2)
# 
# h <- hist(genes.zscore[,1], breaks = breaks, plot=F) 
# cols <- rainbow(ncol(genes.zscore))
# plot(x=h$mids, y = h$density, type = "l", xaxt="n", ylim = c(0,1), xlim = c(-4,6), col = cols[1], 
#      ylab = "density", xlab = "zscore of gene counts")
# axis(1,at=h$breaks,labels=h$breaks)
# text(x = 4.4, y = 1, col = cols[1], labels = colnames(genes.zscore)[1])
# 
# for(i in 2:ncol(genes.zscore)) { 
#   h <- hist(genes.zscore[,i], breaks = breaks, plot=F) 
#   lines(x=h$mids, y = h$density, type = "l", xaxt="n", col = cols[i])
#   text(x = 4.4, y = 1 - (i -1)*0.1, col = cols[i], labels = colnames(genes.zscore)[i])
# } 
# title(main="Fratta Adult F210I")
# 
# dev.off() 

######################################################
## using sva to find latent variables 
library(sva) 

cleanY = function(y, mod, svs) {
    X = cbind(mod, svs)
    Hat = solve(t(X) %*% X) %*% t(X)
    beta = (Hat %*% t(y))
    rm(Hat)
    gc()
    P = ncol(mod)
    return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)

# Assume that we only have one condition here
if (length(list.conditions) != 1) { 
   stop("Can only handle one condition") 
} 

groups = as.factor(support[,list.conditions[1]]) 

mod1 = model.matrix(~groups) 
mod0 = cbind(mod1[,1]) 

#n.sv specifies number of surrogate variables to be removed
svseq = svaseq(genes.counts[-low.counts,], mod1, mod0, n.sv = sv_num)$sv 

# log.gene.counts <- log(genes.counts + 1) 
# 
# log.clean.gene.counts <- cleanY(log.gene.counts, mod1, svseq) 
# clean.genes.counts <- exp(log.clean.gene.counts) - 1
# 
# 
# #plot.fname <- paste0(fig.folder, "zscores_distribution_after_sva.pdf")
# 
# #pdf(plot.fname)
# 
# total.per.sample <- colSums(clean.genes.counts) 
# genes.ratios <- 1e6*clean.genes.counts/total.per.sample
# 
# mean.ratios <- rowMeans(genes.ratios) 
# low.counts <- which(mean.ratios < 0.5) 
# 
# mean.ratios <- rowMeans(genes.ratios[-low.counts,]) 
# sd.ratios <- apply(genes.ratios[-low.counts,], 1, FUN = sd)
# 
# genes.ratios <- genes.ratios[-low.counts,]
# genes.zscore <- genes.ratios 
# for(i in 1:nrow(genes.ratios) ) { 
#    genes.zscore[i,] <- (genes.ratios[i,]-mean.ratios[i])/sd.ratios[i]
# }
#  
# breaks <- seq(min(genes.zscore)-0.2, max(genes.zscore)+0.5, 0.2)
# h <- hist(genes.zscore[,1], breaks = breaks, plot=F)
# cols <- rainbow(ncol(genes.zscore))
# plot(x=h$mids, y = h$density, type = "l", xaxt="n", ylim = c(0,1), xlim = c(-4,6), col = cols[1],
#      ylab = "density", xlab = "zscore of gene counts")
# text(x = 4.4, y = 1, col = cols[1], labels = colnames(genes.zscore)[1])
# axis(1,at=h$breaks,labels=h$breaks)
# 
# for(i in 2:ncol(genes.zscore)) {
#   h <- hist(genes.zscore[,i], breaks = breaks, plot=F)
#   lines(x=h$mids, y = h$density, type = "l", xaxt="n", col = cols[i])
#   text(x = 4.4, y = 1 - (i -1)*0.1, col = cols[i], labels = colnames(genes.zscore)[i])
# }
# 
# title(main= "Fratta F210I Adult - corrected for 1 surrogate variable")
# dev.off()
if(!is.null(dim(svseq)) ){
  for(i in 1:sv_num){
    support <- cbind(support,svseq[,i])
    col.name <- paste0("covariate.",i)
    names(support)[i + 2] <- col.name
    }}else{support$covariate <- svseq}
#write out support table

outfile <- "project/Humphrey_C9orf72_longGenes/long_genes/fratta_adult_f210i/wt3_removed/final_go/design"
if(remove_WT3 == T){outfile <- paste0(outfile, "_wt3_removed") 
  }else{outfile <- paste0(outfile, "_all_samples")
  }
if(sv_num > 0){outfile <- paste0(outfile,"_",sv_num,"sva")}
outfile <- paste0(outfile,".tab")
write.table(support,outfile)

