#Code for making fancy graphs with base plotting. I've since moved on to GGplot2
#Example: C9 vs Tau
#colour genes >500kb as red
#colour genes >100kb as orange
#colour SNORD genes in blue
plot(isaacs_c9_z_score,isaacs_tau_z_score, xlim = c(-15,15),ylim = c(-15,15),xlab = "C9 Z-score",ylab = "Tau Z-score",
     main = "Comparing Isaacs", 
     col=ifelse(length.x > 500000, "red", ifelse(length.x < 5000, "blue",ifelse(length.x >250000,"green","black"))))
     col=ifelse(length.x > 500000, "red", ifelse(length.x < 100, "blue",ifelse(length.x >250000,"green","black"))), pch = ifelse(grepl("SNORD",external_gene_id.x),3,1))

     col=ifelse(length.x > 500000, "red", ifelse(length.x > 250000, "orange",ifelse(length.x >100000,"blue","grey"))),

     col = ifelse(length.x > 500000, "red",ifelse(grepl("SNORD",external_gene_id.x),3,"grey"))
abline(h=c(-1.96,1.96),v=c(-1.96,1.96), col = "black",lty=2) #plot lines at +/-1.96

#add text in the right place
text(x=4,y=-4, label = "length > 100kb",col="blue")
text(x=4,y=-5, label = "length > 500kb",col="red")
text(x=4,y=-4.5, label = "length > 250kb",col="orange")


text(x=10,y=-11, label = "SNORD",col="green")

attach(isaacs_merge)
plot(isaacs_c9_z_score,isaacs_tau_z_score, 
     xlim = c(-15,15),ylim = c(-15,15),
     xlab = "C9 Z-score",ylab = "Tau Z-score",
     main = "Comparing Isaacs",
     col = ifelse(length.x > 500000, "red",ifelse(length.x > 250000,"orange",ifelse(length.x > 100000, "blue",ifelse(grepl("SNORD",external_gene_id.x),3,"grey")))))
abline(h=c(-1.96,1.96),v=c(-1.96,1.96), col = "black",lty=2) #plot lines at +/-1.96
text(x=10,y=-9.25, label='length > 100kb',col="blue")
text(x=10,y=-10.5, label='length > 250kb',col="orange")
text(x=10,y=-11.75, label='length > 500kb',col="red")
text(x=10,y=-13, label='SNORD',col="green")
setwd('long_genes/')

attach(cleve_merge)
pdf(cleveland_comparison_length.pdf)
plot(x=cleveland_tdp_z_score,y=cleveland_fus_z_score,
     xlim=c(-15,15),ylim=c(-15,15),
     xlab="TDP KD Z-score",ylab = "FUS KD Z-score", 
     main = "Comparing Clevelands",
     col=ifelse(length.x > 500000, "red", ifelse(length.x > 250000, "orange",ifelse(length.x >100000,"blue","grey"))))
text(x=-3.5,y=-14,label="FUS",cex = 0.6,col="black")
text(x=-15,y=0.8,label="TDP",cex = 0.6,col="black")
abline(h=c(-1.96,1.96),v=c(-1.96,1.96), col = "black",lty=2) #plot lines at +/-1.96
text(x=10,y=-9.25, label='length > 100kb',col="blue")
text(x=10,y=-10.5, label='length > 250kb',col="orange")
text(x=10,y=-11.75, label='length > 500kb',col="red")
dev.off()

attach(fratta_merge)
plot(x=fratta_merge$fratta_embryo_f210i_z_score,y=fratta_merge$fratta_adult_f210i_z_score,
    xlim =c(-15,15),ylim = c(-15,15),
    col=ifelse(length.x <1000 & chromosome_name.y == 7, "red", "grey"))

# Creating Cumulative Frequency plots of all gene lengths AND the up/downregulated fraction
cleveland_tdp_downregulated <- subset(cleveland_tdp, cleveland_tdp_z_score < -1.96)
fratta_adult_upregulated <- subset(fratta, fratta_adult_f210i_z_score > 1.96)
fratta_embryo_upregulated <- subset(fratta_emb, fratta_embryo_f210i_z_score > 1.96)
plot(c(2,7),c(0.0,1.0), type ="n",xlab = "log10(Gene length)",ylab="Cumulative fraction")
legend(2,1,c("All genes","Downregulated in TDP KD","Upregulated in TDP F210I(embryo)"),lty = c(1,1,1,1),lwd = c(2.5,2.5,2.5),col = c("black","blue","red"),cex = 0.8)
legend(2,1,c("All genes","Downregulated in TDP KD","Upregulated in TDP F210I (embryo)"),lty = c(1,1,1,1),lwd = c(2.5,2.5,2.5),col = c("black","blue","red"),cex = 0.8)
lines(ecdf(log10(fratta_embryo_upregulated$length)),lwd=2.5,col="red")
lines(ecdf(log10(cleveland_tdp_downregulated$length)),lwd=2.5,col="blue")
lines(ecdf(log10(cleveland_tdp$length)),lwd=2.5,col="black")
title(main = "Cumulative Distribution of Gene Lengths")
abline(h=0.5,lty=2)

