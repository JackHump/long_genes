##Code for making fancy graphs
#Example: C9 vs Tau
#colour genes >500kb as red
#colour genes >100kb as orange
#colour SNORD genes in blue
plot(isaacs_c9_z_score,isaacs_tau_z_score, xlim = c(-15,15),ylim = c(-15,15),xlab = "C9 Z-score",ylab = "Tau Z-score",
     main = "Comparing Isaacs", 
     col=ifelse(length.x > 500000, "red", ifelse(length.x < 5000, "blue",ifelse(length.x >250000,"green","black"))))
     col=ifelse(length.x > 500000, "red", ifelse(length.x < 100, "blue",ifelse(length.x >250000,"green","black"))), pch = ifelse(grepl("SNORD",external_gene_id.x),3,1))

     col=ifelse(length.x > 500000, "red", ifelse(length.x > 250000, "orange",ifelse(length.x >100000,"blue","grey"))),

     col = ifelse(length.x > 500000, "red",ifelse(grepl("SNORD",external_gene_id.x),3,1))
abline(h=c(-1.96,1.96),v=c(-1.96,1.96), col = "black",lty=2) #plot lines at +/-1.96

#add text in the right place
text(x=4,y=-4, label = "length > 100kb",col="blue")
text(x=4,y=-5, label = "length > 500kb",col="red")
text(x=4,y=-4.5, label = "length > 250kb",col="orange")


text(x=10,y=-11, label = "SNORD",col="green")
