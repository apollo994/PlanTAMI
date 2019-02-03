#scipt che prende come input il una lista di geneID di ath, sly, hvu e osa e partendo dai conti normalizzati 
#tra 0 e 1 ne genera un grafico
#rispetto alla versione 1 in questa versione 2, gli input sono 2 liste di geni:
# -tutti i geni DE della famiglia (line sottile)
# -i geni DE che risultano piu interessanti
library(scales)

#######################################################################################################################

#INPUT DA INSERIRE

sample<-"bHLH selected candidate"

gene_list<-readLines("fam_gene_list/HB_bHLH_all_DE.txt")

TOP_gene_list<-readLines("fam_gene_list/HB_bHLH_in_situ.txt")

notes=""

#######################################################################################################################

#defaul input

sp1<-"ath"
color_sp1<-"black"

sp2<-"osa"
color_sp2<-"blue"

sp3<-"sly"
color_sp3<-"red"

sp4<-"hvu"
color_sp4<-"green3"



count_sp1<-read.table("dataset/count/ath_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp1<-as.matrix(count_sp1)

count_sp2<-read.table("dataset/count/osa_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp2["na_column"]<-NA
count_sp2<-as.matrix(count_sp2)

count_sp3<-read.table("dataset/count/sly_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp3<-as.matrix(count_sp3)

count_sp4<-read.table("dataset/count/hvu_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp4<-as.matrix(count_sp4)

########################################################################################################################

#genero subset delle matrici di conti con solo i geni di gene list, tenendo comunque separate le specie

#lista DE
sub_count_sp1<-subset(count_sp1, rownames(count_sp1) %in% gene_list)
sub_count_sp2<-subset(count_sp2, rownames(count_sp2) %in% gene_list)
sub_count_sp3<-subset(count_sp3, rownames(count_sp3) %in% gene_list)
sub_count_sp4<-subset(count_sp4, rownames(count_sp4) %in% gene_list)

#lista TOP
TOP_sub_count_sp1<-subset(count_sp1, rownames(count_sp1) %in% TOP_gene_list)
TOP_sub_count_sp2<-subset(count_sp2, rownames(count_sp2) %in% TOP_gene_list)
TOP_sub_count_sp3<-subset(count_sp3, rownames(count_sp3) %in% TOP_gene_list)
TOP_sub_count_sp4<-subset(count_sp4, rownames(count_sp4) %in% TOP_gene_list)

#genero la il grafico

pdf(paste("results/",sample,"_trend.pdf", sep = ""))

plot(c(1,2,3),c(0,0.5,1), main = sample, sub = notes , col="white" , xlab = "", ylab="relative expression", axes=FALSE, cex.sub=0.8, font.sub=2 )
title(xlab="Stages", line=2, cex.lab=1 )
axis(1, at=c(1,2,3), labels=c("EARLY","MID","LATE"))
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
legend(1.25,0.95,legend = c(sp1,sp2,sp3,sp4), fill = c(color_sp1,color_sp2,color_sp3,color_sp4))

#stampo tutti i TOP

for (gene in rownames(TOP_sub_count_sp1)){
  lines(TOP_sub_count_sp1[gene,], col=color_sp1, lwd=1)
  text(2.9-(3-ncol(TOP_sub_count_sp1)),TOP_sub_count_sp1[gene,ncol(count_sp1)] ,labels=gene, col = color_sp1, cex = 0.6)
} 

for (gene in rownames(TOP_sub_count_sp2)){
  lines(TOP_sub_count_sp2[gene,], col=color_sp2, lwd=1)
  text(1.9-(3-ncol(TOP_sub_count_sp2)),TOP_sub_count_sp2[gene,ncol(count_sp2)-1] ,labels=gene, col = color_sp2, cex = 0.6)
} 

for (gene in rownames(TOP_sub_count_sp3)){
  lines(TOP_sub_count_sp3[gene,], col=color_sp3, lwd=1)
  text(2.9-(3-ncol(TOP_sub_count_sp3)),TOP_sub_count_sp3[gene,ncol(count_sp3)] ,labels=gene, col = color_sp3, cex = 0.6)
} 

for (gene in rownames(TOP_sub_count_sp4)){
  lines(TOP_sub_count_sp4[gene,], col=color_sp4, lwd=1)
  text(2.9-(3-ncol(TOP_sub_count_sp4)),TOP_sub_count_sp4[gene,ncol(count_sp4)] ,labels=gene, col = color_sp4, cex = 0.6)
} 


#stampo tutti DE

for (gene in rownames(sub_count_sp1)){
  lines(sub_count_sp1[gene,], col=alpha(color_sp1, 0.1), lty=1)
} 

for (gene in rownames(sub_count_sp2)){
  lines(sub_count_sp2[gene,], col=alpha(color_sp2, 0.1), lty=1)
} 

for (gene in rownames(sub_count_sp3)){
  lines(sub_count_sp3[gene,], col=alpha(color_sp3, 0.1), lty=1)
} 

for (gene in rownames(sub_count_sp4)){
  lines(sub_count_sp4[gene,], col=alpha(color_sp4, 0.1), lty=1)
} 

dev.off()

