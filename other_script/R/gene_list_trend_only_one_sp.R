#scipt che prende come input il una lista di geneID di ath, sly, hvu e osa e partendo dai conti normalizzati 
#tra 0 e 1 ne genera un grafico

#######################################################################################################################

#LISTA INPUT DA INSERIRE

sample<-"LOC_Os09g29360 expression trend"

gene_list<-readLines("fam_gene_list/HB_bHLH_osa.txt")

notes=""

#lista degli input di default 

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

#sub_count_sp1<-count_sp1[rownames(count_sp1) %in% gene_list,]
#sub_count_sp2<-count_sp2[rownames(count_sp2) %in% gene_list,]
#sub_count_sp3<-count_sp3[rownames(count_sp3) %in% gene_list,]
#sub_count_sp4<-count_sp4[rownames(count_sp4) %in% gene_list,]
sub_count_sp1<-subset(count_sp1, rownames(count_sp1) %in% gene_list)
sub_count_sp2<-subset(count_sp2, rownames(count_sp2) %in% gene_list)
sub_count_sp3<-subset(count_sp3, rownames(count_sp3) %in% gene_list)
sub_count_sp4<-subset(count_sp4, rownames(count_sp4) %in% gene_list)


#genero la il grafico

pdf(paste("results/",sample,"_trend.pdf", sep = ""))

plot(c(1,2,3),c(0,0.5,1), main = sample, sub = notes , col="white" , xlab = "", ylab="relative expression", axes=FALSE, cex.sub=0.8, font.sub=2 )
title(xlab="Stages", line=2, cex.lab=1 )
axis(1, at=c(1,2,3), labels=c("EARLY","MID","LATE"))
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
legend(1.25,0.95,legend = c(sp1,sp2,sp3,sp4), fill = c(color_sp1,color_sp2,color_sp3,color_sp4))

for (gene in rownames(sub_count_sp1)){
  lines(sub_count_sp1[gene,], col=color_sp1)
  text(2.9-(3-ncol(sub_count_sp1)),sub_count_sp1[gene,ncol(count_sp1)] ,labels=gene, col = color_sp1, cex = 0.6)
} 

for (gene in rownames(sub_count_sp2)){
  lines(sub_count_sp2[gene,], col=color_sp2)
  text(1.9-(3-ncol(sub_count_sp2)),sub_count_sp2[gene,ncol(count_sp2)-1] ,labels=gene, col = color_sp2, cex = 0.6)
} 

for (gene in rownames(sub_count_sp3)){
  lines(sub_count_sp3[gene,], col=color_sp3)
  text(2.9-(3-ncol(sub_count_sp3)),sub_count_sp3[gene,ncol(count_sp3)] ,labels=gene, col = color_sp3, cex = 0.6)
} 

for (gene in rownames(sub_count_sp4)){
  lines(sub_count_sp4[gene,], col=color_sp4)
  text(2.9-(3-ncol(sub_count_sp4)),sub_count_sp4[gene,ncol(count_sp4)] ,labels=gene, col = color_sp4, cex = 0.6)
} 

dev.off()

