#questo script prende come input liste di geni, i conti noralzzati a uno uniti tutti insieme e la tabella delle famiglie. 
#ne restituisce una heatmap per ogni famiglia

#NON FUNZIONA CON RISO MA SOLO CON CAMPIONI CON 3 CONFRONTI

#input da inserire################################################################################################

plaza_tab<-read.table("dataset/genefamily_data_ath_sly_osa_hvu_plaza.txt", row.names = 3)

sample<-"mid_to_late"

sp1_genes<-readLines("dataset/mid_to_late/ath_FM_ST3.txt")

sp2_genes<-readLines("dataset/mid_to_late/hvu_2.0_3.5.txt")

sp3_genes<-readLines("dataset/mid_to_late/sly_SIM_FM.txt")


count_sp1<-read.table("lib_size_norm_count/norm_to_one_count/ath_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp1<-as.matrix(count_sp1)

count_sp2<-read.table("lib_size_norm_count/norm_to_one_count/hvu_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp2<-as.matrix(count_sp2)

count_sp3<-read.table("lib_size_norm_count/norm_to_one_count/sly_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp3<-as.matrix(count_sp3)

##################################################################################################################

#tabella di plaza ma con solo i geni delle tre liste

plaza_only_en<-as.matrix(plaza_tab[c(sp1_genes,sp2_genes,sp3_genes),])
plaza_only_en<-na.omit(plaza_only_en)

remove(plaza_tab)

#lista di famiglie (unique) presenti nelle tre liste di geni 

en_fam<-as.character(unique(plaza_only_en[,1]))

#tabella con i conti di tutte e 3 le specie

all_count<-rbind(count_sp1,count_sp3,count_sp4)

#ciclo per fare le heatmap

file_name<-paste("results/",sample,"_heatmap.pdf", sep = "")

pdf(file_name)

for (fam in en_fam){
  #print (fam)
  fam_genes<-rownames(plaza_only_en)[plaza_only_en[,1]==fam]
  if (length(fam_genes)>10){
    heatmap(as.matrix(all_count[rownames(all_count) %in% fam_genes,]), main = fam, margins = c(10,10)) 
  }
}

dev.off()

