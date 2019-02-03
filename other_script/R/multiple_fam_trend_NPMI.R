#questo script rispetto a fam_trend_NPMI da la possibilità di vedere l'andamento contemporaneamenti di 4 specie 

#input da inserire################################################################################################

plaza_tab<-read.table("dataset/genefamily_data_ath_sly_osa_hvu_plaza.txt", row.names = 3)

sample<-"early_to_mid"

sp1<-"ath"
color_sp1<-"black"
sp1_genes<-readLines("dataset/early_to_mid/ath_IM_FM.txt")

sp2<-"osa"
color_sp2<-"green"
sp2_genes<-readLines("dataset/early_to_late/osa_ind_det.txt")

sp3<-"sly"
color_sp3<-"red"
sp3_genes<-readLines("dataset/early_to_mid/sly_TM_SIM.txt")

sp4<-"hvu"
color_sp4<-"purple"
sp4_genes<-readLines("dataset/early_to_mid/hvu_1.0_2.0.txt")



count_sp1<-read.table("lib_size_norm_count/norm_to_one_count/ath_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp1<-as.matrix(count_sp1)

count_sp2<-read.table("lib_size_norm_count/norm_to_one_count/osa_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp2<-as.matrix(count_sp2)

count_sp3<-read.table("lib_size_norm_count/norm_to_one_count/sly_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp3<-as.matrix(count_sp3)

count_sp4<-read.table("lib_size_norm_count/norm_to_one_count/hvu_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp4<-as.matrix(count_sp4)

##################################################################################################################

#tabella di plaza ma con solo i geni delle quattro liste

plaza_only_en<-as.matrix(plaza_tab[c(sp1_genes,sp2_genes,sp3_genes,sp4_genes),])
plaza_only_en<-na.omit(plaza_only_en)

#lista di famiglie (unique) presenti nelle quattro liste di geni 

en_fam<-as.character(unique(plaza_only_en[,1]))

number_of_genes<-c(1:nrow(plaza_only_en))

#genero la lista my_data in cui tramite il ciclo annoto famiglia e relativi geni delle quattro specie 

my_data<-data.frame(matrix( ncol = 4))
colnames(my_data)<-c("sp1","sp2","sp3","sp4")

#indice per aggiungere la famiglia e i geni a my data
index<-1

for (family in en_fam){
  
  #liste vuote i cui aggiungere geni di sp1 e sp2 per ogni famiglia
  sp1_g<-character()
  sp2_g<-character()
  sp3_g<-character()
  sp4_g<-character()
  
  #ciclo sulla tabella di plaza con solo i geni delle due liste
  for (i in number_of_genes){
    
    #print (plaza_only_en[i,])
    
    if (plaza_only_en[i,1]==family){
      if (plaza_only_en[i,2]==sp1){
        sp1_g<-append(sp1_g,rownames(plaza_only_en)[i])
        #print(rownames(plaza_only_en)[i])
      }
      if (plaza_only_en[i,2]==sp2){
        sp2_g<-append(sp2_g,rownames(plaza_only_en)[i])
        #print(rownames(plaza_only_en)[i])
      }
      if (plaza_only_en[i,2]==sp3){
        sp3_g<-append(sp3_g,rownames(plaza_only_en)[i])
        #print(rownames(plaza_only_en)[i])
      }
      if (plaza_only_en[i,2]==sp4){
        sp4_g<-append(sp4_g,rownames(plaza_only_en)[i])
        #print(rownames(plaza_only_en)[i])
      }
    }
  }
  
  if (length(sp1_g)!=0 && length(sp2_g)!=0 && length(sp3_g)!=0 && length(sp4_g)!=0){
    my_data[family,1]<-paste(sp1_g, collapse = ',')
    my_data[family,2]<-paste(sp2_g, collapse = ',')
    my_data[family,3]<-paste(sp3_g, collapse = ',')
    my_data[family,4]<-paste(sp4_g, collapse = ',')
  }
}

my_data<-na.omit(my_data)

rm(plaza_tab)

#generazione dei grafici#######################################################################################################

file_name<-paste("results/",sample,"_fam_trend.pdf", sep = "")

pdf(file_name)

for (fam in rownames(my_data)){
  sp1_g<-strsplit(my_data[fam,1], ",")
  sp2_g<-strsplit(my_data[fam,2], ",")
  sp3_g<-strsplit(my_data[fam,3], ",")
  sp4_g<-strsplit(my_data[fam,4], ",")
  
  plot(c(1,2,3),c(0,0.5,1), main = fam, col="white", xlab = "stages" , ylab="relative expression", axes=FALSE )
  axis(1, at=c(1,2,3), labels=colnames(count_sp1))
  axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
  legend(1.25,0.95,legend = c(sp1,sp2,sp3,sp4), fill = c(color_sp1,color_sp2,color_sp3,color_sp4))
  
  for (genes in sp1_g){
    for (gene in genes){
      if (gene %in% rownames(count_sp1)){
        lines(count_sp1[gene,], col=color_sp1)
        text(2.9-(3-ncol(count_sp1)),count_sp1[gene,ncol(count_sp1)] ,labels=gene, col = color_sp1, cex = 0.5) 
      }
    }
  }
  
  for (genes in sp2_g){
    for (gene in genes){
      if (gene %in% rownames(count_sp2)){
        lines(count_sp2[gene,], col=color_sp2)
        text(2.9-(3-ncol(count_sp2)),count_sp2[gene,ncol(count_sp2)] ,labels=gene, col = color_sp2, cex = 0.5)
      }
    }
  }
  
  for (genes in sp3_g){
    for (gene in genes){
      if (gene %in% rownames(count_sp3)){
        lines(count_sp3[gene,], col=color_sp3)
        text(2.9-(3-ncol(count_sp3)),count_sp3[gene,ncol(count_sp3)] ,labels=gene, col = color_sp3, cex = 0.5)
      }
    }
  }
  
  for (genes in sp4_g){
    for (gene in genes){
      if (gene %in% rownames(count_sp4)){
        lines(count_sp4[gene,], col=color_sp4)
        text(2.9-(3-ncol(count_sp4)),count_sp4[gene,ncol(count_sp4)] ,labels=gene, col = color_sp4, cex = 0.5)
      }
    }
  }
}




dev.off()

