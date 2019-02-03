#script che prende come input la tabella di plaza con le famiglie di ortho, le liste di geni significtivi risultanti da NPMI e 
#le tabelle di conti normalizai in modo da avere somma della riga 1
#lo script restituisce un grafico per ogni famiglia con all'interno l'andamento dei geni arricchiti nei tre stadi per le due specie 
#e una heatmap che clasterizza sui trend dei geni rrichiti in ogni famiglia. 

#input da inserire################################################################################################

plaza_tab<-read.table("dataset/genefamily_data_ath_sly_osa_hvu_plaza.txt", row.names = 3)

sample<-"IM_ST3_vs_ind_det"

sp1<-"ath"
color_sp1<-"black"
sp1_genes<-readLines("dataset/IM_ST3/IM_ST3_vs_IND_DET_all_significant_genes_sp1.txt")

sp2<-"osa"
color_sp2<-"green"
sp2_genes<-readLines("dataset/IM_ST3/IM_ST3_vs_IND_DET_all_significant_genes_sp2.txt")

count_sp1<-read.table("lib_size_norm_count/norm_to_one_count/ath_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp1<-as.matrix(count_sp1)
count_sp2<-read.table("lib_size_norm_count/norm_to_one_count/osa_norm_to_one_count.txt", header = T, row.names = 1, sep = "\t")
count_sp2<-as.matrix(count_sp2)

##################################################################################################################

#tabella di plaza ma con solo i geni geni delle due liste

plaza_only_en<-as.matrix(plaza_tab[c(sp1_genes,sp2_genes),])

#lista di famiglie (unique) presenti nelle due liste di geni 

en_fam<-as.character(unique(plaza_only_en[,1]))

number_of_genes<-c(1:nrow(plaza_only_en))

#genero la lista my_data in cui tramite il ciclo annoto famiglia e relativi geni delle due specie 

my_data<-data.frame(matrix( ncol = 2))
colnames(my_data)<-c("sp1","sp2")

#indice per aggiungere la famiglia e i geni a my data
index<-1

for (family in en_fam){
  
  #liste vuote i cui aggiungere geni di sp1 e sp2 per ogni famiglia
  sp1_g<-character()
  sp2_g<-character()
  
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
    }
  }
  
  if (length(sp1_g)!=0 && length(sp2_g)!=0){
    my_data[family,1]<-paste(sp1_g, collapse = ',')
    my_data[family,2]<-paste(sp2_g, collapse = ',')
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
  
  plot(c(1,2,3),c(0,0.5,1), main = fam, col="white", xlab = "stages" , ylab="relative expression", axes=FALSE )
  axis(1, at=c(1,2,3), labels=colnames(count_sp1))
  axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
  legend(1.25,0.95,legend = c(sp1,sp2), fill = c(color_sp1,color_sp2))
  
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
}




dev.off()
