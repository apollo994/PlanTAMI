library(edgeR)
#carica la tabella con i conti di tutti i campioni
data<-read.table("all_count_ath.txt",row.names=1,sep="\t",header=TRUE)
#carica la tabella con i nomi dei campioni e la relativa condizione
design<-read.table("coldata_ara.txt",row.names=1)
#carica una matrice di fattori con campioni e condizioni e relativi 0 e 1
EXPdesign<-read.table("EXP_design.txt",header = T,row.names = 1,sep="\t")
#solglia con cui chiamare i significativi
th=0.05

#elenco dei campioni specificando il nome della condzione
y <- DGEList(counts=data,group=c("IM","IM","FM","FM","ST3","ST3"))
y=calcNormFactors(y)
#filtro sul numero di conti
keep <- rowSums(cpm(y)>2) >= 10
y <- y[keep, keep.lib.sizes=FALSE]
y <- estimateDisp(y,EXPdesign)
fit <- glmQLFit(y, EXPdesign)


contrasts<-makeContrasts(IMvsFM=IM-FM, IMvsST3=IM-ST3, FMvsST3=FM-ST3, levels=EXPdesign)




for (i in 1:ncol(contrasts))
{
  res<-glmQLFTest(fit, contrast=contrasts[,i])
  tableR<-res$table
  FDR<-p.adjust(tableR[,4],method="BH")
  tableR<-cbind(tableR,FDR)
  tableR<-tableR[order(tableR$FDR),]
  
  #stapo la tabella con tutti i risultati
  name<-(colnames(contrasts))[i]
  file<-paste(name,"res.txt",sep=".")
  write.table(tableR,file=file, sep="\t", quote = F)
  
  #stampo la tabella con solo i risultati significativi (th specificato all'inizio)
  soglia<-paste("sig",th,sep = "_")
  file_sig<-paste(soglia,file,sep = "_")
  to_keep<-tableR$FDR<=th
  tableR_sig<-tableR[to_keep,]
  tableR_sig<-tableR_sig[order(tableR_sig$FDR),]
  write.table(tableR_sig,file=file_sig, sep="\t", quote = F)
}


