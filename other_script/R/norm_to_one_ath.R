#breve script per normalizzare i conti in modo da avere le media delle repliche e il peso di ogni stadio (somma della riga = 1)

#tabella con conti normalizzati
norm_cts<-read.table("norm_count/ath_norm_count.txt", row.names = 1, header = T, sep="\t")

#tabella in cui annotare i risultati
norm_to_one<-data.frame(row.names = row.names(norm_cts))

#media delle repliche
norm_to_one$IM<-rowMeans(norm_cts[,c(1,2)])
norm_to_one$FM<-rowMeans(norm_cts[,c(3,4)])
norm_to_one$ST3<-rowMeans(norm_cts[,c(4,6)])

#divisione pper soma riga 
mat<-as.matrix(norm_to_one)
matnorm<-mat/rowSums(mat) 
norm_to_one<-as.data.frame(matnorm) 

write.table(norm_to_one , "norm_to_one_count/ath_norm_to_one_count.txt", sep = "\t", quote = F)