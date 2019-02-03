#breve script per normalizzare i conti
library("DESeq2")

cts<-read.table("NOT_norm/osa_count.txt", row.names = 1, header = T, sep="\t")
coldata<-read.table("NOT_norm/osa_coldata.txt", row.names = 1, header = T, sep="\t")

cts<-round(cts,digits = 0)
cts<-cts[rowMeans(cts)>10,]

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~Condition)
dds <- estimateSizeFactors(dds)

normCounts <- counts(dds, normalized=TRUE)
normCounts <- round(normCounts, digits = 2)

write.table(normCounts, "norm_count/osa_norm_count.txt", sep = "\t", quote = F)