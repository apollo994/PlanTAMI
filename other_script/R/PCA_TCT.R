#dati su cui calcolare la PCA
data<-read.table("norm_count_all.txt", header = T, row.names = 1)
#colori relativi ai campioni
color<-c("red","yellow","green")


project.pca <- prcomp(t(data), scale = T)
summary(project.pca)

#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

pdf("PCA_TCT.pdf")

barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,1:5], col=color, main="Principal components analysis bi-plot\nPCs 1-5", pch=16)


par(xpd = T, mar = par()$mar + c(0,0,0,7))
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col=color, pch=16, cex=1)
legend(45,30,legend = colnames(data) , fill = color)
par(mar=c(5, 4, 4, 2) + 0.1)

dev.off()

