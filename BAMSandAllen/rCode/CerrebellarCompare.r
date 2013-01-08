x <- read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Residuals(LogEuclidianDistance).Correlation.txt")
x <- as.matrix(x)

png("cerebellarCortex2.png", width=640, height=640) 
#boxplot(x["Cerebellar cortex",], x,names=c("Cerebellar cortex", "All regions"))
boxplot(x["Cochlear nuclei",],x["Cerebellar cortex",], x,names=c("Cochlear.nuclei","Cerebellar cortex", "All regions"))
dev.off()

x <- read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Residuals(LogEuclidianDistance).Correlation.out.txt")


rows <- apply(x, 2, mean)
min(rows)
rows["Cerebellar.cortex"]

rows <- apply(x, 2, quantile)

rows <- apply(x, 2, quantile, probs=0.5)
min(rows)
which( min(rows) == rows)
rows["Cerebellar.cortex"]

#ranks, relative
ranks <- apply(x,2,rank)
ranks["Cerebellar.cortex"]

ranks["Cerebellar cortex",] < 112/4

#how many regions have the cerebellar cortex in it's bottom quartile of regions of similar connectivity
length(which(ranks["Cerebellar cortex",] < 112/4))

length(which(ranks["Cerebellar cortex",] < 112/10))


