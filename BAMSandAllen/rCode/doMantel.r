
expressionData <- as.matrix(read.table("ConnectivityAndAllenMinerExpressionPair.MatrixB.txt"))
connectionData <- as.matrix(read.table("ConnectivityAndAllenMinerExpressionPair.MatrixA.txt"))

expressionData <- as.matrix(read.table("ConnectivityAndAllenExpressionPair.MatrixB.txt"))
connectionData <- as.matrix(read.table("ConnectivityAndAllenExpressionPair.MatrixA.txt"))

rexpressionCorrelations <- cor(expressionData, method = "pearson")
rconnectionCorrelations<- cor(connectionData, method = "pearson")

rexpressionCorrelations <- cor(expressionData, method = "spearman")
#don't use spearman on connections - its binary
rconnectionCorrelations<- cor(connectionData, method = "pearson")

dim(rexpressionCorrelations)
max(rexpressionCorrelations)
min(rexpressionCorrelations)
dim(rconnectionCorrelations)
max(rconnectionCorrelations)
min(rconnectionCorrelations)
hist(rexpressionCorrelations, breaks=500)
hist(rconnectionCorrelations, breaks=500)

SE0900022

library("vegan")
mantel(rexpressionCorrelations, rconnectionCorrelations)


#part two

library("vegan")
expressionCorrelations <- as.matrix(read.table("ConnectivityAndAllenExpressionPair.B.Correlation.txt"))
connectionCorrelations <- as.matrix(read.table("ConnectivityAndAllenExpressionPair.A.Correlation.txt"))

hist(expressionCorrelations, breaks=500)
hist(connectionCorrelations, breaks=500)

mantel(expressionCorrelations, connectionCorrelations)
mantel(connectionCorrelations, expressionCorrelations)
mantel(expressionCorrelations, connectionCorrelations, method="spearman")
mantel(expressionCorrelations, connectionCorrelations, permutations=10000)
