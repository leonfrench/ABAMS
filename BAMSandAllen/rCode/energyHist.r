energyTable <- read.table("tempEnergyMatrix.txt")
#energyTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.txt")
#energyTable <- read.table("ConnectivityAndAllenExpressionPair.MatrixB.txt")
energyTable <- energyTable[order(row.names(energyTable)),]
energyTable <- energyTable[,order(colnames(energyTable))]
energyMatrix <- as.matrix(energyTable)
energyMatrix <- subset(energyMatrix, select=c(-Posterodorsal.preoptic.nucleus,-Nucleus.y))
dim(energyMatrix)
length(energyMatrix[energyMatrix==0])
hist(energyMatrix, breaks= 500, xlab="Energy")

png("energyHist.png", width=640, height=640) 
hist(energyMatrix, breaks= 500, xlab="Energy")
dev.off()

png("energyZeroHist.png", width=640, height=640) 
hist(energyMatrix, breaks= 150000, xlab="Energy", xlim=c(0,1))
dev.off()

png("logEnergyHist.png", width=640, height=640) 
hist(log(energyMatrix), breaks= 500, xlab="log(Energy)")
dev.off()

