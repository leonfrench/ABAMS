energyTable <- read.table("AllenMatrixPair.ExpressionEnergy.txt")

energyTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.txt")

energyMatrix <- as.matrix(energyTable)
dim(energyMatrix)
length(energyMatrix[energyMatrix==0])

png("energyHist.png", width=640, height=640) 
hist(energyMatrix, xlab="Energy")
dev.off()


energyMatrix[is.nan(energyMatrix)] = 0

#assume it was logged, undo
energyMatrix <- exp(energyMatrix)
sdev <- sd(t(energyMatrix), na.rm=TRUE)
rowMeans <- rowMeans(energyMatrix, na.rm = TRUE)
nrmSdev <- sdev/rowMeans
hist(nrmSdev, xlab="Normalized Energy Variance")
plot(sdev, rowMeans)

png("energyNormSd.png", width=640, height=640) 
hist(nrmSdev, xlab="Normalized Energy Variance")
dev.off()

png("energySD.png", width=640, height=640) 
hist(sdev, breaks= 100, xlab="Energy Variance")
dev.off()

# Abtb1[1660]

pr<-prcomp(t(energyMatrix))


png("AllenMatrixPair.ExpressionEnergy.Distance.forPlotting.png", width=640, height=640) 
plot(x)
dev.off()
