#run in data/correlation study
#energy <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.txt")
energy <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.out.topList.txt")
energy <- as.matrix(energy)
energy <- energy[-197,]


#png("energy.partialcon.png", width=4640, height=4640) 
#heatmap(energy, margins=c(14,14), scale="none", cexRow = .75, cexCol = .75)
#dev.off()

#png("energy.scaled.partialcon.png", width=4640, height=4640) 
#heatmap(energy, margins=c(14,14), scale="row", cexRow = .75, cexCol = .75, na.rm=TRUE)
#dev.off()


xx <= energy[rowSums(!is.na(energy))!=0, colSums(!is.na(energy))!=0] 

ht <- heatmap(energy, method = "average", keep.dendro = TRUE)
orderedRows <- rownames(energy)[ht$rowInd]
#write.table(orderedRows, "ConnectivityAndAllenExpressionMatrixPair.NewEnergies.rOrder", quote=FALSE, row.names = FALSE, col.names = FALSE)
energy <- energy[ht$rowInd, ht$colInd]
write.table(energy, "Rclustered.matrix.txt", quote=FALSE, sep="\t")


#region classification
regionClassification <- read.csv("RegionClassification.csv", row.names=1)
#make it the variable name setup - a mapping from a region to its parent
row.names(regionClassification) <- make.names(row.names(regionClassification))
regionMatrix <- as.matrix(regionClassification)
#as.matrix(regionMatrix[colnames(energyMatrix),])

numColours <- length(unique(regionMatrix))
colour <- matrix(nrow=numColours)
rownames(colour) <- unique(regionMatrix)[,1]
colour[,1] <- rainbow(numColours)
#original scheme
rownames(colour) <- c("Hindbrain", "Endbrain", "Interbrain", "Midbrain")
#parents of all the regions
parents <- as.matrix(regionMatrix[colnames(energy),])
#convert to colours for all regions
parentCols <- colour[parents,]

ht <- heatmap(energy, method = "average"), ColSideColors=parentCols, margins=c(14,14), scale="row", cexRow = .75, cexCol = .75, na.rm=TRUE))

png("outgoing with divisions.png", width=4640, height=4640) 
heatmap(energy, method = "average",ColSideColors=parentCols, margins=c(14,14), scale="row", cexRow = .75, cexCol = .75)
dev.off()

orderedRows <- rownames(energy)[ht$rowInd]
#write.table(orderedRows, "ConnectivityAndAllenExpressionMatrixPair.NewEnergies.rOrder", quote=FALSE, row.names = FALSE, col.names = FALSE)
energy <- energy[ht$rowInd, ht$colInd]
write.table(energy, "Rclustered.matrix.txt", quote=FALSE, sep="\t")


map <- heatmap(cor(energy), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none", ColSideColors=parentCols, margins=c(10,10))

which(row.names(energy)=="Runx2[1278]")

