# run in data/correlation study

energyOptTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.out.topList.txt")
energy <- as.matrix(energyOptTable)
dim(energy)
length(energy[energy==0])

#region classification
regionClassification <- read.csv("RegionClassificationABAMSOut.csv", row.names=1)
#make it the variable name setup - a mapping from a region to its parent
row.names(regionClassification) <- make.names(row.names(regionClassification))
regionMatrix <- as.matrix(regionClassification)
#as.matrix(regionMatrix[colnames(energyMatrix),])

map <- heatmap(cor(energy, use="pairwise"), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none")


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

map <- heatmap(cor(energy, use="pairwise"), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none", ColSideColors=parentCols, margins=c(10,10))

