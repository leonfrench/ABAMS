#Read in
densityMatrix <- as.matrix(read.table("tempDensityMatrix.txt"))
levelMatrix <- as.matrix(read.table("tempLevelMatrix.txt"))
#convert large level values to 1500
levelMatrix[levelMatrix==1e+37] <- 1500

#what is the energy of the large levels?? fix
#energyMatrix[levelMatrix == 1e+37] 

energyTable <- read.table("tempEnergyMatrix.txt")
#energyTable <- read.table("tempEnergyMatrixReduced.txt")

#region classification
regionClassification <- read.csv("RegionClassification.csv", row.names=1)
#make it the variable name setup - a mapping from a region to its parent
row.names(regionClassification) <- make.names(row.names(regionClassification))
regionMatrix <- as.matrix(regionClassification)
#as.matrix(regionMatrix[colnames(energyMatrix),])


#sort
energyTable <- energyTable[order(row.names(energyTable)),]
energyTable <- energyTable[,order(colnames(energyTable))]

#convert to matrix
energyMatrix <- as.matrix(energyTable)

#remove two brain regions
#Posterodorsal.preoptic.nucleus is all zeroes
#Nucleus.y is all zeroes
energyMatrix <- subset(energyMatrix, select=c(-Posterodorsal.preoptic.nucleus,-Nucleus.y))


#convert zeroes to the roughly lowest number thats not zero
energyMatrix[energyMatrix==0] <- .00000000001
energyMatrix <- log(energyMatrix)

hist(energyMatrix, breaks= 500, xlab="log(energy)")

map <- heatmap(cor(energyMatrix), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none")

numColours <- length(unique(regionMatrix))
colour <- matrix(nrow=numColours)
rownames(colour) <- unique(regionMatrix)[,1]
colour[,1] <- rainbow(numColours)
#original scheme
rownames(colour) <- c("Hindbrain", "Endbrain", "Interbrain", "Midbrain")
#parents of all the regions
parents <- as.matrix(regionMatrix[colnames(energyMatrix),])
#convert to colours for all regions
parentCols <- colour[parents,]

map <- heatmap(cor(energyMatrix), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none", ColSideColors=parentCols, margins=c(10,10))

#write heatmap to PNG
png("heatMapReduced.png", width=1640, height=1640) 
#heatmap(cor(energyMatrix), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE)
#show the one with the bar on top
map <- heatmap(cor(energyMatrix), symm = TRUE, distfun = function(c) as.dist(1 - abs(c)), keep.dendro = TRUE, scale = "none", ColSideColors=parentCols, margins=c(10,10))
dev.off()


#reorganize according to region classification
seenRegions <- intersect(rownames(regionMatrix), colnames(energyMatrix))
regionalized <-energyMatrix[,seenRegions]
#colours
parents <- as.matrix(regionMatrix[seenRegions,])
#convert to colours for all regions
parentCols <- colour[parents,]
heatmap(cor(regionalized), scale="no", symm = TRUE, Rowv=NA, Colv=NA, ColSideColors=parentCols,margins=c(10,10))





#show the dendrogram
par(mar=c(17,4,4,2))
plot(map$Colv)

#write histogram to PNG
png("histogram.png") 
hist(energyMatrix, breaks= 500, xlab="log(energy)")
dev.off()




#write dendro to PNG
png("dendro.png", width=3240, height=740) 
par(mar=c(17,4,4,2))
plot(map$Rowv)
dev.off()


#write level histogram to PNG
png("levelHistogram.png") 
hist(log(levelMatrix), breaks= 1500, xlab="log(level)")
dev.off()

#write density histogram to PNG
png("densityHistogram.png") 
hist(log(densityMatrix), breaks= 1500, xlab="log(density)")
dev.off()

#Low correlation with raphe magnus
#which(c==min(c), arr.ind=TRUE)
#plot(energyMatrix[,"Nucleus.raphé.magnus"])
# it has alot of zeroes
#hist(energyMatrix[,"Nucleus.raphé.magnus"])



#percent of data is NaN
length(expressionData[is.nan(expressionData)])/length(expressionData)
length(expressionData[expressionData==0])/length(expressionData)
expressionData[is.nan(expressionData)] <- 0

