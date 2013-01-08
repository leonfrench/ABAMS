x <- read.table("ConnectivityAndAllenExpressionMatrixPair.GeneToRegion.jacknife.txt")
dim(x)
matX <- as.matrix(x)
hist(matX)
density(matX, na.rm=TRUE)

plot(density(matX, na.rm=TRUE))
max(matX)
max<-max(matX, na.rm=TRUE)
min<-min(matX, na.rm=TRUE)
which(matX == min, arr.ind=TRUE)
which(matX == max, arr.ind=TRUE)

which(is.nan(matX), arr.ind=TRUE)


cords <- matrix(0, 2, dim(matX)[2])
cords[1,] <- matX[10178,]
cords[2,] <- 0
points(t(cords), col="blue")
