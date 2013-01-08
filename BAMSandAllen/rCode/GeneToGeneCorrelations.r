energyOptTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.out.topList.txt")
energyOptMatrixOut <- as.matrix(energyOptTable)
dim(energyOptMatrixOut)
length(energyOptMatrixOut[energyOptMatrixOut==0])
energyOptMatrixOutGoingT <- as.matrix(t(energyOptMatrixOut))
OptCorrelations <- cor(energyOptMatrixOutGoingT, method = "pearson", use = "pairwise")
OptCorrelations[lower.tri(OptCorrelations, diag=TRUE)] <- NA
mean(OptCorrelations, na.rm=TRUE)
plot(density(OptCorrelations, na.rm=TRUE))

energyOptTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.in.toplist.txt")
energyOptMatrixIn <- as.matrix(energyOptTable)
dim(energyOptMatrixIn)
length(energyOptMatrixIn[energyOptMatrixIn==0])
energyOptMatrixInT <- as.matrix(t(energyOptMatrixIn))
OptCorrelationsIn <- cor(energyOptMatrixInT, method = "pearson", use = "pairwise")
OptCorrelationsIn[lower.tri(OptCorrelationsIn, diag=TRUE)] <- NA
mean(OptCorrelationsIn, na.rm=TRUE)

energyOptTable <- read.table("AllenMatrixPair.NewEnergies.space.txt")
energyOptMatrixSpace <- as.matrix(energyOptTable)
dim(energyOptMatrixSpace)
length(energyOptMatrixSpace[energyOptMatrixSpace==0])
energyOptMatrixSpaceT <- as.matrix(t(energyOptMatrixSpace))
OptCorrelationsSpace <- cor(energyOptMatrixSpaceT, method = "pearson", use = "pairwise")
OptCorrelationsSpace[lower.tri(OptCorrelationsSpace, diag=TRUE)] <- NA
mean(OptCorrelationsSpace, na.rm=TRUE)



energyTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.out.txt")
energyMatrixOut <- as.matrix(energyTable)
dim(energyMatrixOut)
length(energyMatrixOut[energyMatrixOut==0])
energyMatrixOutT <- as.matrix(t(energyMatrixOut))
Correlations <- cor(energyMatrixOutT, method = "pearson", use = "pairwise")
Correlations[lower.tri(Correlations)] <- NA
mean(Correlations, na.rm=TRUE)

energyTable <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.in.txt")
energyMatrixIn <- as.matrix(energyTable)
dim(energyMatrixIn)
length(energyMatrixIn[energyMatrixIn==0])
energyMatrixInT <- as.matrix(t(energyMatrixIn))
CorrelationsIn <- cor(energyMatrixInT, method = "pearson", use = "pairwise")
CorrelationsIn[lower.tri(CorrelationsIn)] <- NA
mean(CorrelationsIn, na.rm=TRUE)








# plot
postscript("AllGeneToGeneOut.ps")
png("AllGeneToGeneOut.png", width=800, height=800) 
plot(density(Correlations, na.rm=TRUE, kernel = "rectangular"), lty="solid", main="Density Plot of Expression Correlation Between Image Series Pairs", xlab="Correlation")
lines(density(OptCorrelations, na.rm=TRUE, kernel = "rectangular"), lty="dashed")
legend("topleft", c("All genes", "Optimized outgoing genes"), lty=c("solid", "dashed"))
dev.off()

#ks.test(Correlations, OptCorrelations)
#t.test(Correlations, OptCorrelations)


corVector <- as.vector(Correlations)
corVector <- corVector[!is.na(corVector)]

mean(Correlations, na.rm=TRUE)
sd(as.vector(Correlations), na.rm=TRUE)
var(Correlations, na.rm=TRUE)

optCorVector<-as.vector(OptCorrelations)
optCorVector<- optCorVector[!is.na(optCorVector)]
mean(OptCorrelations, na.rm=TRUE)
sd(optCorVector)

t.test(corVector, optCorVector)
#ks.test(corVector, optCorVector)




# plot
postscript("AllGeneToGeneIn.ps")
png("AllGeneToGeneInSpace.png", width=800, height=800) 
plot(density(CorrelationsIn, na.rm=TRUE, kernel = "rectangular"), lty="solid", main="Density Plot of Expression Correlation Between Image Series Pairs", xlab="Correlation")
lines(density(OptCorrelationsIn, na.rm=TRUE, kernel = "rectangular"), lty="dashed")
lines(density(OptCorrelationsSpace, na.rm=TRUE, kernel = "rectangular"), lty="dotted")
legend("topleft", c("All genes", "Optimized incoming genes", "Proximity genes"), lty=c("solid", "dashed", "dotted"))
dev.off()


