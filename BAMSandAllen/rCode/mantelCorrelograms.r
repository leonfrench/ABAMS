library(mantel.correlog)

####
energies <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.Correlation.txt"))
energiesCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.Residuals(EuclidianDistance).Correlation.txt"))
connections <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Correlation.txt"))
connectionsCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Residuals(EuclidianDistance).Correlation.txt"))
distance <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.EuclidianDistance.Explain.Adjacency.txt"))

###log stuff
#logDistance <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.LogEuclidianDistance.Explain.Correlation.txt"))
energiesLogCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.Residuals(LogEuclidianDistance).Correlation.txt"))
connectionsLogCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Residuals(LogEuclidianDistance).Correlation.txt"))

energies <- -1 * energies
energiesCorrected <- -1 * energiesCorrected
energiesLogCorrected <- -1 * energiesLogCorrected

connections <- -1 * connections
connectionsCorrected <- -1 * connectionsCorrected
connectionsLogCorrected <- -1 * connectionsLogCorrected

plot(mantel.correlog(connections, distance))
plot(mantel.correlog(connectionsLogCorrected, distance))

plot(mantel.correlog(energies, distance))
plot(mantel.correlog(energiesLogCorrected, distance), nperm=100)

x<-mantel.correlog(connections, distance)
plot(x$mantel.res[,3])
x<-mantel.correlog(connectionsLogCorrected, distance)
lines(x$mantel.res[,3])


png("energyMantelCorrelog.png", width=640, height=640) 
plot(mantel.correlog(energies, distance))
dev.off()

png("energyMantelCorrectedCorrelog.png", width=640, height=640) 
plot(mantel.correlog(energiesCorrected, distance))
dev.off()

png("energyMantelLogCorrectedCorrelog.png", width=640, height=640) 
plot(mantel.correlog(energiesLogCorrected, distance))
dev.off()


png("connectionsMantelCorrelog.png", width=640, height=640) 
plot(mantel.correlog(connections, distance))
dev.off()


png("connectionsMantelCorrectedCorrelog.png", width=640, height=640) 
plot(mantel.correlog(connectionsCorrected, distance))
dev.off()

png("connectionsMantelLogCorrectedCorrelog.png", width=640, height=640) 
plot(mantel.correlog(connectionsLogCorrected, distance))
dev.off()


plot(mantel.correlog(connections, distance))
plot(mantel.correlog(energies, distance))
plot(mantel.correlog(energiesCorrected, distance))
plot(mantel.correlog(connectionsCorrected, distance))

###########################################################
#overlay code
#From http://www.bio.umontreal.ca/legendre/indexEn.html
addMantel <- function(x, color="black", alpha=0.05, ...)
{
lim = max(x$n.tests)
if(x$mult=="none") {
	signif = which((x$mantel.res[1:lim,4] <= alpha))
	} else {
	signif = which((x$mantel.res[1:lim,5] <= alpha))
	}
lines(x$mantel.res[1:lim,1], x$mantel.res[1:lim,3], col=color)
points(x$mantel.res[1:lim,1], x$mantel.res[1:lim,3], pch=22, col=color, bg="white")
points(x$mantel.res[signif,1], x$mantel.res[signif,3], pch=22, col=color, bg=color)
}
###########################################################

png("connectionsMantelCorrelog.png", width=640, height=640) 
plot(mantel.correlog(connections, distance))
title("Connectivity Mantel Correlogram")
addMantel(mantel.correlog(connectionsCorrected, distance), color="blue")
addMantel(mantel.correlog(connectionsLogCorrected, distance), color="red")
legend("topright", c("Normal", "Linear Corrected", "Log Transform Corrected"), fill=c("black", "blue", "red"))
dev.off()

png("energiesMantelCorrelog.png", width=640, height=640) 
plot(mantel.correlog(energies, distance))
title("Expression Mantel Correlogram")
addMantel(mantel.correlog(energiesCorrected, distance), color="blue")
addMantel(mantel.correlog(energiesLogCorrected, distance), color="red")
legend("topright", c("Normal", "Linear Corrected", "Log Transform Corrected"), fill=c("black", "blue", "red"))
dev.off()
