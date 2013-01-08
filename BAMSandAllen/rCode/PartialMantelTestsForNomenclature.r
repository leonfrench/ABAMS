library(vegan)
library(mantel.correlog)

#### last ran from /grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/mantel.correlog
energies <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.Correlation.txt"))
connections <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Correlation.txt"))
distance <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.EuclidianDistance.Explain.Adjacency.txt"))
nomen <- as.matrix(read.table("ConnectivityAndAllenNomenclaturePair.Nomenclature.Correlation.txt"))

energiesLogCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.Residuals(LogEuclidianDistance).Correlation.txt"))
connectionsLogCorrected <- as.matrix(read.table("ConnectivityAndAllenPartialExpressionMatrixPair.Connectivity.Residuals(LogEuclidianDistance).Correlation.txt"))


nnames <- row.names(nomen)
dnames <-row.names(distance)
cnames <- row.names(connections)
enames <- row.names(energies)
cor(rank(nnames), rank(dnames))
cor(rank(nnames), rank(cnames))
cor(rank(nnames), rank(enames))

mantel.partial(energies, connections, log(distance))
mantel.partial(energies, nomen, log(distance))
mantel.partial(energies, -1*nomen, log(distance))
mantel.partial(connections, nomen, log(distance))



mantel(energiesLogCorrected, -1*nomen)
mantel(connectionsLogCorrected, nomen)