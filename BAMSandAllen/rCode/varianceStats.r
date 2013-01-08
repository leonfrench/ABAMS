outgoingAll <- read.table("/grp/java/workspace/BAMSandAllen/data/topten/AllGenesOutgoingInfo.geneInfo.csv", header = TRUE)
incomingAll <- read.table("/grp/java/workspace/BAMSandAllen/data/topten/AllGenesIncomingInfo.geneInfo.csv", header = TRUE)
hist(incomingAll[, "sampleStandardDeviation"])

inTopGenes <- read.table("/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final nobed/LOOGenesInOrder.in.partialcon.txt.410.0.018005.topGenes.txt.geneInfo.IN.csv", header = TRUE)
outTopGenes <- read.table("/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final nobed/LOOGenesInOrder.out.partialcon.txt.329.0.014448.topGenes.txt.geneInfo.OUT.csv", header = TRUE)
shuffledOutTopGenes <- read.table("/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final nobed random exp/LOOGenesInOrder.random.out.partialcon.expShuffle.2.txt.256.0.011242.topGenes.txt.geneInfo.csv", header = TRUE)

#in vrs out
t.test(outTopGenes[, "sampleStandardDeviation"], inTopGenes[, "sampleStandardDeviation"])
t.test(outTopGenes[, "expSum"], inTopGenes[, "expSum"])


x <- setdiff(outgoingAll[,"Name"], outTopGenes[,"Name"])
outgoingDiff <- outgoingAll[incomingAll[,"Name"] %in% x, ]

#outgoing
t.test(outTopGenes[, "sampleStandardDeviation"], outgoingDiff[, "sampleStandardDeviation"])
wilcox.test(outTopGenes[, "sampleStandardDeviation"], outgoingDiff[, "sampleStandardDeviation"])
t.test(outTopGenes[, "expSum"], outgoingDiff[, "expSum"])
wilcox.test(outTopGenes[, "expSum"], outgoingDiff[, "expSum"])

x <- setdiff(incomingAll[,"Name"], inTopGenes[,"Name"])
incomingDiff <- incomingAll[incomingAll[,"Name"] %in% x, ]

#incoming
t.test(inTopGenes[, "sampleStandardDeviation"], incomingDiff[, "sampleStandardDeviation"])
wilcox.test(inTopGenes[, "sampleStandardDeviation"], incomingDiff[, "sampleStandardDeviation"])
t.test(inTopGenes[, "expSum"], incomingDiff[, "expSum"])
wilcox.test(inTopGenes[, "expSum"], incomingDiff[, "expSum"])


#shuffled outgoing
t.test(shuffledOutTopGenes[, "sampleStandardDeviation"], outgoingAll[, "sampleStandardDeviation"])
wilcox.test(shuffledOutTopGenes[, "sampleStandardDeviation"], outgoingAll[, "sampleStandardDeviation"])
t.test(shuffledOutTopGenes[, "expSum"], outgoingAll[, "expSum"])
wilcox.test(shuffledOutTopGenes[, "expSum"], outgoingAll[, "expSum"])


