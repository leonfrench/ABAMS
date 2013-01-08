aList <- read.table("LOOGenesInOrder.outgoing.partialcon.A.txt")
bList <- read.table("LOOGenesInOrder.outgoing.partialcon.962.txt")


a <- rank(aList)
b <- rank(bList)
cor(a,b, method="spearman")

#aList <- a
#bList <- b


listSize <- 10

size <- dim(aList)[1]
#size <- length(aList)
aTop <- aList[(size-listSize):size,]
bTop <- bList[(size-listSize):size,]
inter <- intersect(aTop, bTop)
interSize <- length(inter)
interSize/listSize

