decreasing <- read.table("LOOGenesInOrder.incoming.decreasing.txt")
increasing <- read.table("LOOGenesInOrder.incoming.txt")

a <- c("a", "b", "c", "d", "e")
b <- c("e", "c", "b", "d", "a")
a <- rank(a)
b <- rank(b)
cor(a,b, method="spearman")

a <- rank(decreasing)
b <- rev(rank(increasing))
cor(a,b, method="spearman")

cor(decreasing, rev(increasing), method="spearman")
cor(rev(decreasing), increasing, method="spearman")


listSize <- 550

size <- dim(decreasing)[1]
decTop <- decreasing[(size-listSize):size,]
incTop <- increasing[(size-listSize):size,]
intersect(decTop, incTop)

decBot <- decreasing[1:listSize,]
incBot <- increasing[1:listSize,]
intersect(decBot, incBot)

intersect(decBot, incTop)
intersect(decTop, incBot)

