energyTable <- read.table("tempEnergyMatrix.txt")

#sort
energyTable <- energyTable[order(row.names(energyTable)),]
energyTable <- energyTable[,order(colnames(energyTable))]

#convert to matrix
energyMatrix <- as.matrix(energyTable)

#remove two brain regions
#Posterodorsal.preoptic.nucleus is all zeroes
#Nucleus.y is all zeroes
energyMatrix <- subset(energyMatrix, select=c(-Posterodorsal.preoptic.nucleus,-Nucleus.y))

z <- energyMatrix
connected <- scan(file = "SN Efferent Regions in ABA.txt", what = 'character')
genes <- scan(file = "dopamine receptor ABA.txt", what = 'character')
controls <- scan(file = "control genes.txt", what = 'character')
#genes <- controls

#pull out the stuff we want
x <- z[genes, connected]
notconnected <- setdiff(colnames(z) , connected)
y <- z[genes,notconnected]
#par(mar=c(37,24,24,22))
heatmap(x, scale="no", Rowv=NA, Colv=NA)
heatmap(log(x), scale="no", Rowv=NA, Colv=NA)
wilcox.test(x,y)
heatmap(cbind(x,y), scale="no", Rowv=NA, Colv=NA)
#drd1a receptors:
wilcox.test(x[1,],y[1,])$p.value
wilcox.test(x[2,],y[2,])$p.value
#drd 2 receptors:
wilcox.test(x[3,],y[3,])$p.value
wilcox.test(x[4,],y[4,])$p.value
#drd3
wilcox.test(x[5,],y[5,])$p.value
#drd4
wilcox.test(x[6,],y[6,])$p.value
#drd5
wilcox.test(x[7,],y[7,])$p.value

#drd 2 plots
boxplot(x[3,], y[3,])
boxplot(log(x[3,]), log(y[3,]))

#test using emperic null distribution

#all regions (first 24 are connected)
both <- cbind(x,y)

#just drd2 image series
both <- both[3,]

mean(both[1:24])
#[1] 3.15
mean(both[25:207])

samples <- matrix()
#resample brain regions (columns)
#compute average of first 24
for(i in 1:1000) {     
  both.sample <- both[sample(1:length(both))]
  samples[i] <- mean(both.sample[1:24])
}

mean(both[1:24])
samples[samples>mean(both[1:24])]
length(samples[samples>mean(both[1:24])])

