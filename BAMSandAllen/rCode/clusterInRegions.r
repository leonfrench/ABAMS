#run Preliminary R before this

cor <- cor(energyMatrix)
#reorder based on dendrogram
cor <- cor[map$rowInd, map$colInd]

seenRegions <- intersect(rownames(regionMatrix), colnames(energyMatrix))
seenRegions2 <- regionMatrix[seenRegions,]
seenRegions2 <- as.matrix(seenRegions2)
#mapping from region to parent, secondary order is from the dendrogram
seenRegions2[rownames(cor),1]
order(seenRegions2[rownames(cor),1])
sorted <- (cor[order(seenRegions2[rownames(cor),1]),order(seenRegions2[colnames(cor),1])])

#Colours, do the dendrogram ordering
parentColsOrdered <- parentCols[map$colInd]
#do the region sort ording
parentColsOrdered <- parentColsOrdered[order(seenRegions2[rownames(cor),1])]

heatmap(sorted, Rowv=NA, Colv=NA, scale="no", ColSideColors=parentColsOrdered, RowSideColors=parentColsOrdered)

start <- 1
end <- 47
heatmap(sorted[start:end,start:end], Rowv=NA, Colv=NA, scale="no", margins=c(15,15), ColSideColors=parentColsOrdered[start:end])


#write heatmap to PNG
png("heatMapbyRegions.png", width=1640, height=1640) 
heatmap(sorted, Rowv=NA, Colv=NA, scale="no", ColSideColors=parentColsOrdered, RowSideColors=parentColsOrdered)
dev.off()

