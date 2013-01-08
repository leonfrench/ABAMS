#run Preliminary R before this

cor <- cor(energyMatrix)
#reorder based on dendrogram
cor <- cor[map$rowInd, map$colInd]
parentColsOrdered <- parentCols[map$colInd]
heatmap(cor, Rowv=NA, Colv=NA, scale="no")

heatmap(cor[1:27,1:27], Rowv=NA, Colv=NA, scale="no", margins=c(15,15))
png("heatMapReducedbyRegions 1to26.png", width=740, height=640) 
heatmap(cor[1:27,1:27], Rowv=NA, Colv=NA, scale="no", margins=c(15,15))
dev.off()


heatmap(cor[181:207,181:207], Rowv=NA, Colv=NA, scale="no", margins=c(15,15))
png("heatMapReducedbyRegions 181to207.png", width=740, height=640) 
heatmap(cor[181:207,181:207], Rowv=NA, Colv=NA, scale="no", margins=c(15,15))
dev.off()

start <- 187
end <- 207
heatmap(cor[start:end,start:end], Rowv=NA, Colv=NA, scale="no", margins=c(15,15), ColSideColors=parentColsOrdered[start:end])


