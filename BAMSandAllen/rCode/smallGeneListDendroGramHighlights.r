#ran in data/Gene to Gene correlation study
require(graphics)

energy <- read.table("ConnectivityAndAllenExpressionMatrixPair.NewEnergies.out.topList.txt")
energy <- as.matrix(energy)
energy <- energy[-197,]
dim(energy)
length(energy[energy==0])

ht <- heatmap(energy, method = "average", keep.dendro = TRUE)

rowD <- ht$Rowv


coolGenes <- c("Ephb1[79677367]","Epha7[402]","L1cam[80342072]","L1cam[79591595]","Epha7[402]","Gpc3[69173280]","Gpc3[71020431]","Nrp1[79913171]","Bcl11b[74990505]","Cd24a[79591541]","Cit[75079801]","Nrtn[69117402]","Prkg1[73521817]","Apc[74881517]","Nr2e1[77464920]","Ptprz1[71487157]","Reln[890]","Ptk2[75774672]","Sema6a[68844851]")


     ## toy example to set colored leaf labels :
     local({
       colLab <<- function(n) {
           if(is.leaf(n)) {
             a <- attributes(n)
             i <<- i+1
             #check value of label
             if (a$label %in% coolGenes) {
                attr(n, "nodePar") <-
                    c(a$nodePar, list(lab.col = "red"))#, lab.cex=0.67)
             }
           }
           n
       }
       i <- 0
      })
      
#      png("dendroRows.png", width=3740, height=540, res=500) 
      postscript("dendroRows.ps", width=55, height=9, paper="special", horizontal = FALSE) 
      par(mar=c(22,4,4,2))
     dL <- dendrapply(rowD, colLab)
     plot(dL) #, nodePar=list(pch = c(NA),cex=0.67, lab.cex=0.67)) ## --> colored labels!
dev.off()

