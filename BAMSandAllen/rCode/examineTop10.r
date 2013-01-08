
x<- read.table("ConnectivityAndAllenPartialExpressionMatrixPair.NewEnergies.txt")
x<-as.matrix(x)

means <- mean(as.data.frame(x), na.rm=TRUE)

oneGene <- matrix(0, 10, dim(x)[2])
oneGene[1,] <- x["Lhfp[73769323]", ]
oneGene[2,] <- x["Pgrmc1[797]", ]
oneGene[3,] <- x["Esr1[80342169]", ]
oneGene[4,] <- x["Zfyve26[71212229]", ]
oneGene[5,] <- x["Ephb1[79677367]", ]
oneGene[6,] <- x["Nrp2[80514091]", ]
oneGene[7,] <- x["Slc37a2[72001996]", ]
oneGene[8,] <- x["Gm597[70295620]", ]
oneGene[9,] <- x["Alpk3[71574473]", ]
oneGene[10,] <- x["Mm.371780[71210025]", ]

plot(oneGene[3,])
lines(means)

boxplot(as.data.frame(x))
points(oneGene[9,], col="red")
points(log(oneGene[9,]), col="blue")


length(which(oneGene[8,] > means))

for(i in 1:10) {
   print(length(which(is.na(oneGene[i,]))))
}

length(which(is.na(x))) / length(x)
