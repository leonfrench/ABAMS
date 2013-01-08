x<-read.table("ConnectivityAndAllenExpressionMatrixPair.forPlotting.txt")
x<-read.table("AllenMatrixPair.Distance.NomenclatureDistance.forPlotting.txt")
x<-read.table("AllenMatrixPair.Dimensions.NomenclatureDistance.forPlotting.txt")
x<-read.table("AllenMatrixPair.Distance.NomenclatureDistance.forPlotting.txt")
x<-read.table("AllenMatrixPair.ExpressionEnergy.Distance.forPlotting.txt")
x<-read.table("ConnectivityAndAllenSpacePair.Connectivity.EuclidianDistance.forPlotting.txt")
#x<-read.table("ConnectivityAndAllenSpacePair.forPlotting.txt")
colnames(x)
dim(x)
cor(x)
cor(x, use = "complete.obs")
cor(abs(x), use = "complete.obs")
cor(x, method = "spearman")
range(x[,"Connectivity"])
range(x[,"NewEnergies"])

plot(x)

png("ConnectivityAndAllenExpressionMatrixPair.forPlotting.3.png", pointsize=18, width=800, height=800, res=380) 
N <- dim(x)[1]
plot(x, xlab="Connectivity Correlation", ylab="Expression Correlation", type="n")
abline(lm(x[,"NewEnergies"] ~ x[,"Connectivity"]), col="red")
symbols(x, circles = rep(1,N), inches = 0.016, bg=rep(1,N), add=TRUE )
dev.off()



png("ConnectivityAndAllenExpressionMatrixPair.forPlotting.1.png", width=800, height=800) 
plot(x)
abline(lm(x[,"NewEnergies"] ~ x[,"Connectivity"]))
dev.off()

#"Connectivity" "NewEnergies"

oneConn <-  row.names(x)[which(x[,"Connectivity"]==1)]
oneExp <-row.names(x)[which(x[,"NewEnergies"]==1)]

intersect(oneConn, oneExp)

hist(x[,"Connectivity"], breaks=500)

hist(x[,"NewEnergies"], breaks=500)


#regression
y <- x[,"NewEnergies"]
x <- x[,"Connectivity"]
plot(x,y)         # make a plot 
abline(lm(y ~ x)) # plot the regression line 
lm(y ~ x)         # the basic values of the regression analysis 
lm.result = lm(y ~ x) # prepare for graphing regression 
summary(lm.result) # look at the summary 
plot(lm.result) # plot the regression (will produce four graphs) 
par(mfrow=c(2,2));plot(lm.result) # plots all four graphs in a single graph. 



#show the fancy graph

def.par <- par(no.readonly = TRUE) # save default, for resetting...

y <- x[,"NewEnergies"]
x <- x[,"Connectivity"]

png("ConnectivityAndAllenExpressionMatrixPair.forPlotting.2.png", width=800, height=800) 
xhist <- hist(x,  breaks=100, plot=FALSE)
yhist <- hist(y,  breaks=100, plot=FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- c(min(x),max(x))
yrange <- c(min(y),max(y))
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
#layout.show(nf)

par(mar=c(3,3,1,1))
plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
abline(lm(y ~ x)) # plot the regression line 
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)

par(def.par)

dev.off()

