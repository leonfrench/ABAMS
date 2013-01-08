#input provided by RegressionVector.java


#all <- read.table("Expression.Distance.Connectivity.forR.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.forR.txt")

all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pgrmc1.Connections.forR.out.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pgrmc1.Connections.forR.in.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pgrmc1.Connections.forR.in.add.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pgrmc1.Connections.forR.out.add.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pgrmc1.Connections.forR.bi.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pcp2.Pgrmc1.Connections.forR.txt")

all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pcp2.Connections.any.forR.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pcp2.Pgrmc1.Rgs3.Connections.forR.txt")
all <- read.table("Expression.Distance.Connectivity.Nomenclature.Pcp2.Pgrmc1.Rgs3.Connections.TopExpression.forR.txt")

#all <- read.table("Expression.Distance.Connectivity.Nomenclature.forR.Prop.txt")
#all <- read.table("Expression.Distance.Connectivity.Nomenclature.Dimensions.Connections.forR.txt")
#all <- read.table("Target.Explain.Residuals.forR.txt")
dim(all)

rownames(all)
exp <- as.numeric(as.list(all["Expression",]))
space <- as.numeric(as.list(all["Distance",]))
conCor <- as.numeric(as.list(all["Connectivity",]))
nomen <- as.numeric(as.list(all["Nomenclature",]))
pgrmc1 <- as.numeric(as.list(all["Pgrmc1",]))
connections <- as.numeric(as.list(all["Connections",]))
pcp2 <- as.numeric(as.list(all["Pcp2",]))
rgs3 <- as.numeric(as.list(all["Rgs3",]))
topExp <- as.numeric(as.list(all["TopExpression",]))

#dim <- as.numeric(as.list(all[5,]))
#connections <- as.numeric(as.list(all[6,]))

range(exp)
range(conCor)
range(space)
range(nomen)
range(dim)
range(connections)



plot(space, exp)
plot(expression, conCor)

#regress expression using space
fit <- lm(exp ~ space)
#get the residuals
expRes <- residuals(fit)
plot(space, expRes)

#test if it fits anymore
fitAfter <- lm(expRes ~ space)
summary(fitAfter)

#before and after correlation
cor(conCor, exp)
cor(conCor, expRes)

#before and after nomenclature correlation
cor(exp, nomen)
cor(expRes, nomen)
fitNomen <- lm(expRes ~ nomen)
summary(fitNomen)

fitCon <- lm(expRes ~ conCor)
summary(fitCon)


plot(density(connected), col="red")
lines(density(disconnected), col="blue")


#regress connectivty using space
fit <- lm(conCor ~ space)
#get the residuals
conCorRes <- residuals(fit)
plot(space, conCorRes)

#before and after correlation
cor(conCor, exp)
cor(conCorRes, exp)
cor(conCorRes, expRes)

cor(conCor, nomen)
cor(conCorRes, nomen)
fitNomen <- lm(conCorRes ~ nomen)
summary(fitNomen)


#regress nomen using space
fit <- lm(nomen ~ space)
#get the residuals
nomenRes <- residuals(fit)
plot(space, nomenRes)
#before and after correlation
cor(space, nomen)
cor(space, nomenRes)

cor(conCor, nomen)
cor(conCorRes, nomenRes)
cor(nomen, exp)
cor(nomenRes, expRes)


#old hypothesis
length(connections[connections == 1])
length(connections[connections == 0])
connectedInd <- which(connections == 1, arr.ind=TRUE)
disconnectedInd <- which(connections == 0, arr.ind=TRUE)
connected = expRes[connectedInd]
disconnected = expRes[disconnectedInd]
ks.test(connected, disconnected)
t.test(connected, y=disconnected)


#pgrmc1
length(pgrmc1[connections ==1])
length(pgrmc1[connections ==0])

boxplot(pgrmc1[connections ==1], pgrmc1[connections==0])

mean(pgrmc1[connections ==0], na.rm=TRUE)
mean(pgrmc1[connections ==1], na.rm=TRUE)

t.test(pgrmc1[connections ==1], pgrmc1[connections==0])

t.test(pcp2[connections ==1], pcp2[connections==0])
png("pcp2.any.add.png")
plot(density(pcp2[connections ==1], na.rm=TRUE, kernel = "rectangular"), col="red", main="Pcp2 expression levels versus connectivity", xlab="Sum of Pcp2 Expression", ylab="Region Pair Density")
lines(density(pcp2[connections ==0], na.rm=TRUE, kernel = "rectangular"), col="black")
legend("topright", c("Not connected", "Connected"), fill=c("black", "red"))
dev.off()


postscript("pgrmc1.any.add.ps")
png("pgrmc1.any.add.png")
plot(density(pgrmc1[connections ==1], na.rm=TRUE, kernel = "rectangular"), lty="solid", main="Pgrmc1 expression levels versus connectivity", xlab="Sum of Pgrmc1 Expression", ylab="Region Pair Density")
lines(density(pgrmc1[connections ==0], na.rm=TRUE, kernel = "rectangular"), lty="dashed")
legend("topleft", c("Connected", "Not connected"), lty=c("solid", "dashed"))
dev.off()

plot(density(rgs3[connections ==0], na.rm=TRUE, kernel = "rectangular"), col="black", main="Rgs3 expression levels versus connectivity", xlab="Sum of Pgrmc1 Expression", ylab="Region Pair Density")
lines(density(rgs3[connections ==1], na.rm=TRUE, kernel = "rectangular"), col="red")
legend("topleft", c("Not connected", "Connected"), fill=c("black", "red"))

plot(density(rgs3[connections ==0], na.rm=TRUE), col="black", main="Rgs3 expression levels versus connectivity", xlab="Sum of Pgrmc1 Expression", ylab="Region Pair Density")
lines(density(rgs3[connections ==1], na.rm=TRUE), col="red")
legend("topleft", c("Not connected", "Connected"), fill=c("black", "red"))


png("pgrmc1.any.add.hist.png")
hist(pgrmc1[connections ==0], main="Pgrmc1 expression levels versus connectivity", xlab="Sum of Pgrmc1 Expression", ylab="Region Pair Frequency")
hist(pgrmc1[connections ==1], add=TRUE, border="red")
legend("topleft", c("Not connected", "Connected"), fill=c("black", "red"))
dev.off()

plot(density(conCor[connections ==1], na.rm=TRUE, kernel = "rectangular"), col="red", main="Sum of conCor expression levels for connected pairs")
lines(density(conCor[connections ==0], na.rm=TRUE, kernel = "rectangular"), col="black")

hist(conCor[connections ==0])
hist(conCor[connections ==1], add=TRUE, border="red")




require(hdrcde)

x <- topExp
y <- conCor
par(mfrow=c(1,2))
plot(x,y, pch="+", cex=.5)
hdr.boxplot.2d(x,y)
hdr.boxplot.2d(x,y, prob=c(0.01,0.05,0.1,0.50, 0.9))

