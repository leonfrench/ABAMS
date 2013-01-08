#look in the ABACellEnriched/patterns folder for the input
aAll <- read.table("flattened.ap..43exp.txt")
bAll<- read.table("flattened.ap..68exp.txt")

aAll <- read.table("flattened.43exp.txt")
bAll<- read.table("flattened.68exp.txt")

rownames(aAll)
xA <- as.numeric(as.list(aAll["x",]))
degreeA <- as.numeric(as.list(aAll["degree",]))
expA <- as.numeric(as.list(aAll["exp",]))
degreeConA <- as.numeric(as.list(aAll["degreeCon",]))
expConA <- as.numeric(as.list(aAll["expCon",]))

x <- as.numeric(as.list(bAll["x",]))
degreeB <- as.numeric(as.list(bAll["degree",]))
expB <- as.numeric(as.list(bAll["exp",]))
degreeConB <- as.numeric(as.list(bAll["degreeCon",]))
expConB <- as.numeric(as.list(bAll["expCon",]))

paraA <- aAll["exp","Parafascicular.nucleus"]/43
paraA
paraB <- bAll["exp","Parafascicular.nucleus"]/68
paraB 
paraA/paraB

complexA <- aAll["exp","Ventral.posterior.complex.of.the.thalamus"]/43
complexA
complexB <- bAll["exp","Ventral.posterior.complex.of.the.thalamus"]/68
complexB
complexA/complexB

aAvg <- mean(expA)/43
bAvg <- mean(expB)/68

# is para above average for pattern A?
paraA > aAvg
# is complex above average for pattern A?
complexA > aAvg

# is para above average for pattern B?
paraB > bAvg
# is complex above average for pattern B?
complexB > bAvg



bAll["degree","Parafascicular.nucleus"]
bAll["degree","Ventral.posterior.complex.of.the.thalamus"]

cor(expB, expConB)
cor(degreeB, degreeConB)


cor(expA, expConA)
cor(degreeA, degreeConA)


#regress expression using space
fit <- lm(expA ~ x)
#get the residuals
expRes <- residuals(fit)
plot(degreeA, expRes)
cor(degreeA, expRes)
abline(lm(expRes ~ degreeA))

#regress expression using space
fit <- lm(expB ~ x)
#get the residuals
expRes <- residuals(fit)
plot(degreeA, expRes)
cor(degreeA, expRes)

cor(degreeConB, expRes)

abline(lm(expRes ~ degreeA))

plot(expB, x)
abline(lm(x ~ expB))

plot(expA, x)
abline(lm(x ~ expA))


plot(degreeA, x)
abline(lm(x ~ degreeA))

plot(degreeB, x)
abline(lm(x ~ degreeB))


