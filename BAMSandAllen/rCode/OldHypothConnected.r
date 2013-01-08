con <- read.csv("connect.txt", header=FALSE, sep=",")
dis <- read.csv("disconnect.txt", header=FALSE, sep=",")
con<-as.matrix(con)[1,]
dis<-as.matrix(dis)[1,]
t.test(con, y=dis)
ks.test(con, dis)

png("Figure S1z.png", width=840, height=840) 
plot(density(dis), col="red", main="Density Plot of Expression Correlation Between Region Pairs", xlab="Expression correlation")
lines(density(con), col="blue")
legend("topleft", c("Disconnected", "Connected"), fill=c("red", "blue"))
dev.off()


png("Figure S1.png")
postscript("Figure S1.ps")
plot(density(dis, kernel = "rectangular"), lty="dashed", main="Expression Correlation versus connectivity", xlab="Expression correlation", ylab="Region Pair Density")
lines(density(con, kernel = "rectangular"), lty="solid")
legend("topleft", c("Disconnected", "Connected"), lty=c("dashed", "solid"))
dev.off()

superhist2pdf <- function(x, filename = "super_histograms.pdf",
dev = "pdf", title = "Superimposed Histograms", nbreaks ="Sturges") {
junk = NULL
grouping = NULL
for(i in 1:length(x)) {
junk = c(junk,x[[i]])
grouping <- c(grouping, rep(i,length(x[[i]]))) }
grouping <- factor(grouping)
n.gr <- length(table(grouping))
xr <- range(junk)
histL <- tapply(junk, grouping, hist, breaks=nbreaks, plot = FALSE)
maxC <- max(sapply(lapply(histL, "[[", "counts"), max))
if(dev == "pdf") { pdf(filename, version = "1.4") } else{}
if((TC <- transparent.cols <- .Device %in% c("pdf", "png"))) {
cols <- hcl(h = seq(30, by=360 / n.gr, length = n.gr), l = 65, alpha = 0.5) }
else {
h.den <- c(10, 15, 20)
h.ang <- c(45, 15, -30) }
if(TC) {
plot(histL[[1]], xlim = xr, ylim= c(0, maxC), col = cols[1], xlab = "x", main = title) }
else { plot(histL[[1]], xlim = xr, ylim= c(0, maxC), density = h.den[1], angle = h.ang[1], xlab = "x") }
if(!transparent.cols) {
for(j in 2:n.gr) plot(histL[[j]], add = TRUE, density = h.den[j], angle = h.ang[j]) } else {
for(j in 2:n.gr) plot(histL[[j]], add = TRUE, col = cols[j]) }
invisible()
if( dev == "pdf") {
dev.off() }
}

# How to use the function:
l1 = list(dis,con)
superhist2pdf(l1, nbreaks=50)
