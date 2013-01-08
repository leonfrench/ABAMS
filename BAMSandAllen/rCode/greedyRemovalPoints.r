incoming <- read.table("LOOGenesInOrder.in.txt.forPlotting");
outgoing <- read.table("LOOGenesInOrder.out.txt.forPlotting");
bi <- read.table("LOOGenesInOrder.bi.txt.forPlotting");
plot(incoming)
points(outgoing, col="red")
points(bi, col="blue")

plot(incoming, type="l")
lines(outgoing, col="red")
lines(bi, col="blue")


