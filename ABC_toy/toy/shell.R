##################################################################################
# shell.R
#
# 
##################################################################################
library(sfsmisc)
library(magicaxis)
library(MASS)
library(lattice)
library(plotrix)

##################################################################################

# plotting transparency
alphap <- 0.0
 
# for colors, see:
# http://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=4
col1 <- "#e31a1c"
col2 <- "#fd8d3c"
col3 <- "#fecc5c"
#col4 <- "#ffffb2"

# symbols
# character
ch1 <- 19
ch2 <- 19
ch3 <- 19
# size
sz1 <- 1.0
sz2 <- 1.0
sz3 <- 1.0

##################################################################################

roundDown <- function(x) 10^floor(log10(x))
roundUp <- function(x) 10^ceiling(log10(x))

##################################################################################
# DISPLAYED RESULTS ARE FOR ALL SOLUTIONS 
##################################################################################
## make multi-panel graph;
## set margins in order south, west, north, east
## oma is "outer margin" of entire figure
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
## mfcol=c(nrows, ncols) fills in the matrix by columns 
## tck sets tick mark lengths; negative value makes them point
##    outward
## las=1 shows tick mark labels in horizontal orientation
## mgp sets axis label locations relative to edge of inner plot window;
###   first value represents label locations [xlab, ylab], the second
###   the tick mark labels, the third the tick marks; default is c(3,1,0)

# name of file with simulations
simu <- read.table("shell.out", header=TRUE, skip=0, na.strings = "NA")

## redirect graphic output to pdf file  
pdf(file="shell.pdf", width=8, height=5, onefile=F)

par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

######## PLOT lt ##################################################################
# Tpeak versus rhopeak
plot( 1, type="n", lwd=2, xlim=c(1.0e-1,1.0e3), ylim=c(1.0e-1,1.0e3), 
       main="", xlab = "", ylab = "", axes=FALSE,
       cex=1.0, cex.lab=1.0, cex.axis=1.2, cex.main=1.0, log="xy" )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.2,0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste("par1")), line=2.0, cex.lab=1.5)  
title(ylab=expression(paste("par2")), line=1.7, cex.lab=1.5)

# add simulations
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.5 & simu[i,5] <= 2.0) {
      points(simu[i,1], simu[i,2], col=col3, pch=ch3, cex=sz3)
    }
}
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.2 & simu[i,5] <= 1.5){
      points(simu[i,1], simu[i,2], col=col2, pch=ch2, cex=sz2)
    }
}
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.0 & simu[i,5] <= 1.2) {
      points(simu[i,1], simu[i,2], col=col1, pch=ch1, cex=sz1)
    }
}
#for (i in 1:length(simu$ratio_max)) {
#    if(simu[i,9] <= 1.5) {
#      points(simu[i,1]/1e9, simu[i,2], col=col1, pch=ch1, cex=sz1)
#    }
#}

abline(h=12)
abline(v=34)

######## PLOT lb ##################################################################
## timeT versus timerho
plot( 1, type="n", lwd=2, xlim=c(1.0e-1,1.0e3), ylim=c(1.0e-1,1.0e3),
       main="", xlab = "", ylab = "", 
       cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.0, log="xy", axes=FALSE )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.4,0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste("par3")), line=1.9, cex.lab=1.5)
title(ylab=expression(paste("par4")), line=1.8, cex.lab=1.5)

# add simulations
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.5 & simu[i,5] <= 2.0) {
      points(simu[i,3], simu[i,4], col=col3, pch=ch3, cex=sz3)
    }
}
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.2 & simu[i,5] <= 1.5){
      points(simu[i,3], simu[i,4], col=col2, pch=ch2, cex=sz2)
    }
}
for (i in 1:length(simu$ratio_max)) {
    if(simu[i,5] >= 1.0 & simu[i,5] <= 1.2) {
      points(simu[i,3], simu[i,4], col=col1, pch=ch1, cex=sz1)
    }
}
#for (i in 1:length(simu$ratio_max)) {
#    if(simu[i,9] <= 1.5) {
#      points(simu[i,4], simu[i,5], col=col1, pch=ch1, cex=sz1)
#    }
#}

abline(h=78)
abline(v=456)

# plot legend
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", c("1.0-1.2 ", "1.2-1.5 ", "1.5-2.0 "), xpd = TRUE, horiz = TRUE,
     inset = c(0, 0), bty = "n", pch = c(ch1,ch2,ch3), col = c(col1, col2,
     col3), cex = 1.2)

dev.off( )


