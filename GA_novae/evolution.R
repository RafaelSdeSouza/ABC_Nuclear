###################################################################
# evolution.R
#
# - display output results from GA
###################################################################
#
# for adding minor tick marks:
library(Hmisc)
library(plotrix)
library(emdbook)
library(magicaxis)
library(data.table)

###################################################################
# read data input
###################################################################
# read parameter evolution
mydat <- read.table("main.dat", header=FALSE)
# read fittest parameters
simu <- read.table("main.out", header=FALSE, skip=3, nrows=7, 
           na.strings = "NA")
# read parameter bounds
bounds <- read.table("main.in", header=FALSE,  sep=",", skip=16, nrows=7,
             colClasses = c(rep("numeric", 2)), comment.char="!")
# read number of generations and digits
# ngen[1,1]: number of generations
# ngen[2,1]: number of digits
ngen <- read.table("main.in", header=FALSE,  sep="", skip=3, nrows=2,
             colClasses = c(rep("numeric", 1)), comment.char="!")

bnd_l <- vector()
bnd_h <- vector()

# T in MK
bnd_l[1] <- bounds[1,1]/1e6
bnd_h[1] <- bounds[1,2]/1e6

bnd_l[2] <- bounds[2,1]
bnd_h[2] <- bounds[2,2]

bnd_l[3] <- bounds[3,1]
bnd_h[3] <- bounds[3,2]

bnd_l[4] <- bounds[4,1]
bnd_h[4] <- bounds[4,2]

bnd_l[5] <- bounds[5,1]
bnd_h[5] <- bounds[5,2]

bnd_l[6] <- bounds[6,1]
bnd_h[6] <- bounds[6,2]

bnd_l[7] <- bounds[7,1]
bnd_h[7] <- bounds[7,2]

###################################################################
### USER INPUT

# parameter values vs generations
xlim1 = c(0, ngen[1,1])
ylim1 = c(1, 1e4)



######################################################################
# PLOT 1: parameters
######################################################################
## set margins in order south, west, north, east
## oma is "outer margin" of entire figure
## mar is the margin of individual panels, which sets margin sizes in
##    order bottom, left, top, right
## mfcol=c(nrows, ncols) fills in the matrix by columns 
## tck sets tick mark lengths; negative value makes them point
##    outward
## las=1 shows tick mark labels in horizontal orientation
## mgp sets axis label locations relative to edge of inner plot window;
##   first value represents label locations [xlab, ylab], the second
##   the tick mark labels, the third the tick marks; default is c(3,1,0)
## cex controls symbol size
## cex.yy controls label sizes
## pch is the symbol id
## xlim describes plot limits
## main adds a title
pdf(file="evolution1.pdf", width=8, height=5, onefile=F)
par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

################ PLOT lt 
#plot( 1, type="n", lwd=2, xlim=c(200, 300), ylim=c(500, 1000), 
plot( 1, type="n", lwd=2, xlim=c(bnd_l[1], bnd_h[1]), ylim=c(bnd_l[2], bnd_h[2]), 
       main="", xlab = "", ylab = "", axes=FALSE,
       cex=1.0, cex.lab=1.0, cex.axis=1.2, cex.main=1.0, log="xy" )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.2,0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste(T (MK))), line=2.0, cex.lab=1.5)  
title(ylab=expression(paste(rho, " ",(g/cm^3))), line=2.5, cex.lab=1.5)

# add simulations
points(simu[1,]/1e6, simu[2,], col="red", pch=16, cex=2)

abline(v=250)
abline(h=700)

################ PLOT lb 
#plot( 1, type="n", lwd=2, xlim=c(100, 500), ylim=c(10, 200),
plot( 1, type="n", lwd=2, xlim=c(bnd_l[3], bnd_h[3]), ylim=c(bnd_l[4], bnd_h[4]),
       main="", xlab = "", ylab = "", 
       cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.0, log="xy", axes=FALSE )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0,0.4,0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste(tau[T], " ", (s))), line=1.9, cex.lab=1.5)
title(ylab=expression(paste(tau[rho], " ", (s))), line=2.5, cex.lab=1.5)

# add simulations
points(simu[3,], simu[4,], col="red", pch=16, cex=2)

abline(v=200)
abline(h=50)

################ PLOT rt 
## ratio of WD 12C/16O versus mixing fraction f
#plot( 1, type="n", lwd=2, xlim=c(1, 10), ylim=c(0.2, 0.7),
plot( 1, type="n", lwd=2, xlim=c(bnd_l[7], bnd_h[7]), ylim=c(bnd_l[5], bnd_h[5]),
       main="", xlab = "", ylab = "", 
       cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.0, log="x", axes=FALSE )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0.0,0.3,0.0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste(Mixing," ",f)), line=2.0, cex.lab=1.5)
title(ylab=expression(paste("X"[WD],"(",NULL^"12","C",")")), 
        line=1.9, cex.lab=1.5)

wdc12 <- simu[5,]
points(simu[7,], wdc12, col="red", pch=16, cex=2)

abline(v=3)
abline(h=0.5)

################ PLOT rb 
plot( 1, type="n", lwd=2, xlim=c(bnd_l[5], bnd_h[5]), ylim=c(bnd_l[6], bnd_h[6]),
       main="", xlab = "", ylab = "", 
       cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.0, log="", axes=FALSE )
# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number
magaxis(mgp=c(0.0,0.3,0.0), minorn=0, cex=1.0, cex.lab=1.0, cex.main=1.0, cex.axis=1.2)
box()

# control distance between axis and label [line=...]
title(xlab=expression(paste("X"[WD],"(",NULL^"12","C",")")), 
        line=1.9, cex.lab=1.5)
title(ylab=expression(paste("X"[WD],"(",NULL^"22","Ne",")")), 
        line=2.8, cex.lab=1.5)
        
wdc12 <- simu[5,]
wdne22 <- simu[6,]

points(wdc12, wdne22, col="red", pch=16, cex=2)

#abline(v=3)
#abline(h=0)

dev.off()

######################################################################
# PLOT 2: evolution
######################################################################
pdf(file="evolution2.pdf", width=5, height=5, onefile=F)
par(mfcol=c(2,1), mar=c(2.5,5,0,1), oma=c(1.0,2.0,1.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

################ PLOT top 
plot( 1, type="n", lwd=2 , col="black" , xlim=xlim1, ylim=ylim1, 
       axes=FALSE, main="", xlab = "", ylab = "", log="y", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
#title(xlab="Generation #", line=2.0, cex.lab=1.5)
title(ylab="Parameter value", line=3, cex.lab=1.5)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.2,0), cex.axis=1.1)
box()

# plot lines
lines(mydat[,1], (mydat[,3]/(10^ngen[2,1]))*(bnd_h[1] - bnd_l[1]) + bnd_l[1], col="red", lw=1)
lines(mydat[,1], (mydat[,4]/(10^ngen[2,1]))*(bnd_h[2] - bnd_l[2]) + bnd_l[2], col="blue", lw=1)
lines(mydat[,1], (mydat[,5]/(10^ngen[2,1]))*(bnd_h[3] - bnd_l[3]) + bnd_l[3], col="chartreuse3", lw=1)
lines(mydat[,1], (mydat[,6]/(10^ngen[2,1]))*(bnd_h[4] - bnd_l[4]) + bnd_l[4], col="black", lw=1)

legend("bottomright", inset=.01, legend=c(
        expression(paste(T6)), 
        expression(paste(rho)), 
        expression(paste(tau[T])), 
        expression(paste(tau[rho]))), 
        bty = "n",
        pch=c("-", "-", "-", "-"), horiz=TRUE,
        col=c("red", "blue", "chartreuse3", "black"), 
        pt.cex=c(1.5, 1.5, 1.5, 1.5),
        cex=c(1.2, 1.2, 1.2, 1.2))

################ PLOT bottom
# fitness called "loss" in loss.f
fitness <- mydat[,2]
maximum <- max(fitness)

ylim2 <- c(0, 1.5 * maximum)

plot( 1, type="n", lwd=2 , col="black" , xlim=xlim1, ylim=ylim2, 
       axes=FALSE, main="", xlab = "", ylab = "", log="", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
title(xlab="Generation #", line=1.5, cex.lab=1.5)
title(ylab="Best fitness", line=3, cex.lab=1.5)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.2,0), cex.axis=1.1)
box()

lines(mydat[,1], fitness, col="purple", lw=1)
legend("topleft", legend=format(maximum, nsmall=4))

dev.off()





