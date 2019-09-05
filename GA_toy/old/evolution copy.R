
# for adding minor tick marks:
library(Hmisc)
library(plotrix)
library(emdbook)
library(magicaxis)
###########################################################################
# USER INPUT

# axes for plot 1 and 2
xlim1 <- c(1, 10000)
ylim1 <- c(1e-1, 1e5)

# axes for plot 2 - top panel
xlim2 <- c(1.0,5.0e2)
ylim2 <- c(3.0e-1,3.0e2)

# axes for plot 2 - bottom panel
xlim3 <- c(1.0e2,1.5e3)
ylim3 <- c(1.0,1.0e5)

# print errors yes/no: flag=true/false
flag <- TRUE
#flag <- FALSE

###########################################################################
# read data from files
mydat <- read.table("main.dat", header=FALSE)
simu <- read.table("main.out", header=FALSE, skip=3, nrows=4, na.strings = "NA")


###########################################################################
# PLOT: parameter values vs generations
pdf("evolution1.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xlim1, ylim=ylim1, 
       axes=FALSE, main="", xlab = "", ylab = "", log="y", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
title(xlab="Generation #", line=3, cex.lab=2.0)
title(ylab="Parameter value", line=3, cex.lab=2.0)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
#text(1e-7, 2, labels=expression(paste("xxx")), cex=2.0)

yy1 <- 1e4 * mydat[,3]/1e6
yy2 <- 1e4 * mydat[,4]/1e6
yy3 <- 1e4 * mydat[,5]/1e6
yy4 <- 1e5 * mydat[,6]/1e6

lines(mydat[,1], yy1, col="red", lw=1)
lines(mydat[,1], yy2, col="blue", lw=1)
lines(mydat[,1], yy3, col="chartreuse3", lw=1)
lines(mydat[,1], yy4, col="black", lw=1)

legend("bottomleft", inset=.01, legend=c(
        expression(paste(P1)), 
        expression(paste(P2)), 
        expression(paste(P3)), 
        expression(paste(P4))), 
        bty = "n",
        pch=c("-", "-", "-", "-"), horiz=TRUE,
        col=c("red", "blue", "chartreuse3", "black"), 
        pt.cex=c(1.5, 1.5, 1.5, 1.5),
        cex=c(1.2, 1.2, 1.2, 1.2))

dev.off()

###########################################################################
# PLOT: fitness vs generations
pdf("evolution2.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xlim1, ylim=ylim1, 
       axes=FALSE, main="", xlab = "", ylab = "", log="y", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
title(xlab="Generation #", line=3, cex.lab=2.0)
title(ylab="Best Fitness", line=3, cex.lab=2.0)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
#text(1e-7, 2, labels=expression(paste("xxx")), cex=2.0)

fitness <- mydat[,2]

lines(mydat[,1], fitness, col="purple", lw=1)

dev.off()

###########################################################################
# PLOT: parameter values 
pdf(file="evolution3.pdf", width=8, height=5, onefile=F)
par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

############# TOP
plot( 1, type="n", lwd=2, xlim=xlim2, ylim=ylim2, 
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
for(i in 1:length(mydat[,2])){
if( (mydat[i,2] > 0.5 * max(mydat[,2])) && flag ){
     points(yy1[i], yy2[i], col="blue", pch=1, cex=1)
     }
}

points(simu[1,], simu[2,], col="red", pch=16, cex=2)

abline(h=12)
abline(v=34)

############# BOTTOM
plot( 1, type="n", lwd=2, xlim=xlim3, ylim=ylim3,
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
for(i in 1:length(mydat[,2])){
if( (mydat[i,2] > 0.5 * max(mydat[,2])) && flag ){
     points(yy3[i], yy4[i], col="blue", pch=1, cex=1)
     }
}

points(simu[3,], simu[4,], col="red", pch=16, cex=2)

abline(h=78)
abline(v=456)

dev.off()

###########################################################################
# PLOT: parameter uncertainties 

meanyy1 <- mean(yy1)
meanyy2 <- mean(yy2)
meanyy3 <- mean(yy3)
meanyy4 <- mean(yy4)


pdf(file="evolution4.pdf", width=8, height=5, onefile=F)
par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

##
plot(density(yy1), main="", xlab="", ylab="", xlim=c(30,38),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=2.5, cex.lab=1.5)
title(xlab="Parameter 1", line=1.6, cex.lab=1.5)
polygon(density(yy1), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=meanyy1)
##
plot(density(yy2), main="", xlab="", ylab="", xlim=c(0,150),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=2.5, cex.lab=1.5)
title(xlab="Parameter 2", line=1.6, cex.lab=1.5)
polygon(density(yy2), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=meanyy2)

##
plot(density(yy3), main="", xlab="", ylab="", xlim=c(350,550),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
#title(ylab="Probability density", line=1.8, cex.lab=1.5)
title(xlab="Parameter 3", line=1.6, cex.lab=1.5)
polygon(density(yy3), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=meanyy3)

##
plot(density(yy4), main="", xlab="", ylab="", xlim=c(0,2e4),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
#title(ylab="Probability density", line=1.8, cex.lab=1.5)
title(xlab="Parameter 4", line=1.6, cex.lab=1.5)
polygon(density(yy4), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=meanyy4)

dev.off()


######################################################################
# CORRELATION PLOT
######################################################################
pdf("evolution5.pdf", width=6, height=6, onefile=F)

samplesmat <- cbind(yy1, yy2, yy3, yy4)
pairs(~yy1+yy2+yy3+yy4, col=adjustcolor("red", alpha=0.05),  
   data=samplesmat[sample(nrow(samplesmat), size=4000, replace=FALSE),], 
   main="Simple Scatterplot Matrix")


dev.off()
