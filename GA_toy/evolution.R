
# for adding minor tick marks:
library(Hmisc)
library(plotrix)
library(emdbook)
library(magicaxis)
###########################################################################
# read data from files
mydat <- read.table("main.dat", header=FALSE)
simu <- read.table("main.out", header=FALSE, skip=3, nrows=4, na.strings = "NA")
# read number of generations
ngen <- read.table("main.in", header=FALSE,  sep="", skip=3, nrows=2,
             colClasses = c(rep("numeric", 1)), comment.char="!")

###########################################################################
# USER INPUT

# PLOT 1
xlim1 <- c(1, 300)
ylim1 <- c(2e-1, 200)
#ylim1 <- c(4e2, 5e2)

# PLOT 2
xlim2 <- c(100, 2000)
ylim2 <- c(1, 1e5)

# PLOT 3
xlim3 <- c(0, ngen[1,1])
ylim3 <- c(3.0e-1,1.0e4)

# PLOT 4
xlim4 <- c(0, ngen[1,1])

# use exactly the same scaling of parameters as in loss.f
yy1 <- 1e4 * mydat[,3]/(10^ngen[2,1])
yy2 <- 1e4 * mydat[,4]/(10^ngen[2,1])
yy3 <- 1e4 * mydat[,5]/(10^ngen[2,1])
yy4 <- 1e5 * mydat[,6]/(10^ngen[2,1])

######################################################################
# PLOT 1: summary
######################################################################
pdf(file="evolution1.pdf", width=8, height=5, onefile=F)
par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

############# TOP LEFT
plot( 1, type="n", lwd=2, xlim=xlim1, ylim=ylim1, 
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

points(simu[1,], simu[2,], col="red", pch=16, cex=2)

abline(h=12)
abline(v=34)

############# BOTTOM LEFT
plot( 1, type="n", lwd=2, xlim=xlim2, ylim=ylim2,
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
#
points(simu[3,], simu[4,], col="red", pch=16, cex=2)

abline(h=78)
abline(v=456)

############# TOP RIGHT
plot( 1, type="n", lwd=2 , col="black" , xlim=xlim3, ylim=ylim3, 
       axes=FALSE, main="", xlab = "", ylab = "", log="y", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
#title(xlab="Generation #", line=3, cex.lab=1.5)
title(ylab="Parameter value", line=3, cex.lab=1.5)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()

# plot reaction label
#text(1e-7, 2, labels=expression(paste("xxx")), cex=2.0)

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

############# BOTTOM RIGHT
# fitness called "loss" in loss.f
fitness <- mydat[,2]
maximum <- max(fitness)

ylim4 <- c(0, 1.5 * maximum)

plot( 1, type="n", lwd=2 , col="black" , xlim=xlim4, ylim=ylim4, 
       axes=FALSE, main="", xlab = "", ylab = "", log="", yaxs='i', xaxs='i' )

# control distance between axis and label [line=...]
title(xlab="Generation #", line=2, cex.lab=1.5)
title(ylab="Best Fitness", line=3, cex.lab=1.5)

# control distance tick mark labels and axis
# don't touch first number
# second number controls distance tick mark labels and axis
# don't touch third number

magaxis(mgp=c(0,0.6,0), cex.axis=1.3)
box()


lines(mydat[,1], fitness, col="purple", lw=1)
legend("topright", legend=format(maximum, nsmall=4))

dev.off()

######################################################################
######################################################################
######################################################################



if(FALSE)
{
######################################################################
# PLOT 2: parameter values vs generations
######################################################################
pdf("evolution2.pdf",width=10,height=5,onefile=F)
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

abline(v=gen_min, lty=2)
dev.off()

######################################################################
# PLOT 3: chi^2 vs generations
######################################################################
pdf("evolution3.pdf",width=10,height=5,onefile=F)
par(mfcol=c(1,1), mar=c(4.0,7.0,1.0,6.0), oma=c(0.5,1.0,0.5,1.0), tck=0.02, 
     las=1)

# plot axes only...add lines...then data
plot( 1, type="n", lwd=2 , col="black" , xlim=xlim10, ylim=ylim10, 
       axes=FALSE, main="", xlab = "", ylab = "", log="", yaxs='i', xaxs='i' )

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

fitness <- 1/(mydat[,2]^2)

# minimum chi^2
minimum <- min(fitness)

lines(mydat[,1], fitness, col="purple", lw=1)

abline(v=gen_min, lty=2)

dev.off()

######################################################################
# PLOT 4: parameter uncertainties - NEEDS WORK
######################################################################

# analyze means etc. only beyond a user-defined generation number
newyy1 <- yy1[gen_min:length(yy1)]
newyy2 <- yy2[gen_min:length(yy2)]
newyy3 <- yy3[gen_min:length(yy3)]
newyy4 <- yy4[gen_min:length(yy4)]

meanyy1 <- mean(newyy1)
meanyy2 <- mean(newyy2)
meanyy3 <- mean(newyy3)
meanyy4 <- mean(newyy4)


pdf(file="evolution4.pdf", width=8, height=5, onefile=F)
par(mfcol=c(2,2), mar=c(3,5,0,1), oma=c(1.0,2.0,3.0,2.0), tck=0.05, 
    las=1, mgp=c(0,0.1,0))

##
plot(density(newyy1), main="", xlab="", ylab="", xlim=c(33,35),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=2.5, cex.lab=1.5)
title(xlab="Parameter 1", line=1.6, cex.lab=1.5)
polygon(density(newyy1), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=format(meanyy1, nsmall=4))
##
plot(density(newyy2), main="", xlab="", ylab="", xlim=c(0,150),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
title(ylab="Probability density", line=2.5, cex.lab=1.5)
title(xlab="Parameter 2", line=1.6, cex.lab=1.5)
polygon(density(newyy2), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=format(meanyy2, nsmall=4))

##
plot(density(newyy3), main="", xlab="", ylab="", xlim=c(400,500),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
#title(ylab="Probability density", line=1.8, cex.lab=1.5)
title(xlab="Parameter 3", line=1.6, cex.lab=1.5)
polygon(density(newyy3), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=format(meanyy3, nsmall=4))

##
plot(density(newyy4), main="", xlab="", ylab="", xlim=c(0,2e4),
     cex.axis=1.0, yaxs='i', xaxs='i'
     )
#title(ylab="Probability density", line=1.8, cex.lab=1.5)
title(xlab="Parameter 4", line=1.6, cex.lab=1.5)
polygon(density(newyy4), col=adjustcolor("darkgoldenrod2", alpha=0.5))
legend("topright", legend=format(meanyy4, nsmall=4))

dev.off()


######################################################################
# PLOT 4: correlations - NEEDS WORK
######################################################################
pdf("evolution5.pdf", width=6, height=6, onefile=F)

# xxxxxxxx size????
samplesmat <- cbind(newyy1, newyy2, newyy3, newyy4)
pairs(~newyy1+newyy2+newyy3+newyy4, col=adjustcolor("red", alpha=0.05),  
   data=samplesmat[sample(nrow(samplesmat), size=4000, replace=FALSE),], 
   main="Simple Scatterplot Matrix")


dev.off()
}
