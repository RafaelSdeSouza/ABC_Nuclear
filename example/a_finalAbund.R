## produces final abundance plots for output of nucleo.f:
## 

## uses:     a_finalAbund.out

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

# for adding minor tick marks:
library(Hmisc)

# Load the fortran code needed to analyse the final abundance data
if(!is.loaded("finalAbundSub"))
  dyn.load("a_finalAbund.so") #,Local=FALSE)
.Fortran("finalAbundSub")
#dyn.unload("finalAbund.so")

# read data from file
mydatx <- read.table("a_finalAbund.out", header=FALSE, na.strings = "")
#print(mydatx)
# write output to pdf file
pdf(file="a_finalAbund.pdf",width=8,height=5,onefile=F)

# define plot
par(mfcol=c(2,1), mar=c(2.5,6.5,1.0,3.0), oma=c(0.2,0.2,0.2,0.2), tck=0.05, 
    las=1, mgp=c(2.8,0.2,0))

# determine plotranges
xRang = max(mydatx$V2)-min(mydatx$V2)
yRang = max(mydatx$V3)-min(mydatx$V3)
limMult = 0.8
xLim= c( min(mydatx$V2)*limMult , max(mydatx$V2)/limMult )
yLim= c( 1e-10 , max(mydatx$V3)/limMult )

# plot data
# cex controls symbol size
# pch is the symbol id
plot(mydatx$V2, mydatx$V3, main="", 
    cex.lab=1.1, cex.axis=1.0, cex.main=1.0, 
    xlab="", 
  	ylab="Final mass fraction", 
    pch=NA, col=NA, 
    xlim=xLim, ylim=yLim, 
    log="y", yaxt="n", tck=0.07)       

text(mydatx$V2, mydatx$V3, label=mydatx$V1, col='blue', cex=0.4)

mtext(side = 1, "Mass number A", line = 1.5)

axis(side=3, labels=FALSE)
axis(side=4, labels=FALSE)
minor.tick(nx=10, ny=0, tick.ratio=0.5)

aY <- axTicks(2); axis(2, at=aY, label= axTexpr(2, aY), las=2)

######### from web: to change "1 e xx" tick mark labels to "a x 10^y"
 axTexpr <- function(side, at = axTicks(side, axp=axp, usr=usr, log=log),
                    axp = NULL, usr = NULL, log = NULL)
 {
    ## Purpose: Do "a 10^k" labeling instead of "a e<k>"
    ##        this auxiliary should return 'at' and 'label' (expression)
    ## ----------------------------------------------------------------------
    ## Arguments: as for axTicks()
    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  7 May 2004, 18:01
    eT <- floor(log10(abs(at)))# at == 0 case is dealt with below
    mT <- at / 10^eT
    ss <- lapply(seq(along = at),
                 function(i) if(at[i] == 0) quote(0) else
                 substitute(10^E, list(A=mT[i], E=eT[i])))
    do.call("expression", ss)
 }
aY <- axTicks(2); axis(2, at=aY, label= axTexpr(2, aY), las=2)
########


## return output to the terminal
dev.off( )


