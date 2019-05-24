read.fortran("syndat2.i3e")

zzfil <- tempfile("syndat2.i3e")
zz <- file(zzfil, "rb")
readBin(zz,"numeric")