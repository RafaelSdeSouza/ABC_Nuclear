# Fortran 
load("CN3shell.dat")
     


  dyn.load("CN3shell.so") 
  .Fortran("cn3shellSub")


dyn.load("nuclear.so")
.Fortran("nuclear")


is.loaded("multiply")