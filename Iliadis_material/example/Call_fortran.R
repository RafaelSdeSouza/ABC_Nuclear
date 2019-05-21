# Fortran 
load("CN3shell.dat")
     

  dyn.load("cn3shell.so") #,Local=FALSE)
  is.loaded("cn3shellSub")
 .Fortran("cn3shellSub")


  dyn.load("cn3shell.so") 
  .Fortran("cn3shell",10,-1)

  myfacto <- function() {
    dyn.load("cn3shell.so")
    retvals <- .Fortran("cn3shell")
    return(retvals)
  }
  
  
  
dyn.load("nuclear.so")
.Fortran("nuclear")


is.loaded("cn3shell_")