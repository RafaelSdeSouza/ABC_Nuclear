# Fortran 
load("CN3shell.dat")
     

  dyn.load("cn3shell.so")

 .Fortran("cn3shell")


 
 if(!is.loaded("cn3shell_Sub"))
   dyn.load("cn3shell.so") #,Local=FALSE)
 .Fortran("cn3shell")
 
 
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