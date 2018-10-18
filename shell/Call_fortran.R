# Fortran 

# gfortran -shared -m64 -o cn3shell.so cn3shell.f

# Load the fortran code 
if(!is.loaded("cn3shellSub"))
   dyn.load("cn3shell.so") #,Local=FALSE)
   .Fortran("cn3shell.so")


 # Load the fortran code needed to analyse the final abundance data
 if(!is.loaded("finalAbundSub"))
   dyn.load("a_finalAbund.so") #,Local=FALSE)
 .Fortran("finalAbundSub")
 
 
 system("./cn3shell.so")