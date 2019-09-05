      include "loss.f"
      include "pikaia.f" 
      include "toy.f"
c==================================================================  
      program main
c==================================================================
      implicit none

c     External dependencies      

      integer n, seed, i, status, mdata
      parameter (n=4)    ! number of parameters
      real*8 ctrl(12), x(n), f, loss, fit(n), obs(100), pred(100)
      real*8 obserr(100)

      COMMON/DATA/obs,obserr,mdata
      
c     next line needed to pass a user-defined function name as an argument
      external loss

c
c==================================================================
c     read input
c==================================================================
      open(8, file="main.in", status="unknown")
      read(8,*) seed
      read(8,*)
      do 10 i=1,12
         read(8,*) ctrl(i)
   10 continue
      read(8,*)

      do 11 i=1,500
         read(8,*,end=999) obs(i),obserr(i)
   11 continue
  999 continue
      mdata=i-1      ! number of observation
      
c     initialize the random-number generator
      call rninit(seed)

c==================================================================  
c      - Array  x(1:n)  is the "fittest" (optimal) solution found,
c        i.e., the solution which maximizes fitness function ff
c      - Scalar  f  is the value of the fitness function at x
c      - Integer  status  is an indicator of the success or failure
c        of the call to pikaia (0=success; non-zero=failure)

c     Now call pikaia
      call pikaia(loss,n,ctrl,x,f,status)      
      
c==================================================================  
c     Scale results back to original scale
      fit(1) = x(1) * 1.d4
      fit(2) = x(2) * 1.d4
      fit(3) = x(3) * 1.d4
      fit(4) = x(4) * 1.d5
   
c     calculate prediced values from best parameters
      call toy(n,fit,pred)

c==================================================================  
c     output
      open(9,file="main.out",status="unknown")

      write(9,*) 'status: ',status
      write(9,*)
      write(9,*) 'fittest parameter set: '
      write(9,*) fit(1)
      write(9,*) fit(2)
      write(9,*) fit(3)
      write(9,*) fit(4)
      write(9,*)
      write(9,*) 'function f: '
      write(9,'(1pe11.3e3)') f
      write(9,*)
      write(9,*) 'Observations..Error...............Predictions'
      do 12 i=1,mdata
         write(9,9001) obs(i),obserr(i),pred(i)
   12 continue
    
c==================================================================  
 9001 format(1pe12.4e2,2x,1pe12.4e2,8x,1pe12.4e2)
      close(8)
      close(9)
      stop
      end
