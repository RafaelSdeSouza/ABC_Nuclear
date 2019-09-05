      include "loss.f"
      include "pikaia.f" 
      include "readinput.f"
c==================================================================  
      program main
c==================================================================  
c     sets number of parameters, n
c
c     - reads parameters for: GA, nucleo 
c     - outputs fittest solution overall
c
c     note: output of fittest solution of each generation is
c           directly done from pikaia.f
c
c     n:      number of parameters
c     x:      vector of parameters
c     y:      vector of predicted values
c     obs:    vector with observations
c     obserr: vector with observational uncertainties
c     zz:     vector with observed abundance isotope labels
c     mdata:  number of observational data points
c     nis:    number of isotopes in network
c==================================================================
c
c     do not use index "n" for any other quantity
c
      implicit none

      integer seed
      integer n,i,status,mdata,iflag,j
      integer zz(100),nis
      character*5 ciso(1000)
       
c     =============== 
      parameter (n=7)    ! number of parameters
c     =============== 

      real*8 ctrl(12), x(n), f, loss, fit(n) 
      real*8 obs(100),obserr(100),bnd_l(50),bnd_h(50)
      real*8 xf(1000),xele(100)
           
c     pass parameters to other units:
      COMMON/PIKA/ctrl,seed
      COMMON/DATA/obs,obserr,zz,mdata,iflag
c     parameter bounds
      COMMON/PARA/bnd_l,bnd_h,nis
     
c     External dependencies      
c     next line needed to pass a user-defined function name as an argument
      external loss

c==================================================================  
c     read input
      call readinput

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
c     Scale results back to original scale [use here same scaling as in loss.f]
      fit(1) = x(1) * (bnd_h(1) - bnd_l(1)) + bnd_l(1)  ! T (K)
      fit(2) = x(2) * (bnd_h(2) - bnd_l(2)) + bnd_l(2)  ! rho (g/cm3)
      fit(3) = x(3) * (bnd_h(3) - bnd_l(3)) + bnd_l(3)  ! 1/e time for T
      fit(4) = x(4) * (bnd_h(4) - bnd_l(4)) + bnd_l(4)  ! 1/e time for rho
      fit(5) = x(5) * (bnd_h(5) - bnd_l(5)) + bnd_l(5)  ! xwd_c12_min 
      fit(6) = x(6) * (bnd_h(6) - bnd_l(6)) + bnd_l(6)  ! xwd_ne22_min 
      fit(7) = x(7) * (bnd_h(7) - bnd_l(7)) + bnd_l(7)  ! fmix; max=10 means 9% WD admixture
   
c==================================================================  
c     output to main.out
      open(91,file="main.out",status="unknown")

      write(91,*) 'status: ',status
      write(91,*)
      write(91,*) 'fittest parameter set: '
      write(91,9002) fit(1)
      write(91,9002) fit(2)
      write(91,9002) fit(3)
      write(91,9002) fit(4)
      write(91,9002) fit(5)
      write(91,9002) fit(6)
      write(91,9002) fit(7)
      write(91,*)
      write(91,*) 'function f: '
      write(91,'(1pe11.3e3)') f
      write(91,*)
    
c==================================================================  
c     calculate observed and fittest predicted abundances    
     
      call nucleop(fit,xf,ciso)
      call CNshell_ele(nis,ciso,xf,zz,mdata,iflag,xele)   
    
      write(91,*) 'Observations..Error...............Predictions'
      do 12 i=1,mdata
         write(91,9001) obs(i),obserr(i),xele(i)
   12 continue

c==================================================================  

 9001 format(1pe12.4e2,2x,1pe12.4e2,8x,1pe12.4e2)
 9002 format(1pe12.4e2) 
      stop
      end
