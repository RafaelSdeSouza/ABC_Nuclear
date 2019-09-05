      include "loss.f"
      include "pikaia.f" 
      include "toy.f"
c==================================================================  
      program shell
c==================================================================
c     External dependencies      

      implicit  none
      integer    n, seed, i, status
	     parameter (n=4)
	     real*8      ctrl(12), x(n), f, loss, fit(n)
      external loss

c     First, initialize the random-number generator
c
      seed = 123456
      call rninit(seed)
c
c     Set control variables (evolve 50 individuals over 100
c     generations, use defaults values for other input parameters)
c
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=128
      ctrl(2)= 50000
c	    ctrl(9)= 0.0
c	    ctrl(10)= 3

c     Now call pikaia
      call pikaia(loss,n,ctrl,x,f,status)      
	  
c     Scale results back to original scale
	     fit(1) = x(1)*1000.
	     fit(2) = x(2)*1000.
	     fit(3) = x(3)*1000
	     fit(4) = x(4)*1000.
	  
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      fit: ',fit
      write(*,*) '      f: ',f
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)
c***************************************************************

      end
