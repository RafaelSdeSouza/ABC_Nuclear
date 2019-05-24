c     External dependencies      
	  include "toy.f"
	  include "loss.f"
c==================================================================		  
	  program main
c==================================================================
      implicit  none
      integer    n,M
	real*8     x(4)
	real*8     y(5), loss
	external loss

	  x(1) = 34
	  x(2) = 12
	  x(3) = 500
	  x(4) = 7


c     Now call toy
	  call toy(n,x,y)

c     Now call pikaia
      call pikaia(loss,n,ctrl,x,f,status)

c      

c     Print the results
	write(*,*) '      out: ',y
	write(*,*) '      loss: ',loss(y)
c
      end
c***************************************************************

