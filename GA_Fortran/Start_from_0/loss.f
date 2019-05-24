c     External dependencies      
	  include "toy.f"
c==========================================================	  
      function loss(n,x)
c==========================================================
c     Loss function for toy model
c==========================================================
      implicit none
      integer      n,i
      real*8       x(n),sum,obs(n),loss, par(n), y(5)

	  obs(1) =  279565.6000
	  obs(2) = 58891.5000
	  obs(3) = 10187.9800 
	  obs(4) = 434.5424
      obs(5) = 126.6128
c---------- 1. rescale input variables:
	  par(1)=x(1)*1000.
	  par(2)=x(2)*1000.
	  par(3)=x(3)*1000.
	  par(4)=x(4)*1000.



c----------     Now call toy
	  call toy(n,par,y)
	  
	  			
c---------- 2. compute loss
      sum = 0.
      do 1 i = 1,5

         sum = sum + (log(y(i))-log(obs(i)))**2
    1 continue
c---------- 3. define fitness
      loss = -sqrt(sum)

      return
      end
