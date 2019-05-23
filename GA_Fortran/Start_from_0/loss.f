
      function loss(y)
c==========================================================
c     Loss function for toy model
c==========================================================
      implicit none
      integer      n,i, M
      real*8       y(5),sum,obs(5),loss

	  obs(1) =  279565.6000
	  obs(2) = 58891.5000
	  obs(3) = 10187.9800 
	  obs(4) = 434.5424
      obs(5) = 126.6128

c---------- 2. compute loss
      sum=0.
      do 1 i=1,5

         sum = sum + (log(y(i))-log(obs(i)))**2
    1 continue
c---------- 3. define fitness
      loss = sqrt(sum)

      return
      end
