c==========================================================  
      function loss(n,x)
c==========================================================
c     Loss function for toy model
c     n:      number of parameters
c     x:      vector of parameters
c     bnd:    bounds for pikaia
c     y:      vector of predicted values
c     obs:    vector with observations
c     mdata:  number of observational data points
c==========================================================
      implicit none
      integer n,i,mdata
      integer zz(100)
      real*8 x(n), sum, loss, par(n), obs(100), y(500)
      real*8 bnd_l(50),bnd_h(50)

      COMMON/DATA/bnd_l,bnd_h,obs,zz,mdata
      
c---------- **** rescale input variables ****
c     scaling between lower and upper bounds given in input file
c     see PIAKIA user's guide, p. 49
      par(1) = x(1) * (bnd_h(1) - bnd_l(1)) + bnd_l(1)  ! T (K)
      par(2) = x(2) * (bnd_h(2) - bnd_l(2)) + bnd_l(2)  ! rho (g/cm3)
      par(3) = x(3) * (bnd_h(3) - bnd_l(3)) + bnd_l(3)  ! 1/e time for T
      par(4) = x(4) * (bnd_h(4) - bnd_l(4)) + bnd_l(4)  ! 1/e time for rho

c---------- call nucleo
      call nucleo(par,y)

c---------- analyze nucleo output
      
      stop
      
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
      
      
  
c---------- compute loss
      sum = 0.d0
      do 1 i = 1,mdata
         sum = sum + ((y(i))-(obs(i)))**2    ! square residual
    1 continue

c---------- determine fitness
      loss = 1/dsqrt(sum)

      return
      end


