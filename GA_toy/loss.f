c==========================================================  
      function loss(n,x)
c==========================================================
c     Loss function for toy model
c     n:     number of parameters
c     x:     vector of parameters
c     y:     vector of predicted values
c     obs:   vector with observations
c     mdata: number of observational data points
c==========================================================
      implicit none
      integer n,i,mdata,nflag
      real*8 x(n), sum, loss, par(n) 
      real*8 obs(100), y(100), obserr(100)

      COMMON/DATA/obs,obserr,mdata
      
c---------- rescale input variables ****
c           the multiplication factor represents
c           the upper bound of the parameter
         par(1) = x(1) * 1.d4
         par(2) = x(2) * 1.d4
         par(3) = x(3) * 1.d4
         par(4) = x(4) * 1.d5

c---------- call toy
      call toy(n,par,y)
  
c---------- compute loss
      sum = 0.d0
      nflag = 1
c     check if any observational uncertainty is zero; if so,
c     use a different fitness function
      do 2 i=1,mdata
         if(obserr(i).eq.0.d0)then
            nflag=0
         endif
    2 continue

c fitness functions: square to avoid negative differences
ccc   sum = sum + ((y(i))-(obs(i)))**2.d0
ccc   sum = sum + ( ((y(i))-(obs(i)))**2.d0 )/y(i)
ccc   sum = sum + (dlog10(y(i))-dlog10(obs(i)))**2.d0  
ccc   sum = sum + ( (dlog10(y(i))-dlog10(obs(i)))**2.d0 )/dlog10(y(i))   

      if(nflag.eq.1)then
         do 1 i = 1,mdata
            sum = sum + ( ( y(i)-obs(i) )/obserr(i) )**2.d0
c            sum = sum + 
c     &        ( ( dlog(y(i))-dlog(obs(i)) )/dlog(obserr(i)) )**2.d0
    1    continue
       elseif(nflag.eq.0)then
         do 3 i = 1,mdata
            sum = sum + (dlog10(y(i))-dlog10(obs(i)))**2.d0
    3    continue       
      endif

c---------- evaluate fitness function
c      loss = 1.d0/(sum**2.d0)
c      loss = 1.d0/dsqrt(sum)
      loss = 1.d0/sum

      return
      end

