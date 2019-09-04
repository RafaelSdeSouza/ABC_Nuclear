      include "CNshell_ele.f"
c==========================================================  
      function loss(n,x)
c==========================================================
c     calls: 
c     - nucleop to produce network output [iso abundances]
c     - CNshell_ele to use nucleo output to find ele abundances
c
c     n:      number of parameters
c     x:      vector of parameters
c     bnd:    parameter bounds for pikaia
c     y:      vector of predicted isotopic abundances
c     obs:    vector with observed elemental abundances
c     obserr: vector with observational uncertainties
c     zz:     vector with observed abundance isotope labels
c     mdata:  number of observational data points
c     nis:    number of isotopes in network
c     ciso:   vector with nuclide labels
c==========================================================

      implicit none

      integer n,i,mdata,nis,iflag,j,nflag
      integer zz(100)
      real*8 x(n),sum,loss,par(n),obs(100),obserr(100)
      real*8 y(nis+1)
      real*8 bnd_l(50),bnd_h(50),xele(mdata)
      character*5 ciso(nis+1)

      COMMON/DATA/obs,obserr,zz,mdata,iflag
ccc   parameter bounds
      COMMON/PARA/bnd_l,bnd_h,nis
      
c     **************** rescale input variables *****************
c     scaling between lower and upper bounds given in input file
c     see PIAKIA user's guide, p. 49;
c     x() has a value between 0.0 and 1.0
      par(1) = x(1) * (bnd_h(1) - bnd_l(1)) + bnd_l(1)  ! T (K)
      par(2) = x(2) * (bnd_h(2) - bnd_l(2)) + bnd_l(2)  ! rho (g/cm3)
      par(3) = x(3) * (bnd_h(3) - bnd_l(3)) + bnd_l(3)  ! 1/e time for T
      par(4) = x(4) * (bnd_h(4) - bnd_l(4)) + bnd_l(4)  ! 1/e time for rho
      par(5) = x(5) * (bnd_h(5) - bnd_l(5)) + bnd_l(5)  ! xwd_c12 
      par(6) = x(6) * (bnd_h(6) - bnd_l(6)) + bnd_l(6)  ! xwd_ne22 
      par(7) = x(7) * (bnd_h(7) - bnd_l(7)) + bnd_l(7)  ! fmix
c     **************** rescale input variables *****************

c---- call nucleo
c     input:  par      [parameters: T, rho,...]
c     output: y; ciso  [final iso abundances; iso labels]
      call nucleop(par,y,ciso)

c---- analyze nucleo output
c     input:  nis, ciso, y, obs, zz, mdata, iflag 
c     output: xele [predicted ele abundaces for all observed elements]
      call CNshell_ele(nis,ciso,y,zz,mdata,iflag,xele)
        
c---- compute loss
      sum = 0.d0
      nflag = 1
c     check if any observational uncertainty is zero; if so,
c     use a different fitness function
      do 200 i=1,mdata
         if(obserr(i).eq.0.d0)then
            nflag=0
         endif
 200  continue

c fitness functions: square to avoid negative differences
ccc   sum = sum + ((y(i))-(obs(i)))**2.d0
ccc   sum = sum + ( ((y(i))-(obs(i)))**2.d0 )/y(i)
ccc   sum = sum + (dlog10(y(i))-dlog10(obs(i)))**2.d0  
ccc   sum = sum + ( (dlog10(y(i))-dlog10(obs(i)))**2.d0 )/dlog10(y(i))   
      
      if(nflag.eq.1)then
         do 210 i = 1,mdata
            sum = sum + ( ( xele(i)-obs(i) )/obserr(i) )**2.d0  
 210     continue
       elseif(nflag.eq.0)then
         do 220 i = 1,mdata
            sum = sum + ( dlog10(xele(i))-dlog10(obs(i)) )**2.d0
 220     continue       
      endif
      
c---- evaluate fitness function
ccc      loss = 1.d0/dsqrt(sum)
      loss = 1.d0/sum

      return
      end


