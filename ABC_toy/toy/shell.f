c======================================================================
c     ===================================================             |
      program shell
c
c
c
c======================================================================
c   VARIABLES
      implicit double precision (a-h, o-z), integer (i-n)

      integer nsamp,nseed,k,ndat,npar
      real*8 par_min(100),par_max(100)
      real*8 lgpar_min(100),lgpar_max(100)
      real*8 par(100),lgpar(100)
      real*8 dat(100),sim(100)     ! dat: oberved; sim: simulated 
      real*8 ratio_max,ratio_new
 
      character*100 string


c======================================================================
c   DATA INPUT

      npar=4      ! number of parameters

      open(601, file='shell.dat',status='unknown')

      read(601,6000) string
      read(601,6000) string
      read(601,*) par_min(1),par_max(1)
      read(601,*) par_min(2),par_max(2)
      read(601,*) par_min(3),par_max(3)
      read(601,*) par_min(4),par_max(4)
      read(601,6000) string
      read(601,*) nsamp,nseed
      read(601,6000) string
      read(601,6000) string
 
      j=1 
 104  read(601,*,end=111) dat(j)
      j=j+1
      goto 104
 111  continue
      ndat=j-1    ! number of data points or simulation points
c      write(6,*) ndat
 
      close(601)     
      
 6000 format(a100)
      
c======================================================================  

c     for sampling on log scale 
      lgpar_min(1)=dlog10(par_min(1))
      lgpar_max(1)=dlog10(par_max(1))

      lgpar_min(2)=dlog10(par_min(2))
      lgpar_max(2)=dlog10(par_max(2))

      lgpar_min(3)=dlog10(par_min(3))
      lgpar_max(3)=dlog10(par_max(3))

      lgpar_min(4)=dlog10(par_min(4))
      lgpar_max(4)=dlog10(par_max(4))
      
c======================================================================
c   OUTPUT FILES

c     for simulated  parameters
      open(710,file='shell.out',status='unknown') 
      write(710,6003)

 6003 format(5x,'par(1)',7x,'par(2)',7x,'par(3)',7x,'par(4)',5x,
     &   'ratio_max')

c======================================================================
c   LOOP OVER EACH NETWORK CALCULATION   --->
c======================================================================
      do 100 k=1,nsamp
c        print counter
c         write(6,*) ' Run # ',k    ! do not use symbol 'k' anywhere 
c                                     else in this loop !
c        ==============================================================
c        randomize profile parameters
         lgpar(1)=lgpar_min(1)+(lgpar_max(1)-lgpar_min(1))*ran1(nseed)
         par(1)=10.d0**lgpar(1)

         lgpar(2)=lgpar_min(2)+(lgpar_max(2)-lgpar_min(2))*ran1(nseed)
         par(2)=10.d0**lgpar(2)

         lgpar(3)=lgpar_min(3)+(lgpar_max(3)-lgpar_min(3))*ran1(nseed)
         par(3)=10.d0**lgpar(3)

         lgpar(4)=lgpar_min(4)+(lgpar_max(4)-lgpar_min(4))*ran1(nseed)
         par(4)=10.d0**lgpar(4)

c        ==============================================================

c        enter sampled profile parameters into toy.in input file 
         call SYSTEM('cp toy.in toy_help.in')

         open(801,file='toy.in',status='unknown')
 
c        construct new toy.in
         write(801,6005) npar
         do 103 j=1,npar
            write(801,6002) par(j)  
 103     continue
 
         close(801)

 6002    format(1pe10.2)  
 6005    format(i3)

c        =============================================================
c        run toy model     
         call SYSTEM('./toy')
c        =============================================================
         open(707,file='toy.out',status='unknown')       
c        read results
         do 117 j=1,ndat
            read(707,*) sim(j)
            if(sim(j).le.0.d0)then
               write(6,*) 'negative simulated values. stop.'
               stop
            endif
 117     continue
 
         close(707)

c        =============================================================
c        here is where we compared data, dat(i), with simulations, 
c        sim(i), this can be changed...
c        =============================================================
         ratio_max=1.0d0

         do 113 j=1,ndat
            ratio_new=dat(j)/sim(j)
            if(ratio_new.lt.1.0d0)then
c              inverse to compare factors >1
               ratio_new=1.d0/ratio_new     
            endif
c            write(6,*) dat(j),sim(j),ratio_new
 113     continue
         if(ratio_new.ge.ratio_max)then
            ratio_max=ratio_new          ! maximum deviation
         endif

c        output sampled abundances and parameters if agreement is 
c        achieved for all abundances within a factor of xx
         if(ratio_max.le.5.0d0)then
c           parameters
            write(710,6033) (par(j),j=1,npar),ratio_max
         endif

 6033    format(2x,25(1pe10.2e3,3x))

c        =============================================================
c        replace temporary toy.in by help file
         call SYSTEM('rm toy.in')
         call SYSTEM('mv toy_help.in toy.in')

c======================================================================
 100  continue     ! <---- NEXT NETWORK SAMPLE
c======================================================================      
      close(710)
            
c======================================================================
      end

C     ====================================================================
      FUNCTION ran1(idum)
C     ====================================================================
      integer idum,ia,im,iq,ir,ntab,ndiv
      real*8 ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1.d0/im,iq=127773,ir=2836,
     &    ntab=32,ndiv=1+(im-1)/ntab,eps=1.2d-7,rnmx=1.d0-eps)
C         "Minimal" random number generator of Park and Miller with Bays-
C         Durham shuffle and added safeguards. Returns a uniform random
C         deviate between 0.0 and 1.0 (exclusive of the endpoint values).
C         Call with idum a negative integer to initialize; thereafter, do
C         not alter idum between successive deviates in a sequence. RNMX
C         should approximate the largest floating value that is less than 1.
C         (routine taken from NUMERICAL RECIPES); do not call more frequently 
C         than 5% of its period (i.e., number of calls<100,000,000) 
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do 11 j=ntab+8,1,-1
              k=idum/iq
              idum=ia*(idum-k*iq)-ir*k
              if (idum.lt.0) idum=idum+im
              if (j.le.ntab) iv(j)=idum
  11     continue
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)
      
      return
      end


