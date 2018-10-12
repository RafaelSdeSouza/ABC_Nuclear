       PROGRAM a_finalAbund

C     ----------------------------------------------------------
C     produces a table with final abundances and overabundances,
C     this code will be called by R routine
C     ----------------------------------------------------------

      integer ret
      ret=0
      call finalAbundSub(ret)
      end
      
c     -----------------------------------------------------------
      subroutine finalAbundSub(ret)

c     declarations

      INTEGER K,i,nis,ret
      CHARACTER*(132) INSTRING
      real*8 X(1000), XX(1000)
      real*8 xelem(1000)
      character*(5) ON(1000),onhelp
      character*(2) celem(1000)
      character*(3) cmass(1000)
      character*40 netoutfile,netinfile

c     -----------------------------------------------------------
c     I/O

      netoutfile='nucleo.out'
      netinfile='nucleo.dat'

      OPEN (UNIT=9,FILE=netinfile,STATUS='unknown')
      OPEN (UNIT=10,FILE=netoutfile,STATUS='unknown')
      OPEN (UNIT=11,FILE='a_finalAbund.out',STATUS='UNKNOWN')

c     ------------------------------------------------------------
c     read number of species in network

      read(9,'(i6)') nis
      close(9)
c     ------------------------------------------------------------
c     read final abundances from nucleo.out next; they are in the
c     same order as in nucleo.dat

   20 READ(10,1000,END=25) INSTRING
 1000 FORMAT(A144)

      K = INDEX(INSTRING,'FINAL')
      IF (K .NE. 0) THEN 
         read(10,1001) (ON(i),X(i),i=1,NIS) 
       else
         goto 20
      END IF
 25   continue

   21 READ(10,1000,END=26) INSTRING

      K = INDEX(INSTRING,'OVERAB')
      IF (K .NE. 0) THEN 
         read(10,1001) (ON(i),XX(i),i=1,NIS) 
       else
         goto 21
      END IF
 26   continue

      close(10)

      do 867 i=1,nis
         onhelp=on(i)
c        take care of al-6, al*6, al01, al02, al03
         if(onhelp(4:4).eq.('0').or.onhelp(4:4).eq.('*')
     &                          .or.onhelp(4:4).eq.('-'))then
            onhelp(4:5)='26'
         endif
         celem(i)=onhelp(1:2)
         cmass(i)=onhelp(3:5)  
            
         if(X(i).lt.0.d0)then ! avoid negative abundances [from Gear's method]
            X(i)=0.d0
            XX(i)=0.d0
         endif
 
         if(X(i).eq.0.d0)then ! set to very small value for plotting on a log scale
            X(i)=1.0d-99
         endif

         if(XX(i).eq.0.d0)then ! set to very small value for plotting on a log scale
            XX(i)=1.0d-99
         endif
         
 867  continue      

 1001 format(9(1x,a5,1x,1pe10.2E3))
      
c     ----------------------------------------------
c     output

      do 450 i=1,nis
         write(11,9999) celem(i),cmass(i),x(i),xx(i)
 450  continue

 
 9999 format(1x,a2,3x,a3,3x,1pe10.2E3,3x,1pe10.2E3)
 
c     ----------------------------------------------
      close(11)
c      close(9)
c     ----------------------------------------------

      END













