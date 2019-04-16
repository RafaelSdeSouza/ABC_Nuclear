c======================================================================
c     ===================================================             |
      program toy
c
c     this taked the place of nucleo.f. It reads parameters from an
c     input file and outputs simulated values from these parameters   
c
c     with the following 4 parameters:
c
c     34         ! par_1
c     12         ! par_2
c     456        ! par_3
c     78         ! par_4
c
c     we obtain the 5 simulated values:
c
c     1.67E+03   ! val1
c     1.67E+02   ! val2
c     1.67E+01   ! val3
c     1.67E+00   ! val4
c     1.67E-01   ! val5
c
c======================================================================
c   VARIABLES
      integer npar
      real*8 par(100)
      real*8 val1,val2,val3,val4,val5
      
c======================================================================
c   PARAMETER INPUT
      open(601, file='toy.in',status='unknown')
      read(601,1200) npar       ! number of parameters
c      write(6,*) npar
      do 301 i=1,npar
         read(601,*) par(i)
 301  continue      
      close(601)
      
c======================================================================
c   COMPUTE PREDICTED VALUES: BLACK BOX
c======================================================================
c     the numbers in the following equations are "hidden"
      val1=1.d-1*par(1)**4.0d0+
     &     2.d-1*par(2)**3.0d0+
     &     7.d-1*par(3)**2.0d0+
     &     4.d-1*par(4)
      val2=3.d-2*par(1)**4.0d0+
     &     5.d-2*par(2)**3.0d0+
     &     9.d-2*par(3)**2.0d0+
     &     1.d-2*par(4)
      val3=7.d-3*par(1)**4.0d0+
     &     1.d-3*par(2)**3.0d0+
     &     4.d-3*par(3)**2.0d0+
     &     2.d-3*par(4)
      val4=2.d-4*par(1)**4.0d0+
     &     5.d-4*par(2)**3.0d0+
     &     8.d-4*par(3)**2.0d0+
     &     8.d-4*par(4)
      val5=9.d-5*par(1)**4.0d0+
     &     6.d-5*par(2)**3.0d0+
     &     3.d-5*par(3)**2.0d0+
     &     1.d-5*par(4)
      
c======================================================================
c   OUTPUT
     
      open(602, file='toy.out',status='unknown')
         write(602,1100) val1 
         write(602,1100) val2 
         write(602,1100) val3 
         write(602,1100) val4 
         write(602,1100) val5      
      close(602)

 1100 format(1pe9.2)
 1200 format(i3)
      end
