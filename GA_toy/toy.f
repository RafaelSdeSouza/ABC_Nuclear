      subroutine toy(n,par,y)
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
c     2.80E+05   ! val1
c     5.89E+04   ! val2
c     1.02E+04   ! val3
c     4.35E+02   ! val4
c     1.27E+02   ! val5
c
c======================================================================
c   VARIABLES
      integer n,m
      real*8 par(n)
      real*8 y(m)
      
c   COMPUTE PREDICTED VALUES: BLACK BOX
c======================================================================
c     the numbers in the following equations are "hidden"
      y(1)=1.d-1*par(1)**4.0d0+
     &     2.d-1*par(2)**3.0d0+
     &     7.d-1*par(3)**2.0d0+
     &     4.d-1*par(4)
      y(2)=3.d-2*par(1)**4.0d0+
     &     5.d-2*par(2)**3.0d0+
     &     9.d-2*par(3)**2.0d0+
     &     1.d-2*par(4)
      y(3)=7.d-3*par(1)**4.0d0+
     &     1.d-3*par(2)**3.0d0+
     &     4.d-3*par(3)**2.0d0+
     &     2.d-3*par(4)
      y(4)=2.d-4*par(1)**4.0d0+
     &     5.d-4*par(2)**3.0d0+
     &     8.d-4*par(3)**2.0d0+
     &     8.d-4*par(4)
      y(5)=9.d-5*par(1)**4.0d0+
     &     6.d-5*par(2)**3.0d0+
     &     3.d-5*par(3)**2.0d0+
     &     1.d-5*par(4)
      
      return
      end

