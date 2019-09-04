c==================================================================  
      subroutine readinput
c==================================================================  
c     reads parameters for: GA, nucleo 
c     outputs fittest solution overall
c
c     note: output of fittest solution of each generation is
c           directly done from pikaia.f
c
c     n:      number of parameters
c     x:      vector of parameters
c     y:      vector of predicted values
c     obs:    vector with observations
c     obserr: vector with observational uncertainties
c     zz:     vector with observed abundance atomic number
c     mdata:  number of observational data points
c     nis:    number of isotopes in network
c==================================================================

      implicit none

      integer i, mdata, iflag, j, nn, kkk, L, jj
      integer zz(100),nnn,na,nb,nc,nd,ispe,seed
      integer ia,ib,ic,id
      integer kk1(10000),kk2(10000),kk3(10000),kk4(10000)
      integer kk5(10000),kk6(10000),kk7(10000),kk8(10000)
      integer nucin(50)
      integer nre,nis,nsl
      integer nweakr
      
      real*8 ctrl(12)
      real*8 obs(100),obserr(100),bnd_l(50),bnd_h(50)
      real*8 dummy
      real*8 ratw(1000,60,12)
      real*8 rat(10000,60)
      real*8 annn(1000),znnn(1000),aa1,aa2
      real*8 ANOOO(1000),ZNOOO(1000),XINITTT(1000)
      real*8 nurin(50)
      real*8 xx_norm,xx_norm_err,help
      
      character*2 cha_ec(10000)     
      character*1 cha_w(10000)      
      character*3 charz
      character*80 cdummy
      character*5 xyzzz(1000),onnn(1000)
      character*5 SPE1,SPE2,SPE3,SPE4,emitspeee(10000)
      character*25 nuchar(50)
      character*30 reacstrw(1000)
      character*70 reacstr(10000)
      
c     pass parameters to other units:
      COMMON/PIKA/ctrl,seed
      COMMON/DATA/obs,obserr,zz,mdata,iflag
ccc   parameter bounds
      COMMON/PARA/bnd_l,bnd_h,nis
ccc   running parameters from nuclei.in            
      COMMON/NUCINP/nurin,nucin,nuchar
ccc   rate input
      COMMON/RATES1/rat
      COMMON/RATES2/reacstr,cha_ec,cha_w
      COMMON/WEAK1/ratw
      COMMON/WEAK2/reacstrw,nweakr
ccc   nucleo.dat input
      COMMON/NDAT1/onnn,annn,znnn,
     &   kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8
      COMMON/NDAT2/xyzzz,ANOOO,ZNOOO,XINITTT,emitspeee

c==================================================================
c     read input from main.in
c==================================================================
      open(81, file="main.in", status="unknown")
      read(81,*) seed
      read(81,*)
      do 10 i=1,12
         read(81,*) ctrl(i)
   10 continue
      read(81,*)
      read(81,*)
      read(81,*) bnd_l(1),bnd_h(1)
      read(81,*) bnd_l(2),bnd_h(2)
      read(81,*) bnd_l(3),bnd_h(3)
      read(81,*) bnd_l(4),bnd_h(4)
      read(81,*) bnd_l(5),bnd_h(5)
      read(81,*) bnd_l(6),bnd_h(6)
      read(81,*) bnd_l(7),bnd_h(7)
      read(81,*)
      read(81,*)
      read(81,*) iflag
      read(81,*)
      read(81,*)     
 
c     read isotope name and experimental mass fraction
      i=1    ! data label
 71   read(81,*,end=70) charz,obs(i),obserr(i)
      if(charz(1:1).eq.'!')then
          goto 71
       else
         read(charz(1:3),'(i3)') zz(i)
         i=i+1
         goto 71
      endif

 70   continue
      mdata=i-1  ! number of data points
 
c     if iflag.ne.0, mass fraction ratios are analyzed; find
c     normalizing element first 
 
c     find normalizing element Z and abundance 
      xx_norm=0.d0
      if(iflag.ne.0)then
         do 115 jj=1,mdata 
            if(iflag.eq.zz(jj))then 
               xx_norm=obs(jj)
               xx_norm_err=obserr(jj)   
            endif           
 115     continue

c        check if normalizing element is found; if not, stop run       
         if(xx_norm.eq.0.d0)then
            write(6,*) 'Normalizing element not found. Run stop!'
            stop
          elseif(xx_norm.lt.0.d0)then
            write(6,*) 'Cannot normalize to upper limit. Run stop!'
            stop                    
         endif  

c        if normalizing element is found, normalize
         do 112 j=1,mdata
            help=(obserr(j)/obs(j))**2.d0 +
     &           (xx_norm_err/xx_norm)**2.d0 
            obs(j)=obs(j)/xx_norm
            obserr(j)= obs(j) * dsqrt(help)            
 112     continue
      endif
   
c==================================================================
c     read input from nucleo.in
c==================================================================
      open(98,file='nucleo.in',status='unknown') 
      
      read(98,*)
      read(98,*) nucin(1) 
      read(98,*) nucin(2)
      read(98,*)
      read(98,*)
      read(98,*) nuchar(1)
      read(98,*) nuchar(2)
      read(98,*) nuchar(3)
      read(98,*)
      read(98,*) nucin(3)
      read(98,*) nurin(1),nurin(2),nurin(3)
      read(98,*) nuchar(4),nurin(4),nurin(5),nurin(6),nurin(7)
      read(98,*) nurin(8),nurin(9),nurin(10),nurin(11),nurin(12)
      read(98,*) nuchar(5),nurin(13),nurin(14)
      read(98,*) nuchar(6),nurin(15),nurin(16)
      read(98,*)
      read(98,*) nucin(4)
      read(98,*) nurin(17)
      read(98,*) nurin(18)
      read(98,*) nucin(5)
      read(98,*) nurin(19),nucin(6)
      read(98,*) nurin(20)
      read(98,*) nurin(21)
      read(98,*) nurin(22)
      read(98,*)
      read(98,*) nucin(7)
      read(98,*)
      read(98,*) nucin(8)
      read(98,*) nucin(9),nucin(10)
      read(98,*) (nurin(i),i=23,29)
      read(98,*)
      read(98,*) nucin(11)
      
      close(98)

c==================================================================
c     read nucleo.dat input 
c==================================================================
      OPEN(93,FILE=nuchar(1),status='unknown') 

      READ(93,8070) nis
      READ(93,8100) cdummy

c---- read nuclides in network
      DO 40 I=1,nis
         READ(93,9025) KKK,ONNN(I),ANNN(I),ZNNN(I)
   40 continue
      ONNN(nis+1)='    '
      annn(nis+1)=1.d0
      znnn(nis+1)=0.d0

c---- read initial abundances
      READ(93,8070) nsl
      DO 80 I=1,NSL
         READ(93,8050) xyzzz(i),ANOOO(I),ZNOOO(I),XINITTT(I)
c         write(6,*) xyzzz(i)
   80 CONTINUE

c---- read interaction links
      READ(93,8070) nre
      DO 93 L=1,nre
         READ(93,8060) NNN,NA,SPE1,NB,SPE2,ND,SPE4,NC,SPE3
c         write(6,*) spe1
         IA=ISPE(nis,onnn,SPE1)
         IB=ISPE(nis,onnn,SPE2)
         IC=ISPE(nis,onnn,SPE3)
         ID=ISPE(nis,onnn,SPE4)
         KK1(L)=IA
         KK2(L)=NA
         KK3(L)=IB
         KK4(L)=NB
         KK5(L)=IC
         KK6(L)=NC
         KK7(L)=ID
         KK8(L)=ND
         AA1=NA*ANNN(KK1(L)) + NB*ANNN(KK3(L))
         AA2=NC*ANNN(KK5(L)) + ND*ANNN(KK7(L))
         if(AA1.NE.AA2) then
            write(6,*) L
            write(6,*) 'WARNING: BARYON NUMBER NOT CONSERVED!'
            stop
         endif

c        store identity of emitted light particle for each reaction;
c        useful for identifying links that give rise to, e.g., neutron
c        emission, etc.

         emitspeee(l)=spe4

   93 CONTINUE

      close(93)
      
 8100 FORMAT(A80)
 8050 FORMAT(5X,A5,1X,F3.0,1x,F3.0,1x,E30.3)
 8060 FORMAT(1X,I4,3X,I2,1X,A5,2X,I2,1X,A5,1X,I2,1X,A5,
     &  2X,I2,1X,A5)
 8070 format(i7)
 9025 FORMAT(I4,1x,A5,1x,F4.0,F4.0)

c==================================================================
c     read input from reactions.dat
c==================================================================
      OPEN(23,FILE=nuchar(2),STATUS='unknown')

      do 530 j=1,nre

c        read first line with format and label
         READ(23,9001) REACSTR(j)
c        reacstr is needed for a number of things
         cha_ec(j)=REACSTR(j)(46:47)  ! is link labeled "ec"? 
         cha_w(j)=REACSTR(j)(48:48)   ! is link a weak or reverse interaction?
                           
         do 540 i=1,60
c           j: reaction; i: temperature
c           rates of given reaction are multiplied by same factor at all T
            read(23,*) dummy,rat(j,i),dummy
c           if rate=0.0: set equal to 1.d-99 for logarithmic interpolation
            if(rat(j,i).eq.0.d0) then
               rat(j,i)=1.0d-99    
            endif   
 540     continue

 530  continue

 9001 format(a70)

c==================================================================
c     read weak rates input 
c==================================================================
      open(27,file=nuchar(3),status='unknown')
      j=1
c     j: reaction; i: temperature; nn: density      
 62   read(27,9002,end=60) reacstrw(j)
      do 63 i=1,60
         read(27,9003) dummy,(ratw(j,i,nn),nn=1,12)
 63   continue
      j=j+1
      goto 62
      
 60   continue
      nweakr=j-1    ! nrw=number of weak rates read
c     write(6,*) nrw
      close(27)
      
 9002 format(a30)
 9003 format(1pe8.2,12(5x,1pe9.3))

      RETURN
      END

C======================================================================
      FUNCTION ISPE(nis,on,NAME)
C======================================================================
      integer nis,ispe,i
      character*5 name,on(1000)
      
      DO 1 I=1,nis+1
         IF(NAME.EQ.ON(I))THEN
            ISPE=I
            RETURN
         ENDIF
 1    CONTINUE
      
      write(6,*) 'WARNING: UNIDENTIFIED ISOTOPE'
      write(6,*) name
      stop

      RETURN
      END
