C======================================================================
      subroutine nucleop(parm,x,ciso) 
C     versionnum = "2.8"
C     versiondate='August 17, 2017'
C
C     X: mass fractions of all nuclides in network
C
ccc   all adaptations/changes are labeled with "ccc"
C======================================================================
ccc   this script needs to be adapted to the astrophysical site
ccc   under study; search for "change below part"
ccc
ccc   SETUP HERE IS FOR NOVAE
ccc

      include "comgear.inc"
      
      DIMENSION X(NSP),xold(nsp),XS(NSP),XST(KP,NSP)

      integer stepstatus
      integer nucin(50)
      real*8 nurin(50)      
      real*8 parm(7)       ! model parameters
      character*25 nuchar(50)
      character*5 ciso(nsp)

ccc   running parameters from nuclei.in            
      COMMON/NUCINP/nurin,nucin,nuchar

c            
      ishrun=1			! shell number [for profile]

ccc   pass parameters to START
      CALL START(nucin,nurin,nuchar,parm,X)  

   11 continue
     
      T0=0.d0
      T=0.d0
      NPRINT=0
      NT=0
      stepstatus=0
      delta=delta0
      ncurrtime=0

      do i=1,nsp
         xold(i)=0.d0
      enddo

C     Setup the MA48 matrix etc
      IF((INMODE.EQ.2).or.(INMODE.EQ.3).or.(INMODE.EQ.5)) CALL PROFIL(T)
      call cpu_time(stime)
      call srate
      call cpu_time(ftime)
      ratetime=ratetime+(ftime-stime)
      extime=ftime-stime

c     Initialize the sparse matrix
      call initMA48(x)

      nstep=1
      call cpu_time(amstime)

C=============== LOOP FOR EACH RUN ==========>
 10   CONTINUE

      IF(DEBUG .eq. 1)then
         WRITE(6,'(/,"*****************************",/,a,i6)')
     &        'Step:  ',nstep
         WRITE(6,'("Time = ",es11.3," dT = ",es11.3)')T+DELTA,DELTA
         WRITE(6,'("Tau = ",es11.3," Rho = ",es11.3)')tau,rho
         WRITE(6,'("Resteps = ",i4)')nrestep
      ENDIF

C     Make sure the profile is up to date
      IF((nsolver .ne. 0) .and. ((INMODE.ne.1).or.(INMODE.ne.3))) THEN
         if(inmode .eq. 5)then
            call locate(tmodel,numod,T+delta,ncurrtime)
            if(DEBUG .eq. 1)
     &           print*,'Time Counter = ',ncurrtime

         else
            call locate(tprof,numod,T+delta,ncurrtime)
            if(DEBUG .eq. 1)
     &           print*,'Time Counter = ',ncurrtime
         endif

         if(inmode .ne. 4)then
            CALL PROFIL(T+delta)
            call cpu_time(stime)
            CALL SRATE          ! calculate reaction rates at Tau-Rho...
            call cpu_time(ftime)
         else
            call cpu_time(stime)
            call searchmodel(T+delta)
            call instmix(X)
            call cpu_time(ftime)
         endif
         ratetime=ratetime+(ftime-stime)
      ENDIF

C     TAKE A STEP!
      if(nsolver .eq. 0) then
c        print*,delt,dtold,delta
         CALL SNUC(X,SUMX,DTOLD) ! nucleosynthesis for Delta (time)
         IF(DEBUG .eq. 1)then
             WRITE(6,'(a,f12.3)')'New dT = ',DTOLD
             WRITE(6,'(a,i3)') 'Step status: ',stepstatus
         ENDIF        
      else if(nsolver .eq. 1) then
         IF(DEBUG .eq. 1)then
            WRITE(6,'(i3,a,i6)') nordsteps+1,' steps of order: ',
     &           nordord
         ENDIF
         CALL gearstep(gres,X,SUMX,DTOLD,T,stepstatus)
      endif

c     if the step failed, return for another try
      if(stepstatus .eq. -1)goto 10

      T=T+DTOLD                                 ! new time = old time + old timestep

      if(nsolver .eq. 0)then
         IF((T+DELTA).GE.TLAST) DELTA=(TLAST-T)/1.1d0
         DELTA=DMIN1(DELTA,TLAST/DELFAC)        ! new timestep; not too large...

         DO 1273 I=1,NR
            TOTFLUX(I)=TOTFLUX(I)+DTOLD*FLUX(I)	! FLUX(I) are reaction fluxes
            tenergy(i)=tenergy(i)+dtold*ener(i)
 1273    CONTINUE 
         
         TOTENER=0.D0
         DO 1275 I=1,NR
            TOTENER=TOTENER+ENER(I)
 1275    CONTINUE
         TOTENERGY=TOTENERGY+DTOLD*TOTENER

      else
         call fluxenergy(X,dtold)
      endif
     
      if(inmode.eq.5) then
         tenersh(ishrun)=totenergy
      endif

      NPRINT=NPRINT+1
      IF(NPRINT.EQ.1) THEN                      ! stocking results every Nprint timesteps
         NPRINT=-(NLAST/KP)-1                   ! other choices possible (e.g. every dt sec..)
         NT=NT+1                                ! next printing step
         CALL STOCK(NT,T,SUMX,X,XST)
      ENDIF

      IF((NSTEP.GE.NLAST).OR.(T.GE.TLAST).or.(X(ISF).LE.XLAST)) THEN        
         NT=NT+1                                ! end of run

c------- exit program for following conditions
         if(inmode.ne.5) then
ccc         do not use "stop", but exit with return
            goto 12
          elseif((inmode.eq.5).and.(ishrun.eq.nshmax)) then
            call meanx
            stop
          elseif(inmode.eq.5) then   ! next shell
             CALL CLEARSHELL(X)
            ishrun=ishrun+1 
            goto 11
         endif                    
      ENDIF
            
c---- continue present run, new time step
      IF(nsolver .eq. 0)then
         if(((INMODE.EQ.2).or.(INMODE.EQ.3).or.(INMODE.EQ.5))) then   
            call PROFIL(T)      ! calculate new Tau-Rho from profile
         elseif(inmode.eq.4) then
            call searchmodel(t) ! calculate new model
         endif
      endif
      
      nstep = nstep+1

      GO TO 10                                  ! next time step      
C <============= LOOP FOR EACH CALCULATION ========     
 9109 format(15x,'Shell #',i4)
    
ccc   exit nucleo
   12 continue
   
ccc   copy nuclide labels into vector ciso
      do 66 i=1,nsp
         ciso(i)=on(i)
   66 continue
   
      return
      END

C=======================================================================
ccc   pass parameters to START 
      SUBROUTINE START(nucin,nurin,nuchar,parm,X)
C-----------------------------------------------------------------------
C     Reading of input parameters and network
C-----------------------------------------------------------------------
C     ITEST=0: running
C     ITEST=1: print network
C     ILAST=0: single run
C     ILAST=3: Monte Carlo sampling
C     INMODE=1: const T,rho
C     INMODE=2: T-rho profile
C     INMODE=3: T-rho parametrization
C     INMODE=4: average rates over mass shells [instant mising]
C     INMODE=5: sequential runs over shells [no mixing]
C
C     only abundance evolutions with Y>YTMIN can be trusted
C======================================================================

      include 'comnuc.inc'
      dimension X(NSP)
      character*25 xx4,xx11,xx13
      character*25 profilefile
      character*25 shellinputfile
      character*5 versionnum
      character*35 versiondate

      integer nucin(50)
      character*25 nuchar(50)
      real*8 nurin(50)

ccc   define vector of parameter values
      real*8 parm(7)

      real*8 dconv
                  
      itest = nucin(1)
      ilast = nucin(2)
      networkfile = nuchar(1)
      reactionfile = nuchar(2)
      weakinputfile = nuchar(3)
      inmode = nucin(3)
      xx1 = nurin(1)
      xx2 = nurin(2)
      xx3 = nurin(3)
      xx4 = nuchar(4)
      xx5 = nurin(4)
      xx19 = nurin(5)
      xx20 = nurin(6)
      xx21 = nurin(7)
      xx6 = nurin(8)
      xx7 = nurin(9)
      xx8 = nurin(10)
      xx9a = nurin(11)
      xx9b = nurin(12)
      xx11 = nuchar(5)
      xx12 = nurin(13)
      xx16 = nurin(14)
      xx13 = nuchar(6)
      xx14 = nurin(15)
      xx15 = nurin(16)
      nsolver = nucin(4)
      eps = nurin(17)
      ytmin = nurin(18)
      nlast = nucin(5)
      xlast = nurin(19)
      isf = nucin(6)
      delta0 = nurin(20)
      gres = nurin(21)
      delfac = nurin(22)
      iweak = nucin(7)
      iflux = nucin(8)
      isop1 = nurin(9)
      isop2 = nucin(10)
      numab(1)= nurin(23)
      numab(2)= nurin(24) 
      numab(3)= nurin(25)
      numab(4)= nurin(26)
      numab(5)= nurin(27)
      numab(6)= nurin(28)
      numab(7)= nurin(29)
      xx10 = nucin(11)

c --- check some input
      if((itest.ne.0).and.(itest.ne.1))then
         write(6,*) 'Check input: ITEST. Run stop!'
         stop
      endif
ccc
      if(ilast.ne.0)then
         write(6,*) 'Check input: ILAST. Run stop!'
         stop
      endif
ccc
      if(inmode.ne.3)then
         write(6,*) 'Check input: INMODE. Run stop!'
         stop
      endif
      if((iweak.ne.0).and.(iweak.ne.1))then
         write(6,*) 'Check input: IWEAK. Run stop!'
         stop
      endif
ccc
      if(iflux.ne.0)then
         write(6,*) 'Check input: IFLUX. Run stop!'
         stop
      endif
      if((xx10.ne.0).and.(xx10.ne.1))then
         write(6,*) 'Check input: xx10. Run stop!'
         stop
      endif
c -------------------------------------

      nyhandler=1
            
      IF(ITEST.EQ.0) THEN
      
C------- constant tau-rho  -----------------------------
         IF(INMODE.EQ.1) THEN
            TAU0=xx1
            RHO0=xx2
			         TLAST=xx3
         ENDIF

C------- reads tau-rho profile from external file --------
         IF(INMODE.EQ.2) THEN
            profilefile=xx4
            TLAST=xx5
            
            tauprofmax=0.0d0
            tauprocalmax=0.0d0
            i=1
 5543       read(24,*,end=3399) tprof(i),tauprof(i),rhoprof(i)
            if (tauprof(i).ge.tauprofmax) then
               tauprofmax=tauprof(i)
            endif         
            i=i+1
            goto 5543
 3399       continue
            numod=i-1    ! numod: # of profile time grid points
            
c           check if time is increasing with each step 
            do 3398 j=2,numod
               if(tprof(j).le.tprof(j-1))then
                  write(6,*) 'PROFILE TIME STEPS NEED TO INCREASE!!'
                  write(6,*) '  CHECK LINE:',j-1
               stop
               endif
 3398       continue  

c           make sure profile starts at time t=0 
            if(tprof(1).ne.0.0d0)then 
               write(6,*) 'PROFILE MUST START AT ZERO TIME!!'
               stop
            endif
                     
            if(tlast.gt.tprof(numod)) then
c               write(6,*) tlast,tprof(numod)
               write(6,*) 'TLAST EXTENDS BEYOND PROFILE!'
               stop
            endif

c           stretch/compress profile if desired by user:
c            xx19: multiplication of T 
c            xx20: multiplication of rho
c            xx21: multiplication of time           
            if(xx19.eq.0.d0)then
                write(6,*) 'CHECK INPUT (T)'
                stop
              elseif(xx20.eq.0.d0)then
                write(6,*) 'CHECK INPUT (rho)'
                stop
              elseif(xx21.eq.0.d0)then
                write(6,*) 'CHECK INPUT (time)'
                stop                  
            endif
            do 3478 j=1,numod
               tauprof(j)=tauprof(j)*xx19
c              if xx20<0: set density equal to constant xx20 value
               if(xx20.lt.0.d0)then
                  rhoprof(j)=dabs(xx20)
                else
                  rhoprof(j)=rhoprof(j)*xx20
               endif
               tprof(j)=tprof(j)*xx21
 3478       continue  
            tlast=tlast*xx21          
            tauprofmax=tauprofmax*xx19
         ENDIF 

C------  parametrized tau-rho -------------------------------------
         IF(inmode.eq.3) THEN
ccc         insert parameters here, but keep running time from nucleo.in
            TAU0=parm(1)
	           RHO0=parm(2)
	           TLAST=xx8
            scale1=parm(3)
            scale2=parm(4)
         ENDIF

C------  multi-model input -----------------------------------------
         IF((inmode.eq.4).or.(inmode.eq.5))THEN
            if(inmode.eq.4)then
               shellinputfile=xx11
               TLAST=xx12
               dconv=xx16
             elseif(inmode.eq.5)then
               shellinputfile=xx13
               TLAST=xx14
               dconv=xx15
            endif
            open(25,file=shellinputfile,status='unknown')
            open(26,file='profile.out',status='unknown')
ccc            write(26,9119)
            open(28,file='instmixtest.out',status='unknown')

ccc            write(28,'(a)') 'Decay constants in 1/s'
ccc            write(28,'(a)') '...output may have been commented out...' 
c           i denotes model [at given time], j denotes shells in given model        
            i=1
            j=1
            
 5547       read(25,*,end=3622) nsh0(i),tmodel(i)
            do 453 j=1,nsh0(i)
               read(25,*) taush0(i,j),rhosh0(i,j),
     &           vconsh0(i,j),clensh0(i,j),dmsh0(i,j)
 453        continue
            
            i=i+1
            goto 5547
 3622       continue 
  
            numod=i-1   ! number of models [# of profile time grid points]

c           find max. number of shells and corresponding model number
            modmaxnum=1
            nshmax=nsh0(1)

            do 599 i=2,numod
               if(nsh0(i).gt.nsh0(i-1))then
                  nshmax=nsh0(i)   ! max. number of shells
                  modmaxnum=i      ! model number with max. shells
               endif
c              check that time profile points are increasing; if not stop run
               if(tmodel(i).le.tmodel(i-1))then
                  write(6,*) 'PROFILE TIME STEPS NEED TO INCREASE!!'
                  write(6,'(a,i4,a,i6)') 
     &                 '  CHECK TIME:',i,' ON LINE:',(i-1)*(1+nsh0(1))+1
                  stop
               endif
 599        continue
 
c           make sure models start at time t=0 
            if(tmodel(1).ne.0.0d0)then
               write(6,*) 'PROFILE MUST START AT ZERO TIME!!'
               stop
            endif

            if(tlast.gt.tmodel(numod)) then
               write(6,*) 'TLAST EXTENDS BEYOND MODELS!'
               stop
            endif
  
c           re-size input shells  
            call newshells(dconv)

         ENDIF

      ENDIF
       
C---- initialize vector 'mult' - the adjustment factors for each rate
      do 333 i=1,nr
        vmult(i)=1.d0
 333  continue

c---- read multiplication factors for each rate
      if (xx10.eq.1) then
         read(98,*) nmult
         do 222 i=1,nmult
		    read(98,*) nreac,vmult(nreac)
 222     continue
      endif
      close(98)

c--------- temperature grid ---------------------------
      temp(1)=0.001d0
      temp(2)=0.002d0
      temp(3)=0.003d0
      temp(4)=0.004d0
      temp(5)=0.005d0
      temp(6)=0.006d0
      temp(7)=0.007d0
      temp(8)=0.008d0
      temp(9)=0.009d0
      temp(10)=0.010d0
      temp(11)=0.011d0
      temp(12)=0.012d0
      temp(13)=0.013d0
      temp(14)=0.014d0
      temp(15)=0.015d0
      temp(16)=0.016d0
      temp(17)=0.018d0
      temp(18)=0.020d0
      temp(19)=0.025d0
      temp(20)=0.030d0
      temp(21)=0.040d0
      temp(22)=0.050d0
      temp(23)=0.060d0
      temp(24)=0.070d0
      temp(25)=0.080d0
      temp(26)=0.09d0
      temp(27)=0.10d0
      temp(28)=0.11d0
      temp(29)=0.12d0
      temp(30)=0.13d0
      temp(31)=0.14d0
      temp(32)=0.15d0
      temp(33)=0.16d0
      temp(34)=0.18d0
      temp(35)=0.20d0
      temp(36)=0.25d0
      temp(37)=0.30d0
      temp(38)=0.35d0
      temp(39)=0.40d0
      temp(40)=0.45d0
      temp(41)=0.50d0
      temp(42)=0.60d0
      temp(43)=0.70d0
      temp(44)=0.80d0
      temp(45)=0.90d0
      temp(46)=1.0d0
      temp(47)=1.25d0
      temp(48)=1.5d0
      temp(49)=1.75d0
      temp(50)=2.0d0
      temp(51)=2.5d0
      temp(52)=3.0d0
      temp(53)=3.5d0
      temp(54)=4.0d0
      temp(55)=5.0d0
      temp(56)=6.0d0
      temp(57)=7.0d0
      temp(58)=8.0d0
      temp(59)=9.0d0
      temp(60)=10.0d0
   
      do 550 i=1,60
         templog(i)=dlog10(temp(i))
 550  continue

      if(iweak.eq.1)then
         rhoyelog(1)=0.0d0
         rhoyelog(2)=1.0d0
         rhoyelog(3)=2.0d0
         rhoyelog(4)=3.0d0
         rhoyelog(5)=4.0d0
         rhoyelog(6)=5.0d0
         rhoyelog(7)=6.0d0
         rhoyelog(8)=7.0d0
         rhoyelog(9)=8.0d0
         rhoyelog(10)=9.0d0
         rhoyelog(11)=10.0d0
         rhoyelog(12)=11.0d0
      endif
C---------------------------------------------------------------              
 
      CALL RPARAM       ! Parameters for nucleosynthesis
      CALL RNETWORK(parm,X)  ! Reading of Network
      
C---------------- INITIALIZATION -------------------------------
      NSTEP=0
      NREAL=0
      njacb=0   ! Jacobian build and inversions
      njace=0   ! Jacobian evaluations
      nfunc=0   ! RHS function evaluations
      nfail=0   ! Failed steps
      ncor=0    ! Newtonian corrector steps
      nmakejac=jacrate ! only invert jacobian every jacrate steps
      EN=0.d0
      ER=0.d0
      oldtau=0.d0
      oldrho=0.d0
      oldye=0.d0

      ncurrtime=1
      
      do i=1,nis
         TotErr(i)=0.d0
      enddo

      TOTENERGY=0.D0
      DO 1264 I=1,NR
 1264 TOTFLUX(I)=0.D0
C------------------------------------
      IF(INMODE.EQ.1) THEN
         TAU=TAU0
         RHO=RHO0
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
ccc         write(9,9952) tau0,rho0
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
      ENDIF
C------------------------------------
      IF(INMODE.EQ.2) THEN
         CALL PROFIL(T)   ! calculate initial TAU, RHO from profile
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif	 
ccc         WRITE(9,9953) profilefile
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
      ENDIF
C------------------------------------            
      if(INMODE.eq.3) then
		       tau=tau0
	        rho=rho0
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
ccc         WRITE(9,9954) TAU0,RHO0,scale1,scale2 
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
      endif
C------------------------------------
      IF(INMODE.EQ.4) THEN
         call searchmodel(t)  ! calculate initial shell properties
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif	 
ccc         WRITE(9,9957) shellinputfile
ccc         if(ilast.eq.0) then
ccc            write(9,9300)
ccc         endif
      ENDIF
C------------------------------------
      IF(INMODE.EQ.5) THEN
         CALL PROFIL(T)      ! calculate initial TAU, RHO from profile
ccc         write(9,9300)
ccc         WRITE(9,9958) shellinputfile
ccc         write(9,9300)
      ENDIF
C------------------------------------           
ccc      if((ilast.eq.0).and.(inmode.ne.5)) then
ccc         write(9,9962) eps,ytmin,nlast
ccc      endif

ccc      write(9,9955) networkfile
ccc      write(9,9956) reactionfile

ccc      if((ilast.eq.0).and.(inmode.ne.5)) then
ccc         write(9,9322) nis,nr
ccc	     write(9,9300)
ccc      endif
C------------------------------------
      do 676 i=1,nr
ccc         if ((vmult(i).ne.1.d0).and.(ilast.eq.0)) 
ccc     &     write(9,9971) i,vmult(i)
 676  continue

ccc      if((ilast.eq.0).and.(inmode.ne.5)) then 
ccc         write(9,*)
ccc       elseif((ilast.eq.0).and.(inmode.eq.5)) then
ccc         write(9,9300)
ccc      endif

ccc      if((ilast.eq.0).and.(inmode.ne.5)) then
ccc         WRITE(9,9321) NLAST/KP+1
ccc      endif
C------------------------------------

 9106 format(a) 
 9119 format('   nstep',3x,'Time (s)',2x,'Temp (K)',2x,' rho (g/cm3)')
 9300 format(153('$'))	 
 9321 FORMAT(/,20X,'PRINTING EVERY NLAST/KP=',I7,'  TIME STEPS',/)
 9322 format(1x,'    number of nuclides=',i5,
     &   '        number of reactions=',i6)
 9952 FORMAT(1X,'RUN AT CONSTANT TEMPERATURE=',1PE9.2,
     & ' K  AND DENSITY=',1PE9.2,' GR/CM3')
 9953 FORMAT(1X,'RUN WITH TAU-RHO PROFILE FROM FILE=',a)
 9954 FORMAT(1X,'RUN AT EXPONENTIAL TEMPERATURE=',
     & 1PE9.2,' K, DENSITY=',1PE9.2,' GR/CM3, SCALE FACTOR 1=',
     & 1PE7.1,' AND SCALE FACTOR 2=',1PE7.1)
 9955 format(1x,'    networkfile= ',a)
 9956 format(1x,'    reactionfile= ',a)
 9957 FORMAT(1X,'RUN WITH "INSTANT MIXING" FROM FILE=',a)
 9958 FORMAT(1X,'RUN WITH "NO MIXING" FROM FILE=',a)
 9962 format(5x,'EPS=',1PD8.1,',  YTMIN=',1PD8.1,',  NLAST=',I7)   
 9971 format(' ',' Reaction ',i4,' changed by a factor of ',
     &   1PD10.3)
      
      RETURN
      END
C======================================================================
      SUBROUTINE RPARAM
C======================================================================
      include 'comnuc.inc'
      dimension ff(0:6)
      data ff/1.d0,1.d0,2.d0,6.d0,2.4d1,1.2d2,7.2d2/

      PAR=0.25d0
      PAMIN=5.d0
      ATEST=5.d0
      TES=20.d0
      DYMIN=1.D-40
      DEMIN=1.D-40
      YMIN=1.D-35
      YEMINM=YTMIN

      do 10 i=1,7
  10  f(i-1)=ff(i-1)

      RETURN
      END
C======================================================================
      SUBROUTINE RNETWORK(parm,X)
C---------------------------------------------------------------------
C     READS ISOTOPES AND NETWORK
C======================================================================
      include 'comnuc.inc'
      character*5 xyz(nsl),onnn(1000),xyzzz(1000),emitspeee(nr)
      character*7 rspec1(nr),rspec2(nr),rspec3(nr),rspec4(nr)
      character*7 tspec1(nr),tspec2(nr),tspec3(nr),tspec4(nr)
      character*80 dummystr
      character*70 reacstr(10000)
      character*2 cha_ec(10000)     
      character*1 cha_w(10000)      
      real*8 rat(10000,60)
      real*8 annn(1000),znnn(1000)
      real*8 ANOOO(1000),ZNOOO(1000),XINITTT(1000)
      character*30 reacstrw(1000)
      real*8 ratw(1000,60,12),xwd(nsl)
      integer nweakr
      integer kk1(10000),kk2(10000),kk3(10000),kk4(10000)
      integer kk5(10000),kk6(10000),kk7(10000),kk8(10000)
ccc   define vector of parameter values
      real*8 parm(7)

      DIMENSION X(NSP)
      DIMENSION KNET(0:100,0:150)
      DIMENSION ANO(NSL),ZNO(NSL),XNO(NSL),XINIT(NSL)
      DIMENSION VVX(NR,NKT),TTT(NKT)
      INTEGER n,i,j,k,kk

ccc   rate input
      COMMON/RATES1/rat
      COMMON/RATES2/reacstr,cha_ec,cha_w
      COMMON/WEAK1/ratw
      COMMON/WEAK2/reacstrw,nweakr
ccc   nucleo.dat input
      COMMON/NDAT1/onnn,annn,znnn,
     &   kk1,kk2,kk3,kk4,kk5,kk6,kk7,kk8
      COMMON/NDAT2/xyzzz,ANOOO,ZNOOO,XINITTT,emitspeee

c     calculates NA<sv> at these T if network is printed
      DATA TTT/0.0030d0,0.05d0,0.100d0,0.15d0,0.20d0,5.0d0,10.0d0/  ! T in 10**9 K

c     ----------------------------------------------------------------------
c     read reaction rates from REACTIONS.DAT; reactions are ordered exactly 
c     the same as in NUCLEO.DAT
c     ----------------------------------------------------------------------
ccc
      do 530 j=1,nr

         REACSTRING(j) = reacstr(j)
c        reacstring is needed for a number of things
         char_ec(j)=cha_ec(j)  ! is link labeled "ec"? 
         char_w(j)=cha_w(j)    ! is link a weak or reverse interaction?
         prob(j)=0.d0      
                           
         do 540 i=1,60
c           j: reaction; i: temperature
c           rates of given reaction are multiplied by same factor at all T
            rate(j,i)=rat(j,i)
 540     continue

 530  continue
 
c     this is the "flat" parameterization: for each run, a probability
c     factor p is sampled for each reaction; for a given reaction,
c     p has the same value at all temperatures;
c     for INMODE=5, sample only between runs, not between shells; randomly
c     sample all values of p, including for reverse reactions; will be taken
c     care of later
ccc      if(ilast.eq.3)then        
ccc         open(36,file='prob.dat',status='unknown')
ccc         do 206 j=1,nr
ccc            read(36,1011) prob(j) 
ccc 206     continue
ccc         close(36)
ccc      endif
      
 1011 format(7x,1pe20.12)
 
      close(23)
    
c     ----------------------------------------------------------------
c     read weak interaction rates if needed
c     ----------------------------------------------------------------
ccc
      if(iweak.eq.1)then
         nrw = nweakr

         do 60 j=1,nrw
c           j: reaction; i: temperature; n: density      
            reacstrweak(j) = reacstrw(j)
            do 63 i=1,60
               do 64 n=1,12
                  ratew(j,i,n) = ratw(j,i,n)
                  ratelogw(j,i,n)=dlog10(ratew(j,i,n))               
 64            continue
 63         continue
 60      continue
c        write(6,*) nrw
      endif

c     ----------------------------------------------------------------
c     read isotopes and network from NUCLEO.DAT
c     ----------------------------------------------------------------
ccc      OPEN(8,FILE=networkfile,status='unknown') 

c---- nuclides 

ccc      READ(8,8101) ABCD
ccc      READ(8,8100) ABCD
      DO 40 I=1,NIS+1
         ON(I) = onnn(i)
         AN(I) = annn(i)
         ZN(I) = znnn(i)
   40 continue
ccc      ON(NIS+1)='    '
ccc      AN(NIS+1)=1.d0
ccc      ZN(NIS+1)=0.d0

c---- initial abundances 
         do 141 jj=1,nsl
            xwd(jj)=0.0d0
 141     continue

c========================================================
c change below part depending on astrophysical scenario;
c this is for mixing in novae
c========================================================
c        first set WD abundances
         xwd(9)=parm(5)     ! 12C
         xwd(19)=parm(6)    ! 22Ne
         xwd(13)=1.d0-xwd(9)-xwd(19)   ! 16O

ccc      READ(8,8101) ABCD
c     assign initial abundances by mixing with WD matter
      DO 80 I=1,NSL
         xyz(i) = xyzzz(i)
         ANO(I) = ANOOO(i)
         ZNO(I) = ZNOOO(i)
c         XINIT(I) = XINITTT(i)
         XINIT(I) = (xwd(i)+parm(7)*XINITTT(i))/(1.d0+parm(7))
c        calculate initial eta and Ye
ccc         ETA=ETA+(ANO(I)-2.d0*ZNO(I))*(XINIT(I)/ANO(I))
   80 CONTINUE
c========================================================
c     change above part
c========================================================

ccc      YE=0.5d0*(1.d0-ETA)

c     xyz(i): element names + mass number
      DO 90 I=1,NIS+1
        X(I)=0.
        XINI(I)=0.
        DO 85 K=1,NSL  ! setting initial abundances 
          IF((ZNO(K).EQ.ZN(I)).AND.(ANO(K).EQ.AN(I))
     &      .and.(on(i).eq.xyz(k))) THEN
            X(I)=XINIT(K) 
	           XINI(I)=XINIT(K)
          ENDIF
   85   CONTINUE
   90 CONTINUE
      SUMX=0.d0

c     renormalize abundances so that mass fraction sum is 1         
      DO 91 I=1,NSP
         SUMX=SUMX+X(I)
   91 CONTINUE
   
ccc      if((ilast.eq.0).and.(inmode.ne.5))then  
ccc         write(6,*) 'X_i_sum=',sumx
ccc      endif
      
      DO 92 I=1,NSP
         X(I)=X(I)/SUMX           ! Normalisation of Sum X = 1
   92 CONTINUE
c     find max. initial abundance; calculation stops if this
c     abundance is < XLAST
c      isf=1
c      do 84 i=1,nis
c         if(x(i).gt.x(isf)) then
c            isf=i
c         endif
c   84 continue
      
c      if((ilast.eq.0).and.(itest.ne.1).and.(inmode.ne.5)) then      
c         WRITE(9,9978) ON(ISF)
c      endif

c---- interaction links 
      N=0
      AAAAA=' ===>'
      BBBBB='Q MEV'
      DO 93 L=1,NR
         K1(L)=KK1(L)
         K2(L)=KK2(L)
         K3(L)=KK3(L)
         K4(L)=KK4(L)
         K5(L)=KK5(L)
         K6(L)=KK6(L)
         K7(L)=KK7(L)
         K8(L)=KK8(L)
c        store identity of emitted light particle for each reaction;
c        useful for identifying links that give rise to, e.g., neutron
c        emission, etc.

         emitspe(L)=emitspeee(L)

   93 CONTINUE

c     ------------------------------------------------------------------     
c     reverse rates must be multiplied by the same probability factor p;
c     thus find corresponding forward rate for given reverse rate
c     ------------------------------------------------------------------     

ccc      if(ilast.eq.3)then
ccc
ccc         OPEN(8,FILE=networkfile,status='old')
ccc         do 533 j=1,2+nis+1+nsl+1
ccc            read(8,8101) dummystr
ccc 533     continue
ccc         do 534 j=1,nr
ccc            read(8,8101) dummystr
ccc            rspec1(j)=dummystr(10:16)
ccc            rspec2(j)=dummystr(20:26)
ccc            rspec3(j)=dummystr(29:35)
ccc            rspec4(j)=dummystr(39:45)
ccc            
ccc            tspec1(j)=rspec1(j)
ccc            tspec2(j)=rspec2(j)
ccc            tspec3(j)=rspec3(j)
ccc            tspec4(j)=rspec4(j)
ccc 534     continue
ccc         close(8)

c        find corresponding forward and reverse rates
ccc         j=0
ccc 699     j=j+1
ccc         if(char_w(j).eq.'v')then
ccc            kk=0
ccc 698        kk=kk+1
ccc            if((rspec1(j).eq.tspec4(kk)).and.
ccc     &         (rspec2(j).eq.tspec3(kk)).and.
ccc     &         (rspec3(j).eq.tspec2(kk)).and.
ccc     &         (rspec4(j).eq.tspec1(kk)))then
ccc               prob(j)=prob(kk)
ccc               if(j.lt.nr)then 
ccc                  goto 699
ccc               else
ccc                  goto 697
ccc               endif
ccc            else
ccc               if(kk.lt.nr)then 
ccc                  goto 698
ccc               else
ccc                  goto 694
ccc               endif
ccc            endif
ccc 694        write(6,*) ' Forward reaction not found for:'
ccc            write(6,9003) rspec1(j),rspec2(j),rspec3(j),rspec4(j)
ccc            stop
ccc         else
ccc            if(j.lt.nr)then 
ccc               goto 699
ccc            else
ccc               goto 697
ccc            endif
ccc         endif
ccc 697     continue
ccc      endif    

c     ----------------------------------------------------------------
c     write new prob.dat file
c     ----------------------------------------------------------------
ccc      if(ilast.eq.3)then        
ccc         open(37,file='prob.dat',status='unknown')
ccc         do 209 j=1,nr
ccc            write(37,1012) j,prob(j)
c           to account for unknown systematic uncertainty factor, vmult
ccc            vmult(j)=vmult(j)**prob(j)    
ccc 209     continue
ccc         close(37)
ccc      endif
ccc 1012 format(1x,i5,'=',1pe20.12)
c     ----------------------------------------------------------------
c     compute rates
c     ----------------------------------------------------------------
      do 531 j=1,nr
         do 532 i=1,60
ccc            if(ilast.eq.3) then
ccc                rate(j,i)=rate(j,i)*(fu(j,i)**prob(j))
ccc            endif
            ratelog(j,i)=dlog10(rate(j,i))
 532     continue
 531  continue
      
c----------------------------------------------------------------------
c---- print network ---------------------------------------------------
c----------------------------------------------------------------------
      if(itest.eq.1)then
ccc        write(9,*) 'Network taken from ',networkfile
ccc        write(9,*) 'Rates taken from file ',reactionfile
ccc        write(9,*)
        do 696 i=1,nr
ccc          if ((vmult(i).ne.1.d0).and.(ilast.eq.0)) 
ccc     &      write(9,315) i,vmult(i)
 696    continue
 315    format(' ',' Reaction ',i4,' changed by a factor of ',1PD10.3)
ccc        write(9,*)
      endif

      DO 103 I=1,NIS
ccc         IF(ITEST.EQ.1) WRITE(9,9965) I,ON(I),ZN(I),AN(I),X(I)
 103  CONTINUE

      IF(ITEST.EQ.1) THEN
         RHO=2.0d0     ! any value will do because we will divide it out
         DO 115 K=1,NKT
            TAU=TTT(K)*1.D9
            CALL SRATE
c           for EC list decay constant without rho*Ye
            DO 112 N=1,NR
               if(char_ec(N).eq.'ec')then
                  VVX(N,K)=V(N)/(RHO**(TOTPARTIN(N)-1))
                  VVX(N,K)=VVX(N,K)/(RHO*YE)
               elseif(char_w(N).eq.'w')then
                  VVX(N,K)=V(N)
               else
c                  write(6,*) rho
                  VVX(N,K)=V(N)/(RHO**(TOTPARTIN(N)-1))
               endif
  112       CONTINUE
  115    CONTINUE
ccc         write(9,*)
ccc         if(iweak.eq.1)then
ccc            write(9,8899) rho*YE
ccc         endif
ccc         write(9,*)
ccc         WRITE(9,9810) (TTT(K),K=1,NKT)
         DO 118 L=1,NR
ccc            WRITE(9,9071) L,K2(L),ON(K1(L)),K4(L),ON(K3(L)),
ccc     &      AAAAA,K6(L),ON(K5(L)),K8(L),ON(K7(L)),BBBBB,Q(L),
ccc     &      (VVX(L,K),K=1,NKT)
  118    CONTINUE
c----------------------------------------------------------------------
c------- end of run after printing network ----------------------------
c----------------------------------------------------------------------
         stop
      ENDIF

c----------------------------------------------------------------------
 8050 FORMAT(5X,A5,1X,F3.0,1x,F3.0,1x,E30.3)
 8060 FORMAT(1X,I4,3X,I2,1X,A5,2X,I2,1X,A5,1X,I2,1X,A5,
     &2X,I2,1X,A5,3X,F7.3)
 8100 FORMAT(A80)
 8101 FORMAT(A80)
 8899 FORMAT(13x,'STELLAR WEAK RATES ARE CALCULATED FOR rho*Ye=',E10.3)
 9000 format(a70)
 9001 format(a30)
 9002 format(1pe8.2,12(5x,1pe9.3))
 9003 format(2x,i6,2x,a7,' ( ',a7,', ',a7,')  ',a7,4x,1PE10.2) 
 9025 FORMAT(I4,1x,A5,1x,F4.0,F4.0)
 9071 FORMAT(I5,'-',2(1X,I2,1X,A5,' + ',I2,1X,A5,2X,A5),1X,F7.3,
     &7(1PE9.1))
 9111 format(15x,'Network OK, start Run #',i4,' Shell #',i4)
 9810 FORMAT(13X,'REACTIONS AND REACTION RATES NA<SV> AT TEMP(10**9 K):'
     &,1x,7(F9.4),/,153('-'))
 9899 Format(15X,'Network OK, start Run #',i4)
 9965 FORMAT(I4,1X,A5,2F5.0,1PE10.2,12(3X,A5))
c 9978 FORMAT(30X,'NUCLEOSYNTHESIS DURING  ',A5,'  BURNING',/,/)
      RETURN
      END
C======================================================================
      SUBROUTINE STOCK(NT,T,SUMX,X,XST)
C======================================================================
      include 'comnuc.inc'
      dimension X(NSP),XST(KP,NSP)

      NNSTEP(NT)=NSTEP
      TIME(NT)=T
      TEMM(NT)=TAU
      DENN(NT)=RHO
      ERATEY(NT)=EB
      ERATEQ(NT)=EQ
      ENTOTY(NT)=EQQ
      ENTOTQ(NT)=ER
      do 10 i=1,NSP
  10  XST(NT,I)=X(I)

ccc      if((ilast.eq.0).and.(inmode.ne.5)) then
ccc         WRITE(9,9020) ! Intermediate printing of some results
ccc         WRITE(9,9200) T,T/SECY,DELTA,NSTEP
ccc     &     ,TAU,LOG10(TAU),RHO,SUMX,ETA,YE
         
ccc         if(iflux.eq.1) then
ccc   	        WRITE(22,9200) T,T/SECY,DELTA,NSTEP
ccc     &        ,TAU,LOG10(TAU),RHO,SUMX
ccc            WRITE(22,9993)
ccc            DO 1324 I=1,NR
ccc              WRITE(22,1325) I,FLUX(I)
ccc 1324       continue
ccc         endif

ccc         WRITE(9,9995)
ccc      endif

ccc      if((inmode.ne.5).and.(ilast.eq.0))then  
ccc         write(6,9744) t,on(isop1),x(isop1),on(isop2),x(isop2),
ccc     &     sumx,eta,tau
ccc      endif

      SUMENER=0.D0
      DO 433 I=1,NR
  433 SUMENER=SUMENER+ENER(I)

ccc      if((ilast.eq.0).and.(inmode.ne.5)) then
ccc         WRITE(9,2223) SUMENER
ccc         WRITE(9,9997)
ccc         WRITE(9,9998) (ON(I),X(I),I=1,NIS)
ccc         WRITE(9,9020)
ccc      endif

 1325 FORMAT(I5,1X,D22.15)
 2223 FORMAT(1X,1PE9.2)
 9020 FORMAT(153('*'))
 9200 FORMAT(/,1X,'T =',1PE17.10,' S =',E10.3,' YR',2X,
     &'NEW TS =',E10.3,' S',2X,' STEP=',I7,
     &/,1X,'TEMP.=',1PE11.4,' K',2X,'LOG(TEMP)=',0PF7.3,3X,
     &'DENS. =',1PE11.4,' GR CM-3',3X,'SUMX=',0PF16.14,2x,
     &' ETA =',E10.3,2x,' YE =',E10.3)
 9744 FORMAT(1X,'TIME(s)=',1pe8.1,2X,'X(',A5,')=',1pe8.1,
     &2X,'X(',A5,')=',1pe8.1,2x,'SUMX=',0PF16.14,2x,
     &'ETA=',1pe9.2,2x,'T(K)=',1pe8.1)
 9993 FORMAT(/,1X,'REACTION FLUXES  ')
 9994 FORMAT(10(1X,I4,'=',1PE7.1))
 9995 FORMAT(/,1X,'SUM OF REACTION ENERGIES (ERG GR-1 S-1)')
 9996 FORMAT(9(1X,I4,'=',1PE8.1))
 9997 FORMAT(/,1X,'ISOTOPIC ABUNDANCES')
 9998 FORMAT(9(1X,A5,'=',1PE10.2E3))
      RETURN
      END
C======================================================================
      SUBROUTINE SEXIT(NT,T,X,XST)
C======================================================================
      include 'comnuc.inc'
      real*8 deltaprof
      dimension X(NSP),XS(NSP),XST(KP,NSP)

ccc      if(inmode.ne.5) then
ccc         WRITE(9,9270)
ccc      endif

ccc      if((inmode.ne.5).or.(ilast.ne.3))then
ccc         WRITE(9,9271) NSTEP,NLAST,T,TLAST,ON(ISF),X(ISF),XLAST
ccc      endif
         
      if(inmode.eq.2) then
         deltaprof=(tauprofmax-tauprocalmax)/tauprofmax
ccc         if(deltaprof.le.1.d-2) then
ccc            write(9,9007) 'PROFILE CHECK: PASS WITH',deltaprof
ccc          else         
ccc            write(6,9005) 'PROFILE CHECK: WARNING',deltaprof,tauprofmax,
ccc     &          tauprocalmax
ccc            write(9,9005) 'PROFILE CHECK: WARNING',deltaprof,tauprofmax,
ccc     &          tauprocalmax            
ccc         endif
      endif

ccc      if(inmode.eq.4) then
ccc         write(9,*) 'PROFILE CHECK: none'
ccc      endif

      if(inmode.eq.5)then
         trun(ishrun)=t   ! store all running times; they need to be the same for averaging x
         do 777 i=1,nsp
            xshells(i,ishrun)=x(i)
 777     continue
      endif

ccc      if(inmode.eq.5)then
ccc         deltaprof=(tmaxsh(ishrun)-tmaxcode(ishrun))/tmaxsh(ishrun)
ccc         if(deltaprof.ge.1.d-2)then
ccc            write(6,9005) 'PROFILE CHECK: WARNING',deltaprof,
ccc     &          tmaxsh(ishrun),tmaxcode(ishrun)
ccc         endif
ccc      endif      

ccc      if((inmode.eq.5).and.(ilast.eq.0))then
c        ishrun labels the shells for sequential runs if INMODE=5
ccc         write(9,9305) ishrun         
ccc      endif

ccc   ilast=3 in our case and makes first condition TRUE;
ccc   enter if statement
ccc      if((inmode.ne.5).or.(ilast.ne.3))then
ccc         if(inmode.ne.5)then
ccc             WRITE(9,9343) 
ccc           elseif(inmode.eq.5)then
ccc             WRITE(9,9345) 
ccc         endif 
ccc         write(9,9344) TOTENERGY,ETA,YE
ccc         if(inmode.ne.5)then
ccc             WRITE(9,9997)
ccc           elseif(inmode.eq.5)then
ccc             WRITE(9,9991)
ccc         endif 
ccc         WRITE(9,9999) (ON(I),X(I),I=1,NSP)
ccc      endif
      
C     X(I) are mass fractions
      DO 484 I=1,NSP          ! overabundances (X/Xinitial)
         IF(XINI(I).EQ.0.d0) THEN
            XS(I)=0.
          ELSE
            XS(I)=X(I)/XINI(I)
         ENDIF      
C     careful with overabundances:
C      - X(I) are normalized to 1 (see below)
C      - XINI(I) are not normalized            
 484  CONTINUE

ccc      if((inmode.ne.5).and.(ilast.eq.0)) then
ccc         WRITE(9,9998)
ccc         WRITE(9,9999) (ON(I),XS(I),I=1,NSP)
ccc         write(9,9432)
ccc         WRITE(9,9996) (I,tenergy(I),I=1,NR)
ccc         write(9,9221)
ccc         WRITE(9,9220) (I,100*tenergy(I)/TOTENERGY,I=1,NR)

ccc         WRITE(9,9992) (ON(numab(i)),i=1,7)
ccc         DO  486 I=1,NT
ccc           WRITE(9,9988) NNSTEP(I),TIME(I),TEMM(I),DENN(I),
ccc     &       (XST(I,numab(k)),k=1,7)
ccc  486    CONTINUE
ccc         WRITE(9,9270)
ccc      endif

ccc      if(ilast.eq.0)then
ccc         write(6,9376)
ccc      endif
      
ccc      if(inmode.ne.5) then 
ccc         DO 1324 I=1,NR
ccc           if(totflux(i).le.1.d-99) then          
ccc              totflux(i)=0.d0
ccc           endif
ccc           WRITE(21,1325) I, TOTFLUX(I)
ccc 1324    continue
ccc      endif            
            
ccc      close(21)
ccc      close(22)
ccc      close(24)
ccc      close(25)
ccc      close(28)

 1325 FORMAT(I5,1X,D22.15)
 9005 FORMAT(1x,a,1pe8.1,2x,'TauProfile_max (K)=',1PE11.4,2x,
     &  'TauNucleo_max (K)=',1PE11.4)
 9007 FORMAT(1x,a,1pe8.1)
 9020 FORMAT(153('*'))
 9220 FORMAT(9(1X,I5,'=',f4.0))
 9221 format(/,1x,'TOTAL REACTION ENERGIES (%)')
 9270 FORMAT(/,153('$'),/,5('  **** RUN STOP OKAY **** '),/,153('$'),/)
 9271 FORMAT(1X,'NSTEP=',I7,' (NLAST=',I7,')'
     &        ,2X,'TIME=',1pe10.3,' s (TLAST=',1pe10.3,' s)',
     & 2X,'X(',A5,')=',1pe10.3,' (XLAST=',1pe10.3,')')
 9305 Format(1x,'SHELL# ',i4)
 9343 FORMAT(1X,'TOTAL ENERGY (ERG/G)',6x,'ETA',10x,'YE')
 9345 FORMAT(1X,'SHELL ENERGY (ERG/G)',6x,'ETA',10x,'YE')
 9344 FORMAT(3X,1PE10.3,10x,1PE10.3,3x,1PE10.3)
 9376 FORMAT('     **** RUN STOP OKAY **** ')
 9432 format(/,1x,'TOTAL REACTION ENERGIES (erg/g)')
 9987 FORMAT(I4,1X,1PE8.2,4(E8.1),E9.2,9(E8.1))
 9988 FORMAT(1X,I7,1X,1PE8.2,2(E10.3),7(E15.8))
 9992 FORMAT(/,153('-'),/,1X,'NSTEP ','   TIME  ',
     &'  TEMPER  ','  DENSITY ',7('     ',A5,'     '),/,153('-'))
 9996 FORMAT(9(1X,I5,'=',1PE8.1))
 9997 FORMAT(1X,'FINAL ABUNDANCES')
 9991 FORMAT(1X,'ABUNDANCES')
 9998 FORMAT(1X,'FINAL OVERABUNDANCES')
 9999 FORMAT(9(1X,A5,'=',1PE10.2E3))

      return
      END
C======================================================================
      SUBROUTINE SNUC(X,SUMX,DTOLD)
C----------------------------------------------------------------------
C     NUCLEOSYNTHESIS FOR TIME-STEP DELTA
C======================================================================
      include 'comnuc.inc'
      DIMENSION Y0(NSP),YB(NSP),Y(NSP),YT(NSP),X(NSP)
      LOGICAL C
      
      O=1.D-100
      DELT=DELTA
      T=0.d0
   3  DO 5 I=1,NSP       ! Assign Y values...
c        convert X into mole fractions
         Y0(I)=X(I)/AN(I)
         YB(I)=0.d0
         YT(I)=0.d0
   5  CONTINUE
   
      if(inmode.ne.4)then
         call srate
      elseif(inmode.eq.4)then
         call instmix(x)
      endif
      
      C=.FALSE.
      NNNN=0
      NL=0
      PRO=TES

  110 DO 120 I=1,NSP
         Y(I)=Y0(I)              
  120 CONTINUE
C=======================================================================
C     HERE IS THE HARD CORE OF THE PROGRAM: DON'T TOUCH IT!!! 
C     CONSTRUCTION AND SOLUTION OF THE LINEARIZED SYSTEM OF EQUATIONS
C     BY THE 2-STEP METHOD OF WAGONER(1969:AP.J.SUP.18,P.247)
C     CALCULATION OF THE MATRIX ELEMENTS (SUBR MATRIX)
C     SOLUTION OF THE LINEARIZED SYSTEM : A(I,J)*Y(J) = Y(I)(SUBR EIGEN)
C=======================================================================
  130 DO 230 ISTEP=1,2
         IF(C) THEN
              DO 150 I=1,NIS
  150         Y(I)=Y0(I)
         ENDIF
         CALL MATRIX(DELT,Y)
         DO 170 I=1,NIS
  170    YT(I)=Y0(I)
         CALL EIGEN(YT,YB)
         IF(ISTEP.LT.2) THEN
             DO 220 I=1,NIS
                DYT=(YB(I)-Y0(I))/DELT
                IF(ABS(DYT).LT.DYMIN) DYT=DYMIN
                Y(I)=Y0(I)+DYT*DELT
  220        CONTINUE
             C=.FALSE.
         ENDIF
  230 CONTINUE

      continue
      IF(NL.NE.1) PRO1=PRO
      CALL STEST(YB,Y0,PRO,KTEST,DEY)  ! Determine largest DY/Y to affect DT
      AK=PRO1/PRO
      IF(AK-ATEST) 340,340,330

  330 DELT=DELT/PAMIN ! DY/Y too large: results cancelled, Dt=Dtold/pamin..
      NL=1
      NNNN=NNNN+1
      IF(NNNN.GT.50) GO TO 555
      
C     commented out in order to save disk space
c      WRITE(9,9021)
c      WRITE(9,9230) NSTEP,DELT,PRO,AK,ON(KTEST),Y0(KTEST),YB(KTEST),DEY

      C=.TRUE.
      GO TO 110       ! and back to 110 for new calculation

  340 NL=0            ! results acceptable...
      NREAL=NREAL+1

      EQ=0.d0
      EQQ=0.d0
      DO 350 I=1,NR   ! calculate reaction fluxes and energies...
      I1=K1(I)
      I2=K2(I)
      I3=K3(I)
      I4=K4(I)
      I21=I2-1
      I41=I4-1
      A1A=((Y(I1)+O)**I21)*((Y(I3)+O)**I4)*YB(I1)*I2
      A2A=((Y(I3)+O)**I41)*((Y(I1)+O)**I2)*YB(I3)*I4
      FLUX(I)=V(I)*(A1A+A2A)/(F(I2)*F(I4)*(I2+I4))

C     FLUX(I)=V(I) x [rho x NA<sv>]

      AE=Q(I)*FLUX(I)
      AEE=Q(I)*V(I)*(YB(I1)+O)**I2*(YB(I3)+O)**I4/(F(I2)*F(I4))
      ENER(I)=AE*E
      EQ=EQ+AE
      EQQ=EQQ+AEE
  350 CONTINUE
c      IIIII=K1(912)
      IF(NREAL.GE.2) ER=ER+EQ*DELT*E
      EQQ=EQQ*E
      EQ=EQ*E

      SUMX=0.d0
      EB=0.d0
      ETA=0.d0
      DO 360 I=1,NSP
         YDD=Y0(I)
         Y0(I)=YB(I)          ! new Y.....
         X(I)=AN(I)*Y0(I)     ! ... and X values
         ETA=ETA+(AN(I)-2.d0*ZN(I))*Y0(I)
         SUMX=SUMX+X(I)*1.d0
  360 CONTINUE
      YE=0.5d0*(1.d0-ETA)
      EB=E*EB/DELT

C---------   COMPUTATION OF N,P,A EXPOSURES     -----------

      IF(NREAL.EQ.1) SUMX=1.d0
      SE=SUMX-1.d0

C     Check of SumX conservation (Epsilon~1.d-6 OK)
      IF(ABS(SE).GT.EPS)GO TO 430 

      EN=EN+EB*DELT
      T=T+DELT
      DTOLD=DELT
      DELTA=PAR*PRO*DELT
      RETURN
  430 WRITE(9,9260) EPS,SUMX,NSTEP
      STOP
  555 WRITE(9,9998) NNNN
      STOP

 9021 FORMAT(153('='))
 9230 FORMAT(2X,'NS :',I7,4X,'NEW TS =',1PE9.2,' SEC',2X,'PRO =',
     &F10.5,2X,'AK =',D10.3,2X,'RAPID =',2X,A5,2X,'Y0=',1PE8.2,2X,
     &'YB=',1PE8.2,2X,'DEY/Y = ',0PE9.2)
 9260 FORMAT(///,5X,'SUMX ERROR GREATER THAN :',1PE8.2,5X,'SUMX=',
     &0PF20.16,5X,'STEP NO :',I7,///,
     &5X,'CHECK :   1) BARYON NUMBER CONSERVATION IN THE NETWORK',/,
     &5X,'          2) NON-NEGATIVE REACTION RATES IN ROUTINE SRATE',/,
     &5X,'          3) MATRIX INVERSION ROUTINE (EIGEN OR OTHER)',///,
     &2X,6('*- PROGRAM ABORTED -*'))
 9998 FORMAT(///,153('='),///,2X,'NO CONVERGENCE',I7,' ITERATIONS')

      RETURN
      END
C======================================================================
      SUBROUTINE MATRIX(DELT,Y)
C----------------------------------------------------------------------
C     COMPUTES MATRIX ELEMENTS FROM THE NETWORK
C======================================================================
      include 'comnuc.inc'
      DIMENSION Y(NSP)

      NDIP=NDI+1
      AIN=0.D0
      O=1.D-100
C
C     FOR EVERY ISOTOPE I,CALCULATE THE MATRIX LINE(A(I,J),J=I,NISP)
C     FOR EVERY REACTION M : N(II)*II+N(IJ)*IJ---N(IK)*IK+N(IL)*IL
C     THERE IS A CONTRIBUTION :
C      --TO AB(I, I),AB(I, J)  :  IF  I=II  OR  I=IJ
C      --TO AB(I,II),AB(I,IJ)  :  IF  I=IK  OR  I=IL
C      --0.                       IF  I . NE . II,IJ,IK,IL
C
      DO 15 I2 = 1 , NSP
    	 DO 10 I1 = 1 , IJD-1
 10   AV(I1,I2) = 0.D0
 15   CONTINUE

      DO 25 I2 = 1 , IJD-1
    	 DO 20 I1 = 1 , NSP
 20   AH(I1,I2) = 0.D0
 25   CONTINUE

      DO 35 I2 = 1 , NSP
	     DO 30 I1 = 1 , NDI + NDS +1
 30   AD(I1,I2) = 0.D0
 35   CONTINUE

      DO 50 M=1,NR
      II   = K1(M)
      NII  = K2(M)
      IJ   = K3(M)
      NIJ  = K4(M)
      IK   = K5(M)
      NIK  = K6(M)
      IL   = K7(M)
      NIL  = K8(M)

      AMT=V(M) * DELT / (F(NII)*F(NIJ)*(NII+NIJ))
      NI1=NII-1
      NJ1=NIJ-1
      AM=AMT*NII
      A1(M)=AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      A2(M)=AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      AM=AMT*NIJ
      A3(M)=AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      A4(M)=AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      AM=AMT*NIK
      A5(M)=-AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      A6(M)=-AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      AM=AMT*NIL
      A7(M)=-AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      A8(M)=-AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
 50   CONTINUE

      DO 100 M=1,NR
      CALL MATFILL(K1(M),K1(M),A1(M))
      CALL MATFILL(K3(M),K1(M),A2(M))
      CALL MATFILL(K3(M),K3(M),A3(M))
      CALL MATFILL(K1(M),K3(M),A4(M))
      CALL MATFILL(K1(M),K5(M),A5(M))
      CALL MATFILL(K3(M),K5(M),A6(M))
      CALL MATFILL(K1(M),K7(M),A7(M))
      CALL MATFILL(K3(M),K7(M),A8(M))
 100  CONTINUE

      DO 150 I=1,IJD-1
 150  AH(I,I)=AH(I,I)+1.D0
      DO 160 I=IJD,NIS
 160  AD(NDI+1,I)=AD(NDI+1,I)+1.D0

      RETURN
      END
C======================================================================
      SUBROUTINE MATFILL(I,J,AX)
C======================================================================
      include 'comnuc.inc'

      NDIP=NDI+1
      IF(I.EQ.NSP) RETURN
      IF(J.EQ.NSP) RETURN

      IF(J.LT.IJD)THEN
          AH(I,J)=AH(I,J)+AX
      ELSE
          IF(I.LT.IJD)THEN
	         AV(I,J)=AV(I,J)+AX
           ELSE
	         AD(I-J+NDI+1,J)=AD(I-J+NDI+1,J)+AX
           ENDIF
      ENDIF

      RETURN
      END
C======================================================================
      SUBROUTINE EIGEN(C,X)
C======================================================================
      include 'comnuc.inc'
      DIMENSION C(NSP),X(NSP),DIV(NSP),CT(NSP),XSOM(NSP),C0(NSP),
     &          XJAC(NSP),CTT(NSP),RR(NIS)

      NDIP=NDI+1
      N=NIS
      IJDM=IJD-1

      SOMME=0.d0
      DO 6543 I=1,N
      RR(I)=C(I)
      XSOM(I)=C(I)
 6543 CONTINUE

      DO 1 J=1,N
      CTT(J)=C(J)
  1   C0(J)=0.D0

      DO 2 J=1,IJDM
      DO 3 I=1,N
      ATH(I,J)=AH(I,J)
  3   C0(J)=C0(J)+AH(I,J)*C(I)
  2   C0(J)=C(J)-C0(J)

      DO 4 J=IJD,N
      DO 5 I=1,IDEL
      ATV(I,J)=AV(I,J)
  5   C0(J)=C0(J)+AV(I,J)*C(I)
      IMIN=J-NDI
      IMAX=J+NDS
      IF(IMIN.LT.IJD)IMIN=IJD
      IF(IMAX.GT.N)IMAX=N
      JMNDIP=NDIP-J
      DO 6 I=IMIN,IMAX
      ATD(I+JMNDIP,J)=AD(I+JMNDIP,J)
  6   C0(J)=C0(J)+AD(I+JMNDIP,J)*C(I)
  4   C0(J)=C(J)-C0(J)

      DO 7 J=1,N
      C(J)=C0(J)
      CT(J)=C(J)
  7   XSOM(J)=0.D0

C------------ BEGINNING OF LOOP FOR 2 STEP ITERATION ------------------

      DO 10000 ITER=1,2

C---BEGINNING OF THE GAUSSIAN ELIMINATION PROCEDURE
C-(1)-ELIMINATION OF THE LOWER DIAGONALS

      DO 1000 JBAL=IJD,N-1
      DIV(JBAL)=-1.D0/AD(NDIP,JBAL)
      JMAX=JBAL+NDI
      IF(JMAX.GT.N)JMAX=N
      IMAX=JBAL+NDS
      IF(IMAX.GT.N)IMAX=N
      JBNDIP=JBAL+NDIP
      JBMNDIP=NDIP-JBAL
      DO 1000 J=JBAL+1,JMAX
      IF(AD(JBNDIP-J,J).EQ.0.D0)GOTO 1000
      DIVJ=DIV(JBAL)*AD(JBNDIP-J,J)
      DO 10 I=1,IDEL
  10  AV(I,J)=DIVJ*AV(I,JBAL)+AV(I,J)
      JMNDIP=NDIP-J

      DO 20 I=JBAL+1,IMAX
  20  AD(I+JMNDIP,J)=DIVJ*AD(I+JBMNDIP,JBAL)+AD(I+JMNDIP,J)
      C(J)=DIVJ*C(JBAL)+C(J)
 1000 CONTINUE
      DIV(N)=-1.D0/AD(NDIP,N)

C-(2)-ELIMINATION OF THE UPPER DIAGONALS AND OF THE HORIZONTAL BAND

      DO 2000 JBAL=N,IJD+1,-1
      JMIN=JBAL-NDI
      IF(JMIN.LT.IJD)JMIN=IJD
      JBNDIP=JBAL+NDIP
      DO 200 J=JMIN,JBAL-1
      IF(AD(JBNDIP-J,J).EQ.0.D0)GOTO 200
      DIVJ=DIV(JBAL)*AD(JBNDIP-J,J)
      DO 30 I=1,IDEL
 30   AV(I,J)=DIVJ*AV(I,JBAL)+AV(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
 200  CONTINUE
      DO 300 J=1,JDEL
      IF(AH(JBAL,J).EQ.0.D0)GOTO 300
      DIVJ=DIV(JBAL)*AH(JBAL,J)
      DO 40 I=1,IDEL
  40  AH(I,J)=DIVJ*AV(I,JBAL)+AH(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
 300  CONTINUE
 2000 CONTINUE
      DO 400 J=1,JDEL
      IF(AH(IJD,J).EQ.0.D0)GOTO 400
      DIVJ=DIV(IJD)*AH(IJD,J)
      DO 50 I=1,IDEL
 50   AH(I,J)=DIVJ*AV(I,IJD)+AH(I,J)
      C(J)=DIVJ*C(IJD)+C(J)
 400  CONTINUE

C-(3)-GAUSSIAN ELIMINATION OF THE UPPER LEFT SQUARE MATRIX

      DO 3000 JBAL=1,IJD-2
      DIV(JBAL)=-1.D0/AH(JBAL,JBAL)
      DO 3000 J=JBAL+1,IJD-1
      IF(AH(JBAL,J).EQ.0.D0)GOTO 3000
      DIVJ=DIV(JBAL)*AH(JBAL,J)
      DO 60 I=JBAL+1,IJD-1
  60  AH(I,J)=DIVJ*AH(I,JBAL)+AH(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
 3000 CONTINUE
      DIV(IJDM)=-1.D0/AH(IJDM,IJDM)
      X(IJDM)=-DIV(IJDM)*C(IJDM)
      DO 4000 JBAL=IJD-2,1,-1
      SOM=0.D0
      DO 70 I=IJD-1,JBAL+1,-1
  70  SOM=SOM+AH(I,JBAL)*X(I)
 4000 X(JBAL)=DIV(JBAL)*(SOM-C(JBAL))
      SOM=0.D0
      DO 80 I=1,IDEL
  80  SOM=SOM+AV(I,N)*X(I)
      X(N)=DIV(N)*(SOM-C(N))
      DO 5000 JBAL=N-1,IJD,-1
      SOM=0.D0
      DO 90 I=1,IDEL
  90  SOM=SOM+AV(I,JBAL)*X(I)
 5000 X(JBAL)=DIV(JBAL)*(SOM-C(JBAL))

C  MODIFICATION OF XSOM
      DO 12 J=1,N
      XSOM(J)=XSOM(J)+X(J)
 12   C(J)=0.D0

      IF(ITER.EQ.2)GOTO 10000

C RE-INITIALIZATION OF THE MATRIX A (CONSERVED IN AT)
      DO 13 J=1,IJDM
      DO 14 I=1,N
      AH(I,J)=ATH(I,J)
 14   C(J)=C(J)+AH(I,J)*X(I)
 13   C(J)=CT(J)-C(J)
      DO 15 J=IJD,N
      DO 16 I=1,IDEL
      AV(I,J)=ATV(I,J)
 16   C(J)=C(J)+AV(I,J)*X(I)
      IMIN=J-NDI
      IMAX=J+NDS
      IF(IMIN.LT.IJD)IMIN=IJD
      IF(IMAX.GT.N)IMAX=N
      JMNDIP=NDIP-J
      DO 17 I=IMIN,IMAX
      AD(I+JMNDIP,J)=ATD(I+JMNDIP,J)
 17   C(J)=C(J)+AD(I+JMNDIP,J)*X(I)
 15   C(J)=CT(J)-C(J)

10000 CONTINUE

C     ADD THE VALUE OF THE VECTOR "C" THAT WAS INITIALLY SUBSTRACTED
      DO 950 I=1,N
 950  X(I)=XSOM(I)+CTT(I)

      somme=0.d0
      do 7777 i=1,n
         if(x(i).lt.0.d0) x(i)=RR(i)*0.5d0
 7777 somme=somme + x(i)*an(i)

      RETURN
      END
C======================================================================
      SUBROUTINE STEST(YN,YO,PROV,KMM,DDDD)
C-----------------------------------------------------------------------
C     CHOICE OF THE LARGEST ABUNDANCE VARIATION DDDD= DY/Y = (YN-Y0)/Y0
C     OR, EQUIVALENTLY,OF THE SMALLEST PRO=Y(I)/DELY(I).
C     IMPORTANT: ONLY SIGNIFICANT ISOTOPES (Y.GT.YTMIN) CONSIDERED.
C======================================================================
      include 'comnuc.inc'
      DIMENSION YN(NSP),YO(NSP),TEST(NSP)

      KMM=1
      DO 10 I=1,NIS
      DELY=YN(I)-YO(I)
      IF(ABS(DELY).LT.DEMIN) DELY=DEMIN
   4  IF(ABS(YO(I))-YTMIN) 8,7,7
   7  TEST1= ABS(YN(I)/DELY)
      TEST2= ABS(YO(I)/DELY)
      TEST(I)=MIN(TEST1,TEST2)
      GO TO 10
   8  TEST(I)=TES
  10  CONTINUE
      PRO=TEST(1)
      KMM=1
      DO 20 I=1,NIS
      IF(PRO-TEST(I)) 20,20,15
  15  PRO=TEST(I)
      KMM=I
  20  CONTINUE
      PROV=PRO
      YYYY=YO(KMM)
      IF(YYYY.EQ.0.d0) GO TO 999
      DDDD=(YN(KMM)-YYYY)/YYYY
      RETURN
  999 DDDD=0.d0

      RETURN
      END
C=======================================================================
      SUBROUTINE FORMATR
C======================================================================
      include 'comnuc.inc'
      dimension IT(NSP,NSP), JP(NSP)

      kkk = 0
      do 90 i = 1, NSP
         kkk = kkk + 1
         if(kkk .eq. 10) kkk = 0
         JP(i) = kkk
         do 90 j = 1, NSP
   90       IT(j,i) = '.'
      do 50 m = 1, NR
         ii = K1(m)
         ij = K3(m)
         ik = K5(m)
         il = K7(m)
         IT(ii,ii) = '+'
         IT(ij,ii) = '+'
         IT(ij,ij) = '+'
         IT(ii,ij) = '+'
         IT(ii,ik) = '+'
         IT(ij,ik) = '+'
         IT(ii,il) = '+'
         IT(ij,il) = '+'
   50 continue
      do 150 i = 1, NSP
         IT(i,i) = '*'
  150 continue
      kprin = NSP
      nprin = NSP/128
      do 300 kkp = 1, nprin+1
         lk1 = (kkp - 1)*128 + 1
         lk2 = kkp * 128
         if(lk2 .GT. NSP) lk2 = NSP
         if(NSP .GT. 128) kprin = 128
ccc         write(9,9996) (JP(i), i = lk1, lk2)
ccc         do 200 i = lk1, lk2
ccc         write(9,9998) i, (IT(j,i), j = lk1, lk2)
ccc  200    continue
  300 continue

 9996 format(/,4x,128i1,/)
 9998 format(i3,1x,130a1)

      return
      end 
C=========================================================================
      SUBROUTINE SRATE
C-------------------------------------------------------------------------
C     NUCLEAR REACTION RATES V(...)= RHO * NA * <SIGMA*V>  (SEC^-1);
C     reactions in vector V are ordered the same as in REACTIONS.DAT
C======================================================================

      include 'comnuc.inc'
                                                                                
      INTEGER REACTION,TOTIN,TOTOUT,ctemp,k,crho 
      double precision t9log,extraplograte,rate2log(60),rhoye2log
      double precision deltatau,deltarho,deltaye

C  calculate rates only if tau, rho or Ye change significantly      
C  -----------------------------------------------------------
      
      deltatau=abs(tau-oldtau)/tau
      deltarho=abs(rho-oldrho)/rho
      if(iweak.eq.1)then
         deltaye=abs(ye-oldye)/ye
      else
         deltaye=0.d0
      endif
       
      if((deltatau.gt.1.d-5).or.(deltarho.gt.1.d-5).or.
     &    (deltaye.gt.1.d-5))then      

         t9log=dlog10(tau*1.d-09)  
C        find where we are on T grid; T is between grid points ctemp and ctemp+1
         call locate(templog,60,t9log,ctemp)
         if(ctemp.eq.60)then
            write(6,*) 'Temperature > 10 GK! Run stop!'
            write(6,*) 'TAU=',tau*1.d-09,' GK'
            stop
         endif

         do 16 REACTION=1,nr  

c---------- see how many particles we have in the entrance/exit channels
            TOTIN=K2(REACTION)+K4(REACTION)
            totpartin(REACTION)=TOTIN
            TOTOUT=K6(REACTION)+K8(REACTION)

            IF(ctemp.ne.0)THEN			! otherwise T is below 1 MK
            
c---------- convert 2D into 1D array of table rates; only 2 temperature grid points needed
            do 17 i=ctemp,ctemp+1
               rate2log(i)=ratelog(REACTION,i)    
   17       continue 
         
c           for iweak=1 and char_w(REACTION).eq.'w', replace rate2log with the one for the
c           proper rho*Ye
            if((iweak.eq.1).and.(char_w(REACTION).eq.'w'))then
             
               rhoye2log=dlog10(rho*YE)               
               do 20 k=1,nrw 
                  if(reacstring(REACTION)(6:30).eq.
     &               reacstrweak(k)(6:30))then
c                    find where we are on rho grid  
                   
                     call locate(rhoyelog,12,rhoye2log,crho)
                     if(crho.eq.12)then
                        write(6,*) 'rho*Ye > 10**11! Run stop!'
                        stop
                     endif

                     IF(crho.ne.0)THEN

                     rate2log(ctemp)=
     &                 (ratelogw(k,ctemp,crho+1)-
     &                 ratelogw(k,ctemp,crho))*
     &                 (rhoye2log-rhoyelog(crho))/
     &                 (rhoyelog(crho+1)-rhoyelog(crho))+
     &                 ratelogw(k,ctemp,crho)
     
                     rate2log(ctemp+1)=
     &                 (ratelogw(k,ctemp+1,crho+1)-
     &                 ratelogw(k,ctemp+1,crho))*
     &                 (rhoye2log-rhoyelog(crho))/
     &                 (rhoyelog(crho+1)-rhoyelog(crho))+
     &                 ratelogw(k,ctemp+1,crho)
                  
c                    include probability factor; factor uncertainty=2 is hardwired!
                     rate2log(ctemp)=
     &                 dlog10((10.d0**rate2log(ctemp))*
     &                 (2.d0**prob(REACTION)))
                     rate2log(ctemp+1)=
     &                 dlog10((10.d0**rate2log(ctemp+1))*
     &                 (2.d0**prob(REACTION)))

                     ENDIF
                     
                     goto 21
                  endif
   20          continue                                                                                
            endif
 21         continue
         
            ENDIF 
c--------------------------------------------------------------------            
c---------- find interpolated log rate ----------------------------->

            if((char_ec(REACTION).ne.'ec').and.
     &           (char_w(REACTION).ne.'w'))then
               if(tau.ge.1.d6)then 
                  extraplograte=(rate2log(ctemp+1)-
     &              rate2log(ctemp))*
     &              (t9log-templog(ctemp))/
     &              (templog(ctemp+1)-templog(ctemp))
     &              +rate2log(ctemp)
                  V(REACTION)=vmult(reaction)*
     &              RHO**(TOTIN-1)*
     &              (10.d0**extraplograte)
               else
                  V(REACTION)=0.d0   ! rate zero for T < 1 MK 
               endif

            elseif(char_ec(REACTION).eq.'ec')then
               if(tau.ge.1.d6)then
                  extraplograte=(rate2log(ctemp+1)-
     &              rate2log(ctemp))*
     &              (t9log-templog(ctemp))/
     &              (templog(ctemp+1)-templog(ctemp))
     &              +rate2log(ctemp)
c-------------    multiply all links labeled 'ec' by rho*Ye=rho*(1-eta)/2; 
c                 this only applies to the 2 rates: pep and 3He decay,
c                 both adopted from CF88; none of these are contained
c                 in stellar weak interaction file!!!                      
                  V(REACTION)=vmult(reaction)*RHO**(TOTIN-1)*
     &             (10.d0**extraplograte)*RHO*YE
               else
                  V(REACTION)=0.d0
               endif

            elseif(char_w(REACTION).eq.'w')then
               if(iweak.eq.0)then
c                 if iweak=0: chose lab rates [here chosen as lab rate at
c                 10 MK; doesn't matter which T value is picked, because
c                 rate values are the same at all T]
                  V(REACTION)=vmult(reaction)*rate(reaction,10)
               elseif(iweak.eq.1)then
                  if((ctemp.ne.0).and.(crho.ne.0))then

                     extraplograte=(rate2log(ctemp+1)-
     &                rate2log(ctemp))*
     &                (t9log-templog(ctemp))/
     &                (templog(ctemp+1)-templog(ctemp))
     &                +rate2log(ctemp)
                     V(REACTION)=vmult(reaction)*
     &                (10.d0**extraplograte)
                  else      
c                    for T < 1 MK or rho*Ye < 1 g/cm3 use lab rate 
                     V(REACTION)=vmult(reaction)*rate(reaction,10)
                  endif
               endif
            endif
c<--------- find interpolated log rate ------------------------------
c--------------------------------------------------------------------
            
            if(V(REACTION).lt.0.0) then
               write(6,*) 'WARNING: RATE IS NEGATIVE!'
               stop
            endif
c           -------------------------------------------------
c           this part allows user to set blocks of reaction
c           rates to zero; useful to find the reaction that
c           gives rise to a certain effect, for example, neutron
c           production, etc.
c           -------------------------------------------------
c            if((reaction.ge.221).and.(reaction.le.225))then
c               if(emitspe(reaction).eq.'NE  1')then
c                  v(reaction)=0.d0
c               endif
c            endif
c           -------------------------------------------------
   16    continue      
         oldtau=tau
         oldrho=rho
         oldye=ye

      endif

      RETURN
      END
C===========================================================================
      SUBROUTINE PROFIL(T)
C---------------------------------------------------------------------------
C     finds temperature (tau) and density (rho) for given time (t) 
C     value;
C     for INMODE=2: T, rho from external profile via interpolation in t
C     for INMODE=3: T, rho from analytical expression
C     for INMODE=5: T, rho from external shells input via interpolation in t
C===========================================================================
      include 'comnuc.inc'

c---- interpolation for external profile
      if(inmode.eq.2) then      
         i=1

         call locate(tprof,numod,t,i) ! t between grid points i and i+1
         
         tau=((t-tprof(i))/(tprof(i+1)-tprof(i)))*(tauprof(i+1)-
     &          tauprof(i))+tauprof(i)
         rho=((t-tprof(i))/(tprof(i+1)-tprof(i)))*(rhoprof(i+1)-
     &          rhoprof(i))+rhoprof(i)
         if(tau.ge.tauprocalmax) then
            tauprocalmax=tau
         endif
	 
c        monitor profile
cc         if(ilast.eq.0) then
ccc            write(26,5068) t,tau,rho
c            write(26,5065) nstep,t,tau,rho
cc		 endif
      endif

c---- calculation from analytical function     
      if(inmode.eq.3) then
         tau=tau0*dexp(-t/scale1)
         rho=rho0*dexp(-t/scale2)
         if(ilast.eq.0) then	 
ccc	           write(26,5065) nstep,t,tau,rho
	        endif
      endif

c---- interpolation for sequential runs over shells
      if(inmode.eq.5) then    
         i=1      ! i denotes model
         j=ishrun ! j denotes shell

         call locate(tmodel,numod,t,i) ! t between grid points i and i+1

         tau=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (taushell(i+1,j)-taushell(i,j))+taushell(i,j)
         rho=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (rhoshell(i+1,j)-rhoshell(i,j))+rhoshell(i,j)

         if(tau.ge.tmaxcode(ishrun)) then
            tmaxcode(ishrun)=tau
         endif
	 
c        monitor profile
ccc         write(26,5066) ishrun,nstep,t,tau,rho
         
      endif

 5065 format(1x,i7,2x,1pe22.15,e10.3,e10.3)
 5068 format(1x,1pe22.14e3,e10.3,e10.3)
 5066 format(1x,'shell#',i7,2x,i7,2x,1pe22.15,e10.3,e10.3)
            
      RETURN
      END
C========================================================================
      SUBROUTINE NEWSHELLS(dconv)     
C------------------------------------------------------------------------
C     routine checks if each model has the same number of shells [only 
C     for INMODE=4,5]; if not, shells need to be re-sized using the model 
C     with the largest number of shells as template; a test is also 
C     performed to check if the total mass in each model is approximately
C     constant; if not, the run stops
C========================================================================
      include 'comnuc.inc'
      character*1 cdiffer
      real*8 totmash(numod)    
      real*8 mgridhelp,mgrid(0:nshmax),msh0(nmod,0:nshmax)
      real*8 mhelp,help1,help2,help3,dconv
      real*8 cumshell(0:nshmax)
      integer ndep
      
      open(16,file='newgrid.out',status='unknown')

c---- i denotes model (at given time), j denotes shells in given model        
c---- numod: number of models 
c---- nshmax: maximum number of shells
c---- modmaxnum: model number with max. number of shells
c---- totmash: total mass of each model
c---- mgrid: cumulative mass for template used to re-size other shells
c---- msh0: cumulative mass for all models
c---- cumshell: cumulative mass distribution

c---- test: do we have the same number of shells in each model? 
      cdiffer='y'
      do 200 i=1,numod
         if(nsh0(i).ne.nshmax)then
            cdiffer='n'
         endif
 200  continue

      if(cdiffer.eq.'y')then
c     ...yes; shells do not need re-sizing
         do 300 i=1,numod
            do 400 j=1,nshmax
               taushell(i,j)=taush0(i,j)
               rhoshell(i,j)=rhosh0(i,j)
               vconshell(i,j)=vconsh0(i,j)
               clenshell(i,j)=clensh0(i,j)
               dmshell(i,j)=dmsh0(i,j)
 400        continue
 300     continue
       elseif(cdiffer.eq.'n')then      
c      ...no; shells need re-sizing

c----    first, find total mass for each model
         do 20 i=1,numod
            totmash(i)=0.d0
            do 10 j=1,nsh0(i)
               totmash(i)=totmash(i)+dmsh0(i,j)
 10         continue
 20      continue

c----    is total mass in each model the same (within 0.1%, say)?          
         do 30 i=1,numod
           if(dabs((totmash(modmaxnum)-totmash(i))/totmash(modmaxnum))
     &         .ge.0.001) then
              write(6,*) 'Masses in models are not the same! Run stop!'
              stop
           endif
 30      continue

c----    model with maximum number of shells serves as template grid; 
c----    establish a cumulative mass coordinate for template 
         mgridhelp=0.d0
         mgrid(0)=0.d0
         do 600 j=1,nshmax
            mgridhelp=mgridhelp+dmsh0(modmaxnum,j)
            mgrid(j)=mgridhelp
 600     continue
c----    mgrid(nshmax) is the total mass of template grid
         write(16,*) 'Total mass of template model=',mgrid(nshmax)

c----    find cumulative mass coordinates for all models
         do 900 i=1,numod
            mhelp=0.d0
            msh0(i,0)=0.d0
            do 800 j=1,nsh0(i)
               mhelp=mhelp+dmsh0(i,j)
               msh0(i,j)=mhelp
 800        continue

c----       set total mass of model equal to total template grid mass 
            msh0(i,j)=mgrid(nshmax)
 900     continue
    
c----    apply template grid to all models; find T, rho etc. 

         do 1000 i=1,numod          
            do 2000 m=0,nshmax-1

               k=-1
 2500          k=k+1
               if(msh0(i,k).gt.mgrid(m))then
                  kmin=k
                else
                  goto 2500 
               endif
            
               k=-1
 2600          k=k+1
               if(msh0(i,k).ge.mgrid(m+1))then
                  kmax=k
                else
                  goto 2600 
               endif
            
               if(kmin.eq.kmax)then
                  dmshell(i,m+1)=mgrid(m+1)-mgrid(m)
                  taushell(i,m+1)=taush0(i,kmax)
                  rhoshell(i,m+1)=rhosh0(i,kmax)
                  vconshell(i,m+1)=vconsh0(i,kmax)
                  clenshell(i,m+1)=clensh0(i,kmax)
                else
                  
                  dmshell(i,m+1)=mgrid(m+1)-mgrid(m)
                  
c----             temperature
                  help1=(msh0(i,kmin)-mgrid(m))*taush0(i,kmin)
                  help2=(mgrid(m+1)-msh0(i,kmax-1))*taush0(i,kmax)  
                  help3=0.d0
                  do 999 k=kmin+1,kmax-1
                     help3=help3+dmsh0(i,k)*taush0(i,k)
 999              continue
                  taushell(i,m+1)=(help1+help2+help3)/dmshell(i,m+1)
                                    
c----             density
                  help1=(msh0(i,kmin)-mgrid(m))*rhosh0(i,kmin)
                  help2=(mgrid(m+1)-msh0(i,kmax-1))*rhosh0(i,kmax)  
                  help3=0.d0
                  do 998 k=kmin+1,kmax-1
                     help3=help3+dmsh0(i,k)*rhosh0(i,k)
 998              continue
                  rhoshell(i,m+1)=(help1+help2+help3)/dmshell(i,m+1)

c----             convection velocity
                  help1=(msh0(i,kmin)-mgrid(m))*vconsh0(i,kmin)
                  help2=(mgrid(m+1)-msh0(i,kmax-1))*vconsh0(i,kmax)  
                  help3=0.d0
                  do 997 k=kmin+1,kmax-1
                     help3=help3+dmsh0(i,k)*vconsh0(i,k)
 997              continue
                  vconshell(i,m+1)=(help1+help2+help3)/dmshell(i,m+1)

c----             mixing length
                  help1=(msh0(i,kmin)-mgrid(m))*clensh0(i,kmin)
                  help2=(mgrid(m+1)-msh0(i,kmax-1))*clensh0(i,kmax)  
                  help3=0.d0
                  do 996 k=kmin+1,kmax-1
                     help3=help3+dmsh0(i,k)*clensh0(i,k)
 996              continue
                  clenshell(i,m+1)=(help1+help2+help3)/dmshell(i,m+1)                  
               endif

 2000       continue
 1000    continue
      
c----    output of re-sized shells
         write(16,*) 'Total number of models [time steps]=',numod
         write(16,*) 'Original mass shells for each model, temperature'
         write(16,*)         
         do 111 i=1,numod
            write(16,9001) i,tmodel(i)
            do 222 j=1,nsh0(i)
               write(16,9000) dmsh0(i,j),taush0(i,j)
 222        continue
            write(16,*)
 111     continue
         write(16,*)
         write(16,*)
         write(16,*) 'Re-sized mass shells for each model, temperature' 
         write(16,*)        
         do 444 i=1,numod
            write(16,9001) i,tmodel(i)
            do 333 j=1,nshmax
               write(16,9000) dmshell(i,j),taushell(i,j)
 333        continue
            write(16,*)
 444     continue

      endif
                
c     if only outer shells burn, up to a certain depth, and inner shells burn
c     but are not mixed to surface: disregard all inner shells [INMODE=5 only]   
      if(dconv.ne.0.d0)then
         cumshell(0)=0.d0
c        find cumulative shell mass distribution from first model
         do 251 j=1,nshmax  
            cumshell(j)=dmshell(1,j)+cumshell(j-1)
            write(16,*)
            write(16,*) 'Cumulative mass of shells'
            write(16,*) j,cumshell(j)
 251     continue

         call locate(cumshell,nshmax,dconv,ndep)
c        depth of convection is located between grid points ndep and ndep+1
         
         if((ndep.eq.0).or.(ndep.eq.nshmax))then
            write(6,*) 'dconv out of range! Run stop!'
            stop
         endif
         
         write(16,*)
         write(16,*) 'dconv=',dconv    ! set in nucleo.in
         write(16,*) 'ndep=',ndep	   ! found by code
         
c        take burning of all shells with *** index .gt. ndep *** into account   
c        [excluding ndep; from ndep+1 to nshmax]      
                                  
         do 302 i=1,numod
            do 402 j=ndep+1,nshmax
               taushell(i,j-ndep)=taushell(i,j)
               rhoshell(i,j-ndep)=rhoshell(i,j)
               vconshell(i,j-ndep)=vconshell(i,j)
               clenshell(i,j-ndep)=clenshell(i,j)
               dmshell(i,j-ndep)=dmshell(i,j)
 402        continue
 302     continue

         nshmax=nshmax-ndep     ! new maximum number of shells

c        set values of all other [inner] shells equal to zero
         do 303 i=1,numod
            do 403 j=nshmax+1,nshmax+ndep
               taushell(i,j)=0.d0
               rhoshell(i,j)=0.d0
               vconshell(i,j)=0.d0
               clenshell(i,j)=0.d0
               dmshell(i,j)=0.d0
 403        continue
 303     continue
        
      endif
      
c     find maximum temperature for each shell

      do 456 j=1,nshmax
         tmaxsh(j)=0.d0
         tmaxcode(j)=0.d0
 456  continue

      do 454 i=1,numod
         do 453 j=1,nshmax
            if(taushell(i,j).ge.tmaxsh(j))then
               tmaxsh(j)=taushell(i,j)
            endif
 453     continue
 454  continue    
      
      close(16)

 9000 format(2x,1PE13.6,3x,1PE13.6)
 9001 format(1x,'model=',i5,5x,'time=',1PE13.6)
 
      RETURN
      END
C=========================================================================
      SUBROUTINE INSTMIX(x)   
C-----------------------------------------------------------------------
C     averages rates over mass shells [INMODE=4]; info on temperature, 
C     density, mass of shells is given in file shellinputfile
C-----------------------------------------------------------------------
C     NUCLEAR REACTION RATES V(...)= RHO * NA<SIGMA*V>  (SEC-1);
C     reactions in vector V are ordered the same as in REACTIONS.DAT
C======================================================================

      include 'comnuc.inc'
	  
	  DIMENSION X(NSP)
      INTEGER REACTION,TOTIN,totout,ctemp,k,crho  
      double precision t9logshell,rate2log(60),rhoye2log
      double precision ysecond(60),extraplograte
	  double precision vinspect,timeconvshell

      do 1001 REACTION=1,nr

         v(reaction)=0.d0

c------- see how many particles we have in the entrance/exit channels
         TOTIN=K2(REACTION)+K4(REACTION)
         TOTOUT=K6(REACTION)+K8(REACTION)
         
c------- j labels shells
         do 20 j=1,nshmax
            t9logshell=dlog10(taus(j)*1.d-09)

c           find two T grid points
            call locate(templog,60,t9logshell,ctemp)            
            if(ctemp.eq.60)then
               write(6,*) 'Temperature > 10 GK! Run stop!'
               stop
            endif

            IF(ctemp.ne.0)THEN			! otherwise T is below 1 MK
            
c           convert 2D into 1D array of table rates; only two T values are needed
            do 17 i=ctemp,ctemp+1
               rate2log(i)=ratelog(REACTION,i)              
  17        continue 

            timeconvshell=clens(j)/vcons(j)

c           for iweak=1: replace rate2log with the one for the proper rho*Ye
            if((iweak.eq.1).and.(char_w(REACTION).eq.'w'))then

               rhoye2log=dlog10(rhos(j)*YE)               
               do 30 k=1,nrw            
                  if(reacstring(REACTION)(6:30)
     &               .eq.reacstrweak(k)(6:30))then
c                 find where we are on rho grid                     
                     call locate(rhoyelog,12,rhoye2log,crho)
                     if(crho.eq.12)then
                        write(6,*) 'rho*Ye > 10**11! Run stop!'
                        stop
                     endif

                     IF(crho.ne.0)THEN

                     rate2log(ctemp)=
     &                 (ratelogw(k,ctemp,crho+1)-
     &                 ratelogw(k,ctemp,crho))*
     &                 (rhoye2log-rhoyelog(crho))/
     &                 (rhoyelog(crho+1)-rhoyelog(crho))+
     &                 ratelogw(k,ctemp,crho)
     
                     rate2log(ctemp+1)=
     &                 (ratelogw(k,ctemp+1,crho+1)-
     &                 ratelogw(k,ctemp+1,crho))*
     &                 (rhoye2log-rhoyelog(crho))/
     &                 (rhoyelog(crho+1)-rhoyelog(crho))+
     &                 ratelogw(k,ctemp+1,crho)
                  
c                    include probability factor; factor uncertainty=2 is hardwired!
                     rate2log(ctemp)=
     &                 dlog10((10.d0**rate2log(ctemp))*
     &                 (2.d0**prob(REACTION)))
                     rate2log(ctemp+1)=
     &                 dlog10((10.d0**rate2log(ctemp+1))*
     &                 (2.d0**prob(REACTION)))

                     ENDIF

                     goto 31
                  endif
 30            continue                                                                    
            endif
 31         continue

            ENDIF

c--------------------------------------------------------------------            
c---------- find interpolated log rate ------------------------------

            if((char_ec(REACTION).ne.'ec').and.
     &           (char_w(REACTION).ne.'w'))then
               if(taus(j).ge.1.d6)then 
                  extraplograte=(rate2log(ctemp+1)-rate2log(ctemp))*
     &              (t9logshell-templog(ctemp))/
     &              (templog(ctemp+1)-templog(ctemp))
     &              +rate2log(ctemp)
                  vinspect=vmult(reaction)*(rhos(j)**
     &              (TOTIN-1))*(10.d0**extraplograte)
               else
                  vinspect=0.d0               
               endif
               V(REACTION)=V(REACTION)+vinspect*dms(j)

            elseif(char_ec(REACTION).eq.'ec')then
               if(taus(j).ge.1.d6)then
                  extraplograte=(rate2log(ctemp+1)-rate2log(ctemp))*
     &              (t9logshell-templog(ctemp))/
     &              (templog(ctemp+1)-templog(ctemp))
     &              +rate2log(ctemp)
c-------------    multiply all links labeled 'ec' by rho*Ye=rho*(1-eta)/2; 
c                 this only applies to the 3 rates: pep, 3He decay, and 7Be
c                 decay, all adopted from CF88; none of these are contained
c                 in stellar weak interaction file!!!                      
                  vinspect=vmult(reaction)*rhos(j)**(TOTIN-1)*
     &             (10.d0**extraplograte)*rhos(j)*YE
               else
                  vinspect=0.d0
               endif
               V(REACTION)=V(REACTION)+vinspect*dms(j)

            elseif(char_w(REACTION).eq.'w')then
               if(iweak.eq.0)then
                  vinspect=vmult(reaction)*rate(reaction,10)
               elseif(iweak.eq.1)then
                  if((ctemp.ne.0).and.(crho.ne.0))then
                     extraplograte=(rate2log(ctemp+1)-rate2log(ctemp))*
     &                 (t9logshell-templog(ctemp))/
     &                 (templog(ctemp+1)-templog(ctemp))
     &                 +rate2log(ctemp)
                     vinspect=vmult(reaction)*
     &                 (10.d0**extraplograte)
                  else
c                    for T < 1 MK or rho*Ye < 1 g/cm3 use lab rate
                     vinspect=vmult(reaction)*rate(reaction,10)
                  endif
               endif
               V(REACTION)=V(REACTION)+vinspect*dms(j)
            endif
c---------- find interpolated log rate ------------------------------
c--------------------------------------------------------------------
c---------- use the following if you are interested in this output:            
cc---------- output lifetimes if 1/tau_nucl > 1/tau_conv
c            if((vinspect*X(K3(reaction))/AN(K3(reaction))).ge.
c     &         1.d0/timeconvshell) then
c               write(28,9071) reaction,K2(reaction),
c     &            ON(K1(reaction)),K4(reaction),
c     &            ON(K3(reaction)),K6(reaction),ON(K5(reaction)), 
c     &            K8(reaction),ON(K7(reaction)),
c     &            (vinspect*X(K3(reaction)))/AN(K3(reaction)),
c     &            1.d0/timeconvshell
c             elseif((vinspect*X(K1(reaction))/AN(K1(reaction))).ge.
c     &            1.d0/timeconvshell) then
c			   WRITE(28,9072) reaction,K2(reaction),
c     &            ON(K1(reaction)),K4(reaction),
c     &            ON(K3(reaction)),K6(reaction),ON(K5(reaction)), 
c     &            K8(reaction),ON(K7(reaction)),
c     &            (vinspect*X(K1(reaction)))/AN(K1(reaction)),
c     &            1.d0/timeconvshell
c            endif 
c--------------------------------------------------------------------
            
   20    continue     ! loop over shells

         v(reaction)=v(reaction)/dmst

         if(V(REACTION).lt.0.d0)then
            write(6,*) 'WARNING: RATE IS NEGATIVE!'
            stop
         endif

 1001 continue        ! loop over reactions    

 9071 FORMAT(I5,2X,I2,1x,A5,' + ',2X,I2,1x,A5,2x,'===>',
     &   2x,I2,1x,A5,' + ',2X,I2,1x,A5,4x,'target:    ',1PE9.1,3x,
     &   'conv:',1PE9.1)
 9072 FORMAT(I5,2X,I2,1x,A5,' + ',2X,I2,1x,A5,2x,'===>',
     &   2x,I2,1x,A5,' + ',2X,I2,1x,A5,4x,'projectile:',1PE9.1,3x,
     &   'conv:',1PE9.1)

      return
      end
C=======================================================================
      SUBROUTINE CLEARSHELL(X)
C-----------------------------------------------------------------------
C     initialization after each run for sequential runs over shells;
C     for INMODE=5
C======================================================================
      include 'comnuc.inc'

      dimension X(NSP)

      CALL RPARAM      

c     set mass fractions to normalized initial abundances
      do 344 i=1,nsp
         x(i)=xini(i)
 344  continue
      
C---------------- INITIALIZATION -------------------------------
      NSTEP=0
      NREAL=0
      njacb=0   ! Jacobian build and inversions
      njace=0   ! Jacobian evaluations
      nfunc=0   ! RHS function evaluations
      nfail=0   ! Failed steps
      ncor=0    ! Newtonian corrector steps
      nmakejac=jacrate ! only invert jacobian every jacrate steps
      EN=0.d0
      ER=0.d0
      oldtau=0.d0
      oldrho=0.d0
      oldye=0.d0

      ncurrtime=1

      TOTENERGY=0.D0
      DO 1264 I=1,NR
 1264 TOTFLUX(I)=0.D0

C     Zero abundance errors
      do i=1,nis
         TotErr(i)=0.d0
      enddo

      CALL PROFIL(T)   ! calculate initial TAU, RHO from profile

ccc      if(ilast.eq.0)then
ccc         write(9,9300)
ccc      endif
      
 9300 format(/,153('$'),/)	
        
      RETURN
      END
C======================================================================
      SUBROUTINE MEANX
C----------------------------------------------------------------------
C     for sequential runs over shells, the final mass fractions are 
C     weighted at the end; need to make sure that nucleosynthesis for
C     every shell was stopped at the same time [for INMODE=5]
C======================================================================
      include 'comnuc.inc'
      dimension xav(nsp),xs(nsp)
      real*8 xavsum,totenerav
      character*1 cdiff

c     check first if all running times are the same before averaging x
      cdiff='n'
      do 200 i=1,nshmax
c         print*,i,trun(i),tlast
         if(trun(i).lt.tlast) then
            cdiff='y'
         endif
 200  continue

c     each shell must run until TLAST, otherwise stop run
      if(cdiff.eq.'y') then
         write(6,*) 'SHELLS DID NOT RUN UNTIL TLAST! CHANGE TLAST 
     *     IN INPUT!'
         stop
      endif
      
c     get masses of all shells at TLAST: dms(j)
      call searchmodel(tlast)
      
c     weight mass fractions by mass of shell; i: isotopes; j: shells      
      do 400 i=1,nsp
         xav(i)=0.d0
         do 300 j=1,nshmax
            xshells(i,j)=xshells(i,j)*dms(j)/dmst
            xav(i)=xav(i)+xshells(i,j)
 300     continue
 400  continue

c     check mass fraction sum
      xavsum=0.d0
      do 500 i=1,nsp
         xavsum=xavsum+xav(i)
 500  continue
  
      totenerav=0.d0
c     find total energy  
      do 700 j=1,nshmax
         totenerav=totenerav+tenersh(j)*dms(j)/dmst
 700  continue
 
c     find eta
      ETA=0.d0
      DO 800 I=1,NSP
         ETA=ETA+(AN(I)-2.d0*ZN(I))*(xav(I)/AN(I))
 800  CONTINUE
      YE=0.5d0*(1.d0-ETA)
 
c     output results
ccc      write(9,9300)

ccc      WRITE(9,9271) tlast,xavsum,dmst
ccc      write(9,9272) 
ccc      write(9,9273) totenerav,eta,ye

ccc      WRITE(9,9997)
ccc      WRITE(9,9999) (ON(I),XAV(I),I=1,NSP)

      if((abs(xavsum-1.d0)).gt.0.001) then
         write(6,*) 'WARNING: AVERAGED MASS FRACTION SUM IS NOT 1!'
         stop
      endif

      DO 600 I=1,NSP          ! overabundances (X/Xinitial)
         IF(XINI(I).EQ.0.d0) THEN
            XS(I)=0.d0
          ELSE
            XS(I)=XAV(I)/XINI(I)
         ENDIF
c        careful with overabundances:
c        - X(I) are normalized to 1 (see below)
c        - XINI(I) are not normalized            
 600  CONTINUE

ccc      if(ilast.eq.0)then
ccc         write(9,*)
ccc         WRITE(9,9998)
ccc         WRITE(9,9999) (ON(I),XS(I),I=1,NSP)
ccc      endif
      
 9271 FORMAT(1x,'TLAST=',1pe10.3,' s','  AVER. MASS FRACTION SUM:',
     &  1pe14.7,'  MASS OF ALL SHELLS:',
     &  1pe11.4)
 9272 format(1X,'TOTAL ENERGY (ERG/G)',6x,'ETA',10x,'YE')
 9273 format(3X,1PE10.3,10x,E10.3,3x,E10.3)           
 9300 FORMAT(/,153('$'),/,5('  **** RUN STOP OKAY **** '),/,153('$'),/)
 9997 FORMAT(1X,'FINAL AVERAGED ABUNDANCES')
 9998 FORMAT(1X,'FINAL OVERABUNDANCES')
 9999 FORMAT(9(1X,A5,'=',1PE10.2E3)) 
 
      return
      END
C=========================================================================
      SUBROUTINE SEARCHMODEL(t)     
C-------------------------------------------------------------------------
C     uses time to interpolate between models when rates are averaged over 
C     shells [INMODE=4]
C=========================================================================
      include 'comnuc.inc'

      dmst=0.d0

      i=1

      call locate(tmodel,numod,t,i) ! t between grid points i and i+1
      if((i.eq.0).or.(i.eq.numod))then
         write(6,*) 'searchmodel out of range! Run stop!'
         stop
      endif

c---- interpolation
      do 300 j=1,nshmax
            taus(j)=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (taushell(i+1,j)-taushell(i,j))+taushell(i,j)
            rhos(j)=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (rhoshell(i+1,j)-rhoshell(i,j))+rhoshell(i,j)
            vcons(j)=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (vconshell(i+1,j)-vconshell(i,j))+vconshell(i,j)
            clens(j)=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (clenshell(i+1,j)-clenshell(i,j))+clenshell(i,j)
            dms(j)=((t-tmodel(i))/(tmodel(i+1)-tmodel(i)))*
     &        (dmshell(i+1,j)-dmshell(i,j))+dmshell(i,j)
c           write(6,*) taus(j),rhos(j)
            dmst=dmst+dms(j)    ! total mass of shells in interpolated model
 300  continue

      RETURN
      END
C======================================================================
      SUBROUTINE steplimiter(eata,t,delt)
C======================================================================
C     Subroutine to limit step sizes
C======================================================================

      include 'comgear.inc'

      PARAMETER(dtaurec  =2.0d7)

      double precision newtime,nexttime,dtau,otau
      integer idmax
      logical hitgrids,tcheck,extralimit
      hitgrids = .false.
      tcheck = .true.
      extralimit = .true.

      if(DEBUG .eq. 1)
     &     print*,"start",eata,eata*delt


      if(extralimit)then
C        Limit step size changes:
C        Step 1-10: eta <= 100
C        Step 10-end: eta <= 10
         if(nstep .lt. 5)then
            eata = dmin1(eata,100.d0)
         else if(nstep .lt. 10) then
            eata = dmin1(eata,10.d0)
            eata = dmax1(eata,0.1d0)
         else
c           if(delt .lt. 1.0d-5)then
c           eata = dmin1(eata,1.d2)
c           else if(delt .lt. 1.0d-3)then
c           eata = dmin1(eata,10.d0)
c           else
            eata = dmin1(eata,2.d0)
c           endif
            eata = dmax1(eata,0.5d0)
         endif
         if(DEBUG .eq. 1)
     &        print*,"first limit",eata,eata*delt
      endif
c      print*,'eata (simple limit) = ',eata

c     Step size limit from input file
c      eata = dmin1(eata,tlast/(sscale*delt))
c      eata = dmin1(eata,hmax/delt)
c      if(DEBUG .eq. 1)
c     &     print*,"Input limitor",eata,eata*delt,hmax,hmax/delt

C     It's important that we do not step past large profile
C     changes. Sometimes we will miss some nucleosynthesis on the
C     up-slope of a profile spike.
      if(tcheck .and. (.not. (inmode .eq. 1)) .and. 
     &     (.not. (inmode.eq.4))) then
C        Call the profile
c         print*,"in eata tcheck"
         tt = t+delt+eata*delt
c         i=max(1,ncurrtime-1)
c        i=1
         otau=tau
         
        if(DEBUG .eq. 1)then ! .or. ishrun.eq.13)then
           print*,'******'
           print*,'ncurrtime = ',ncurrtime,numod
           print*,'future time = ',tt+6.125869143945095d10
        endif
         dtau=0.d0
         if(inmode.eq.5)then
            j=ishrun
            dtau=0.d0
            ttau=0.d0
c            print*,ncurrtime,numod
            do i=ncurrtime+2,numod-2
               if(tmodel(i) .gt. tt)goto 8923
c               if(j.eq.13)print*,"   ",i,tmodel(i)+6.125869143945095d10,
c     &              ",   ",tt+6.125869143945095d10
               dtau=dmax1(dtau,dabs(otau-taushell(i,j)))
               if(taushell(i,j).gt.ttau)idmax=i
               ttau=dmax1(ttau,taushell(i,j))
            enddo
         else
            dtau=0.d0
            ttau=0.d0
            do i=ncurrtime+1,numod-1
               if(tprof(i) .gt. tt)goto 8923
               dtau=dmax1(dtau,dabs(otau-tauprof(i)))
               if(tauprof(i).gt.ttau)idmax=i
               ttau=dmax1(ttau,tauprof(i))
            enddo
         endif
 8923    continue
         
         if(DEBUG .eq. 1)then ! .or. ishrun.eq.13)then
            print*,"current grid point",tmodel(ncurrtime)+
     &           6.125869143945095d10
            write(6,'(a,f16.0,1pe10.3,1pe10.3,1pe10.3,i8)')
     &           'time,ttau,otau,dtau,idmax',
     &           tt-90018360555.d0,ttau,otau,dtau,idmax
         endif

         if(dtau .gt. dtaurec)then
            if(DEBUG .eq. 1) ! .or. ishrun.eq.13)
     &           print*,"Large temperature change with eata=",eata,dtau
c            eata=dtaurec*eata/dtau
c            print*,"calculating new eata"
            if(inmode.eq.5)then
               j=ishrun
               eata = (dtaurec/delt)*(tmodel(idmax)-(t+delt))/
     &              dabs((taushell(idmax,j)-tau))
               if(DEBUG .eq. 1)! .or. ishrun.eq.13)
     &              print*,(tmodel(idmax)-(t+delt)),
     &              dabs((taushell(idmax,j)-tau))
            else
               eata = (dtaurec/delt)*(tprof(idmax)-(t+delt))/
     &              dabs((tauprof(idmax)-tau))
            endif
c            print*,"done calculating new eata"

            if(DEBUG .eq. 1)! .or. ishrun.eq.13)
c            if(ishrun .gt. 13)
     &           print*,"new eata, new timestep",eata,eata*delt
         endif
c         call profil(t)
      endif
 8925 continue
c      if(DEBUG .eq. 1)
c          print*,"Temperature Changes",eata,eata*delt
      
C     Need to make sure we land on profile grid points
      if(hitgrids .and. (INMODE .gt. 1)) then
         newtime = t+delt
         nexttime = newtime+eata*delt
         iproflimit=0
         i=ncurrtime+1
c         if(i .gt. nprof)print*,i
         if(nexttime .gt. tprof(i)) then
            eata=(tprof(i)-(t+delt))/delt
            changdir=0
c           Flag for grid limit
            iproflimit=1        
            tprof(i)=tprof(i)
            if(DEBUG .eq. 1) then
c            tprof(i)=tprof(i)
               print*,"hitting point at ",tprof(i)
               print*,'newtime,eata*delt,nexttime=',
     &              newtime,eata*delt,newtime+eata*delt
            endif
         endif
      endif
c      print*,"Grid points",eata,eata*delt

C     Make sure we're not at the end
      if((t+delt)+eata*delt .gt. tlast)then
         changdir=0
         eata=(tlast-(t+delt))/delt
      endif

C     delta has to be larger than 0
      eata=dmax1(eata,deltamin/delt)

      return
      end
C======================================================================
      SUBROUTINE locate(xx,n,x,j)
C======================================================================
c     Given an array xx(1:n), and given a value x, returns a value j such that
c     x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either
c     increasing or decreasing. j=0 or j=n is returned to indicate that x is 
c     out of range [from Numerical Recipes]
C======================================================================

      integer j,n
      real*8 x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
 10   if(ju-jl.gt.1)then
          jm=(ju+jl)/2
	  if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
	      jl=jm
	  else
	      ju=jm
	  endif
      goto 10
      endif
      if(x.eq.xx(1))then
           j=1
      else if(x.eq.xx(n))then
           j=n-1
      else
           j=jl
      endif
      return
      end
C======================================================================
      SUBROUTINE welcomescreen(versionnum,versiondate)
C======================================================================
C     Subroutine for welcome screen
C======================================================================

      character*5 versionnum
      character*35 versiondate

      write(6,*)
      write(6,*)" *****************************************************"
      write(6,*)" *                Welcome to nucleo                  *"
      write(6,'(a,a5,a,a17,a)')"  *            V. ",versionnum,"  ",
     $     versiondate,"            *"
      write(6,*)" *                                                   *"
      write(6,*)" *            C. Iliadis and R. Longland             *"
      write(6,*)" *        Univ. North Carolina at Chapel Hill        *"
      write(6,*)" *****************************************************"
      write(6,*)

      return
      end
