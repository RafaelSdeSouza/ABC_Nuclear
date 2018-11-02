c*********************************************************************|
      program CN3shell
c*********************************************************************|
c
c     CLASSICAL NOVA NUCLEOSYNTHESIS
c
c     10/16/2018; author: Christian Iliadis
c
c     shell assumes INMODE=3 [exponentially decaying profile];
c     it samples:
c
c     - T_peak, rho_peak, tau_T, tau_rho
c     - WD outer core composition [CO or ONe]; notice that the 16O
c       abundance is not listed in the input, but is calculated in the
c       shell as: X(16O) = 1 - X(12C) - X(20Ne) - ...
c     - mixing between outer WD core and accreted solar matter;
c       mixing factor defined as: 1 part WD matter, ff parts solar matter 
c
c     USER INPUT 
c     - sampling ranges of these parameters [now listed in input file]
c     - uncertainty factor, uf, for observed elemental mass fractions
c       [for output of rate variation factors only]
c
c     the shell script checks if ILAST=3 [MC option] and INMODE=3 
c     [exponentially decaying profile] are set, otherwise run is
c     aborted
c  
c     EXPERIMENTAL INPUT 
c     atomic number, Z, and mass fractions of all observed elements;
c     negative mass fractions are interpreted as upper limits; input 
c     in CN3shell.dat;
c     integer flag in this file let's user input either mass fractions  
c     [0] or mass fractions divided by the hydrogen mass fraction [1] 
c
c     OUTPUT
c     * only if simulated and observed abundances agree for all observed
c     * elements within a factor of 10:
c
c     CN3shell_a.out: simulated elemental mass fractions or mass fraction
c                     relative to X_H, depending on flag in input 
c     CN3shell_b.out: sampled parameters 
c
c     * only if simulated and observed abundances agree for all observed
c     * elements within user-defined factor, uf:
c
c     CN3shell_c.out: sampled isotopic abundances and reaction rate
c                     variation factors
c
c     to find simulated elemental abundances from isotopic abundances,
c     the script adds the mass fractions of several isotopes that
c     contribute, after decay, to the abundance of a daughter element 
c     [the way these are added is site-dependent]
c
c
c=======================================================================
c   VARIABLES
      implicit double precision (a-h, o-z), integer (i-n)

      integer i,j,k,lll,nexp,nsamp,nseed,niso,nini,nreac,iflag
      integer ilast,inmode,zexp(100)

c     xexp: oberved abundance; xx: simulated abundance        
      real*8 xexp(100),xx(100)
                   
      real*8 trun,uf
      real*8 temp_min,temp_max,rhop_min,rhop_max
      real*8 timetem_min,timetem_max,timerho_min,timerho_max
      real*8 lgtemp_min,lgtemp_max,lgrhop_min,lgrhop_max
      real*8 lgtimetem_min,lgtimetem_max,lgtimerho_min,lgtimerho_max
      real*8 prob(100000)
      real*8 pro1,pro2,pro3,pro4,pro5
      real*8 xf(5000)
      real*8 x_h,x_he,x_li,x_be,x_b,x_c,x_n,x_o,x_f,x_ne,x_na,x_mg
      real*8 x_al,x_si,x_p,x_s,x_cl,x_ar,x_k,x_ca,x_sc,x_ti,x_v,x_cr

      real*8 lgtemp,lgrhop,lgtimetem,lgtimerho
      real*8 temp,rhop,timetem,timerho

      real*8 xwd_c12,xwd_o16,xwd_ne20,xwd_ne21,xwd_ne22,xwd_na23
      real*8 xwd_mg24,xwd_mg25,xwd_mg26,xwd_al27
      real*8 xwd_o16_min,xwd_o16_max          
      real*8 xwd_ne20_min,xwd_ne20_max,xwd_ne21_min,xwd_ne21_max
      real*8 xwd_ne22_min,xwd_ne22_max,xwd_na23_min,xwd_na23_max
      real*8 xwd_mg24_min,xwd_mg24_max,xwd_mg25_min,xwd_mg25_max
      real*8 xwd_mg26_min,xwd_mg26_max,xwd_al27_min,xwd_al27_max
      
      real*8 f_min,f_max,lgf_min,lgf_max,lgf,ff
      real*8 xwd(500),xsum1,xsum2,xsolar(500),xini
      real*8 ratio_max,ratio_new
      
      character*100 string,string2,string3,cdum
      character*1 choice
      character*28 nucleofile
      character*5 on(5000)
      character*13 chelp(500)
       
      logical lel(100)

c=======================================================================
c   DATA INPUT      
      open(601, file='CN3shell.dat',status='unknown')

      read(601,6000) string
      read(601,6000) string
      read(601,*) temp_min,temp_max
      read(601,*) rhop_min,rhop_max
      read(601,*) timetem_min,timetem_max
      read(601,*) timerho_min,timerho_max
      read(601,*) xwd_c12_min,xwd_c12_max
      read(601,*) xwd_ne20_min,xwd_ne20_max
      read(601,*) xwd_ne21_min,xwd_ne21_max
      read(601,*) xwd_ne22_min,xwd_ne22_max
      read(601,*) xwd_na23_min,xwd_na23_max
      read(601,*) xwd_mg24_min,xwd_mg24_max
      read(601,*) xwd_mg25_min,xwd_mg25_max
      read(601,*) xwd_mg26_min,xwd_mg26_max
      read(601,*) xwd_al27_min,xwd_al27_max
      read(601,*) f_min,f_max
      read(601,6000) string
      read(601,6000) string3      
c     iflag=0: input is in mass fractions; iflag=1: input is mass
c     fractions relative to X(H)
      read(601,'(i2)') iflag
      read(601,6000) string
      i=1                                    ! data label
c     read isotope name and experimental mass fraction
 61   read(601,*,end=60) zexp(i),xexp(i)
      i=i+1
      goto 61
 60   continue
      nexp=i-1                               ! number of data points
      close(601)

c=======================================================================
c     for reaction rate variations: 
c     output results if agreement between observed and simulated elemental 
c     abundances is better than factor uf 
      uf=1.5d0
c=======================================================================

c     for sampling, a log scale is of advantage
      lgtemp_min=dlog10(temp_min)
      lgtemp_max=dlog10(temp_max)
      
      lgrhop_min=dlog10(rhop_min)
      lgrhop_max=dlog10(rhop_max)

      lgtimetem_min=dlog10(timetem_min)
      lgtimetem_max=dlog10(timetem_max)

      lgtimerho_min=dlog10(timerho_min) 
      lgtimerho_max=dlog10(timerho_max)

      lgf_min=dlog10(f_min)
      lgf_max=dlog10(f_max)

c=======================================================================
c   USER CONSOLE INPUT
      write(6,*) 'Enter number of network samples and seed:'
      read(5,*) nsamp,nseed
       
      if(nsamp.lt.1)then
         write(6,*) 'Pick number of samples >1. Run stop.'
         stop
      endif
         
      if(nseed.gt.0)then
         write(6,*) 'Pick random seed <0. Run stop.'
         stop
      endif       
       
       
c     always no       
      write(6,*) 'Random reaction rates? (y,n)'
      read(5,*) choice

c=======================================================================
c   OUTPUT FILES

c     for simulated elemental abundances [all solutions]
      open(725,file='CN3shell_a.out',status='unknown') 
      write(725,6000) string3
      write(725,*) 'total # of samples:',nsamp
      write(725,*) '       random seed:',nseed
      write(725,*) 'sample thermonuclear rates: ',choice
      write(725,6002) (zexp(j),j=1,nexp)     

c     for simulated profile parameters [all solutions]
      open(726,file='CN3shell_b.out',status='unknown') 
      write(726,6000) string3
      write(726,*) 'total # of samples:',nsamp
      write(726,*) '       random seed:',nseed
      write(726,*) 'sample thermonuclear rates: ',choice
      write(726,6004)

c     for reaction rate variation factors [choice='y' only]
      open(720,file='CN3shell_c.out',status='unknown') 
      write(720,6000) string3
      write(720,*) 'total # of samples:',nsamp
      write(720,*) '       random seed:',nseed
      write(720,*) 'sample thermonuclear rates: ',choice

c=======================================================================
c   READ NETWORK INFO FROM nucleo.in AND nucleo.dat
c     check if ILAST=3 in nucleo.in [requirement for MC run] and
c     INMODE=3 for run with exponential profile
      open(702,file='nucleo.in',status='unknown')
      do 71 i=1,2
         read(702,6000) string
 71   continue
      read(702,*) ilast
      if(ilast.ne.3)then
         write(6,*) 'ILAST must be =3. Run stop!'
         stop
      endif
      do 77 i=1,6
         read(702,6000) string
 77   continue
      read(702,*) inmode
      if(inmode.ne.3)then
         write(6,*) 'INMODE must be =3. Run stop!'
         stop
      endif
      close(702) 

c     get # of reactions and # of initial abundances
      open(703,file='nucleo.in',status='unknown')
      do 73 i=1,5
         read(703,6000) string
 73   continue
      read(703,*) nucleofile
      close(703) 
          
      open(704,file=nucleofile,status='unknown')
      read(704,*) niso  
      read(704,6000) string           
      do 74 i=1,niso
         read(704,6000) string
 74   continue
      read(704,*) nini
      do 75 i=1,nini
         read(704,6000) string
 75   continue
      read(704,*) nreac             ! nreac: number of reactions
      close(704) 

c=======================================================================
c   CHECK: IS SAMPLE SIZE OF RANDOM NUMBERS SUFFICIENTLY SMALL?
c     sample: reactions [optional], profile parameters, etc.
      if(nseed*(nreac+4).gt.1.d8) then
         write(6,*) 'Too many calls to random number generator!'
         stop
      endif

c=======================================================================
c   SET ALL RATE PROBABILITIES EQUAL TO ZERO AT THIS STAGE
      open(705,file='prob.dat',status='unknown')          
      do 76 j=1,nreac
         prob(j)=0.d0
         write(705,1000) j,prob(j)
 76   continue
      close(705)

c=======================================================================
c   LOOP OVER EACH NETWORK CALCULATION   --->
c=======================================================================
      do 100 k=1,nsamp
c        print counter
         write(6,*) ' Run # ',k    ! do not use symbol 'k' anywhere else
c                                    in this loop !
c        ===============================================================
c        randomize profile parameters
         lgtemp=lgtemp_min+(lgtemp_max-lgtemp_min)*ran1(nseed)
         temp=10.d0**lgtemp

         lgrhop=lgrhop_min+(lgrhop_max-lgrhop_min)*ran1(nseed)
         rhop=10.d0**lgrhop

         lgtimetem=lgtimetem_min+(lgtimetem_max-lgtimetem_min)
     &             *ran1(nseed)
         timetem=10.d0**lgtimetem

         lgtimerho=lgtimerho_min+(lgtimerho_max-lgtimerho_min)
     &             *ran1(nseed)
         timerho=10.d0**lgtimerho

c        ===============================================================
c        enter sampled profile parameters into nucleo input file 
         call SYSTEM('cp nucleo.in nucleo_shell.in')

         open(801,file='nucleo.in',status='unknown')
         open(802,file='nucleo_shell.in',status='unknown')

c        read nucleo_shell.in and construct new nucleo.in
         do 102 j=1,12
            read(802,6000) string
            write(801,6000) string 
 102     continue
         read(802,*) pro1,pro2,pro3,pro4,pro5   ! pro3: running time
         write(801,1100) temp,rhop,pro3,timetem,timerho

         do 103 j=1,500
            read(802,6000,end=105) string
            write(801,6000) string
 103     continue
 105     continue         
          
         close(801)
         close(802)

c        ==============================================================
c        sample WD outer core mass fractions and construct new 
c        nucleo.dat
         xwd_c12=xwd_c12_min+(xwd_c12_max-xwd_c12_min)
     &                *ran1(nseed)
         xwd_ne20=xwd_ne20_min+(xwd_ne20_max-xwd_ne20_min)
     &                *ran1(nseed)
         xwd_ne21=xwd_ne21_min+(xwd_ne21_max-xwd_ne21_min)
     &                *ran1(nseed)
         xwd_ne22=xwd_ne22_min+(xwd_ne22_max-xwd_ne22_min)
     &                *ran1(nseed)
         xwd_na23=xwd_na23_min+(xwd_na23_max-xwd_na23_min)
     &                *ran1(nseed)
         xwd_mg24=xwd_mg24_min+(xwd_mg24_max-xwd_mg24_min)
     &                *ran1(nseed)
         xwd_mg25=xwd_mg25_min+(xwd_mg25_max-xwd_mg25_min)
     &                *ran1(nseed)
         xwd_mg26=xwd_mg26_min+(xwd_mg26_max-xwd_mg26_min)
     &                *ran1(nseed)
         xwd_al27=xwd_al27_min+(xwd_al27_max-xwd_al27_min)
     &                *ran1(nseed)

c        oxygen abundance given by 1 minus the other abundances
         xwd_o16=1.d0-xwd_c12-xwd_ne20-xwd_ne21-xwd_ne22-xwd_na23
     &               -xwd_mg24-xwd_mg25-xwd_mg26-xwd_al27

         call SYSTEM('cp nucleo.dat nucleo_shell.dat')

         open(223,file='nucleo.dat',status='unknown')
         open(224,file='nucleo_shell.dat',status='unknown')

c        !!!!!!!!!!!
c        following index values apply to INITIAL abundances;
c        they do not label the isotope number
c        !!!!!!!!!!!
                   
c        abundances of white dwarf
c        ...first, set everything to zero                
         do 141 jj=1,nini
            xwd(jj)=0.0d0
 141     continue
c        ...then input WD abundances
         xwd(9)=xwd_c12      ! 12C
         xwd(13)=xwd_o16     ! 16O
         xwd(17)=xwd_ne20    ! 20Ne
         xwd(18)=xwd_ne21    ! 21Ne
         xwd(19)=xwd_ne22    ! 22Ne
         xwd(21)=xwd_na23    ! 23Na
         xwd(22)=xwd_mg24    ! 24Mg
         xwd(23)=xwd_mg25    ! 25Mg
         xwd(24)=xwd_mg26    ! 26Mg
         xwd(26)=xwd_al27    ! 27Al

c        ...draw random mixing fraction [1 part WD with ff parts solar];
c        defined as: X_mix = (X_WD + ff*X_solar)/(1+ff)

         lgf=lgf_min+(lgf_max-lgf_min)*ran1(nseed)
         ff=10.d0**lgf
             
         do 901 jj=1,niso+3
            read(224,6000) cdum
            write(223,6000) cdum
 901     continue     

c        test to see if "solar" abundances are used in input; if
c        not, stop run
         if(index(cdum,'solar').eq.0)then
            write(6,*) 'nucleo.dat must list solar initial abundances!'
c           replace temporary nucleo.in by help file
            call SYSTEM('rm nucleo.in')
            call SYSTEM('mv nucleo_shell.in nucleo.in')

c           replace temporary nucleo.dat by help file
            call SYSTEM('rm nucleo.dat')
            call SYSTEM('mv nucleo_shell.dat nucleo.dat')
            stop        
         endif

c        solar abundances read from nucleo.dat stored in xsolar;
c        mix with WD matter 
         xsum1=0.d0
         xsum2=0.d0
         do 143 jj=1,nini
              read(224,8050) chelp(jj),xsolar(jj)
              xsum1=xsum1+xsolar(jj)
 143     continue
c        normalize all solar mass fractions from nucleo.dat
         do 147 jj=1,nini
            xsolar(jj)=xsolar(jj)/xsum1
            xini=(xwd(jj)+ff*xsolar(jj))/(1.d0+ff)     
            write(223,8051) chelp(jj),xini
            xsum2=xsum2+xini            ! check mass fraction sum
 147     continue
 
c        output on screen
         write(6,9555) temp,rhop,ff,xsum1,xsum2 

         do 148 jj=1,nreac+1
             read(224,6000) cdum
             write(223,6000) cdum
 148     continue          
          
         close(223)
         close(224)
                       
c        ==============================================================
c        randomize rates (choice.eq.'y')
         if(choice.eq.'y')then
            open(706,file='prob.dat',status='unknown')            
            do 104 j=1,nreac
               prob(j)=gasdev(nseed)
               write(706,1000) j,prob(j)
 104        continue
            close(706)
         endif

c        ==============================================================
c        run nucleo     
         call SYSTEM('./nucleo')
c        ==============================================================

c        see if "RUN STOP OKAY"          
         lll=0
         open(707,file='nucleo.out',status='unknown')
         do 106 j=1,6
            read(707,6000,end=200) string
 106     continue
          
         lll=index(string,'RUN STOP OKAY')
        
         if(lll.eq.0)then                ! network does not complete
            goto 200
           else                          ! read iso abundances for later output
            do 107 j=1,2
               read(707,6000) string     ! dummy
 107        continue
            read(707,1099) string2,trun  ! read running time
            do 108 j=1,3
               read(707,6000) string     ! dummy
 108        continue
            read(707,1001) (on(j),xf(j),j=1,niso)
         endif
         close(707)

c        ==============================================================
c        analyze results          

c        first task is to add abundances of stable isotopes of a given
c        element [summation of stable and very long-lived species in rows]

         x_h=xf(1)+xf(3)       ! H = 1H + 2H,
         x_he=xf(5)+xf(6)      ! He = 3He + 4He, etc.
         x_li=xf(10)+xf(12)
         x_be=xf(16)
         x_b=xf(19)+xf(24)
         x_c=xf(7)+xf(30)
         x_n=xf(34)+xf(36)
         x_o=xf(8)+xf(44)+xf(48)    
         x_f=xf(53)         
         x_ne=xf(56)+xf(62)+xf(65)
         x_na=xf(69)
         x_mg=xf(74)+xf(79)+xf(82)
         x_al=xf(92)
         x_si=xf(96)+xf(101)+xf(105)
         x_p=xf(110)
         x_s=xf(114)+xf(119)+xf(120)+xf(130)
         x_cl=xf(128)+xf(138)
         x_ar=xf(135)+xf(144)+xf(152)
         x_k=xf(146)+xf(153)+xf(154)
         x_ca=xf(149)+xf(161)+xf(166)+xf(170)+xf(178)+xf(190)
         x_sc=xf(175)
         x_ti=xf(179)+xf(187)+xf(193)+xf(197)+xf(200)
         x_v=xf(203)+xf(204)
         x_cr=xf(201)+xf(207)+xf(210)+xf(212)
     
c        next task: add abundances of radioactive nuclides that
c        contribute to abundance of a given element; 
c        this part may need to be modified, depending on the 
c        scenario; see also comments in the header;
c        index of xx() is equal to atomic number!

         xx(1)=x_h+xf(2)+xf(4)    ! H = 1H + 2H + 3H + n
         xx(2)=x_he               ! He = 3He + 4He
         xx(3)=x_li               ! Li = 6Li + 7Li 
         xx(4)=x_be+xf(11)+xf(18) ! Be = 7Be + 9Be + 10Be
         xx(5)=x_b+xf(22)         ! B = 10B + 11B + 11C
         xx(6)=x_c+xf(33)+xf(27)  ! C = 12C + 13C + 14C + 13N
         xx(7)=x_n+xf(31)+xf(35)  ! N = 14N + 15N + 14O + 15O
         xx(8)=x_o+xf(50)         ! O = 16O + 17O + 18O + 18F    
         xx(9)=x_f         
         xx(10)=x_ne
         xx(11)=x_na+xf(66)       ! Na = 23Na + 22Na
         xx(12)=x_mg
         xx(13)=x_al+xf(88)       ! Al = 27Al + 26Al
         xx(14)=x_si
c        P = 31P + 31Si + 23P + 33P
         xx(15)=x_p+xf(109)+xf(115)+xf(116)       
c        S = 23S + 33S + 34S + 36S + 35S 
         xx(16)=x_s+xf(129)       
c        Cl = 35Cl + 37Cl + 36Cl
         xx(17)=x_cl+xf(132)      
c        Ar = 36Ar + 38Ar + 40Ar + 37Ar + 39Ar
         xx(18)=x_ar+xf(139)+xf(148)              
c        K = 39K + 40K +41K + 41Ar
         xx(19)=x_k+xf(156)       
c        Ca = stable Ca + 41Ca + 45Ca + 47Ca + 43Sc + 44Sc + 42K + 43K
         xx(20)=x_ca+xf(155)+xf(174)+xf(185)+xf(163)+xf(169)
     &      +xf(162)+xf(165)      
c        Sc = 45Sc + 45Ti  + 46Sc
         xx(21)=x_sc+xf(173)+xf(183)      
c        Ti = stable Ti + 44Ti + 47Sc + 48Sc + 47V
         xx(22)=x_ti+xf(167)+xf(184)+xf(194)+xf(188)      
c        V = 50V + 51V + 48V + 49V + 48Cr + 49Cr
         xx(23)=x_v+xf(195)+xf(198)+xf(191)+xf(196)         
         xx(24)=x_cr+xf(206)    ! Cr = stable Cr + 51Cr
         
c        ==============================================================
c        if iflag=1: 
c        divide all computed elemental abundances by hydrogen abundance
         if(iflag.eq.1)then
            if(xf(1).lt.0.d0) goto 200      ! stop run if X(1H)<0.0
            do 111 jj=1,24
               xx(jj)=xx(jj)/x_h
c               write(6,*) xx(jj)
 111        continue
         endif
         
c        ==============================================================
c        determine maximum ratio of observed and computed abundances
         ratio_max=1.0d0
         do 113 j=1,nexp   
            if(xexp(j).gt.0.d0)then         ! no upper limit input
               ratio_new=xx(zexp(j))/xexp(j)
               if(ratio_new.le.1.0d0)then       
                  ratio_new=1/ratio_new     ! inverse to compare factors >1
               endif
             else                           ! upper limit input [negative value]
               ratio_new=xx(zexp(j))/dabs(xexp(j))
               if(ratio_new.le.1.0d0)then   ! simulation agrees with upper limit;
                  ratio_new=1.0d0           ! nothing else can be concluded
               endif
            endif
            if(ratio_new.ge.ratio_max)then
               ratio_max=ratio_new
            endif
 113     continue

c        output sampled abundances and parameters if agreement is achieved
c        for all abundances within a factor of 10
         if(ratio_max.le.10.0d0)then
c           elemental abundances 
            write(725,6001) (xx(zexp(j)),j=1,nexp)

c           profile parameters
            write(726,1103) temp,rhop,pro3,timetem,timerho,
     &                   xwd_c12,xwd_o16,ff,ratio_max
         endif

c        ==============================================================              
c        output sampled reaction rate variation factors if choice='y' 
c        and all solutions are within user-defined factor, uf    
         if((ratio_max.le.uf).and.(choice.eq.'y'))then
            write(720,9211)
            write(720,9214) k
            write(720,9999) (on(j),xf(j),j=1,niso)
            write(720,9216)

c           read new prob.dat file created by nucleo [with proper
c           forward/reverse factor sampling]
            open(730,file='prob.dat',status='unknown')
            do 219 j=1,nreac
               read(730,9218) prob(j)
 219        continue
            close(730)
c           output probabilities
            write(720,9217) (j,prob(j),j=1,nreac)
         endif
         
c        ==============================================================
c        if network does not complete, no solution is obtained, or 
c        hydrogen mass fraction becomes negative:
 200     continue      
         if((lll.eq.0).or.(xh.lt.0.d0))then      
            close(707)
         endif                  
                                                    
c        replace temporary nucleo.in by help file
         call SYSTEM('rm nucleo.in')
         call SYSTEM('mv nucleo_shell.in nucleo.in')

c        replace temporary nucleo.dat by help file
         call SYSTEM('rm nucleo.dat')
         call SYSTEM('mv nucleo_shell.dat nucleo.dat')

c=======================================================================
 100  continue     ! <---- NEXT NETWORK SAMPLE
c=======================================================================
      close(720)
      close(725)
      close(726)
      
c=======================================================================
c   FORMATS

 1000 format(1x,i5,'=',1pe20.12)          
 1001 format(9(1x,a5,1x,1pe10.2E3))
 1099 format(a37,1pe10.3)
 1100 format(4(1pe9.2e2,','),1pe9.2e2)
 1103 format(9(1pe9.2e2,2x))
 6000 format(a100)
 6001 format(25(1pe10.2e3,2x))
 6002 format(5x,25(i2,10x))
 6004 format(2x,'T_peak',4x,'rho_peak',3x,'time_run',4x,'tau_T',
     &       5x,'tau_rho',4x,'xwd_c12',4x,'xwd_o16',6x,'ff',7x,
     &       'ratio_max')
 8050 format(5X,A13,1x,E30.3)   
 8051 format(5X,A13,7x,1pe10.3)
 9211 format(153('$'))
 9214 format(1x,'FINAL ABUNDANCES','       RUN: ',i6)
 9216 format(1x,'PROBABILITIES')
 9217 format(9(1x,i5,'=',1pe10.2e3))
 9218 format(7x,1pe20.12) 
 9555     format(5x,'T_peak: ',E10.3,5x,'rho_peak: ',
     &       E10.3,5x,
     &       'pre-mixing: ',E10.3,5x,'Xsum_sol: ',E17.10,5x,
     &       'Xsum_premix: ',E17.10)
 9999 format(9(1x,a5,'=',1pe10.2e3))  

      end

c*********************************************************************|
c*********************************************************************|
      FUNCTION gasdev(idum)
c*********************************************************************|
      integer idum
      real*8 gasdev
C     uses ran1
C         Returns a normally distributed deviate with zero mean and unit 
C         variance, using ran1(idum) as the source of uniform deviates
      integer iset
      real*8 fac,gset,rsq,v1,v2
      real*8 ran1
      save iset,gset
      data iset/0/
      if(iset.eq.0) then
  1      v1=2.d0*ran1(idum)-1.d0
         v2=2.d0*ran1(idum)-1.d0
         rsq=v1**2.d0+v2**2.d0
         if(rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
         fac=sqrt(-2.d0*dlog(rsq)/rsq)
         gset=v1*fac
         gasdev=v2*fac
         iset=1
       else
         gasdev=gset
         iset=0
      endif
      Return
      end
c*********************************************************************|
      FUNCTION ran1(idum)
c*********************************************************************|
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


