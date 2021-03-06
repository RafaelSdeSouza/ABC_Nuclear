C     ************************************************************
C     *                                                          *
C     *      N U C L E O S Y N T H E S I S         C O D E       *
C     *   (from Nikos Prantzos, received at TRIUMF Apr. 1995)    *
C     *                                                          *
C     ************************************************************
C
C    SUMMER '95
C    ````````````````````````````````
C  * Reaction rates are no longer in the subroutine SRATE.
C    SRATE now opens the file REACTIONS.DAT and reads in the
C    parameters for each reaction from the Thielemann library.
C    REACTIONS.DAT is made from a subset of REACLIB.DAT using
C    the SEARCH program.  For further information, see the 
C    documentation for the PROCDATA program.
C
C  * When running with constant temperature/density the reaction
C    rates are only calculated once, instead of at each time step.
C
C  * Subroutines COUFRA, DESCOUP, SMOK have been removed as they
C    are no longer used in rate calculation.
C
C  * Each reaction flux is time integrated as the program runs and
C    the total integral is output into TOTFLUX.OUT.
C
C  * The program now asks the user if she/he wants to change the
C    value of EPS,YTMIN,NLAST, or the file consulted for the 
C    reaction rates.  Unless a 'y' or 'Y' is entered in either of
C    these cases, the program proceeds using the defaults.
C
C    -Paul Tupper
C
C
C    SUMMER '97
C    ````````````````````````````````   
C  * The program outputs now the total energy generation (erg s-1 g-1) 
C    when the results are stocked every NLAST/KP time steps in NUCLEO.OUT.
C
C  * The program outputs now the total time-integrated energy generation 
C    (erg g-1) after the end of the calculation in NUCLEO.OUT.
C
C    -CI
C
C
C    SPRING '98
C    ````````````````````````````````   
C  * If INMODE=2 an external temperature-density profile is read from a
C    file. The format has to be TPROF(s), TAUPROF(K), RHOPROF(g/cm3). The
C    subroutine PROFIL finds the temperature TAU and density RHO for the next
C    time-step integration by interpolation. The profile MUST start with 
C    TPROF=0 s. Reaction rates are calculated only if TAU or RHO change
C    significantly (say, by 0.1%).
C       Interpolated values for TAU and RHO are printed for each timestep
C    NSTEP at time T to an output file PROFILE.OUT in order to check if 
C    the program used the profile properly.
C
C    (Remark: in nucleo.out T is the new time calculated with the old
C             values for TAU and RHO).
C
C    -CI
C
C
C    SPRING 2000
C    `````````````````````````````````
C  * Multiplicative factors for rate changes can now be entered in an
C    input file.
C
C  * Corrected a small bug in code. X(ISF) is now the largest initial
C    abundance (the code stops if X(ISF).le.XLAST). For cases where H
C    was not the most abundant isotope, program would default to 3He
C    burning. It now will properly track the most abundant species.
C
C    -AEC
C
C
C    FALL 2000
C    ```````````````````````````````
C  * If INMODE=3 a parametrized temperature-density profile is used
C    (an exponentially decaying tau-rho profile). Input values are 
C    the initial temperature and density, the time TLAST, and a scale 
C    factor that determines the expansion time. The subroutine PROFIL 
C    finds the temperature TAU and density RHO for the next time-step
C    by using the relations (see, for example, Arnett: Supernovae and
C    Nulceosynthesis, p. 252/253):
C
C            tau=tau0*exp(-t/(3*thd))
C            rho=rho0*exp(-t/thd)
C            thd=446.0*scale/sqrt(rho0)
C
C    Reaction rates are calculated only if TAU or RHO change significantly 
C    (say, by 0.1%).
C       Parametrized values for TAU and RHO are printed for each timestep
C    NSTEP at time T to an output file PROFILE.OUT in order to check if 
C    the program calculated the profile properly.
C
C  * Program now stops if the temperature drops below T9=0.001.
C
C    -CI
C
C
C    SPRING 2001
C    `````````````````````````````````
C  * Program is now reading input from file instead of terminal
C
C  * The code allows now for multiple calculations for different
C    sets of conditions: 
C
C    - for a single calculation, ILAST=0; in this case, a long 
C      output is written to NUCLEO.OUT. 
C
C    - for multiple calculations, ILAST=1 for all runs except the 
C      last run (ILAST=2); in this case, NUCLEO.OUT contains for 
C      each calculation only info on reaction rate changes, total 
C      energy production and final abundances.
C
C    -CI
C
C
C    SPRING 2004
C    `````````````````````````````````
C  * In the input file, 6 isotopes can be specified for output of
C    abundances versus time(step) at the end of NUCLEO.OUT; they must be 
C    labeled by their isotope number (given in the .dat input file); 
C    this output is much more precise than the one given in the block
C    output (that is analyzed by ISOBAUND) and is suitable, for example, 
C    for plots of abundances vs. consumed hydrogen.
C
C    -CI
C
C
C    FALL 2004
C    `````````````````````````````````
C  * Program now outputs total reaction energy for every link at end of
C    calculation; output is in erg/g and in %  
C
C    -CI
C
C
C    SUMMER 2005
C    `````````````````````````````````
C  * Commented out condition for initial time step DELTA 
C    ("DELTA=TLAST*1.d-10"); previously, 1.d-6 was used and the program 
C    frequently found no convergence (if TLAST was incidentally too large); 
C    DELTA (in s) as well as XLAST are now specified in NUCLEO.IN 
C
C  * Neutron excess parameter ETA is printed in output
C
C    -CI
C
C
C    SPRING 2009
C    `````````````````````````````````
C  * Changed reaction rate input from fitting parameters to tabular format
C    in REACTIONS.DAT. The rates are interpolated using a cubic spline 
C    [routine is adopted from Numerical Recipes]. The rates are now read
C    once in the beginning into arrays [i.e., the rates are stored in RAM 
C    instead of an external file on the disk]; this significntly speeds 
C    up program execution.
C
C    -CI   
C
C
C    SUMMER 2009
C    `````````````````````````````````
C  * All links labeled 'ec' in REACTIONS.DAT are now multiplied by rho*Ye.
C
C    -CI
C
C
C    WINTER 2009
C    `````````````````````````````````
C  * Changed maximum number of iterations from 20 to 40.
C
C    -CI
C
C
C    SPRING 2010
C    `````````````````````````````````
C  * Sometimes, when code runs with inmode=2, the external profile is 
C    not exactly mapped; this is especially important for running near
C    peak temperature. Now code determines max. T for external profile
C    and maximum mapped T, and outputs their difference at end of
C    calculation.    
C
C    -CI
C
C
C    WINTER 2011
C    `````````````````````````````````
C  * Added two running options:
C
C    INMODE=4
C    If a region is fully convective AND if the convective timescale in
C    a given shell, given by the ratio of mixing length and convection
C    velocity, is shorter than the lifetimes of nuclear reactions, then 
C    it can be shown that the full hydro run should be close to a post-
C    processing network calculation in which the rates are being mass-
C    averaged over all shells of a given model (i.e., at given time).
C    This reflects the "instantaneous mixing approximation". 
C
C    Each time the subroutine "instmix" is called, a test is performed:
C    if the convection time in a given shell is longer than a nuclear
C    mean lifetime, then both values are written to an output file 
C    [instmixtest.out] for closer inspection.
C
C    INMODE=5
C    Postprocessing calculations are performed sequentially over each
C    shell and at the end the mass fractions are weighted by the mass 
C    of each shell before being summed to find the average mass fraction
C    for each species. This reflects the "no mixing approximation". With
C    this option, the ILAST entry in the input file is meaningless. 
C
C    All individual runs must be performed for the same time, which is 
C    chosen to be TLAST. If some runs terminate earlier, a warning appears. 
C    In that case, reduce TLAST and redo the run.
C       
C
C    Input for both modes is given in shellinputfile, listing:
C
C    #of shells, time for model 1
C    T, rho, mix.length, conv. velocity, mass of shell 1
C    T, rho, mix.length, conv. velocity, mass of shell 2
C    T, rho, mix.length, conv. velocity, mass of shell 3
C    etc.
C    #of shells, time for model 2
C    T, rho, conv. velocity, mix. length, mass of shell 1
C    T, rho, conv. velocity, mix. length, mass of shell 2
C    T, rho, conv. velocity, mix. length, mass of shell 3
C    T, rho, conv. velocity, mix. length, mass of shell 4
C    etc.
C    etc.
C
C    [Units:
C    ---, seconds
C    K, g/cm3, cm/s, cm, rel. mass]
C    Shells should be ordered from interior to exterior [usually meaning
C    from hot to cold].
C    
C    First "time" value in input file must be t=0; and running time 
C    (TLAST) should not extend beyond last time value given in input.
C
C    When the input stellar models [i.e., input for different times] do
C    not have the same number of shells, then the shells need to be re-
C    sized so that the code runs always with the same number of shells at
C    each given time. The standard mass grid is constructed in the following
C    manner: it is adopted from the model with the largest number of
C    shells.
C
C    The original and re-sized grids (together with temperature) are output 
C    in NEWGRID.OUT. The code stops if the original models differ in TOTAL
C    mass by more than 5%.
C 
C    -CI
C
C
C    SPRING 2011
C    `````````````````````````````````
C  * Interpolation (cubic spline) gave wrong results near 10 MK [the boundary
C    below which most rates are set equal to zero. Changed algorithm to simple
C    linear interpolation between logarithmic values. 
C
C    -CI
C
C
C    SPRING 2012
C    `````````````````````````````````
C  * Implemented Monte Carlo sampling of reaction rates (ILAST=3); this applies
C    only to INMODE=1,2,3,4, but not to INMODE=5; for each run [constant T, or
C    T-rho profile, or shells] rates are sampled only once and network is evolved
C    with this fixed set of samples; note, however, that the factor uncertainty  
C    is *not* constant, rate = (e^mu) * (e^sigma)**p, because sigma is changing
C    with temperature [e^mu: median rate; e^sigma: factor uncertainty];
C        we chose a "flat" parameterization: for each run, a probability
C        factor p is sampled for each reaction; for a given reaction,
C        p has the same value at all temperatures
C    
C  * Rates [and probability samples, if MC option is selected] have been moved
C    to subroutine RNETWORK
C
C  * Final abundances are now given with more significant figures
C
C
C    SUMMER 2012
C    `````````````````````````````````
C  * Implemented Monte Carlo sampling for INMODE=5; the same random numbers are 
C    used from shell to shell [i.e., probabilities p_i are only sampled once for
C    shells], but rates are sampled randomly from run sample to run sample
C
C
C    WINTER 2012
C    `````````````````````````````````
C  * Warnings are now printed on screen if, for INMODE=2 or 5, the code does not
C    track peak temperature properly; peak temperature values of external profile
C    and of the calculation are also printed on screen  
C
C
C    WINTER 2013
C    `````````````````````````````````
C
C  * If the time, T, for the previous step plus DELTA exceeds the user specified 
C    running time, TLAST, then code picked DELTA as the difference between
C    TLAST and T:
C    
C    IF((T+DELTA).GE.TLAST) DELTA=TLAST-T
C 
C    This means that, if DELTA does not change in the next iteration, then 
C    T'=T+DELTA=TLAST and the code terminates. If there is a huge rise in T during 
C    this large time step, the code will not find it. 
C
C    Thus, above line was replaced by:
C
C    IF((T+DELTA).GE.TLAST) DELTA=(TLAST-T)/1.1
C
C    Now we are forcing DELTA to be less than simply TLAST-T, and when it will be 
C    added to the time, T, the sum will not simply equal TLAST, meaning the code 
C    should not simply jump over that last step.
C
C  * Replaced line in code by:
C
C    DELTA=MIN(DELTA,TLAST/DELFAC)
C
C    Parameter Delfac is now set in nucleo.in.   
C
C    
C    SPRING 2014
C    `````````````````````````````````
C    versionnum = "2.2"
C    versiondate='February 26, 2014'
C
C  * Richard Longland implemented Gear's method as network solver. As a result, 
C    format of input file has changed slightly. The parameters EPS and YTMIN have
C    a different meaning with Gear's method. The Gear parameters are:
C
C    errorscale: entered instead of EPS in nucleo.in; convergence criterium on 
C                Newton-Raphson iteration [controls step size]
C
C    esc:        entered instead of YTMIN in nucleo.in; set value between 1e-8 and
C                1e-12; species with Y>YTMIN are weighted equally in convergence 
C                procedure; very small values will cause it to fail 
C
C    etascale:   set near 0.25, but can be changed in nucleo.in
C
C    thresh:     hardwired to 1e-3
CC
C  - ma48d.f and ddeps.f are linear algebra solvers
C  - better profile tracking to make sure no temperature spikes are missed
C  - IJD, NDS are not needed for Gear's method runs [still needed for Wagoner runs]
C  - NLAST, XLAST, DELTA0 work same for both methods
C  - DELFAC is ignored with Gear's method
C
C  * Welcome screen has been changed. We now keep track of version numbers. 
C
C  * Ye is now listed in nucleo.out
C
C  * Code has been changed to allow for input of *stellar* weak rates (IWEAK=1);
C    this is not without problems and care should be taken when using this option;
C
C    main problems at this point:
C
C    (i)   for T < 1 MK, the 7Be decay rate is set equal to zero; this is not
C          correct but derives from a consistent treatment of rates labeled 'ec'
C          [the other ones are pep and 3He-->t decay; in both of these cases the
C          above assumption is correct]
C
C    (ii)  when *stellar* weak rates are used and either T < 1 MK or rho*Ye <
C          10 g/cm3, then the lab rates are used; this is not correct either
C          since (a) the lab rates do correspond to some bound electron density 
C          near the nucleus [would be interesting to find out what the values
C          are], and (b) it is not clear if at very low density the atoms are
C          completely ionized [no electron capture] or are neutral [lab rate]
C
C    (iii) the laboratory beta-delayed particle decay rates are not replaced by
C          IWEAK=1 since the stellar rates for these links are not known; in
C          reality these will also be affected by T and rho in a plasma; perhaps
C          one could estimate those for the stellar case based on the pure 
C          decay/delayed particle rate ratio in the lab, but this is not clear
C
C    again: careful when stellar rates are used with shells where T or
C           rho*Ye become very small 
C
C    note: *** the factor uncertainty of the stellar weak rates has been
C          hard-wired to a factor of 2! ***    
C
C  * Identical grid points are identified and corrected by adding a 10%
C    DeltaT
C
C
C    SPRING 2014
C    `````````````````````````````````
C
C    versionnum = "2.3"
C    versiondate='April 19, 2014'
C
C  * For ILAST=3 [MC option], forward and reverse rates are not sampled
C    independently anymore; a given reverse rate [labeled by "v"] has now the 
C    same probability factor p as the corresponding forward rate
C
C
C    SUMMER 2014
C    `````````````````````````````````
C
C    versionnum = "2.4"
C    versiondate='July 18, 2014'
C
C  * A number of table searches that were performed "manually" [for profile and 
C    shell input] are now done using subroutine 'locate'
C
C  * Variable 'numprof' seems to duplicate 'numod'; former variable has been 
C    removed
C
C  * INMODE=4,5 is modified to allow for just the outer shells to burn; only shells
C    with cumulative mass between new input variable 'dconv' [in relative mass
C    units, same as profile input] and the surface burn;
C
C     0.1  0.2  0.1  0.3  0.2   0.4
C    |   |     |   |     |   |      |
C    | 1 |  2  | 3 |  4  | 5 |   6  |
C    |   |     |   |     |   |      |
C      0.1   0.3 0.4   0.7 0.9    1.3   [cumulative mass]
C
C    Example: dconv=0.71; then shell counting variable ndep=4 and only active shells
C    are **exterior** to ndep, i.e., shells #5 and #6
C
C    dconv=0.d0 in nucleo.in means that all shells will burn [new option disabled]
C
C
C    SPRING 2015
C    `````````````````````````````````
C
C    versionnum = "2.5"
C    versiondate='April 28, 2015'
C
C  * Removed multi-run options ILAST=1 and =2; sequential runs are best done
C    using a shell script that wraps around nucleo; cleaned up logic and loops
C    of conditions
C
C  * Removed Monte Carlo loop from nucleo; this is also best done in a shell
C    script; removed IPRINT variable as well; new organization of output:
C    ILAST=0: long output
C    ILAST=3: short output in nucleo.out and on screen; in addition, the
C             following actions:
C             - "prob.dat" is opened and probability factors are read for all 
C               reactions
C             - rates are multiplied by: rate(j,i)=rate(j,i)*(fu(j,i)**prob(j))
C               [j: reaction; i: temperature]
C             - new "prob.dat" is written with proper forward/reverse values
C
C     *** ILAST=3 *** must be set for the MC shell to work properly ***
C
C    [both options initiate a single network run when ./nucleo is run; i.e., 
C     setting ILAST=3 does not by itself result in a MC run; this needs to be 
C     initiated in shell script]
C
C
C    WINTER 2016
C    `````````````````````````````````
C
C    versionnum = "2.6"
C    versiondate='December 10, 2016'
C
C *  Changed slightly the output so that energy and neutron excess can be
C    analyzed by MCshell.f
C
C *  ISF was the isotope number of the most abundant species and the code
C    stopped running if its abundance dropped below XLAST; this is a problem,
C    for example, in carbon burning, where **oxygen**, in fact, is the most 
C    abundant species; ISF [the number of the nuclide whose abundance falling
C    below XLAST will stop the code] is now entered in nucleo.in
C
C
C    SPRING 2017
C    `````````````````````````````````
C
C    versionnum = "2.7"
C    versiondate='February 4, 2017'
C
C *  BUG ALERT! multi-zone running, INMODE=4 or 5: the cumulative mass of all
C    zones was only calculated correctly if the mass zones had to be resized;
C    in Jordi's nova profiles, for example, the mass of the zones do not change
C    from model to model [i.e., for different times]; in that case, the code
C    still used the cumulative mass distribution of the "templet" that, however,
C    is only calculated for multi-zone profiles with zone mass distributions 
C    that change with time; error was corrected.
C
C *  suppose you want to inflate the rate uncertainty factor, f.u., by an
C    "unknown systematic effect" uncertainty in a MC network study; included 
C    an option where vmult [previously the rate multiplication factor listed  
C    at the end of nucleo.in] becomes this inflation factor in a MC network 
C    calculation if INMODE=3 is selected; normally we get [see Longland, AA 548, 
C    A30 (2012)]:
C
C      x(T) = exp[mu(T)] * exp[sig(T)]**p(T) = x_med * [f.u.**p(T)]
C
C     inflating the factor uncertainty, f.u., by an additional factor, f_a,
C     corresponds then to:
C
C      x_a(T) = x_med * [f_a * f.u.]**p(T) = x_med * [f_a**p(T)] * [f.u.**p(T)]
C
C    when this option is used, forward and reverse rate should always be
C    entered with the same factor f_a in nucleo.in
C
C *  changed INMODE=3 option: up to now, the adiabatic expansion was hard-wired,
C    Eq. (5.144) in NPoS, with the hydrodynamic free-fall time given by Eq.
C    (5.145); the problem is that T and rho are correlated, which does not leave
C    enough freedom to vary the expansion time scales for T and rho independently;
C
C    the code has been changed so that the user can enter independent time scales
C    for T and rho, according to:
C
C    T = T0 exp(-t/scale1)    and     rho = rho0 exp(-t/scale2)
C
C    where scale 1 and scale2 are the times at which T and rho, respectively, have 
C    fallen to 1/e of the peak values;
C
C    the old options can still be used, e.g.:
C
C    - for SNIa, where we assumed a "fixed" adiabatic expansion [i.e., without 
C      assuming a hydrodynamic free fall], use: scale1=3*[desired expansion time] 
C      and scale2=[desired expansion time] 
C    - for SNII, where we assumed both an adiabatic expansion and a hydrodynamic
C      expansion time, use: scale1=3*scale2 and scale2=446.0/sqrt(rho0)
C
C *  changed INMODE=2 option: besides profile filename and running time, three
C    more parameters can be set by the user: these are the factors by which 
C    temperature, density and the time axis is multiplied [for stretching or 
C    compressing external profile]; if a run without stretching/compressing the 
C    external profile is desired, set these three input parameters to 1.0
C
C
C    SUMMER 2017
C    `````````````````````````````````
C
C    versionnum = "2.8"
C    versiondate='August 17, 2017'
C
C *  changes have been made to convertWeak.out, some of which required input format
C    modifications in nucleo; in convertWeak.out:
C
C    - previously, we had set all stellar weak decay rates below T=10 MK equal to 
C      1.00e-99; this leads to discontinuities; therefore, we now set at those low 
C      temperatures the rates equal to the rate at the lowest T grid point of the 
C      old grid [i.e., 10 MK]; this *extrapolation* should be more reliable than 
C      simply setting the decay rates to essentially zero  
C    - another density column is added for rho*Ye=1 g/cm3; adopted values are those 
C      for 10 g/cm3 [this extrapolation should be a reasonable approximation] 
C    - stellar 7Be decay is now implemented [i.e., for IWEAK=0, the 7Be laboratory
C      decay rate will be used by nucleo]
C
C      to summarize: convertWeak.out can now be used for 1 MK < T < 10 GK and
C      1 g/cm3 < rho*Ye < 10^11 g/cm3.
C
C
C======================================================================
       
      Program nucleo 

      include 'comgear.inc'

      DIMENSION X(NSP),xold(nsp),XS(NSP),XST(KP,NSP)
      integer stepstatus
      
      inputfile='nucleo.in'
      open(98,file=inputfile,status='unknown')
      
 9106 format(a)   
      
      ishrun=1			! shell number [for profile]

      CALL START(X)  
C=============== NEXT RUN ===================>

   11 continue
   
      if((inmode.eq.5).and.(ilast.eq.0))then
         WRITE(6,9109) ishrun
      endif
  
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
C            do while(T+delta-tmodel(ncurrtime+1) .gt. -1.d-12)
C               if(DEBUG .eq. 1)
C     &              print*,'Incrementing Time Counter to ',ncurrtime+1
C               ncurrtime=ncurrtime+1
C            end do
C            do while(T+delta-tmodel(ncurrtime) .lt. 1.d-12 .and.
C     &           ncurrtime .gt. 1)
C               if(DEBUG .eq. 1)
C     &              print*,'Decreasing Time Counter to',ncurrtime-1
C               ncurrtime=ncurrtime-1
C            end do
            call locate(tmodel,numod,T+delta,ncurrtime)
            if(DEBUG .eq. 1)
     &           print*,'Time Counter = ',ncurrtime

         else
C            do while(T+delta-tprof(ncurrtime+1) .gt. -1.d-12)
C               if(DEBUG .eq. 1)
C     &              print*,'Incrementing Time Counter to ',ncurrtime+1
C               ncurrtime=ncurrtime+1
C            end do
C            do while(T+delta-tprof(ncurrtime) .lt. 1.d-12 .and.
C     &           ncurrtime .gt. 1)
C               if(DEBUG .eq. 1)
C     &              print*,'Decreasing Time Counter'
C               ncurrtime=ncurrtime-1
C            end do
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
         CALL STOCK(NT,T,SUMX,X,XST)            ! stocking the last values
         CALL SEXIT(NT,T,X,XST)                 ! exit run

c------- exit program for following conditions
         if(inmode.ne.5) then
            stop
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

      STOP
      END
C=======================================================================
      SUBROUTINE START(X)
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

      real*8 dconv
            
      OPEN(9,FILE='nucleo.out',status='unknown')     
      OPEN(21,FILE='totflux.out',STATUS='unknown')
      OPEN(22,FILE='instantflux.out',STATUS='unknown')
      
      read(98,*)
      read(98,*) itest 
      read(98,*) ilast
      read(98,*)
      read(98,*)
      read(98,*) networkfile
      read(98,*) reactionfile
      read(98,*) weakinputfile
      read(98,*)
      read(98,*) inmode
      read(98,*) xx1,xx2,xx3
      read(98,*) xx4,xx5,xx19,xx20,xx21
      read(98,*) xx6,xx7,xx8,xx9a,xx9b
      read(98,*) xx11,xx12,xx16
      read(98,*) xx13,xx14,xx15
      read(98,*)
      read(98,*) nsolver
      read(98,*) eps
      read(98,*) ytmin
      read(98,*) nlast
      read(98,*) xlast,isf
      read(98,*) delta0
      read(98,*) gres
      read(98,*) delfac
      read(98,*)
      read(98,*) iweak
      read(98,*)
      read(98,*) iflux
      read(98,*) isop1,isop2
      read(98,*) (numab(i),i=1,7)
      read(98,*)
      read(98,*) xx10

c --- welcome screen  

      versionnum = "2.8"
      versiondate='August 17, 2017'
    
      if(ilast.eq.0)then 
         call welcomescreen(versionnum,versiondate)
      endif

c --- check some input
      if((itest.ne.0).and.(itest.ne.1))then
         write(6,*) 'Check input: ITEST. Run stop!'
         stop
      endif
      if((ilast.ne.0).and.(ilast.ne.3))then
         write(6,*) 'Check input: ILAST. Run stop!'
         stop
      endif
      if((inmode.lt.1).or.(inmode.gt.5))then
         write(6,*) 'Check input: INMODE. Run stop!'
         stop
      endif
      if((iweak.ne.0).and.(iweak.ne.1))then
         write(6,*) 'Check input: IWEAK. Run stop!'
         stop
      endif
      if((iflux.ne.0).and.(iflux.ne.1))then
         write(6,*) 'Check input: IFLUX. Run stop!'
         stop
      endif
      if((xx10.ne.0).and.(xx10.ne.1))then
         write(6,*) 'Check input: xx10. Run stop!'
         stop
      endif
c -------------------------------------

      nyhandler=1
      
c     check combinations of parameters, disallow certain combinations,
c     e.g., (Monte Carlo) AND (xx10=1) [manual variation of rates]     
      
      if(((itest.eq.1).and.(ilast.eq.3)))then
         write(6,*) 'unacceptable running parameters:'
         write(6,*) 'check ILAST, ITEST etc.'
         write(6,*)
         stop
      endif
      
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

            open(24,file=profilefile,status='unknown')
	        open(26,file='profile.out',status='unknown')
            write(26,9119) 
            
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

C------  parametrized tau-rho (adiabatic expansion) ----------------
         IF(inmode.eq.3) THEN
            TAU0=xx6
	        RHO0=xx7
	        TLAST=xx8
	        scale1=xx9a
            scale2=xx9b
            open(26,file='profile.out',status='unknown')
            write(26,9119)
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
            write(26,9119)
            open(28,file='instmixtest.out',status='unknown')

            write(28,'(a)') 'Decay constants in 1/s'
            write(28,'(a)') '...output may have been commented out...' 
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
      CALL RNETWORK(X)  ! Reading of Network
      
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
         if(ilast.eq.0) then
            write(9,9300)
         endif
         write(9,9952) tau0,rho0
         if(ilast.eq.0) then
            write(9,9300)
         endif
      ENDIF
C------------------------------------
      IF(INMODE.EQ.2) THEN
         CALL PROFIL(T)   ! calculate initial TAU, RHO from profile
         if(ilast.eq.0) then
            write(9,9300)
         endif	 
         WRITE(9,9953) profilefile
         if(ilast.eq.0) then
            write(9,9300)
         endif
      ENDIF
C------------------------------------            
      if(INMODE.eq.3) then
		 tau=tau0
	     rho=rho0
         if(ilast.eq.0) then
            write(9,9300)
         endif
         WRITE(9,9954) TAU0,RHO0,scale1,scale2 
         if(ilast.eq.0) then
            write(9,9300)
         endif
      endif
C------------------------------------
      IF(INMODE.EQ.4) THEN
         call searchmodel(t)  ! calculate initial shell properties
         if(ilast.eq.0) then
            write(9,9300)
         endif	 
         WRITE(9,9957) shellinputfile
         if(ilast.eq.0) then
            write(9,9300)
         endif
      ENDIF
C------------------------------------
      IF(INMODE.EQ.5) THEN
         CALL PROFIL(T)      ! calculate initial TAU, RHO from profile
         write(9,9300)
         WRITE(9,9958) shellinputfile
         write(9,9300)
      ENDIF
C------------------------------------           
      if((ilast.eq.0).and.(inmode.ne.5)) then
         write(9,9962) eps,ytmin,nlast
      endif

      write(9,9955) networkfile
      write(9,9956) reactionfile

      if((ilast.eq.0).and.(inmode.ne.5)) then
         write(9,9322) nis,nr
	     write(9,9300)
      endif
C------------------------------------
      do 676 i=1,nr
         if ((vmult(i).ne.1.d0).and.(ilast.eq.0)) 
     &     write(9,9971) i,vmult(i)
 676  continue

      if((ilast.eq.0).and.(inmode.ne.5)) then 
         write(9,*)
       elseif((ilast.eq.0).and.(inmode.eq.5)) then
         write(9,9300)
      endif

      if((ilast.eq.0).and.(inmode.ne.5)) then
         WRITE(9,9321) NLAST/KP+1
      endif
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
      SUBROUTINE RNETWORK(X)
C---------------------------------------------------------------------
C     READS ISOTOPES AND NETWORK
C======================================================================
      include 'comnuc.inc'
      character*5 xyz(nsl)
      character*7 rspec1(nr),rspec2(nr),rspec3(nr),rspec4(nr)
      character*7 tspec1(nr),tspec2(nr),tspec3(nr),tspec4(nr)
      character*80 dummystr
      DIMENSION X(NSP)
      DIMENSION KNET(0:100,0:150)
      DIMENSION ANO(NSL),ZNO(NSL),XNO(NSL),XINIT(NSL)
      DIMENSION VVX(NR,NKT),TTT(NKT)
      INTEGER n,i,j,k,kk

c     calculates NA<sv> at these T if network is printed
      DATA TTT/0.0030d0,0.05d0,0.100d0,0.15d0,0.20d0,5.0d0,10.0d0/  ! T in 10**9 K

c     ----------------------------------------------------------------------
c     read reaction rates from REACTIONS.DAT; reactions are ordered exactly 
c     the same as in NUCLEO.DAT
c     ----------------------------------------------------------------------
      OPEN(23,FILE=reactionfile,STATUS='unknown')

      do 530 j=1,nr

c        read first line with format and label
         READ(23,9000) REACSTRING(j)
c        reacstring is needed for a number of things
         char_ec(j)=REACSTRING(j)(46:47)  ! is link labeled "ec"? 
         char_w(j)=REACSTRING(j)(48:48)   ! is link a weak or reverse interaction?
         prob(j)=0.d0      
                           
         do 540 i=1,60
c           j: reaction; i: temperature
c           rates of given reaction are multiplied by same factor at all T
            read(23,*) dummy,rate(j,i),fu(j,i)
c           if rate=0.0: set equal to 1.d-99 for logarithmic interpolation
            if(rate(j,i).eq.0.d0) then
               rate(j,i)=1.0d-99    
            endif   
 540     continue

 530  continue
 
c     this is the "flat" parameterization: for each run, a probability
c     factor p is sampled for each reaction; for a given reaction,
c     p has the same value at all temperatures;
c     for INMODE=5, sample only between runs, not between shells; randomly
c     sample all values of p, including for reverse reactions; will be taken
c     care of later
      if(ilast.eq.3)then        
         open(36,file='prob.dat',status='unknown')
         do 206 j=1,nr
            read(36,1011) prob(j) 
 206     continue
         close(36)
      endif
      
 1011 format(7x,1pe20.12)
 
      close(23)
    
c     ----------------------------------------------------------------
c     read weak interaction rates if needed
c     ----------------------------------------------------------------
      if(iweak.eq.1)then

         open(27,file=weakinputfile,status='unknown')
         j=1
c        j: reaction; i: temperature; n: density      
 62      read(27,9001,end=60) reacstrweak(j)
         do 63 i=1,60
            read(27,9002) dummy,(ratew(j,i,n),n=1,12)
            do 64 n=1,12
               ratelogw(j,i,n)=dlog10(ratew(j,i,n))
 64         continue
 63      continue
         j=j+1
         goto 62
 60      continue
         nrw=j-1		! nrw=number of weak rates read
c        write(6,*) nrw
         close(27)
      endif
c     ----------------------------------------------------------------
c     read isotopes and network from NUCLEO.DAT
c     ----------------------------------------------------------------
      OPEN(8,FILE=networkfile,status='unknown') 

c---- nuclides 

      READ(8,8101) ABCD
      READ(8,8100) ABCD
      DO 40 I=1,NIS
   40 READ(8,9025) K,ON(I),AN(I),ZN(I)
      ON(NIS+1)='    '
      AN(NIS+1)=1.d0
      ZN(NIS+1)=0.d0

c---- initial abundances 

      READ(8,8101) ABCD
      DO 80 I=1,NSL
         READ(8,8050) xyz(i),ANO(I),ZNO(I),XINIT(I)
c        calculate initial eta and Ye
         ETA=ETA+(ANO(I)-2.d0*ZNO(I))*(XINIT(I)/ANO(I))
   80 CONTINUE
      YE=0.5d0*(1.d0-ETA)
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
   
      if((ilast.eq.0).and.(inmode.ne.5))then  
         write(6,*) 'X_i_sum=',sumx
      endif
      
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
      READ(8,8101) ABCD
      DO 93 L=1,NR
         READ(8,8060) NNNR,NA,SPE1,NB,SPE2,ND,SPE4,NC,SPE3,Q(L)
         IA=ISPE(SPE1)
         IB=ISPE(SPE2)
         IC=ISPE(SPE3)
         ID=ISPE(SPE4)
         K1(L)=IA
         K2(L)=NA
         K3(L)=IB
         K4(L)=NB
         K5(L)=IC
         K6(L)=NC
         K7(L)=ID
         K8(L)=ND
         AA1=NA*AN(K1(L)) + NB*AN(K3(L))
         AA2=NC*AN(K5(L)) + ND*AN(K7(L))
         if(AA1.NE.AA2) then
            write(6,*) L
            write(6,*) 'WARNING: BARYON NUMBER NOT CONSERVED!'
            stop
         endif

c        store identity of emitted light particle for each reaction;
c        useful for identifying links that give rise to, e.g., neutron
c        emission, etc.

         emitspe(l)=spe4

   93 CONTINUE
      close(8)

c     ------------------------------------------------------------------     
c     reverse rates must be multiplied by the same probability factor p;
c     thus find corresponding forward rate for given reverse rate
c     ------------------------------------------------------------------     

      if(ilast.eq.3)then

         OPEN(8,FILE=networkfile,status='old')
         do 533 j=1,2+nis+1+nsl+1
            read(8,8101) dummystr
 533     continue
         do 534 j=1,nr
            read(8,8101) dummystr
            rspec1(j)=dummystr(10:16)
            rspec2(j)=dummystr(20:26)
            rspec3(j)=dummystr(29:35)
            rspec4(j)=dummystr(39:45)
            
            tspec1(j)=rspec1(j)
            tspec2(j)=rspec2(j)
            tspec3(j)=rspec3(j)
            tspec4(j)=rspec4(j)
 534     continue
         close(8)

c        find corresponding forward and reverse rates
         j=0
 699     j=j+1
         if(char_w(j).eq.'v')then
            kk=0
 698        kk=kk+1
            if((rspec1(j).eq.tspec4(kk)).and.
     &         (rspec2(j).eq.tspec3(kk)).and.
     &         (rspec3(j).eq.tspec2(kk)).and.
     &         (rspec4(j).eq.tspec1(kk)))then
               prob(j)=prob(kk)
               if(j.lt.nr)then 
                  goto 699
               else
                  goto 697
               endif
            else
               if(kk.lt.nr)then 
                  goto 698
               else
                  goto 694
               endif
            endif
 694        write(6,*) ' Forward reaction not found for:'
            write(6,9003) rspec1(j),rspec2(j),rspec3(j),rspec4(j)
            stop
         else
            if(j.lt.nr)then 
               goto 699
            else
               goto 697
            endif
         endif
 697     continue
      endif    

c     ----------------------------------------------------------------
c     write new prob.dat file
c     ----------------------------------------------------------------
      if(ilast.eq.3)then        
         open(37,file='prob.dat',status='unknown')
         do 209 j=1,nr
            write(37,1012) j,prob(j)
c           to account for unknown systematic uncertainty factor, vmult
            vmult(j)=vmult(j)**prob(j)    
 209     continue
         close(37)
      endif
 1012 format(1x,i5,'=',1pe20.12)
c     ----------------------------------------------------------------
c     compute rates
c     ----------------------------------------------------------------
      do 531 j=1,nr
         do 532 i=1,60
            if(ilast.eq.3) then
                rate(j,i)=rate(j,i)*(fu(j,i)**prob(j))
            endif
            ratelog(j,i)=dlog10(rate(j,i))
 532     continue
 531  continue
      
c----------------------------------------------------------------------
c---- print network ---------------------------------------------------
c----------------------------------------------------------------------
      if(itest.eq.1)then
        write(9,*) 'Network taken from ',networkfile
        write(9,*) 'Rates taken from file ',reactionfile
        write(9,*)
        do 696 i=1,nr
          if ((vmult(i).ne.1.d0).and.(ilast.eq.0)) 
     &      write(9,315) i,vmult(i)
 696    continue
 315    format(' ',' Reaction ',i4,' changed by a factor of ',1PD10.3)
        write(9,*)
      endif

      DO 103 I=1,NIS
         IF(ITEST.EQ.1) WRITE(9,9965) I,ON(I),ZN(I),AN(I),X(I)
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
         write(9,*)
         if(iweak.eq.1)then
            write(9,8899) rho*YE
         endif
         write(9,*)
         WRITE(9,9810) (TTT(K),K=1,NKT)
         DO 118 L=1,NR
            WRITE(9,9071) L,K2(L),ON(K1(L)),K4(L),ON(K3(L)),
     &      AAAAA,K6(L),ON(K5(L)),K8(L),ON(K7(L)),BBBBB,Q(L),
     &      (VVX(L,K),K=1,NKT)
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

      if((ilast.eq.0).and.(inmode.ne.5)) then
         WRITE(9,9020) ! Intermediate printing of some results
         WRITE(9,9200) T,T/SECY,DELTA,NSTEP
     &     ,TAU,LOG10(TAU),RHO,SUMX,ETA,YE
         
         if(iflux.eq.1) then
	        WRITE(22,9200) T,T/SECY,DELTA,NSTEP
     &        ,TAU,LOG10(TAU),RHO,SUMX
            WRITE(22,9993)
            DO 1324 I=1,NR
              WRITE(22,1325) I,FLUX(I)
 1324       continue
         endif

         WRITE(9,9995)
      endif

      if((inmode.ne.5).and.(ilast.eq.0))then  
         write(6,9744) t,on(isop1),x(isop1),on(isop2),x(isop2),
     &     sumx,eta,tau
      endif

      SUMENER=0.D0
      DO 433 I=1,NR
  433 SUMENER=SUMENER+ENER(I)

      if((ilast.eq.0).and.(inmode.ne.5)) then
         WRITE(9,2223) SUMENER
         WRITE(9,9997)
         WRITE(9,9998) (ON(I),X(I),I=1,NIS)
         WRITE(9,9020)
      endif

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

      if(inmode.ne.5) then
         WRITE(9,9270)
      endif

      if((inmode.ne.5).or.(ilast.ne.3))then
         WRITE(9,9271) NSTEP,NLAST,T,TLAST,ON(ISF),X(ISF),XLAST
      endif
         
      if(inmode.eq.2) then
         deltaprof=(tauprofmax-tauprocalmax)/tauprofmax
         if(deltaprof.le.1.d-2) then
            write(9,9007) 'PROFILE CHECK: PASS WITH',deltaprof
         else         
            write(6,9005) 'PROFILE CHECK: WARNING',deltaprof,tauprofmax,
     &          tauprocalmax
            write(9,9005) 'PROFILE CHECK: WARNING',deltaprof,tauprofmax,
     &          tauprocalmax            
         endif
      endif

      if(inmode.eq.4) then
         write(9,*) 'PROFILE CHECK: none'
      endif

      if(inmode.eq.5)then
         trun(ishrun)=t   ! store all running times; they need to be the same for averaging x
         do 777 i=1,nsp
            xshells(i,ishrun)=x(i)
 777     continue
      endif

      if(inmode.eq.5)then
         deltaprof=(tmaxsh(ishrun)-tmaxcode(ishrun))/tmaxsh(ishrun)
         if(deltaprof.ge.1.d-2)then
            write(6,9005) 'PROFILE CHECK: WARNING',deltaprof,
     &          tmaxsh(ishrun),tmaxcode(ishrun)
         endif
      endif      

      if((inmode.eq.5).and.(ilast.eq.0))then
c        ishrun labels the shells for sequential runs if INMODE=5
         write(9,9305) ishrun         
      endif

      if((inmode.ne.5).or.(ilast.ne.3))then
         if(inmode.ne.5)then
             WRITE(9,9343) 
           elseif(inmode.eq.5)then
             WRITE(9,9345) 
         endif 
         write(9,9344) TOTENERGY,ETA,YE
         if(inmode.ne.5)then
             WRITE(9,9997)
           elseif(inmode.eq.5)then
             WRITE(9,9991)
         endif 
         WRITE(9,9999) (ON(I),X(I),I=1,NSP)
      endif
      
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

      if((inmode.ne.5).and.(ilast.eq.0)) then
         WRITE(9,9998)
         WRITE(9,9999) (ON(I),XS(I),I=1,NSP)
         write(9,9432)
         WRITE(9,9996) (I,tenergy(I),I=1,NR)
         write(9,9221)
         WRITE(9,9220) (I,100*tenergy(I)/TOTENERGY,I=1,NR)

         WRITE(9,9992) (ON(numab(i)),i=1,7)
         DO  486 I=1,NT
           WRITE(9,9988) NNSTEP(I),TIME(I),TEMM(I),DENN(I),
     &       (XST(I,numab(k)),k=1,7)
  486    CONTINUE
         WRITE(9,9270)
      endif

      if(ilast.eq.0)then
         write(6,9376)
      endif
      
      if(inmode.ne.5) then 
         DO 1324 I=1,NR
           if(totflux(i).le.1.d-99) then          
              totflux(i)=0.d0
           endif
           WRITE(21,1325) I, TOTFLUX(I)
 1324    continue
      endif            
            
      close(21)
      close(22)
      close(24)
      close(25)
      close(28)

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
      FUNCTION ISPE(NAME)
C======================================================================
      include 'comnuc.inc'

      DO 1 I=1,NSP
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
         write(9,9996) (JP(i), i = lk1, lk2)
         do 200 i = lk1, lk2
         write(9,9998) i, (IT(j,i), j = lk1, lk2)
  200    continue
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
            write(26,5068) t,tau,rho
c            write(26,5065) nstep,t,tau,rho
cc		 endif
      endif

c---- calculation from analytical function     
      if(inmode.eq.3) then
         tau=tau0*dexp(-t/scale1)
         rho=rho0*dexp(-t/scale2)
         if(ilast.eq.0) then	 
	        write(26,5065) nstep,t,tau,rho
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
         write(26,5066) ishrun,nstep,t,tau,rho
         
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

      if(ilast.eq.0)then
         write(9,9300)
      endif
      
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
      write(9,9300)

      WRITE(9,9271) tlast,xavsum,dmst
      write(9,9272) 
      write(9,9273) totenerav,eta,ye

      WRITE(9,9997)
      WRITE(9,9999) (ON(I),XAV(I),I=1,NSP)

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

      if(ilast.eq.0)then
         write(9,*)
         WRITE(9,9998)
         WRITE(9,9999) (ON(I),XS(I),I=1,NSP)
      endif
      
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
