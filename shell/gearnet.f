C------------------------------------------------------------
C     Subroutine to initialise the MA48 routines etc.
      subroutine initma48(x)

      include 'comgear.inc'

      dimension x(nsp),y(nsp)

C     Get the sparse matrix pattern
      call getjacpat

      do 5 i=1,nis
         y(i) = x(i)/an(i)
 5    continue

c     Make the jacobian form the starting parameters. The numbers will
c     change througout the calculation, but the matrix pattern will not.
      call makejac(delta,y)

C     Set default controls
      call ma48id(cntl,icntl)

         

C     Drop tolerances
c      cntl(3)=dymin
c      cntl(4)=dymin

C     Find pivot order, etc. Restrict to the diagonal by setting job=3
      job = 3
      call ma48ad(nis,nis,ne,job,la,amat,irow,jcol,keep,cntl,icntl,
     $     iw,info,rinfo)
      if (info(1).lt.0) then
        write (6,'(a,i3)') 'error stop from ma48a/ad with info(1) ='
     +    ,info(1)
        stop
      end if
c     Print some info about the matrix call
c      print*,'Done with ma48ad'
c      write(6,*)(info(i),i=1,12)
c      print*,(irow(i),jcol(i),i=1,ne)

C     Factorize matrix initially. This should speed it up later
      job=1
c      icntl(10)=1
      call ma48bd(nis,nis,ne,job,la,amat,irow,jcol,keep,cntl,icntl,
     $     w,iw,info,rinfo)
      if (info(1).ne.0) then
        write (6,fmt='(a,i3/a)') 'stop from ma48b/bd with info(1) =',
     +    info(1),'solution not possible'
        stop
      end if
c     Print some info about the matrix call
c      print*,'Done with ma48bd'
c      write(6,*)(info(i),i=1,12)

      if(nsolver .eq. 1) then
C        Initialize the Nordsieck array. 
         do i=1,nis
            vnord(i,1)=X(i)/AN(I)
         enddo
         do j=2,7
            ptau(j)=0.d0        ! Previous values for h
            do i=1,NIS
               vnord(i,j)=0.d0
            enddo
         enddo
         nordord = MINORD
         
         do i=1,7
            ptau(i)=0.d0
            el(i)=0.d0
         enddo
         el(8)=0.d0
         nordsteps=0
         nrestep=0
      endif  


C     Finished setting it all up. Now return to the burner and start
C     doing some actual nucleosynthesis!
      return
      end

C---------------------------------------------------------------
C     Subroutine to calculate the right hand side of the evolution
C     equation:
C               (I/h - J)Delta = dydt
      subroutine rhs(y,dydt,dt)

      include 'comgear.inc'
      dimension y(nsp),dydt(nsp)
      integer r,IA,IAA,A,AA,IB,IBB,B,BB

      O=1.d-100

      do 100 i=1,nis
         dydt(i) = 0.d0
 100  continue

C     For each reaction, add or subtract the dydt for each part
C     Note: For each reaction, the nuclei are labelled like:
C              iA(ja,kb)lB
C              k1 is index of A, k2 is number involved
C              k3    index of a, k4 "   "
C              k5             b, k6 "   "
C              k7             B, k8 "   "
      do 200 r=1,nr
         IA=k2(r)
         A=k1(r)
         IAA=k4(r)
         AA=k3(r)
         IB=k6(r)
         B=k5(r)
         IBB=k8(r)
         BB=k7(r)

         AFactor = -v(r)/(F(IA)*F(IAA))
c         AFactor = -v(r)/(F(IA)*F(IAA)*(IA+IAA))
c         if(Y(A) .lt. YMIN)Y(A)=0.d0
c         if(Y(AA) .lt. YMIN)Y(AA)=0.d0

         dydt(A) = dydt(A)   + AFactor*IA*((Y(A))**(IA))*
     $        ((Y(AA))**IAA)
         if(IAA .gt. 0)
     $        dydt(AA) = dydt(AA) + AFactor*IAA*((Y(A))**(IA))*
     $        ((Y(AA))**IAA)
         dydt(B) = dydt(B)   - AFactor*IB*((Y(A))**(IA))*
     $        ((Y(AA))**IAA)
         if(IBB .gt. 0)
     $        dydt(BB) = dydt(BB) - AFactor*IBB*((Y(A))**(IA))*
     $        ((Y(AA))**IAA)

 200  continue
      
C     Multiply by the step length
      do 300 i=1,nis
c         if((dydt(i) .lt. dymin) .and. (dydt(i) .gt. 0.d0))
c     $        dydt(i) = dymin
c         if((dabs(dydt(i)) .lt. dymin) .and. (dydt(i) .lt. 0.d0))
c     $        dydt(i) = -dymin

         dydt(i) = dydt(i)*dt
c         if(dabs(dydt(i)) .lt. dymin)dydt(i)=0.d0
c         print*,dydt(i)
 300  continue

      return
      end

C======================================================================
C     Subroutine to go through all reaction rates and pack them into the
C     Jacobian. Then multiply by h, subtract from unity matrix to get
C     the A matrix in Timmes, Astro. J. Supp. 124 (1999) Page 12
      subroutine getjacpat

      include 'comgear.inc'

      integer pat(2,8),nuc(4),count,kmax,reac
      logical diag

C     Maximum ne can be is unique entries for every reaction. In reality
C     this will be a lot smaller, which we intend to find here!
      ne=nr*8

C     Zero everything
      do 100 i=1,la
         amat(i)=0.d0
         irow(i)=0
         jcol(i)=0
 100  continue

c     pattern
C        Each reaction has 6 entries in the jacobian
C           (A,A),(A,a),(a,A),(a,a),(B,A),(B,a),(b,A),(b,a)
      pat = reshape((/ 1,1,1,2,2,1,2,2,3,1,3,2,4,1,4,2 /),shape(pat))

      count=0
      do 200 reac=1,nr
C        Note: For each reaction, the nuclei are labelled like:
C              iA(ja,kb)lB
C              k1 is index of A, k2 is number involved
C              k3    index of a, k4 "   "
C              k5             B, k6 "   "
C              k7             b, k8 "   "
         nuc(1)=k1(reac)
         nuc(2)=k3(reac)
         nuc(3)=k5(reac)
         nuc(4)=k7(reac)

C        Does this reaction produce gammas? If so, ignore them
         kmax=8
c         if(k8(reac) .eq. 0)kmax=6

         do k=1,kmax
c           only fill matrix if it's not a reaction with gamma
            if((nuc(pat(1,k)) .ne. nsp) .and. 
     $           (nuc(pat(2,k)) .ne. nsp)) then
               count = count+1
               ent(k,reac)=count
               
C              Check if J(x,y) is already filled.
               do i=1,max0(0,count-1)
                  if((irow(i) .eq. nuc(pat(1,k)))
     $                 .and. (jcol(i) .eq. nuc(pat(2,k)))) then
                     ent(k,reac)=i
                     count = count-1
                     exit
                  endif
               enddo
C              Assign the jacobian element positions. These will be
C              re-assigned with idential numbers sometimes, but that's
C              ok...
               irow(ent(k,reac))=nuc(pat(1,k))
               jcol(ent(k,reac))=nuc(pat(2,k))
                  
            endif
         enddo

 200  continue

C     So we now have, for each reaction, the 6 jacobian elements that it
C     affects stored in the ent(k,nr) array. An example of how to access
C     it is (psuedo code):
C     Entry = ent(5,reac)
C     A(Entry) += v(reac)Y
C     irow(Entry) = nuc(5,reac)
C     jcol(Entry) = nuc(5,reac)


C     The diagonal elements must also be added if not there
      do i=1,nis
         diag = .false.
         do j=1,count
            if((irow(j) .eq. jcol(j)) .and. (irow(j) .eq. i))then
               diag = .true.
C              Let the routines know which entry the diagonal is
               diagindex(i)=j
            endif
         enddo
         if(.not. diag)then  ! if not assigned, make the diagonal point
            count = count+1
            irow(count) = i
            jcol(count) = i
            diagindex(i)=count
         endif
         
      enddo

C     Also know the number of jacobian elements now!
      ne = count
c      print*,'diagindex(1) = ',diagindex(1)
      return
      end

C======================================================================
C     Subroutine to go through all reaction rates and pack them into the
C     Jacobian. Then multiply by h, subtract from unity matrix to get
C     the A matrix in Timmes, Astro. J. Supp. 124 (1999) Page 12
      subroutine makejac(dt,y)

      include 'comgear.inc'

      integer reac,j,IA,IAA,A,AA,IB,IBB,B,BB
      dimension y(nsp)
      
      O=1.D-100
C     Zero everything
      do 100 i=1,ne
         amat(i)=0.d0
 100  continue

      do 200 reac=1,nr
C        Note: For each reaction, the nuclei are labelled like:
C              iA(ja,kb)lB
C              k1 is index of A, k2 is number involved
C              k3    index of a, k4 "   "
C              k5             b, k6 "   "
C              k7             B, k8 "   "
         
         IA=k2(reac)    ! number of A in reaction
         A=k1(reac)
         IAA=k4(reac)   ! number of a     "
         AA=k3(reac)
         IB=k6(reac)    ! number of B     "
         B=k5(reac)
         IBB=k8(reac)   ! number of b     "
         BB=k7(reac)

         if(nsolver .eq. 0)then
            AFactor = v(reac)*dt/((F(IA)*F(IAA)*(IA+IAA)))
         else 
            AFactor = v(reac)*dt/((F(IA)*F(IAA)))
         endif
c         if(dabs(AFactor) .gt. 1.d-100) then
C           (A,A),(A,a),(a,A),(a,a),(B,A),(B,a),(b,A),(b,a)
            
            i=ent(1,reac)       !(A,A)
            amat(i) = amat(i)+AFactor*IA*IA*((Y(A))**(IA-1))*
     $           ((Y(AA))**(IAA))

            i=ent(5,reac)       !(B,A)
            amat(i) = amat(i)-(AFactor*IB*IA*((Y(A))**(IA-1))*
     $           ((Y(AA))**(IAA)))

            if(IBB .gt. 0)then
               i=ent(7,reac)    !(b,A)
               amat(i) = amat(i)-(AFactor*IBB*IA*((Y(A))**(IA-1))*
     $              ((Y(AA))**(IAA)))

            endif
            if(IAA .gt. 0)then
               i=ent(2,reac)    !(A,a)
               amat(i) = amat(i)+(AFactor*IAA*IA*
     $              ((Y(AA))**(IAA-1))*((Y(A))**(IA)))

               i=ent(3,reac)    !(a,A)
               amat(i) = amat(i)+(AFactor*IAA*IA*
     $              ((Y(A))**(IA-1))*((Y(AA))**(IAA)))

               i=ent(4,reac)    !(a,a)
               amat(i) = amat(i)+(AFactor*IAA*IAA*
     $              ((Y(AA))**(IAA-1))*((Y(A))**(IA)))

               i=ent(6,reac)    !(B,a)
               amat(i) = amat(i)-(AFactor*IAA*IB*
     $              ((Y(AA))**(IAA-1))*((Y(A))**(IA)))

               if(IBB .gt. 0)then
                  i=ent(8,reac) !(b,a)
                  amat(i) = amat(i)-(AFactor*IAA*IBB*
     $                 ((Y(AA))**(IAA-1))*((Y(A))**(IA)))
               endif
            endif
c         endif

 200  continue


C     Add the unity matrix
      do 300 i=1,nis
         j=diagindex(i)
         amat(j) = amat(j)+1.d0
 300  continue


      return
      end

      SUBROUTINE gearstep(eatascale,X,SUMX,DTOLD,T,stepstatus)
C----------------------------------------------------------------------
C     NUCLEOSYNTHESIS FOR TIME-STEP DELTA
C----------------------------------------------------------------------
      include 'comgear.inc'
 
c      PARAMETER(eatascale=0.5d0)
c      PARAMETER(eatascale=0.25d0,thresh=1.0d-3)
      PARAMETER(thresh=1.0d-3)

      DIMENSION Y0(NSP),YB(NSP),YT(NSP),X(NSP)
      dimension dydt(nsp),dydtb(nsp),deltay(nsp),
     $     rh(nsp),yscal(nsp)
      integer step,changedir,negy,k,stepstatus
      double precision T,stepError(3),crate(nsp)
      double precision ar(nsp)
      double precision astime,rstime,aftime,rftime
      logical takeanotherstep,ordchange,trans
            
      O=1.D-100
      DELT=DELTA
c      print*,"delt = ",delt
      flotnis = dble(nis)
      ordchange=.false.    ! Control if order change will be attempted

C     Error scale to check errors
c     errorscale = EPS*delt/tlast
c      errorscale = EPS*EPS          !*delt
      errorscale = EPS

C     Don't increment order step counter unless it's a real step
      if(stepstatus .eq. 0)then
c        Count the steps taken at the current order
         nordsteps = nordsteps+1
         nmakejac = nmakejac-1
         nrestep=0
      endif

      stepstatus=0

      if(nordsteps > nordord+1)ordchange=.true.

   3  DO 5 I=1,NSP       ! Assign Y values...
c        convert X into number abundance
         Y0(I)=X(I)/AN(I)
         YB(I)=0.d0
         YT(I)=0.d0
         yscal(i)=max(y0(i),ytmin)
c         YS(I)=0.d0
   5  CONTINUE

C     in order to avoid that NA<sv> is calculated at each time step
C     for const. T,rho
C     -----------------------------------------------------
C     This is where to return to if a step fails
 10   CONTINUE

C     make the right hand side (eqn 2.23 in Byrne)
C     Note that dydt is actually h*dydt
c      call rhs(y0,dydtb,delt)
c      do i=1,nis
c         vnord(i,2)=dydtb(i)       ! The second column: rates of change
c      enddo

      takeanotherstep=.false.   ! Control if we need to try again

C     Set the first 2 columns of the Nordsiek array by hand. The rest is
C     automatically updated as we go
      IF(DEBUG .eq. 1)then
         print*,"The initial y and dydt values are"
         do i=1,10
            WRITE(6,'(i4,7es11.3)')i,(vnord(i,j),j=1,7)
         enddo
         do i=nis-10,nis
            WRITE(6,'(i4,7es11.3)')i,(vnord(i,j),j=1,7)
         enddo
      endif

C     ---- PREDICTION ----
C     Do prediction by multiplying the nordsiek array by the pascal
C     matrix
      call gearpredict
      IF(DEBUG .eq. 1)
     $     print*,"PREDICTRED Y and DYDT VALUES"

C     Set up the abundances to evolve:
C        y0 = start value
C        yt = predicted value
C        yb = corrected value
      do i=1,nis 
         yt(i) = vnord(i,1)
         dydt(i) = vnord(i,2)
         yb(i) = yt(i)          ! yt continues to be the predictor
         dydtb(i) = dydt(i)     ! while the yb and dydtb are corrected
         deltay(i)=0.d0
c         crate(i)=1.d0
         crate(i)=0.d0
      enddo         
      IF(DEBUG .eq. 1)then
         do i=1,10
            WRITE(6,'(i4,7es11.3)')i,(vnord(i,j),j=1,7)
         enddo
         do i=nis-10,nis
            WRITE(6,'(i4,7es11.3)')i,(vnord(i,j),j=1,7)
         enddo
      endif

C     ---- ----
C     And then calculate the el parameter
      call calcel(delt,ordchange)
      rel1=1.d0/el(2)

C     Only make jacobian on first pass or if the step failed more than 3
C     times
      if((nrestep .eq. 0) .or. (nrestep .ge. 2))then
         nmakejac = jacrate     ! reset the counter 
C        Pack the rates into
C        the (h/l1)*jacobian+I = A matrix and factorise it
         deltrel1 = delt*rel1
         call cpu_time(astime)
         call makejac(deltrel1,yt)
         call cpu_time(aftime)
         ratetime=ratetime+(aftime-astime)

         call cpu_time(astime)
         job=1   ! 2=fast factorisation of Jacobian on second step since
                 ! we've already done the hard part in the initialisation
         call ma48bd(nis,nis,ne,job,la,amat,irow,jcol,keep,cntl,icntl,
     $        w,iw,info,rinfo)
         njacb=njacb+1          ! count one jacobian evaluation
         call cpu_time(aftime)
         algebratime=algebratime+(aftime-astime)

      endif

C     ---- CORRECTOR STEPS ----
C     No do the corrector steps as many times as neccessary (up to 3)
      do step=1,4
         ncor=ncor+1
         negy=0 
         IF(DEBUG .eq. 1)
     $        write(6,'("Corrector Step: ",i3)')step

C        make the right hand side (eqn 2.23 in Byrne)
C        Note that dydt is actually h*dydt
         call cpu_time(astime)
         call rhs(yb,dydtb,delt)
         call cpu_time(aftime)
         algebratime=algebratime+(aftime-astime)

         nfunc=nfunc+1
         do i=1,nis
            rh(i) = -(yb(i)-yt(i)) + rel1*(dydtb(i) - dydt(i))
         enddo
      
C        Solve the matrix equation to get DeltaY
         trans = .false.        ! don't solve the transpose
         job=1                  ! no iterative refinement
         call cpu_time(astime)
         call ma48cd(nis,nis,trans,job,la,amat,irow,keep,cntl,icntl,
     $        rh,deltay,error,w,iw,info)
         call cpu_time(aftime)
         algebratime=algebratime+(aftime-astime)
         njace=njace+1

C        Corrected Y(m+1) = Y(m)+deltaY  (m=0,1,...)
         err=0.d0
         do i=1,nis
            yb(i) = yb(i) + deltay(i) ! Actually correct Y
c            if(yb(i) .lt. 0.d0-dymin) then
            if(yb(i) .lt. 0.d0-dymin)then
               select case (nyhandler)
               case (0)  ! ignore
                  continue
               case (1)  ! don't allow
                  negy=1        ! negatives are baaaad
               case (2)  ! fix
                  yb(i) = 0.d0
               end select
            endif

C           Calculate the convergence rate
c            ar(i) = (yb(i)-yt(i))/yb(i)
c            crate(i) = yb(i)*deltay(i)/(yscal(i)*crate(i))
            crate(i) = deltay(i)/yscal(i) !- crate(i)
c            if(yb(i) .lt. ytmin .or. dabs(deltay(i)).lt.dymin)then
c               crate(i) = 1.d0
c            else
            err = err + crate(i)*crate(i)
c            endif

         enddo
         err = sqrt(err)
         IF(DEBUG .eq. 1)then
            print*,"Yb, DeltaY, ar, and crate:"
            do i=1,10
               write(6,'(i4,4es11.3)')i,yb(i),deltay(i),ar(i),crate(i)
            enddo
            do i=nis-10,nis
               write(6,'(i4,4es11.3)')i,yb(i),deltay(i),ar(i),crate(i)
            enddo
C           Convergence?
            write(6,'("Error and scaled eps: ",2es12.3)')err,errorscale
         endif

C        Exit this loop if we have good convergence
         if(((negy .eq. 0) .and. (err .lt. errorscale)) .and. 
     $        (step .gt. 1)) goto 100

      enddo !step=1,3

C     ---- END OF CORRECTOR STEPS ----
C     Check. If code made it here, no convergence occured
C     and we need to reduce step to try again.
      takeanotherstep=.true.
      if(debug .eq. 1)then
         if(negy .eq. 1)print*,'There was a negative abundance'
         print*,'Did not converge. Taking another step'
      endif
         
C     ---- CORRECT THE NORDSIEK VECTOR (vnord) ----
C     Calculate the correction amount (e_n in Byrne) and apply it to the
C     Nordsiek array (See Eq. 2.19 in Byrne)
 100  maxcorr=0.d0
      do i=1,nis
c         corr(i) = 0.d0
c         if(dabs(yb(i)) .gt. ytmin)
         corr(i)=yb(i)-yt(i)
      enddo
      call calcerrors(ordchange,stepError,yt,yb)

C     Change stepsize to retry step if needed
 101  if((stepError(1) .gt. errorscale) .or. takeanotherstep)then
         nfail=nfail+1
         nrestep=nrestep+1
         if(nrestep .ge. 10)goto 555

C        Correct the nordsiek vector to before the predictor
         do j1 = 1,nordord
            do j2 = j1,nordord
               j = (nordord + j1) - j2
               do i = 1,nis
                  vnord(i,j) = vnord(i,j) - vnord(i,j+1)
               enddo
            enddo
 230     enddo

C        Have we hit a tough spot? If so, change the order
         if(nrestep .eq. 3)then
            if(DEBUG .eq. 1)then
               WRITE(6,'("Changing order from",i3," to 1")')nordord
            endif
            k=nordord
c            print*,k
            do i=k,2,-1         ! reduce the order down to 2
               call reduceord(delt)
               nordord=nordord-1
            enddo
            nordord=1
            nordsteps=1
            ordchange=.false.
            eata=1.d0
         elseif((nrestep .eq. 1) .and. takeanotherstep) then
C           Sometimes we just need to try again
            goto 10
         else                   ! if not, try a smaller step
C           Calculate the new step size to try
            pow=1.d0/(dble(nordord)+1.d0)
            eata=0.5d0*(errorscale/stepError(1))**pow
            eata = dmax1(dmin1(eata,5.0d-1),1.0d-4)
c            eata = dmax1(eata,1.0d-4)
            delt = delt*eata
C           Then use it to correct the Nordsiek Vector
            do i=1,nis
               do j=2,nordord+1
                  vnord(i,j)=vnord(i,j)*eata**(j-1)
               enddo
               vnord(i,1)=y0(i)
            enddo
         endif

         if(debug .eq. 1)then
            write(6,'("Step",i4," Failed",i4," Times!")')nstep,nrestep
            write(6,'("Correcting time step by eata =",es12.3)')eata
            write(6,'("Scaled EPS = ",es9.1," stepError = ",
     $           es11.3)')errorscale,stepError(1)
            write(6,'("New time step = ",es12.3)')delt
            write(6,'("New order = ",i4)')nordord
         endif

C        Go back and try again
c         goto 10
         dtold=delt
         delta=delt
         stepstatus = -1
         return
      endif

C     If error tests are passed, finally correct the Nordsiek vector
      do i=1,nis
         do j=1,nordord+1
            vnord(i,j)=vnord(i,j)+corr(i)*el(j)
c            if(dabs(vnord(i,j)) .lt. dymin)vnord(i,j)=0.d0
         enddo
c         if(dabs(vnord(i,1)) .lt. dymin)vnord(i,1)=0.d0
         pcorr(i)=corr(i)    ! Save the previous values
C        And the total errors
         TotErr(i)=TotErr(i)+Efact(1)*corr(i)
      enddo

      if(DEBUG .eq. 1) then
         print*,'Final Nordsiek vector'
         testsum=-1.d0
         do i=1,nis
            testsum = testsum+an(i)*yb(i)
         enddo
         do i=1,10
            WRITE(6,'(i4,7es12.3)')i,(vnord(i,j),j=1,7)
         enddo
         do i=nis-10,nis
            WRITE(6,'(i4,7es12.3)')i,(vnord(i,j),j=1,7)
         enddo                  
         print*,'X-1=',testsum
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Update X and calculate SUMX and eta
      SUMX=0.d0
      ETA=0.d0
      DO 360 I=1,NSP
         Y0(I)=YB(I)          ! New Y.....
         X(I)=AN(I)*Y0(I)     ! .. and X values
         ETA=ETA+(AN(I)-2.*ZN(I))*Y0(I)
         SUMX=SUMX+X(I)*1.d0
  360 CONTINUE
      YE=0.5d0*(1-ETA)

      SE=SUMX-1.d0
C     Check of SumX conservation (Epsilon~1.d-6 OK)
      IF(dABS(SE).GT.EPS)GO TO 430 

C     If we're at a good order change point, see if it's worthwhile
      eata = 0.d0
      changedir=0
      if(ordchange)then
C        nordsteps=0
         if(nordord .lt. maxord)then   ! not alowed above q=5
            changedir=1
            pow=0.5d0/(dble(nordord)+2.d0)
            eata= 1.d0*eatascale*(errorscale/stepError(3))**pow
         endif
         if(nordord .gt. minord)then   ! not allowed below q=1
            pow=0.5d0/(dble(nordord))
            eatat = 1.d0*eatascale*(errorscale/stepError(2))**pow
            if(DEBUG .EQ. 1)
     $           write(6,'(a,es12.3,a,es12.3)')"Reduction eata",eatat,
     $           " Increase eata",eata
            if((eatat-eata) .gt. thresh)then
               eata=eatat
               changedir=-1
            endif
         endif
      endif

C     Now calculate the new time step, and see if an order change is
C     worth it.
      pow=0.5d0/(dble(nordord+1.d0))
      eatat=eatascale*(1.d0*errorscale/stepError(1))**pow
      if(DEBUG .EQ. 1)
     $           write(6,'(a,es12.3)')"Constant eata",eatat
c      eata=eatat
c      if((eatat- eata)/eatat .gt. -1.d0)then   
      if((eatat- eata) .gt. thresh)then   
c        now, an order change is not worth it
         eata = eatat
         changedir=0
      endif

      DTOLD=DELT      ! The time step that was actually used is DELT

C     Call the step limitor, which applies all kinds of checks on the
C     stepsize and limits them
      if(debug .eq. 1)then
         write(6,'("Before Step Limitor, eata = ",es12.3)')eata
      endif      

      call steplimiter(eata,T,delt)

      if(debug .eq. 1)then
         write(6,'("Calculated that the order should change by",
     $              i2," with eata = ",es12.3)')changedir,eata
      endif

C     At this point, we know if the order will be changed. This depends
C     on changedir whether or not ordchange was true.
      if(changedir .eq. -1)then   ! For reducing the order
         call reduceord(delt)
         nordord=nordord-1
         nordsteps=0
      else if(changedir .eq. 1)then  ! increase the order
         nordord=nordord+1
         nordsteps=0
C        Make sure the next column is zero
         do i=1,nis
            vnord(i,nordord+1)=0.d0
         enddo
      endif


C     Write to the previous time vector (ptau)
      do i=1,6
         j=8-i
         ptau(j) = ptau(j-1)
      enddo
      ptau(1) = dtold

C     Finall, no matter what heppened, we need to change the time step
c     eata=eata*0.5d0
      do j=2,nordord+1  ! Correct the nordsiek vector for new time step
         do i=1,nis
            vnord(i,j)=vnord(i,j)*eata**(j-1)
         enddo
      enddo
      
C     Change the time step
      delta=delt*eata
      delt=delta


      return
  430 WRITE(900000+NREAC,9260) EPS,SUMX,NSTEP
      STOP
  555 WRITE(900000+NREAC,9998) nrestep
      STOP

 9021 FORMAT(132('='))
 9230 FORMAT(2X,'NS :',I7,4X,'NEW TS =',1PE9.2,' SEC',2X,'PRO =',
     1F10.5,2X,'AK =',D10.3,2X,'RAPID =',2X,A5,2X,'Y0=',1PE9.2,2X,
     2'YB=',1PE9.2,2X,'DEY/Y = ',0PE9.2)
 9260 FORMAT(///,5X,'SUMX ERROR GREATER THAN :',1PE9.2,5X,'SUMX=',
     10PF20.16,5X,'STEP NO :',I7,///,
     25X,'CHECK :   1) BARYON NUMBER CONSERVATION IN THE NETWORK',/,
     35X,'          2) NON-NEGATIVE REACTION RATES IN ROUTINE SRATE',/,
     45X,'          3) MATRIX INVERSION ROUTINE (EIGEN OR OTHER)',///,
     52X,6('*- PROGRAM ABORTED -*'))
 9998 FORMAT(///,132('='),///,2X,'NO CONVERGENCE',I7,' ITERATIONS')
      RETURN
      END

C======================================================================
C     Subroutine to calculate the first predictor step using the pascal
C     triagle matrix.
C     This is done by summing
      subroutine gearpredict

      include 'comgear.inc'

      integer ord

      ord=nordord

C     Copied from EPSODE code (See Byrne, acm transactions on
c       mathematical software, 1 (1975), pp. 71-96
      do j1 = 1,ord
         do j2 = j1,ord
            j = (ord + j1) - j2
            do i = 1,NIS
 210           vnord(i,j) = vnord(i,j) + vnord(i,j+1)
            enddo
         enddo
      enddo

      return
      end

C======================================================================
C     Subroutine to calculate the el values and other things related to
C     error calculations
      subroutine calcel(delt,ordchange)

      include 'comgear.inc'
      dimension ff(0:7)
      data ff/1.d0,1.d0,2.d0,6.d0,2.4d1,1.2d2,7.2d2,5.04d3/
      double precision rxi,rxiprodmo,rq,tauprod,rxiprod
      integer iq
      logical ordchange

      l=nordord+1

C     First get the xi values (ratios of times to current stepsize
      h=delt
      rq=dble(nordord)
      rxi=1.d0
      hsum=h
      hsum1=0.d0
      el(1)=1.d0
      el(2)=1.d0
      tauprod=1.d0
      rxiprod=1.d0
      do i=3,l     ! zero everything
         el(i)=0.d0
      enddo
C     Use reccursion relationship to calculate el
      do i=1,nordord-1
         hsum=hsum+ptau(i)    ! note that these are i+1 , so this loop
         hsum1=hsum1+ptau(i)  ! goes over 2->q for tauprod
c         print*,'hsum, hsum1',hsum,hsum1
         rxi=h/hsum1        ! 1/xi
         tauprod = tauprod*hsum/hsum1
         rxiprod = rxiprod*rxi
         do j=1,i+1
            iq = (i + 3)-j     ! eg. for i=3 we loop: 5,4,3,2
            el(iq) = el(iq) + el(iq-1)*rxi
         enddo
      enddo
c      print*,'tauprod',tauprod
      tauprod=1.d0+tauprod
c      rxiprod=rxiprodmo !*h/(hsum+ptau(nordord))
      rxiprodmo=rxiprod/rxi

C     The multiplacation factors for error calculations
      Efact(1)=-1.d0/(el(2)*tauprod)  ! For En(iq)
      Efact(2)=0.d0
      Efact(3)=0.d0
      Ceenmo = Ceen
      Ceen = tauprod/(ff(nordord+1)*rxiprod)
c      print*,"rxiprod, tauprod, ff(nordord+1)"
c      print*,rxiprod,tauprod,ff(nordord+1)

C     only calculate Qn, En(iq-1), En(iq+1) if signalled to do so
      if(ordchange)then
         xiprodmo=1.d0/rxiprodmo
c         xiprod=rxiprod*(hsum+ptau(nordord))/h
         xiprod=1.d0/rxiprod
         xipo=(hsum+ptau(nordord))/h
c         print*,'xipo',xipo
c         print*,'hsum, ptau(nordord),h',hsum,ptau(nordord),h
c         print*,'xiprod = ',xiprod,rq
         if(nordord .gt. 1)
     $        Efact(2)=-xiprodmo/(el(2)*(rq-1.d0)) ! En(q-1)
         Efact(3)=-xipo/((rq+1.d0)*(rq+2.d0)*el(2)*tauprod) ! En(q+1)
         Quen = (Ceen/Ceenmo)*(h/ptau(1))**(nordord+1)
         if(debug .eq. 1)then
            print*,"Calculate Order changes"
            print*,"  Ceen, Ceenmo,h/ptau(1),Quen"
            write(6,'(1x,5es11.3)')Ceen,Ceenmo,h/ptau(1),Quen
            print*,"  Efact(1),  Efact(2),  Efact(3)"
            write(6,'(1x,3es11.3)') Efact(1),Efact(2),Efact(3)
         endif
      endif

      if(DEBUG .eq. 1)then
         print*,"The el values are (0,1,...q+1)"
         WRITE(6,'(8f6.3)')(el(i),i=1,maxord+3)
         print*,"The error factors are"
         WRITE(6,'(3f9.3)')(Efact(i),i=1,3)
      endif

      return
      end

C======================================================================
C     Calculate the truncation errors that will be used for stepsize
C     changing as well as order changing if ordchange is set to true
      subroutine calcerrors(ordchange,stepError,yi,yf)

      include 'comgear.inc'
      double precision stepError(3),yi(nsp),yf(nsp),yscale(nsp),
     $     Earray(nis,2),maxcorr,terr,ar
      integer scalebyyi
      logical ordchange

C     Scale by yi?
      scalebyyi = 1
      if(scalebyyi .eq. 1)then
         do i=1,nis
c            yscale(i)=yi(i)
            yscale(i)=max(yi(i),ytmin)
         enddo
      else
         do i=1,nis
            yscale(i)=1.d0
         enddo
      endif
      

      do i=1,3
         stepError(i)=0.d0
      enddo

C     Current errors
      maxcorr = 0.d0
      terr = 0.d0
      nc=0
      do i=1,nis
         ar = (yf(i)-yi(i))/yscale(i)  !/vnord(i,1)
c         if(dabs(yi(i)) .lt. ytmin .or. dabs(yf(i)-yi(i)) .lt. dymin) 
c     $        ar = 0.d0
         if(dabs(ar) .gt. maxcorr) maxcorr=dabs(ar)
         
         terr = terr + ar*ar
      enddo
      stepError(1) = -Efact(1)*terr

C     For stepsize changes
      if(ordchange)then
         terrmo = 0.d0
         terrpo = 0.d0
         nc=0
         do i=1,nis
c            ar = 1.d0*(vnord(i,nordord))/yscale(i) 
            ar = (vnord(i,nordord+1)+corr(i)*el(nordord+1))/
     $           yscale(i) 
c            if(yi(i) .lt. ytmin) ar = 0.d0
            terrmo = terrmo + ar*ar

            if(DEBUG .eq. 1)then
               write(6,'(i4,3es12.3)')i,yf(i)-yi(i),ar,
     $              ((yf(i)-yi(i))-Quen*pcorr(i))/yscale(i)
            endif            

            ar = ((yf(i)-yi(i)) - Quen*pcorr(i))/yscale(i)
c            if(yi(i) .lt. ytmin) ar = 0.d0
            terrpo = terrpo + ar*ar

c            if(DEBUG .eq. 1)then
c               print*,i,ar
c            endif    

            Earray(i,1)=Efact(2)*vnord(i,nordord+1)
            Earray(i,2)=Efact(3)*(corr(i)-Quen*pcorr(i))
c            print*,Efact(3),Earray(i,2),Quen
c            if(dabs(Earray(i,1)) .gt. stepError(2))
c     $           stepError(2)=dabs(Earray(i,1))
c            if(dabs(Earray(i,2)) .gt. stepError(3))
c     $           stepError(3)=dabs(Earray(i,2))
         enddo
         stepError(2) = -Efact(2)*terrmo
         stepError(3) = -Efact(3)*terrpo
      endif

      if(debug .eq. 1)then
         print*,'The stepErrors are:'
         write(6,'(3es12.3)')(stepError(i),i=1,3)
      endif

      return
      end

C======================================================================
C     Reduce the order of the integration if needed.
      subroutine reduceord(delt)

      include 'comgear.inc'
      
      double precision scaler(7)
      integer iq

C     Return without doing anything if the order is 2
      if(nordord .eq. 2)return

      hsum=0.d0
      h=delt
      do i=1,7
         scaler(i)=0.d0
      enddo
      scaler(3)=1.d0
      hsum=0.d0
      do i=1,nordord-2
         hsum=hsum+ptau(i)
         xi=hsum/h
         do j=1,i+1
            iq = (i + 4)-j
            scaler(iq) = scaler(iq)*xi + scaler(iq-1)
         enddo
      enddo
      if(debug .eq. 1)then
         print*,"Reducing order with scalers:"
         write(6,'(7es12.3)')(scaler(i),i=1,7)
      endif
      do j=3,nordord
         do i=1,nis
c            print*,vnord(i,j),'-',vnord(i,nordord+1)*scaler(j)
            vnord(i,j)=vnord(i,j) - 
     $           vnord(i,nordord+1)*scaler(j)
         enddo
      enddo
      
      return
      end

C-------------------------------------------------------
C Subroutine to calculate the fluxes, energies, etc.
      subroutine fluxenergy(X,dtold)

      include 'comnuc.inc'

      DIMENSION Y(NSP),X(NSP)

      real*8 O

      O=1.d-100

      do i=1,nis
         y(i) = x(i)/an(i)
      enddo

      EQ=0.
      EQQ=0.
      DO 350 I=1,NR    ! calculate reaction fluxes and energies...
         I1=K1(I)
         I2=K2(I)
         I3=K3(I)
         I4=K4(I)
         I21=I2-1
         I41=I4-1
         A1A=((Y(I1)+O)**I21)*((Y(I3)+O)**I4)*I2*Y(I1)
         A2A=((Y(I3)+O)**I41)*((Y(I1)+O)**I2)*I4*Y(I3)
         FLUX(I)=V(I)*(A1A+A2A)/(F(I2)*F(I4)*(I2+I4))
C        FLUX(I)=V(I) x [rho x NA<sv>]
         
         AE=Q(I)*FLUX(I)
         AEE=Q(I)*V(I)*(Y(I1)+O)**I2*(Y(I3)+O)**I4/(F(I2)*F(I4))
         ENER(I)=AE*E
         EQ=EQ+AE
         EQQ=EQQ+AEE
  350 CONTINUE

C---------COMPUTATION OF N,P,A EXPOSURES     -----------
c     print*,INEUT,IPROT,IALPH
c     SYN=SYN+RHO*DELT*Y0(INEUT)
c     SYP=SYP+RHO*DELT*Y0(IPROT)
c     SYA=SYA+RHO*DELT*Y0(IALPH)


C     Write some stepping info to file
      sumen = 0.d0
      do i=1,nr
         sumen=sumen+ener(i)
      enddo

      DO 1273 I=1,NR
        TOTFLUX(I)=TOTFLUX(I)+DTOLD*FLUX(I)
        tenergy(i)=tenergy(i)+dtold*ener(i)
 1273 continue
C     FLUX(I) are reaction fluxes

      TOTENER=0.D0
      DO I=1,NR
 1275    TOTENER=TOTENER+ENER(I)
      enddo
      TOTENERGY=TOTENERGY+DTOLD*TOTENER
      
      return
      end



