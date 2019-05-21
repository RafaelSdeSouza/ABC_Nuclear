c     For Gear's method and MA48 matrix solver
      include 'comnuc.inc'

      integer la,ne,alength
      parameter (la=nr*40) 
      PARAMETER(MINORD=1,MAXORD=5,jacrate=1)
     
      integer irow(la),jcol(la),ent(8,nr)
      integer icntl(20),info(20),iw(9*nsp),job,
     *     keep(10*nsp+7),diagindex(nsp)

      COMMON/MA48/amat(la),cntl(10),rinfo(10),error(3),
     *     w(4*nsp),rscale(nis),cscale(nis),
     *     irow,jcol,ent,icntl,iw,job,keep,ne,
     *     diagindex
      COMMON/GEAR/ptau(7),el(8),pcorr(NSP),corr(NSP),
     *	   vnord(NSP,7)
      COMMON/GEARFACT/tauprod,xiprod,Ceen,Quen,Efact(3),factor
      COMMON/GEARINT/nordord,nordsteps,nmakejac,
     *	   nrestep

      PARAMETER (nmax=3500,kmaxx=7,imax=kmaxx+1,safe1=0.25d0,
     *     safe2=0.7d0,redmax=1.d-5,redmin=0.7d0,tiny=1.d-30,
     *     scalmx=0.1d0)
      COMMON/BADERREAL/cost(imax),alf(kmaxx,kmaxx),tnew,kmax
      COMMON/BADERINT/nseq(imax)
