      implicit double precision (a-h, o-z), integer (i-n)

      PARAMETER(NIS=213,NSP=NIS+1,NR=2373,NSL=60,NKT=7,KP=50000)
      PARAMETER(NDS=56,NDI=56,IJD=9,IDEL=IJD-1,JDEL=IJD-1)
      PARAMETER(SECY=3.1536D+7,XNAV=6.023D+23,E=9.636D+17)
      PARAMETER(NPROF=70000,nmod=5000,ntsh=5000)
      PARAMETER(DEBUG=0)
      PARAMETER(NWEAK=1000)
      
      character*5  ONET,ON,ONO,AAAAA,AAA,ZZZZZZ,OPR,BBBBB
      character*5  SPE1,SPE2,SPE3,SPE4,NAME,EMITSPE
      character    KNET*1,IT*1,ABCD*80,char_ec*2,char_w*1
      character*25 reactionfile,networkfile,weakinputfile,inputfile
      character*70 reacstring
      character*30 reacstrweak

      COMMON/CEMIT/emitspe(NR)
      COMMON/CINTER/NLAST,reactionfile,networkfile,weakinputfile,
     * inputfile
      COMMON/CONTROL/INMODE,ITEST,ilast,iswitch,iweak,isop1,isop2
      COMMON/CMODE/ISF,KTEST,DEY,IFLUX,nsolver,nrw
      COMMON/CSHELL/ishrun,nshmax,numod,nsh0(nmod),ishell,modmaxnum,
     * numprof

      COMMON/RRATES/temp(60),prob(NR),rate(NR,60),fu(NR,60),templog(60),
     * ratelog(NR,60),rate2(NR,60),char_ec(NR),reacstring(NR),char_w(NR)
      COMMON/WRATES/ratew(nweak,60,12),ratelogw(nweak,60,12),
     * reacstrweak(nweak),rhoyelog(12)

      COMMON/CM/K1(NR),K2(NR),K3(NR),K4(NR),K5(NR),K6(NR),K7(NR),K8(NR),
     * totpartin(NR)
      COMMON/CMM/numab(7)
      COMMON/CNETW/AN(NSP),ZN(NSP),XINI(NSP),Q(NR),ENER(NR),ON(NSP)
      COMMON/CPARM/PAR,PAMIN,DYMIN,DEMIN,YMIN,YTMIN,EPS,ATEST,TES,F(0:6)
      COMMON/CENER/EB,EN,EQ,ER,EQQ,FLUX(NR),TOTFLUX(NR),V(0:NR),SUMENER,
     * TOTENER,TOTENERGY,vmult(0:nr),tenergy(nr)
      COMMON/CEXPL/TAU0,RHO0,TAU,RHO,oldtau,oldrho,oldye,
     * TPROF(NPROF),TAUPROF(NPROF),RHOPROF(NPROF),TotErr(nsp),
     * scale1,scale2,tauprofmax,tauprocalmax,ncurrtime
      COMMON/SHELL0/taush0(nmod,ntsh),tmaxsh(ntsh),rhosh0(nmod,ntsh),
     * vconsh0(nmod,ntsh),clensh0(nmod,ntsh),dmsh0(nmod,ntsh),
     * tmaxcode(ntsh)
      COMMON/SHELL/tmodel(nmod),taushell(nmod,ntsh),
     * rhoshell(nmod,ntsh),vconshell(nmod,ntsh),clenshell(nmod,ntsh),
     * dmshell(nmod,ntsh),xshells(nsp,ntsh),trun(ntsh),tenersh(ntsh),
     * taus(ntsh),rhos(ntsh),vcons(ntsh),clens(ntsh),dms(ntsh),dmst
      COMMON/CLAST/TLAST,XLAST,DELTA,delta0,delfac,NSTEP,NREAL,gres
      COMMON/CET1/ETA,YE
      COMMON/CET2/TIME(KP),DENN(KP),TEMM(KP),ERATEY(KP),ERATEQ(KP),
     * ENTOTY(KP),ENTOTQ(KP),NNSTEP(KP)
      COMMON /MATAT2/AV(JDEL,NSP),AH(NSP,IDEL),AD(NDI+NDS+1,NSP)
      COMMON/AA/A1(NR),A2(NR),A3(NR),A4(NR),A5(NR),A6(NR),A7(NR),A8(NR)
      COMMON/AVHD/ATV(JDEL,NSP),ATH(NSP,IDEL),ATD(NDI+NDS+1,NSP)
