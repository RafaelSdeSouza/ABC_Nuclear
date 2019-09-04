c     ==============================================================
      subroutine CNshell_ele(nis,on,xf,zexp,nexp,iflag,xele)
c     ==============================================================
c     this script needs to be adapted to the astrophysical site
c     under study     
c
c     SETUP HERE IS FOR NOVAE
c
c     nis:    number of isotopes in network
c     on:     vector with nuclide labels
c     xf:     vector of predicted isotopic abundances
c     zexp:   vector with observed abundance atomic number
c     nexp:   number of observational data points
c     iflag:  flag for absolute vs ratio abundance input

      integer lll,j,i,nis,iflag,zexp(nexp),nxx,nexp
      real*8 trun
      real*8 xf(nis+1),xele(nexp)
      real*8 x_h,x_he,x_li,x_be,x_b,x_c,x_n,x_o,x_f,x_ne,x_na,x_mg
      real*8 x_al,x_si,x_p,x_s,x_cl,x_ar,x_k,x_ca,x_sc,x_ti,x_v,x_cr
      real*8 xx(100)
      real*8 xx_norm
      character*100 string,string2
      character*5 on(nis+1)

c     ==============================================================
c     analyze results          

c     first task is to add abundances of stable isotopes of a given
c     element [summation of stable and very long-lived species in rows]

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
     
c     next task: add abundances of radioactive nuclides that
c     contribute to abundance of a given element; 
c     this part may need to be modified, depending on the 
c     scenario; see also comments in the header;

c     !!! index of xx() is equal to atomic number !!!

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
c     P = 31P + 31Si + 23P + 33P
      xx(15)=x_p+xf(109)+xf(115)+xf(116)       
c     S = 23S + 33S + 34S + 36S + 35S 
      xx(16)=x_s+xf(129)       
c     Cl = 35Cl + 37Cl + 36Cl
      xx(17)=x_cl+xf(132)      
c     Ar = 36Ar + 38Ar + 40Ar + 37Ar + 39Ar
      xx(18)=x_ar+xf(139)+xf(148)              
c     K = 39K + 40K +41K + 41Ar
      xx(19)=x_k+xf(156)       
c     Ca = stable Ca + 41Ca + 45Ca + 47Ca + 43Sc + 44Sc + 42K + 43K
      xx(20)=x_ca+xf(155)+xf(174)+xf(185)+xf(163)+xf(169)
     &   +xf(162)+xf(165)      
c     Sc = 45Sc + 45Ti  + 46Sc
      xx(21)=x_sc+xf(173)+xf(183)      
c     Ti = stable Ti + 44Ti + 47Sc + 48Sc + 47V
      xx(22)=x_ti+xf(167)+xf(184)+xf(194)+xf(188)      
c     V = 50V + 51V + 48V + 49V + 48Cr + 49Cr
      xx(23)=x_v+xf(195)+xf(198)+xf(191)+xf(196)         
      xx(24)=x_cr+xf(206)    ! Cr = stable Cr + 51Cr

      nxx=24   ! nxx: number of element species analyzed 
         
c     =============================================================
c     if iflag.ne.0: predicted mass fraction ratios [i.e., relative  
c     to X(iflag)] will be analyzed
c     =============================================================
      if(iflag.ne.0)then
            
c        divide all elemental abundances by xx(iflag)
         xx_norm=xx(iflag) 
c        if norm element abundance is negative, set all abundances
c        to a very small value, so that run is discarded because of 
c        poor fitness
         if(xx_norm.lt.0.d0)then    
c            write(6,*) 'xx_norm is negative. Run stop.'            
            write(6,*) 'xx_norm is negative'            
            xx_norm=1.d0
            do 112 jj=1,nxx
               xx(jj)=1.d-99
 112        continue
c            stop
         endif
         do 111 jj=1,nxx
            xx(jj)=xx(jj)/xx_norm
 111     continue

      endif
         
c     =============================================================
c     compute predicted elemental abundances for elements that
c     have been observed         
         
      do 91 j=1,nexp
         xele(j)=xx(zexp(j))
   91 continue
                  
c     =============================================================
 6000 format(a100)
 6001 format(25(1pe10.2e3,2x))
     
      return
      end


