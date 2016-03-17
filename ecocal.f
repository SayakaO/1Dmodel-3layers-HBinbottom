c-----------------------------------------------------------------------
c
c     -- Physical Properties and State Variables --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine ecocal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      if(nt.ne.1) goto 200
c
c     -- density --
c
      do 10 k=1,nz
       dens(k)=1028.14d0-.0735d0*tmp(k)-.00469d0*tmp(k)*tmp(k)
     &                  -35.d0*(.802d0-.002d0*tmp(k))
   10 continue
c
  200 continue
c
c     -- vertical eddy diffusivity --
c
	call khcal
c
c     -- source terms --
c
      call pro
c
c     -- water temperature --
c
	ww=0.d0
	call adcal(ww,tmp,tmpn,qtmp,fkh)
c
c     -- state variables in the ecosystem --
c
c	if(nsw(6).eq.1) then
	if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
       do 20 m=1,np
        ww=wphy(m)
        call adcalf(m,ww,phy,phyn,qphy,fkh)
   20  continue
c
       do 22 m=1,nzp
        ww=0.d0
        call adcalf(m,ww,zoo,zoon,qzoo,fkh)
   22  continue
c
c       ww=0.d0
c       call adcal(ww,bac,bacn,qbac,fkh)
       ww=wpoc
       call adcal(ww,poc,pocn,qpoc,fkh)
       ww=0.d0
       call adcal(ww,doc,docn,qdoc,fkh)
       call adcal(ww,dip,dipn,qdip,fkh)
       call adcal(ww,dinh,dinhn,qdinh,fkh)
c      call adcal(ww,dis,disn,qdis,fkh)
       call adcal(ww,dox,doxn,qdox,fkh)
       call adcalb(ww,hb,hbn,qhb,fkh)
c
	endif
c
c     -- update water temperature and density --
c
      do 30 k=0,nz
       tmp(k)=tmpn(k)
       dens(k)=1028.14d0-.0735d0*tmp(k)-.00469d0*tmp(k)*tmp(k)
     &                  -35.d0*(.802d0-.002d0*tmp(k))
   30 continue
c
c     -- update state variables in the ecosystem --
c
c	if(nsw(6).eq.1) then
	if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
       do 32 k=1,nz
        do 34 m=1,np
         phy(m,k)=phyn(m,k)
         if(phy(m,k).le.1.d-8) phy(m,k)=1.d-8
   34   continue
        do 36 m=1,nzp
         zoo(m,k)=zoon(m,k)
         if(zoo(m,k).le.1.d-8) zoo(m,k)=1.d-8
   36   continue
c	  bac(k)=bacn(k)
        poc(k)=pocn(k)
        doc(k)=docn(k)
        dip(k)=dipn(k)
        dinh(k)=dinhn(k)
c        dis(k)=disn(k)
        dox(k)=doxn(k)
        hb(k)=hbn(k)
c        if(bac(k).le.1.d-8) bac(i,j,k)=1.d-8
        if(poc(k).le.1.d-8) poc(k)=1.d-8
        if(doc(k).le.1.d-8) doc(k)=1.d-8
        if(dip(k).le.1.d-8) dip(k)=1.d-8
        if(dinh(k).le.1.d-8) dinh(k)=1.d-8
c        if(dis(k).le.1.d-8) dis(k)=1.d-8
        if(dox(k).le.1.d-8) dox(k)=1.d-8
        if(hb(k).le.1.d-8) hb(k)=1.d-8
   32  continue
c
	endif
c
c     -- vertical mixing --
c
      if(nsw(7).ne.2) call mixcal
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Vertical Eddy Diffusivity Coefficients --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine khcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- Kh is constant --
c
	if(nsw(7).eq.0) then
c
	 if(nsw(2).eq.1) then
	  do 10 k=0,nz
	   fkh(k)=fkh0
   10   continue
       endif
c
c	-- Kh depends on stratification function --
c
	elseif(nsw(7).eq.1) then
c
	 if(nsw(2).eq.1) then
c
        do 20 k=1,nz-1
c
         rhodz=(dens(k+1)-dens(k))/dzh(k)
	   if(rhodz.lt.0.d0) rhodz=0.d0
         rhom=(dens(k)*ddz(k+1)+dens(k+1)*ddz(k))*.5d0/dzh(k)
         udz=(ui(k+1)-ui(k))*.5d0/dzh(k)
         vdz=(vj(k+1)-vj(k))*.5d0/dzh(k)
         uvmdz2=udz*udz+vdz*vdz
         epsi=btkc*g*rhodz/(rhom*((fkhmin/fkh0)**(1.d0/alkc)-1.d0))
c
         if(rhodz.gt.1.0d-10) then
          if(uvmdz2.le.epsi) then
           fkh(k)=fkhmin
          else
           fkh(k)=fkh0*(1.d0+btkc*g*rhodz/rhom/uvmdz2)**alkc
          endif
         endif
c
         if(fkh(k).lt.fkhmin) fkh(k)=fkhmin
c
   20   continue
c
	 endif
c
	endif
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Production Terms --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine pro
      implicit double precision (a-h,o-z)
      include 'param.h'
      dimension b1(npmax),b2(npmax),b3(npmax),b4(npmax)
      dimension b6(npmax),b7(npmax),b8(npmax),b9(npmax)
      dimension d1(npmax),d2(npmax),d3(npmax)
c    エラー対策↓-----------param.h に追加-------------------------------------------------------*解決？
c      common /pol/qpoly(nzmax),poly(nzmax),bother(nzmax),qbother(nzmax)
c
      do 10 k=nz,1,-1
c
	 qtmp(k)=0.d0
c
c       if(nsw(6).eq.0) ec=0.2d0
       ec=0.2d0
c
c     -- bulk formula --
c
       if(nsw(4).eq.0) then 
c
        q1=qs0*(1.d0-.71d0*cloud)*(1.d0-ref)
c
c        q2=s_s*sigma*(tmpaq+273.15d0)**4.d0
c     &     *(.39d0-5.8d-2*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
c     &     +4.d0*s_s*sigma*(tmpaq+273.15d0)**3.d0
c     &     *(tmp(i,j,k)-tmpaq)
        q2=s_s*sigma*(tmp(k)+273.15d0)**4.d0
     &     *(.49d0-.066d0*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
     &     +4.d0*s_s*sigma*(tmp(k)+273.15d0)**3.d0
     &     *(tmp(k)-tmpaq)
c
        ww=sqrt(wx*wx+wy*wy)
        if(ww.ge.1.d0)then
         q3=f_i*densa*ce*(q_s-q_e)*ww
         q4=c_a*densa*ch*(tmp(k)-tmpa)*ww
        else
         tmpfc=(tmp(k)-tmpa)+.61d0*(tmpa+273.15d0)*(q_s-q_e)
         if(tmpfc.lt.0.d0)then
          q3=-f_i*densa*ce2*(q_s-q_e)*(abs(tmpfc))**(1.d0/3.d0)
          q4=-c_a*densa*ch2*(tmp(k)-tmpa)
         elseif(tmpfc.ge.0.d0)then
          q3=f_i*densa*ce2*(q_s-q_e)*tmpfc**(1.d0/3.d0)
          q4=c_a*densa*ch2*(tmp(k)-tmpa)*tmpfc**(1.d0/3.d0)
         end if
        endif
c
       elseif(nsw(4).eq.1) then
c
        e_s=6.1078d0*(10.d0**(7.5d0*tmpaq/(237.3d0+tmpaq)))
        e_a=e_s*humq
        q_s=.622d0*e_s/(presq-.378d0*e_s)
        q_e=.622d0*e_a/(presq-.378d0*e_a)
        densa=1.293d0*273.15d0*(presq-.378d0*e_a)
     &       /(273.15d0+tmpaq)/1013.25d0
c
        q1=sunnq*(1.d0-ref)
c
c        q2=s_s*sigma*(tmpaq+273.15d0)**4.d0
c     &     *(.39d0-5.8d-2*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
c     &     +4.d0*s_s*sigma*(tmpaq+273.15d0)**3.d0
c     &     *(tmp(i,j,k)-tmpaq)
        q2=s_s*sigma*(tmp(k)+273.15d0)**4.d0
     &     *(.49d0-.066d0*sqrt(e_a))*(1.d0-clat*(cldq**2.d0))
     &     +4.d0*s_s*sigma*(tmp(k)+273.15d0)**3.d0
     &     *(tmp(k)-tmpaq)
c
        ww=sqrt(wx*wx+wy*wy)
        if(ww.ge.1.d0)then
         q3=f_i*densa*ce*(q_s-q_e)*ww
         q4=c_a*densa*ch*(tmp(k)-tmpaq)*ww
        else
         tmpfc=(tmp(k)-tmpaq)
     &         +.61d0*(tmpaq+273.15d0)*(q_s-q_e)
         if(tmpfc.lt.0.d0)then
          q3=-f_i*densa*ce2*(q_s-q_e)*(abs(tmpfc))**(1.d0/3.d0)
          q4=-c_a*densa*ch2*(tmp(k)-tmpaq)
     &       *(abs(tmpfc))**(1.d0/3.d0)
         elseif(tmpfc.ge.0.d0)then
          q3=f_i*densa*ce2*(q_s-q_e)*tmpfc**(1.d0/3.d0)
          q4=c_a*densa*ch2*(tmp(k)-tmpaq)
     &           *tmpfc**(1.d0/3.d0)
         endif
        endif
	 endif
c
       if(k.eq.1) then
        fu=q1*(exp(-ec*zb(k-1))-exp(-ec*zb(k)))/(dens0*c_0)
     &    -(q2+q3+q4)/(dens0*c_0)
       else
        if(k.eq.nz) then
         fu=q1*(exp(-ec*zb(k-1))-0.d0)/(dens0*c_0)
        else
         fu=q1*(exp(-ec*zb(k-1))-exp(-ec*zb(k)))/(dens0*c_0)
        endif
       endif
c
       qtmp(k)=fu/ddz(k)
c
c     -- river --
c
       if(k.eq.1) then
        do 20 n=1,nriv
         qtmp(k)=qtmp(k)+qriv(n)*(triv(n)-tmp(k))
     &                  /(ddx*ddy*ddz(k)+qriv(n))
         goto 22
   20   continue
       endif
   22  continue
c
c	 if(nsw(6).eq.1) then
	 if(nsw(6).eq.1.and.mod(nt,intbal).eq.0) then
c
c     -- extinction coefficient --
c
        physum=0.d0
        phys1=0.d0
        chlaave=0.d0
        do 30 m=1,np
         physum=physum+phy(m,k)
         phys1=phys1+phy(m,1)
         chlaave=chlaave+chlac(m)
   30   continue
        chlaave=chlaave/np
        ec=ec0+ec1*chlaave*phys1
c
c     -- phytoplankton --
c
        do 40 m=1,np
         b1(m)=gp(m)*thep(m)**(tmp(k)-20.d0)
     &        *dmin1(dip(k)/(dip(k)+hsp(m)),
     &               dinh(k)/(dinh(k)+hsn(m)))
     &        *(dexp(1.d0-q1*dexp(-ec*zb(k))/aip(m))
     &         -dexp(1.d0-q1*dexp(-ec*zb(k-1))/aip(m)))
     &         /ec/ddz(k)
     &        *phy(m,k)
	   if(k.eq.nz) then
          b1(m)=gp(m)*thep(m)**(tmp(k)-20.d0)
     &         *dmin1(dip(k)/(dip(k)+hsp(m)),
     &          dinh(k)/(dinh(k)+hsn(m)))
     &         *(dexp(1.d0-q1*0.d0/aip(m))
     &          -dexp(1.d0-q1*dexp(-ec*zb(k-1))/aip(m)))
     &          /ec/ddz(k)
     &         *phy(m,k)
	   endif
         b2(m)=rp(m)*thep(m)**(tmp(k)-20.d0)*phy(m,k)
         b3(m)=erp(m)*dexp(gamp(m)*chlac(m)*phy(m,k))*b1(m)
         b4(m)=dmp(m)*phy(m,k)*phy(m,k)
   40   continue
c
c     -- zooplankton --
c   
        zoosum=0.d0
        do m2=1,nzp
         zoosum=zoosum+zoo(m,k)
         end do 
        do 42 m=1,nzp
c         b6(m)=az(m)*cz(m)*hsz(m)*(physum+poc(i,j,k))
c     &        /(hsz(m)+physum+poc(i,j,k))*zoo(m,i,j,k)
         b6(m)=cz(m)*thez(m)**(tmp(k)-20.d0)*(1.d0-
     &         dexp(.007d0*(hsz(m)-physum-poc(k))))*zoo(m,k)
         b7(m)=rz(m)*thez(m)**(tmp(k)-20.d0)*zoo(m,k)
         b8(m)=(1.d0-az(m))*b6(m)
         b9(m)=dmz(m)*zoo(m,k)*zoo(m,k)
   42   continue
c
c     -- particulate organic carbon --
c
        b10=rpoc*thepoc**(tmp(k)-20.d0)*poc(k)
        b11=ft*b13
c
c     -- dissolved organic carbon --
c
        b13=rdoc*thedoc**(tmp(k)-20.d0)*doc(k)
c
c     -- dissolved oxygen --
c
        if(k.eq.1) then
         dosu=14.161d0-.3943d0*tmp(k)+7.714d-3*tmp(k)**2.d0
     &       -6.46d-5*tmp(k)**3.d0
     &       -0.d0*(.1519d0-4.62d-3*tmp(k)+6.76d-5*tmp(k)**2.d0)
         b17=rea*(dosu-dox(k))/ddz(k)
        else
         b17=0.d0
        endif
c   -- hervibous benthos -- 
c feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
         pq4hb=0.101
         hbsr2phy=0.133
         vulhbandphy=1.0
      b18=pq4hb*(hbsr2phy*vulhbandphy*physum*hb(k)/2*vulhbandphy
     &    +hbsr2phy*hb(k))
c
c  other mortality = (1-EE)*PB*B
        ee4hb=0.383
        pb4hb=1.875d-09
      b20 = (1-ee4hb)* pb4hb * hb(k)
c
c   -- polychaeta -- 
c         
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c         pq4poly=0.205
c         polysr2det=0.0222
c         vulpolyanddet=1.03
c      b31=pq4poly*(polysr2det*vulpolyanddet*poc(k)*poly(k)/2
c      & *vulpolyanddet+polysr2det*poly(k))
c
c  other mortality = (1-EE)*PB*B
c        ee4poly=0.010
c        pb4poly=2.708d-09
c      b32 = (1-ee4poly)* pb4poly * poly(k)
c
c   -- other benthos -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c      other benthos VS detritus         
c          pq4bother=0.185
c          othersr2det=0.0089
c          vulotheranddet=1.76
c        bofd1=pq4bother*othersr2det*vulotheranddet*poc(k)*bother(k)
c        bofd2=2*vulotheranddet+othersr2det*bother(k)
c    エラー箇所↓-----------------------------------------------------------------------------*解決？
c       obvsdet=bofd1/bofd2
c       other benthos VS  zooplankton         
c          othersr2zoop=0.035
c          vulotherandzoop=1.9d+12
c        bofz1=pq4bother*othersr2zoop*vulotherandzoop*zoosum*bother(k)
c        bofz2=2*vulotherandzoop+othersr2zoop*bother(k)
c    エラー箇所↓-----------------------------------------------------------------------------*解決？
c      obvszp=bofz1/bofz2
c      b41=obvsdet + obvszp
c
c  other mortality = (1-EE)*PB*B
c        ee4bother=0.125
c        pb4bother=1.957d-09
c      b42=(1-ee4bother)* pb4bother * bother(k)
c
c   -- small pelagic fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c      small pelagic fish VS phtoplankton        
c          pq4spf=0.113
c          spfsr2phy=0.05733
c          vulspfandphy=3.41
c        spffphy1=pq4spf*spfsr2phy*vulspfandphy*physum*spfish(k)
c        spffphy2=2*vulspfandphy+spfsr2phy*spfish(k)
c       spfvsphy=spffphy1/spffphy2
c      small pelagic fish VS  zooplankton  
c          spfsr2zoop=0.02266
c          vulspfandzoop=1.00
c        spffz1=pq4bother*spfsr2zoop*vulspfandzoop*zoosum*spfish(k)
c        spffz2=2*vulspfandzoop+spfsr2zoop*spfish(k)
c 
c      spfvszp=spffz1/spffz2
c      b51=obvsdet + obvszp
c
c  other mortality = (1-EE)*PB*B
c        ee4spf=0.294
c        pb4spf=1.283d-09
c      b52=(1-ee4spf)* pb4spf * spfish(k)
c
c   -- small demersal fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c      small demersal fish VS   hervibous benthos      
c          pq4sdf=0.106
c          sdfsr2hb=0.0055
c          vulsdfandhb=1.00
c        sdffhb1=pq4sdf*sdfsr2hb*vulsdfandhb*hb(k)*sdfish(k)
c        sdffhb2=2*vulsdfandhb+sdfsr2hb*sdfish(k)
c       sdfvshb=sdffhb1/sdffhb2
c      small demersal fish VS   polychaeta      
c          pq4sdf=0.106
c          sdfsr2poly=0.000639
c          vulsdfandpoly=2.06
c        sdffpoly1=pq4sdf*sdfsr2poly*vulsdfandpoly*poly(k)*sdfish(k)
c        sdffpoly2=2*vulsdfandpoly+sdfsr2poly*sdfish(k)
c       sdfvspoly=sdffpoly1/sdffpoly2
c      spfvspoly=spffz1/spffz2
c      small demersal fish VS   other benthos      
c          pq4sdf=0.106
c          sdfsr2other=0.0063
c          vulsdfandother=37.6
c        sdffother1=pq4sdf*sdfsr2other*vulsdfandother*bother(k)*sdfish(k)
c        sdffother2=2*vulsdfandother+sdfsr2other*sdfish(k)
c       sdfvsother=sdffother1/sdffother2
c      small demersal fish VS   detritus      
c          pq4sdf=0.106
c          sdfsr2det=0.00737
c          vulsdfanddet=1.0202
c        sdffdet1=pq4sdf*sdfsr2det*vulsdfanddet*poc(k)*sdfish(k)
c        sdffdet2=2*vulsdfanddet+sdfsr2det*sdfish(k)
c       sdfvsdet=sdffdet1/sdffdet2
c
c      b61= sdfvshb + sdfvsother + sdfvsdet + sdfvspoly
c
c  other mortality = (1-EE)*PB*B
c        ee4sdf=0.248
c        pb4sdf=1.030d-09
c      b62=(1-ee4sdf)* pb4sdf * sdfish(k)
c
c   -- piscivorous fish -- 
c  feeding = (P/Q)*(search rate*vulnerability*B prey*B predator/2*vul+search rate*B predator)
c  piscivorous fish VS  hervibous benthos      
c          pq4pfish=0.091
c          pfishsr2hb=0.003357
c          vulpfishandhb=1.00
c        pfishfhb1=pq4pfish*pfishsr2hb*vulpfishandhb*hb(k)*pfish(k)
c        pfishfhb2=2*vulpfishandhb+pfishsr2hb*pfish(k)
c       pfishvshb=pfishfhb1/pfishfhb2
c  piscivorous fish VS   polychaeta      
c          pq4pfish=0.091
c          pfishsr2ply=0.00056
c          vulpfishandply=1.09
c       pfishfply1=pq4pfish*pfishsr2ply*vulpfishandply*poly(k)*pfish(k)
c       pfishfply2=2*vulpfishandply+pfishsr2ply*pfish(k)
c       pfishvspoly=pfishfply1/pfishfply2
c  piscivorous fish VS   other benthos      
c          pq4pfish=0.091
c          pfishsr2other=0.003817
c          vulpfishandother=1.01
c       pfishfother1=pq4pfish*pfishsr2other*vulpfishandother*bother(k)*pfish(k)
c       pfishfother2=2*vulpfishandother+pfishsr2other*pfish(k)
c       pfishvsother=pfishfother1/pfishfother2
c  piscivorous fish VS  small pelagic fish      
c          pq4pfish=0.091
c          pfishsr2spf=0.041447
c          vulpfishandspf=1.65
c       pfishfspf1=pq4pfish*pfishsr2spf*vulpfishandspf*spfish(k)*pfish(k)
c       pfishfspf2=2*vulpfishandspf+pfishsr2spf*pfish(k)
c       pfishvsspf=pfishfspf1/pfishfspf2
c  piscivorous fish VS  small demersal fish      
c          pq4pfish=0.091
c          pfishsr2sdf=0.018909
c          vulpfishandsdf=1.15
c       pfishfsdf1=pq4pfish*pfishsr2sdf*vulpfishandsdf*sdfish(k)*pfish(k)
c       pfishfsdf2=2*vulpfishandsdf+pfishsr2sdf*pfish(k)
c       pfishvssdf=pfishfsdf1/pfishfsdf2
c  piscivorous fish VS  zooplankton      
c          pq4pfish=0.091
c          pfishsr2zoop=0.0034
c          vulpfishandzoop=1.7
c       pfishfzoop1=pq4pfish*pfishsr2zoop*vulpfishandzoop*zoosum*pfish(k)
c       pfishfzoop2=2*vulpfishandzoop+pfishsr2zoop*pfish(k)
c       pfishvszoop=pfishfzoop1/pfishfzoop2
c      b71= pfishvshb + pfishvsother + pfishvspoly + pfishvssdf + pfishvsspf + pfishvszoop
c
c  other mortality = (1-EE)*PB*B
c        ee4pfish=0.373
c        pb4pfish=5.798d-10
c      b72=(1-ee4pfish)* pb4pfish * pfish(k)
c
c     -- river --
c
        if(k.eq.1) then
         if(nriv.ge.1) then
          do 50 n=1,nriv
           qph=qriv(n)*0.d0/ddx/ddy/ddz(k)
           qzo=qriv(n)*0.d0/ddx/ddy/ddz(k)
           qpc=qriv(n)*(pcriv(n)-poc(k))/ddx/ddy/ddz(k)
           qdc=qriv(n)*(dcriv(n)-doc(k))/ddx/ddy/ddz(k)
           qdp=qriv(n)*(dpriv(n)-dip(k))/ddx/ddy/ddz(k)
           qdnh=qriv(n)*(dnhriv(n)-dinh(k))/ddx/ddy/ddz(k)
           qdo=qriv(n)*(dxriv(n)-dox(k))/ddx/ddy/ddz(k)
           goto 52
   50     continue
         endif
        else
         qph=0.d0
         qzo=0.d0
         qpc=0.d0
         qdc=0.d0
         qdp=0.d0
         qdnh=0.d0
         qds=0.d0
         qdo=0.d0
         endif 
   52   continue
c
c     -- prediction of state variables --
c
        b1psum=0.d0
        b1nsum=0.d0
        b1ssum=0.d0
        b2psum=0.d0
        b2nsum=0.d0
        b2ssum=0.d0
        b3sum=0.d0
        b4sum=0.d0
        b6sum=0.d0
        b7psum=0.d0
        b7nsum=0.d0
        b7ssum=0.d0
        b8sum=0.d0
        b9sum=0.d0
        d1sum=0.d0
        d2sum=0.d0
        d3sum=0.d0
c
        do 60 m=1,np
         b1psum=b1psum+dipph(m)*b1(m)
         b1nsum=b1nsum+dinph(m)*b1(m)
         b1ssum=b1ssum+disph(m)*b1(m)
         b2psum=b2psum+dipph(m)*b2(m)
         b2nsum=b2nsum+dinph(m)*b2(m)
         b2ssum=b2ssum+disph(m)*b2(m)
         b3sum=b3sum+b3(m)
         b4sum=b4sum+b4(m)
   60   continue
        do 62 m=1,nzp
         b6sum=b6sum+b6(m) 
         b7psum=b7psum+dipzo(m)*b7(m) 
         b7nsum=b7nsum+dinzo(m)*b7(m) 
         b7ssum=b7ssum+diszo(m)*b7(m) 
         b8sum=b8sum+b8(m) 
         b9sum=b9sum+b9(m) 
   62   continue
        do 64 m=1,np
         d1sum=d1sum+todph(m)*b1(m)
         d2sum=d2sum+todph(m)*b2(m)
   64   continue
        do 66 m=1,nzp
         d3sum=d3sum+todzo(m)*b7(m)
   66   continue
c
        do 70 m=1,np
         qphy(m,k)=b1(m)-b2(m)-b3(m)-b4(m)
     &            -b6sum*phy(m,k)/(physum+poc(k))+qph
   70   continue
        do 72 m=1,nzp
         qzoo(m,k)=b6(m)-b7(m)-b8(m)-b9(m)+qzo
   72   continue
        qpoc(k)=b4sum-b6sum*poc(k)/(physum+poc(k))
     &         +b8sum+b9sum-b10-b11+qpc
        qdoc(k)=b3sum+b11-b13+qdc
c
        dipzod=(dipph(1)*(b7(1)+b8(1)+b9(1))-dipzo(1)*b9(1))
     &        /(b7(1)+b8(1))
        dinzod=(dinph(1)*(b7(1)+b8(1)+b9(1))-dinzo(1)*b9(1))
     &        /(b7(1)+b8(1))
        b8psud=dipzod*b8(1)
        b8nsud=dinzod*b8(1)
c
        dppomd=(dipph(1)*b4(1)+dipzod*b8(1)+dipzo(1)*b9(1))
     &        /(b4(1)+b8(1)+b9(1))
	  dnpomd=(dinph(1)*b4(1)+dinzod*b8(1)+dinzo(1)*b9(1))
     &        /(b4(1)+b8(1)+b9(1))
c
        dpdomd=(dipph(1)*b3(1)+dppomd*b11)
     &        /(b3(1)+b11)
	  dndomd=(dinph(1)*b3(1)+dnpomd*b11)
     &        /(b3(1)+b11)
c
       qdip(k)=-b1psum+b2psum+b7psum+dppomd*b10+dpdomd*b13+qdp
       qdinh(k)=-b1nsum+b2nsum+b7nsum+dnpomd*b10+dndomd*b13+qdnh
       qdox(k)=d1sum-d2sum-d3sum-topom*b10-todom*b13+b17+qdo
c
       qhb(k)=b18-b20
c       hbn(k)=hb(k)+tbal*dt*dhb(k)
         write (*,*)  hb(k)
c       qpoly(k)=(b31-b32-spfvspoly-pfishvspoly)*dt+qpoly(k-1)
c       qbother(k)=(b41-b42-sdfvsother-pfishvsother)*dt+qbother(k-1)
c       qspfish(k)=(b51-b52-pfishvsspf)*dt+qspfish(k-1)
c       qsdfish(k)=(b61-b62-pfishvssdf)*dt+qsdfish(k-1)
c       qpfish(k)=(b71-b72)*dt+qpfish(k-1)
c
c     -- exchange of nutrients and oxygen --
c
        if(k.eq.nz) then
c
         wphsum=0.
         do 78 m=1,np
          wphsum=wphsum+wphy(m)*phy(m,k)
   78    continue
         wpcsum=wpoc*poc(k)
c
         qdip(k)=qdip(k)+.0d0*(dipph(1)*wphsum+dppomd*wpcsum)/ddz(k)
         qdinh(k)=qdinh(k)+.4d0*(dinph(1)*wphsum+dnpomd*wpcsum)/ddz(k)
	 qdox(k)=qdox(k)-topom*(wphsum+wpcsum)/ddz(k)
c
        endif
c
	 endif
c
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Vertical Mixing under Unstable Condition --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine mixcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- mean values of physical properties and state variables --
c
      do 10 k=2,nz-1
c
       if(dens(k-1).ge.dens(k)) then
c
        dzu=ddz(k-1)
        dzl=ddz(k)
c
        tmpm=(tmp(k)*dzl+tmp(k-1)*dzu)/(dzu+dzl)
        tmp(k)=tmpm
        tmp(k-1)=tmpm
c
        densm=1028.14-.0735d0*tmpm-.00469*tmpm*tmpm
     &               -35.*(.802-0.002*tmpm)
        dens(k)=densm
        dens(k-1)=densm
c
        do 12 m=1,np
         phym=(phy(m,k)*dzl+phy(m,k-1)*dzu)/(dzu+dzl)
         phy(m,k)=phym
         phy(m,k-1)=phym
   12   continue
c
        do 14 m=1,nzp
         zoom=(zoo(m,k)*dzl+zoo(m,k-1)*dzu)/(dzu+dzl)
         zoo(m,k)=zoom
         zoo(m,k-1)=zoom
   14   continue
c
c        bacm=(bac(k)*dzl+bac(k-1)*dzu)/(dzu+dzl)
c        bac(k)=bacm
c        bac(k-1)=bacm
c
        pocm=(poc(k)*dzl+poc(k-1)*dzu)/(dzu+dzl)
        poc(k)=pocm
        poc(k-1)=pocm
c
        docm=(doc(k)*dzl+doc(k-1)*dzu)/(dzu+dzl)
        doc(k)=docm
        doc(k-1)=docm
c
        dipm=(dip(k)*dzl+dip(k-1)*dzu)/(dzu+dzl)
        dip(k)=dipm
        dip(k-1)=dipm
c
        dinhm=(dinh(k)*dzl+dinh(k-1)*dzu)/(dzu+dzl)
        dinh(k)=dinhm
        dinh(k-1)=dinhm
c
c        dism=(dis(k)*dzl+dis(k-1)*dzu)/(dzu+dzl)
c        dis(k)=dism
c        dis(k-1)=dism
c
        doxm=(dox(k)*dzl+dox(k-1)*dzu)/(dzu+dzl)
        dox(k)=doxm
        dox(k-1)=doxm
c
        hbm=(hb(k)*dzl+hb(k-1)*dzu)/(dzu+dzl)
        hb(k)=hbm
        hb(k-1)=hbm
c
       endif
c
   10 continue
c
      return
      end
