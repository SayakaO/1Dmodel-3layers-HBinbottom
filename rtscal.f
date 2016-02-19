c-----------------------------------------------------------------------
c
c     -- Boundary Conditions for Real Time Simulation --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine rtscal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      if(time.gt.tat(nct(7))) nct(7)=nct(7)+1
c
c     -- wind velocity --
c
      wx=wxt(nct(7)-1)+(time-tat(nct(7)-1))
     &  *((wxt(nct(7))-wxt(nct(7)-1))
     &  /(tat(nct(7))-tat(nct(7)-1)))
      wy=wyt(nct(7)-1)+(time-tat(nct(7)-1))
     &  *((wyt(nct(7))-wyt(nct(7)-1))
     &  /(tat(nct(7))-tat(nct(7)-1)))
c
c     -- wind stress --
c
      uwx=wx-ui(1)
      vwx=wy-vj(1)
      wstx=densa*cd*uwx*dsqrt(uwx*uwx+vwx*vwx)
c
      uwy=wx-ui(1)
      vwy=wy-vj(1)
      wsty=densa*cd*vwy*dsqrt(uwy*uwy+vwy*vwy)
c
c     -- meteorological conditions --
c
      tmpaq=tempa(nct(7)-1)+(time-tat(nct(7)-1))
     &   *((tempa(nct(7))-tempa(nct(7)-1))
     &    /(tat(nct(7))-tat(nct(7)-1)))
      presq=pres(nct(7)-1)+(time-tat(nct(7)-1))
     &   *((pres(nct(7))-pres(nct(7)-1))
     &    /(tat(nct(7))-tat(nct(7)-1)))
      sunnq=sunn(nct(7)-1)+(time-tat(nct(7)-1))
     &   *((sunn(nct(7))-sunn(nct(7)-1))
     &    /(tat(nct(7))-tat(nct(7)-1)))
      cldq=cld(nct(7)-1)+(time-tat(nct(7)-1))
     &  *((cld(nct(7))-cld(nct(7)-1))
     &   /(tat(nct(7))-tat(nct(7)-1)))
      humq=hum(nct(7)-1)+(time-tat(nct(7)-1))
     &  *((hum(nct(7))-hum(nct(7)-1))
     &   /(tat(nct(7))-tat(nct(7)-1)))
      rainq=raind(nct(7)-1)+(time-tat(nct(7)-1))
     &   *((raind(nct(7))-raind(nct(7)-1))
     &    /(tat(nct(7))-tat(nct(7)-1)))
c
c     -- river conditions --
c
      if(nsw(5).eq.1) then
       if(time.gt.qrvr(1,nct(6))) nct(6)=nct(6)+1
       if(nra.ge.1) then
        do 40 n=1,nra
         qriv(n)=qrvr(n*7-5,nct(6)-1)
     &         +(time-qrvr(1,nct(6)-1))
     &        *((qrvr(n*7-5,nct(6))-qrvr(n*7-5,nct(6)-1))
     &         /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         triv(n)=qrvr(n*7-4,nct(6)-1)
     &         +(time-qrvr(1,nct(6)-1))
     &        *((qrvr(n*7-4,nct(6))-qrvr(n*7-4,nct(6)-1))
     &         /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         pcriv(n)=qrvr(n*7-3,nct(6)-1)
     &          +(time-qrvr(1,nct(6)-1))
     &         *((qrvr(n*7-3,nct(6))-qrvr(n*7-3,nct(6)-1))
     &          /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         dcriv(n)=qrvr(n*7-2,nct(6)-1)
     &          +(time-qrvr(1,nct(6)-1))
     &         *((qrvr(n*7-2,nct(6))-qrvr(n*7-2,nct(6)-1))
     &          /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         dpriv(n)=qrvr(n*7-1,nct(6)-1)
     &          +(time-qrvr(1,nct(6)-1))
     &         *((qrvr(n*7-1,nct(6))-qrvr(n*7-1,nct(6)-1))
     &          /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         dnhriv(n)=qrvr(n*7,nct(6)-1)
     &           +(time-qrvr(1,nct(6)-1))
     &          *((qrvr(n*7,nct(6))-qrvr(n*7,nct(6)-1))
     &           /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
         dxriv(n)=qrvr(n*7+1,nct(6)-1)
     &          +(time-qrvr(1,nct(6)-1))
     &         *((qrvr(n*7+1,nct(6))-qrvr(n*7+1,nct(6)-1))
     &          /(qrvr(1,nct(6))-qrvr(1,nct(6)-1)))
   40   continue
       endif
      endif
c
c     -- output of boundary conditions for RTS --
c
      if(mod(nt,100).eq.0) then
       write(61,100) time,wx,wy,tmpaq,presq,sunnq,cldq,humq,rainq
  100  format(f15.1,7f12.2,f15.5)
      endif
c
      if(mod(nt,100).eq.0) then
       write(62,200) time,(qriv(i),triv(i),pcriv(i),dcriv(i),
     &                     dpriv(i),dnhriv(i),
c     &                     dsriv(i),
     &                     dxriv(i),i=1,nra)
  200  format(f15.1,300f12.2)
      endif
c
      return
      end
