c-----------------------------------------------------------------------
c
c     -- Advection-Diffusion Equation --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine adcal(ww,eco,econ,qeco,fkhq)
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      dimension eco(0:nzmax+1)
      dimension econ(0:nzmax+1)
      dimension qeco(nzmax)
      dimension fkhq(0:nzmax)
c
      do 10 k=nz,1,-1
c
c	-- advection term in the z direction --
c
	 if(k.eq.1) then
	  advez=ww*eco(k)
	 else
        advez=ww*(eco(k)-eco(k-1))/ddz(k)
	 endif
c
c	-- diffusion term --
c
       ep=eco(k)
       eu=eco(k-1)
       el=eco(k+1)
c
	 if(k.eq.1) eu=ep
	 if(k.eq.nz) el=ep
c
	 fkhqu=fkhq(k-1)
	 fkhql=fkhq(k)
c
       if(k.eq.1) then
        fu=0.d0
       else
        fu=fkhqu*(ep-eu)/dzh(k-1)
       endif
c
	 if(k.eq.nz) then
	  fl=0.d0
	 else
        fl=fkhql*(el-ep)/dzh(k)
	 endif
c
       difez=(fl-fu)/ddz(k)
c
       econ(k)=eco(k)+tbal*dt*(-advez+difez+qeco(k))
c
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Advection-Diffusion Equation for Plankton --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine adcalf(m,ww,eco,econ,qeco,fkhq)
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      dimension eco(npmax,0:nzmax+1)
      dimension econ(npmax,0:nzmax+1)
      dimension qeco(npmax,nzmax)
	dimension fkhq(0:nzmax)
c
      do 10 k=nz,1,-1
c
c	-- advection term in the z direction --
c
	 if(k.eq.1) then
	  advez=ww*eco(m,k)
	 elseif(k.eq.nz) then
	  advez=ww*(eco(m,k)*.5d0-eco(m,k-1))/ddz(k)
	 else
        advez=ww*(eco(m,k)-eco(m,k-1))/ddz(k)
	 endif
c
c	-- diffusion term --
c
       ep=eco(m,k)
       eu=eco(m,k-1)
       el=eco(m,k+1)
c
	 if(k.eq.1) eu=ep
	 if(k.eq.nz) el=ep
c
	 fkhqu=fkhq(k-1)
	 fkhql=fkhq(k)
c
       if(k.eq.1) then
        fu=0.d0
       else
        fu=fkhqu*(ep-eu)/dzh(k-1)
       endif
c
	 if(k.eq.nz) then
	  fl=0.d0
	 else
        fl=fkhql*(el-ep)/dzh(k)
	 endif
c
       difez=(fl-fu)/ddz(k)
c
       econ(m,k)=eco(m,k)+tbal*dt*(-advez+difez+qeco(m,k))
c
   10 continue
      return
      end
c
c-----------------------------------------------------------------------
c     -- dynamics for HTL --
c                                               Author: Ayaka SAKAMOTO
c                                               Update: 2016.2.18
c-----------------------------------------------------------------------
      subroutine adcalb(ww,eco,econ,qeco,fkhq)
      implicit double precision (a-h,o-z)
      include 'param.h'

      dimension eco(0:nzmax+1)
      dimension econ(0:nzmax+1)
      dimension deco(0:nzmax)
      dimension fkhq(0:nzmax)
c
c      do 10 k=nz,1,-1
c
c     if(k.eq.nz) then
c
c      econ(k)=eco(k)+tbal*dt*deco(k)
c     end if 
c     10 continue
c      
c      return
c      end