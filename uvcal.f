c-----------------------------------------------------------------------
c
c     -- Horizontal Flow Velocities --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine uvcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- vertical eddy viscosity --
c
	call kmcal
c
c     -- flow velocity in the x direction --
c
      call ucal
c
c     -- flow velocity in the y direction --
c
      call vcal
c
c     -- update horizontal flow velocities --
c
      do 10 k=1,nz
       ui(k)=uin(k)
       vj(k)=vjn(k)
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Vertical Eddy Viscosity Coefficients --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine kmcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- Km is constant --
c
	if(nsw(7).eq.0) then
c
c     -- vertical eddy viscosity coefficient in the x direction --
c
	 if(nsw(1).eq.1) then
	  do 10 k=0,nz
	   fkmx(k)=fkm0
   10   continue
       endif
c
c	-- Km depends on stratification function --
c
	elseif(nsw(7).eq.1) then
c
	 if(nsw(1).eq.1) then
c
        do 20 k=1,nz-1
c
c     -- vertical eddy viscosity coefficient in the x direction --
c
	   if(k.eq.0.or.k.eq.nz) then
	    fkmx(k)=0.d0
	    goto 20
	   endif
c
         rhodz=(dens(k+1)-dens(k))*.5d0/dzh(k)
         rho=(dens(k)*ddz(k+1)+dens(k+1)*ddz(k))/dzh(k)*.25d0
         udz=(ui(k+1)-ui(k))/dzh(k)
         vdz=(vj(k+1)-vj(k))/dzh(k)
         uvdz2=udz*udz+vdz*vdz
         epsi=btkm*g*rhodz/(rho*((fkmmin/fkm0)**(1.d0/alkm)-1.d0))
c
         if(rhodz.gt.1.0d-10) then
          if(uvdz2.le.epsi) then
           fkmx(k)=fkmmin
          else
           fkmx(k)=fkm0*(1.d0+btkm*g*rhodz/rho/uvdz2)**alkm
          endif
         endif
c
         if(fkmx(k).le.fkmmin) fkmx(k)=fkmmin
c
   20   continue
	 endif
	endif
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Flow Velocity in the X Direction --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine ucal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      do 10 k=nz,1,-1
c
c	-- Coriolis term --
c
       colx=f*vj(k)
c
c	-- eddy viscosity term --
c
       fkxu=fkmx(k-1)
       fkxl=fkmx(k)
c
	 up=ui(k)
	 vp=vj(k)
       uu=ui(k-1)
	 ul=ui(k+1)
	 if(k.eq.1) uu=up
	 if(k.eq.nz) ul=up
c
       quu=0.d0
       qul=0.d0
c
       if(k.eq.1) then
        quu=wstx/dens0
	  qul=fkxl*(ul-up)/dzh(k)
       elseif(k.eq.nz) then
	  quu=fkxu*(up-uu)/dzh(k-1)
        btmstx=-gamma*up*dens0*dsqrt(up*up+vp*vp)
	  qul=btmstx/dens0
	 else
	  quu=fkxu*(up-uu)/dzh(k-1)
	  qul=fkxl*(ul-up)/dzh(k)
       endif
c
       viscxz=(qul-quu)/ddz(k)
c
       uin(k)=up+dt*(colx+viscxz-up/864000.d0)
c
   10 continue
c
      return
      end
c-----------------------------------------------------------------------
c
c     -- Flow Velocity in the Y Direction --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine vcal
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      do 10 k=nz,1,-1
c
c	-- Coriolis term --
c
       coly=f*ui(k)
c
c	-- eddy viscosity term --
c
       fkyu=fkmx(k-1)
       fkyl=fkmx(k)
c
	 up=ui(k)
	 vp=vj(k)
       vu=vj(k-1)
	 vl=vj(k+1)
	 if(k.eq.1) vu=vp
	 if(k.eq.nz) vl=vp
c
       qvu=0.d0
       qvl=0.d0
c
       if(k.eq.1) then
        qvu=wsty/dens0
	  qvl=fkyl*(vl-vp)/dzh(k)
       elseif(k.eq.nz) then
	  qvu=fkyu*(vp-vu)/dzh(k-1)
        btmsty=-gamma*vp*dens0*dsqrt(vp*vp+up*up)
	  qvl=btmsty/dens0
	 else
	  qvu=fkyu*(vp-vu)/dzh(k-1)
	  qvl=fkyl*(vl-vp)/dzh(k)
       endif
c
       viscyz=(qvl-qvu)/ddz(k)
c
       vjn(k)=vp+dt*(-coly+viscyz-vp/864000.d0)
c
   10 continue
c
      return
      end
