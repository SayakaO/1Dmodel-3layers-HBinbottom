c-----------------------------------------------------------------------
c
c     -- Initial Conditions --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine init
      implicit double precision (a-h,o-z)
      include 'param.h'
c
      ntmax=idint(tmax/dt)
      ntis=idint(tis/dt)
      ntia=idint(tia/dt)
      nstrt=idint(strt/dt)
c
c     -- initial values of flow velocities --
c
      do 20 k=0,nz
       ui(k)=0.d0
       uin(k)=0.d0
   20 continue
c
      do 30 k=0,nz
       vj(k)=0.d0
       vjn(k)=0.d0
   30 continue
c
c     -- initial values of physical properties and state variables --
c
      do 60 k=0,nz+1
       dens(k)=dens0
       tmp(k)=tmp0
       tmpn(k)=0.d0
       do 62 m=1,np
        phy(m,k)=phy0(m)
        phyn(m,k)=0.d0
   62  continue
       do 66 m=1,nzp
        zoo(m,k)=zoo0(m)
        zoon(m,k)=0.d0
   66  continue
       bac(k)=bac0
       bacn(k)=0.d0
       poc(k)=poc0
       pocn(k)=0.d0
       doc(k)=doc0
       docn(k)=0.d0
       dip(k)=dip0
       dipn(k)=0.d0
       dinh(k)=dinh0
       dinhn(k)=0.d0
       dis(k)=dis0
       disn(k)=0.d0
       dox(k)=dox0
       doxn(k)=0.d0
       hb(k)=hb0
       hbn(k)=0.d0
   60 continue
c
      wx=wx0
      wy=wy0
c
c     -- initial values of eddy viscosity coefficients --
c
	do 120 k=0,nz
	 fkmx(k)=fkm0
	 fkh(k)=fkh0
  120 continue
c
      return
      end
