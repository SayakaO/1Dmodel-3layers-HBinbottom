c-----------------------------------------------------------------------
c
c     -- Computational Conditions and Parameter Values --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine input
      implicit double precision (a-h,o-z)
      include 'param.h'
      character title*72
c
c     -- input and output title of numerical simulation --
c
      read(10,500) title
  500 format(a72)
c
c     -- computational condition --
c
      read(10,*) (nsw(i),i=1,7)
c
c     -- time step --
c
      read(10,*) dt,tmax
      read(10,*) intbal
c
c     -- interval of time histories output --
c
      read(10,*) tis
c
c     -- mesh --
c
      read(10,*) nz
      read(10,*) ddx,ddy
      read(10,*) (zb(k),k=1,nz)
c
c     -- river --
c
      read(10,*) nriv
      if(nriv.gt.0) then
       do 10 i=1,nriv
        read(10,*) iriv(i),jriv(i),qriv(i),triv(i),pcriv(i),
     &             dcriv(i),dpriv(i),dnhriv(i),dxriv(i)
   10  continue
      endif
c
c     -- meteorological boundary conditions for RTS --
c
      if(nsw(4).eq.1)then
       read(12,*) nta
       do 30 i=1,nta
        read(12,*) tat(i),tempa(i),pres(i),sunn(i),cld(i),hum(i),
     &             raind(i),wxt(i),wyt(i)
   30  continue
      endif
c
c     -- river boundary conditions for RTS --
c
      if(nsw(5).eq.1) then
      read(13,*) nra,nrb
      if(nra.ge.1) then
       do 40 j=1,nrb
        read(13,*) (qrvr(i,j),i=1,nra*7+1)
   40  continue
      endif
      endif
c
c     -- parameter values for momentum equations --
c
      dens0=1.d3
      pa0=1.01d5
      f=8.42d-5
      cp=3.93d6
      gamma=2.6d-3
      cd=1.5d-3
      densa=1.226d0
      am0=1.d0
      fkm0=1.d-2
      ac0=1.d0
      fkh0=1.d-2
c
c     -- parameter values for bulk formula --
c
      tmpa=27.d0
      qr0=341.8d0
      cloud=.5d0
      e_a=27.71d0
      q_s=2.23d-2
      q_e=1.73d-2
      rain=0.d0
      c_0=3.93d3
      c_a=1.006d3
      ref=9.d-2
      s_s=.96d0
      sigma=5.67051d-8
      clat=.65d0
      f_i=2.45d6
      ce=1.2d-3
      ch=1.2d-3
      ce2=1.2d-3
      ch2=1.2d-3
      wx0=0.d0
      wy0=0.d0
c
c     -- parameter values for mellor-yamada closure model --
c
	a1my=.92d0
	a2my=.74d0
	b1my=16.6d0
	b2my=10.1d0
	c1my=.08d0
	c2my=.7d0
	c3my=.2d0
	e1my=1.8d0
	e2my=1.33d0
	e3my=5.093d0
	ckar=.4d0
	fghc=-6.d0
c
c     -- initial values of physical properties  --
c
      tmp0=10.d0
	qs0=1.d-4
	qsl0=1.d-4
c
c     -- parameter values for species of plankton --
c
      np=1
      nzp=1
c
c     -- parameter values for state variables in the ecosystem --
c
      phy0(1)=20.d0
      zoo0(1)=2.d0
	  bac0=1000.d0
      poc0=100.d0
      doc0=1000.d0
      dip0=50.d0
      dinh0=400.d0
      dis0=300.d0
      dox0=10.d0
      hb0=20.d0
      poly0=20.d0
c      bother(1)=10000.d0
c      spfish(1)=10000.d0
c      sdfish(1)=10000.d0
c      pfish(1)=10000.d0
c
c     -- parameter values for phytoplankton --
c
      gp(1)=2.d0/86400.d0
      thep(1)=1.05d0
      aip(1)=1.d2
      ec0=.3d0
      ec1=.02d0
      hsp(1)=2.d0
      hsn(1)=25.d0
      erp(1)=.13d0
      gamp(1)=-8.44d-4
      rp(1)=.03d0/86400.d0
      dmp(1)=.0005d0/86400.d0
      wphy(1)=.1d0/86400.d0
      chlac(1)=.05d0
      dipph(1)=.05d0
      dinph(1)=.5d0
c      dipph(1)=.0244d0
c      dinph(1)=.176d0
      disph(1)=.15d0
      todph(1)=.00349d0
c
c     -- parameter values for zooplankton --
c
      cz(1)=.65d0/86400.d0
      thez(1)=1.05d0
      hsz(1)=0.d0
      rz(1)=.15d0/86400.d0
      az(1)=.6d0
      dmz(1)=.005d0/86400.d0
      dipzo(1)=.026d0
      dinzo(1)=.2d0
c      dipzo(1)=.0244d0
c      dinzo(1)=.176d0
      diszo(1)=.15d0
      todzo(1)=.00349d0
c
c     -- parameter values for particulate organic carbon --
c
      rpoc=.01d0/86400.d0
      thepoc=1.05d0
      ft=.25d0
      wpoc=.5d0/86400.d0
      dppom=.026d0
      dnpom=.2d0
c      dppom=.0244d0
c      dnpom=.176d0
      dspom=.15d0
      topom=.00349d0
c
c     -- parameter values for dissolved organic carbon --
c
      rdoc=.001d0/86400.d0
      thedoc=1.05d0
      dpdom=.025d0
      dndom=.2d0
c      dpdom=.0244d0
c      dndom=.176d0
      dsdom=.15d0
      todom=.00349d0
c
c     -- parameter values for aeration --
c
      rea=3.d0/86400.d0
c
      return
      end
