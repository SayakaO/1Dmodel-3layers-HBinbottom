c-----------------------------------------------------------------------
c
c     -- Numerical Simulation Program of the Ecosystem in Lake Biwa --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      program main
      implicit double precision (a-h,o-z)
      include 'param.h'
      character title*72
c
c     -- open files for computational conditions --
c
      open(10,file='input/eco.dat',form='formatted',status='unknown')
      open(12,file='input/rts.dat',form='formatted',status='unknown')
      open(13,file='input/riv.dat',form='formatted',status='unknown')
c
c     -- input computational conditions and parameter values --
c
      call input
c
      write(*,*) 'Numerical simulation was started successfully!!' 
c
c     -- close files for computational conditions --
c
      close(10)
      close(12)
      close(13)
c
c     -- grid systems and flags --
c
      call grid
c
c     -- initial conditions --
c
      call init
c
c     -- open files for output of time histories --
c
      open(21,file='output/tu.dat',form='formatted',
     &        status='unknown')
      open(22,file='output/tv.dat',form='formatted',
     &        status='unknown')
      if(nsw(2).eq.1) then
       open(23,file='output/tt.dat',form='formatted',
     &         status='unknown')
       open(24,file='output/td.dat',form='formatted',
     &         status='unknown')
       open(25,file='output/tkm.dat',form='formatted',
     &         status='unknown')
       open(26,file='output/tkh.dat',form='formatted',
     &         status='unknown')
       if(nsw(6).eq.1) then
        open(27,file='output/tph.dat',form='formatted',
     &          status='unknown')
        open(28,file='output/tzo.dat',form='formatted',
     &          status='unknown')
        open(29,file='output/tpc.dat',form='formatted',
     &          status='unknown')
        open(30,file='output/tdc.dat',form='formatted',
     &          status='unknown')
        open(31,file='output/tdp.dat',form='formatted',
     &          status='unknown')
        open(32,file='output/tdn.dat',form='formatted',
     &          status='unknown')
        open(33,file='output/tdo.dat',form='formatted',
     &          status='unknown')
        open(34,file='output/thb.dat',form='formatted',
     &          status='unknown')
        open(35,file='output/tpoly.dat',form='formatted',
     &          status='unknown')
        open(50,file='output/tother.dat',form='formatted',
     &          status='unknown')
        open(70,file='output/tspfish.dat',form='formatted',
     &          status='unknown')  
        open(71,file='output/tsdfish.dat',form='formatted',
     &          status='unknown') 
        open(72,file='output/tpfish.dat',form='formatted',
     &          status='unknown') 
        open(73,file='output/tqpoc(pos).dat',form='formatted',
     &          status='unknown') 
        open(74,file='output/tqpoc(neg).dat',form='formatted',
     &          status='unknown') 
       endif
      endif
c
c     -- open files for output of computational conditions for RTS --
c
      if(nsw(4).eq.1) then
       open(61,file='output/rtst.dat',status='unknown')
       open(62,file='output/rivt.dat',status='unknown')
      endif
c
c     -- initial values of counter --
c
      nt=0
      time=0.d0
      tbal=real(intbal)
      do 10 i=1,3
       nct(i)=0
   10 continue
      do 20 i=4,12
       nct(i)=1
   20 continue
c
c-----------------------------------------------------------------------
c
  100 continue
c
      nt=nt+1
      time=time+dt
c
      if(mod(nt,2880).eq.0) then
       write(*,200) nt,'(step)',time/86400.d0,'(day)'
  200  format(i10,a7,f10.2,a6)
      endif
c
c     -- boundary conditions for RTS --
c
      if(nsw(4).eq.1) call rtscal
c
      if(nsw(1).eq.1) then
c
c     -- horizontal flow velocities --
c
       call uvcal
c
      endif
c
c     -- physical properties and state variables in the ecosystem --
c
      if(nsw(2).eq.1) then
c       if(mod(nt,intbal).eq.0) call ecocal
       call ecocal
      endif
c
c     -- output time histories --
c
      if(ntis.ge.1) then
       if(mod(nt,ntis).eq.0) then
        call outprt
       endif
      endif
c
      if(nt.lt.ntmax) goto 100
c
c-----------------------------------------------------------------------
c
c
c     -- output time histories at the last time step --
c
      if(mod(nt,ntis).ne.0) then
       call outprt
      endif
c
c     -- close files --
c
      close(31)
      close(36)
      if(nsw(2).eq.1) then
       close(41)
       if(nsw(6).eq.1) then
        close(46)
       endif
      endif
c
      if(nsw(4).eq.1) then
       close(61)
       close(62)
	endif
c
      close(63)
      close(64)
      close(65)
c
      write(*,*) 'Numerical simulation was completed successfully!!'
c
  999 stop 
      end
