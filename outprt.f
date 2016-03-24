c-----------------------------------------------------------------------
c
c     -- Output Time Histories --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine outprt
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- flow velocities --
c
      write(21,100) time,(ui(k),k=1,nz)
      write(22,100) time,(vj(k),k=1,nz)
c
c     -- physical properties --
c
      if(nsw(2).eq.1) then
       write(23,100) time,(tmp(k),k=1,nz)
       write(24,100) time,(dens(k),k=1,nz)
       write(25,100) time,(fkmx(k),k=1,nz)
       write(26,100) time,(fkh(k),k=1,nz)
c
c     -- state variables in the ecosystem --
c
       if(nsw(6).eq.1) then
        write(27,100) time,((phy(m,k),k=1,nz),m=1,np)
        write(28,100) time,((zoo(m,k),k=1,nz),m=1,nzp)
        write(29,100) time,(poc(k),k=1,nz)
        write(30,100) time,(doc(k),k=1,nz)
        write(31,100) time,(dip(k),k=1,nz)
        write(32,100) time,(dinh(k),k=1,nz)
        write(33,100) time,(dox(k),k=1,nz)
        write(34,100) time,(hb(k),k=1,nz)
        write(35,100) time,(qphy(1,k),k=1,nz)
        write(50,100) time,(qbother(k),k=1,nz)
        write(70,100) time,(qspfish(k),k=1,nz)
        write(71,100) time,(qsdfish(k),k=1,nz)
        write(72,100) time,(qpfish(k),k=1,nz)
       endif
      endif
c
  100 format(f15.1,200e14.4)
c
      return
      end
