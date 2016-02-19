c-----------------------------------------------------------------------
c
c     -- Grid Systems and Flags --
c                                               Author: Daisuke Kitazawa
c                                               Last Update: 2005.7.29
c
c-----------------------------------------------------------------------
      subroutine grid
      implicit double precision (a-h,o-z)
      include 'param.h'
c
c     -- grid system in the veritical direction --
c
      zb(0)=0.d0
      z(0)=0.d0
      do 10 k=1,nz
       z(k)=zb(k-1)+(zb(k)-zb(k-1))*.5d0
   10 continue
c
      do 12 k=1,nz
       ddz(k)=zb(k)-zb(k-1)
   12 continue
c
      dzh(0)=z(1)
      do 14 k=1,nz-1
       dzh(k)=z(k+1)-z(k)
   14 continue
      dzh(nz)=ddz(nz)*.5d0
c
      return
      end
