      parameter(pi=3.14159265d0)
      parameter(g=9.80665d0)
      parameter(nzmax=100)
      parameter(nrmax=30000)
      parameter(lbmax=10)
      parameter(npmax=1)
      parameter(nrvmax=1)
      parameter(mmmp=11302)
      parameter(mmmb=738)
      parameter(alkm=-1.,btkm=5.2)
      parameter(fkmmin=1.049d-6)
      parameter(fkmmax=1.d-2)
      parameter(btkc=10./3.,alkc=-3./2.)
      parameter(fkhmin=1.4d-7,fkqmin=1.4d-7)
      parameter(fkhmax=1.d-2,fkqmax=1.d-2)
      parameter(rkap=.5)
c
      common /times/dt,tmax,time,tbal,intbal
      common /count/nt,ntmax,nout,nrf,nct(12),nsw(7)
      common /coord/nx,ny,nz
      common /grids/ddx,ddy,ddz(nzmax),
     &              z(0:nzmax),zb(0:nzmax),dzh(0:nzmax)
      common /lb/nxlb,nylb,ixlb(lbmax),jxlb(lbmax),
     &           iylb(lbmax),jylb(lbmax)
      common /vphy/f,cp,pa0,densa,wx0,wy0,
     &             am0,ac0,fkm0,fkh0,fkq0,gamma,cd
      common /bulk/tmpa,qr0,cloud,e_a,q_s,q_e,rain,
     &             c_0,c_a,ref,s_s,sigma,clat,f_i,ce,ch,ce2,ch2
      common /mymodel/a1my,a2my,b1my,b2my,c1my,c2my,c3my,e1my,e2my,e3my,
     &                ckar,fghc
      common /vnos/np,nzp
      common /vphy/gp(npmax),thep(npmax),aip(npmax),ec,ec0,ec1,
     &             hsp(npmax),hsn(npmax),erp(npmax),gamp(npmax),
     &             rp(npmax),dmp(npmax),wphy(npmax),chlac(npmax),
     &             dipph(npmax),dinph(npmax),disph(npmax),todph(npmax)
      common /vzoo/cz(npmax),thez(npmax),hsz(npmax),rz(npmax),
     &             az(npmax),dmz(npmax),dipzo(npmax),
     &             dinzo(npmax),diszo(npmax),todzo(npmax)
      common /vpoc/rpoc,thepoc,ft,wpoc,dppom,dnpom,dspom,topom
      common /vdoc/rdoc,thedoc,dpdom,dndom,dsdom,todom
      common /vdox/rea
      common /rivn/nriv,iriv(nrvmax),jriv(nrvmax)
      common /rivp/qriv(nrvmax),triv(nrvmax)
      common /rive/pcriv(nrvmax),dcriv(nrvmax),dpriv(nrvmax),
     &             dnhriv(nrvmax),dsriv(nrvmax),dxriv(nrvmax),
     &             nra,nrb,qrvr(10,150)
      common /windt/wxt(nrmax),wyt(nrmax),nw
      common /wind/wx,wy,wstx,wsty
      common /ini/dens0,tmp0,qs0,qsl0,phy0(npmax),zoo0(npmax),
     &            bac0,poc0,doc0,dip0,dinh0,dis0,dox0
      common /rts/nta,nba,
     &            tat(nrmax),tempa(nrmax),
     &            pres(nrmax),sunn(nrmax),cld(nrmax),
     &            hum(nrmax),raind(nrmax),
     &            tmpaq,presq,sunnq,cldq,humq,rainq
      common /vel/ui(0:nzmax+1),uin(0:nzmax+1),
     &            vj(0:nzmax+1),vjn(0:nzmax+1)
      common /kmhq/fkmx(0:nzmax),fkh(0:nzmax)
      common /dens/dens(0:nzmax+1)
      common /tmp/tmp(0:nzmax+1),tmpn(0:nzmax+1),qtmp(nzmax)
      common /phy/phy(npmax,0:nzmax+1),phyn(npmax,0:nzmax+1),
     &            qphy(npmax,nzmax)
      common /zoo/zoo(npmax,0:nzmax+1),zoon(npmax,0:nzmax+1),
     &            qzoo(npmax,nzmax)
      common /bac/bac(0:nzmax+1),bacn(0:nzmax+1),qbac(nzmax)
      common /poc/poc(0:nzmax+1),pocn(0:nzmax+1),qpoc(nzmax)
      common /doc/doc(0:nzmax+1),docn(0:nzmax+1),qdoc(nzmax)
      common /dip/dip(0:nzmax+1),dipn(0:nzmax+1),qdip(nzmax)
      common /dinh/dinh(0:nzmax+1),dinhn(0:nzmax+1),qdinh(nzmax)
      common /dis/dis(0:nzmax+1),disn(0:nzmax+1),qdis(nzmax)
      common /dox/dox(0:nzmax+1),doxn(0:nzmax+1),qdox(nzmax)
      common /outv/nv,iv(20),jv(20)
      common /outz/nzt,iz(20),jz(20)
      common /outd/nd,id(20),jd(20)
      common /outnt/ntout(15),ntis,ntia,nstrt
      common /outt/tout(15),tis,tia,strt,nav
      common /outnrf/ntrf(500),nprd(500)
      common /outrf/trf(500),prd(500)
      common /outnav/ntav(500),npav(500)
      common /outav/tav(500),pav(500)
      common /hb/hb(0:nzmax+1),hbn(0:nzmax+1),qhb(0:nzmax)
      common /poly/poly(0:nzmax+1),polyn(0:nzmax+1),qpoly(nzmax)
      common /bother/bother(0:nzmax+1),bothern(0:nzmax+1),qbother(nzmax)
      common /spfish/spfish(0:nzmax+1),spfishn(0:nzmax+1),qspfish(nzmax)
      common /sdfish/sdfish(0:nzmax+1),sdfishn(0:nzmax+1),qsdfish(nzmax)
      common /pfish/pfish(0:nzmax+1),pfishn(0:nzmax+1),qpfish(nzmax)