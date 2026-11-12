      block data ZZpolarization_data
c--- defaults for the WBFZZ helicity selection: off
      implicit none
      integer ZZpol_mode,ZZpol_hel(2),ZZpol_frame
      common/ZZpolarization/ZZpol_mode,ZZpol_hel,ZZpol_frame
      data ZZpol_mode/0/
      data ZZpol_hel/0,0/
      data ZZpol_frame/0/
      end


      subroutine getVVZZamps_pol(p,amp,ampa,ampb,za,zb,zab,zba,
     & j1,j2,j3,j4,j5,j6,j7,j8,doHO,doBO)
c--- Wrapper around getVVZZamps. With ZZpol_mode=0 it only calls
c--- getVVZZamps. With ZZpol_mode=1 it returns the amplitudes with the
c--- Z/gamma* -> (j3,j4) and Z/gamma* -> (j5,j6) projected onto the
c--- helicities chosen in ZZpol_hel (see ZZpolarization.f).
c---
c--- Method. In a diagram where the pair (j3,j4) comes from a single
c--- Z/gamma* propagator, the amplitude is M.L with the lepton current
c--- L^mu = <j3|gamma^mu|j4] and M depending on j3,j4 only through
c--- P = p(j3)+p(j4). Since P.L = 0,
c---     L^mu = - sum_lam eps^mu(lam) (eps*(lam).L),
c--- so the helicity-lam part is -(eps*(lam).L) M.eps(lam). M.eps(lam)
c--- is found by writing eps(lam) as a combination of the currents of
c--- three reference lepton configurations with the same P and calling
c--- getVVZZamps for those configurations. The same is done for
c--- (j5,j6), so up to 3x3 calls are needed. Diagrams where a lepton
c--- pair does not come from its own propagator are switched off via
c--- ZZpol_resonly while the polarized amplitudes are computed.
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'ZZpolarization.f'
      integer nmax
      parameter(nmax=10)
      integer j1,j2,j3,j4,j5,j6,j7,j8
      double precision p(mxpart,4)
      double complex zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & amp(nmax,2,2,2,2),ampa(nmax,2,2,2,2),ampb(nmax,2,2,2,2)
      logical doHO,doBO
      double precision prot(mxpart,4),kref(4),PA(4),PB(4),
     & qA(4,2,3),qB(4,2,3)
      double complex zaR(mxpart,mxpart),zbR(mxpart,mxpart),
     & zabR(mxpart,4,mxpart),zbaR(mxpart,4,mxpart),
     & ampR(nmax,2,2,2,2),ampaR(nmax,2,2,2,2),ampbR(nmax,2,2,2,2),
     & LactA(4,2),LactB(4,2),wA(3,2),wB(3,2),wt
      integer ncfgA,ncfgB,ia,ib,h3,h5,mu,nwarn
      logical okA,okB
      external ZZpolarization_data
      data nwarn/0/
      save nwarn

      ZZpol_resonly=.false.

      if (ZZpol_mode .eq. 0) then
        call getVVZZamps(amp,ampa,ampb,za,zb,zab,zba,
     &   j1,j2,j3,j4,j5,j6,j7,j8,doHO,doBO)
        return
      endif

      if ((ZZpol_mode .ne. 1)
     & .or. (ZZpol_hel(1) .lt. 0) .or. (ZZpol_hel(1) .gt. 5)
     & .or. (ZZpol_hel(2) .lt. 0) .or. (ZZpol_hel(2) .gt. 5)
     & .or. (ZZpol_frame .lt. 0) .or. (ZZpol_frame .gt. 2)) then
        write(6,*) 'getVVZZamps_pol: invalid ZZpol settings ',
     &   ZZpol_mode,ZZpol_hel(1),ZZpol_hel(2),ZZpol_frame
        stop
      endif

c--- pair momenta and the frame in which the helicity axes are defined
      do mu=1,4
        PA(mu)=p(j3,mu)+p(j4,mu)
        PB(mu)=p(j5,mu)+p(j6,mu)
        kref(mu)=0d0
      enddo
      if (ZZpol_frame .eq. 1) then
        do mu=1,4
          kref(mu)=PA(mu)+PB(mu)
        enddo
      elseif (ZZpol_frame .eq. 2) then
        do mu=1,4
          kref(mu)=-p(1,mu)-p(2,mu)
        enddo
      else
        kref(4)=1d0
      endif

c--- lepton currents of the input momenta;
c--- h=1: <j3|gamma|j4], h=2: <j4|gamma|j3] (likewise for j5,j6)
      do mu=1,4
        LactA(mu,1)=zab(j3,mu,j4)
        LactA(mu,2)=zab(j4,mu,j3)
        LactB(mu,1)=zab(j5,mu,j6)
        LactB(mu,2)=zab(j6,mu,j5)
      enddo

c--- reference configurations and projection weights for each pair
      call ZZpol_setup(p,j3,j4,ZZpol_hel(1),kref,LactA,
     & ncfgA,qA,wA,okA)
      call ZZpol_setup(p,j5,j6,ZZpol_hel(2),kref,LactB,
     & ncfgB,qB,wB,okB)

      amp(:,:,:,:,:)=czip
      ampa(:,:,:,:,:)=czip
      ampb(:,:,:,:,:)=czip

      if (okA .and. okB) then
        ZZpol_resonly=.true.
        do ia=1,ncfgA
        do ib=1,ncfgB
          prot(:,:)=p(:,:)
          do mu=1,4
            prot(j3,mu)=qA(mu,1,ia)
            prot(j4,mu)=qA(mu,2,ia)
            prot(j5,mu)=qB(mu,1,ib)
            prot(j6,mu)=qB(mu,2,ib)
          enddo
          call spinorcurr(8,prot,zaR,zbR,zabR,zbaR)
          call getVVZZamps(ampR,ampaR,ampbR,zaR,zbR,zabR,zbaR,
     &     j1,j2,j3,j4,j5,j6,j7,j8,doHO,doBO)
          do h3=1,2
          do h5=1,2
            wt=wA(ia,h3)*wB(ib,h5)
            amp(:,:,:,h3,h5)=amp(:,:,:,h3,h5)+wt*ampR(:,:,:,h3,h5)
            ampa(:,:,:,h3,h5)=ampa(:,:,:,h3,h5)+wt*ampaR(:,:,:,h3,h5)
            ampb(:,:,:,h3,h5)=ampb(:,:,:,h3,h5)+wt*ampbR(:,:,:,h3,h5)
          enddo
          enddo
        enddo
        enddo
        ZZpol_resonly=.false.
      else
        nwarn=nwarn+1
        if (nwarn .le. 10) then
          write(6,*) 'getVVZZamps_pol: no usable reference '//
     &     'configuration, amplitude set to zero for this point'
        endif
      endif

c--- restore the invariants s(i,j) of the input momenta, which
c--- spinorcurr has overwritten above
      call spinorcurr(8,p,zaR,zbR,zabR,zbaR)

      return
      end


      subroutine ZZpol_setup(p,n3,n4,isel,kref,Lact,ncfg,q,w,ok)
c--- For the lepton pair (n3,n4): reference configurations
c--- q(:,1:2,1:ncfg) (momenta of n3 and n4) and weights w(icfg,h)
c--- that project the pair onto the helicities selected by isel.
c--- h=1 refers to the current <n3|gamma|n4], h=2 to <n4|gamma|n3].
c--- isel=0 needs no projection: one configuration, the input one.
      implicit none
      include 'constants.f'
      integer n3,n4,isel,ncfg
      double precision p(mxpart,4),kref(4),q(4,2,3)
      double complex Lact(4,2),w(3,2)
      logical ok
      double precision ptmp(mxpart,4),Ppair(4),reldet,thr
      double complex zaT(mxpart,mxpart),zbT(mxpart,mxpart),
     & zabT(mxpart,4,mxpart),zbaT(mxpart,4,mxpart),
     & eps(4,3),Lcfg(4,3,2)
      integer iset,icfg,mu,k
      logical goodmom

      ok=.true.
      if (isel .eq. 0) then
        ncfg=1
        do mu=1,4
          q(mu,1,1)=p(n3,mu)
          q(mu,2,1)=p(n4,mu)
        enddo
        w(1,1)=cone
        w(1,2)=cone
        return
      endif

      ncfg=3
      do mu=1,4
        Ppair(mu)=p(n3,mu)+p(n4,mu)
      enddo
      call ZZpol_epsvec(Ppair,kref,eps)

c--- spinorcurr is singular for a momentum along -x, i.e. E+px=0.
c--- Reference leptons must not be much closer to -x than the input
c--- ones: (E+px)/E >= 1d-3 * min(1, value for the input leptons)
      thr=min(1d0,(p(n3,4)+p(n3,1))/p(n3,4),(p(n4,4)+p(n4,1))/p(n4,4))
      thr=1d-3*thr

c--- four sets of reference directions; a later set is used only if
c--- the previous one fails the test above or gives a badly
c--- conditioned system
      do iset=1,4
        goodmom=.true.
        do icfg=1,3
          call ZZpol_cfgmom(Ppair,iset,icfg,q(1,1,icfg),q(1,2,icfg))
          do k=1,2
            if (q(4,k,icfg)+q(1,k,icfg) .lt. thr*q(4,k,icfg)) then
              goodmom=.false.
            endif
          enddo
        enddo
        if (goodmom) then
          do icfg=1,3
            ptmp(:,:)=p(:,:)
            do mu=1,4
              ptmp(n3,mu)=q(mu,1,icfg)
              ptmp(n4,mu)=q(mu,2,icfg)
            enddo
            call spinorcurr(8,ptmp,zaT,zbT,zabT,zbaT)
            do mu=1,4
              Lcfg(mu,icfg,1)=zabT(n3,mu,n4)
              Lcfg(mu,icfg,2)=zabT(n4,mu,n3)
            enddo
          enddo
          call ZZpol_weights(Lact,Lcfg,eps,isel,w,reldet)
          if (reldet .gt. 1d-6) return
        endif
      enddo

      ok=.false.
      return
      end


      subroutine ZZpol_cfgmom(P,iset,icfg,q3,q4)
c--- massless q3, q4 with q3+q4 = P, back to back along +n and -n in
c--- the rest frame of P; n is direction icfg of the orthonormal set
c--- iset
      implicit none
      integer iset,icfg,mu
      double precision P(4),q3(4),q4(4),rm,dirs(3,3,4),r3(4),r4(4)
      data dirs/
     & 0.36d0,0.48d0,0.80d0,  -0.80d0,0.60d0,0.00d0,
     & -0.48d0,-0.64d0,0.60d0,
     & 0.60d0,0.00d0,0.80d0,  0.00d0,1.00d0,0.00d0,
     & -0.80d0,0.00d0,0.60d0,
     & 0.6666666666666667d0,0.6666666666666667d0,0.3333333333333333d0,
     & -0.6666666666666667d0,0.3333333333333333d0,0.6666666666666667d0,
     & 0.3333333333333333d0,-0.6666666666666667d0,0.6666666666666667d0,
     & 0.00d0,0.60d0,0.80d0,  0.80d0,-0.48d0,0.36d0,
     & 0.60d0,0.64d0,-0.48d0/
      save dirs

      rm=dsqrt(P(4)**2-P(1)**2-P(2)**2-P(3)**2)
      do mu=1,3
        r3(mu)=0.5d0*rm*dirs(mu,icfg,iset)
        r4(mu)=-r3(mu)
      enddo
      r3(4)=0.5d0*rm
      r4(4)=0.5d0*rm
      call ZZpol_rboost(P,r3,q3,-1)
      call ZZpol_rboost(P,r4,q4,-1)
      return
      end


      subroutine ZZpol_epsvec(P,kref,eps)
c--- polarization vectors of a vector boson with momentum P (P.P>0),
c--- helicity axis = direction of P in the rest frame of kref.
c--- eps(:,1): lam=-1, eps(:,2): lam=0, eps(:,3): lam=+1.
c--- Components (x,y,z,E), in the frame of the input momenta.
      implicit none
      include 'constants.f'
      double precision P(4),kref(4),PF(4),zp(3),xp(3),yp(3),ax(3),
     & pabs,rm,adotz,xnorm
      double complex eps(4,3),eF(4)
      integer mu,il,lam

      call ZZpol_rboost(kref,P,PF,1)
      rm=dsqrt(P(4)**2-P(1)**2-P(2)**2-P(3)**2)
      pabs=dsqrt(PF(1)**2+PF(2)**2+PF(3)**2)
      if (pabs .gt. 1d-9*dabs(PF(4))) then
        do mu=1,3
          zp(mu)=PF(mu)/pabs
        enddo
      else
        zp(1)=0d0
        zp(2)=0d0
        zp(3)=1d0
      endif

c--- transverse axes: xp in the plane of zp and the beam (z) axis
      ax(1)=0d0
      ax(2)=0d0
      ax(3)=1d0
      if (dabs(zp(3)) .ge. 0.99d0) then
        ax(1)=1d0
        ax(3)=0d0
      endif
      adotz=ax(1)*zp(1)+ax(2)*zp(2)+ax(3)*zp(3)
      do mu=1,3
        xp(mu)=ax(mu)-adotz*zp(mu)
      enddo
      xnorm=dsqrt(xp(1)**2+xp(2)**2+xp(3)**2)
      do mu=1,3
        xp(mu)=xp(mu)/xnorm
      enddo
      yp(1)=zp(2)*xp(3)-zp(3)*xp(2)
      yp(2)=zp(3)*xp(1)-zp(1)*xp(3)
      yp(3)=zp(1)*xp(2)-zp(2)*xp(1)

c--- longitudinal
      do mu=1,3
        eF(mu)=dcmplx(PF(4)/rm*zp(mu),0d0)
      enddo
      eF(4)=dcmplx(pabs/rm,0d0)
      call ZZpol_cboost(kref,eF,eps(1,2),-1)

c--- transverse: eps(lam) = -lam (xp + i lam yp)/sqrt(2)
      do il=1,3,2
        lam=il-2
        do mu=1,3
          eF(mu)=-dble(lam)*dcmplx(xp(mu),dble(lam)*yp(mu))/rt2
        enddo
        eF(4)=czip
        call ZZpol_cboost(kref,eF,eps(1,il),-1)
      enddo

      return
      end


      subroutine ZZpol_weights(Lact,Lcfg,eps,isel,w,reldet)
c--- For each lepton helicity h, writes eps(lam) = sum_a c(a) Lcfg(a)
c--- (solving the Gram system G c = r with G(a,b) = -Lcfg(a)*.Lcfg(b))
c--- and returns w(a,h) = sum_{lam kept} d(lam) c(a), with
c--- d(lam) = -eps*(lam).Lact(h). reldet is the smallest relative
c--- determinant of G met (1 = perfectly conditioned).
      implicit none
      include 'constants.f'
      integer isel,h,ia,ib,il
      double complex Lact(4,2),Lcfg(4,3,2),eps(4,3),w(3,2),
     & G(3,3),r(3),c(3),d,ZZpol_hdot
      double precision reldet,rd
      logical keep(3)
      external ZZpol_hdot

c--- keep(1): lam=-1, keep(2): lam=0, keep(3): lam=+1
      keep(1)=(isel.eq.2).or.(isel.eq.4).or.(isel.eq.5)
      keep(2)=(isel.eq.1).or.(isel.eq.5)
      keep(3)=(isel.eq.2).or.(isel.eq.3).or.(isel.eq.5)

      reldet=1d0
      do h=1,2
        do ia=1,3
          do ib=1,3
            G(ia,ib)=-ZZpol_hdot(Lcfg(1,ia,h),Lcfg(1,ib,h))
          enddo
          w(ia,h)=czip
        enddo
        do il=1,3
          if (keep(il)) then
            do ia=1,3
              r(ia)=-ZZpol_hdot(Lcfg(1,ia,h),eps(1,il))
            enddo
            call ZZpol_solve3(G,r,c,rd)
            d=-ZZpol_hdot(eps(1,il),Lact(1,h))
            do ia=1,3
              w(ia,h)=w(ia,h)+d*c(ia)
            enddo
            reldet=min(reldet,rd)
          endif
        enddo
      enddo

      return
      end


      subroutine ZZpol_solve3(G,r,c,rd)
c--- solves G c = r (3x3, complex) by Cramer's rule;
c--- rd = |det G|/|G(1,1) G(2,2) G(3,3)|, 0 if G is singular
      implicit none
      include 'constants.f'
      double complex G(3,3),r(3),c(3),M(3,3),det,ZZpol_det3
      double precision rd
      integer i,k
      external ZZpol_det3

      det=ZZpol_det3(G)
      if (abs(det) .eq. 0d0) then
        rd=0d0
        do k=1,3
          c(k)=czip
        enddo
        return
      endif
      rd=abs(det)/abs(G(1,1)*G(2,2)*G(3,3))
      do k=1,3
        do i=1,3
          M(i,1)=G(i,1)
          M(i,2)=G(i,2)
          M(i,3)=G(i,3)
        enddo
        do i=1,3
          M(i,k)=r(i)
        enddo
        c(k)=ZZpol_det3(M)/det
      enddo
      return
      end


      double complex function ZZpol_det3(M)
      implicit none
      double complex M(3,3)
      ZZpol_det3=M(1,1)*(M(2,2)*M(3,3)-M(2,3)*M(3,2))
     &          -M(1,2)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))
     &          +M(1,3)*(M(2,1)*M(3,2)-M(2,2)*M(3,1))
      return
      end


      double complex function ZZpol_hdot(a,b)
c--- conjg(a).b with metric (+,-,-,-), components (x,y,z,E)
      implicit none
      double complex a(4),b(4)
      ZZpol_hdot=dconjg(a(4))*b(4)-dconjg(a(1))*b(1)
     &          -dconjg(a(2))*b(2)-dconjg(a(3))*b(3)
      return
      end


      subroutine ZZpol_rboost(k,v,vout,idir)
c--- idir=+1: v (frame of the input momenta) -> rest frame of k;
c--- idir=-1: the inverse. v and vout must be different arrays.
      implicit none
      integer idir,mu
      double precision k(4),v(4),vout(4),bet(3),b2,gam,bv,km

      km=dsqrt(k(4)**2-k(1)**2-k(2)**2-k(3)**2)
      b2=0d0
      bv=0d0
      do mu=1,3
        bet(mu)=dble(idir)*k(mu)/k(4)
        b2=b2+bet(mu)**2
        bv=bv+bet(mu)*v(mu)
      enddo
      if (b2 .lt. 1d-30) then
        do mu=1,4
          vout(mu)=v(mu)
        enddo
        return
      endif
      gam=k(4)/km
      vout(4)=gam*(v(4)-bv)
      do mu=1,3
        vout(mu)=v(mu)+((gam-1d0)*bv/b2-gam*v(4))*bet(mu)
      enddo
      return
      end


      subroutine ZZpol_cboost(k,v,vout,idir)
c--- as ZZpol_rboost, for a complex vector v
      implicit none
      integer idir,mu
      double precision k(4),bet(3),b2,gam,km
      double complex v(4),vout(4),bv

      km=dsqrt(k(4)**2-k(1)**2-k(2)**2-k(3)**2)
      b2=0d0
      bv=dcmplx(0d0,0d0)
      do mu=1,3
        bet(mu)=dble(idir)*k(mu)/k(4)
        b2=b2+bet(mu)**2
        bv=bv+bet(mu)*v(mu)
      enddo
      if (b2 .lt. 1d-30) then
        do mu=1,4
          vout(mu)=v(mu)
        enddo
        return
      endif
      gam=k(4)/km
      vout(4)=gam*(v(4)-bv)
      do mu=1,3
        vout(mu)=v(mu)+((gam-1d0)*bv/b2-gam*v(4))*bet(mu)
      enddo
      return
      end
