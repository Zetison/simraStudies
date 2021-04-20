!///////////////////////////////////////////////////////////////////
!  ~/Simra/Prepro/ ibc.f90
!
!    This program is composed of ibc_UM2Simra & setIBC, 
!    to calculate simplified IBC's for SIMRA
!...................................................................
!                                                                  
!     Purpose:                                                      
!                 to calculate IBC data for SIMRA                  
!                 for test case purposes
!
!     Input :     Data from SIMRA:  parm.dat, mesh.dat             
!                                                                  
!     Output:     init-data to 'init.dat'                          
!                 BC  data to  'boun.dat'                          
!                                                                  
!     NOTES:   
!              *  All input data in SI units                       
!              *  Automatic change of Dirichlet BCs at (E,W,S,N)           
!                   according to the actual flow situation           
!              *  Structured mesh numbering is assumed
!              *  Initial turbulence is calculated from eddy model,        
!                   use fixed BCs for turbulence along inflow sides             
!                                                                  
!///////////////////////////////////////////////////////////////////

      module linalg

          implicit none
          
          contains
          
          FUNCTION cross(a, b)
            real, DIMENSION(3) :: cross
            real, DIMENSION(3), INTENT(IN) :: a, b
          
            cross(1) = a(2) * b(3) - a(3) * b(2)
            cross(2) = a(3) * b(1) - a(1) * b(3)
            cross(3) = a(1) * b(2) - a(2) * b(1)
          END FUNCTION cross
        
      end module linalg
      use linalg
      implicit none
      integer nstep,iout,istrat,niter,npoin,nelem,imax,jmax,kmax,   &
             nelM,i,ip,ie,numIC,j,itop,k,i3,ib,ib2,iu,iv,iw,ipp,ik, &
             ilog,istr,icont,iew,new,i3w,i3n,nlog,itest_S,itest_N,  &
             i3p,i3m,itest_E,itest_W,nfixp,nfixu,nfixv,nfixw,nfixk, &
             nfixe,ikke,i1,kont,isolve,iturb,mp,mb,m1
      integer ipN,ipS,ipW,ipE,mxy
      real dt,eps,omega,visc,rho0,hscal,uinf,cappa,z0,delta,       &
           ustar,u_top,tkmin,z_bakke,z,eta,rlmix,ww,r2,rad,        &
           uinf2,uinf3,time
      logical isinflow

      integer, allocatable:: kon_we(:),icod_b(:),icod_n(:),icu(:),  &
                             icv(:),icw(:),ick(:),icp(:),ice(:)
      real, allocatable:: pd(:),ps(:),rhos(:),rho(:),u1(:),u2(:),    &
       u3(:),tk(:),td(:),vtef(:),pt(:),pts(:),coord(:,:),fixu(:),    &
       fixv(:),fixw(:),fixk(:),fixd(:),fixe(:),u1m(:),u2m(:),u3m(:), &
       tleng(:),s1(:),s2(:),s3(:),dimas(:),detj(:),epk(:),dist(:),z0L(:)
      real, allocatable:: z1(:),rho1(:),tk1(:),ep1(:),av1(:),u(:),v(:)
      real, DIMENSION(3,4) :: normals
      real, DIMENSION(3) :: PswBot,PswTop,PseBot,PseTop,            &
                            PnwBot,PnwTop,PneBot,PneTop
!-------------------------------------------------------------------
!
! ----- Read parameters for SIMRA, SI units
!
      open (unit=9,file='parm.dat',status='old')
      read(9,*) isolve,kont,iout,istrat
      read(9,*) nstep,dt
      read(9,*) niter,eps,omega
      read(9,*) visc,rho0
      read(9,*) hscal,uinf
!
! ----- Read grid data for SIMRA, SI units
!
      open (unit=10,file='mesh.dat',form='unformatted')
      read(10) npoin,nelem,imax,jmax,kmax,nelM
      allocate (coord(3,npoin))
      read(10) ((coord(i,ip),i=1,3),ip=1,npoin)
      !read(10) ((lnode(i,ie),i=1,8),ie=1,nelem)
      close(10)
      print*,' mesh: (imax,jmax,kmax)=',imax,jmax,kmax

      mp=npoin
      mb=imax*jmax*2+imax*kmax*2+jmax*kmax*2
      m1=max(imax,jmax)
      allocate (kon_we(mb),icod_b(mb),icod_n(mb),icu(mb),icv(mb),  &
                icw(mb),ick(mb),icp(mb),ice(mb))
      allocate (pd(mp),ps(mp),rhos(mp),rho(mp),u1(mp),u2(mp),u3(mp), &
        tk(mp),td(mp),vtef(mp),pt(mp),pts(mp),dimas(mp), &
        detj(mp),epk(mp),dist(mp),u1m(mp),u2m(mp),u3m(mp),tleng(mp))
      allocate (fixu(mb),fixv(mb),fixw(mb),fixk(mb),fixd(mb), &
                fixe(mb),s1(mb),s2(mb),s3(mb))
      allocate (z1(m1),rho1(m1),tk1(m1),ep1(m1),av1(m1),u(m1),v(m1))
      mxy=imax*jmax
      allocate (z0L(mxy))
!-------------------------------------------------------------------
!
! ----- Compute vertical profile for (ui, pts) and distribute
!       as initial values
!
      cappa=0.42
      z0=0.2 !0.3              ! (m) ad hoc      
      delta=1500. !2000.       ! (m) assumed BL height 
      delta=delta+z0

      ! Wind direction is set in subroutine profile
      call profile(coord,imax,jmax,kmax,delta,ustar,z0,u1,u2,u3,pts, &
                   istrat,mp,m1)
!
! ---- calculate hydrostatic pressure & density from State Equation
!
      call hstat (pts,ps,rhos,coord,imax,jmax,kmax,mp,mb)
!     pts=state eq! must be included here
!      ...and set initial density equal to hydrostatic value
!         same with pot.temp
      do ip=1,npoin
        rho(ip)=rhos(ip)
        pt(ip)=pts(ip)   
      enddo
!-------------------------------------------------------------------
!
! ----- Compute initial turbulence from simplified, 1D eddy model
!       along each quasi-vertical line (subroutine eddy)
!

      do i=1,imax
        do j=1,jmax
          do k=1,kmax
              i3=k+(i-1)*kmax+(j-1)*imax*kmax
              z1(k)=coord(3,i3)
              rho1(k)=rho(i3)
              u(k)=u1(i3)
              v(k)=u2(i3)
          enddo

          CALL eddy (kmax,delta,ustar,z1,rho1,u,v,tk1,ep1,av1,z0,m1)
          !if (i==1.and.j==1) then   !!
          !   do k=1,kmax
          !     write(31,*) tk1(k),z1(k)   !!
          !     write(32,*) sqrt(u(k)**2+v(k)**2),z1(k)   !!
          !   enddo
          !endif
          do k=1,kmax
            i3=k+(i-1)*kmax+(j-1)*imax*kmax
            vtef(i3)=av1(k)
            tk(i3)=tk1(k)
            td(i3)=ep1(k)
          enddo            
        enddo
      enddo
!-------------------------------------------------------------------
!
! ----- Smoothing of potential temperature along ground 
!        (usually unneccessary, only if strange data'holes')
!
      call mean (imax,jmax,kmax,pt,pt,1,mp)
!********************************************************************
!
! ----- BOUNDARY CONDITIONS              
!
!********************************************************************
!
      iu=0
      iv=0
      iw=0
      ipp=0
      ik=0
      ilog=0
      istr=0

      !print*,'  iturb = 1: inflow Dirichlet BC for (K,eps)'
      !print*,'  iturb = 0: complete free side BCs for (K,eps)'
      !read*, iturb
      iturb=1
!--------------------------------------------------------------
!
!  GROUND LEVEL (k=1) :  LOG-BC's along ground/terrain surface
!
!
!  ----- wall elements:
!
      icont=0
      do j=1,jmax-1
        do i=1,imax-1
          icont=icont+1
          iew=1+(kmax-1)*(i-1)+(imax-1)*(kmax-1)*(j-1)
          kon_we(icont)=iew
        enddo
      enddo
      new=icont
!
!  ----- wall & near-wall nodes:
!
        z0L=z0   ! land roughness

      do j=1,jmax           
      do i=1,imax           
        i3w=1+(i-1)*kmax+(j-1)*imax*kmax 
        i3n=2+(i-1)*kmax+(j-1)*imax*kmax 
        ilog=ilog+1
        icod_b(ilog)=i3w  !  i.e. point 1 -> on ground
        icod_n(ilog)=i3n  !  i.e. point 2 -> near-ground
          if (coord(3,i3w) <= 10.0) z0L(ilog)=0.001  ! sea roughness
      enddo
      enddo
      nlog=ilog
!
!  ----- Include LOG-points as Dirichlet conditions:
!
      do ib=1,nlog
        i3=icod_b(ib)
        iu=iu+1
        iv=iv+1
        iw=iw+1
        ik=ik+1
        istr=istr+1
        icu(iu)=i3
        icv(iv)=i3
        icw(iw)=i3
        ick(ik)=icod_n(ib)  ! point 2, above ground for K,eps !
        ice(istr)=i3
        fixu(iu)=0.0        ! no-slip
        fixv(iv)=0.0        !  "
        fixw(iw)=0.0        !  "
        fixk(ik)=0.0        ! computed in SIMRA 
        fixd(ik)=0.0        !    "
        fixe(istr)=pt(i3)   ! temp. along ground, from UM 
      end do
!--------------------------------------------------------------
! ---- Dirichlet BC if inflow: dot(u,n) < 0 
! The four boundaries are assumed to be planar    
      PswBot=coord(:,1)
      PswTop=coord(:,kmax)
      PseBot=coord(:,(imax-1)*kmax+1)
      PseTop=coord(:,kmax*imax)
      PnwBot=coord(:,(jmax-1)*imax*kmax+1)
      PnwTop=coord(:,(jmax-1)*imax*kmax+kmax)
      PneBot=coord(:,(jmax-1)*imax*kmax+(imax-1)*kmax+1)
      PneTop=coord(:,jmax*imax*kmax)
      normals(:,1)=cross(PseBot-PswBot,PswTop-PswBot)
      normals(:,2)=cross(PnwBot-PneBot,PneTop-PneBot)
      normals(:,3)=cross(PneBot-PseBot,PseTop-PseBot)
      normals(:,4)=cross(PswBot-PnwBot,PnwTop-PnwBot)
!SOU THERN BOUNDARY (j=1):  conditional
      itest_S=1
      if (itest_S.eq.1) then
          j=1
          do i=1,imax
              do k=1,kmax
                  i3=k+(i-1)*kmax+(j-1)*imax*kmax
                  call checkFlowDir (isInflow,u1(i3),u2(i3),u3(i3), &
                      normals(:,1))
                  if (isInflow) then
                       if (k.eq.kmax.and.i.eq.1) & 
                           print*, 'wind from bottom boundary'
                       iu=iu+1
                       iv=iv+1
                       iw=iw+1
                       istr=istr+1
                       icu(iu)=i3
                       icv(iv)=i3
                       icw(iw)=i3
                       ice(istr)=i3
                       fixu(iu)=u1(i3)
                       fixv(iv)=u2(i3)
                       fixw(iw)=u3(i3)
                       fixe(istr)=pt(i3)
                       if (iturb.eq.1) then
                           ik=ik+1
                           ick(ik)=i3
                           fixk(ik)=tk(i3)
                           fixd(ik)=td(i3)
                       endif
                  endif
              enddo
          enddo
      endif
!NOR THERN BOUNDARY (j=jmax):  conditional 
      itest_N=1
      if (itest_N.eq.1) then
          j=jmax
          do i=1,imax
              do k=1,kmax
                  i3=k+(i-1)*kmax+(j-1)*imax*kmax
                  call checkFlowDir (isInflow,u1(i3),u2(i3),u3(i3), &
                      normals(:,2))
                  if (isInflow) then
                      if (k.eq.kmax.and.i.eq.1) &
                          print*, 'wind from top boundary'
                      iu=iu+1
                      iv=iv+1
                      iw=iw+1
                      istr=istr+1
                      icu(iu)=i3
                      icv(iv)=i3
                      icw(iw)=i3
                      ice(istr)=i3
                      fixu(iu)=u1(i3)
                      fixv(iv)=u2(i3)
                      fixw(iw)=u3(i3)
                      fixe(istr)=pt(i3)
                      if (iturb.eq.1) then
                          ik=ik+1
                          ick(ik)=i3
                          fixk(ik)=tk(i3)
                          fixd(ik)=td(i3)
                      endif
                  endif
              enddo
          enddo
      endif
! EA STERN BOUNDARY (i=imax):  conditional
      itest_E=1  !
      if (itest_E.eq.1) then
          i=imax
          do j=1,jmax
              do k=1,kmax
                  i3=k+kmax*(i-1)+kmax*imax*(j-1)
                  call checkFlowDir (isInflow,u1(i3),u2(i3),u3(i3), &
                      normals(:,3))
                  if (isInflow) then
                      if (k.eq.kmax.and.j.eq.1) & 
                          print*, 'wind from right boundary'
                      iu=iu+1
                      iv=iv+1
                      iw=iw+1
                      istr=istr+1
                      icu(iu)=i3
                      icv(iv)=i3
                      icw(iw)=i3
                      ice(istr)=i3
                      fixu(iu)=u1(i3)
                      fixv(iv)=u2(i3)
                      fixw(iw)=u3(i3)
                      fixe(istr)=pt(i3)
                      if (iturb.eq.1) then
                          ik=ik+1
                          ick(ik)=i3
                          fixk(ik)=tk(i3)
                          fixd(ik)=td(i3)
                      endif
                  endif
              enddo
          enddo
      endif
!  W ESTERN BOUNDARY (i=1):   Conditional     
      itest_W=1
      if (itest_W.eq.1) then
          do j=1,jmax
              do k=1,kmax
                  i3=k+kmax*imax*(j-1)
                  call checkFlowDir (isInflow,u1(i3),u2(i3),u3(i3), &
                      normals(:,4))
                  if (isInflow) then
                      if (k.eq.kmax.and.j.eq.1) &
                          print*, 'wind from left boundary'
                      iu=iu+1
                      iv=iv+1
                      iw=iw+1
                      istr=istr+1
                      icu(iu)=i3
                      icv(iv)=i3
                      icw(iw)=i3
                      ice(istr)=i3
                      fixu(iu)=u1(i3)
                      fixv(iv)=u2(i3)
                      fixw(iw)=u3(i3)
                      fixe(istr)=pt(i3)
                      if (iturb.eq.1) then
                          ik=ik+1
                          ick(ik)=i3
                          fixk(ik)=tk(i3)
                          fixd(ik)=td(i3)
                      endif
                  endif
              enddo
          enddo
      endif
!!--------------------------------------------------------------
!!
!!  SOUTHERN BOUNDARY (j=1):  conditional
!!
!      itest_S=1
!      if (itest_S.eq.1) then
!      j=1
!      do i=1,imax
!      do k=1,kmax
!        i3=k+(i-1)*kmax+(j-1)*imax*kmax
!!
!! ---- Dirichlet BC if inflow: u2 > 0 
!!
!        if (u2(i3).gt.0.0) then
!          if (k.eq.1.and.i.eq.1) print*, 'wind from S'
!          iu=iu+1
!          iv=iv+1
!          iw=iw+1
!          istr=istr+1
!          icu(iu)=i3
!          icv(iv)=i3
!          icw(iw)=i3
!          ice(istr)=i3
!          fixu(iu)=u1(i3)   
!          fixv(iv)=u2(i3)   
!          fixw(iw)=u3(i3)   
!          fixe(istr)=pt(i3)  
!          if (iturb.eq.1) then
!            ik=ik+1
!            ick(ik)=i3
!            fixk(ik)=tk(i3)
!            fixd(ik)=td(i3)
!          endif
!        endif
!      enddo
!      enddo
!
!      endif
!!--------------------------------------------------------------
!!
!!  NORTHERN BOUNDARY (j=jmax):  conditional 
!!
!      itest_N=1  
!      if (itest_N.eq.1) then
!
!      j=jmax
!      do i=1,imax
!      do k=1,kmax
!        i3=k+(i-1)*kmax+(j-1)*imax*kmax
!!
!! ---- Dirichlet BC if inflow (u2 < 0 along N boundary)
!!
!      if (u2(i3).lt.0.0) then
!          if (k.eq.1.and.i.eq.1) print*, 'wind from N'
!        iu=iu+1
!        iv=iv+1
!        iw=iw+1
!        istr=istr+1
!        icu(iu)=i3
!        icv(iv)=i3
!        icw(iw)=i3
!        ice(istr)=i3
!        fixu(iu)=u1(i3)   
!        fixv(iv)=u2(i3)   
!        fixw(iw)=u3(i3)   
!        fixe(istr)=pt(i3)  
!          if (iturb.eq.1) then
!            ik=ik+1
!            ick(ik)=i3
!            fixk(ik)=tk(i3)
!            fixd(ik)=td(i3)
!          endif
!      endif
!      enddo
!      enddo
!
!      endif
!!-------------------------------------------------------------------
!!
!!  EASTERN BOUNDARY (i=imax):  conditional
!!
!      itest_E=1  ! 
!      if (itest_E.eq.1) then
!
!      i=imax
!      do k=1,kmax 
!      do j=1,jmax
!       i3=k+kmax*(i-1)+kmax*imax*(j-1)
!!
!! ---- Dirichlet BC if inflow (u1 < 0 along E boundary)
!!
!      if (u1(i3).lt.0.0) then
!         if (k.eq.1.and.j.eq.1) print*, 'wind from E'
!      
!        iu=iu+1
!        iv=iv+1
!        iw=iw+1
!        istr=istr+1
!        icu(iu)=i3
!        icv(iv)=i3
!        icw(iw)=i3
!        ice(istr)=i3
!        fixu(iu)=u1(i3)
!        fixv(iv)=u2(i3)
!        fixw(iw)=u3(i3)
!        fixe(istr)=pt(i3)    
!          if (iturb.eq.1) then
!            ik=ik+1
!            ick(ik)=i3
!            fixk(ik)=tk(i3)
!            fixd(ik)=td(i3)
!          endif
!       endif
!      enddo
!      enddo
!      
!      endif
!!-------------------------------------------------------------------
!!
!!  WESTERN BOUNDARY (i=1):   Conditional     
!!
!      itest_W=1 
!      if (itest_W.eq.1) then
!
!      do k=1,kmax
!      do j=1,jmax
!       i3=k+kmax*imax*(j-1)
!!
!! ---- Dirichlet BC if inflow (u1 > 0 along W boundary)
!!
!      if (u1(i3).gt.0.0) then
!          if (k.eq.1.and.j.eq.1) print*, 'wind from W'
!        iu=iu+1
!        iv=iv+1
!        iw=iw+1
!        istr=istr+1
!        icu(iu)=i3
!        icv(iv)=i3
!        icw(iw)=i3
!        ice(istr)=i3
!        fixu(iu)=u1(i3)
!        fixv(iv)=u2(i3)
!        fixw(iw)=u3(i3)
!        fixe(istr)=pt(i3)  
!          if (iturb.eq.1) then
!            ik=ik+1
!            ick(ik)=i3
!            fixk(ik)=tk(i3)
!            fixd(ik)=td(i3)
!          endif
!       endif
!      enddo
!      enddo
!
!      endif
!-------------------------------------------------------------------
!
!  UPPER BOUNDARY (K=KMAX):
!
      do i=1,imax
      do j=1,jmax
        i3=kmax+(i-1)*kmax+(j-1)*imax*kmax
!
! ---- Dirichlet BC for (K, eps)
!
        ik=ik+1
        ick(ik)=i3
        fixk(ik)=tk(i3)
        fixd(ik)=td(i3)
!
! ---- Dirichlet BC for W 
!
        iw=iw+1
        icw(iw)=i3
        fixw(iw)=u3(i3)    

! ---- Dirichlet BC for temp
        istr=istr+1
        ice(istr)=i3
        fixe(istr)=pt(i3)

        ! extra:
!        iu=iu+1
!        icu(iu)=i3
!        fixu(iu)=u1(i3)
!        iv=iv+1
!        icv(iv)=i3
!        fixv(iv)=u2(i3)
      enddo
      enddo
!-------------------------------------------------------------------
!
! Secure safe bottom boundaries
!
      do i=1,imax
        ipN=1+(i-1)*kmax+(jmax-1)*imax*kmax
        iu=iu+1; iv=iv+1; iw=iw+1
        icu(iu)=ipN; icv(iv)=ipN; icw(iw)=ipN
        fixu(iu)=0.0; fixv(iv)=0.0; fixw(iw)=0.0

        ipS=1+(i-1)*kmax
        iu=iu+1; iv=iv+1; iw=iw+1
        icu(iu)=ipS; icv(iv)=ipS; icw(iw)=ipS
        fixu(iu)=0.0; fixv(iv)=0.0; fixw(iw)=0.0
      enddo

      do j=1,jmax
        ipE=1+(imax-1)*kmax+(j-1)*imax*kmax
        iu=iu+1; iv=iv+1; iw=iw+1
        icu(iu)=ipE; icv(iv)=ipE; icw(iw)=ipE
        fixu(iu)=0.0; fixv(iv)=0.0; fixw(iw)=0.0

        ipW=1+(j-1)*imax*kmax
        iu=iu+1; iv=iv+1; iw=iw+1
        icu(iu)=ipW; icv(iv)=ipW; icw(iw)=ipW
        fixu(iu)=0.0; fixv(iv)=0.0; fixw(iw)=0.0
      enddo
!-------------------------------------------------------------------
!
! ----- Estimate other parameters
!
      do ie=1,nelem
        pd(ie)=0.0
      enddo

      do i=1,mb
        icp(i)=0
      enddo
!-------------------------------------------------------------------
!
!  ----- Summation: number of BC points
!
      nfixp=0
      nfixu=iu
      nfixv=iv
      nfixw=iw
      nfixk=ik
      nfixe=istr
!-------------------------------------------------------------------
!
! ----- Calculate distance from ground to field points
!       (for low-Rk turbulence model)
!
      do i=1,imax
        do j=1,jmax
          i1=1+kmax*(i-1)+imax*kmax*(j-1)   ! wall point
          do k=1,kmax
            ip=k+kmax*(i-1)+imax*kmax*(j-1)
            r2=(coord(1,ip)-coord(1,i1))**2 &
     &      +(coord(2,ip)-coord(2,i1))**2 &
     &      +(coord(3,ip)-coord(3,i1))**2
            dist(ip)=sqrt(r2)+z0
          enddo
        enddo
      enddo
!********************************************************************
!
! ----- INITIAL CONDITIONS
!
!********************************************************************
!
! ----- Output of initial BC's to 'init.dat'
!       ... data written in physical dimensions
!
      open (unit=14,file='init.dat',form='unformatted')
      time=0.0
      write(14) time,(u1(i),u2(i),u3(i),ps(i),tk(i),td(i), &
     &          vtef(i),pt(i),pts(i),rho(i),rhos(i),i=1,npoin)
      write(14) (pd(i),i=1,nelem)
      close(14)
      print*,'  init.dat completed'
!-------------------------------------------------------------------
!
! ----- Output of BC's to 'boun.dat', data in SI units
!
      open (unit=16,file='boun.dat',status='unknown')
      write(16,900) 'Boundary conditions'
      write(16,910)  nfixu,nfixv,nfixw,nfixp,nfixe,nfixk,new,nlog,z0
      write(16,940) (z0l(i),i=1,nlog)
      write(16,900) 'U-component :'
      write(16,920) (icu(i),fixu(i),i=1,nfixu)
      write(16,900) 'V-component :'
      write(16,920) (icv(i),fixv(i),i=1,nfixv)
      write(16,900) 'W-component :'
      write(16,920) (icw(i),fixw(i),i=1,nfixw)
      write(16,900) 'P-points :'
      write(16,930) (icp(i),i=1,nfixp)
      write(16,900) 'Wall/Log-points:'
      write(16,930) (kon_we(i),i=1,new)
      write(16,930) (icod_b(i),icod_n(i),i=1,nlog)
      write(16,900) '(K,e)-values:'
      write(16,930) (ick(i),i=1,nfixk)
      write(16,940) (fixk(i),fixd(i),i=1,nfixk)
      write(16,940) (dist(ip),ip=1,npoin)
      write(16,900) 'Pot.temp :'
      write(16,930) (ice(i),i=1,nfixe)
      write(16,940) (fixe(i),i=1,nfixe)
      print*,'  boun.dat completed'
!------------------------------------------------------------------
      print*,' ...ibc_UM2Simra finished'
!------------------------------------------------------------------
 900  format(a)
 910  format(8i7,e12.4)
 920  format(3(i8,e12.4))
 930  format(9i8)
 940  format(6e12.4)
!------------------------------------------------------------------
      END

!###################################################################
      subroutine checkFlowDir (isInflow,u,v,w,normal)                                                                   
      !    Check if flow is inflow
      logical isInflow
      real u,v,w,normal(3)
      
      isInflow = (u*normal(1)+v*normal(2)+w*normal(3)).lt.0.0

      END subroutine checkFlowDir
!*******************************************************************
!
      subroutine hstat (pts,ps,rhos,coord,imax,jmax,kmax,mp,mb)
!===================================================================
!  Purpose:
!     Computation of hydrostatic pressure & density (ps, rhos).
!
!  Method:
!     Hydrostatic pressure (ps) and density (rhos) are computed
!     from the hydrostatic equation and the state equation,
!     given hydrostatic potential temperature
!     at a horizontal reference position (ihs, jhs).
!     These values are then interpolated to other positions.
!
!  Input:
!     pts  = hydrostatic potential temperature at reference point
!     constants in state equation (se below)
!
!  Output:
!     (ps, rhos) = computed hydrostatic pressure and density
!------------------------------------------------------------------
      integer i,j,ie,ip,npoin,nelem,kmax,k,mp,mb
      real coord(3,mp),ps(mp),pts(mp),rhos(mp),p_ref,R_gas,gamma
      real t_s(mb),z_ref(mb),xa(mb),yap(mb),yar(mb)
!
! ----- Constants
!
      grav=9.81
      p_ref = 100000.   ! static pressure at sea level [Pa]
      R_gas = 287.
      gamma = 2./7.
!-----------------------------------------------------------------
!
! ----- For given horizontal point (i,j), find (rhos,ps) along z-axis
!       assuming given pts-profile
!
      z_ref(1)=0.0 ! at sea level
      i=1; j=1
      do k=2,kmax
        ip=k+kmax*(i-1)+imax*kmax*(j-1)
        z_ref(k)=coord(3,ip)       ! measured from sea level
      enddo

! ..... Niveau 1:  sea level, k = 1
      k=1
      ip1=1+kmax*(i-1)+imax*kmax*(j-1)
      ps(ip1)=p_ref
      Ts1=pts(ip1)
      rhos(ip1)=p_ref/(R_gas*pts(ip1))                ! pts given
      xa(1)=z_ref(1)
      yap(1)=ps(ip1)
      yar(1)=rhos(ip1)

! .... Other levels:
      do it=1,2

        do k=2,kmax
          ip=k+kmax*(i-1)+imax*kmax*(j-1)
          ipm=k-1+kmax*(i-1)+imax*kmax*(j-1)
          dz=z_ref(k)-z_ref(k-1)
          if (it.eq.1) then
            rhos_avr=rhos(ipm)
          else
            rhos_avr=0.5*(rhos(ipm)+rhos(ip))
          endif
!      ... Hydrostatic Equation
          ps(ip)=ps(ipm)-grav*rhos_avr*dz            ! hydrostat. eq.
!      ... State Equation
          rhos(ip)=ps(ip)/(R_gas*pts(ip))*(p_ref/ps(ip))**gamma

! ..... Storing for interpolation:
          yap(k)=ps(ip)
          yar(k)=rhos(ip)
          xa(k)=z_ref(k)
        enddo

      enddo
!-----------------------------------------------------------------
!
! ----- Printing of density profile at horizontal point
!
      print *, ' Hydrostatic density profile at point (1,1):'
      print *,'   Z coord [m] above sea level:'
      write(*,900) (coord(3,k),k=1,kmax)
      print *,'   rho_s [kg/m^3] :'
      write(*,910) (rhos(k),k=1,kmax)
      print *,'     p_s [N/m^2] :'
      write(*,920) (ps(k),k=1,kmax)
 900  format(8f9.2)
 910  format(8f9.3)
 920  format(8e10.2)
!----------------------------------------------------------------
!
! ----- INTERPOLATE the hydrostatic fields (ps, rhos)
!         from reference positions to actual domain points
!
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
            ip=k+kmax*(i-1)+imax*kmax*(j-1)
            z=coord(3,ip)       ! measured from sea level
            call polint (xa,yar,kmax,z,y,mb)
            rhos(ip)=y          ! same as in non-dim.form
            call polint (xa,yap,kmax,z,y,mb)  ! interpolation
            ps(ip)=y            ! dimensional
          enddo
        enddo
      enddo
!----------------------------------------------------------------
      END subroutine hstat

!******************************************************************
      subroutine polint (xa,ya,n,x,y,mb)
!==================================================================
!
!  Purpose:
!     Linear interpolation
!  Method:
!      Given:     ya(i) = f{xa(i)}, i = 1, n
!      Find :     y     = f(x) for given x inside xa(i)
!  Input:
!      1-D coordinate xa(i)
!      Corresponding function value ya(i)
!      Coordinate value x
!  Output:
!      Interpolated function value y
!      corresponding to coordinate x
!  Restrictions:
!      If x < min (xa(i)) then  y = ya (xa_min)
!      If x > max (xa(i)) then  y = ya (xa_max)
!
!  22/5/01    T. Utnes
!------------------------------------------------------------------
      implicit none
      integer i,n,mb
      real x,y,den,xa(mb),ya(mb)
!------------------------------------------------------------------

      if (x.le.xa(1)) then    ! xa(1) = xa_min
        y=ya(1)               ! lower level restriction
        RETURN
      endif

      do i=1,n-1
        if (x.ge.xa(i).and.x.le.xa(i+1)) then
          den=xa(i+1)-xa(i)
          y=ya(i)+(ya(i+1)-ya(i))*(x-xa(i))/den
          RETURN
        endif
      enddo
      if (x.ge.xa(n)) y=ya(n)   ! upper level restriction

      END subroutine polint


!*******************************************************************
      subroutine mean (imax,jmax,kmax,x,xx,k,mp)
!=================================================================== 
!    Purpose:                                                      
!            to average a variable x() at given level k
!-------------------------------------------------------------------
      implicit none
      integer imax,jmax,kmax,k,i,j,ij,imjm,ipjm,ipjp,imjp,ijm,ipj, &
     &        ijp,imj,mp
      real x(mp),xx(mp)
!-------------------------------------------------------------------
!
! ---- Interpolation at level k
!
      do j=2,jmax-2
        do i=2,imax-2
          ij=k+kmax*(i-1)+imax*kmax*(j-1)
          imjm=k+kmax*(i-2)+imax*kmax*(j-2)
          ipjm=k+kmax*(i)+imax*kmax*(j-2)
          ipjp=k+kmax*(i)+imax*kmax*(j)
          imjp=k+kmax*(i-2)+imax*kmax*(j)
          ijm=k+kmax*(i-1)+imax*kmax*(j-2)
          ipj=k+kmax*(i)+imax*kmax*(j-1)
          ijp=k+kmax*(i-1)+imax*kmax*(j)
          imj=k+kmax*(i-2)+imax*kmax*(j-1)
          xx(ij)=0.0625*(x(imjm)+x(ipjm)+x(ipjp)+x(imjp)) &
     &       +0.125*(x(ijm)+x(ipj)+x(ijp)+x(imj))+0.25*x(ij)
        enddo
      enddo
      RETURN
      END

!***********************************************************************
      subroutine profile(coord,imax,jmax,kmax,delta,ustar,z0,u1,u2,u3, &
     &                   pts,istrat,mp,m1)
!=======================================================================
      real coord(3,mp),u1(mp),u2(mp),u3(mp),pts(mp)
      real pro_u(m1), pro_v(m1), pro_w(m1), pro_t(m1),delta,ustar
!-----------------------------------------------------------------------

         grav=9.81
         delta=delta+z0
         u_a=SED_U_A  ! wind at delta   

         dz=abs(coord(3,2)-coord(3,1))
         cappa=0.42
!
! ---- Geostrophic values (upper boundary)    
!
         eta=1.0
         wake=3*eta-2*eta**2
         !wake=0.
         ustar=u_a*cappa/(alog(delta/z0)+wake)
!
! ---- Wind direction relative to x-axis = alpha[deg]
!
         alpha=180. !   E wind
         alpha=90. !   S wind
         alpha=67.5 !   SSW wind (202.5 grader)
         alpha=315. !  NW 

         alpha=135. !  SE wind
         alpha=90. !   S wind

         alpha=alpha/57.295
         
         phi=SED_PHI
         alpha=270.0-phi
         alpha=alpha/57.295
!
! ---- Stratification:
!        T0 = temp. at ground, BVN = Brunt-Vaisala freq.
!
         T0 = 290.
         BVN=1.0   ! = N.h/U, dvs F = 1.0
         BVN=0.0   ! dvs. neutral
!-------------------------------------------------------------------
!
!  ----- Inflow profile
!
       DO k=1,kmax
          ib=1
          z_bakke=coord(3,ib)
          z=(coord(3,k)-z_bakke)+z0
          eta=min(z/delta,1.0)
          wake=3*eta-2*eta**2
          !wake=0.

!   ----- velocity
          profu=min((ustar/cappa)*(alog(z/z0)+wake),u_a) 
          pro_u(k)=profu*cos(alpha)
          pro_v(k)=profu*sin(alpha)
          pro_w(k)=0.0

!   ----- potensial temp.
          z1=z-z0
          !!pro_t(k)=T0*exp((BVN**2)*z1/grav)
          if (istrat.eq.1) then
            call teta_profile(z1,t1) ! see subroutine
            pro_t(k)=t1
          elseif (istrat.eq.0) then
            pro_t(k)=T0  ! constant
          endif
       ENDDO
!--------------------------------------------------------------------
!
! ----- INITIAL CONDITIONS
!         Inflow profiles used as initialization
!
      do i = 1, imax
        do j = 1, jmax
          do k = 1, kmax
            ip=k+kmax*(i-1)+imax*kmax*(j-1)
            u1(ip)   = pro_u(k)   ! given inflow profile, see above
            u2(ip)   = pro_v(k)
            u3(ip)   = pro_w(k)
            pts(ip)  = pro_t(k)
          enddo
        enddo
      enddo
!-------------------------------------------------------------------
      END subroutine profile


!*******************************************************************
      SUBROUTINE eddy (kmax,delta,ustar,z,rho,u,v,tk,td,av,z0,m1)
!=================================================================== 
!                                                                   
!    1D (vertical profile) eddy viscosity model 
!
!    - given vertical geometry, velocity and density
!    - output: 
!       eddy viscosity (av), turb kin energy (tk), dissipation (td)
!                                                                   
!    Programmed: 
!       July 2004: T. Utnes                                          
!    Revised:
!       Oct. 2004: TU
!-------------------------------------------------------------------
      integer kmax,m1
      real tk(m1),td(m1),u(m1),v(m1),av(m1),z(m1),rho(m1),delta,ustar
      real z0
      av0=1.0e-5; rho0=1.3
      cm34=0.09**0.75
      u_top=sqrt(u(kmax)**2+v(kmax)**2)
      ub=sqrt(u(2)**2+v(2)**2)
      tkmin=0.001*u_top*u_top
      tkmin=0.0005*u_top*u_top
      !ustar=0.4*u_top/(alog(delta/z0)) ! from 'profile' via parameter

      do k=2,kmax   
        kp=min(k+1,kmax)    
        km=max(k-1,1)

        ! Neutral length scale
        z1=z(k)-z(1)+z0   
        al0=0.4*z1/(1.0+0.4*z1/40.0) ! cf Lock, p.10, Note, eq(8)

        ! Ri-number
        dz=z(kp)-z(km)
        drdz=(rho(kp)-rho(km))/dz
        any2=-(9.81/rho0)*drdz    
        dudz=(u(kp)-u(km))/dz
        dvdz=(v(kp)-v(km))/dz
        prod=(dudz)**2+(dvdz)**2 + 1.0e-8
        Ri=any2/prod

        ! Eddy viscosity & turbulence parameters
        fact=max(1.0-0.5*Ri, 0.00001)
        av(k)=al0*al0*sqrt(prod)*sqrt(fact)        ! Note, eq(7)
        av(k)=max(av(k),av0)

        tk(k)=av(k)**2/(0.3*al0*al0)               ! Note, eq(9)

        !.....simple alternative........
        rm=max(1.0-z1/delta,0.0)
        tk(k)=3.3333*ustar*ustar*rm*rm
        ! tk(k)=1.2*tk(k)
        ! Alternative >>
         al0=0.4*z1/(1.0+4*z1/delta)
        !...............................

        tk(k)=max(tk(k),tkmin)

        td(k)=cm34*tk(k)**1.5/al0                ! standard
      enddo

      tk(2)=(0.4*ub)**2/(0.3*alog(z(2)/z0)**2)  ! near-bottom BC

      tk(1)=tk(2); td(1)=td(2)  ! not used at bottom level

      tk = max(tk,tkmin)
!-----------------------------------------------------------------
      END SUBROUTINE eddy 

!******************************************************************
      subroutine teta_profile(z,t)
!
!     Potential temperature profile,
!               similar to case Hammerfest 1/5/05
!------------------------------------------------------------------

      t0=272.3; z1=800.

        a1=0.0007
        if (z.le.z1) t=t0+a1*z
        t1=t0+a1*z1
        dz12=1100.; a12=0.009
        if (z.gt.z1.and.z.le.z1+dz12) t=t1+a12*(z-z1)
        z2=z1+dz12; t2=t1+a12*dz12 
        dz23=700.
        z3=z2+dz23; a23=a12
        if (z.gt.z2.and.z.le.z3) t=t2+a23*(z-z2)
        t3=t2+a23*dz23
        a3=a1
        if (z.gt.z3) t=t3+a3*(z-z3)

      end subroutine teta_profile

