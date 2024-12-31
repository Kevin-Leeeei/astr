!+---------------------------------------------------------------------+
!| This module contains user defined subroutines to interfere  program |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 18-08-2023  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module userdefine
  !
  implicit none
  !
  real(8) :: hsource
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to set flow environment, such as, incoming     |
  !| free stream variables.                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_setflowenv
    !
!     use commvar,  only: roinf,uinf,vinf,winf,pinf,tinf,spcinf,num_species
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     real(8) :: specr(num_species)
!     ! 
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     ! pinf=5.d0*pinf
!     uinf=0.97d0
!     vinf=0.d0
!     winf=0.d0
!     tinf=300.d0
!     spcinf(:)=specr(:)
!     roinf=thermal(pressure=pinf,temperature=tinf,species=spcinf(:))
!     !
! #endif
    !
  end subroutine udf_setflowenv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_setflowenv.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate fluctuations for inflow            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Oct-2023: Created by by Jian Fang @ Daresbury                  |
  !+-------------------------------------------------------------------+
  subroutine udf_inflow_fluc(umean,uinst)
    !
    use commvar, only : jm,km
    !
    real(8),intent(in) ::  umean(0:jm,1:3)  ! inflow mean velocity
    real(8),intent(out) :: uinst(0:jm,0:km,1:3)  ! velocity with fluctuations
    !
  end subroutine udf_inflow_fluc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_inflow_fluc.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise flowfield by a user              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-May-2023: Created by Yifan Xu @ Peking University              |
  !| 18-Aug-2023: Rename and relocated by Jian Fang @ Daresbury        |
  !+-------------------------------------------------------------------+
  subroutine udf_flowinit
    !
!     use commvar,  only: im,jm,km,ndims,roinf,uinf,nondimen,xmax,pinf,  &
!                         ia,num_species
!     use commarray,only: x,vel,rho,prs,spc,tmp,q
!     use parallel, only: lio
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     !
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     ! local data
!     integer :: i,j,k
!     real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specr(num_species),  &
!       specp(num_species),arg,prgvar,masflx,specx(num_species)
!     real(8) :: pthick
!     !
!     tmpr=300.d0
!     xloc=3.d0*xmax/4.d0
!     xwid=xmax/(12.d0*5.3d0*2.d0)
!     !
!     !reactants
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     !products
!     tmpp=1814.32d0
!     !
!     ! pthick=1.d-4
!     !
!     do k=0,km
!     do j=0,jm
!     do i=0,im
!       !
!       xc=x(i,j,k,1)
!       !
!       !prgvar=0.5d0*(1.d0+tanh(10.d0*(xc-xloc)/xloc))
!       ! if(xc-xloc<xwid*0.5d0*1.2d0) then 
!       !   prgvar=0.d0
!       !   if(xc-xloc>xwid*0.5d0) &
!       !   prgvar=1.d0-(xc-xloc-(xwid*0.5d0))/(xwid*0.5d0*0.2d0)
!       ! else
!       !   prgvar=1.d0
!       ! endif
!       !
!       prgvar=1.d0*exp(-0.5d0*((xc-xloc)/xwid)**2)
!       !
!       spc(i,j,k,:)=specr(:)
!       !
!       vel(i,j,k,1)=uinf
!       !
!       vel(i,j,k,2)=0.d0
!       vel(i,j,k,3)=0.d0
!       !
!       tmp(i,j,k)=tmpr+prgvar*(tmpp-tmpr)
!       !
!       prs(i,j,k)=pinf
!       !
!       rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
!                           species=spc(i,j,k,:))
!     enddo
!     enddo
!     enddo
!     !
!     !
!     if(lio)  write(*,'(A,I1,A)')'  ** HIT flame initialised.'
!     !
! #endif
    !
  end subroutine udf_flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_flowinit.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate grid.                              | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_grid
  end subroutine udf_grid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_grid.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to list something during a computation.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !| 12-Oct-2023: adapted by Chensheng Luo @ Beihang University        |
  !+-------------------------------------------------------------------+
  subroutine udf_stalist
    !
    use constdef
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km,roinf
    use commarray,only : vel,rho,tmp,dvel,q,vorbis,dvor,prs
    use fludyna,  only : miucal,sos
    use comsolver,only : solvrinit,grad
    use utility,  only : listinit,listwrite
    use parallel, only : dataswap,psum,lio,pmax,pmin
    use constdef, only : pi
    !
    integer :: i,j,k,ns
    real(8) :: s11,s12,s13,s23,s22,s33,div,miu,dissa,dissloc,kolmloc
    real(8) :: du11,du12,du21,du22,du33,du23,du13,du32,du31
    real(8) :: m11m11,m22m22,m11m22,m12m21,m12m12,m21m21,m33m33
    real(8) :: m13m13,m23m23,m31m31,m32m32
    real(8) :: m11m11m11,m22m22m22,m11m11m22,m22m22m11,m11m12m12,&
              m22m21m21,m22m21m12,m11m12m21,m11m21m21,m22m12m12,&
              m33m33m33
    real(8) :: dPhi, dPsi,comlen
    real(8) :: divU,Omega,Psi,Phi
    real(8) :: divUdivU,dPsidPsi,dPhidPhi,OmegadivU,OmegaOmega
    real(8) :: dPsidPsidivU,dPhidPhidivU,divUdivUdivU,OmegaOmegadivU,traceSS,traceSSth,wsw,s3
    real(8) :: rsamples,miudrho,dudx2,csavg,v2,cs,ufluc,ens,omegaz,omegax,omegay
    real(8) :: rhoavg,rho2nd,w2drho
    real(8) :: urms,energy,taylorlength,kolmoglength,Retaylor,machrms,macht,rhoe,skewness,du2,du3,ReL
    real(8) :: mfpath,R ! mfpath = mean free path
    real(8) :: rhoprime,Kw,Krho,ptheta,uutheta,rhouugrad
    !
    logical :: fex
    logical,save :: linit=.true.
    integer,save :: hand_fs,hand_mom2nd,hand_mom3rd,hand_skew,hand_en,hand_weakly
    !
    if(ka==0) then 
      ! 2D part
      R = 8.31446261815324d0
      !
      if(linit) then
        !
        if(lio) then
          call listinit(filename='log/fturbstats2d.dat',handle=hand_fs, &
             firstline='nstep time urms ens talen 2Dkolm kolmloc Rel ReL Ensdis mfpath machrms macht Tavg skewness')
          call listinit(filename='log/mom2nd.dat',handle=hand_mom2nd, &
                        firstline='nstep time m11m11 m22m22 m11m22 m12m21 m12m12 m21m21')
          call listinit(filename='log/mom3rd.dat',handle=hand_mom3rd, &
                        firstline='nstep time A111111 A222222 A111122 A222211 A111212 A221212 A112121 A222121 A222112 A111221')
          call listinit(filename='log/skewness.dat',handle=hand_skew, &
                        firstline='ns time th o ps f th2 ps2 f2 oth o2 ps2th f2th th3 o2th m11s m22s m11c m22c')
          call listinit(filename='log/fenergy.dat',handle=hand_en, &
                        firstline='ns time rhoavg rho2nd w2drho')
          call listinit(filename='log/weakly.dat',handle=hand_weakly, &
                        firstline='ns time Kw Krho ptheta uutheta rhouugrad')
        endif
        !
        linit=.false.
        !
      endif
      !
      rsamples=dble(ia*ja)
      !
      urms=0.d0
      energy=0.d0
      !
      divU=0.d0
      comlen=2*pi
      kolmloc=2*pi
      mfpath=0.d0
      Omega=0.d0
      Psi=0.d0
      Phi=0.d0
      !
      OmegaOmega=0.d0
      OmegadivU=0.d0
      divUdivU=0.d0
      dPsidPsi=0.d0
      dPhidPhi=0.d0
      !
      OmegaOmegadivU=0.d0
      divUdivUdivU=0.d0
      dPsidPsidivU=0.d0
      dPhidPhidivU=0.d0
      !
      m11m11=0.d0
      m22m22=0.d0
      m11m22=0.d0
      m12m21=0.d0
      m12m12=0.d0
      m21m21=0.d0
      m11m11m11=0.d0
      m22m22m22=0.d0
      m11m11m22=0.d0
      m22m22m11=0.d0
      m11m12m12=0.d0
      m22m12m12=0.d0
      m11m21m21=0.d0
      m22m21m21=0.d0
      m11m12m21=0.d0
      m22m21m12=0.d0
      !
      machrms=0.d0
      !
      miudrho=0.d0
      dudx2=0.d0
      csavg=0.d0
      dissa=0.d0
      ens=0.d0
      du3=0.d0
      du2=0.d0
      rhoavg=0.d0
      rho2nd=0.d0
      w2drho=0.d0
      Kw=0.d0
      Krho=0.d0
      ptheta=0.d0
      uutheta=0.d0
      rhouugrad=0.d0
      !
      k=0
      do j=0,jm
      do i=0,im
        vorbis(i,j,k)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      enddo
      enddo
      !
      call dataswap(vorbis)
      !
      !
      dvor=grad(vorbis)
      !
      k=0
      do j=1,jm
      do i=1,im
        !
        ! This point values
        !
        miu=miucal(tmp(i,j,k))/Reynolds
        rhoprime = rho(i,j,k) - roinf
        !
        du11=dvel(i,j,k,1,1)
        du12=dvel(i,j,k,1,2)
        du21=dvel(i,j,k,2,1)
        du22=dvel(i,j,k,2,2)
        !
        s11=du11
        s12=0.5d0*(du12+du21)
        s22=du22
        !
        div=du11+du22
        !
        omegaz=du21-du12
        !
        dPsi = du12+du21
        dPhi = du11-du22
        !
        v2=vel(i,j,k,1)**2+vel(i,j,k,2)**2
        !
        cs=sos(tmp(i,j,k))
        !
        ! Volume average values 
        !
        urms = urms+v2
        energy = energy + 0.5d0*rho(i,j,k)*v2
        ens = ens + 0.5d0*(omegaz*omegaz)
        !
        divU = divU+div
        comlen = min(comlen,sqrt(miu/rho(i,j,k)/abs(div)))
        Omega = Omega + omegaz
        Psi = Psi + dPsi
        Phi = Phi + dPhi
        !
        OmegaOmega = OmegaOmega + omegaz*omegaz
        OmegadivU = OmegadivU + omegaz*div
        divUdivU = divUdivU + div*div
        dPsidPsi = dPsidPsi + dPsi*dPsi
        dPhidPhi = dPhidPhi + dPhi*dPhi
        !
        OmegaOmegadivU = OmegaOmegadivU + omegaz*omegaz*div
        divUdivUdivU = divUdivUdivU + div*div*div
        dPsidPsidivU = dPsidPsidivU + dPsi*dPsi*div
        dPhidPhidivU = dPhidPhidivU + dPhi*dPhi*div
        !
        m11m11=m11m11+du11*du11
        m22m22=m22m22+du22*du22
        m11m22=m11m22+du11*du22
        m12m21=m12m21+du12*du21
        m12m12=m12m12+du12*du12
        m21m21=m21m21+du21*du21
        !
        m11m11m11=m11m11m11+du11*du11*du11
        m22m22m22=m22m22m22+du22*du22*du22
        m11m11m22=m11m11m22+du11*du11*du22
        m22m22m11=m22m22m11+du22*du22*du11
        m11m12m12=m11m12m12+du11*du12*du12
        m22m12m12=m22m12m12+du22*du12*du12
        m11m21m21=m11m21m21+du11*du21*du21
        m22m21m21=m22m21m21+du22*du21*du21
        m11m12m21=m11m12m21+du11*du12*du21
        m22m21m12=m22m21m12+du22*du21*du12
        !
        dudx2=dudx2+du11**2+du22**2
        !
        miudrho=miudrho+miu/rho(i,j,k)
        !
        csavg=csavg+cs
        !
        machrms=machrms+v2/(cs*cs)
        !
        rhoe=rhoe+tmp(i,j,k)
        !
        ! This is the enstrophy dissipation!
        dissa = dissa+miu/rho(i,j,k)*(dvor(i,j,k,1)**2+dvor(i,j,k,2)**2)
        ! This is the energy dissipation!
        dissloc = 2.d0*miu/rho(i,j,k)*(s11**2+s22**2+2.d0*(s12**2)-0.5d0*div**2)
        kolmloc = min(kolmloc,sqrt(sqrt((miu/rho(i,j,k))**3/dissloc)))
        !
        !
        du3=du3+(du11*du11*du11+du22*du22*du22)
        du2=du2+(du11*du11+du22*du22)
        !
        rhoavg = rhoavg + rho(i,j,k)
        rho2nd = rho2nd + rho(i,j,k)*rho(i,j,k)
        w2drho = w2drho + (omegaz*omegaz)/rho(i,j,k)
        !
        mfpath = max(mfpath,2.d0*miu/rho(i,j,k)/0.921/sqrt(3.d0*R*tmp(i,j,k)))
        !
        Kw       = Kw+ 0.5d0*roinf*(vel(i,j,k,1)**2+vel(i,j,k,2)**2)
        Krho     = Krho + 0.5d0*rhoprime*(vel(i,j,k,1)**2+vel(i,j,k,2)**2)
        ptheta   = ptheta + prs(i,j,k)*div
        uutheta  = uutheta + 0.5d0*roinf*(vel(i,j,k,1)**2+vel(i,j,k,2)**2)*div
        rhouugrad= rhouugrad + rhoprime * (vel(i,j,k,1)*vel(i,j,k,1)*du11 + &
                                           vel(i,j,k,1)*vel(i,j,k,2)*du12 + &
                                           vel(i,j,k,2)*vel(i,j,k,1)*du21 + &
                                           vel(i,j,k,2)*vel(i,j,k,2)*du22)
        !
      enddo
      enddo
      !
      urms= sqrt(psum(urms)/rsamples)
      energy= psum(energy)/rsamples
      ens=psum(ens)/rsamples
      !
      divU  = psum(divU)/rsamples
      Omega  = psum(Omega)/rsamples
      Psi  = psum(Psi)/rsamples
      Phi  = psum(Phi)/rsamples
      !
      OmegaOmega  = psum(OmegaOmega)/rsamples
      OmegadivU  = psum(OmegadivU)/rsamples
      divUdivU  = psum(divUdivU)/rsamples
      dPsidPsi  = psum(dPsidPsi)/rsamples
      dPhidPhi  = psum(dPhidPhi)/rsamples
      !
      OmegaOmegadivU  = psum(OmegaOmegadivU)/rsamples
      divUdivUdivU  = psum(divUdivUdivU)/rsamples
      dPsidPsidivU  = psum(dPsidPsidivU)/rsamples
      dPhidPhidivU  = psum(dPhidPhidivU)/rsamples
      !
      m11m11  = psum(m11m11)/rsamples
      m22m22  = psum(m22m22)/rsamples
      m11m22  = psum(m11m22)/rsamples
      m12m21  = psum(m12m21)/rsamples
      m12m12  = psum(m12m12)/rsamples
      m21m21  = psum(m21m21)/rsamples
      !
      m11m11m11  = psum(m11m11m11)/rsamples
      m22m22m22  = psum(m22m22m22)/rsamples
      m11m11m22  = psum(m11m11m22)/rsamples
      m22m22m11  = psum(m22m22m11)/rsamples
      m11m12m12  = psum(m11m12m12)/rsamples
      m22m12m12  = psum(m22m12m12)/rsamples
      m11m21m21  = psum(m11m21m21)/rsamples
      m22m21m21  = psum(m22m21m21)/rsamples
      m11m12m21  = psum(m11m12m21)/rsamples
      m22m21m12  = psum(m22m21m12)/rsamples
      !
      dudx2      = num1d3*psum(dudx2)/rsamples
      miudrho    = psum(miudrho)/rsamples
      csavg      = psum(csavg)/rsamples
      dissa      = psum(dissa)/rsamples
      !
      machrms=sqrt(psum(machrms)/rsamples)
      !
      rhoe=psum(rhoe)/rsamples
      !
      ufluc=urms/sqrt(2.d0)
      !
      macht         = urms/csavg
      taylorlength  = sqrt(miudrho*ens/dissa) ! enstrophy based
      retaylor      = ens**(3.d0/2.d0)/dissa ! enstrophy based
      ReL           = energy/(miudrho*(dissa**(1.d0/3.d0))) ! enstrophy based
      kolmoglength  = (miudrho**3/dissa)**(1.d0/6.d0)
      kolmloc       = pmin(kolmloc)
      comlen        = pmin(comlen) ! Compressible length
      mfpath        = pmax(mfpath)
      !
      skewness      = psum(du3)/(2.d0*rsamples)/sqrt((psum(du2)/(2.d0*rsamples))**3)
      ! kolmogvelocity= sqrt(sqrt(dissipation*miudrho))
      ! kolmogtime    = sqrt(miudrho/dissipation)
      rhoavg = psum(rhoavg)/rsamples
      rho2nd = psum(rho2nd)/rsamples
      w2drho = psum(w2drho)/rsamples
      !
      Kw       = psum(Kw)/rsamples
      Krho     = psum(Krho)/rsamples
      ptheta   = psum(ptheta)/rsamples
      uutheta  = psum(uutheta)/rsamples
      rhouugrad= psum(rhouugrad)/rsamples
      !
      if(lio) then 
        call listwrite(hand_fs,urms,ens,taylorlength,kolmoglength, &
                        kolmloc, Retaylor,ReL,dissa,mfpath,machrms, &
                        macht, rhoe,skewness)
        call listwrite(hand_mom2nd,m11m11,m22m22,m11m22,m12m21,m12m12,m21m21)
        call listwrite(hand_mom3rd,m11m11m11,m22m22m22,m11m11m22,m22m22m11,&
                      m11m12m12,m22m12m12,m11m21m21,m22m21m21,m22m21m12,m11m12m21)
        call listwrite(hand_skew,divU,Omega,Psi,Phi,divUdivU,dPsidPsi,&
                      dPhidPhi,OmegadivU,OmegaOmega,dPsidPsidivU,dPhidPhidivU,&
                      divUdivUdivU,OmegaOmegadivU,m11m11,m22m22,m11m11m11,m22m22m22)
        call listwrite(hand_en,rhoavg,rho2nd,w2drho)
        call listwrite(hand_weakly,Kw,Krho,ptheta,uutheta,rhouugrad)
      endif
    else
      ! 3D part
      R = 8.31446261815324d0
      !
      if(linit) then
        !
        if(lio) then
          call listinit(filename='log/fturbstats.dat',handle=hand_fs, &
              firstline='nstep time urms ens talen kolmavg kolmloc Reta comlen mfpath machrms macht Tavg hsource skewness')
          call listinit(filename='log/skewness.dat',handle=hand_skew, &
                        firstline='ns time m11s m22s m33s m11c m22c m33c th2 o2 s2 th3 o2th s2th wsw s3')
          call listinit(filename='log/mom2nd.dat',handle=hand_mom2nd, &
                        firstline='nstep time m12m12 m13m13 m21m21 m23m23 m31m31 m32m32')
          call listinit(filename='log/mom3rd.dat',handle=hand_mom3rd, &
                        firstline='nstep time m11m12m12 m11m12m21 m11m11m22')
        endif
        !
        linit=.false.
        !
      endif
      !
      rsamples=dble(ia*ja*ka)
      !
      urms=0.d0
      energy=0.d0
      !
      comlen=2*pi
      kolmloc=2*pi
      mfpath=0.d0
      !
      OmegaOmega=0.d0
      divUdivU=0.d0
      traceSS=0.d0
      !
      OmegaOmegadivU=0.d0
      divUdivUdivU=0.d0
      traceSSth=0.d0
      wsw = 0.d0
      s3 = 0.d0
      !
      m11m11  = 0.d0
      m22m22  = 0.d0
      m33m33  = 0.d0
      m11m11m11  = 0.d0
      m22m22m22  = 0.d0
      m33m33m33  = 0.d0
      !
      m12m12  = 0.d0
      m13m13  = 0.d0
      m21m21  = 0.d0
      m23m23  = 0.d0
      m31m31  = 0.d0
      m32m32  = 0.d0
      m11m12m12 = 0.d0
      m11m12m21 = 0.d0
      m11m11m22 = 0.d0
      !
      machrms=0.d0
      !
      miudrho=0.d0
      dudx2=0.d0
      csavg=0.d0
      dissa=0.d0
      ens=0.d0
      du3=0.d0
      du2=0.d0
      rhoavg=0.d0
      rho2nd=0.d0
      w2drho=0.d0
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        ! This point values
        !
        miu=miucal(tmp(i,j,k))/Reynolds
        !
        du11=dvel(i,j,k,1,1)
        du12=dvel(i,j,k,1,2)
        du13=dvel(i,j,k,1,3)
        du21=dvel(i,j,k,2,1)
        du22=dvel(i,j,k,2,2)
        du23=dvel(i,j,k,2,3)
        du31=dvel(i,j,k,3,1)
        du32=dvel(i,j,k,3,2)
        du33=dvel(i,j,k,3,3)
        !
        s11=du11
        s12=du12+du21
        s13=du13+du31
        s22=du22
        s23=du23+du32
        s33=du33
        !
        div=du11+du22+du33
        !
        omegax=du32-du23
        omegay=du13-du31
        omegaz=du21-du12
        !
        dPsi = du11-du22
        dPhi = du11-du33
        !
        v2=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
        !
        cs=sos(tmp(i,j,k))
        !
        ! Volume average values 
        !
        urms=urms+v2
        energy = energy + 0.5d0*rho(i,j,k)*v2 !
        !
        comlen = min(comlen,sqrt(4.d0/3.d0*miu/rho(i,j,k)/abs(div)))
        !
        OmegaOmega = OmegaOmega + omegax*omegax  + omegay*omegay + omegaz*omegaz
        divUdivU = divUdivU + div*div
        traceSS = traceSS + s12*s12/2.d0 + s13*s13/2.d0 + s23*s23/2.d0 + 2.d0*dPsi*dPsi/3.d0  &
                - 2.d0*dPsi*dPhi/3.d0 + 2.d0*dPhi*dPhi/3.d0
        !
        OmegaOmegadivU = OmegaOmegadivU + (omegax**2  + omegay**2 + omegaz**2)*div
        divUdivUdivU = divUdivUdivU + div*div*div
        traceSSth = traceSSth + (s12*s12/2.d0 + s13*s13/2.d0 + s23*s23/2.d0 + 2.d0*dPsi*dPsi/3.d0 &
                  - 2*dPsi*dPhi/3.d0 + 2*dPhi*dPhi/3.d0)*div
        wsw = wsw + s12*omegay*omegax/4.d0 + s13*omegaz*omegax/4.d0 + s23*omegay*omegaz/4.d0   &
              + dPsi*omegax*omegax/12.d0 - dPsi*omegay*omegay/6.d0  + dPsi*omegaz*omegaz/12.d0 &
              + dPhi*omegax*omegax/12.d0 + dPhi*omegay*omegay/12.d0 - dPhi*omegaz*omegaz/6.d0 
        s3 = s3 - s12*s12*dPsi/4.d0 + s13*s13*dPsi/2.d0 - s23*s23*dPsi/4.d0 +3.d0/4.d0*s12*s13*s23 &
            + s12*s12*dPhi/2.d0 - s13*s13*dPhi/4.d0 - s23*s23*dPhi/4.d0 &
            - 2.d0/9.d0 * dPsi**3 + dPsi**2*dPhi /3.d0 + dPsi*dPhi**2/3.d0 - 2.d0/9.d0 * dPhi**3
        !
        m11m11  = m11m11 + du11**2
        m22m22  = m22m22 + du22**2
        m33m33  = m33m33 + du33**2
        m11m11m11  = m11m11m11 + du11**3
        m22m22m22  = m22m22m22 + du22**3
        m33m33m33  = m33m33m33 + du33**3
        !
        m12m12  = m12m12 + du12*du12
        m13m13  = m13m13 + du13*du13
        m21m21  = m21m21 + du21*du21
        m23m23  = m23m23 + du23*du23
        m31m31  = m31m31 + du31*du31
        m32m32  = m32m32 + du32*du32
        !
        m11m12m12=m11m12m12+du11*du12*du12
        m11m12m21=m11m12m21+du11*du12*du21
        m11m11m22=m11m11m22+du11*du11*du22
        !
        dudx2=dudx2+du11**2+du22**2+du33**2
        !
        miudrho=miudrho+miu/rho(i,j,k)
        !
        csavg=csavg+cs
        !
        machrms=machrms+v2/(cs*cs)
        !
        rhoe=rhoe+tmp(i,j,k)
        !
        dissloc = miu/rho(i,j,k)*(du11**2 + du12**2 + du13**2 + du21**2 + du22**2 + du23**2 + &
                  du31**2 + du32**2+ du33**2)
        !dissloc = 2.d0*miu/rho(i,j,k)*(s11**2+s22**2+s33**2+(s12**2+s13**2+s23**2)/2.d0-div**2/3.d0)
        !
        dissa=dissa+dissloc
        !
        kolmloc = min(kolmloc,sqrt(sqrt((miu/rho(i,j,k))**3/dissloc)))
        !
        ens=ens+(omegaz*omegaz)
        !
        du3=du3+(du11*du11*du11+du22*du22*du22+du33*du33*du33)
        du2=du2+(du11*du11+du22*du22+du33*du33)
        !
        rhoavg = rhoavg + rho(i,j,k)
        rho2nd = rho2nd + rho(i,j,k)*rho(i,j,k)
        w2drho = w2drho + (omegaz*omegaz)/rho(i,j,k)
        !
        mfpath = max(mfpath,2.d0*miu/rho(i,j,k)/0.921/sqrt(3.d0*R*tmp(i,j,k)))
        !
      enddo
      enddo
      enddo
      !
      urms  = sqrt(psum(urms)/rsamples)
      energy  = psum(energy)/rsamples
      !
      OmegaOmega  = psum(OmegaOmega)/rsamples
      divUdivU  = psum(divUdivU)/rsamples
      traceSS  = psum(traceSS)/rsamples
      !
      OmegaOmegadivU  = psum(OmegaOmegadivU)/rsamples
      divUdivUdivU  = psum(divUdivUdivU)/rsamples
      traceSSth  = psum(traceSSth)/rsamples
      wsw = psum(wsw)/rsamples
      s3 = psum(s3)/rsamples
      !
      m11m11  = psum(m11m11)/rsamples
      m22m22  = psum(m22m22)/rsamples
      m33m33  = psum(m33m33)/rsamples
      m11m11m11  = psum(m11m11m11)/rsamples
      m22m22m22  = psum(m22m22m22)/rsamples
      m33m33m33  = psum(m33m33m33)/rsamples
      !
      m12m12  = psum(m12m12)/rsamples
      m13m13  = psum(m13m13)/rsamples
      m21m21  = psum(m21m21)/rsamples
      m23m23  = psum(m23m23)/rsamples
      m31m31  = psum(m31m31)/rsamples
      m32m32  = psum(m32m32)/rsamples
      !
      m11m11m22  = psum(m11m11m22)/rsamples
      m11m12m12  = psum(m11m12m12)/rsamples
      m11m12m21  = psum(m11m12m21)/rsamples
      !
      dudx2      = num1d3*psum(dudx2)/rsamples
      miudrho    = psum(miudrho)/rsamples
      csavg      = psum(csavg)/rsamples
      dissa      = psum(dissa)/rsamples
      !
      machrms=sqrt(psum(machrms)/rsamples)
      !
      rhoe=psum(rhoe)/rsamples
      !
      ens=0.5d0*psum(ens)/rsamples
      !
      ufluc=urms/sqrt(2.d0)
      !
      macht         = urms/csavg
      taylorlength  = ufluc/sqrt(dudx2)
      retaylor      = ufluc*taylorlength/miudrho
      kolmoglength  = sqrt(sqrt(miudrho**3/dissa))
      kolmloc       = pmin(kolmloc)
      comlen        = pmin(comlen) ! Compressible length
      mfpath        = pmax(mfpath)
      !
      skewness      = psum(du3)/(3.d0*rsamples)/sqrt((psum(du2)/(3.d0*rsamples))**3)
      ! kolmogvelocity= sqrt(sqrt(dissipation*miudrho))
      ! kolmogtime    = sqrt(miudrho/dissipation)
      rhoavg = psum(rhoavg)/rsamples
      rho2nd = psum(rho2nd)/rsamples
      w2drho = psum(w2drho)/rsamples
      !
      if(lio) then 
        call listwrite(hand_fs,urms,ens,taylorlength,kolmoglength, &
                        kolmloc, Retaylor,comlen,mfpath,machrms, &
                        macht, rhoe,hsource,skewness)
        call listwrite(hand_skew,m11m11,m22m22,m33m33,m11m11m11,m22m22m22, &
                        m33m33m33,divUdivU,OmegaOmega,traceSS,divUdivUdivU, &
                        OmegaOmegadivU,traceSSth,wsw,s3)
        call listwrite(hand_mom2nd,m12m12,m13m13,m21m21,m23m23,m31m31,m32m32)
        call listwrite(hand_mom3rd,m11m12m12,m11m12m21,m11m11m22)
      endif
    endif
    !
  end subroutine udf_stalist
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_stalist.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add vortical fluctuations to initial field  | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  ! subroutine addvortex(xc,yc,radius,amp)
    !
    ! use commvar,  only: im,jm,km,ndims,roinf,uinf
    ! use parallel, only: lio
    ! use commarray,only: x,vel,rho,prs,spc,tmp,q
    ! use fludyna,  only: thermal
    ! !
    ! ! local data
    ! real(8),intent(in) :: xc,yc,radius,amp
    ! !
    ! integer :: i,j,k
    ! real(8) :: var1,radi2,cvor
    ! !
    ! cvor=amp*uinf*radius
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/radius/radius
    !   var1=cvor/radius/radius*exp(-0.5d0*radi2)
    !   !
    !   vel(i,j,k,1)=vel(i,j,k,1)-var1*(x(i,j,k,2)-yc)
    !   if(ndims>=2) vel(i,j,k,2)=vel(i,j,k,2)+var1*(x(i,j,k,1)-xc)
    !   if(ndims==3) vel(i,j,k,3)=0.d0
    !   prs(i,j,k)  =prs(i,j,k)-0.5d0*roinf*cvor*cvor/radi2/radi2*exp(-radi2)
    !   !
    !   tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
    !   !
    ! enddo
    ! enddo
    ! enddo
    !
  ! end subroutine addvortex
  !+-------------------------------------------------------------------+
  !| The end of the subroutine addvortex.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| hit flame.                                                        |
  !| a random force acting like fans to input energy at largest scale  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-06-2023: Created by Yifan Xu @ Peking University               |
  !+-------------------------------------------------------------------+
  subroutine udf_src
    !
    use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,rkstep,xmax,ymax,zmax,&
                         lforce,nstep, forcenum
    use parallel, only : lio,psum,bcast
    use commarray,only : rho,tmp,vel,qrhs,x,jacob,forcep,forcek
    use utility,  only : listinit,listwrite
    use constdef, only : pi
    !
    logical,save :: linit=.true.
    integer,save :: hand_force
    integer, allocatable, dimension(:), save :: hand_forcea
    integer,save :: last_step = -1
    ! Random iniforce generation
    integer :: NumTheta, n, i,j,k,t
    real(8) :: theta
    real(8) :: power,rsamples,Tpower
    real(8) :: forcekT
    real(8), allocatable,dimension(:) :: alphas, alphad
    character(len=4) :: forcename
    !
    if(.not. allocated(hand_forcea))then
      allocate(hand_forcea(1:forcenum))
    else
      if(size(hand_forcea) .ne. forcenum) stop 'Error in hand_forcea @ udf_src'
    endif
    !
    allocate(alphas(1:forcenum),alphad(1:forcenum))
    !
    if(lforce) then
      !
      if(linit) then
        !
        if(lio) then
          !
          call listinit(filename='log/forcestat.dat',handle =hand_force, firstline='nstep time forcekT power')
          !
          do t=1,forcenum
            write(forcename,'(i4.4)')t
            call listinit(filename='log/forcestat'//forcename//'.dat',handle =hand_forcea(t), &
             firstline='nstep time forcek alphas alphad')
          enddo
          !
        endif
        !
        linit = .false.
        !
      endif
      !
      if(last_step .ne. nstep) then
        !
        call udf_generate_force(alphas,alphad)
        ! Calculate power
        power = 0.0d0
        Tpower = 0.0d0
        !
        if(ndims == 2)then
          do j=1,jm
          do i=1,im
            power  = power  + rho(i,j,0)*(forcep(i,j,0,1)*vel(i,j,0,1) + forcep(i,j,0,2)*vel(i,j,0,2))
            Tpower = Tpower + tmp(i,j,0)**4
          enddo
          enddo
          !
          rsamples=dble(ia*ja)
        endif
        !
        power = psum(power)/rsamples
        Tpower = psum(Tpower)/rsamples
        forcekT = power/Tpower
        !
        if(lio) then 
          call listwrite(hand_force,forcekT,power)
          !
          do t=1,forcenum
            call listwrite(hand_forcea(t),real(forcek(t),8),alphas(t),alphad(t))
          enddo
        endif
        !
        last_step = nstep
        !
      endif
      !
      ! Add in qrhs and calculate power
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        !
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+rho(i,j,k)*forcep(i,j,k,1)*jacob(i,j,k)
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*forcep(i,j,k,2)*jacob(i,j,k)
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+rho(i,j,k)*forcep(i,j,k,3)*jacob(i,j,k)
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+rho(i,j,k)*(forcep(i,j,k,1)*vel(i,j,k,1) + &
                                                forcep(i,j,k,2)*vel(i,j,k,2) + &
                                                forcep(i,j,k,3)*vel(i,j,k,3) )*jacob(i,j,k)
        !
        !
        ! temperation dissipation
        qrhs(i,j,k,5)=qrhs(i,j,k,5)-forcekT*(tmp(i,j,k)**4)*jacob(i,j,k)
        !
      end do
      end do
      end do
      !
      !
    endif
    !
  end subroutine udf_src
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_src.                                |
  !+-------------------------------------------------------------------+
  !
  subroutine udf_generate_force(alphas, alphad)
    !
    use commvar, only: ndims
    real(8), allocatable, dimension(:), intent(out) :: alphas, alphad
    !   
    if(ndims == 2) then
      call udf_generate_force_2D(alphas, alphad)
    else
      stop "Not implemented error! udf_generate_force3D"
    endif
    !
  end subroutine udf_generate_force
    !
  subroutine udf_generate_force_2D(alphas, alphad)
    !
    use, intrinsic :: iso_c_binding
    use commvar,        only : forcenum,hypervisk,hypervismiu,im,jm,ia,ja
    use commarray,      only : vel, forcep,forcek,forcespes,forcesped
    use fftwlink,       only : jmf, alloc_local, iafftw, jmfftw, fftw_grid_fence, fftw_fence_grid, jafftw, j0f
    use parallel,       only : MPI_COMM_WORLD, psum, mpiright, mpiup, mpitag, mpileft, mpidown
    use udf_tool,       only : GenerateWave, kint
    use mpi
    include 'fftw3-mpi.f03'
    !
    real(8), allocatable, dimension(:), intent(out):: alphas, alphad
    real(8), allocatable, dimension(:,:) :: localvel1t, localvel2t, force1t, force2t
    real(8), allocatable, dimension(:,:) :: fftvel1, fftvel2, fftforce1, fftforce2
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_fftforce1, c_fftforce2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8), allocatable, dimension(:,:) :: usspe,udspe,u1s,u2s,u1d,u2d
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: fftforce1c,fftforce2c
    real(8) ::  dk, kk
    integer :: i,j,ierr,allkmax,kOrdinal,t
    real(8), allocatable, dimension(:) :: Ed, Es
    real(8), allocatable, dimension(:) :: sendim,recvim, sendjm,recvjm
    integer :: status(mpi_status_size) 
    !
    dk = 1.d0
    !
    allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
    !
    if(.not. allocated(alphas)) allocate(alphas(1:forcenum))
    if(.not. allocated(alphad)) allocate(alphad(1:forcenum))
    allocate(localvel1t(1:jm,1:im),localvel2t(1:jm,1:im))
    allocate(force1t(1:jm,1:im),force2t(1:jm,1:im))
    allocate(fftvel1(1:ia,1:jmf),fftvel2(1:ia,1:jmf))
    allocate(fftforce1(1:ia,1:jmf),fftforce2(1:ia,1:jmf))
    !
    do j=1,jm
    do i=1,im
      !
      localvel1t(j,i)=vel(i,j,0,1)
      localvel2t(j,i)=vel(i,j,0,2)
      !
    enddo
    enddo
    !
    !
    fftvel1 = fftw_grid_fence(localvel1t)
    fftvel2 = fftw_grid_fence(localvel2t)
    !
    ! Begin FFTW
      !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [iafftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [iafftw,jmfftw])
    !
    !!!! Do S-C decomposition
    c_fftforce1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_fftforce1, fftforce1c, [iafftw,jmfftw])
    c_fftforce2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_fftforce2, fftforce2c, [iafftw,jmfftw])
    !
    ! planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !
    do j=1,jmf
    do i=1,ia
        !
        u1spe(i,j)=CMPLX(fftvel1(i,j),0.d0,C_INTPTR_T);
        u2spe(i,j)=CMPLX(fftvel2(i,j),0.d0,C_INTPTR_T);
        !
    end do
    end do
    !
    !
    !!!! Do 2d FFT
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jmf
    do i=1,ia
        !
        u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
        u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
        !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:ia,1:jmf),k2(1:ia,1:jmf))
    call GenerateWave(ia,jmf,ia,ja,j0f,k1,k2)
    !
    allocate(usspe(1:ia,1:jmf),udspe(1:ia,1:jmf))
    allocate(u1s(1:ia,1:jmf),u1d(1:ia,1:jmf),u2s(1:ia,1:jmf),u2d(1:ia,1:jmf))
    allocate(Es(0:allkmax),Ed(0:allkmax))
    Ed = 0.d0
    Es = 0.d0
    !
    do j=1,jmf
    do i=1,ia
        kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
        udspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
        u1d(i,j)=  udspe(i,j)*k1(i,j)/kk
        u2d(i,j)=  udspe(i,j)*k2(i,j)/kk
        u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
        kOrdinal = kint(kk,dk,2,1.d0)
        if(kOrdinal <= allkmax)then
          Es(kOrdinal) = Es(kOrdinal) + usspe(i,j)*conjg(usspe(i,j))/2
          Ed(kOrdinal) = Ed(kOrdinal) + udspe(i,j)*conjg(udspe(i,j))/2
        endif
        !
    end do
    end do
    !
    !
    do i=0,allkmax
      Es(i) = psum(Es(i))
      Ed(i) = psum(Ed(i))
    enddo
    !
    do t=1,forcenum
      alphas(t) = min(max(forcespes(t)/Es(forcek(t))-1.d0, 0.d0),10.d0)
      alphad(t) = min(max(forcesped(t)/Ed(forcek(t))-1.d0, 0.d0),10.d0)
    enddo
    !
    !
    fftforce1c = 0.d0
    fftforce2c = 0.d0
    do j=1,jmf
    do i=1,ia
      kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
      do t=1,forcenum
        if(kint(kk,dk,2,1.d0)==forcek(t))then
          fftforce1c(i,j) = fftforce1c(i,j) + alphas(t) * u1s(i,j) + alphad(t) * u1d(i,j)
          fftforce2c(i,j) = fftforce2c(i,j) + alphas(t) * u2s(i,j) + alphad(t) * u2d(i,j)
        endif
      enddo
      !
      if(kk>hypervisk)then
        fftforce1c(i,j) = fftforce1c(i,j) - hypervismiu * kk**4 * u1spe(i,j)
        fftforce2c(i,j) = fftforce2c(i,j) - hypervismiu * kk**4 * u2spe(i,j)
      endif
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(backward_plan,fftforce1c,fftforce1c)
    call fftw_mpi_execute_dft(backward_plan,fftforce2c,fftforce2c)
    !
    do j=1,jmf
    do i=1,ia
      fftforce1(i,j) = real(fftforce1c(i,j))
      fftforce2(i,j) = real(fftforce2c(i,j))
    enddo
    enddo
    !
    force1t = fftw_fence_grid(fftforce1)
    force2t = fftw_fence_grid(fftforce2)
    !
    forcep = 0.d0
    !
    !
    do j=1,jm
    do i=1,im
      forcep(i,j,0,1) = force1t(j,i)
      forcep(i,j,0,2) = force2t(j,i)
    enddo
    enddo
    !
    call mpi_sendrecv(forcep(im,jm,0,1),1,mpi_real8,mpiright,mpitag, &
                      forcep(0,jm,0,1),1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,2),1,mpi_real8,mpiright,mpitag, &
                      forcep(0,jm,0,2),1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,1),1,mpi_real8,mpiup,mpitag,    &
                      forcep(im,0,0,1),1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,2),1,mpi_real8,mpiup,mpitag,    &
                      forcep(im,0,0,2),1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    allocate(sendim(0:im),recvim(0:im),sendjm(0:jm),recvjm(0:jm))
    !
    sendim(0:im) = forcep(0:im,jm,0,1)
    call mpi_sendrecv(sendim,im+1,mpi_real8,mpiup,mpitag,     &
                      recvim,im+1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    forcep(0:im,0,0,1) = recvim(0:im)
    mpitag=mpitag+1
    !
    sendim(0:im) = forcep(0:im,jm,0,2)
    call mpi_sendrecv(sendim,im+1,mpi_real8,mpiup,mpitag,     &
                      recvim,im+1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    forcep(0:im,0,0,2) = recvim(0:im)
    !
    sendjm(0:jm) = forcep(im,0:jm,0,1)
    call mpi_sendrecv(sendjm,jm+1,mpi_real8,mpiright,mpitag,  &
                      recvjm,jm+1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    forcep(0,0:jm,0,1) = recvjm(0:jm)
    !
    sendjm(0:jm) = forcep(im,0:jm,0,2)
    call mpi_sendrecv(sendjm,jm+1,mpi_real8,mpiright,mpitag,  &
                      recvjm,jm+1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    forcep(0,0:jm,0,2) = recvjm(0:jm)
    mpitag=mpitag+1
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_fftforce1)
    call fftw_free(c_fftforce2)
    !
    deallocate(localvel1t, localvel2t, force1t, force2t)
    deallocate(fftvel1, fftvel2, fftforce1, fftforce2)
    deallocate(k1,k2,usspe,udspe, u1d, u1s, u2d, u2s)
    deallocate(sendim,sendjm,recvim,recvjm)
    !
  end subroutine udf_generate_force_2D
  !+-------------------------------------------------------------------+
  !| This subroutine is to defined an output by a user.                | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_write
    !
    
    !
  end subroutine udf_write
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_write.                              |
  !+-------------------------------------------------------------------+
  !
    !+-------------------------------------------------------------------+
  !| This subroutine is to manipulate data solver as one likes at the  |
  !| end of each loop.                                                 | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-Oct-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_eom_set
    !
  end subroutine udf_eom_set
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_eom_set.                            |
  !+-------------------------------------------------------------------+
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
