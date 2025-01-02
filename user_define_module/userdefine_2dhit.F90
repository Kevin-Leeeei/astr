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
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km
    use commarray,only : vel,rho,tmp,dvel,q
    use fludyna,  only : miucal,sos
    use utility,  only : listinit,listwrite
    use parallel, only : psum,lio,pmax,pmin
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
    real(8) :: urms,energy,taylorlength,kolmoglength,Retaylor,machrms,macht,rhoe,skewness,du2,du3
    real(8) :: mfpath,R ! mfpath = mean free path
    !
    logical :: fex
    logical,save :: linit=.true.
    integer,save :: hand_fs,hand_mom2nd,hand_mom3rd,hand_skew,hand_en
    !
    if(ka==0) then 
      R = 8.31446261815324d0
      !
      if(linit) then
        !
        if(lio) then
          call listinit(filename='log/fturbstats.dat',handle=hand_fs, &
             firstline='nstep time urms ens talen kolmavg kolmloc Reta comlen mfpath machrms macht Tavg hsource skewness')
          call listinit(filename='log/mom2nd.dat',handle=hand_mom2nd, &
                        firstline='nstep time m11m11 m22m22 m11m22 m12m21 m12m12 m21m21')
          call listinit(filename='log/mom3rd.dat',handle=hand_mom3rd, &
                        firstline='nstep time A111111 A222222 A111122 A222211 A111212 A221212 A112121 A222121 A222112 A111221')
          call listinit(filename='log/skewness.dat',handle=hand_skew, &
                        firstline='ns time th o ps f th2 ps2 f2 oth o2 ps2th f2th th3 o2th m11s m22s m11c m22c')
          call listinit(filename='log/fenergy.dat',handle=hand_en, &
                        firstline='ns time rhoavg rho2nd w2drho')
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
      !
      k=0
      ! do k=1,km
      do j=1,jm
      do i=1,im
        !
        ! This point values
        !
        miu=miucal(tmp(i,j,k))/Reynolds
        !
        du11=dvel(i,j,k,1,1)
        du12=dvel(i,j,k,1,2)
        du21=dvel(i,j,k,2,1)
        du22=dvel(i,j,k,2,2)
        !
        s11=du11
        s12=0.5d0*(du12+du21)
        !s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
        s22=du22
        !s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
        !s33=dvel(i,j,k,3,3)
        !
        div=du11+du22
        !
        !omegax=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
        !omegay=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
        omegaz=du21-du12
        !
        dPsi = du12+du21
        dPhi = du11-du22
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
        divU = divU+div
        comlen = min(comlen,sqrt(2.d0*4.d0/3.d0*miu/rho(i,j,k)/abs(div)))
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
        dissa=dissa+2.d0*miu/rho(i,j,k)*(s11**2+s22**2+2.d0*(s12**2)-0.5d0*div**2)
        !
        dissloc = 2.d0*miu/rho(i,j,k)*(s11**2+s22**2+2.d0*(s12**2)-0.5d0*div**2)
        kolmloc = min(kolmloc,sqrt(sqrt((miu/rho(i,j,k))**3/dissloc)))
        !
        ens=ens+(omegaz*omegaz)
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
      enddo
      enddo
      ! enddo
      urms  = sqrt(psum(urms)/rsamples)
      energy  = psum(energy)/rsamples
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
      ens=0.5d0*psum(ens)/rsamples
      !
      ufluc=urms/sqrt(2.d0)
      !
      macht         = urms/csavg!wrong
      taylorlength  = ufluc/sqrt(dudx2)
      retaylor      = ufluc*taylorlength/miudrho
      kolmoglength  = sqrt(sqrt(miudrho**3/dissa))
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
      if(lio) then 
        call listwrite(hand_fs,urms,ens,taylorlength,kolmoglength, &
                        kolmloc, Retaylor,comlen,mfpath,machrms, &
                        macht, rhoe,hsource,skewness)
        call listwrite(hand_mom2nd,m11m11,m22m22,m11m22,m12m21,m12m12,m21m21)
        call listwrite(hand_mom3rd,m11m11m11,m22m22m22,m11m11m22,m22m22m11,&
                      m11m12m12,m22m12m12,m11m21m21,m22m21m21,m22m21m12,m11m12m21)
        call listwrite(hand_skew,divU,Omega,Psi,Phi,divUdivU,dPsidPsi,&
                      dPhidPhi,OmegadivU,OmegaOmega,dPsidPsidivU,dPhidPhidivU,&
                      divUdivUdivU,OmegaOmegadivU,m11m11,m22m22,m11m11m11,m22m22m22)
        call listwrite(hand_en,rhoavg,rho2nd,w2drho)
      endif
    else
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
                        lforce,nstep,forcek,forceamp,numftheta,forcekT
    use parallel, only : lio,mpirank,psum,bcast
    use commarray,only : rho,tmp,vel,qrhs,x,jacob,Ak,phik 
    use utility,  only : listinit,listwrite
    use constdef,only: pi
    !
    logical,save :: linit=.true.
    integer,save :: hand_force
    integer,save :: last_step
    ! Random iniforce generation
    integer :: NumTheta, n, i,j,k
    real(8) :: theta, k1, k2
    real(8) :: power,rsamples
    real(8),allocatable,dimension(:,:,:,:) :: force
    !
    NumTheta = 200
    last_step = -1
    !
    allocate(force(0:im,0:jm,0:0,1:2))
    !
    !
    if(lforce) then
      !
      if(linit) then
        !
        if(lio) then
          !
          call listinit(filename='log/forcestat.dat',handle =hand_force, firstline='nstep time power')
          !
        endif
        !
        linit = .false.
        !
      endif
      !
      !      
      ! Generate random force
      do j=1,jm
      do i=1,im
        force(i,j,0,1) = 0.d0
        force(i,j,0,2) = 0.d0
      end do
      end do
      !
      !
      !
      do n=0, (NumTheta-1)
        theta = n * 2 * pi / NumTheta
        ! 
        k1 = forcek * dcos(theta)
        k2 = forcek * dsin(theta)
        do j=1,jm
        do i=1,im
          force(i,j,0,1) = force(i,j,0,1) + forceamp * Ak(n) * k2 * dcos(k1 * x(i,j,0,1) + k2 * x(i,j,0,2) + phik(n))
          force(i,j,0,2) = force(i,j,0,2) - forceamp * Ak(n) * k1 * dcos(k1 * x(i,j,0,1) + k2 * x(i,j,0,2) + phik(n))
        end do
        end do
      end do
      !
      force(0,:,0,1) = force(im,:,0,1)
      force(0,:,0,2) = force(im,:,0,2)
      force(:,0,0,1) = force(:,jm,0,1)
      force(:,0,0,2) = force(:,jm,0,2)
        !
      if(last_step .ne. nstep) then
        ! Calculate power
        power = 0.0d0
        !
        do j=0,jm
        do i=0,im
          power=power + force(i,j,0,1)*vel(i,j,0,1) + force(i,j,0,2)*vel(i,j,0,2)
        enddo
        enddo
        !
        rsamples=dble(ia*ja)
        power= psum(power)/rsamples
        !
        if(lio) then 
          call listwrite(hand_force,abs(power))
        endif
        !
        last_step = nstep
        !
      endif
      !
      ! Inverse if the power is negative
      if (power<0) then
        force(:,:,:,:) = - force(:,:,:,:)
      endif
      !
      ! Add in qrhs and calculate power
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        !
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+rho(i,j,k)*force(i,j,k,1)*jacob(i,j,k)
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*force(i,j,k,2)*jacob(i,j,k)
        !
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+rho(i,j,k)*( force(i,j,k,1)*vel(i,j,k,1) + &
                                                force(i,j,k,2)*vel(i,j,k,2) )*jacob(i,j,k)
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
      deallocate(force)
      !
    endif
    !
  end subroutine udf_src
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_src.                                |
  !+-------------------------------------------------------------------+
  !
  subroutine udf_generate_Aphi
    !
    use commvar, only: numftheta,lforce
    use commarray, only: Ak, phik 
    use parallel, only : lio,mpirank,psum,bcast
    use constdef,only: pi
    !
    real(8) :: ranAk, ranphik
    integer :: n,seedsize,clock,irandom
    integer,allocatable :: seed(:)
    !
    if(lforce) then
      if(allocated(Ak)) then
        deallocate(Ak)
      endif

      if(allocated(phik)) then
        deallocate(phik)
      endif
      !
      allocate(Ak(0:(numftheta-1)),phik(0:(numftheta-1)))
      !
      if(mpirank==0) then
        !
        call random_seed(size=seedsize) ! find out size of seed
        allocate(seed(seedsize))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, seedsize) /)
        call random_seed(put=seed) ! set current seed
        deallocate(seed)           ! safe
        !
        do n=0, (numftheta-1)
          ! Generate random amplitude and phase
          !
          call random_number(ranAk)
          call random_number(ranphik)
          Ak(n) = 1 + (ranAk-0.5)*0.1
          phik(n) = 2*pi*ranphik
        end do
      end if
      !
      call bcast(Ak)
      call bcast(phik)
    endif
    !
  end subroutine udf_generate_Aphi
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
