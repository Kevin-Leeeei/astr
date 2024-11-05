!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| physical SGS method for energy flux calculation.                    |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_hitgen
    !
    !
    use constdef
    use stlaio,  only: get_unit
    !
    implicit none
    !
    contains
    !
    subroutine ppHitgenentrance
        !
        use cmdefne
        use parallel,        only : mpirank,bcast,mpisize,lio
        !
        ! local data
        character(len=64) :: casefolder,inputfile,outputfile,viewmode, &
                             flowfieldfile, readmode
        integer :: filenumb
        !
        !
        if(mpirank == 0) then
          call readkeyboad(readmode) 
        endif
        call bcast(readmode)
        !
        if(trim(readmode)=='3DIC4') then
          !
          if(mpisize .ne. 1) then
            if(mpirank == 0)   print *, ' Parallel used, generate IC4'
            call hitgen_ic4_parallel
          else
            print *, 'Serial mode, only use 1 proc, generate IC4'
            call hitgen_ic4
          endif
          !
        elseif(trim(readmode)=='pic') then
          !
          ! TODO: implement
          stop 'NOT IMPLEMENTED ERROR'
          !
          !
        elseif(trim(readmode)=='2DIC4') then
          !
          if(mpisize .ne. 1) then
            if(mpirank == 0)   print *, ' Parallel used, generate IC4'
            call hitgen2d_ic4_parallel
          else
            print *, 'Serial mode, only use 1 proc, generate IC4'
            call hitgen2d_ic4
          endif
          !
        elseif(trim(readmode)=='2Dpic') then
          !
          if(mpirank == 0)   print *, 'Generate one pic field'
          call hitgen2d_pic
          !
        elseif(trim(readmode)=='2Dherring') then
          !
          if(mpirank == 0)   print *, 'Generate Herrings field'
          call hitgen2d_herring
          !
        else
          print* ,"Readmode is not defined!", readmode
        endif
    end subroutine ppHitgenentrance
    !
    !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a random field for homogeneous|
  !| isotropic turbulence.                                             |
  !+-------------------------------------------------------------------+
  !| Ref: Blaisdell, G. A., Numerical simulation of compressible       |
  !|      homogeneous turbulence, Phd, 1991, Stanford University       |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-04-2023: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine hitgen_ic4
    !
    use cmdefne,   only : readkeyboad
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen
    use bc,        only : twall
    use readwrite, only : readgrid, readic, readinput
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : mpisizedis,parapp,parallelini
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io,    only : h5srite,h5sread
    use gridgeneration
    use tecio
    !
    real(8) :: urms,kenergy,ufmx,roav,uav,vav,wav,tav,pav
    integer :: i,j,k
    character(len=4) :: genmethod
    !
    print*,' ** cmd to generate a box turbulence: astr pp hitgen <input> . gen/read'
    !
    !
    call readinput
    !
    call readic
    !
    im=ka ! ensure a cubic box
    !
    jm=ka
    !
    km=ka
    !
    call mpisizedis
    !
    call parapp
    !
    call parallelini
    !
    call refcal
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(rho(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km),prs(0:im,0:jm,0:km))
#ifdef COMB
    allocate(spc(0:im,0:jm,0:km,1:num_species))
#endif
    !
    ! call gridhitflame(mode='cubic')
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    ! call readgrid(trim(gridfile))
    !
    call geomcal
    !
    call readkeyboad(genmethod)
    !
    if(trim(genmethod)=='gen') then
      !
      call div_free_gen(im,jm,km,vel(0:im,0:jm,0:km,1),   &
                                 vel(0:im,0:jm,0:km,2),   &
                                 vel(0:im,0:jm,0:km,3) )
      !
      urms=1.d0 
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        rho(i,j,k)  = roinf
        tmp(i,j,k)  = tinf 
        if(nondimen) then
            prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        else
            spc(i,j,k,:)= spcinf(:)
            prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
        endif
        vel(i,j,k,1)= urms*vel(i,j,k,1)
        vel(i,j,k,2)= urms*vel(i,j,k,2)
        vel(i,j,k,3)= urms*vel(i,j,k,3)
        !
      enddo
      enddo
      enddo
      !
      call div_test(vel,dvel)
      !
      call hitsta
      !
    elseif(trim(genmethod)=='read') then
      !
      call h5sread(vel(0:im,0:jm,0:km,1),'u1',im,jm,km,'outdat/flowfield.h5')
      call h5sread(vel(0:im,0:jm,0:km,2),'u2',im,jm,km,'outdat/flowfield.h5')
      call h5sread(vel(0:im,0:jm,0:km,3),'u3',im,jm,km,'outdat/flowfield.h5')
      call h5sread(rho(0:im,0:jm,0:km),  'ro',im,jm,km,'outdat/flowfield.h5')
      call h5sread(tmp(0:im,0:jm,0:km),  't', im,jm,km,'outdat/flowfield.h5')
      !
      call div_test(vel,dvel)
      !
      call hitsta
      !
    else
      print*,genmethod
      stop ' !! genmethod error '
    endif
    ! 
    call h5srite(var=rho,                  varname='ro',filename='flowini3d.h5',explicit=.true.,newfile=.true.) 
    call h5srite(var=vel(0:im,0:jm,0:km,1),varname='u1',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,2),varname='u2',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,3),varname='u3',filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=prs,                  varname='p', filename='flowini3d.h5',explicit=.true.)
    call h5srite(var=tmp,                  varname='t', filename='flowini3d.h5',explicit=.true.)
    !
    roav=0.d0
    uav=0.d0
    vav=0.d0
    wav=0.d0
    tav=0.d0
    pav=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im 
      if(nondimen) then
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      else
        spc(i,j,k,:)= spcinf(:)
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
      endif
      !
      if(i==0 .or. j==0 .or. k==0) cycle
      !
      roav=roav+rho(i,j,k)
      uav=uav+vel(i,j,k,1)
      vav=vav+vel(i,j,k,2)
      wav=wav+vel(i,j,k,3)
      tav=tav+tmp(i,j,k)
      pav=pav+prs(i,j,k)
      !
    enddo
    enddo
    enddo
    !
    roav=roav/dble(im*jm*km)
    uav = uav/dble(im*jm*km)
    vav = vav/dble(im*jm*km)
    wav = wav/dble(im*jm*km)
    tav = tav/dble(im*jm*km)
    pav = pav/dble(im*jm*km)
    !
    print*, '** mean density    :',roav
    print*, '** mean velocity   :',uav,vav,wav
    print*, '** mean temperature:',tav
    print*, '** mean pressure   :',pav
    !
    rho(:,:,:)   = rho(:,:,:)  -roav
    vel(:,:,:,1) = vel(:,:,:,1)-uav
    vel(:,:,:,2) = vel(:,:,:,2)-vav
    vel(:,:,:,3) = vel(:,:,:,3)-wav
    tmp(:,:,:)   = tmp(:,:,:)  -tav
    prs(:,:,:)   = prs(:,:,:)  -pav
    !
    call h5srite(var=x(0:im,0,0,1),        varname='x', filename='datin/flowin.h5',explicit=.true.,newfile=.true.) 
    call h5srite(var=rho,                  varname='ro',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,1),varname='u1',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,2),varname='u2',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=vel(0:im,0:jm,0:km,3),varname='u3',filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=prs,                  varname='p', filename='datin/flowin.h5',explicit=.true.)
    call h5srite(var=tmp,                  varname='t', filename='datin/flowin.h5',explicit=.true.)
    
    ! call tecbin('techit.plt',x(0:im,0:jm,0:km,1),'x', &
    !                          x(0:im,0:jm,0:km,2),'y', &
    !                          x(0:im,0:jm,0:km,3),'z', &
    !                        vel(0:im,0:jm,0:km,1),'u', &
    !                        vel(0:im,0:jm,0:km,2),'v', &
    !                        vel(0:im,0:jm,0:km,3),'w' )
    !
  end subroutine hitgen_ic4
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen.                                 |
  !+-------------------------------------------------------------------+
  !
  subroutine hitgen_ic4_parallel
    !
    use, intrinsic :: iso_c_binding
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use readwrite, only : readgrid, readic, readinput
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    real(8) :: urms,kenergy,ufmx,roav,uav,vav,wav,tav,pav,vmax
    integer :: i,j,k,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2,wn3, wn12, wna, var1, var2
    real(8) :: ran1, ran2,ran3,ran4,rn3
    complex(8) :: vac1, vac2, vac3, crn1, crn2, crn4
    real(8) :: ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1c,u2c,u3c
    real(C_DOUBLE), pointer, dimension(:,:,:) :: u1r,u2r,u3r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u3c, c_u1r, c_u2r, c_u3r
    character(len=1) :: modeio
    modeio = 'h'
    !
    !
    call readinput
    !
    call readic
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:km),tmp(0:(2*im),0:jm,0:km),prs(0:(2*im),0:jm,0:km))
    !
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j,k) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j,k) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if(j <= (ja/2+1)) then
        k2(i,j,k) = real(j-1,8)
      else if(j<=ja) then
        k2(i,j,k) = real(j-ja-1,8)
      else
        print *,"Error, no wave number possible, j must smaller than ja-1 !"
      end if
      !
      !
      if((k+k0) <= (ka/2+1)) then
        k3(i,j,k) = real(k+k0-1,8)
      else if((k+k0)<=(ka)) then
        k3(i,j,k) = real(k+k0-ka-1,8)
      else
        print *,"Error, no wave number possible, (k+k0) must smaller than ka-1 !"
      end if
      !
      !
    end do
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw,kmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw,kmfftw])
    c_u3r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u3r, u3r, [2*imfftw,jmfftw,kmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_3d(kafftw,jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
        u1c(i,j,k)=0.d0
        u2c(i,j,k)=0.d0
        u3c(i,j,k)=0.d0
      else
        !
        call random_number(ran1)
        call random_number(ran2)
        call random_number(ran3)
        call random_number(ran4)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        rn3 =ran3*2.d0*pi
        crn4=ran4*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j,k))
        wn2=real(k2(i,j,k))
        wn3=real(k3(i,j,k))
        wn12=sqrt(wn1**2+wn2**2)
        wna=sqrt(wn1**2+wn2**2+wn3**2)
        !
        var1=IniEnergDis(ISEA*2.d0,ickmax*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        !
        vac1=var2*cdexp(crn1)*dcos(rn3)
        vac2=var2*cdexp(crn2)*dsin(rn3)
        vac3=var2*cdexp(crn4)
        !
        if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
          u1c(i,j,k)=vac1*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=vac2*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=vac3*wn3/wna*icdilatational
        else
          u1c(i,j,k)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=-vac2*wn12/wna*icsolenoidal + vac3*wn3/wna*icdilatational
        end if
        !
      endif
    enddo
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u3c,u3r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do k=1,km
    do j=1,jm
    do i=1,im
        ! 
        vel(i,j,k,1)=u1r(i,j,k)
        vel(i,j,k,2)=u2r(i,j,k)
        vel(i,j,k,3)=u3r(i,j,k)
        !
    end do
    end do
    end do
    !
    vel(0,1:jm,1:km,1)=vel(im,1:jm,1:km,1)
    vel(0,1:jm,1:km,2)=vel(im,1:jm,1:km,2)
    vel(0,1:jm,1:km,3)=vel(im,1:jm,1:km,3)
    !
    vel(0:im,0,1:km,1)=vel(0:im,jm,1:km,1)
    vel(0:im,0,1:km,2)=vel(0:im,jm,1:km,2)
    vel(0:im,0,1:km,3)=vel(0:im,jm,1:km,3)
    ! !
    vel(0:im,0:jm,0:(km-1),1)=vel(0:im,0:jm,1:km,1)
    vel(0:im,0:jm,0:(km-1),2)=vel(0:im,0:jm,1:km,2)
    vel(0:im,0:jm,0:(km-1),3)=vel(0:im,0:jm,1:km,3)
    !
    call dataswap(vel)
    !
    roav=0.d0
    uav=0.d0
    vav=0.d0
    wav=0.d0
    tav=0.d0
    pav=0.d0
    vmax=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im 
      rho(i,j,k)  = roinf
      tmp(i,j,k)  = tinf 
      if(nondimen) then
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      else
        spc(i,j,k,:)= spcinf(:)
        prs(i,j,k)  = thermal(density=rho(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
      endif
      !
      vel(i,j,k,1)= icamplitude*vel(i,j,k,1)
      vel(i,j,k,2)= icamplitude*vel(i,j,k,2)
      vel(i,j,k,3)= icamplitude*vel(i,j,k,3)
      !
      if(i==0 .or. j==0 .or. k==0) cycle
      !
      roav=roav+rho(i,j,k)
      uav=uav+vel(i,j,k,1)
      vav=vav+vel(i,j,k,2)
      wav=wav+vel(i,j,k,3)
      tav=tav+tmp(i,j,k)
      pav=pav+prs(i,j,k)
      vmax = max(vmax,sqrt(vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2))
      !
    enddo
    enddo
    enddo
    !
    roav= psum(roav)/dble(ia*ja*ka)
    uav = psum(uav) /dble(ia*ja*ka)
    vav = psum(vav) /dble(ia*ja*ka)
    wav = psum(wav) /dble(ia*ja*ka)
    tav = psum(tav) /dble(ia*ja*ka)
    pav = psum(pav) /dble(ia*ja*ka)
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** mean density    :',roav
      print*, '** mean velocity   :',uav,vav,wav
      print*, '** mean temperature:',tav
      print*, '** mean pressure   :',pav
      print*, '** max velocity    :',vmax
    endif
    !
    vel(:,:,:,1) = vel(:,:,:,1)-uav
    vel(:,:,:,2) = vel(:,:,:,2)-vav
    vel(:,:,:,3) = vel(:,:,:,3)-wav
    ! 
    call mpi_barrier(mpi_comm_world,ierr)
    !
    call h5io_init(trim('datin/flowini3d.h5'),mode='write')
    !
    if((k0+km)==ka)then
      call h5write(var=rho(0:im,0:jm,0:km),  varname='ro', mode = modeio) 
      call h5write(var=vel(0:im,0:jm,0:km,1),varname='u1', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:km,2),varname='u2', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:km,3),varname='u3', mode = modeio)
      call h5write(var=prs(0:im,0:jm,0:km),  varname='p', mode = modeio)
      call h5write(var=tmp(0:im,0:jm,0:km),  varname='t', mode = modeio)
    else
      call h5write(var=rho(0:im,0:jm,0:(km-1)),  varname='ro', mode = modeio) 
      call h5write(var=vel(0:im,0:jm,0:(km-1),1),varname='u1', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:(km-1),2),varname='u2', mode = modeio)
      call h5write(var=vel(0:im,0:jm,0:(km-1),3),varname='u3', mode = modeio)
      call h5write(var=prs(0:im,0:jm,0:(km-1)),  varname='p', mode = modeio)
      call h5write(var=tmp(0:im,0:jm,0:(km-1)),  varname='t', mode = modeio)
    end if
    call h5io_end
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    !
    call geomcal
    !
    call div_test(vel,dvel)
    ! !
    ! call hitsta!
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call fftw_free(c_u3r)
    call mpistop
    !
  end subroutine hitgen_ic4_parallel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  !| Ref: Samtaney, R., Pullin, D. I., & Kosović, B. (2001). Direct    |
  !|      numerical simulation of decaying compressible turbulence and |
  !|      shocklet statistics. Physics of Fluids, 13(5), 1415-1430.    |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17 Nov. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine hitgen2d_ic4
    !
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use readwrite, only : readgrid,readinput, readic
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : mpisizedis,parapp,parallelini
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    !
    integer :: i,j
    !
    call readinput
    call readic
    !
    call mpisizedis
    print*, '** mpisizedis done!'
    !
    call parapp
    print*, '** parapp done!'
    !
    call parallelini
    print*, '** parallelini done!'
    !
    call refcal
    print*, '** refcal done!'
    !
    allocate(   x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:im,0:jm,0:0),tmp(0:im,0:jm,0:0),prs(0:im,0:jm,0:0))
    print*, '** allocation finished!'
    !
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    ! call readgrid(trim(gridfile))
    !
    ! generate field
    !call div_free_2d_gen(im,jm,ickmax,vel(0:im,0:jm,0,1), vel(0:im,0:jm,0,2))
    call solenoidal_dilatational_2d_gen(im,jm,ickmax,icsolenoidal,icdilatational,vel(0:im,0:jm,0,1),vel(0:im,0:jm,0,2))
    print*, '** field generated!'
    !
    vel(0:im,0:jm,0,3)=0.d0
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
    enddo
    enddo
    !
    !
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
    !
    call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
    call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
    call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    call h5io_end
    !
    !
    ! test
    print*, '-- Test!'
    !
    call geomcal
    ! 
    call div_test_2d(vel,dvel)
    !
    call hitsta2d
    !
    call mpistop
    !
  end subroutine hitgen2d_ic4
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d.                               |
  !+-------------------------------------------------------------------+
  !
  !
  subroutine hitgen2d_ic4_parallel
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput, readgrid, readic
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    integer :: i,j,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2, wna, var1, var2, ran1, ran2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0,vmax
    complex(8) :: vac1, vac2, crn1, crn2
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c
    real(C_DOUBLE), pointer, dimension(:,:) :: u1r,u2r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u1r, c_u2r
    !
    call readinput
    call readic
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:0),tmp(0:(2*im),0:jm,0:0),prs(0:(2*im),0:jm,0:0))
    if(mpirank==0)  print*, '** allocation finished!'
    !
    ! call readgrid(trim(gridfile))
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_2d(jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do j=1,jm
    do i=1,im
      if(k1(i,j)==0 .and. k2(i,j)==0) then
        u1c(i,j)=0.d0
        u2c(i,j)=0.d0
      else
        call random_number(ran1)
        call random_number(ran2)
        ! ran1: random number distributied in (0,1)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j))
        wn2=real(k2(i,j))
        wna=sqrt(wn1**2+wn2**2)
        !
        var1=IniEnergDis(ISEA*2.d0,ickmax*1.d0,wna)
        var2=sqrt(var1/2.d0/pi/wna)
        !
        vac1=var2*cdexp(crn1)
        vac2=var2*cdexp(crn2)
        !
        u1c(i,j)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(i,j)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do j=1,jm
    do i=1,im
      ! 
      vel(i,j,0,1)=u1r(i,j)
      vel(i,j,0,2)=u2r(i,j)
      !
    end do
    end do
    !
    vel(0,1:jm,0,1)=vel(im,1:jm,0,1)
    vel(0,1:jm,0,2)=vel(im,1:jm,0,2)
    !
    vel(0:im,0:(jm-1),0,1)=vel(0:im,1:jm,0,1)
    vel(0:im,0:(jm-1),0,2)=vel(0:im,1:jm,0,2)
    !
    ! Prepare other field
    vel(0:im,0:jm,0,3)=0.d0
    !
    call dataswap(vel)
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
      vmax = max(vmax,sqrt(vel(i,j,0,1)**2+vel(i,j,0,2)**2))
      !
    enddo
    enddo
    !
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** max velocity    :',vmax
    endif
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    if((j0+jm)==ja)then
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    else
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:(jm-1),0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:(jm-1),0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:(jm-1),0),  dir='k')
    endif
    !
    call h5io_end
    !
    ! Param
    !
    ! if(mpirank == 0) then
    !   ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(ickmax**5)
    !   en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(ickmax**7)
    !   lint=sqrt(2.d0*pi)/ke0
    !   tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(ickmax**7))
    !   eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !   !
    !   print*,' ---------------------------------------------------------------'
    !   print*,'        statistics according to the initial energy spectrum     '
    !   print*,' --------------------------+------------------------------------'
    !   print*,'                   kenergy |',ke0
    !   print*,'                 enstrophy |',en0
    !   print*,'           integral length |',lint
    !   print*,'  large-eddy-turnover time |',tau
    !   print*,'         kolmogorov length |',eta0
    !   print*,' --------------------------+------------------------------------'
    ! endif
    !
    ! Test
    if(mpirank==0)  print*, '-- Test!'
    ! 
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    call geomcal
    !
    call div_test_2d(vel,dvel)
    !
    !call hitsta2d
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call mpistop
    ! 
    !
  end subroutine hitgen2d_ic4_parallel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d_parallel.                               |
  !+-------------------------------------------------------------------+
  !
  !
  subroutine hitgen2d_pic
    ! This code is directly the parallel version
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput, readgrid, readic
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    integer :: i,j,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2, wna, var1, var2, ran1, ran2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0,vmax
    complex(8) :: vac1, vac2, crn1, crn2
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c
    real(C_DOUBLE), pointer, dimension(:,:) :: u1r,u2r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u1r, c_u2r
    !
    call readinput
    call readic
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:0),tmp(0:(2*im),0:jm,0:0),prs(0:(2*im),0:jm,0:0))
    if(mpirank==0)  print*, '** allocation finished!'
    !
    ! call readgrid(trim(gridfile))
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_2d(jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do j=1,jm
    do i=1,im
      if(k1(i,j)==0 .and. k2(i,j)==0) then
        u1c(i,j)=0.d0
        u2c(i,j)=0.d0
      else
        call random_number(ran1)
        call random_number(ran2)
        ! ran1: random number distributied in (0,1)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j))
        wn2=real(k2(i,j))
        wna=sqrt(wn1**2+wn2**2)
        !
        if( wna> ickmax-0.5d0 .and. wna < ickmax+0.5d0) then
          var1=ISEA*2.d0
        else
          var1 = 0.d0
        endif
        var2=sqrt(var1/2.d0/pi/wna)
        !
        vac1=var2*cdexp(crn1)
        vac2=var2*cdexp(crn2)
        !
        u1c(i,j)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(i,j)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do j=1,jm
    do i=1,im
      ! 
      vel(i,j,0,1)=u1r(i,j)
      vel(i,j,0,2)=u2r(i,j)
      !
    end do
    end do
    !
    vel(0,1:jm,0,1)=vel(im,1:jm,0,1)
    vel(0,1:jm,0,2)=vel(im,1:jm,0,2)
    !
    vel(0:im,0:(jm-1),0,1)=vel(0:im,1:jm,0,1)
    vel(0:im,0:(jm-1),0,2)=vel(0:im,1:jm,0,2)
    !
    ! Prepare other field
    vel(0:im,0:jm,0,3)=0.d0
    !
    call dataswap(vel)
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
      vmax = max(vmax,sqrt(vel(i,j,0,1)**2+vel(i,j,0,2)**2))
      !
    enddo
    enddo
    !
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** max velocity    :',vmax
    endif
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    if((j0+jm)==ja)then
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    else
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:(jm-1),0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:(jm-1),0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:(jm-1),0),  dir='k')
    endif
    !
    call h5io_end
    !
    ! Param
    !
    ! if(mpirank == 0) then
    !   ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(ickmax**5)
    !   en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(ickmax**7)
    !   lint=sqrt(2.d0*pi)/ke0
    !   tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(ickmax**7))
    !   eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !   !
    !   print*,' ---------------------------------------------------------------'
    !   print*,'        statistics according to the initial energy spectrum     '
    !   print*,' --------------------------+------------------------------------'
    !   print*,'                   kenergy |',ke0
    !   print*,'                 enstrophy |',en0
    !   print*,'           integral length |',lint
    !   print*,'  large-eddy-turnover time |',tau
    !   print*,'         kolmogorov length |',eta0
    !   print*,' --------------------------+------------------------------------'
    ! endif
    !
    ! Test
    if(mpirank==0)  print*, '-- Test!'
    ! 
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    call geomcal
    !
    call div_test_2d(vel,dvel)
    !
    !call hitsta2d
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call mpistop
    ! 
    !
  end subroutine hitgen2d_pic
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d_pic.                               |
  !+-------------------------------------------------------------------+
  !
  subroutine hitgen2d_herring
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput, readgrid, readic
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          tinf,roinf,spcinf,num_species,nondimen,&
                          ickmax,iomode,icamplitude,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : x,vel,rho,tmp,prs,dvel,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpistop,mpi_sizeof,mpirank, psum, pmax, dataswap
    use geom,      only : geomcal
    use fludyna,   only : thermal
    use hdf5io
    use gridgeneration
    use tecio
    include 'fftw3-mpi.f03'
    !
    integer :: i,j,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2, wna, var1, var2, ran1, ran2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0,vmax
    complex(8) :: vac1, vac2, crn1, crn2
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c
    real(C_DOUBLE), pointer, dimension(:,:) :: u1r,u2r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u1r, c_u2r
    !
    call readinput
    call readic
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:0),tmp(0:(2*im),0:jm,0:0),prs(0:(2*im),0:jm,0:0))
    if(mpirank==0)  print*, '** allocation finished!'
    !
    ! call readgrid(trim(gridfile))
    !
    ! Generate field
    ISEA=1.d0/224.7699d0
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    do j=1,jm
    do i=1,im
      !
      ! if(jm .ne. ja)then
      !   stop "error! jm /= ja"
      ! endif
      !
      if(i <= (ia/2+1)) then
        k1(i,j) = real(i-1,8)
      else if(i<=(ia)) then
        k1(i,j) = real(i-ia-1,8)
      else
        print *,"Error, no wave number possible, i must smaller than ia-1 !"
      end if
      !
      if((j+j0) <= (ja/2+1)) then
        k2(i,j) = real(j+j0-1,8)
      else if((j+j0)<=(ja)) then
        k2(i,j) = real(j+j0-ja-1,8)
      else
        print *,"Error, no wave number possible, (j+jm) must smaller than ja-1 !"
      end if
      !
    end do
    end do
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_2d(jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do j=1,jm
    do i=1,im
      if(k1(i,j)==0 .and. k2(i,j)==0) then
        u1c(i,j)=0.d0
        u2c(i,j)=0.d0
      else
        call random_number(ran1)
        call random_number(ran2)
        ! ran1: random number distributied in (0,1)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=real(k1(i,j))
        wn2=real(k2(i,j))
        wna=sqrt(wn1**2+wn2**2)
        !
        var1=ISEA*wna/(1.d0*ickmax)*exp(-wna/(1.d0*ickmax))
        var2=sqrt(var1/2.d0/pi/wna)
        !
        vac1=var2*cdexp(crn1)
        vac2=var2*cdexp(crn2)
        !
        u1c(i,j)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(i,j)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, '** field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    !
    if(mpirank==0)  print*,' ** project to physical space. '
    !
    im = im*2-2
    !
    do j=1,jm
    do i=1,im
      ! 
      vel(i,j,0,1)=u1r(i,j)
      vel(i,j,0,2)=u2r(i,j)
      !
    end do
    end do
    !
    vel(0,1:jm,0,1)=vel(im,1:jm,0,1)
    vel(0,1:jm,0,2)=vel(im,1:jm,0,2)
    !
    vel(0:im,0:(jm-1),0,1)=vel(0:im,1:jm,0,1)
    vel(0:im,0:(jm-1),0,2)=vel(0:im,1:jm,0,2)
    !
    ! Prepare other field
    vel(0:im,0:jm,0,3)=0.d0
    !
    call dataswap(vel)
    !
    do j=0,jm
    do i=0,im
      !
      rho(i,j,0)  = roinf
      tmp(i,j,0)  = tinf 
      if(nondimen) then
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0))
      else
          spc(i,j,0,:)= spcinf(:)
          prs(i,j,0)  = thermal(density=rho(i,j,0),temperature=tmp(i,j,0),species=spc(i,j,0,:))
      endif
      vel(i,j,0,1)= icamplitude*vel(i,j,0,1)
      vel(i,j,0,2)= icamplitude*vel(i,j,0,2)
      !
      vmax = max(vmax,sqrt(vel(i,j,0,1)**2+vel(i,j,0,2)**2))
      !
    enddo
    enddo
    !
    vmax= pmax(vmax)
    !
    if(mpirank==0)  then
      print*, '** max velocity    :',vmax
    endif
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    if((j0+jm)==ja)then
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    else
      call h5wa2d_r8(varname='ro',var=rho(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='u1',var=vel(0:im,0:(jm-1),0,1),dir='k')
      call h5wa2d_r8(varname='u2',var=vel(0:im,0:(jm-1),0,2),dir='k')
      call h5wa2d_r8(varname='p', var=prs(0:im,0:(jm-1),0),  dir='k')
      call h5wa2d_r8(varname='t', var=tmp(0:im,0:(jm-1),0),  dir='k')
    endif
    !
    call h5io_end
    !
    ! Param
    !
    ! if(mpirank == 0) then
    !   ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(ickmax**5)
    !   en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(ickmax**7)
    !   lint=sqrt(2.d0*pi)/ke0
    !   tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(ickmax**7))
    !   eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !   !
    !   print*,' ---------------------------------------------------------------'
    !   print*,'        statistics according to the initial energy spectrum     '
    !   print*,' --------------------------+------------------------------------'
    !   print*,'                   kenergy |',ke0
    !   print*,'                 enstrophy |',en0
    !   print*,'           integral length |',lint
    !   print*,'  large-eddy-turnover time |',tau
    !   print*,'         kolmogorov length |',eta0
    !   print*,' --------------------------+------------------------------------'
    ! endif
    !
    ! Test
    if(mpirank==0)  print*, '-- Test!'
    ! 
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    call gridcube(2.d0*pi,2.d0*pi,0.d0)
    call geomcal
    !
    call div_test_2d(vel,dvel)
    !
    !call hitsta2d
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call mpistop
    ! 
    !
  end subroutine hitgen2d_herring
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  subroutine hitsta
    !
    use constdef
    use commvar,   only : reynolds,mach,im,jm,km
    use commarray, only : vel,dvel,rho,tmp,x
    use fludyna,   only : miucal
    !
    real(8) :: roav,kenergy,enstrophy,dissp,miuav,urms,kolmlength,ukolm, &
               timekolm,taylorlength,retaylor,tav,machrms
    !
    integer :: i,j,k
    real(8) :: var1,vort1,vort2,vort3,vorts,s11,s22,s33,s12,s13,s23, &
               div,dudx2,miu,ufmx
    !
      roav=0.d0
      tav=0.d0
      !
      kenergy=0.d0
      !
      enstrophy=0.d0
      !
      dissp=0.d0
      !
      miuav=0.d0
      urms=0.d0
      dudx2=0.d0
      !
      ufmx=0.d0
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        roav=roav+rho(i,j,k)
        tav=tav+tmp(i,j,k)
        !
        var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
        kenergy=kenergy+rho(i,j,k)*var1
        !
        urms=urms+var1
        !
        vort1=dvel(i,j,k,2,3)-dvel(i,j,k,3,2)
        vort2=dvel(i,j,k,3,1)-dvel(i,j,k,1,3)
        vort3=dvel(i,j,k,1,2)-dvel(i,j,k,2,1)
        vorts=vort1*vort1+vort2*vort2+vort3*vort3
        !
        enstrophy=enstrophy+rho(i,j,k)*vorts
        !
        miu=miucal(tmp(i,j,k))/reynolds
        miuav=miuav+miu
        !
        s11=dvel(i,j,k,1,1)
        s22=dvel(i,j,k,2,2)
        s33=dvel(i,j,k,3,3)
        div=s11+s22+s33
        s12=0.5d0*(dvel(i,j,k,2,1)+dvel(i,j,k,1,2))
        s13=0.5d0*(dvel(i,j,k,3,1)+dvel(i,j,k,1,3))
        s23=0.5d0*(dvel(i,j,k,3,2)+dvel(i,j,k,2,3))
        dissp=dissp+2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)-num1d3*div**2)
        
        !
        dudx2=dudx2+dvel(i,j,k,1,1)**2+dvel(i,j,k,2,2)**2+dvel(i,j,k,3,3)**2
        !
        ufmx=max(ufmx,abs(vel(i,j,k,1)),abs(vel(i,j,k,2)),abs(vel(i,j,k,3)))
      end do
      end do
      end do
      !
      var1=dble(im*jm*km)

      roav      =roav     /var1   
      tav       =tav     /var1     
      kenergy   =kenergy  /var1*0.5d0
      enstrophy =enstrophy/var1*0.5d0
      dissp     =dissp    /var1
      !
      miuav     =miuav    /var1
      !
      urms      =urms /var1
      dudx2     =dudx2/var1
      !
      kolmlength=sqrt(sqrt((miuav/roav)**3/dissp))
      !
      ukolm=sqrt(sqrt(dissp*miuav/roav))
      !
      timekolm=sqrt(miuav/roav/dissp)
      !
      taylorlength=sqrt(urms/dudx2)
      retaylor=taylorlength*roav*urms/miuav/1.7320508075688773d0
      !
      machrms=urms/sqrt(tav)*mach
      !
      print*,' ---------------------------------------------------------------'
      print*,'              statistics according to actual field              '
      print*,' --------------------------+------------------------------------'
      print*,'                      urms |',urms
      print*,'                   machrms |',machrms
      print*,'                   kenergy |',kenergy
      print*,'           max fluctuation |',ufmx
      print*,'                 enstrophy |',Enstrophy
      print*,'             Kolmlength, η |',kolmlength
      print*,'                       η/Δ |',kolmlength/(x(1,0,0,1)-x(0,0,0,1))
      print*,'                     ukolm |',ukolm
      print*,'                     tkolm |',kolmlength/ukolm
      print*,'              Taylorlength |',taylorlength
      print*,'                  Retaylor |',retaylor
      print*,' --------------------------+------------------------------------'
      !
      !
  end subroutine hitsta
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitsta.                                 |
  !+-------------------------------------------------------------------+
  !
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  subroutine hitsta2d
    !
    use constdef
    use commvar,   only : reynolds,mach,im,jm,ia,ja
    use commarray, only : vel,dvel,rho,tmp,x
    use fludyna,   only : miucal
    use parallel,  only : psum
    !
    real(8) :: roav,kenergy,enstrophy,dissp,miuav,urms,kolmlength,ukolm, &
               timekolm,taylorlength,retaylor,tav,machrms
    !
    integer :: i,j,k
    real(8) :: var1,vort3,vorts,s11,s22,s12,div,dudx2,miu,ufmx
    !
      roav=0.d0
      tav=0.d0
      !
      kenergy=0.d0
      !
      enstrophy=0.d0
      !
      dissp=0.d0
      !
      miuav=0.d0
      urms=0.d0
      dudx2=0.d0
      !
      ufmx=0.d0
      !
      do j=1,jm
      do i=1,im
        !
        roav=roav+rho(i,j,0)
        tav=tav+tmp(i,j,0)
        !
        var1=vel(i,j,0,1)**2+vel(i,j,0,2)**2
        kenergy=kenergy+rho(i,j,0)*var1
        !
        urms=urms+var1
        !
        vort3=dvel(i,j,0,1,2)-dvel(i,j,0,2,1)
        vorts=vort3*vort3
        !
        enstrophy=enstrophy+rho(i,j,0)*vorts
        !
        miu=miucal(tmp(i,j,0))/reynolds
        miuav=miuav+miu
        !
        s11=dvel(i,j,0,1,1)
        s22=dvel(i,j,0,2,2)
        div=s11+s22
        s12=0.5d0*(dvel(i,j,0,2,1)+dvel(i,j,0,1,2))
        dissp=dissp+2.d0*miu*(s11**2+s22**2+2.d0*s12**2-num1d3*div**2)
        
        !
        dudx2=dudx2+dvel(i,j,0,1,1)**2+dvel(i,j,0,2,2)**2+dvel(i,j,0,3,3)**2
        !
        ufmx=max(ufmx,abs(vel(i,j,0,1)),abs(vel(i,j,0,2)))
      end do
      end do
      !
      var1=dble(ia*ja)
      !
      roav      =psum(roav)     /var1   
      tav       =psum(tav)     /var1     
      kenergy   =psum(kenergy)  /var1*0.5d0
      enstrophy =psum(enstrophy)/var1*0.5d0
      dissp     =psum(dissp)   /var1
      !
      miuav     =psum(miuav)    /var1
      !
      urms      =sqrt(psum(urms) /var1)
      dudx2     =psum(dudx2)/var1
      !
      kolmlength=sqrt(sqrt((miuav/roav)**3/dissp))
      !
      ukolm=sqrt(sqrt(dissp*miuav/roav))
      !
      timekolm=sqrt(miuav/roav/dissp)
      !
      taylorlength=urms/sqrt(dudx2)
      retaylor=taylorlength*roav*urms/miuav/1.7320508075688773d0
      !
      machrms=urms/sqrt(tav)*mach
      !
      print*,' ---------------------------------------------------------------'
      print*,'              statistics according to actual field              '
      print*,' --------------------------+------------------------------------'
      print*,'                      urms |',urms
      print*,'                   machrms |',machrms
      print*,'                   kenergy |',kenergy
      print*,'           max fluctuation |',ufmx
      print*,'                 enstrophy |',Enstrophy
      print*,'             Kolmlength, η |',kolmlength
      print*,'                       η/Δ |',kolmlength/(x(1,0,0,1)-x(0,0,0,1))
      print*,'                     ukolm |',ukolm
      print*,'                     tkolm |',kolmlength/ukolm
      print*,'              Taylorlength |',taylorlength
      print*,'                  Retaylor |',retaylor
      print*,' --------------------------+------------------------------------'
      !
      !
  end subroutine hitsta2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitsta2d.                               |
  !+-------------------------------------------------------------------+
  !
  !
    !| This subroutine is used calcualted the divergence.                |
  !+-------------------------------------------------------------------+
  subroutine div_test(u,du)
    !
    use commvar,   only : im,jm,km,hm,ia,ja,ka
    use parallel,  only : dataswap, psum, mpirank
    use comsolver, only : solvrinit,grad
    !
    real(8),intent(inout) :: u(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3)
    !
    real(8),allocatable,intent(out),dimension(:,:,:,:,:) :: du
    !
    integer :: i,j,k
    real(8) :: div,div2
    !
    allocate( du(0:im,0:jm,0:km,1:3,1:3))
    !
    call dataswap(u)
    !
    call solvrinit

    du(:,:,:,:,1)=grad(u(:,:,:,1))
    du(:,:,:,:,2)=grad(u(:,:,:,2))
    du(:,:,:,:,3)=grad(u(:,:,:,3))
    !
    div=0.d0
    div2=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      div =div +du(i,j,k,1,1)+du(i,j,k,2,2)+du(i,j,k,3,3)
      div2=div2+(du(i,j,k,1,1)+du(i,j,k,2,2)+du(i,j,k,3,3))**2
    enddo
    enddo
    enddo
    !
    div = psum(div) /dble(ia*ja*ka)
    div2= psum(div2)/dble(ia*ja*ka)
    !
    if(mpirank == 0) then
      print*,' ** averaged div is:',div
      print*,' ** variance div is:',div2
    endif
    !
  end subroutine div_test
  !+-------------------------------------------------------------------+
  !| The end of the subroutine div_test.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the divergence.                |
  !+-------------------------------------------------------------------+
  subroutine div_test_2d(u,du)
    !
    use commvar,   only : im,jm,km,hm,ia,ja,ka
    use parallel,  only : dataswap, psum, mpirank
    use comsolver, only : solvrinit,grad
    !
    real(8),intent(inout) :: u(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3)
    !
    real(8),allocatable,intent(out),dimension(:,:,:,:,:) :: du
    !
    integer :: i,j,k
    real(8) :: div,div2
    !
    allocate( du(0:im,0:jm,0:0,1:3,1:3))
    !
    call dataswap(u)
    !
    call solvrinit
    !
    du(:,:,:,:,1)=grad(u(:,:,:,1))
    du(:,:,:,:,2)=grad(u(:,:,:,2))
    du(:,:,:,:,3)=0.0d0
    !
    div=0.d0
    div2=0.d0
    do j=1,jm
    do i=1,im
      div =div +du(i,j,0,1,1)+du(i,j,0,2,2)
      div2=div2+(du(i,j,0,1,1)+du(i,j,0,2,2))**2
    enddo
    enddo
    !
    div = psum(div) /dble(ia*ja)
    div2= psum(div2)/dble(ia*ja)
    !
    if(mpirank == 0) then
      print*,' ** averaged div is:',div
      print*,' ** variance div is:',div2
    endif
    !
  end subroutine div_test_2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine div_test_2d.                            |
  !+-------------------------------------------------------------------+
  !
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a divergence free fluctuation.|
  !+-------------------------------------------------------------------+
  !| Ref: Blaisdell, G. A., Numerical simulation of compressible       |
  !|      homogeneous turbulence, Phd, 1991, Stanford University       |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 26-09-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine div_free_gen(idim,jdim,kdim,u1,u2,u3)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,kdim
    real(8),intent(out),dimension(0:idim,0:jdim,0:kdim) :: u1,u2,u3
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wn3,wn12,wna
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wn3: modul of wavenumber in k direction
    ! wn12: wn12=sqrt(wn1**2+wn2**2)
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8), allocatable, dimension(:,:,:) :: u1tp,u2tp,u3tp
    !
    complex(8), allocatable, dimension(:,:,:) :: u1c,u2c,u3c,u1ct,u2ct,u3ct,u4ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,k3,k0,i,j,k
    real(8) :: ran1,ran2,ran3,rn1,rn2,rn3,var1,var2,var3,ISEA
    complex(8) :: vac1,vac2,vac3,vac4,crn1,crn2
    real(8) :: dudi,lambda,ke0,en0,lint,tau,eta0
    !
    kmi=idim/2
    kmj=jdim/2
    kmk=kdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2+kmk**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u2c(-kmi:kmi,-kmj:kmj,-kmk:kmk),                        &
             u3c(-kmi:kmi,-kmj:kmj,-kmk:kmk)                         )
    allocate(u1ct(1:idim,1:jdim,1:kdim),u2ct(1:idim,1:jdim,1:kdim),              &  
             u3ct(1:idim,1:jdim,1:kdim) )
    allocate(u1tp(1:idim,1:jdim,1:kdim),u2tp(1:idim,1:jdim,1:kdim),              &
             u3tp(1:idim,1:jdim,1:kdim),u4ct(1:idim,1:jdim,1:kdim)               )
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    k0=4
    !
    ! Generate the random velocity field according the given energy 
    ! spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefor, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    do k1=0,kmi
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      !
      call random_number(ran1)
      call random_number(ran2)
      call random_number(ran3)
      ! ran1,ran2,ran3: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      rn2=ran2*2.d0*pi
      rn3=ran3*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wn3=real(k3,8)
      wn12=sqrt(wn1**2+wn2**2)
      wna=sqrt(wn1**2+wn2**2+wn3**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0 .and. k3==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      crn1=rn1*(0.d0,1.d0)
      crn2=rn2*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)*dcos(rn3)
      vac2=var2*cdexp(crn2)*dsin(rn3)
      !
      if(k1==0 .and. k2==0 .and. k3==0) then
        u1c(k1,k2,k3)=0.d0
        u2c(k1,k2,k3)=0.d0
        u3c(k1,k2,k3)=0.d0
      elseif(k1==0 .and. k2==0) then
        u1c(k1,k2,k3)=vac1
        u2c(k1,k2,k3)=vac2
        u3c(k1,k2,k3)=0.d0
      else
        u1c(k1,k2,k3)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)
        u2c(k1,k2,k3)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)
        u3c(k1,k2,k3)=-vac2*wn12/wna
      end if
      !
    end do
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
    do k3=-kmk,kmk
      u1c(k1,k2,k3)=conjg(u1c(-k1,-k2,-k3))
      u2c(k1,k2,k3)=conjg(u2c(-k1,-k2,-k3))
      u3c(k1,k2,k3)=conjg(u3c(-k1,-k2,-k3))
    end do
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) fo rthe 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do k=1,kdim
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      if(k<=kdim/2+1) then
        k3=k-1
      else
        k3=k-kdim-1
      end if
      !
      u1ct(i,j,k)=u1c(k1,k2,k3)
      u2ct(i,j,k)=u2c(k1,k2,k3)
      u3ct(i,j,k)=u3c(k1,k2,k3)
    end do
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    u3ct=FFT(u3ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do k=1,kdim
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNyNz) for return standard FFT
      u1(i,j,k)=real(u1ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      u2(i,j,k)=real(u2ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      u3(i,j,k)=real(u3ct(i,j,k),8)*sqrt(real(idim*jdim*kdim,8))
      !
    end do
    end do
    end do
    !
    u1(0,1:jdim,1:kdim)=u1(idim,1:jdim,1:kdim)
    u2(0,1:jdim,1:kdim)=u2(idim,1:jdim,1:kdim)
    u3(0,1:jdim,1:kdim)=u3(idim,1:jdim,1:kdim)
    !
    u1(0:idim,0,1:kdim)=u1(0:idim,jdim,1:kdim)
    u2(0:idim,0,1:kdim)=u2(0:idim,jdim,1:kdim)
    u3(0:idim,0,1:kdim)=u3(0:idim,jdim,1:kdim)
    !
    u1(0:idim,0:jdim,0)=u1(0:idim,0:jdim,kdim)
    u2(0:idim,0:jdim,0)=u2(0:idim,0:jdim,kdim)
    u3(0:idim,0:jdim,0)=u3(0:idim,0:jdim,kdim)
    !
    ! urms=0.d0
    ! ! ufmx=0.d0
    ! ! Kenergy=0.d0
    ! do k=1,kdim
    ! do j=1,jdim
    ! do i=1,idim
    !   Kenergy=Kenergy+0.5d0*(u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2)
    !   urms=urms+u1(i,j,k)**2+u2(i,j,k)**2+u3(i,j,k)**2
    !   ufmx=max(ufmx,dabs(u1(i,j,k)),dabs(u2(i,j,k)),dabs(u3(i,j,k)))
    ! end do
    ! end do
    ! end do
    ! urms=sqrt(urms/real(idim*jdim*kdim,8))
    ! Kenergy=Kenergy/real(idim*jdim*kdim,8)
    ! !
    ! u1=u1/urms
    ! u2=u2/urms
    ! u3=u3/urms
    ! Kenergy=Kenergy/urms/urms
    ! urms=urms/urms
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    ! call h5srite(var=u1,varname='u1',filename='velocity.h5',explicit=.true.,newfile=.true.)
    ! call h5srite(var=u2,varname='u2',filename='velocity.h5',explicit=.true.)
    ! call h5srite(var=u3,varname='u3',filename='velocity.h5',explicit=.true.)
    !
    deallocate(u1c,u2c,u3c)
    deallocate(u1ct,u2ct,u3ct)
    deallocate(u1tp,u2tp,u3tp)
    !
  end subroutine div_free_gen
  !+-------------------------------------------------------------------+
  !| The end of the function div_free_gen.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a divergence free fluctuation |
  !| in 2D.                                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16 Oct. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine div_free_2d_gen(idim,jdim,k0,u1,u2)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,k0
    real(8),intent(out),dimension(0:idim,0:jdim,0:0) :: u1,u2
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wna, Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    complex(8), allocatable, dimension(:,:) :: u1c,u2c,u1ct,u2ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,i,j,n,clock, irandom
    !
    real(8) :: ran1,rn1,var1,var2
    complex(8) :: crn1,vac1
    real(8) :: ke0,en0,lint,tau,eta0
    integer,allocatable :: seed(:)
    !
    kmi=idim/2
    kmj=jdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj),u2c(-kmi:kmi,-kmj:kmj))
    allocate(u1ct(1:idim,1:jdim),u2ct(1:idim,1:jdim))
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    !
    ! Generate the random velocity field according the given energy spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefore, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    !
    ! Insert random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    do k1=0,kmi
    do k2=-kmj,kmj
      !
      !
      call random_number(ran1)
      ! ran1: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wna=sqrt(wn1**2+wn2**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      !
      crn1=rn1*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)
      !
      if(k1==0 .and. k2==0) then
        u1c(k1,k2)=0.d0
        u2c(k1,k2)=0.d0
      else
        u1c(k1,k2)=vac1*wn2/wna
        u2c(k1,k2)=-vac1*wn1/wna
      end if
      !
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
      u1c(k1,k2)=conjg(u1c(-k1,-k2))
      u2c(k1,k2)=conjg(u2c(-k1,-k2))
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) for the 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      !
      u1ct(i,j)=u1c(k1,k2)
      u2ct(i,j)=u2c(k1,k2)
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNy) for return standard FFT
      u1(i,j,0)=real(u1ct(i,j),8)*sqrt(real(idim*jdim,8))
      u2(i,j,0)=real(u2ct(i,j),8)*sqrt(real(idim*jdim,8))
      !
    end do
    end do
    !
    u1(0,1:jdim,0)=u1(idim,1:jdim,0)
    u2(0,1:jdim,0)=u2(idim,1:jdim,0)
    !
    u1(0:idim,0,0)=u1(0:idim,jdim,0)
    u2(0:idim,0,0)=u2(0:idim,jdim,0)
    !
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    !
    deallocate(u1c,u2c)
    deallocate(u1ct,u2ct)
    !
  end subroutine div_free_2d_gen
  !+-------------------------------------------------------------------+
  !| The end of the function div_free_2d_gen.                             |
  !+-------------------------------------------------------------------+
  !!
    !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a solenoidal-dilatational     !
  !| fluctuation in 2D.                                                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17 Nov. 2023: Created by C.S. Luo @ Beihang University            |
  !+-------------------------------------------------------------------+
  subroutine solenoidal_dilatational_2d_gen(idim,jdim,k0,icsolenoidal,icdilatational,u1,u2)
    !
    use singleton
    use commvar,only : Reynolds
    !
    integer,intent(in) :: idim,jdim,k0
    real(8),intent(in) :: icsolenoidal,icdilatational
    real(8),intent(out),dimension(0:idim,0:jdim,0:0) :: u1,u2
    !
    ! local data
    integer :: kmi,kmj,kmk,kmax
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! kmi: maximal wavenumber in i direction
    ! kmj: maximal wavenumber in j direction
    ! kmk: maximal wavenumber in k direction
    ! kmax: maximal wavenumber in all direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: wn1,wn2,wna, Ac
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! wn1: modul of wavenumber in i direction
    ! wn2: modul of wavenumber in j direction
    ! wna: modul of wavenumber in all direction
    ! (k0*1.d0): the wavenumber at maximum given 
    !     spectrum
    ! Ac: the intensity of given spectrum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    complex(8), allocatable, dimension(:,:) :: u1c,u2c,u1ct,u2ct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Egv: the given initial energy spectrum
    ! u1c: the spectral velocity in k1 direction
    ! u2c: the spectral velocity in k2 direction
    ! u3c: the spectral velocity in k3 direction
    ! uct: the spectrl variable in (1~*2km)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(8) :: Kenergy,Enstropy,ITGscale,LETT,KolmLength,urms,ufmx,ISEA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Kenergy: initial volume averaged turbulent 
    !          kinetic energy
    ! Enstropy: initial volume averaged e
    !           nstrophy
    ! ITGscale: initial integral length scale
    ! LETT: initial large-eddy-turnover time
    ! KolmLength: initial Kolmogorov scale
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: k1,k2,i,j,n,clock, irandom
    !
    real(8) :: ran1,ran2,rn1,rn2,var1,var2
    complex(8) :: crn1,crn2,vac1,vac2
    real(8) :: ke0,en0,lint,tau,eta0
    integer,allocatable :: seed(:)
    !
    kmi=idim/2
    kmj=jdim/2
    kmax=idnint(sqrt((kmi**2+kmj**2)*1.d0))+1
    !
    allocate(u1c(-kmi:kmi,-kmj:kmj),u2c(-kmi:kmi,-kmj:kmj))
    allocate(u1ct(1:idim,1:jdim),u2ct(1:idim,1:jdim))
    !
    ! Give the inital energy spectrum.
    ISEA=1.d0/224.7699d0
    !
    ! Generate the random velocity field according the given energy spectrum
    ! Blaisdell, G. A. 1991 took E(k)=Integer(ui*uicoj*dA(k)). 
    ! This program takes E(k)=Integer(0.5*ui*uicoj*dA(k)). 
    ! Therefore, we take the Ek as twice of that from Blaisdell.
    print*,' ** Generate the random velocity field according the given energy spectrum'
    !
    ! Insert random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    do k1=0,kmi
    do k2=-kmj,kmj
      !
      !
      call random_number(ran1)
      call random_number(ran2)
      ! ran1: random number distributied in (0,1)
      !
      rn1=ran1*2.d0*pi
      rn2=ran2*2.d0*pi
      !
      ! Calculate the modul of the wavenumber in each direction
      wn1=real(k1,8)
      wn2=real(k2,8)
      wna=sqrt(wn1**2+wn2**2)
      !
      ! Calculate the initidiml energy spectral
      if(k1==0 .and. k2==0) then
        var1=0.d0
        var2=0.d0
      else
        var1=IniEnergDis(ISEA*2.d0,K0*1.d0,wna)
        var2=sqrt(var1/4.d0/pi/wna**2)
        ! var2=1.d0
      end if
      !
      ! Gererate the velocity spectrum in half-wavenumber space.
      !
      crn1=rn1*(0.d0,1.d0)
      crn2=rn2*(0.d0,1.d0)
      !
      vac1=var2*cdexp(crn1)
      vac2=var2*cdexp(crn2)
      !
      if(k1==0 .and. k2==0) then
        u1c(k1,k2)=0.d0
        u2c(k1,k2)=0.d0
      else
        u1c(k1,k2)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(k1,k2)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
      !
    end do
    end do
    !
    print*,' ** Generate the velocity spectrum in another half-wavenumber space '
    ! Generate the velocity spectrum in another half-wavenumber space
    ! by using conjunction relation
    do k1=-kmi,-1
    do k2=-kmj,kmj
      u1c(k1,k2)=conjg(u1c(-k1,-k2))
      u2c(k1,k2)=conjg(u2c(-k1,-k2))
    end do
    end do
    ! !
    ! Transform the spectrum from (-N/2+1,N/2) to (1,N) for the 
    ! convenience of using external FFT subroutine
    print*,' ** Transform the spectrum from (-N/2+1,N/2) to (1,N)  '
    !
    do j=1,jdim
    do i=1,idim
      if(i<=idim/2+1) then
        k1=i-1
      else
        k1=i-idim-1
      end if
      if(j<=jdim/2+1) then
        k2=j-1
      else
        k2=j-jdim-1
      end if
      !
      u1ct(i,j)=u1c(k1,k2)
      u2ct(i,j)=u2c(k1,k2)
    end do
    end do
    ! !
    u1ct=FFT(u1ct,inv=.true.)
    u2ct=FFT(u2ct,inv=.true.)
    !
    print*,' ** project to physical space. '
    !
    do j=1,jdim
    do i=1,idim
      ! multiply sqrt(NxNy) for return standard FFT
      u1(i,j,0)=real(u1ct(i,j),8)*sqrt(real(idim*jdim,8))
      u2(i,j,0)=real(u2ct(i,j),8)*sqrt(real(idim*jdim,8))
      !
    end do
    end do
    !
    u1(0,1:jdim,0)=u1(idim,1:jdim,0)
    u2(0,1:jdim,0)=u2(idim,1:jdim,0)
    !
    u1(0:idim,0,0)=u1(0:idim,jdim,0)
    u2(0:idim,0,0)=u2(0:idim,jdim,0)
    !
    ! !
    ke0=3.d0*ISEA/64.d0*sqrt(2.d0*pi)*dble(k0**5)
    en0=15.d0*ISEA/256.d0*sqrt(2.d0*pi)*dble(k0**7)
    lint=sqrt(2.d0*pi)/ke0
    tau =sqrt(32.d0/ISEA*sqrt(2.d0*pi))/sqrt(dble(k0**7))
    eta0=1.d0/sqrt(sqrt(2.d0*en0*Reynolds**2))
    !
    print*,' ---------------------------------------------------------------'
    print*,'        statistics according to the initial energy spectrum     '
    print*,' --------------------------+------------------------------------'
    print*,'                   kenergy |',ke0
    print*,'                 enstrophy |',en0
    print*,'           integral length |',lint
    print*,'  large-eddy-turnover time |',tau
    print*,'         kolmogorov length |',eta0
    print*,' --------------------------+------------------------------------'
    ! !
    !
    deallocate(u1c,u2c)
    deallocate(u1ct,u2ct)
    !
  end subroutine solenoidal_dilatational_2d_gen
  !+-------------------------------------------------------------------+
  !| The end of the function solenoidal_dilatational_2d_gen.           |
  !+-------------------------------------------------------------------+
  !!
  !
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calcuate the spectral energy at any 
  ! wavenumber.
  ! Ref: S. JAMME, et al. Direct Numerical Simulation of the 
  ! Interaction between a Shock Wave and Various Types of Isotropic 
  ! Turbulence, Flow, Turbulence and Combustion, 2002, 68:227每268.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function IniEnergDis(Ac,k0,wnb)
    !
    real(8) :: k0,Ac,var1,wnb,IniEnergDis
    !
    var1=-2.d0*(wnb/k0)**2
    IniEnergDis=Ac*wnb**4*exp(var1)
    !IniEnergDis=Ac*wnb**(-5.d0/3.d0)
    !
    return
    !
  end function IniEnergDis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function Ek.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end module udf_pp_hitgen