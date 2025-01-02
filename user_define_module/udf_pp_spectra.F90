!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| spectrual calculations.                                             |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_spectra
    !
    !
    use constdef
    use stlaio,  only: get_unit
    use udf_tool
    !
    implicit none
    !
    contains
    !
    subroutine ppSpectraentrance
        !
        use cmdefne
        use parallel,        only : mpirank,bcast,mpisize,lio
        !
        ! local data
        character(len=64) :: casefolder,inputfile,outputfile,viewmode, &
                             flowfieldfile, readmode, method
        integer :: filenumb, methodnumb
        !
        !
        if(mpirank == 0) then
            call readkeyboad(readmode) 
        endif
        call bcast(readmode)
        !
        if(trim(readmode)=='instant2D') then
          !
          if(mpirank == 0) then
              call readkeyboad(inputfile) 
              read(inputfile,'(i4)') filenumb
              call readkeyboad(method)
              read(method,'(i1)') methodnumb
          endif
          call bcast(filenumb)
          call bcast(methodnumb)
          call instantspectra2D(filenumb,methodnumb)
          !
        elseif(trim(readmode)=='instant3D') then
          !
          if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            call readkeyboad(method)
            read(method,'(i1)') methodnumb
          endif
          call bcast(filenumb)
          call bcast(methodnumb)
          call instantspectra3D(filenumb,methodnumb)
          !
        elseif(trim(readmode)=='skewness2D') then
          !
          if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
          endif
          call bcast(filenumb)
          call instantspectraskewness2D(filenumb)
          !
        elseif(trim(readmode)=='triad2D') then
          !
          if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            call readkeyboad(method)
            read(method,'(i1)') methodnumb
          endif
          call bcast(filenumb)
          call bcast(methodnumb)
          call instanttriad2D(filenumb,methodnumb)
          !
        elseif(trim(readmode)=='triad3D') then
          !
          if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            call readkeyboad(method)
            read(method,'(i1)') methodnumb
          endif
          call bcast(filenumb)
          call bcast(methodnumb)
          call instanttriad3D(filenumb,methodnumb)
          !
        elseif(trim(readmode)=='initparam2D') then
          !
          call initparam2D
          !
        elseif(trim(readmode)=='initparam3D') then
          ! 
          call initparam3D
          !
        else
            print* ,"Readmode is not defined!", readmode
        endif
        !
    end subroutine ppSpectraentrance
    !
  !
  !
  !
  subroutine instantspectra2D(thefilenumb,method)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka,roinf
    use commarray, only: vel, rho, prs
    use hdf5io
    use solver,    only : refcal
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    integer,intent(in) :: method
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe,rou1spe,rou2spe
    real(8), allocatable, dimension(:) :: ES,ED,Ecount,Eweak,Pud,kn,Erho
    complex(8), allocatable, dimension(:,:) :: usspe,udspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1d,u2d,u1s,u2s
    ! real(8), allocatable, dimension(:,:) :: usspeR,usspeI,ucspeR,ucspeI
    real(8), allocatable, dimension(:,:) :: u1dR,u2dR,u1sR,u2sR,udud,udus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer :: allkmax, kOrdinal
    real(8) :: kk,dk,lambda
    real(8) :: Edspe,Esspe,Erhospe,Pudspe,Edphy,Esphy,Erhophy,udusphy,Edmax
    real(8) :: ReUsConjUdmax,ReUsConjUdmin,k2Es,k2Ed,IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b,hand_c
    character(len=1) :: modeio
    integer :: i,j,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, &
     c_u1d, c_u2d, c_u1s, c_u2s, c_rou1spe, c_rou2spe
    !
    ! Initialization
    if((method < 1 ) .or. (method >3))then
      stop "Error! method @ instantspectra2D problem, not among 1, 2, 3"
    else
      if(mpirank==0) print *, "Using method", method
    endif
    !
    call readinput
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    modeio='h'
    !
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    dk = 1.d0
    if(method==1)then
      lambda = 1.21d0
      allkmax = ceiling(log(real(sqrt(2.d0)/3*min(ia,ja))/dk)/log(lambda))
    elseif((method == 2) .or. (method==3))then
      allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
    endif
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"knumber:",allkmax
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "** Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
    do j=1,jm
      u1mean = u1mean + vel(i,j,0,1)
      u2mean = u2mean + vel(i,j,0,2)
      rhomean = rhomean + rho(i,j,0)
      prsmean = prsmean + prs(i,j,0)
    enddo
    enddo
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    !
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    c_rou1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rou1spe, rou1spe, [imfftw,jmfftw])
    c_rou2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rou2spe, rou2spe, [imfftw,jmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1),0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2),0.d0,C_INTPTR_T);
      pspe(i,j)=CMPLX(prs(i,j,0),0.d0,C_INTPTR_T);
      rou1spe(i,j)=CMPLX((rho(i,j,0)-roinf)*vel(i,j,0,1),C_INTPTR_T)
      rou2spe(i,j)=CMPLX((rho(i,j,0)-roinf)*vel(i,j,0,2),C_INTPTR_T)
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    call fftw_mpi_execute_dft(forward_plan,rou1spe,rou1spe)
    call fftw_mpi_execute_dft(forward_plan,rou2spe,rou2spe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      rou1spe(i,j)=rou1spe(i,j)/(1.d0*ia*ja)
      rou2spe(i,j)=rou2spe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
    !
    !!!! Do S-C decomposition
    allocate(usspe(1:im,1:jm),udspe(1:im,1:jm))
    c_u1d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1d, u1d, [imfftw,jmfftw])
    c_u2d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2d, u2d, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    allocate(u1dR(1:im,1:jm),u2dR(1:im,1:jm),u1sR(1:im,1:jm),u2sR(1:im,1:jm),&
            udud(1:im,1:jm),udus(1:im,1:jm),usus(1:im,1:jm))
    allocate(ReUsConjUd(1:im,1:jm))
    !
    ! 
    !
    do j=1,jm
    do i=1,im
        kk=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
        udspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
        !
        u1d(i,j)=  udspe(i,j)*k1(i,j)/kk
        u2d(i,j)=  udspe(i,j)*k2(i,j)/kk
        u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
        !
        ReUsConjUd(i,j) = real(conjg(usspe(i,j))*udspe(i,j))
        !
      end do
    end do
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !
    !!!! Give S-C spectra and spectral energy
    !
    allocate(ES(0:allkmax),ED(0:allkmax),Erho(0:allkmax),Ecount(0:allkmax))
    allocate(Eweak(0:allkmax),Pud(0:allkmax),kn(0:allkmax))
    !
    ES = 0.0d0
    ED = 0.0d0
    Ecount = 0.0d0
    Eweak = 0.0d0
    Erho = 0.0d0
    Pud = 0.0d0
    kn = 0.d0
    Edspe = 0.0d0
    Esspe = 0.0d0
    Pudspe = 0.0d0
    k2Es = 0.d0
    k2Ed = 0.d0
    !
    do j=1,jm
    do i=1,im
        kk=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        kOrdinal = kint(kk,dk,method,lambda)
        if (kOrdinal<=allkmax) then
          Ecount(kOrdinal) = Ecount(kOrdinal) + 1
          if((method == 1) .or. (method == 2))then
            ES(kOrdinal) = ES(kOrdinal) + roinf * usspe(i,j)*conjg(usspe(i,j))*kk/2
            ED(kOrdinal) = ED(kOrdinal) + roinf * udspe(i,j)*conjg(udspe(i,j))*kk/2
            Eweak(kOrdinal) = Eweak(kOrdinal) + roinf * u1spe(i,j)*conjg(u1spe(i,j))*kk/2 + &
                            roinf * u2spe(i,j)*conjg(u2spe(i,j))*kk/2
            Pud(kOrdinal) = Pud(kOrdinal) - dimag(pspe(i,j)*dconjg(udspe(i,j))*kk)*kk
            Erho(kOrdinal) = Erho(kOrdinal) + real(rou1spe(i,j)*conjg(u1spe(i,j))+ &
                              rou2spe(i,j)*conjg(u2spe(i,j)))*kk/2
            kn(kOrdinal) = kn(kOrdinal) + kk
          elseif(method == 3)then
            ES(kOrdinal) = ES(kOrdinal) + roinf * usspe(i,j)*conjg(usspe(i,j))/2
            ED(kOrdinal) = ED(kOrdinal) + roinf * udspe(i,j)*conjg(udspe(i,j))/2
            Eweak(kOrdinal) = Eweak(kOrdinal) + roinf * u1spe(i,j)*conjg(u1spe(i,j))/2 + &
                              roinf * u2spe(i,j)*conjg(u2spe(i,j))/2
            Pud(kOrdinal) = Pud(kOrdinal) - dimag(pspe(i,j)*dconjg(udspe(i,j))*kk)
            Erho(kOrdinal) = Erho(kOrdinal) + real(rou1spe(i,j)*conjg(u1spe(i,j))+ &
                              rou2spe(i,j)*conjg(u2spe(i,j)))/2
          endif
        end if
        Edspe = Edspe + roinf * (udspe(i,j)*dconjg(udspe(i,j)))/2
        Esspe = Esspe + roinf * (usspe(i,j)*dconjg(usspe(i,j)))/2
        Erhospe = Erhospe + real(rou1spe(i,j)*conjg(u1spe(i,j))+ rou2spe(i,j)*conjg(u2spe(i,j)))/2
        IntLengthAbove = IntLengthAbove + roinf * usspe(i,j)*conjg(usspe(i,j))/2/kk + &
                        roinf * (udspe(i,j)*dconjg(udspe(i,j)))/2/kk
        Pudspe = Pudspe + dimag(pspe(i,j)*dconjg(udspe(i,j))*kk)/2
        k2Es = k2Es + kk**2 * (usspe(i,j)*dconjg(usspe(i,j)))/2
        k2Ed = k2Ed + kk**2 * (udspe(i,j)*dconjg(udspe(i,j)))/2
      end do
    end do
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      ED(i) = psum(ED(i))
      Eweak(i) = psum(Eweak(i))
      Erho(i) = psum(Erho(i))
      Pud(i) = psum(Pud(i))
      Ecount(i) = psum(Ecount(i))
      if((method == 1) .or. (method == 2))then
        kn(i) = psum(kn(i))
        if(Ecount(i) .ne. 0)then
          ES(i) = ES(i)/Ecount(i)*2*pi
          ED(i) = ED(i)/Ecount(i)*2*pi
          Eweak(i) = Eweak(i)/Ecount(i)*2*pi
          Pud(i) = Pud(i)/Ecount(i)*2*pi
          kn(i) = kn(i)/Ecount(i)
          Erho(i) = Erho(i)/Ecount(i)*2*pi
        endif
      else
        kn(i) = real(i)
      endif
    end do
    !
    Edspe = psum(Edspe)
    Esspe = psum(Esspe)
    Erhospe = psum(Erhospe)
    Pudspe = psum(Pudspe)
    IntLengthAbove = psum(IntLengthAbove)
    k2Es = psum(k2Es)
    k2Ed = psum(k2Ed)
    !
    if(mpirank==0)  print*, '** Summation & average'
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1d,u1d)
    call fftw_mpi_execute_dft(backward_plan,u2d,u2d)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    !
    !!!! Give S-C physical energy
    Edphy = 0.0d0
    Esphy = 0.0d0
    Erhophy = 0.0d0
    udusphy = 0.0d0
    Edmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    do j=1,jm
    do i=1,im
        u1dR(i,j) = real(u1d(i,j))
        u2dR(i,j) = real(u2d(i,j))
        u1sR(i,j) = real(u1s(i,j))
        u2sR(i,j) = real(u2s(i,j))
        udud(i,j) = u1dR(i,j)*u1dR(i,j)+u2dR(i,j)*u2dR(i,j)
        udus(i,j) = u1dR(i,j)*u1sR(i,j)+u2dR(i,j)*u2sR(i,j)
        usus(i,j) = u1sR(i,j)*u1sR(i,j)+u2sR(i,j)*u2sR(i,j)
        Edmax = max(Edmax,udud(i,j))
        Edphy = Edphy + roinf * udud(i,j) / 2.d0
        Esphy = Esphy + roinf * usus(i,j) / 2.d0
        Erhophy = Erhophy + (rho(i,j,0)-roinf)*(vel(i,j,0,1)**2 + vel(i,j,0,2)**2) / 2.d0
        udusphy = udusphy + roinf * udus(i,j)
        ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j))
        ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j))
      end do
    end do
    !
    Edphy = psum(Edphy)/(1.d0*ia*ja)
    Esphy = psum(Esphy)/(1.d0*ia*ja)
    udusphy = psum(udusphy)/(1.d0*ia*ja)
    Erhophy = psum(Erhophy)/(1.d0*ia*ja)
    Edmax = pmax(Edmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES ED Eweak Erho Pud')
      do i=0,allkmax
        if(Ecount(i)>1e-3) call listwrite(hand_a,kn(i),ES(i),ED(i),Eweak(i),Erho(i),Pud(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time Edphy Esphy udusphy Eweakphy Erhophy Edspe Esspe Pudspe Eweakspe Erhospe Edmax k2Es k2Ed')! IntLen')
      call listwrite(hand_b,Edphy,Esphy,udusphy,Edphy+Esphy+udusphy,Erhophy,&
                    Edspe,Esspe,Pudspe,Edspe+Esspe,Erhospe,Edmax,k2Es,k2Ed)!,IntLengthAbove/(Edspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      ! if (thefilenumb .ne. 0) then
      !   outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      ! else
      !   outfilename = 'pp/Espec_ReUsConjUd.dat'
      ! endif
      ! !
      ! call listinit(filename=outfilename,handle=hand_b, &
      !               firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      ! call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      ! !
      ! print*,' <<< '//outfilename//'... done.'
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_mean'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_mean.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_c, &
                    firstline='nstep time u1mean u2mean rhomean prsmean')
      call listwrite(hand_c,u1mean,u2mean,rhomean,prsmean)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_rou1spe)
    call fftw_free(c_rou2spe)
    call fftw_free(c_u1d)
    call fftw_free(c_u1s)
    call fftw_free(c_u2d)
    call fftw_free(c_u2s)
    call mpistop
    
    deallocate(ES,ED,Ecount,Eweak,Pud,kn,Erho)
    deallocate(usspe,udspe)
    deallocate(u1dR,u2dR,u1sR,u2sR,udud,udus,usus,ReUsConjUd)
    deallocate(k1,k2)
    !
  end subroutine instantspectra2D
  !
  !
  subroutine instantspectra3D(thefilenumb,method)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka,roinf
    use commarray, only: vel, rho, prs
    use hdf5io
    use solver,    only : refcal
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    integer,intent(in) :: method
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe,rou1spe,rou2spe,rou3spe
    real(8), allocatable, dimension(:) :: ES,ED,Ecount,Eweak,Pud,kn,Erho
    complex(8), allocatable, dimension(:,:,:) :: udspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1d,u2d,u3d,u1s,u2s,u3s
    real(8), allocatable, dimension(:,:,:) :: udud,udus,usus,ReUsConjUd
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer :: allkmax, kOrdinal
    real(8) :: kk,dk,lambda
    real(8) :: Edspe,Esspe,Erhospe,Pudspe,Edphy,Esphy,Erhophy,udusphy,Edmax,ReUsConjUdmax,ReUsConjUdmin
    real(8) :: IntLengthAbove
    character(len=128) :: outfilename
    integer :: hand_a,hand_b,hand_c
    character(len=1) :: modeio
    integer :: i,j,k,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe, c_u1d, c_u2d, c_u3d, &
                  c_u1s, c_u2s, c_u3s, c_rou1spe, c_rou2spe, c_rou3spe
    !
    ! Initialization
    if((method < 1 ) .or. (method >3))then
      stop "Error! method @ instantspectra3D problem, not among 1, 2, 3"
    else
      if(mpirank==0) print *, "Using method", method
    endif
    !
    call readinput
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    modeio='h'
    !
    if(ka == 0) stop 'Please use instantspectra2D'
    !
    dk = 1.d0
    if(method==1)then
      lambda = 1.21d0
      allkmax = ceiling(log(real(sqrt(2.d0)/3*min(min(ia,ja),ka))/dk)/log(lambda))
    elseif((method == 2) .or. (method==3))then
      allkmax=ceiling(real(sqrt(2.d0)/3*min(min(ia,ja),ka))/dk)
    endif
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:", ka,"knumber:",allkmax
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    c_rou1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rou1spe, rou1spe, [imfftw,jmfftw,kmfftw])
    c_rou2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rou2spe, rou2spe, [imfftw,jmfftw,kmfftw])
    c_rou3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rou3spe, rou3spe, [imfftw,jmfftw,kmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2)-u2mean,0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3)-u3mean,0.d0,C_INTPTR_T);
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T);
      rou1spe(i,j,k)=CMPLX((rho(i,j,k)-roinf)*vel(i,j,k,1),C_INTPTR_T)
      rou2spe(i,j,k)=CMPLX((rho(i,j,k)-roinf)*vel(i,j,k,2),C_INTPTR_T)
      rou3spe(i,j,k)=CMPLX((rho(i,j,k)-roinf)*vel(i,j,k,3),C_INTPTR_T)
      !
    enddo
    enddo
    enddo
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    call fftw_mpi_execute_dft(forward_plan,rou1spe,rou1spe)
    call fftw_mpi_execute_dft(forward_plan,rou2spe,rou2spe)
    call fftw_mpi_execute_dft(forward_plan,rou3spe,rou3spe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      rou1spe(i,j,k)=rou1spe(i,j,k)/(1.d0*ia*ja*ka)
      rou2spe(i,j,k)=rou2spe(i,j,k)/(1.d0*ia*ja*ka)
      rou3spe(i,j,k)=rou3spe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
    !
    !!!! Do S-C decomposition
    allocate(udspe(1:im,1:jm,1:km))
    c_u1d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1d, u1d, [imfftw,jmfftw,kmfftw])
    c_u2d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2d, u2d, [imfftw,jmfftw,kmfftw])
    c_u3d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3d, u3d, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    allocate(udud(1:im,1:jm,1:km),udus(1:im,1:jm,1:km),usus(1:im,1:jm,1:km))
    allocate(ReUsConjUd(1:im,1:jm,1:km))
    !
    ! 
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15
      !
      udspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1d(i,j,k)=  k1(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u2d(i,j,k)=  k2(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u3d(i,j,k)=  k3(i,j,k)*k1(i,j,k)/kk * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/kk * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/kk * u3spe(i,j,k)
      u1s(i,j,k)=  u1spe(i,j,k) - u1d(i,j,k)
      u2s(i,j,k)=  u2spe(i,j,k) - u2d(i,j,k)
      u3s(i,j,k)=  u3spe(i,j,k) - u3d(i,j,k)
      !
      ReUsConjUd(i,j,k) = real(conjg(u1s(i,j,k))*udspe(i,j,k)) + real(conjg(u2s(i,j,k))*udspe(i,j,k)) + &
                         real(conjg(u3s(i,j,k))*udspe(i,j,k))
      !
    enddo
    enddo
    enddo
    !
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !!!! Give S-C spectra and spectral energy
    !
    allocate(ES(0:allkmax),ED(0:allkmax),Erho(0:allkmax),Ecount(0:allkmax))
    allocate(Eweak(0:allkmax),Pud(0:allkmax),kn(0:allkmax))
    !
    ES = 0.0d0
    ED = 0.0d0
    Ecount = 0.0d0
    Eweak = 0.0d0
    Erho = 0.d0
    Pud = 0.0d0
    kn = 0.d0
    Edspe = 0.0d0
    Esspe = 0.0d0
    Pudspe = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      kOrdinal = kint(kk,dk,method,lambda)
      if (kOrdinal<=allkmax) then
        Ecount(kOrdinal) = Ecount(kOrdinal) + 1
        if((method == 1) .or. (method == 2))then
          ES(kOrdinal) = ES(kOrdinal) + roinf * u1s(i,j,k)*conjg(u1s(i,j,k))*kk*kk/2 + &
                        roinf * u2s(i,j,k)*conjg(u2s(i,j,k))*kk*kk/2 + &
                        roinf * u3s(i,j,k)*conjg(u3s(i,j,k))*kk*kk/2
          ED(kOrdinal) = ED(kOrdinal) + roinf * u1d(i,j,k)*conjg(u1d(i,j,k))*kk*kk/2 + &
                        roinf * u2d(i,j,k)*conjg(u2d(i,j,k))*kk*kk/2 + &
                        roinf * u3d(i,j,k)*conjg(u3d(i,j,k))*kk*kk/2
          Eweak(kOrdinal) = Eweak(kOrdinal) + roinf * u1spe(i,j,k)*conjg(u1spe(i,j,k))*kk*kk/2 + &
                        roinf * u2spe(i,j,k)*conjg(u2spe(i,j,k))*kk*kk/2 + &
                        roinf * u3spe(i,j,k)*conjg(u3spe(i,j,k))*kk*kk/2
          Erho(kOrdinal) = Erho(kOrdinal) + real(rou1spe(i,j,k)*conjg(u1spe(i,j,k))+ &
                        rou2spe(i,j,k)*conjg(u2spe(i,j,k)) + &
                        rou3spe(i,j,k)*conjg(u3spe(i,j,k)))*kk*kk/2
          Pud(kOrdinal) = Pud(kOrdinal) - dimag(pspe(i,j,k)*dconjg(udspe(i,j,k))*kk)*kk*kk
          kn(kOrdinal) = kn(kOrdinal) + kk
        elseif(method == 3)then
          ES(kOrdinal) = ES(kOrdinal) + roinf * u1s(i,j,k)*conjg(u1s(i,j,k))/2 + &
                        roinf * u2s(i,j,k)*conjg(u2s(i,j,k))/2 + &
                        roinf * u3s(i,j,k)*conjg(u3s(i,j,k))/2
          ED(kOrdinal) = ED(kOrdinal) + roinf * u1d(i,j,k)*conjg(u1d(i,j,k))/2 + &
                        roinf * u2d(i,j,k)*conjg(u2d(i,j,k))/2 + &
                        roinf * u3d(i,j,k)*conjg(u3d(i,j,k))/2
          Eweak(kOrdinal) = Eweak(kOrdinal) + roinf * u1spe(i,j,k)*conjg(u1spe(i,j,k))/2 + &
                        roinf * u2spe(i,j,k)*conjg(u2spe(i,j,k))/2 + &
                        roinf * u3spe(i,j,k)*conjg(u3spe(i,j,k))/2
          Erho(kOrdinal) = Erho(kOrdinal) + real(rou1spe(i,j,k)*conjg(u1spe(i,j,k))+ &
                        rou2spe(i,j,k)*conjg(u2spe(i,j,k)) + &
                        rou3spe(i,j,k)*conjg(u3spe(i,j,k)))/2
          Pud(kOrdinal) = Pud(kOrdinal) - dimag(pspe(i,j,k)*dconjg(udspe(i,j,k))*kk)
        endif
      end if
      Edspe = Edspe +  roinf * u1d(i,j,k)*conjg(u1d(i,j,k))/2 + &
              roinf * u2d(i,j,k)*conjg(u2d(i,j,k))/2 + &
              roinf * u3d(i,j,k)*conjg(u3d(i,j,k))/2
      Esspe = Esspe + roinf * u1s(i,j,k)*conjg(u1s(i,j,k))/2 + &
              roinf * u2s(i,j,k)*conjg(u2s(i,j,k))/2 + &
              roinf * u3s(i,j,k)*conjg(u3s(i,j,k))/2
      IntLengthAbove = IntLengthAbove + roinf * u1d(i,j,k)*conjg(u1d(i,j,k))/2/kk + &
                      roinf * u2d(i,j,k)*conjg(u2d(i,j,k))/2/kk + &
                      roinf * u3d(i,j,k)*conjg(u3d(i,j,k))/2/kk + &
                      roinf * u1s(i,j,k)*conjg(u1s(i,j,k))/2/kk + &
                      roinf * u2s(i,j,k)*conjg(u2s(i,j,k))/2/kk + &
                      roinf * u3s(i,j,k)*conjg(u3s(i,j,k))/2/kk
      Pudspe = Pudspe + dimag(pspe(i,j,k)*dconjg(udspe(i,j,k))*(kk**2))/2
      Erhospe = Erhospe + real(rou1spe(i,j,k)*conjg(u1spe(i,j,k))+ &
                rou2spe(i,j,k)*conjg(u2spe(i,j,k)) + &
                rou3spe(i,j,k)*conjg(u3spe(i,j,k)))/2
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      ES(i) = psum(ES(i))
      ED(i) = psum(ED(i))
      Eweak(i) = psum(Eweak(i))
      Erho(i) = psum(Erho(i))
      Pud(i) = psum(Pud(i))
      Ecount(i) = psum(Ecount(i))
      if((method==1) .or. (method==2))then
        kn(i) = psum(kn(i))
        if(Ecount(i) .ne. 0)then
          ES(i) = ES(i)/Ecount(i)*4*pi
          ED(i) = ED(i)/Ecount(i)*4*pi
          Erho(i) = Erho(i)/Ecount(i)*4*pi
          Eweak(i) = Eweak(i)/Ecount(i)*4*pi
          Pud(i) = Pud(i)/Ecount(i)*4*pi
          kn(i) = kn(i)/Ecount(i)
        endif
      else
        kn(i) = real(i)
      endif
    end do
    !
    Edspe = psum(Edspe)
    Esspe = psum(Esspe)
    Erhospe = psum(Erhospe)
    Pudspe = psum(Pudspe)
    IntLengthAbove = psum(IntLengthAbove)
    !
    if(mpirank==0)  print*, '** Summation & average'
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1d,u1d)
    call fftw_mpi_execute_dft(backward_plan,u2d,u2d)
    call fftw_mpi_execute_dft(backward_plan,u3d,u3d)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    !
    !!!! Give S-C physical energy
    Edphy = 0.0d0
    Esphy = 0.0d0
    udusphy = 0.0d0
    Erhophy = 0.0d0
    Edmax = 0.0d0
    ReUsConjUdmax = -100.d0
    ReUsConjUdmin = 100.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      udud(i,j,k) = real(u1d(i,j,k))*real(u1d(i,j,k)) + real(u2d(i,j,k))*real(u2d(i,j,k)) + real(u3d(i,j,k))*real(u3d(i,j,k))
      udus(i,j,k) = real(u1d(i,j,k))*real(u1s(i,j,k)) + real(u2d(i,j,k))*real(u2s(i,j,k)) + real(u3d(i,j,k))*real(u3s(i,j,k))
      usus(i,j,k) = real(u1s(i,j,k))*real(u1s(i,j,k)) + real(u2s(i,j,k))*real(u2s(i,j,k)) + real(u3s(i,j,k))*real(u3s(i,j,k))
      Edmax = max(Edmax,udud(i,j,k))
      Edphy = Edphy + roinf * udud(i,j,k) / 2.d0
      Esphy = Esphy + roinf * usus(i,j,k) / 2.d0
      udusphy = udusphy + roinf * udus(i,j,k)
      Erhophy = Erhophy + (rho(i,j,k)-roinf)*(vel(i,j,k,1)**2 + vel(i,j,k,2)**2+ vel(i,j,k,3)**2) / 2.d0
      ReUsConjUdmax = max(ReUsConjUdmax,ReUsConjUd(i,j,k))
      ReUsConjUdmin = min(ReUsConjUdmin,ReUsConjUd(i,j,k))
    enddo
    enddo
    enddo
    !
    Erhophy = psum(Erhophy)/(1.d0*ia*ja*ka)
    Edphy = psum(Edphy)/(1.d0*ia*ja*ka)
    Esphy = psum(Esphy)/(1.d0*ia*ja*ka)
    udusphy = psum(udusphy)/(1.d0*ia*ja*ka)
    Edmax = pmax(Edmax)
    ReUsConjUdmax = pmax(ReUsConjUdmax)
    ReUsConjUdmin = pmin(ReUsConjUdmin)
    !
    if(mpirank==0)  print*, '** physical calculation finish'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec'//stepname//'.dat'
      else
        outfilename = 'pp/Espec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k ES ED Eweak Erho Pud')
      do i=0,allkmax
        if(Ecount(i)>1e-3) call listwrite(hand_a,kn(i),ES(i),ED(i),Eweak(i),Erho(i),Pud(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time Edphy Esphy udusphy Eweakphy Erhophy Edspe Esspe Pudspe Eweakspe Erhospe Edmax')! IntLen')
      call listwrite(hand_b,Edphy,Esphy,udusphy,Edphy+Esphy+udusphy, Erhophy,&
                    Edspe,Esspe,Pudspe,Edspe+Esspe,Erhospe,Edmax)!IntLengthAbove/(Ecspe+Esspe))
      !
      print*,' <<< '//outfilename//'... done.'
      !
      ! if (thefilenumb .ne. 0) then
      !   outfilename = 'pp/Espec_ReUsConjUd'//stepname//'.dat'
      ! else
      !   outfilename = 'pp/Espec_ReUsConjUd.dat'
      ! endif
      ! !
      ! call listinit(filename=outfilename,handle=hand_b, &
      !               firstline='nstep time ReUsConjUdmax ReUsConjUdmin')
      ! call listwrite(hand_b,ReUsConjUdmax,ReUsConjUdmin)
      ! !
      ! print*,' <<< '//outfilename//'... done.'
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Espec_mean'//stepname//'.dat'
      else
        outfilename = 'pp/Espec_mean.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_c, &
                    firstline='nstep time u1mean u2mean u3mean rhomean prsmean')
      call listwrite(hand_c,u1mean,u2mean,u3mean,rhomean,prsmean)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_rou1spe)
    call fftw_free(c_rou2spe)
    call fftw_free(c_rou3spe)
    call fftw_free(c_u1d)
    call fftw_free(c_u1s)
    call fftw_free(c_u2d)
    call fftw_free(c_u2s)
    call fftw_free(c_u3d)
    call fftw_free(c_u3s)
    call mpistop
    
    deallocate(ES,ED,Ecount,Eweak,Pud,kn,Erho)
    deallocate(udspe)
    deallocate(udud,udus,usus,ReUsConjUd)
    deallocate(k1,k2,k3)
    !
  end subroutine instantspectra3D
  !
  subroutine initparam2D
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput,readic
    use fftwlink
    use commvar,   only : im,jm,km,ia,ja,ka,ickmax,Reynolds
    use commarray, only : vel
    use hdf5io
    use utility,   only : listinit,listwrite
    use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    use fludyna,   only : miucal
    use solver,    only : refcal
    include 'fftw3-mpi.f03'
    !
    ! arguments
    character(len=128) :: infilename
    character(len=4) :: stepname
    integer :: allkmax
    character(len=1) :: modeio
    !
    real(8) :: u1mean,u2mean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    type(C_PTR) :: forward_plan, c_u1spe, c_u2spe
    real(8), allocatable, dimension(:,:) :: k1,k2
    !
    real(8) :: k,dk !wave number
    real(8) :: Espeall,Dissip,urms,tau,L,miu
    integer :: i,j,n
    
    !
    call readinput
    !
    call readic
    !
    call refcal
    if(mpirank==0)  print*, '** Referencecal done'
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:2))
    !
    !
    infilename='datin/flowini2d.h5'
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + vel(i,j,0,1)
        u2mean = u2mean + vel(i,j,0,2)
      enddo
    enddo
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
    !
    !
    dk = 1.d0
    !
    Espeall = 0.d0
    Dissip = 0.d0
    !
    do j=1,jm
    do i=1,im
        k=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        Espeall = Espeall + (u1spe(i,j)*dconjg(u1spe(i,j)))/2 + &
                            (u2spe(i,j)*dconjg(u2spe(i,j)))/2
        Dissip = Dissip + (u1spe(i,j)*dconjg(u1spe(i,j)))*(k**4) + &
                              (u2spe(i,j)*dconjg(u2spe(i,j)))*(k**4)
      end do
    end do
    !
    !
    miu = miucal(1.d0)/reynolds
    Dissip = psum(Dissip)*miu
    Espeall = psum(Espeall)
    !
    urms = sqrt(Espeall)
    L = Espeall**(1.d0/2.d0)/Dissip**(1.d0/3.d0)
    tau = L/urms
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !
    if(mpirank == 0) then
      !
      print *, 'urms:', urms
      print *, 'L:', L
      print *, 'tau:', tau
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call mpistop
    !
    deallocate(k1,k2)
    !
  end subroutine initparam2D
  !
  !
  subroutine initparam3D
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput,readic
    use fftwlink
    use commvar,only : im,jm,km,ia,ja,ka,ickmax
    use commarray, only: vel
    use hdf5io
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    use solver, only: refcal
    include 'fftw3-mpi.f03'
    !
    ! arguments
    character(len=128) :: infilename
    character(len=4) :: stepname
    integer :: allkmax
    character(len=1) :: modeio
    !
    real(8) :: u1mean,u2mean,u3mean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe
    type(C_PTR) :: forward_plan, c_u1spe, c_u2spe, c_u3spe
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    !
    real(8) :: kk,dk !wave number
    real(8) :: Espeall,TauAbove,urms,tau,L
    integer :: i,j,k,n
    !
    call readinput
    !
    call readic
    !
    call refcal
    if(mpirank==0)  print*, '** Referencecal done'
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(min(ia,ja),ka))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:", ka,"allkmax:",allkmax
    if(ka == 0) stop 'Please use instantspectra2D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    !
    allocate(vel(0:im,0:jm,0:km,1:3))
    !
    !
    !
    infilename='datin/flowini3d.h5'
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
    enddo
    enddo
    enddo
    !
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    !
    ! Planning
    forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, u1spe,u1spe, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2)-u2mean,0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3)-u3mean,0.d0,C_INTPTR_T);
      !
    enddo
    enddo
    enddo
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    !
    call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
    !
    dk = 1.0d0
    !
    Espeall = 0.d0
    TauAbove = 0.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      !
      Espeall = Espeall + (u1spe(i,j,k)*dconjg(u1spe(i,j,k)))/2 + &
                            (u2spe(i,j,k)*dconjg(u2spe(i,j,k)))/2+&
                            (u3spe(i,j,k)*dconjg(u3spe(i,j,k)))/2
      TauAbove = TauAbove + (u1spe(i,j,k)*dconjg(u1spe(i,j,k)))/2/kk + &
                            (u2spe(i,j,k)*dconjg(u2spe(i,j,k)))/2/kk + &
                            (u3spe(i,j,k)*dconjg(u3spe(i,j,k)))/2/kk
      !
    enddo
    enddo
    enddo
    !
    !
    if(mpirank==0)  print*, '** spectral decomposition calculation finish'
    !
    !!!! Give S-C spectra and spectral energy
    !
    TauAbove = psum(TauAbove)
    Espeall = psum(Espeall)
    !
    urms = sqrt(2.d0*Espeall/3.d0)
    L = TauAbove/urms/urms*pi/2
    tau = L/urms
    !
    if(mpirank==0)  print*, '** spectra calculation finish'
    !
    !
    if(mpirank == 0) then
      !
      print *, 'urms:', urms
      print *, 'L1:', sqrt(2*pi)/ickmax
      print *, 'L:', L
      print *, 'tau:', tau
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call mpistop
    !
    deallocate(k1,k2,k3)
    !
  end subroutine initparam3D
  !
  subroutine instanttriad2D(thefilenumb,method)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use solver,    only : refcal
    use utility,  only : listinit,listwrite
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    integer,intent(in) :: method
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe,pspe
    complex(8), allocatable, dimension(:,:) :: usspe,ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: theta,u1c,u2c,u1s,u2s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11,T12,T21,T22,Tstheta1,Tstheta2,Tdtheta1,Tdtheta2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: T11s,T12s,T21s,T22s,T11d,T12d,T21d,T22d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tdd,Tcount,Tss,Tdstheta,Tddtheta,kn
    real(8) :: TsO,TsthetaO,TdO,TddO,TssO,TdsthetaO,TddthetaO
    real(8), allocatable, dimension(:,:) :: k1,k2,u1,u2
    integer :: allkmax, kOrdinal
    real(8) :: kk,dk,p,q,kx,ky,lambda !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t
    real(8) :: k2Ts, k2Tc,Tsall,Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_pspe, c_theta,c_u1c,c_u2c,c_u1s,c_u2s
    type(C_PTR) :: c_T11,c_T12,c_T21,c_T22,c_Tstheta1,c_Tstheta2,c_Tdtheta1,c_Tdtheta2
    type(C_PTR) :: c_T11s,c_T12s,c_T21s,c_T22s,c_T11d,c_T12d,c_T21d,c_T22d
    !
    ! Initialization
    if((method < 1 ) .or. (method >3))then
      stop "Error! method @ instantspectra2D problem, not among 1, 2, 3"
    else
      if(mpirank==0) print *, "Using method", method
    endif
    !
    call readinput
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    modeio='h'
    !
    if(ka .ne. 0) stop 'Please use instantspectra3D'
    !
    dk = 1.d0
    !
    if(method==1)then
      lambda = 1.21d0
      allkmax = ceiling(log(real(sqrt(2.d0)/3*min(ia,ja))/dk)/log(lambda))
    elseif((method == 2) .or. (method==3))then
      allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
    endif
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"knumber:",allkmax
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm),u2(1:im,1:jm))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do i=1,im
    do j=1,jm
      u1mean = u1mean + vel(i,j,0,1)
      u2mean = u2mean + vel(i,j,0,2)
      rhomean = rhomean + rho(i,j,0)
      prsmean = prsmean + prs(i,j,0)
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja)
    u2mean = psum(u2mean) / (1.d0*ia*ja)
    prsmean = psum(prsmean) / (1.d0*ia*ja)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do j=1,jm
    do i=1,im
      !
      u1(i,j) = vel(i,j,0,1)-u1mean
      u2(i,j) = vel(i,j,0,2)-u2mean
      u1spe(i,j)=CMPLX(u1(i,j),0.d0,C_INTPTR_T)
      u2spe(i,j)=CMPLX(u2(i,j),0.d0,C_INTPTR_T)
      pspe(i,j)=CMPLX(prs(i,j,0)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      pspe(i,j)=pspe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
    !
    allocate(usspe(1:im,1:jm),ucspe(1:im,1:jm))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    !!!! Do S-C decomposition
    do j=1,jm
      do i=1,im
        !
        kk=dsqrt(k1(i,j)**2+k2(i,j)**2+1.d-15)
        usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
        ucspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
        !
        u1c(i,j)=  ucspe(i,j)*k1(i,j)/kk
        u2c(i,j)=  ucspe(i,j)*k2(i,j)/kk
        u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
        u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
        !
        theta(i,j)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j) * kk)
        !
      end do
    end do
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw])
    !
    ! This part is in physical space
    do j=1,jm
    do i=1,im
      T11(i,j) = u1(i,j)*u1(i,j)
      T12(i,j) = u1(i,j)*u2(i,j)
      T21(i,j) = u2(i,j)*u1(i,j)
      T22(i,j) = u2(i,j)*u2(i,j)
      !
      T11s(i,j) = u1s(i,j)*u1s(i,j)
      T12s(i,j) = u1s(i,j)*u2s(i,j)
      T21s(i,j) = u2s(i,j)*u1s(i,j)
      T22s(i,j) = u2s(i,j)*u2s(i,j)
      !
      T11d(i,j) = u1c(i,j)*u1c(i,j)
      T12d(i,j) = u1c(i,j)*u2c(i,j)
      T21d(i,j) = u2c(i,j)*u1c(i,j)
      T22d(i,j) = u2c(i,j)*u2c(i,j)
      !
      Tstheta1(i,j) = theta(i,j)*u1s(i,j)
      Tstheta2(i,j) = theta(i,j)*u2s(i,j)
      Tdtheta1(i,j) = theta(i,j)*u1c(i,j)
      Tdtheta2(i,j) = theta(i,j)*u2c(i,j)
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    !
    do j=1,jm
    do i=1,im
      T11(i,j)=T11(i,j)/(1.d0*ia*ja)
      T12(i,j)=T12(i,j)/(1.d0*ia*ja)
      T21(i,j)=T21(i,j)/(1.d0*ia*ja)
      T22(i,j)=T22(i,j)/(1.d0*ia*ja)
      !
      T11s(i,j)=T11s(i,j)/(1.d0*ia*ja)
      T12s(i,j)=T12s(i,j)/(1.d0*ia*ja)
      T21s(i,j)=T21s(i,j)/(1.d0*ia*ja)
      T22s(i,j)=T22s(i,j)/(1.d0*ia*ja)
      !
      T11d(i,j)=T11d(i,j)/(1.d0*ia*ja)
      T12d(i,j)=T12d(i,j)/(1.d0*ia*ja)
      T21d(i,j)=T21d(i,j)/(1.d0*ia*ja)
      T22d(i,j)=T22d(i,j)/(1.d0*ia*ja)
      !
      Tstheta1(i,j)=Tstheta1(i,j)/(1.d0*ia*ja)
      Tstheta2(i,j)=Tstheta2(i,j)/(1.d0*ia*ja)
      Tdtheta1(i,j)=Tdtheta1(i,j)/(1.d0*ia*ja)
      Tdtheta2(i,j)=Tdtheta2(i,j)/(1.d0*ia*ja)
      !
      u1s(i,j)=u1s(i,j)/(1.d0*ia*ja)
      u2s(i,j)=u2s(i,j)/(1.d0*ia*ja)
      u1c(i,j)=u1c(i,j)/(1.d0*ia*ja)
      u2c(i,j)=u2c(i,j)/(1.d0*ia*ja)
    enddo
    enddo
    !
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    allocate(kn(0:allkmax))
    !
    Ts= 0.0d0
    Tstheta= 0.0d0
    Td= 0.0d0
    Tddtheta= 0.0d0
    Tdstheta= 0.0d0
    Tcount=0
    Tss=0.d0
    Tdd=0.d0
    kn=0.d0
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j)
      ky = k2(i,j)
      kk=dsqrt(kx**2+ky**2+1.d-15)
      kOrdinal = kint(kk,dk,method,lambda)
      !
      if(kOrdinal==0)then
        print *, u1s(i,j), u2s(i,j), T11s(i,j), T12s(i,j), T21s(i,j), T22s(i,j)
      endif
      !
      TsO       = ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11(i,j))) + &
                  ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11(i,j))) + &
                  ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12(i,j))) + &
                  ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12(i,j))) + &
                  ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21(i,j))) + &
                  ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21(i,j))) + &
                  ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22(i,j))) + &
                  ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22(i,j)))
      !
      TssO       =ProjectP3(1,1,1,kx,ky) * dimag(u1s(i,j)*conjg(T11s(i,j))) + &
                  ProjectP3(2,1,1,kx,ky) * dimag(u2s(i,j)*conjg(T11s(i,j))) + &
                  ProjectP3(1,1,2,kx,ky) * dimag(u1s(i,j)*conjg(T12s(i,j))) + &
                  ProjectP3(2,1,2,kx,ky) * dimag(u2s(i,j)*conjg(T12s(i,j))) + &
                  ProjectP3(1,2,1,kx,ky) * dimag(u1s(i,j)*conjg(T21s(i,j))) + &
                  ProjectP3(2,2,1,kx,ky) * dimag(u2s(i,j)*conjg(T21s(i,j))) + &
                  ProjectP3(1,2,2,kx,ky) * dimag(u1s(i,j)*conjg(T22s(i,j))) + &
                  ProjectP3(2,2,2,kx,ky) * dimag(u2s(i,j)*conjg(T22s(i,j)))
      !
      TsthetaO   =ProjectP2(1,1,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta1(i,j)+Tdtheta1(i,j))) + &
                  ProjectP2(1,2,kx,ky) * dreal(u1s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j))) + &
                  ProjectP2(2,1,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta1(i,j)+Tstheta1(i,j))) + &
                  ProjectP2(2,2,kx,ky) * dreal(u2s(i,j)*conjg(Tstheta2(i,j)+Tdtheta2(i,j)))
      !
      TdO       = ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11(i,j))) + &
                  ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11(i,j))) + &
                  ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12(i,j))) + &
                  ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12(i,j))) + &
                  ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21(i,j))) + &
                  ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21(i,j))) + &
                  ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22(i,j))) + &
                  ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22(i,j)))
      !
      TddO      = ProjectPi3(1,1,1,kx,ky) * dimag(u1c(i,j)*conjg(T11d(i,j))) + &
                  ProjectPi3(2,1,1,kx,ky) * dimag(u2c(i,j)*conjg(T11d(i,j))) + &
                  ProjectPi3(1,1,2,kx,ky) * dimag(u1c(i,j)*conjg(T12d(i,j))) + &
                  ProjectPi3(2,1,2,kx,ky) * dimag(u2c(i,j)*conjg(T12d(i,j))) + &
                  ProjectPi3(1,2,1,kx,ky) * dimag(u1c(i,j)*conjg(T21d(i,j))) + &
                  ProjectPi3(2,2,1,kx,ky) * dimag(u2c(i,j)*conjg(T21d(i,j))) + &
                  ProjectPi3(1,2,2,kx,ky) * dimag(u1c(i,j)*conjg(T22d(i,j))) + &
                  ProjectPi3(2,2,2,kx,ky) * dimag(u2c(i,j)*conjg(T22d(i,j)))
      !
      TdsthetaO = ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta1(i,j))) + &
                  ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tstheta2(i,j))) + &
                  ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta1(i,j))) + &
                  ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tstheta2(i,j)))
      !
      TddthetaO = ProjectPi2(1,1,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta1(i,j))) + &
                  ProjectPi2(1,2,kx,ky) * dreal(u1c(i,j)*conjg(Tdtheta2(i,j))) + &
                  ProjectPi2(2,1,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta1(i,j))) + &
                  ProjectPi2(2,2,kx,ky) * dreal(u2c(i,j)*conjg(Tdtheta2(i,j)))
      !
      if(kOrdinal<= allkmax)then
        Tcount(kOrdinal)   = Tcount(kOrdinal)   + 1
        if((method == 1) .or. (method == 2))then
          Ts(kOrdinal)       = Ts(kOrdinal)       - TsO*kk
          Tstheta(kOrdinal)  = Tstheta(kOrdinal)  + TsthetaO*kk
          Td(kOrdinal)       = Td(kOrdinal)       - TdO*kk
          Tddtheta(kOrdinal) = Tddtheta(kOrdinal) + TddthetaO*kk
          Tdstheta(kOrdinal) = Tdstheta(kOrdinal) + TdsthetaO*kk
          Tss(kOrdinal)      = Tss(kOrdinal)      - TssO*kk
          Tdd(kOrdinal)      = Tdd(kOrdinal)      - TddO*kk
          kn(kOrdinal)       = kn(kOrdinal)       + kk
        elseif(method == 3)then
          Ts(kOrdinal)       = Ts(kOrdinal)       - TsO
          Tstheta(kOrdinal)  = Tstheta(kOrdinal)  + TsthetaO
          Td(kOrdinal)       = Td(kOrdinal)       - TdO
          Tddtheta(kOrdinal) = Tddtheta(kOrdinal) + TddthetaO
          Tdstheta(kOrdinal) = Tdstheta(kOrdinal) + TdsthetaO
          Tss(kOrdinal)      = Tss(kOrdinal)      - TssO
          Tdd(kOrdinal)      = Tdd(kOrdinal)      - TddO
        endif
      end if
      !
      Tsall = Tsall + TsO
      Tssall= Tssall+ TssO
      Tcall = Tcall + TdO
      Tccall= Tccall+ TddO
      k2Ts  = k2Ts  + kk**2 * TsO
      k2Tss = k2Tss + kk**2 * TssO
      k2Tc  = k2Tc  + kk**2 * TdO
      k2Tcc = k2Tcc + kk**2 * TddO
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = psum(Tss(i))
      Tdd(i) = psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
      if((method == 1) .or. (method == 2))then
        kn(i) = psum(kn(i))
        if(Tcount(i) .ne. 0)then
          Ts(i) = Ts(i)/Tcount(i)*2*pi
          Tstheta(i) = Tstheta(i)/Tcount(i)*2*pi
          Td(i) = Td(i)/Tcount(i)*2*pi
          Tdstheta(i) = Tdstheta(i)/Tcount(i)*2*pi
          Tddtheta(i) = Tddtheta(i)/Tcount(i)*2*pi
          Tss(i) = Tss(i)/Tcount(i)*2*pi
          Tdd(i) = Tdd(i)/Tcount(i)*2*pi
          kn(i) = kn(i)/tcount(i)
        endif
      else
        kn(i) = real(i)
      endif
    end do
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank==0)  print*, '** Summation & average'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        if(Tcount(i)>1e-3) call listwrite(hand_a,kn(i),Ts(i),Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call mpistop
    !
  end subroutine instanttriad2D
  !
  subroutine instanttriad3D(thefilenumb,method)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,   only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only : vel, rho, prs
    use hdf5io
    use solver,    only : refcal
    use utility,   only : listinit,listwrite
    use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    integer,intent(in) :: method
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,u3mean,rhomean,prsmean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe
    complex(8), allocatable, dimension(:,:,:) :: ucspe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: theta,u1c,u2c,u3c,u1s,u2s,u3s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11,T12,T13,T21,T22,T23,T31,T32,T33
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Tstheta1,Tstheta2,Tstheta3,Tdtheta1,Tdtheta2,Tdtheta3
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11s,T12s,T13s,T21s,T22s,T23s,T31s,T32s,T33s
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: T11d,T12d,T13d,T21d,T22d,T23d,T31d,T32d,T33d
    real(8), allocatable, dimension(:) :: Ts,Tstheta,Td,Tddtheta,Tdstheta,Tcount,Tss,Tdd,kn
    real(8) :: TsO,TsthetaO,TdO,TddO,TssO,TdsthetaO,TddthetaO
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3,u1,u2,u3
    integer :: allkmax, kOrdinal
    real(8) :: kk,dk,p,q,kx,ky,kz,lambda !wave number
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n,s,t,k
    real(8) :: k2Ts, k2Tc, Tsall, Tcall, Tssall, Tccall, k2Tss, k2Tcc
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe
    type(C_PTR) :: c_theta,c_u1c,c_u2c,c_u3c,c_u1s,c_u2s,c_u3s
    type(C_PTR) :: c_T11,c_T12,c_T13,c_T21,c_T22,c_T23,c_T31,c_T32,c_T33
    type(C_PTR) :: c_Tstheta1,c_Tstheta2,c_Tstheta3,c_Tdtheta1,c_Tdtheta2,c_Tdtheta3
    type(C_PTR) :: c_T11s,c_T12s,c_T13s,c_T21s,c_T22s,c_T23s,c_T31s,c_T32s,c_T33s
    type(C_PTR) :: c_T11d,c_T12d,c_T13d,c_T21d,c_T22d,c_T23d,c_T31d,c_T32d,c_T33d
    !
    ! Initialization
    if((method < 1 ) .or. (method >3))then
      stop "Error! method @ instantspectra2D problem, not among 1, 2, 3"
    else
      if(mpirank==0) print *, "Using method", method
    endif
    !
    call readinput
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    modeio='h'
    !
    if(ka==0) stop 'Please use instantspectra2D'
    !
    dk = 1.0d0
    !
    if(method==1)then
      lambda = 1.21d0
      allkmax = ceiling(log(real(sqrt(2.d0)/3*min(ia,ja,ka))/dk)/log(lambda))
    elseif((method == 2) .or. (method==3))then
      allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja,ka))/dk)
    endif
    !
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka,"knumber:",allkmax
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    allocate(u1(1:im,1:jm,1:km),u2(1:im,1:jm,1:km),u3(1:im,1:jm,1:km))
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    u3mean = 0.0d0
    rhomean = 0.0d0
    prsmean = 0.0d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      u1mean = u1mean + vel(i,j,k,1)
      u2mean = u2mean + vel(i,j,k,2)
      u3mean = u3mean + vel(i,j,k,3)
      rhomean = rhomean + rho(i,j,k)
      prsmean = prsmean + prs(i,j,k)
    enddo
    enddo
    enddo
    !
    rhomean = psum(rhomean) / (1.0d0*ia*ja*ka)
    u1mean = psum(u1mean) / (1.d0*ia*ja*ka)
    u2mean = psum(u2mean) / (1.d0*ia*ja*ka)
    u3mean = psum(u3mean) / (1.d0*ia*ja*ka)
    prsmean = psum(prsmean) / (1.d0*ia*ja*ka)
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean, 'u3mean=',u3mean, 'prsmean=',prsmean
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    !
    forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1(i,j,k) = vel(i,j,k,1)-u1mean
      u2(i,j,k) = vel(i,j,k,2)-u2mean
      u3(i,j,k) = vel(i,j,k,3)-u3mean
      u1spe(i,j,k)=CMPLX(u1(i,j,k),0.d0,C_INTPTR_T)
      u2spe(i,j,k)=CMPLX(u2(i,j,k),0.d0,C_INTPTR_T)
      u3spe(i,j,k)=CMPLX(u3(i,j,k),0.d0,C_INTPTR_T)
      pspe(i,j,k)=CMPLX(prs(i,j,k)-prsmean,0.d0,C_INTPTR_T)
      !
    end do
    end do
    end do
    !
    !!!! Do FFT
    !
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    call fftw_mpi_execute_dft(forward_plan,u3spe,u3spe)
    call fftw_mpi_execute_dft(forward_plan,pspe,pspe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      !
    end do
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
    !
    allocate(ucspe(1:im,1:jm,1:km))
    !
    c_theta = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_theta, theta, [imfftw,jmfftw,kmfftw])
    c_u1s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1s, u1s, [imfftw,jmfftw,kmfftw])
    c_u2s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2s, u2s, [imfftw,jmfftw,kmfftw])
    c_u3s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3s, u3s, [imfftw,jmfftw,kmfftw])
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    !!!! Do S-C decomposition
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
      !
      ucspe(i,j,k) = k1(i,j,k)/kk * u1spe(i,j,k) + k2(i,j,k)/kk * u2spe(i,j,k) + k3(i,j,k)/kk * u3spe(i,j,k)
      u1c(i,j,k)   =  k1(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k1(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k1(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u2c(i,j,k)   =  k2(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k2(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k2(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u3c(i,j,k)   =  k3(i,j,k)*k1(i,j,k)/(kk**2) * u1spe(i,j,k) + k3(i,j,k)*k2(i,j,k)/(kk**2) * u2spe(i,j,k) &
                + k3(i,j,k)*k3(i,j,k)/(kk**2) * u3spe(i,j,k)
      u1s(i,j,k)   =  u1spe(i,j,k) - u1c(i,j,k)
      u2s(i,j,k)   =  u2spe(i,j,k) - u2c(i,j,k)
      u3s(i,j,k)   =  u3spe(i,j,k) - u3c(i,j,k)
      !
      theta(i,j,k)= CMPLX(0.d0,1.d0,C_INTPTR_T) * (ucspe(i,j,k) * kk)
      !
    end do
    end do
    enddo
    !
    !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,theta,theta)
    call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(backward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(backward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(backward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(backward_plan,u3c,u3c)
    !
    c_T11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11,T11, [imfftw,jmfftw,kmfftw])
    c_T12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12,T12, [imfftw,jmfftw,kmfftw])
    c_T13 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13,T13, [imfftw,jmfftw,kmfftw])
    c_T21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21,T21, [imfftw,jmfftw,kmfftw])
    c_T22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22,T22, [imfftw,jmfftw,kmfftw])
    c_T23 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23,T23, [imfftw,jmfftw,kmfftw])
    c_T31 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31,T31, [imfftw,jmfftw,kmfftw])
    c_T32 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32,T32, [imfftw,jmfftw,kmfftw])
    c_T33 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33,T33, [imfftw,jmfftw,kmfftw])
    !
    c_T11s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11s,T11s, [imfftw,jmfftw,kmfftw])
    c_T12s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12s,T12s, [imfftw,jmfftw,kmfftw])
    c_T13s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13s,T13s, [imfftw,jmfftw,kmfftw])
    c_T21s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21s,T21s, [imfftw,jmfftw,kmfftw])
    c_T22s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22s,T22s, [imfftw,jmfftw,kmfftw])
    c_T23s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23s,T23s, [imfftw,jmfftw,kmfftw])
    c_T31s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31s,T31s, [imfftw,jmfftw,kmfftw])
    c_T32s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32s,T32s, [imfftw,jmfftw,kmfftw])
    c_T33s = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33s,T33s, [imfftw,jmfftw,kmfftw])
    !
    c_T11d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T11d,T11d, [imfftw,jmfftw,kmfftw])
    c_T12d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T12d,T12d, [imfftw,jmfftw,kmfftw])
    c_T13d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T13d,T13d, [imfftw,jmfftw,kmfftw])
    c_T21d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T21d,T21d, [imfftw,jmfftw,kmfftw])
    c_T22d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T22d,T22d, [imfftw,jmfftw,kmfftw])
    c_T23d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T23d,T23d, [imfftw,jmfftw,kmfftw])
    c_T31d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T31d,T31d, [imfftw,jmfftw,kmfftw])
    c_T32d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T32d,T32d, [imfftw,jmfftw,kmfftw])
    c_T33d = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_T33d,T33d, [imfftw,jmfftw,kmfftw])
    !
    c_Tstheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta1,Tstheta1, [imfftw,jmfftw,kmfftw])
    c_Tstheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta2,Tstheta2, [imfftw,jmfftw,kmfftw])
    c_Tstheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tstheta3,Tstheta3, [imfftw,jmfftw,kmfftw])
    c_Tdtheta1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta1,Tdtheta1, [imfftw,jmfftw,kmfftw])
    c_Tdtheta2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta2,Tdtheta2, [imfftw,jmfftw,kmfftw])
    c_Tdtheta3 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_Tdtheta3,Tdtheta3, [imfftw,jmfftw,kmfftw])
    !
    ! This part is in physical space
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k) = u1(i,j,k)*u1(i,j,k)
      T12(i,j,k) = u1(i,j,k)*u2(i,j,k)
      T13(i,j,k) = u1(i,j,k)*u3(i,j,k)
      T21(i,j,k) = u2(i,j,k)*u1(i,j,k)
      T22(i,j,k) = u2(i,j,k)*u2(i,j,k)
      T23(i,j,k) = u2(i,j,k)*u3(i,j,k)
      T31(i,j,k) = u3(i,j,k)*u1(i,j,k)
      T32(i,j,k) = u3(i,j,k)*u2(i,j,k)
      T33(i,j,k) = u3(i,j,k)*u3(i,j,k)
      !
      T11s(i,j,k) = u1s(i,j,k)*u1s(i,j,k)
      T12s(i,j,k) = u1s(i,j,k)*u2s(i,j,k)
      T13s(i,j,k) = u1s(i,j,k)*u3s(i,j,k)
      T21s(i,j,k) = u2s(i,j,k)*u1s(i,j,k)
      T22s(i,j,k) = u2s(i,j,k)*u2s(i,j,k)
      T23s(i,j,k) = u2s(i,j,k)*u3s(i,j,k)
      T31s(i,j,k) = u3s(i,j,k)*u1s(i,j,k)
      T32s(i,j,k) = u3s(i,j,k)*u2s(i,j,k)
      T33s(i,j,k) = u3s(i,j,k)*u3s(i,j,k)
      !
      T11d(i,j,k) = u1c(i,j,k)*u1c(i,j,k)
      T12d(i,j,k) = u1c(i,j,k)*u2c(i,j,k)
      T13d(i,j,k) = u1c(i,j,k)*u3c(i,j,k)
      T21d(i,j,k) = u2c(i,j,k)*u1c(i,j,k)
      T22d(i,j,k) = u2c(i,j,k)*u2c(i,j,k)
      T23d(i,j,k) = u2c(i,j,k)*u3c(i,j,k)
      T31d(i,j,k) = u3c(i,j,k)*u1c(i,j,k)
      T32d(i,j,k) = u3c(i,j,k)*u2c(i,j,k)
      T33d(i,j,k) = u3c(i,j,k)*u3c(i,j,k)
      !
      Tstheta1(i,j,k) = theta(i,j,k)*u1s(i,j,k)
      Tstheta2(i,j,k) = theta(i,j,k)*u2s(i,j,k)
      Tstheta3(i,j,k) = theta(i,j,k)*u3s(i,j,k)
      Tdtheta1(i,j,k) = theta(i,j,k)*u1c(i,j,k)
      Tdtheta2(i,j,k) = theta(i,j,k)*u2c(i,j,k)
      Tdtheta3(i,j,k) = theta(i,j,k)*u3c(i,j,k)
    enddo
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(forward_plan,T11,T11)
    call fftw_mpi_execute_dft(forward_plan,T12,T12)
    call fftw_mpi_execute_dft(forward_plan,T13,T13)
    call fftw_mpi_execute_dft(forward_plan,T21,T21)
    call fftw_mpi_execute_dft(forward_plan,T22,T22)
    call fftw_mpi_execute_dft(forward_plan,T23,T23)
    call fftw_mpi_execute_dft(forward_plan,T31,T31)
    call fftw_mpi_execute_dft(forward_plan,T32,T32)
    call fftw_mpi_execute_dft(forward_plan,T33,T33)
    !
    call fftw_mpi_execute_dft(forward_plan,T11s,T11s)
    call fftw_mpi_execute_dft(forward_plan,T12s,T12s)
    call fftw_mpi_execute_dft(forward_plan,T13s,T13s)
    call fftw_mpi_execute_dft(forward_plan,T21s,T21s)
    call fftw_mpi_execute_dft(forward_plan,T22s,T22s)
    call fftw_mpi_execute_dft(forward_plan,T23s,T23s)
    call fftw_mpi_execute_dft(forward_plan,T31s,T31s)
    call fftw_mpi_execute_dft(forward_plan,T32s,T32s)
    call fftw_mpi_execute_dft(forward_plan,T33s,T33s)
    !
    call fftw_mpi_execute_dft(forward_plan,T11d,T11d)
    call fftw_mpi_execute_dft(forward_plan,T12d,T12d)
    call fftw_mpi_execute_dft(forward_plan,T13d,T13d)
    call fftw_mpi_execute_dft(forward_plan,T21d,T21d)
    call fftw_mpi_execute_dft(forward_plan,T22d,T22d)
    call fftw_mpi_execute_dft(forward_plan,T23d,T23d)
    call fftw_mpi_execute_dft(forward_plan,T31d,T31d)
    call fftw_mpi_execute_dft(forward_plan,T32d,T32d)
    call fftw_mpi_execute_dft(forward_plan,T33d,T33d)
    !
    call fftw_mpi_execute_dft(forward_plan,Tstheta1,Tstheta1)
    call fftw_mpi_execute_dft(forward_plan,Tstheta2,Tstheta2)
    call fftw_mpi_execute_dft(forward_plan,Tstheta3,Tstheta3)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta1,Tdtheta1)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta2,Tdtheta2)
    call fftw_mpi_execute_dft(forward_plan,Tdtheta3,Tdtheta3)
    !
    call fftw_mpi_execute_dft(forward_plan,u1s,u1s)
    call fftw_mpi_execute_dft(forward_plan,u2s,u2s)
    call fftw_mpi_execute_dft(forward_plan,u3s,u3s)
    call fftw_mpi_execute_dft(forward_plan,u1c,u1c)
    call fftw_mpi_execute_dft(forward_plan,u2c,u2c)
    call fftw_mpi_execute_dft(forward_plan,u3c,u3c)
    !
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      T11(i,j,k)=T11(i,j,k)/(1.d0*ia*ja*ka)
      T12(i,j,k)=T12(i,j,k)/(1.d0*ia*ja*ka)
      T13(i,j,k)=T13(i,j,k)/(1.d0*ia*ja*ka)
      T21(i,j,k)=T21(i,j,k)/(1.d0*ia*ja*ka)
      T22(i,j,k)=T22(i,j,k)/(1.d0*ia*ja*ka)
      T23(i,j,k)=T23(i,j,k)/(1.d0*ia*ja*ka)
      T31(i,j,k)=T31(i,j,k)/(1.d0*ia*ja*ka)
      T32(i,j,k)=T32(i,j,k)/(1.d0*ia*ja*ka)
      T33(i,j,k)=T33(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11s(i,j,k)=T11s(i,j,k)/(1.d0*ia*ja*ka)
      T12s(i,j,k)=T12s(i,j,k)/(1.d0*ia*ja*ka)
      T13s(i,j,k)=T13s(i,j,k)/(1.d0*ia*ja*ka)
      T21s(i,j,k)=T21s(i,j,k)/(1.d0*ia*ja*ka)
      T22s(i,j,k)=T22s(i,j,k)/(1.d0*ia*ja*ka)
      T23s(i,j,k)=T23s(i,j,k)/(1.d0*ia*ja*ka)
      T31s(i,j,k)=T31s(i,j,k)/(1.d0*ia*ja*ka)
      T32s(i,j,k)=T32s(i,j,k)/(1.d0*ia*ja*ka)
      T33s(i,j,k)=T33s(i,j,k)/(1.d0*ia*ja*ka)
      !
      T11d(i,j,k)=T11d(i,j,k)/(1.d0*ia*ja*ka)
      T12d(i,j,k)=T12d(i,j,k)/(1.d0*ia*ja*ka)
      T13d(i,j,k)=T13d(i,j,k)/(1.d0*ia*ja*ka)
      T21d(i,j,k)=T21d(i,j,k)/(1.d0*ia*ja*ka)
      T22d(i,j,k)=T22d(i,j,k)/(1.d0*ia*ja*ka)
      T23d(i,j,k)=T23d(i,j,k)/(1.d0*ia*ja*ka)
      T31d(i,j,k)=T31d(i,j,k)/(1.d0*ia*ja*ka)
      T32d(i,j,k)=T32d(i,j,k)/(1.d0*ia*ja*ka)
      T33d(i,j,k)=T33d(i,j,k)/(1.d0*ia*ja*ka)
      !
      Tstheta1(i,j,k)=Tstheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta2(i,j,k)=Tstheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tstheta3(i,j,k)=Tstheta3(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta1(i,j,k)=Tdtheta1(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta2(i,j,k)=Tdtheta2(i,j,k)/(1.d0*ia*ja*ka)
      Tdtheta3(i,j,k)=Tdtheta3(i,j,k)/(1.d0*ia*ja*ka)
      !
      u1s(i,j,k)=u1s(i,j,k)/(1.d0*ia*ja*ka)
      u2s(i,j,k)=u2s(i,j,k)/(1.d0*ia*ja*ka)
      u3s(i,j,k)=u3s(i,j,k)/(1.d0*ia*ja*ka)
      u1c(i,j,k)=u1c(i,j,k)/(1.d0*ia*ja*ka)
      u2c(i,j,k)=u2c(i,j,k)/(1.d0*ia*ja*ka)
      u3c(i,j,k)=u3c(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !
    allocate(Ts(0:allkmax),Tstheta(0:allkmax),Td(0:allkmax),Tdstheta(0:allkmax))
    allocate(Tddtheta(0:allkmax),Tcount(0:allkmax),Tss(0:allkmax),Tdd(0:allkmax))
    allocate(kn(0:allkmax))
    !
    Ts= 0.0d0
    Tstheta= 0.0d0
    Td= 0.0d0
    Tddtheta= 0.0d0
    Tdstheta= 0.0d0
    Tcount=0
    Tss=0.d0
    Tdd=0.d0
    kn=0.d0
    k2Ts = 0.d0
    k2Tc = 0.d0
    Tsall = 0.d0
    Tcall = 0.d0
    Tssall = 0.d0
    Tccall = 0.d0
    k2Tss = 0.d0
    k2Tcc = 0.d0
    !
    do k = 1,km
    do j = 1,jm
    do i = 1,im
      kx = k1(i,j,k)
      ky = k2(i,j,k)
      kz = k3(i,j,k)
      kk=dsqrt(kx**2+ky**2+kz**2+1.d-15)
      kOrdinal = kint(kk,dk,method,lambda)
      !
      TsO         = ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33(i,j,k))) + &
                    ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33(i,j,k))) + &
                    ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33(i,j,k)))
        !
      TssO       =  ProjectP3(1,1,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T11s(i,j,k))) + &
                    ProjectP3(2,1,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T11s(i,j,k))) + &
                    ProjectP3(3,1,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T11s(i,j,k))) + &
                    ProjectP3(1,1,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T12s(i,j,k))) + &
                    ProjectP3(2,1,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T12s(i,j,k))) + &
                    ProjectP3(3,1,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T12s(i,j,k))) + &
                    ProjectP3(1,1,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T13s(i,j,k))) + &
                    ProjectP3(2,1,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T13s(i,j,k))) + &
                    ProjectP3(3,1,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T13s(i,j,k))) + &
                    ProjectP3(1,2,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T21s(i,j,k))) + &
                    ProjectP3(2,2,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T21s(i,j,k))) + &
                    ProjectP3(3,2,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T21s(i,j,k))) + &
                    ProjectP3(1,2,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T22s(i,j,k))) + &
                    ProjectP3(2,2,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T22s(i,j,k))) + &
                    ProjectP3(3,2,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T22s(i,j,k))) + &
                    ProjectP3(1,2,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T23s(i,j,k))) + &
                    ProjectP3(2,2,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T23s(i,j,k))) + &
                    ProjectP3(3,2,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T23s(i,j,k))) + &
                    ProjectP3(1,3,1,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T31s(i,j,k))) + &
                    ProjectP3(2,3,1,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T31s(i,j,k))) + &
                    ProjectP3(3,3,1,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T31s(i,j,k))) + &
                    ProjectP3(1,3,2,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T32s(i,j,k))) + &
                    ProjectP3(2,3,2,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T32s(i,j,k))) + &
                    ProjectP3(3,3,2,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T32s(i,j,k))) + &
                    ProjectP3(1,3,3,kx,ky,kz) * dimag(u1s(i,j,k)*conjg(T33s(i,j,k))) + &
                    ProjectP3(2,3,3,kx,ky,kz) * dimag(u2s(i,j,k)*conjg(T33s(i,j,k))) + &
                    ProjectP3(3,3,3,kx,ky,kz) * dimag(u3s(i,j,k)*conjg(T33s(i,j,k)))
        !
      TsthetaO    = ProjectP2(1,1,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                    ProjectP2(1,2,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                    ProjectP2(1,3,kx,ky,kz) * dreal(u1s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
                    ProjectP2(2,1,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                    ProjectP2(2,2,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                    ProjectP2(2,3,kx,ky,kz) * dreal(u2s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k))) + &
                    ProjectP2(3,1,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta1(i,j,k)+Tdtheta1(i,j,k))) + &
                    ProjectP2(3,2,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta2(i,j,k)+Tdtheta2(i,j,k))) + &
                    ProjectP2(3,3,kx,ky,kz) * dreal(u3s(i,j,k)*conjg(Tstheta3(i,j,k)+Tdtheta3(i,j,k)))
        !
      TdO         = ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11(i,j,k))) + &
                    ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12(i,j,k))) + &
                    ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13(i,j,k))) + &
                    ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21(i,j,k))) + &
                    ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22(i,j,k))) + &
                    ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23(i,j,k))) + &
                    ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31(i,j,k))) + &
                    ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32(i,j,k))) + &
                    ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33(i,j,k))) + &
                    ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33(i,j,k))) + &
                    ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33(i,j,k)))
        !
      TddO        = ProjectPi3(1,1,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T11d(i,j,k))) + &
                    ProjectPi3(2,1,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T11d(i,j,k))) + &
                    ProjectPi3(3,1,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T11d(i,j,k))) + &
                    ProjectPi3(1,1,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T12d(i,j,k))) + &
                    ProjectPi3(2,1,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T12d(i,j,k))) + &
                    ProjectPi3(3,1,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T12d(i,j,k))) + &
                    ProjectPi3(1,1,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T13d(i,j,k))) + &
                    ProjectPi3(2,1,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T13d(i,j,k))) + &
                    ProjectPi3(3,1,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T13d(i,j,k))) + &
                    ProjectPi3(1,2,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T21d(i,j,k))) + &
                    ProjectPi3(2,2,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T21d(i,j,k))) + &
                    ProjectPi3(3,2,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T21d(i,j,k))) + &
                    ProjectPi3(1,2,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T22d(i,j,k))) + &
                    ProjectPi3(2,2,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T22d(i,j,k))) + &
                    ProjectPi3(3,2,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T22d(i,j,k))) + &
                    ProjectPi3(1,2,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T23d(i,j,k))) + &
                    ProjectPi3(2,2,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T23d(i,j,k))) + &
                    ProjectPi3(3,2,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T23d(i,j,k))) + &
                    ProjectPi3(1,3,1,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T31d(i,j,k))) + &
                    ProjectPi3(2,3,1,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T31d(i,j,k))) + &
                    ProjectPi3(3,3,1,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T31d(i,j,k))) + &
                    ProjectPi3(1,3,2,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T32d(i,j,k))) + &
                    ProjectPi3(2,3,2,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T32d(i,j,k))) + &
                    ProjectPi3(3,3,2,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T32d(i,j,k))) + &
                    ProjectPi3(1,3,3,kx,ky,kz) * dimag(u1c(i,j,k)*conjg(T33d(i,j,k))) + &
                    ProjectPi3(2,3,3,kx,ky,kz) * dimag(u2c(i,j,k)*conjg(T33d(i,j,k))) + &
                    ProjectPi3(3,3,3,kx,ky,kz) * dimag(u3c(i,j,k)*conjg(T33d(i,j,k)))
        !
      TdsthetaO   = ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                    ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                    ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
                    ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                    ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                    ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tstheta3(i,j,k))) + &
                    ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta1(i,j,k))) + &
                    ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta2(i,j,k))) + &
                    ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tstheta3(i,j,k)))
        !
      TddthetaO   = ProjectPi2(1,1,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                    ProjectPi2(1,2,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                    ProjectPi2(1,3,kx,ky,kz) * dreal(u1c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
                    ProjectPi2(2,1,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                    ProjectPi2(2,2,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                    ProjectPi2(2,3,kx,ky,kz) * dreal(u2c(i,j,k)*conjg(Tdtheta3(i,j,k))) + &
                    ProjectPi2(3,1,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta1(i,j,k))) + &
                    ProjectPi2(3,2,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta2(i,j,k))) + &
                    ProjectPi2(3,3,kx,ky,kz) * dreal(u3c(i,j,k)*conjg(Tdtheta3(i,j,k)))
      if(kOrdinal<= allkmax)then
        Tcount(kOrdinal)   = Tcount(kOrdinal)   + 1
        if((method == 1) .or. (method == 2))then
          Ts(kOrdinal)       = Ts(kOrdinal)       - TsO*kk*kk
          Tstheta(kOrdinal)  = Tstheta(kOrdinal)  + TsthetaO*kk*kk
          Td(kOrdinal)       = Td(kOrdinal)       - TdO*kk*kk
          Tddtheta(kOrdinal) = Tddtheta(kOrdinal) + TddthetaO*kk*kk
          Tdstheta(kOrdinal) = Tdstheta(kOrdinal) + TdsthetaO*kk*kk
          Tss(kOrdinal)      = Tss(kOrdinal)      - TssO*kk*kk
          Tdd(kOrdinal)      = Tdd(kOrdinal)      - TddO*kk*kk
          kn(kOrdinal)       = kn(kOrdinal)       + kk
        elseif(method == 3)then
          Ts(kOrdinal)       = Ts(kOrdinal)       - TsO
          Tstheta(kOrdinal)  = Tstheta(kOrdinal)  + TsthetaO
          Td(kOrdinal)       = Td(kOrdinal)       - TdO
          Tddtheta(kOrdinal) = Tddtheta(kOrdinal) + TddthetaO
          Tdstheta(kOrdinal) = Tdstheta(kOrdinal) + TdsthetaO
          Tss(kOrdinal)      = Tss(kOrdinal)      - TssO
          Tdd(kOrdinal)      = Tdd(kOrdinal)      - TddO
        endif
      end if
      Tsall = Tsall  + TsO
      Tssall= Tssall + TssO
      Tcall = Tcall  + TdO
      Tccall= Tccall + TddO
      k2Ts  = k2Ts   + kk**2 * TsO
      k2Tss = k2Tss  + kk**2 * TssO
      k2Tc  = k2Tc   + kk**2 * TdO
      k2Tcc = k2Tcc  + kk**2 * TddO
    enddo
    enddo
    enddo
    !
    do i=0,allkmax
      Ts(i) = psum(Ts(i))
      Tstheta(i) = psum(Tstheta(i))
      Td(i) = psum(Td(i))
      Tdstheta(i) = psum(Tdstheta(i))
      Tddtheta(i) = psum(Tddtheta(i))
      Tss(i) = psum(Tss(i))
      Tdd(i) = psum(Tdd(i))
      Tcount(i) = psum(Tcount(i))
      if((method == 1) .or. (method == 2))then
        kn(i) = psum(kn(i))
        if(Tcount(i) .ne. 0)then
          Ts(i) = Ts(i)/Tcount(i)*4*pi
          Tstheta(i) = Tstheta(i)/Tcount(i)*4*pi
          Td(i) = Td(i)/Tcount(i)*4*pi
          Tdstheta(i) = Tdstheta(i)/Tcount(i)*4*pi
          Tddtheta(i) = Tddtheta(i)/Tcount(i)*4*pi
          Tss(i) = Tss(i)/Tcount(i)*4*pi
          Tdd(i) = Tdd(i)/Tcount(i)*4*pi
          kn(i) = kn(i)/tcount(i)
        endif
      else
        kn(i) = real(i)
      endif
    end do
    !
    Tsall = psum(Tsall)
    Tcall = psum(Tcall)
    k2Ts = psum(k2Ts)
    k2Tc = psum(k2Tc)
    Tssall = psum(Tssall)
    Tccall = psum(Tccall)
    k2Tss = psum(k2Tss)
    k2Tcc = psum(k2Tcc)
    !
    if(mpirank==0)  print*, '** Summation & average'
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_a, &
                        firstline='nstep time k Ts Tss Tstheta Td Tdd Tdstheta Tddtheta')
      do i=0,allkmax
        if(Tcount(i)>1e-3) call listwrite(hand_a,kn(i),Ts(i),Tss(i),Tstheta(i),Td(i),Tdd(i),Tdstheta(i),Tddtheta(i))
      end do
      !
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/Tspec_aux'//stepname//'.dat'
      else
        outfilename = 'pp/Tspec_aux.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
            firstline='nstep time k2Ts k2Tc Tsall Tcall Tssall Tccall k2Tss k2Tcc')
      call listwrite(hand_b,k2Ts,k2Tc,Tsall,Tcall,Tssall,Tccall,k2Tss,k2Tcc)
      !
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_theta)
    call fftw_free(c_u1s)
    call fftw_free(c_u2s)
    call fftw_free(c_u3s)
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_T11s)
    call fftw_free(c_T12s)
    call fftw_free(c_T13s)
    call fftw_free(c_T21s)
    call fftw_free(c_T22s)
    call fftw_free(c_T23s)
    call fftw_free(c_T31s)
    call fftw_free(c_T32s)
    call fftw_free(c_T33s)
    call fftw_free(c_T11d)
    call fftw_free(c_T12d)
    call fftw_free(c_T13d)
    call fftw_free(c_T21d)
    call fftw_free(c_T22d)
    call fftw_free(c_T23d)
    call fftw_free(c_T31d)
    call fftw_free(c_T32d)
    call fftw_free(c_T33d)
    call fftw_free(c_T11)
    call fftw_free(c_T12)
    call fftw_free(c_T13)
    call fftw_free(c_T21)
    call fftw_free(c_T22)
    call fftw_free(c_T23)
    call fftw_free(c_T31)
    call fftw_free(c_T32)
    call fftw_free(c_T33)
    call fftw_free(c_Tstheta1)
    call fftw_free(c_Tstheta2)
    call fftw_free(c_Tstheta3)
    call fftw_free(c_Tdtheta1)
    call fftw_free(c_Tdtheta2)
    call fftw_free(c_Tdtheta3)
    call mpistop
    !
  end subroutine instanttriad3D
  !
  !
  !
  subroutine instantspectraskewness2D(thefilenumb)
    !
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,   only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only : vel, rho
    use hdf5io
    use solver,    only : refcal
    use utility,   only : listinit,listwrite
    use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer,intent(in) :: thefilenumb
    character(len=128) :: infilename
    character(len=4) :: stepname
    real(8) :: u1mean,u2mean,rhomean
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u11,u12,u21,u22
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer :: allkmax
    real(8) :: u11c,u11s,u22c,u22s,thetac,thetas,omegas,psi2theta,phi2theta,thetamax
    character(len=128) :: outfilename
    integer :: hand_a,hand_b
    character(len=1) :: modeio
    integer :: i,j,n
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u11, c_u12, c_u21, c_u22
    !
    call readinput
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    modeio='h'
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    !
    allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
    if(mpirank==0)  print *, "ia:",ia,",ja:",ja,"allkmax:",allkmax
    if(ka .ne. 0) stop 'Please use skewness3D'
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
    !
    !
    if (thefilenumb .ne. 0) then
      write(stepname,'(i4.4)')thefilenumb
      infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    else
      infilename='outdat/flowfield.'//modeio//'5'
    endif
    !
    !
    call h5io_init(filename=infilename,mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='time',var=time)
    call h5read(varname='nstep',var=nstep)
    !
    call h5io_end
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Calculate favre average
    u1mean = 0.0d0
    u2mean = 0.0d0
    rhomean = 0.0d0
    !
    do i=1,im
      do j=1,jm
        u1mean = u1mean + rho(i,j,0) * vel(i,j,0,1)
        u2mean = u2mean + rho(i,j,0) * vel(i,j,0,2)
        rhomean = rhomean + rho(i,j,0)
      enddo
    enddo
    rhomean = psum(rhomean) / (1.0d0*ia*ja)
    u1mean = psum(u1mean) / (1.d0*ia*ja*rhomean)
    u2mean = psum(u2mean) / (1.d0*ia*ja*rhomean)
    !
    if(mpirank==0) print *, 'u1mean=',u1mean, 'u2mean=',u2mean
    !
    !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw])
    !
    ! planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !!!! Convert to complex
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=CMPLX(vel(i,j,0,1)-u1mean,0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(vel(i,j,0,2)-u2mean,0.d0,C_INTPTR_T);
      !
    end do
    end do
    !
    !!!! Do 2d FFT
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
      u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
      !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
    !
    c_u11 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u11, u11, [imfftw,jmfftw])
    c_u21 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u21, u21, [imfftw,jmfftw])
    c_u12 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u12, u12, [imfftw,jmfftw])
    c_u22 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u22, u22, [imfftw,jmfftw])
    !
    !
    do j=1,jm
      do i=1,im
        u11(i,j) = u1spe(i,j) * CMPLX(0.d0,k1(i,j),C_INTPTR_T)
        u12(i,j) = u1spe(i,j) * CMPLX(0.d0,k2(i,j),C_INTPTR_T)
        u21(i,j) = u2spe(i,j) * CMPLX(0.d0,k1(i,j),C_INTPTR_T)
        u22(i,j) = u2spe(i,j) * CMPLX(0.d0,k2(i,j),C_INTPTR_T)
      end do
    end do
    !
    !!!! Do inverse FFT
    call fftw_mpi_execute_dft(backward_plan,u11,u11)
    call fftw_mpi_execute_dft(backward_plan,u21,u21)
    call fftw_mpi_execute_dft(backward_plan,u12,u12)
    call fftw_mpi_execute_dft(backward_plan,u22,u22)
    !
    u11c = 0.d0
    u11s = 0.d0
    u22c = 0.d0
    u22s = 0.d0
    !
    thetac = 0.d0
    thetas = 0.d0
    omegas = 0.d0
    psi2theta = 0.d0
    phi2theta = 0.d0
    thetamax = 0.d0
    !
    do j=1,jm
      do i=1,im
        u11c = u11c + real(u11(i,j))**3
        u11s = u11s + real(u11(i,j))**2
        u22c = u22c + real(u22(i,j))**3
        u22s = u22s + real(u22(i,j))**2
        !
        thetac = thetac + real(u11(i,j)+u22(i,j))**3
        thetas = thetas + real(u11(i,j)+u22(i,j))**2
        omegas = omegas + real(u21(i,j)-u12(i,j))**2
        thetamax = max(thetamax,abs(u11(i,j)+u22(i,j)))
        !
        psi2theta = psi2theta + real(u12(i,j)+u21(i,j))**2*real(u11(i,j)+u22(i,j))
        phi2theta = phi2theta + real(u11(i,j)-u22(i,j))**2*real(u11(i,j)+u22(i,j))
        !
      end do
    end do
    !
    u11c      = psum(u11c)     /(1.d0*ia*ja)
    u11s      = psum(u11s)     /(1.d0*ia*ja)
    u22c      = psum(u22c)     /(1.d0*ia*ja)
    u22s      = psum(u22s)     /(1.d0*ia*ja)
    thetac    = psum(thetac)   /(1.d0*ia*ja)
    thetas    = psum(thetas)   /(1.d0*ia*ja)
    omegas    = psum(omegas)   /(1.d0*ia*ja)
    psi2theta = psum(psi2theta)/(1.d0*ia*ja)
    phi2theta = psum(phi2theta)/(1.d0*ia*ja)
    thetamax  = pmax(thetamax)
    !
    if(mpirank == 0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SkewnessSpec'//stepname//'.dat'
      else
        outfilename = 'pp/SkewnessSpec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time u11c u11s s1 u22c u22s s2 ullc ulls skew thetamax')
      call listwrite(hand_b,u11c,u11s,u11c/sqrt(u11s**3),&
              u22c,u22s,u22c/sqrt(u22s**3),&
              u11c+u22c,u11s+u22s,(u11c+u22c)/2.d0/sqrt((u11s/2.d0+u22s/2.d0)**3),thetamax)
      print*,' <<< '//outfilename//'... done.'
      !
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SkewmomSpec'//stepname//'.dat'
      else
        outfilename = 'pp/SkewmomSpec.dat'
      endif
      !
      call listinit(filename=outfilename,handle=hand_b, &
                    firstline='nstep time thetac thetas omegas psi2theta phi2theta Sa Sb')
      call listwrite(hand_b,thetac,thetas,omegas,psi2theta,phi2theta,&
                    5.d0/16.d0*thetac/sqrt((omegas/8.d0+3.d0*thetas/8.d0)**3),&
                    (thetac/8.d0+3.d0*psi2theta/16.d0+3.d0*phi2theta/16.d0)/sqrt((omegas/8.d0+3.d0*thetas/8.d0)**3))
      print*,' <<< '//outfilename//'... done.'
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u11)
    call fftw_free(c_u12)
    call fftw_free(c_u21)
    call fftw_free(c_u22)
    call mpistop
    !
    deallocate(k1,k2)
    !
  end subroutine instantspectraskewness2D
  !
  !
end module udf_pp_spectra