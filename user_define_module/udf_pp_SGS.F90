!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| physical SGS method for energy flux calculation.                    |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_SGS
    !
    !
    use constdef
    use stlaio,  only: get_unit
    !
    implicit none
    !
    contains
    !
    subroutine ppSGSentrance
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
            print*,' ** ppSGS readmode command: ',readmode
        endif
        !
        call bcast(readmode)
        !
        if(trim(readmode)=='Pi2Dint') then
        ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dint"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dint(filenumb)
        !
        elseif(trim(readmode)=='Pi2Dtot') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dtot"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dtot(filenumb)
            !
        elseif(trim(readmode)=='PiOmega2D') then
            ! 
            if(mpirank == 0) then
              print* ," ** Use SGSPiOmega2D"
              call readkeyboad(inputfile) 
              read(inputfile,'(i4)') filenumb
              print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPiOmega2D(filenumb)
            !
        elseif(trim(readmode)=='Pi2Dlocal') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dlocal"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dlocal(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dint') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dint"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dint(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dtot') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dtot"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dtot(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dlocal') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dlocal"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dlocal(filenumb)
            !
        elseif(trim(readmode)=='LES3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSLES3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSLES3D(filenumb)
            !
        elseif(trim(readmode)=='T3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGST3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGST3D(filenumb)
            !
        elseif(trim(readmode)=='stress2D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSstress2D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSstress2D(filenumb)
            !
        elseif(trim(readmode)=='stress3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSstress3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSstress3D(filenumb)
            !
        else
            print* ,"Readmode is not defined!", readmode
        endif
    end subroutine ppSGSentrance
    !
    !
  subroutine SGSPiOmega2D(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka,reynolds
      use commarray, only: vel,tmp,rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use fludyna,   only : miucal
      use solver,    only : refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1,u2,w,wu1,wu2
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      integer :: allkmax
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl, beta,roav,miu,miuav,miudrho
      real(8), allocatable, dimension(:) :: Pi_omega
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1_filted,u2_filted,w_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wu1_filted,wu2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1_filted,wx2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1,wx2
      !
      !
      type(C_PTR) :: c_u1,c_u2,c_w,forward_plan,backward_plan
      type(C_PTR) :: c_wx1,c_wx2,c_wu1,c_wu2
      type(C_PTR) :: c_u1_filted,c_u2_filted,c_w_filted
      type(C_PTR) :: c_wu1_filted,c_wu2_filted
      type(C_PTR) :: c_wx1_filted,c_wx2_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja, 'Reynolds:',reynolds
      allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2),rho(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),mode = modeio)
      call h5read(varname='t', var=tmp(0:im,0:jm,0:km),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
          print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_u1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u1, u1, [imfftw,jmfftw])
      c_u2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u2, u2, [imfftw,jmfftw])
      c_w = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w, w, [imfftw,jmfftw])
      c_wu1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu1, wu1, [imfftw,jmfftw])
      c_wu2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu2, wu2, [imfftw,jmfftw])
      c_wx1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx1, wx1, [imfftw,jmfftw])
      c_wx2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx2, wx2, [imfftw,jmfftw])
      !
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        u1(i,j)=CMPLX(vel(i,j,0,1),0.d0,C_INTPTR_T);
        u2(i,j)=CMPLX(vel(i,j,0,2),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, u1,u2,w are in spectral space
      call fftw_mpi_execute_dft(forward_plan,u1,u1)
      call fftw_mpi_execute_dft(forward_plan,u2,u2)
      !
      do j=1,jm
      do i=1,im
        !
        u1(i,j)=u1(i,j)/(1.d0*ia*ja)
        u2(i,j)=u2(i,j)/(1.d0*ia*ja)
        !
        w(i,j)=imag*k1(i,j)*u2(i,j)-imag*k2(i,j)*u1(i,j)
        !
        if(sqrt(k1(i,j)**2+k2(i,j)**2) > allkmax)then
          w(i,j) = 0.d0
        endif
        !
        wx1(i,j)=imag*k1(i,j)*w(i,j)
        wx2(i,j)=imag*k2(i,j)*w(i,j)
        !
      end do
      end do
      !
      !After this bloc,u1,u2,w,wx1,wx2 are in physical space
      call fftw_mpi_execute_dft(backward_plan,u1,u1)
      call fftw_mpi_execute_dft(backward_plan,u2,u2)
      call fftw_mpi_execute_dft(backward_plan,w,w)
      call fftw_mpi_execute_dft(backward_plan,wx1,wx1)
      call fftw_mpi_execute_dft(backward_plan,wx2,wx2)
      !
      beta = 0.d0
      roav = 0.d0
      miu = 0.d0
      miuav = 0.d0
      miudrho = 0.d0
      !!! 
      do j=1,jm
      do i=1,im
        !
        !
        wu1(i,j)=w(i,j)*u1(i,j)
        wu2(i,j)=w(i,j)*u2(i,j)
        !
        !
        miu=miucal(tmp(i,j,0))/reynolds
        beta = beta + miu/rho(i,j,0) * (dreal(wx1(i,j))**2 + dreal(wx2(i,j))**2)
        !
        roav=roav+rho(i,j,0)
        miuav=miuav+miu
        miudrho = miudrho+miu/rho(i,j,0)
        !
      end do
      end do
      !
      beta  = psum(beta) / (ia*ja)
      roav  = psum(roav) / (ia*ja)
      miuav = psum(miuav)/ (ia*ja)
      miudrho = psum(miudrho) / (ia*ja)
      !
      !After this bloc, u1,u2,w,wu1,wu2 are in spectral space
      call fftw_mpi_execute_dft(forward_plan,u1,u1)
      call fftw_mpi_execute_dft(forward_plan,u2,u2)
      call fftw_mpi_execute_dft(forward_plan,w,w)
      call fftw_mpi_execute_dft(forward_plan,wu1,wu1)
      call fftw_mpi_execute_dft(forward_plan,wu2,wu2)
      !
      do j=1,jm
      do i=1,im
        !
        u1(i,j)=u1(i,j)/(1.d0*ia*ja)
        u2(i,j)=u2(i,j)/(1.d0*ia*ja)
        w(i,j)=w(i,j)/(1.d0*ia*ja)
        wu1(i,j)=wu1(i,j)/(1.d0*ia*ja)
        wu2(i,j)=wu2(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_omega(1:num_l))
      !
      c_u1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u1_filted, u1_filted, [imfftw,jmfftw])
      c_u2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u2_filted, u2_filted, [imfftw,jmfftw])
      c_w_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w_filted, w_filted, [imfftw,jmfftw])
      c_wu1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu1_filted, wu1_filted, [imfftw,jmfftw])
      c_wu2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu2_filted, wu2_filted, [imfftw,jmfftw])
      c_wx1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx1_filted, wx1_filted, [imfftw,jmfftw])
      c_wx2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx2_filted, wx2_filted, [imfftw,jmfftw])
      !
      Pi_omega = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, u1_filted is in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
          !
          u1_filted(i,j) = u1(i,j) *Gl
          u2_filted(i,j) = u2(i,j) *Gl
          w_filted(i,j)  = w(i,j)  *Gl
          wu1_filted(i,j)= wu1(i,j)*Gl
          wu2_filted(i,j)= wu2(i,j)*Gl
          wx1_filted(i,j)= w(i,j)  *Gl*imag*k1(i,j)
          wx2_filted(i,j)= w(i,j)  *Gl*imag*k2(i,j)
          !
        enddo
        enddo
        !
        ! After this bloc, u1_filted is in physical space
        call fftw_mpi_execute_dft(backward_plan,u1_filted,u1_filted)
        call fftw_mpi_execute_dft(backward_plan,u2_filted,u2_filted)
        call fftw_mpi_execute_dft(backward_plan,w_filted,w_filted)
        call fftw_mpi_execute_dft(backward_plan,wu1_filted,wu1_filted)
        call fftw_mpi_execute_dft(backward_plan,wu2_filted,wu2_filted)
        call fftw_mpi_execute_dft(backward_plan,wx1_filted,wx1_filted)
        call fftw_mpi_execute_dft(backward_plan,wx2_filted,wx2_filted)
        !
        ! 
        do j=1,jm
        do i=1,im
          !
          Pi_omega(m) = Pi_omega(m) + &
              dreal(wu1_filted(i,j) - w_filted(i,j) * u1_filted(i,j))*dreal(wx1_filted(i,j)) + &
              dreal(wu2_filted(i,j) - w_filted(i,j) * u2_filted(i,j))*dreal(wx2_filted(i,j)) 
          !
        enddo
        enddo
        !
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_omega(m) =	 psum(Pi_omega(m)) / (ia*ja)
        !
        !
      enddo
      !
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Piomega_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Piomega.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell piomega beta miudrho miudrho2 lens')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_omega(m), beta,miuav/roav,miudrho,((miudrho)**3/beta)**(1.d0/6.d0))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_u1)
      call fftw_free(c_u2)
      call fftw_free(c_w)
      call fftw_free(c_wx1)
      call fftw_free(c_wx2)
      call fftw_free(c_wu1)
      call fftw_free(c_wu2)
      call fftw_free(c_u1_filted)
      call fftw_free(c_u2_filted)
      call fftw_free(c_w_filted)
      call fftw_free(c_wu1_filted)
      call fftw_free(c_wu2_filted)
      call fftw_free(c_wx1_filted)
      call fftw_free(c_wx2_filted)
      call mpistop
      deallocate(k1,k2)
      deallocate(l_lim)
      deallocate(Pi_omega)
      !
    end subroutine SGSPiOmega2D
    !
    subroutine SGSPi2Dtot(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pi_tot
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
      real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
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
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
        w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
        w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
        w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
        !
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
          print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      allocate(tau(1:2,1:2,1:im,1:jm))
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_tot(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
      !
      Pi_tot = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
          !
          w1_filted(i,j)    = w1(i,j)   *Gl
          w2_filted(i,j)    = w2(i,j)   *Gl
          !
          w1w1_filted(i,j)  = w1w1(i,j) *Gl
          w1w2_filted(i,j)  = w1w2(i,j) *Gl
          w2w1_filted(i,j)  = w2w1(i,j) *Gl
          w2w2_filted(i,j)  = w2w2(i,j) *Gl
          !
          rho_filted(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          !
          tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
          tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
          tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
          tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
          !
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        !
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          !
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        !
        do j=1,jm
        do i=1,im
          !
          Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j) * A11_filted(i,j) + &
                                  tau(1,2,i,j) * A12_filted(i,j) + &
                                  tau(2,1,i,j) * A21_filted(i,j) + &
                                  tau(2,2,i,j) * A22_filted(i,j)
          !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pitot.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pitot')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_tot(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call mpistop
      deallocate(k1,k2,tau)
      deallocate(l_lim)
      deallocate(Pi_tot)
      !
    end subroutine SGSPi2Dtot
    !
    subroutine SGSPi2Dlocal(thefilenumb)
      !
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1f,w2f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
      !
      complex(8), allocatable, dimension(:,:) :: All
      complex(8), allocatable, dimension(:,:) :: S11,S12,S21,S22
      complex(8), allocatable, dimension(:,:) :: W12,W21
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
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
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm),&
              S11(1:im,1:jm),S12(1:im,1:jm),&
              S21(1:im,1:jm),S22(1:im,1:jm),&
              W12(1:im,1:jm),W21(1:im,1:jm))
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
          print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw])
      !
      Pis1=0.d0
      Pis2=0.d0
      Pim2=0.d0
      Pim3=0.d0
      Pid=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1 is (rho*u1) in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1f(i,j)    = w1(i,j)    *Gl
          w2f(i,j)    = w2(i,j)    *Gl
          !
          rhof(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1 is (rho*u1) in physical space
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        ! After this bloc, w1 is u1 in physical space
        do j=1,jm
        do i=1,im
          w1f(i,j) = w1f(i,j)/rhof(i,j)
          w2f(i,j) = w2f(i,j)/rhof(i,j)
          !
          !
        enddo
        enddo
        !
        ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        !
        do j=1,jm
        do i=1,im
          !
          w1f(i,j)  = w1f(i,j)/(1.d0*ia*ja)
          w2f(i,j)  = w2f(i,j)/(1.d0*ia*ja)
          !
          A11(i,j) = imag*w1f(i,j)*k1(i,j)
          A21(i,j) = imag*w2f(i,j)*k1(i,j)
          A12(i,j) = imag*w1f(i,j)*k2(i,j)
          A22(i,j) = imag*w2f(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        !
        ! After this bloc, A11 is A11 in physical space
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        !
        do j=1,jm
        do i=1,im
          !
          All(i,j) = A11(i,j)+A22(i,j)
          !
          S11(i,j) = A11(i,j) - 1.d0/2.d0 * All(i,j)
          S22(i,j) = A22(i,j) - 1.d0/2.d0 * All(i,j)
          S12(i,j) = (A12(i,j) + A21(i,j))*0.5d0
          S21(i,j) = S12(i,j)
          !
          W12(i,j) = (A12(i,j)-A21(i,j))*0.5d0
          W21(i,j) = -1.d0*W12(i,j)
          !
        end do
        end do
        !
        do j=1,jm
        do i=1,im
          !
          Pis1(m) = Pis1(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*S11(i,j)+ &
                                          S12(i,j)*S12(i,j)*S11(i,j)+ &
                                          S11(i,j)*S21(i,j)*S12(i,j)+ &
                                          S12(i,j)*S22(i,j)*S12(i,j)+ &
                                          S21(i,j)*S11(i,j)*S21(i,j)+ &
                                          S22(i,j)*S12(i,j)*S21(i,j)+ &
                                          S21(i,j)*S21(i,j)*S22(i,j)+ &
                                          S22(i,j)*S22(i,j)*S22(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pis2(m) = Pis2(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*S11(i,j)+ &
                                            W21(i,j)*W12(i,j)*S22(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pim2(m) = Pim2(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*All(i,j)+ &
                                            S12(i,j)*S12(i,j)*All(i,j)+ &
                                            S21(i,j)*S21(i,j)*All(i,j)+ &
                                            S22(i,j)*S22(i,j)*All(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pim3(m) = Pim3(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*All(i,j)+ &
                                            W21(i,j)*W12(i,j)*All(i,j)) * &
                                        l_lim(m) * l_lim(m)
          Pid(m) = Pid(m) + rhof(i,j) * All(i,j) * All(i,j) * All(i,j) * &
                                      l_lim(m) * l_lim(m)
          !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pis1(m) =	 psum(Pis1(m)) / (ia*ja)
        Pis2(m) =	 - psum(Pis2(m)) / (ia*ja)
        Pim2(m) =	 psum(Pim2(m)) / (ia*ja) /2
        Pim3(m) =	 - psum(Pim3(m)) / (ia*ja) /2
        Pid(m) =	 psum(Pid(m)) / (ia*ja) /4
        !
        !
      enddo
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pilocal.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pis1 pis2 pim2 pim3 pid')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call mpistop
      deallocate(k1,k2)
      deallocate(l_lim)
      deallocate(All,S11,S12,S21,S22,W12,W21)
      deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
      !
    end subroutine SGSPi2Dlocal
    !
    subroutine SGSPi2Dint(thefilenumb)
      !
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
      complex(8), allocatable, dimension(:,:) :: All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted
      complex(8), allocatable, dimension(:,:) :: All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: term1_11,term1_12,term1_21,term1_22,term2, &
                                                            term3_11,term3_12,term3_21,term3_22, &
                                                            term4_11,term4_22,term5,term6_11,term6_12, &
                                                            term6_21,term6_22,term7
      real(8) :: vxr_D1,vxr_D2,vxr_D3,vxr_D4,vxr_D5,vxr_D6,vxr_D7
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
      type(C_PTR) :: c_term1_11,c_term1_12,c_term1_21,c_term1_22,c_term2,c_term3_11,c_term3_12,c_term3_21,c_term3_22, &
                                                  c_term4_11,c_term4_22,c_term5,c_term6_11,c_term6_12,c_term6_21,c_term6_22,c_term7
      real(8), allocatable, dimension(:) :: Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
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
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
          print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!! allocation
      allocate(Pi1(1:num_l), Pi2(1:num_l), Pi3(1:num_l), Pi4(1:num_l), Pi5(1:num_l), Pi6(1:num_l), Pi7(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
      !
      allocate(All_filted(1:im,1:jm),S11_filted(1:im,1:jm),S12_filted(1:im,1:jm),S21_filted(1:im,1:jm),S22_filted(1:im,1:jm), &
                W12_filted(1:im,1:jm),W21_filted(1:im,1:jm))
      !
      allocate(All_filted_l(1:im,1:jm),S11_filted_l(1:im,1:jm),S12_filted_l(1:im,1:jm),S21_filted_l(1:im,1:jm),  &
              S22_filted_l(1:im,1:jm))
      !
      c_term1_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_11, term1_11, [imfftw,jmfftw])
      c_term1_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_12, term1_12, [imfftw,jmfftw])
      c_term1_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_21, term1_21, [imfftw,jmfftw])
      c_term1_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_22, term1_22, [imfftw,jmfftw])
      c_term2    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term2, term2, [imfftw,jmfftw])
      c_term3_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_11, term3_11, [imfftw,jmfftw])
      c_term3_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_12, term3_12, [imfftw,jmfftw])
      c_term3_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_21, term3_21, [imfftw,jmfftw])
      c_term3_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_22, term3_22, [imfftw,jmfftw])
      c_term4_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_11, term4_11, [imfftw,jmfftw])
      c_term4_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_22, term4_22, [imfftw,jmfftw])
      c_term5    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term5, term5, [imfftw,jmfftw])
      c_term6_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_11, term6_11, [imfftw,jmfftw])
      c_term6_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_12, term6_12, [imfftw,jmfftw])
      c_term6_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_21, term6_21, [imfftw,jmfftw])
      c_term6_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_22, term6_22, [imfftw,jmfftw])
      c_term7    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term7, term7, [imfftw,jmfftw])
      !
      !
      Pi1 = 0.d0
      Pi2 =	0.d0
      Pi3 =	0.d0
      Pi4 =	0.d0
      Pi5 =	0.d0
      Pi6 =	0.d0
      Pi7 =	0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename2,handle=hand_b, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        endif
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :1
          w1_filted(i,j) = w1(i,j)     *Gl
          w2_filted(i,j) = w2(i,j)     *Gl
          rho_filted(i,j) = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        !
        do j=1,jm
        do i=1,im
          All_filted_l(i,j) = A11_filted(i,j)+A22_filted(i,j)
        !
          S11_filted_l(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
          S12_filted_l(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
          S21_filted_l(i,j) = S12_filted_l(i,j)
          S22_filted_l(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
        !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !!!! Velocity Favre average and density average
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          do i=1,im
          do j=1,jm
            Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j)  = w1(i,j)    *Galpha
            w2_filted(i,j)  = w2(i,j)    *Galpha
            rho_filted(i,j) = rhocom(i,j)*Galpha
          enddo
          enddo
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          do i=1,im
          do j=1,jm
            w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
            w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          enddo
          enddo
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
            w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
            A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
            A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
            A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
            A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
            !
          end do
          end do
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
          call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
          call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
          call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
          !
          !
          do j=1,jm
          do i=1,im
            All_filted(i,j) = A11_filted(i,j)+A22_filted(i,j)
            ! 
            S11_filted(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
            S12_filted(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
            S21_filted(i,j) = S12_filted(i,j)
            S22_filted(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
            !
            W12_filted(i,j) = (A12_filted(i,j)-A21_filted(i,j))*0.5d0
            W21_filted(i,j) = -1.d0*W12_filted(i,j)
          enddo
          enddo
          !
          !!!! Pi terms
          !
          do j=1,jm
          do i=1,im
            !term1_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted
            term1_11(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)
            term1_21(i,j) = rho_filted(i,j)*S21_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S12_filted(i,j)
            term1_12(i,j) = rho_filted(i,j)*S11_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S22_filted(i,j)
            term1_22(i,j) = rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
            ! 
            !
            term2(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j)+rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)+ &
                    rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j)+rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
            !
            ! term3_IJ = rho_filted*All_filted*SIJ_filted
            term3_11(i,j) = rho_filted(i,j)*All_filted(i,j)*S11_filted(i,j)
            term3_12(i,j) = rho_filted(i,j)*All_filted(i,j)*S12_filted(i,j)
            term3_21(i,j) = rho_filted(i,j)*All_filted(i,j)*S21_filted(i,j)
            term3_22(i,j) = rho_filted(i,j)*All_filted(i,j)*S22_filted(i,j)
            !
            !term4_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted
            term4_11(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j)
            term4_22(i,j) = rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
            !
            term5(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j) + rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
            !
            !term6_IJ= rho_filted*(S1J_filted*WI1_filted-SI1_filted*W1J_filted) + rho_filted*(S2J_filted*WI2_filted-SI2_filted*W2J_filted)
            term6_11(i,j)=   rho_filted(i,j)*(S21_filted(i,j)*W12_filted(i,j)-S12_filted(i,j)*W21_filted(i,j))
            term6_12(i,j)= - rho_filted(i,j)*(S11_filted(i,j)*W12_filted(i,j)) + rho_filted(i,j)*(S22_filted(i,j)*W12_filted(i,j))
            term6_21(i,j)=   rho_filted(i,j)*(S11_filted(i,j)*W21_filted(i,j)) - rho_filted(i,j)*(S22_filted(i,j)*W21_filted(i,j))
            term6_22(i,j)=   rho_filted(i,j)*(S12_filted(i,j)*W21_filted(i,j)-S21_filted(i,j)*W12_filted(i,j))
            !
            term7(i,j) = rho_filted(i,j)*All_filted(i,j)*All_filted(i,j)
          enddo
          enddo
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(forward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(forward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(forward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
          call fftw_mpi_execute_dft(forward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(forward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(forward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(forward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(forward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(forward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(forward_plan,term5   ,term5   )
          call fftw_mpi_execute_dft(forward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(forward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(forward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(forward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(forward_plan,term7   ,term7   )
          !
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            term1_11(i,j) = term1_11(i,j)*Gphi/(1.d0*ia*ja)
            term1_12(i,j) = term1_12(i,j)*Gphi/(1.d0*ia*ja)
            term1_21(i,j) = term1_21(i,j)*Gphi/(1.d0*ia*ja)
            term1_22(i,j) = term1_22(i,j)*Gphi/(1.d0*ia*ja)
            term2(i,j)    = term2(i,j)   *Gphi/(1.d0*ia*ja)
            term3_11(i,j) = term3_11(i,j)*Gphi/(1.d0*ia*ja)
            term3_12(i,j) = term3_12(i,j)*Gphi/(1.d0*ia*ja)
            term3_21(i,j) = term3_21(i,j)*Gphi/(1.d0*ia*ja)
            term3_22(i,j) = term3_22(i,j)*Gphi/(1.d0*ia*ja)
            term4_11(i,j) = term4_11(i,j)*Gphi/(1.d0*ia*ja)
            term4_22(i,j) = term4_22(i,j)*Gphi/(1.d0*ia*ja)
            term5(i,j)    = term5(i,j)   *Gphi/(1.d0*ia*ja)
            term6_11(i,j) = term6_11(i,j)*Gphi/(1.d0*ia*ja)
            term6_12(i,j) = term6_12(i,j)*Gphi/(1.d0*ia*ja)
            term6_21(i,j) = term6_21(i,j)*Gphi/(1.d0*ia*ja)
            term6_22(i,j) = term6_22(i,j)*Gphi/(1.d0*ia*ja)
            term7(i,j)    = term7(i,j)   *Gphi/(1.d0*ia*ja)
            !
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(backward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(backward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(backward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(backward_plan,term2   ,term2   )
          call fftw_mpi_execute_dft(backward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(backward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(backward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(backward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(backward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(backward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(backward_plan,term5   ,term5   )
          call fftw_mpi_execute_dft(backward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(backward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(backward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(backward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(backward_plan,term7   ,term7   )
          !
          !
          do j=1,jm
          do i=1,im
            vxr_D1 = term1_11(i,j) * S11_filted_l(i,j) + term1_12(i,j) * S12_filted_l(i,j) + &
                    term1_21(i,j) * S21_filted_l(i,j) + term1_22(i,j) * S22_filted_l(i,j)
            Pi1(m) = Pi1(m) + vxr_D1 * dl_alpha(m,n)
            !
            vxr_D2 = term2(i,j) * All_filted_l(i,j)
            Pi2(m) = Pi2(m) + vxr_D2 * dl_alpha(m,n) * 0.5d0 
            !
            vxr_D3 = term3_11(i,j) * S11_filted_l(i,j) + term3_12(i,j) * S12_filted_l(i,j) + &
                    term3_21(i,j) * S21_filted_l(i,j) + term3_22(i,j) * S22_filted_l(i,j)
            Pi3(m) = Pi3(m) + vxr_D3 * dl_alpha(m,n)
            !
            vxr_D4 = term4_11(i,j) * S11_filted_l(i,j) + term4_22(i,j) * S22_filted_l(i,j)
            Pi4(m) = Pi4(m) - vxr_D4 * dl_alpha(m,n)
            !
            vxr_D5 = term5(i,j) * All_filted_l(i,j)
            Pi5(m) = Pi5(m) - vxr_D5 * dl_alpha(m,n) * 0.5d0
            !
            vxr_D6 = term6_11(i,j) * S11_filted_l(i,j) + term6_12(i,j) * S12_filted_l(i,j) + &
                    term6_21(i,j) * S21_filted_l(i,j) + term6_22(i,j) * S22_filted_l(i,j)
            Pi6(m) = Pi6(m) + vxr_D6 * dl_alpha(m,n)
            !
            vxr_D7 = term7(i,j) * All_filted_l(i,j)
            Pi7(m) = Pi7(m) + vxr_D7 * dl_alpha(m,n) * 0.25d0
          enddo
          enddo
          !
          vxr_D1 = psum(vxr_D1)
          vxr_D2 = psum(vxr_D2)
          vxr_D3 = psum(vxr_D3)
          vxr_D4 = psum(vxr_D4)
          vxr_D5 = psum(vxr_D5)
          vxr_D6 = psum(vxr_D6)
          vxr_D7 = psum(vxr_D7)
          !
          if(mpirank==0) then
            call listwrite(hand_b,l_sqrtalpha(m,n),vxr_D1 * dl_alpha(m,n), &
                          vxr_D2 * dl_alpha(m,n) / 2.d0, &
                          vxr_D3 * dl_alpha(m,n) , &
                          - vxr_D4 * dl_alpha(m,n), &
                          - vxr_D5 * dl_alpha(m,n) * 0.5d0, &
                          vxr_D6 * dl_alpha(m,n), &
                          vxr_D7 * dl_alpha(m,n) * 0.25d0)
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        Pi1(m) =	 psum(Pi1(m)) / (ia*ja)
        Pi2(m) =	 psum(Pi2(m)) / (ia*ja)
        Pi3(m) =	 psum(Pi3(m)) / (ia*ja)
        Pi4(m) =	 psum(Pi4(m)) / (ia*ja)
        Pi5(m) =	 psum(Pi5(m)) / (ia*ja)
        Pi6(m) =	 psum(Pi6(m)) / (ia*ja)
        Pi7(m) =	 psum(Pi7(m)) / (ia*ja)
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pi.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi1(m), Pi2(m), &
                          Pi3(m), Pi4(m), Pi5(m),     &
                          Pi6(m), Pi7(m))
        end do
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_term1_11)
      call fftw_free(c_term1_12)
      call fftw_free(c_term1_21)
      call fftw_free(c_term1_22)
      call fftw_free(c_term2)
      call fftw_free(c_term3_11)
      call fftw_free(c_term3_12)
      call fftw_free(c_term3_21)
      call fftw_free(c_term3_22)
      call fftw_free(c_term4_11)
      call fftw_free(c_term4_22)
      call fftw_free(c_term5)
      call fftw_free(c_term6_11)
      call fftw_free(c_term6_12)
      call fftw_free(c_term6_21)
      call fftw_free(c_term6_22)
      call fftw_free(c_term7)
      call mpistop
      deallocate(All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted)
      deallocate(All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l)
      deallocate(k1,k2)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7)
      !
    end subroutine SGSPi2Dint
    !
    subroutine SGSPi3Dtot(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pi_tot
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
      type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
      type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
      type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
      type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
      c_w1w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
      c_w2w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
      c_w3w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
      c_w3w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
      c_w3w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
      c_w1w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
      c_w2w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
      c_w3w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
      c_w3w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
      call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
      call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
      call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
        w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
        w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
        w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
        w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
        w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
        w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
        w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
        w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      allocate(tau(1:3,1:3,1:im,1:jm,1:km))
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_tot(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      Pi_tot = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          !
          w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
          w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
          w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
          w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
          w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
          w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
          w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
          w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
          w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
          !
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w1_filted(i,j,k))
          tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w2_filted(i,j,k))
          tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w3_filted(i,j,k))
          tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w1_filted(i,j,k))
          tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w2_filted(i,j,k))
          tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w3_filted(i,j,k))
          tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w1_filted(i,j,k))
          tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w2_filted(i,j,k))
          tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w3_filted(i,j,k))
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j,k) * A11_filted(i,j,k) + &
                                  tau(1,2,i,j,k) * A12_filted(i,j,k) + &
                                  tau(1,3,i,j,k) * A13_filted(i,j,k) + &
                                  tau(2,1,i,j,k) * A21_filted(i,j,k) + &
                                  tau(2,2,i,j,k) * A22_filted(i,j,k) + &
                                  tau(2,3,i,j,k) * A23_filted(i,j,k) + &
                                  tau(3,1,i,j,k) * A31_filted(i,j,k) + &
                                  tau(3,2,i,j,k) * A32_filted(i,j,k) + &
                                  tau(3,3,i,j,k) * A33_filted(i,j,k)
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja*ka)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pitot.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pitot')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_tot(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w1w3)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w2w3)
      call fftw_free(c_w3w1)
      call fftw_free(c_w3w2)
      call fftw_free(c_w3w3)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w1w3_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_w2w3_filted)
      call fftw_free(c_w3w1_filted)
      call fftw_free(c_w3w2_filted)
      call fftw_free(c_w3w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call mpistop
      deallocate(k1,k2,k3,tau)
      deallocate(l_lim)
      deallocate(Pi_tot)
      !
    end subroutine SGSPi3Dtot
    !
    subroutine SGSPi3Dlocal(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1f,w2f,w3f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
      !
      complex(8), allocatable, dimension(:,:,:) :: All
      complex(8), allocatable, dimension(:,:,:) :: S11,S12,S13
      complex(8), allocatable, dimension(:,:,:) :: S21,S22,S23
      complex(8), allocatable, dimension(:,:,:) :: S31,S32,S33
      complex(8), allocatable, dimension(:,:,:) :: W12,W21
      complex(8), allocatable, dimension(:,:,:) :: W13,W31
      complex(8), allocatable, dimension(:,:,:) :: W23,W32
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_w3f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A13
      type(C_PTR) :: c_A21,c_A22,c_A23
      type(C_PTR) :: c_A31,c_A32,c_A33
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm,1:km),&
              S11(1:im,1:jm,1:km),S12(1:im,1:jm,1:km),S13(1:im,1:jm,1:km),&
              S21(1:im,1:jm,1:km),S22(1:im,1:jm,1:km),S23(1:im,1:jm,1:km),&
              S31(1:im,1:jm,1:km),S32(1:im,1:jm,1:km),S33(1:im,1:jm,1:km),&
              W12(1:im,1:jm,1:km),W21(1:im,1:jm,1:km),&
              W13(1:im,1:jm,1:km),W31(1:im,1:jm,1:km),&
              W23(1:im,1:jm,1:km),W32(1:im,1:jm,1:km))
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw,kmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw,kmfftw])
      c_w3f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3f, w3f,  [imfftw,jmfftw,kmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw,kmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw,kmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw,kmfftw])
      c_A13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13, A13,[imfftw,jmfftw,kmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw,kmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw,kmfftw])
      c_A23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23, A23,[imfftw,jmfftw,kmfftw])
      c_A31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31, A31,[imfftw,jmfftw,kmfftw])
      c_A32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32, A32,[imfftw,jmfftw,kmfftw])
      c_A33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33, A33,[imfftw,jmfftw,kmfftw])
      !
      Pis1=0.d0
      Pis2=0.d0
      Pim2=0.d0
      Pim3=0.d0
      Pid=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1 is (rho*u1) in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1f(i,j,k)    = w1(i,j,k)    *Gl
          w2f(i,j,k)    = w2(i,j,k)    *Gl
          w3f(i,j,k)    = w3(i,j,k)    *Gl
          !
          rhof(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1 is (rho*u1) in physical space
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,w3f,w3f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        ! After this bloc, w1 is u1 in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1f(i,j,k) = w1f(i,j,k)/rhof(i,j,k)
          w2f(i,j,k) = w2f(i,j,k)/rhof(i,j,k)
          w3f(i,j,k) = w3f(i,j,k)/rhof(i,j,k)
          !
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(forward_plan,w3f,w3f)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1f(i,j,k)  = w1f(i,j,k)/(1.d0*ia*ja*ka)
          w2f(i,j,k)  = w2f(i,j,k)/(1.d0*ia*ja*ka)
          w3f(i,j,k)  = w3f(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11(i,j,k) = imag*w1f(i,j,k)*k1(i,j,k)
          A21(i,j,k) = imag*w2f(i,j,k)*k1(i,j,k)
          A31(i,j,k) = imag*w3f(i,j,k)*k1(i,j,k)
          A12(i,j,k) = imag*w1f(i,j,k)*k2(i,j,k)
          A22(i,j,k) = imag*w2f(i,j,k)*k2(i,j,k)
          A32(i,j,k) = imag*w3f(i,j,k)*k2(i,j,k)
          A13(i,j,k) = imag*w1f(i,j,k)*k3(i,j,k)
          A23(i,j,k) = imag*w2f(i,j,k)*k3(i,j,k)
          A33(i,j,k) = imag*w3f(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11 is A11 in physical space
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A31,A31)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        call fftw_mpi_execute_dft(backward_plan,A32,A32)
        call fftw_mpi_execute_dft(backward_plan,A13,A13)
        call fftw_mpi_execute_dft(backward_plan,A23,A23)
        call fftw_mpi_execute_dft(backward_plan,A33,A33)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All(i,j,k) = A11(i,j,k)+A22(i,j,k)+A33(i,j,k)
          !
          S11(i,j,k) = A11(i,j,k) - 1.d0/3.d0 * All(i,j,k)
          S22(i,j,k) = A22(i,j,k) - 1.d0/3.d0 * All(i,j,k)
          S33(i,j,k) = A33(i,j,k) - 1.d0/3.d0 * All(i,j,k)
          S12(i,j,k) = (A12(i,j,k) + A21(i,j,k))*0.5d0
          S21(i,j,k) = S12(i,j,k)
          S13(i,j,k) = (A13(i,j,k) + A31(i,j,k))*0.5d0
          S31(i,j,k) = S13(i,j,k)
          S23(i,j,k) = (A23(i,j,k) + A32(i,j,k))*0.5d0
          S32(i,j,k) = S23(i,j,k)
          !
          W12(i,j,k) = (A12(i,j,k)-A21(i,j,k))*0.5d0
          W21(i,j,k) = -1.d0*W12(i,j,k)
          W13(i,j,k) = (A13(i,j,k)-A31(i,j,k))*0.5d0
          W31(i,j,k) = -1.d0*W13(i,j,k)
          W23(i,j,k) = (A23(i,j,k)-A32(i,j,k))*0.5d0
          W32(i,j,k) = -1.d0*W23(i,j,k)
          !
        end do
        end do
        end do
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Pis1(m) = Pis1(m) + rhof(i,j,k) *(S11(i,j,k)*S11(i,j,k)*S11(i,j,k)+ &
                                          S12(i,j,k)*S12(i,j,k)*S11(i,j,k)+ &
                                          S13(i,j,k)*S13(i,j,k)*S11(i,j,k)+ &
                                          S11(i,j,k)*S21(i,j,k)*S12(i,j,k)+ &
                                          S12(i,j,k)*S22(i,j,k)*S12(i,j,k)+ &
                                          S13(i,j,k)*S23(i,j,k)*S12(i,j,k)+ &
                                          S11(i,j,k)*S31(i,j,k)*S13(i,j,k)+ &
                                          S12(i,j,k)*S32(i,j,k)*S13(i,j,k)+ &
                                          S13(i,j,k)*S33(i,j,k)*S13(i,j,k)+ &
                                          S21(i,j,k)*S11(i,j,k)*S21(i,j,k)+ &
                                          S22(i,j,k)*S12(i,j,k)*S21(i,j,k)+ &
                                          S23(i,j,k)*S13(i,j,k)*S21(i,j,k)+ &
                                          S21(i,j,k)*S21(i,j,k)*S22(i,j,k)+ &
                                          S22(i,j,k)*S22(i,j,k)*S22(i,j,k)+ &
                                          S23(i,j,k)*S23(i,j,k)*S22(i,j,k)+ &
                                          S21(i,j,k)*S31(i,j,k)*S23(i,j,k)+ &
                                          S22(i,j,k)*S32(i,j,k)*S23(i,j,k)+ &
                                          S23(i,j,k)*S33(i,j,k)*S23(i,j,k)+ &
                                          S31(i,j,k)*S11(i,j,k)*S31(i,j,k)+ &
                                          S32(i,j,k)*S12(i,j,k)*S31(i,j,k)+ &
                                          S33(i,j,k)*S13(i,j,k)*S31(i,j,k)+ &
                                          S31(i,j,k)*S21(i,j,k)*S32(i,j,k)+ &
                                          S32(i,j,k)*S22(i,j,k)*S32(i,j,k)+ &
                                          S33(i,j,k)*S23(i,j,k)*S32(i,j,k)+ &
                                          S31(i,j,k)*S31(i,j,k)*S33(i,j,k)+ &
                                          S32(i,j,k)*S32(i,j,k)*S33(i,j,k)+ &
                                          S33(i,j,k)*S33(i,j,k)*S33(i,j,k)) * &
                                          l_lim(m) * l_lim(m)
          Pis2(m) = Pis2(m) + rhof(i,j,k) *(W12(i,j,k)*W21(i,j,k)*S11(i,j,k)+ &
                                            W13(i,j,k)*W31(i,j,k)*S11(i,j,k)+ &
                                            W13(i,j,k)*W32(i,j,k)*S12(i,j,k)+ &
                                            W12(i,j,k)*W23(i,j,k)*S13(i,j,k)+ &
                                            W23(i,j,k)*W31(i,j,k)*S21(i,j,k)+ &
                                            W21(i,j,k)*W12(i,j,k)*S22(i,j,k)+ &
                                            W23(i,j,k)*W32(i,j,k)*S22(i,j,k)+ &
                                            W21(i,j,k)*W13(i,j,k)*S23(i,j,k)+ &
                                            W32(i,j,k)*W21(i,j,k)*S31(i,j,k)+ &
                                            W31(i,j,k)*W12(i,j,k)*S32(i,j,k)+ &
                                            W31(i,j,k)*W13(i,j,k)*S33(i,j,k)+ &
                                            W32(i,j,k)*W23(i,j,k)*S33(i,j,k)) * &
                                          l_lim(m) * l_lim(m)
          Pim2(m) = Pim2(m) + rhof(i,j,k) *(S11(i,j,k)*S11(i,j,k)*All(i,j,k)+ &
                                            S12(i,j,k)*S12(i,j,k)*All(i,j,k)+ &
                                            S13(i,j,k)*S13(i,j,k)*All(i,j,k)+ &
                                            S21(i,j,k)*S21(i,j,k)*All(i,j,k)+ &
                                            S22(i,j,k)*S22(i,j,k)*All(i,j,k)+ &
                                            S23(i,j,k)*S23(i,j,k)*All(i,j,k)+ &
                                            S31(i,j,k)*S31(i,j,k)*All(i,j,k)+ &
                                            S32(i,j,k)*S32(i,j,k)*All(i,j,k)+ &
                                            S33(i,j,k)*S33(i,j,k)*All(i,j,k)) * &
                                          l_lim(m) * l_lim(m)
          Pim3(m) = Pim3(m) + rhof(i,j,k) *(W12(i,j,k)*W21(i,j,k)*All(i,j,k)+ &
                                            W13(i,j,k)*W31(i,j,k)*All(i,j,k)+ &
                                            W21(i,j,k)*W12(i,j,k)*All(i,j,k)+ &
                                            W23(i,j,k)*W32(i,j,k)*All(i,j,k)+ &
                                            W31(i,j,k)*W13(i,j,k)*All(i,j,k)+ &
                                            W32(i,j,k)*W23(i,j,k)*All(i,j,k)) * &
                                        l_lim(m) * l_lim(m)
          Pid(m) = Pid(m) + rhof(i,j,k) * All(i,j,k) * All(i,j,k) * All(i,j,k) * &
                                      l_lim(m) * l_lim(m)
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pis1(m) =	 psum(Pis1(m)) / (ia*ja*ka)
        Pis2(m) =	 - psum(Pis2(m)) / (ia*ja*ka)
        Pim2(m) =	 psum(Pim2(m)) / (ia*ja*ka) /3
        Pim3(m) =	 - psum(Pim3(m)) / (ia*ja*ka) /3
        Pid(m) =	 psum(Pid(m)) / (ia*ja*ka) /9
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pilocal.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pis1 pis2 pim2 pim3 pid')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_w3f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A13)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_A23)
      call fftw_free(c_A31)
      call fftw_free(c_A32)
      call fftw_free(c_A33)
      call mpistop
      deallocate(k1,k2,k3)
      deallocate(l_lim)
      deallocate(All,S11,S12,S13,S21,S22,S23,S31,S32,S33)
      deallocate(W12,W13,W21,W23,W31,W32)
      deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
      !
    end subroutine SGSPi3Dlocal
    !
    subroutine SGSLES3D(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: AllAll,SijSij
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      complex(8), allocatable, dimension(:,:,:) :: All_filted
      complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
      complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
      complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      !
      allocate(All_filted(1:im,1:jm,1:km),&
              S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
              S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
              S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km))
      !
      allocate(AllAll(1:num_l),SijSij(1:num_l))
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      AllAll = 0.d0
      SijSij = 0.d0
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
          !
          S11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
          S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
          S21_filted(i,j,k) = S12_filted(i,j,k)
          S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
          S31_filted(i,j,k) = S13_filted(i,j,k)
          S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
          S32_filted(i,j,k) = S23_filted(i,j,k)
          !
          AllAll(m) = AllAll(m) + dreal(All_filted(i,j,k))*dreal(All_filted(i,j,k))
          SijSij(m) = SijSij(m) + dreal(S11_filted(i,j,k))* dreal(S11_filted(i,j,k)) + &
                                  dreal(S12_filted(i,j,k))* dreal(S12_filted(i,j,k)) + &
                                  dreal(S13_filted(i,j,k))* dreal(S13_filted(i,j,k)) + &
                                  dreal(S21_filted(i,j,k))* dreal(S21_filted(i,j,k)) + &
                                  dreal(S22_filted(i,j,k))* dreal(S22_filted(i,j,k)) + &
                                  dreal(S23_filted(i,j,k))* dreal(S23_filted(i,j,k)) + &
                                  dreal(S31_filted(i,j,k))* dreal(S31_filted(i,j,k)) + &
                                  dreal(S32_filted(i,j,k))* dreal(S32_filted(i,j,k)) + &
                                  dreal(S33_filted(i,j,k))* dreal(S33_filted(i,j,k))
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        AllAll(m) =	 psum(AllAll(m)) / (ia*ja*ka)
        SijSij(m) =	 psum(SijSij(m)) / (ia*ja*ka)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_LES_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_LES.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell AllAll SijSij')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), AllAll(m),SijSij(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call mpistop
      deallocate(All_filted,S11_filted,S12_filted,S13_filted,&
                  S21_filted,S22_filted,S23_filted,&
                  S31_filted,S32_filted,S33_filted)
      deallocate(k1,k2,k3)
      deallocate(l_lim)
      deallocate(AllAll,SijSij)
      !
    end subroutine SGSLES3D
    !
    subroutine SGSPi3Dint(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      real(8), allocatable, dimension(:) :: Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      complex(8), allocatable, dimension(:,:,:) :: All_filted_l
      complex(8), allocatable, dimension(:,:,:) :: S11_filted_l,S12_filted_l,S13_filted_l
      complex(8), allocatable, dimension(:,:,:) :: S21_filted_l,S22_filted_l,S23_filted_l
      complex(8), allocatable, dimension(:,:,:) :: S31_filted_l,S32_filted_l,S33_filted_l
      complex(8), allocatable, dimension(:,:,:) :: All_filted
      complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
      complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
      complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
      complex(8), allocatable, dimension(:,:,:) :: W12_filted,W21_filted
      complex(8), allocatable, dimension(:,:,:) :: W13_filted,W31_filted
      complex(8), allocatable, dimension(:,:,:) :: W23_filted,W32_filted
      !
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: term1_11,term1_12,term1_13,&
                                                              term1_21,term1_22,term1_23,&
                                                              term1_31,term1_32,term1_33,&
                                                              term2,term5,term7,&
                                                              term3_11,term3_12,term3_13,&
                                                              term3_21,term3_22,term3_23,&
                                                              term3_31,term3_32,term3_33,&
                                                              term4_11,term4_12,term4_13,&
                                                              term4_21,term4_22,term4_23,&
                                                              term4_31,term4_32,term4_33,&
                                                              term6_11,term6_12,term6_13,&
                                                              term6_21,term6_22,term6_23,&
                                                              term6_31,term6_32,term6_33
      real(8) :: vxr_D1,vxr_D2,vxr_D3,vxr_D4,vxr_D5,vxr_D6,vxr_D7
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      type(C_PTR) :: c_term1_11,c_term1_12,c_term1_13,c_term1_21,c_term1_22,c_term1_23,c_term1_31,c_term1_32,c_term1_33
      type(C_PTR) :: c_term2,c_term5,c_term7
      type(C_PTR) :: c_term3_11,c_term3_12,c_term3_13,c_term3_21,c_term3_22,c_term3_23,c_term3_31,c_term3_32,c_term3_33
      type(C_PTR) :: c_term4_11,c_term4_12,c_term4_13,c_term4_21,c_term4_22,c_term4_23,c_term4_31,c_term4_32,c_term4_33
      type(C_PTR) :: c_term6_11,c_term6_12,c_term6_13,c_term6_21,c_term6_22,c_term6_23,c_term6_31,c_term6_32,c_term6_33
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi1(1:num_l), Pi2(1:num_l), Pi3(1:num_l), Pi4(1:num_l), Pi5(1:num_l), Pi6(1:num_l), Pi7(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      allocate(All_filted_l(1:im,1:jm,1:km),&
              S11_filted_l(1:im,1:jm,1:km),S12_filted_l(1:im,1:jm,1:km),S13_filted_l(1:im,1:jm,1:km),&
              S21_filted_l(1:im,1:jm,1:km),S22_filted_l(1:im,1:jm,1:km),S23_filted_l(1:im,1:jm,1:km),&
              S31_filted_l(1:im,1:jm,1:km),S32_filted_l(1:im,1:jm,1:km),S33_filted_l(1:im,1:jm,1:km))
      !
      allocate(All_filted(1:im,1:jm,1:km),&
              S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
              S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
              S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km),&
              W12_filted(1:im,1:jm,1:km),W21_filted(1:im,1:jm,1:km),&
              W13_filted(1:im,1:jm,1:km),W31_filted(1:im,1:jm,1:km),&
              W23_filted(1:im,1:jm,1:km),W32_filted(1:im,1:jm,1:km))
      !
      c_term1_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_11, term1_11, [imfftw,jmfftw,kmfftw])
      c_term1_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_12, term1_12, [imfftw,jmfftw,kmfftw])
      c_term1_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_13, term1_13, [imfftw,jmfftw,kmfftw])
      c_term1_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_21, term1_21, [imfftw,jmfftw,kmfftw])
      c_term1_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_22, term1_22, [imfftw,jmfftw,kmfftw])
      c_term1_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_23, term1_23, [imfftw,jmfftw,kmfftw])
      c_term1_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_31, term1_31, [imfftw,jmfftw,kmfftw])
      c_term1_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_32, term1_32, [imfftw,jmfftw,kmfftw])
      c_term1_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_33, term1_33, [imfftw,jmfftw,kmfftw])
      c_term2    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term2,    term2,    [imfftw,jmfftw,kmfftw])
      c_term3_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_11, term3_11, [imfftw,jmfftw,kmfftw])
      c_term3_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_12, term3_12, [imfftw,jmfftw,kmfftw])
      c_term3_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_13, term3_13, [imfftw,jmfftw,kmfftw])
      c_term3_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_21, term3_21, [imfftw,jmfftw,kmfftw])
      c_term3_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_22, term3_22, [imfftw,jmfftw,kmfftw])
      c_term3_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_23, term3_23, [imfftw,jmfftw,kmfftw])
      c_term3_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_31, term3_31, [imfftw,jmfftw,kmfftw])
      c_term3_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_32, term3_32, [imfftw,jmfftw,kmfftw])
      c_term3_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_33, term3_33, [imfftw,jmfftw,kmfftw])
      c_term4_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_11, term4_11, [imfftw,jmfftw,kmfftw])
      c_term4_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_12, term4_12, [imfftw,jmfftw,kmfftw])
      c_term4_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_13, term4_13, [imfftw,jmfftw,kmfftw])
      c_term4_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_21, term4_21, [imfftw,jmfftw,kmfftw])
      c_term4_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_22, term4_22, [imfftw,jmfftw,kmfftw])
      c_term4_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_23, term4_23, [imfftw,jmfftw,kmfftw])
      c_term4_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_31, term4_31, [imfftw,jmfftw,kmfftw])
      c_term4_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_32, term4_32, [imfftw,jmfftw,kmfftw])
      c_term4_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_33, term4_33, [imfftw,jmfftw,kmfftw])
      c_term5    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term5,    term5,    [imfftw,jmfftw,kmfftw])
      c_term6_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_11, term6_11, [imfftw,jmfftw,kmfftw])
      c_term6_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_12, term6_12, [imfftw,jmfftw,kmfftw])
      c_term6_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_13, term6_13, [imfftw,jmfftw,kmfftw])
      c_term6_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_21, term6_21, [imfftw,jmfftw,kmfftw])
      c_term6_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_22, term6_22, [imfftw,jmfftw,kmfftw])
      c_term6_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_23, term6_23, [imfftw,jmfftw,kmfftw])
      c_term6_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_31, term6_31, [imfftw,jmfftw,kmfftw])
      c_term6_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_32, term6_32, [imfftw,jmfftw,kmfftw])
      c_term6_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_33, term6_33, [imfftw,jmfftw,kmfftw])
      c_term7    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term7,    term7,    [imfftw,jmfftw,kmfftw])
      !
      !
      Pi1 = 0.d0
      Pi2 =	0.d0
      Pi3 =	0.d0
      Pi4 =	0.d0
      Pi5 =	0.d0
      Pi6 =	0.d0
      Pi7 =	0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename2,handle=hand_b, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        endif
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          !
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All_filted_l(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
          !
          S11_filted_l(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
          S22_filted_l(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
          S33_filted_l(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted_l(i,j,k)
          S12_filted_l(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
          S21_filted_l(i,j,k) = S12_filted_l(i,j,k)
          S13_filted_l(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
          S31_filted_l(i,j,k) = S13_filted_l(i,j,k)
          S23_filted_l(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
          S32_filted_l(i,j,k) = S23_filted_l(i,j,k)
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !!!! Velocity Favre average and density average
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          do k=1,km
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
            w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
            w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
            rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
          enddo
          enddo
          enddo
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          do k=1,km
          do j=1,jm
          do i=1,im
            w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
            w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
            w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
            w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
            w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
            !
            A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
            A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
            A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
            A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
            A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
            A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
            A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
            A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
            A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
            !
          end do
          end do
          end do
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
          call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
          call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
          call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
          call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
          call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
          call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
          call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
          call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
          !
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            All_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
            !
            S11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
            S22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
            S33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted(i,j,k)
            S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
            S21_filted(i,j,k) = S12_filted(i,j,k)
            S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
            S31_filted(i,j,k) = S13_filted(i,j,k)
            S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
            S32_filted(i,j,k) = S23_filted(i,j,k)
            !
            W12_filted(i,j,k) = (A12_filted(i,j,k)-A21_filted(i,j,k))*0.5d0
            W21_filted(i,j,k) = -1.d0*W12_filted(i,j,k)
            W13_filted(i,j,k) = (A13_filted(i,j,k)-A31_filted(i,j,k))*0.5d0
            W31_filted(i,j,k) = -1.d0*W13_filted(i,j,k)
            W23_filted(i,j,k) = (A23_filted(i,j,k)-A32_filted(i,j,k))*0.5d0
            W32_filted(i,j,k) = -1.d0*W23_filted(i,j,k)
            !
          end do
          end do
          end do
          !
          !!!! Pi terms
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !term1_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted + rho_filted*SI3_filted*SJ3_filted
            term1_11(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S13_filted(i,j,k)
            term1_12(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S23_filted(i,j,k)
            term1_13(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S33_filted(i,j,k)
            term1_21(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S13_filted(i,j,k)
            term1_22(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S23_filted(i,j,k)
            term1_23(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S33_filted(i,j,k)
            term1_31(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S13_filted(i,j,k)
            term1_32(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S23_filted(i,j,k)
            term1_33(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S33_filted(i,j,k)
            ! 
            ! term2 
            term2(i,j,k) =  rho_filted(i,j,k)*S11_filted(i,j,k)*S11_filted(i,j,k)+&
                            rho_filted(i,j,k)*S12_filted(i,j,k)*S12_filted(i,j,k)+&
                            rho_filted(i,j,k)*S13_filted(i,j,k)*S13_filted(i,j,k)+& ! i=1,k=1,2,3
                            rho_filted(i,j,k)*S21_filted(i,j,k)*S21_filted(i,j,k)+&
                            rho_filted(i,j,k)*S22_filted(i,j,k)*S22_filted(i,j,k)+&
                            rho_filted(i,j,k)*S23_filted(i,j,k)*S23_filted(i,j,k)+& ! i=2,k=1,2,3
                            rho_filted(i,j,k)*S31_filted(i,j,k)*S31_filted(i,j,k)+&
                            rho_filted(i,j,k)*S32_filted(i,j,k)*S32_filted(i,j,k)+&
                            rho_filted(i,j,k)*S33_filted(i,j,k)*S33_filted(i,j,k) ! i=3,k=1,2,3
            !
            ! term3_IJ = rho_filted*All_filted*SIJ_filted
            term3_11(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S11_filted(i,j,k)
            term3_12(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S12_filted(i,j,k)
            term3_13(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S13_filted(i,j,k)
            !
            term3_21(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S21_filted(i,j,k)
            term3_22(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S22_filted(i,j,k)
            term3_23(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S23_filted(i,j,k)
            !
            term3_31(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S31_filted(i,j,k)
            term3_32(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S32_filted(i,j,k)
            term3_33(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S33_filted(i,j,k)
            !
            !term4_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted + &
            !rho_filted*WI3_filted*W3J_filted 
            term4_11(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W21_filted(i,j,k) + &
                              rho_filted(i,j,k)*W13_filted(i,j,k)*W31_filted(i,j,k) 
            term4_21(i,j,k) = rho_filted(i,j,k)*W23_filted(i,j,k)*W31_filted(i,j,k) 
            term4_31(i,j,k) = rho_filted(i,j,k)*W32_filted(i,j,k)*W21_filted(i,j,k)
            !
            term4_12(i,j,k) = rho_filted(i,j,k)*W13_filted(i,j,k)*W32_filted(i,j,k) 
            term4_22(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W12_filted(i,j,k) + &
                              rho_filted(i,j,k)*W23_filted(i,j,k)*W32_filted(i,j,k) 
            term4_32(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W12_filted(i,j,k)
            !
            term4_13(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W23_filted(i,j,k)
            term4_23(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W13_filted(i,j,k)
            term4_33(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W13_filted(i,j,k) + &
                              rho_filted(i,j,k)*W32_filted(i,j,k)*W23_filted(i,j,k)
            !
            ! term5
            term5(i,j,k) = 2.d0*rho_filted(i,j,k)*(W12_filted(i,j,k)*W21_filted(i,j,k) + &
                            W13_filted(i,j,k)*W31_filted(i,j,k) + W23_filted(i,j,k)*W32_filted(i,j,k))
            !
            !term6_IJ= rho_filted*(S1J_filted*WI1_filted-SI1_filted*W1J_filted) + &
            !          rho_filted*(S2J_filted*WI2_filted-SI2_filted*W2J_filted) + &
            !          rho_filted*(S3J_filted*WI3_filted-SI3_filted*W3J_filted)
            term6_11(i,j,k)= rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                            -rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                            +rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k)
            term6_21(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S22_filted(i,j,k)*W21_filted(i,j,k) &
                            +rho_filted(i,j,k)*S31_filted(i,j,k)*W23_filted(i,j,k) &
                            -rho_filted(i,j,k)*S23_filted(i,j,k)*W31_filted(i,j,k)
            term6_31(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W31_filted(i,j,k) &
                            +rho_filted(i,j,k)*S21_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S32_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S33_filted(i,j,k)*W31_filted(i,j,k) 
            !
            term6_12(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S22_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S32_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S13_filted(i,j,k)*W32_filted(i,j,k)
            term6_22(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k) &
                            -rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k)
            term6_32(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W31_filted(i,j,k) &
                            -rho_filted(i,j,k)*S31_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S22_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S33_filted(i,j,k)*W32_filted(i,j,k)
            !
            term6_13(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W13_filted(i,j,k) &
                            +rho_filted(i,j,k)*S23_filted(i,j,k)*W12_filted(i,j,k) &
                            -rho_filted(i,j,k)*S12_filted(i,j,k)*W23_filted(i,j,k) &
                            +rho_filted(i,j,k)*S33_filted(i,j,k)*W13_filted(i,j,k)
            term6_23(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S21_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S22_filted(i,j,k)*W23_filted(i,j,k) &
                            +rho_filted(i,j,k)*S33_filted(i,j,k)*W23_filted(i,j,k)
            term6_33(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k) &
                            -rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                            +rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k)
            !
            ! term7
            term7(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*All_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(forward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(forward_plan,term1_13,term1_13)
          call fftw_mpi_execute_dft(forward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(forward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(forward_plan,term1_23,term1_23)
          call fftw_mpi_execute_dft(forward_plan,term1_31,term1_31)
          call fftw_mpi_execute_dft(forward_plan,term1_32,term1_32)
          call fftw_mpi_execute_dft(forward_plan,term1_33,term1_33)
          !
          call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
          !
          call fftw_mpi_execute_dft(forward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(forward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(forward_plan,term3_13,term3_13)
          call fftw_mpi_execute_dft(forward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(forward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(forward_plan,term3_23,term3_23)
          call fftw_mpi_execute_dft(forward_plan,term3_31,term3_31)
          call fftw_mpi_execute_dft(forward_plan,term3_32,term3_32)
          call fftw_mpi_execute_dft(forward_plan,term3_33,term3_33)
          !
          call fftw_mpi_execute_dft(forward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(forward_plan,term4_12,term4_12)
          call fftw_mpi_execute_dft(forward_plan,term4_13,term4_13)
          call fftw_mpi_execute_dft(forward_plan,term4_21,term4_21)
          call fftw_mpi_execute_dft(forward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(forward_plan,term4_23,term4_23)
          call fftw_mpi_execute_dft(forward_plan,term4_31,term4_31)
          call fftw_mpi_execute_dft(forward_plan,term4_32,term4_32)
          call fftw_mpi_execute_dft(forward_plan,term4_33,term4_33)
          !
          call fftw_mpi_execute_dft(forward_plan,term5   ,term5   )
          !
          call fftw_mpi_execute_dft(forward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(forward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(forward_plan,term6_13,term6_13)
          call fftw_mpi_execute_dft(forward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(forward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(forward_plan,term6_23,term6_23)
          call fftw_mpi_execute_dft(forward_plan,term6_31,term6_31)
          call fftw_mpi_execute_dft(forward_plan,term6_32,term6_32)
          call fftw_mpi_execute_dft(forward_plan,term6_33,term6_33)
          !
          call fftw_mpi_execute_dft(forward_plan,term7   ,term7   )
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            term1_11(i,j,k) = term1_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_12(i,j,k) = term1_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_13(i,j,k) = term1_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_21(i,j,k) = term1_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_22(i,j,k) = term1_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_23(i,j,k) = term1_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_31(i,j,k) = term1_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_32(i,j,k) = term1_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term1_33(i,j,k) = term1_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            term2(i,j,k)    = term2(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
            !
            term3_11(i,j,k) = term3_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_12(i,j,k) = term3_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_13(i,j,k) = term3_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_21(i,j,k) = term3_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_22(i,j,k) = term3_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_23(i,j,k) = term3_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_31(i,j,k) = term3_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_32(i,j,k) = term3_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term3_33(i,j,k) = term3_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            term4_11(i,j,k) = term4_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_12(i,j,k) = term4_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_13(i,j,k) = term4_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_21(i,j,k) = term4_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_22(i,j,k) = term4_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_23(i,j,k) = term4_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_31(i,j,k) = term4_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_32(i,j,k) = term4_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term4_33(i,j,k) = term4_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            term5(i,j,k)    = term5(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
            !
            term6_11(i,j,k) = term6_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_12(i,j,k) = term6_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_13(i,j,k) = term6_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_21(i,j,k) = term6_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_22(i,j,k) = term6_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_23(i,j,k) = term6_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_31(i,j,k) = term6_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_32(i,j,k) = term6_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            term6_33(i,j,k) = term6_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            term7(i,j,k)    = term7(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
            !
          enddo
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(backward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(backward_plan,term1_13,term1_13)
          call fftw_mpi_execute_dft(backward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(backward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(backward_plan,term1_23,term1_23)
          call fftw_mpi_execute_dft(backward_plan,term1_31,term1_31)
          call fftw_mpi_execute_dft(backward_plan,term1_32,term1_32)
          call fftw_mpi_execute_dft(backward_plan,term1_33,term1_33)
          !
          call fftw_mpi_execute_dft(backward_plan,term2   ,term2   )
          !
          call fftw_mpi_execute_dft(backward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(backward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(backward_plan,term3_13,term3_13)
          call fftw_mpi_execute_dft(backward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(backward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(backward_plan,term3_23,term3_23)
          call fftw_mpi_execute_dft(backward_plan,term3_31,term3_31)
          call fftw_mpi_execute_dft(backward_plan,term3_32,term3_32)
          call fftw_mpi_execute_dft(backward_plan,term3_33,term3_33)
          !
          call fftw_mpi_execute_dft(backward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(backward_plan,term4_12,term4_12)
          call fftw_mpi_execute_dft(backward_plan,term4_13,term4_13)
          call fftw_mpi_execute_dft(backward_plan,term4_21,term4_21)
          call fftw_mpi_execute_dft(backward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(backward_plan,term4_23,term4_23)
          call fftw_mpi_execute_dft(backward_plan,term4_31,term4_31)
          call fftw_mpi_execute_dft(backward_plan,term4_32,term4_32)
          call fftw_mpi_execute_dft(backward_plan,term4_33,term4_33)
          !
          call fftw_mpi_execute_dft(backward_plan,term5   ,term5   )
          !
          call fftw_mpi_execute_dft(backward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(backward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(backward_plan,term6_13,term6_13)
          call fftw_mpi_execute_dft(backward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(backward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(backward_plan,term6_23,term6_23)
          call fftw_mpi_execute_dft(backward_plan,term6_31,term6_31)
          call fftw_mpi_execute_dft(backward_plan,term6_32,term6_32)
          call fftw_mpi_execute_dft(backward_plan,term6_33,term6_33)
          !
          call fftw_mpi_execute_dft(backward_plan,term7   ,term7   )
          !
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            vxr_D1 = real(term1_11(i,j,k) * S11_filted_l(i,j,k) + &
                    term1_12(i,j,k) * S12_filted_l(i,j,k) + &
                    term1_13(i,j,k) * S13_filted_l(i,j,k) + &
                    term1_21(i,j,k) * S21_filted_l(i,j,k) + &
                    term1_22(i,j,k) * S22_filted_l(i,j,k) + &
                    term1_23(i,j,k) * S23_filted_l(i,j,k) + &
                    term1_31(i,j,k) * S31_filted_l(i,j,k) + &
                    term1_32(i,j,k) * S32_filted_l(i,j,k) + &
                    term1_33(i,j,k) * S33_filted_l(i,j,k))
            Pi1(m) = Pi1(m) + vxr_D1 * dl_alpha(m,n)
            !
            vxr_D2 = real(term2(i,j,k) * All_filted_l(i,j,k))
            Pi2(m) = Pi2(m) + vxr_D2 * dl_alpha(m,n) / 3.d0
            !
            vxr_D3 = real(term3_11(i,j,k) * S11_filted_l(i,j,k) + &
                    term3_12(i,j,k) * S12_filted_l(i,j,k) + &
                    term3_13(i,j,k) * S13_filted_l(i,j,k) + &
                    term3_21(i,j,k) * S21_filted_l(i,j,k) + &
                    term3_22(i,j,k) * S22_filted_l(i,j,k) + &
                    term3_23(i,j,k) * S23_filted_l(i,j,k) + &
                    term3_31(i,j,k) * S31_filted_l(i,j,k) + &
                    term3_32(i,j,k) * S32_filted_l(i,j,k) + &
                    term3_33(i,j,k) * S33_filted_l(i,j,k))
            Pi3(m) = Pi3(m) + vxr_D3 * dl_alpha(m,n) * 2.d0/3.d0
            !
            vxr_D4 = real(term4_11(i,j,k) * S11_filted_l(i,j,k) + &
                    term4_12(i,j,k) * S12_filted_l(i,j,k) + &
                    term4_13(i,j,k) * S13_filted_l(i,j,k) + &
                    term4_21(i,j,k) * S21_filted_l(i,j,k) + &
                    term4_22(i,j,k) * S22_filted_l(i,j,k) + &
                    term4_23(i,j,k) * S23_filted_l(i,j,k) + &
                    term4_31(i,j,k) * S31_filted_l(i,j,k) + &
                    term4_32(i,j,k) * S32_filted_l(i,j,k) + &
                    term4_33(i,j,k) * S33_filted_l(i,j,k))
            Pi4(m) = Pi4(m) - vxr_D4 * dl_alpha(m,n)
            !
            vxr_D5 = real(term5(i,j,k) * All_filted_l(i,j,k))
            Pi5(m) = Pi5(m) - vxr_D5 * dl_alpha(m,n) / 3.d0
            !
            vxr_D6 = real(term6_11(i,j,k) * S11_filted_l(i,j,k) + &
                    term6_12(i,j,k) * S12_filted_l(i,j,k) + &
                    term6_13(i,j,k) * S13_filted_l(i,j,k) + &
                    term6_21(i,j,k) * S21_filted_l(i,j,k) + &
                    term6_22(i,j,k) * S22_filted_l(i,j,k) + &
                    term6_23(i,j,k) * S23_filted_l(i,j,k) + &
                    term6_31(i,j,k) * S31_filted_l(i,j,k) + &
                    term6_32(i,j,k) * S32_filted_l(i,j,k) + &
                    term6_33(i,j,k) * S33_filted_l(i,j,k))
            Pi6(m) = Pi6(m) + vxr_D6 * dl_alpha(m,n)
            !
            vxr_D7 = term7(i,j,k) * All_filted_l(i,j,k)
            Pi7(m) = Pi7(m) + vxr_D7 * dl_alpha(m,n) / 9.d0
          enddo
          enddo
          enddo
          !
          vxr_D1 = psum(vxr_D1)
          vxr_D2 = psum(vxr_D2)
          vxr_D3 = psum(vxr_D3)
          vxr_D4 = psum(vxr_D4)
          vxr_D5 = psum(vxr_D5)
          vxr_D6 = psum(vxr_D6)
          vxr_D7 = psum(vxr_D7)
          !
          if(mpirank==0) then
            call listwrite(hand_b,l_sqrtalpha(m,n),vxr_D1 * dl_alpha(m,n), &
                          vxr_D2 * dl_alpha(m,n) / 3.d0, &
                          vxr_D3 * dl_alpha(m,n) * 2.d0/3.d0, &
                          - vxr_D4 * dl_alpha(m,n), &
                          - vxr_D5 * dl_alpha(m,n) / 3.d0, &
                          vxr_D6 * dl_alpha(m,n), &
                          vxr_D7 * dl_alpha(m,n) / 9.d0)
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        Pi1(m) =	 psum(Pi1(m)) / (ia*ja*ka)
        Pi2(m) =	 psum(Pi2(m)) / (ia*ja*ka)
        Pi3(m) =	 psum(Pi3(m)) / (ia*ja*ka)
        Pi4(m) =	 psum(Pi4(m)) / (ia*ja*ka)
        Pi5(m) =	 psum(Pi5(m)) / (ia*ja*ka)
        Pi6(m) =	 psum(Pi6(m)) / (ia*ja*ka)
        Pi7(m) =	 psum(Pi7(m)) / (ia*ja*ka)
        !
        if(mpirank==0)then
          print *, '>>>>', outfilename2
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pi.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi1(m), Pi2(m), &
                          Pi3(m), Pi4(m), Pi5(m),     &
                          Pi6(m), Pi7(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call fftw_free(c_term1_11)
      call fftw_free(c_term1_12)
      call fftw_free(c_term1_13)
      call fftw_free(c_term1_21)
      call fftw_free(c_term1_22)
      call fftw_free(c_term1_23)
      call fftw_free(c_term1_31)
      call fftw_free(c_term1_32)
      call fftw_free(c_term1_33)
      call fftw_free(c_term2)
      call fftw_free(c_term5)
      call fftw_free(c_term7)
      call fftw_free(c_term3_11)
      call fftw_free(c_term3_12)
      call fftw_free(c_term3_13)
      call fftw_free(c_term3_21)
      call fftw_free(c_term3_22)
      call fftw_free(c_term3_23)
      call fftw_free(c_term3_31)
      call fftw_free(c_term3_32)
      call fftw_free(c_term3_33)
      call fftw_free(c_term4_11)
      call fftw_free(c_term4_12)
      call fftw_free(c_term4_13)
      call fftw_free(c_term4_21)
      call fftw_free(c_term4_22)
      call fftw_free(c_term4_23)
      call fftw_free(c_term4_31)
      call fftw_free(c_term4_32)
      call fftw_free(c_term4_33)
      call fftw_free(c_term6_11)
      call fftw_free(c_term6_12)
      call fftw_free(c_term6_13)
      call fftw_free(c_term6_21)
      call fftw_free(c_term6_22)
      call fftw_free(c_term6_23)
      call fftw_free(c_term6_31)
      call fftw_free(c_term6_32)
      call fftw_free(c_term6_33)
      call mpistop
      deallocate(All_filted_l,S11_filted_l,S12_filted_l,S13_filted_l)
      deallocate(S21_filted_l,S22_filted_l,S23_filted_l)
      deallocate(S31_filted_l,S32_filted_l,S33_filted_l)
      deallocate(All_filted,S11_filted,S12_filted,S13_filted)
      deallocate(S21_filted,S22_filted,S23_filted)
      deallocate(S31_filted,S32_filted,S33_filted)
      deallocate(W12_filted,W21_filted,W13_filted,W31_filted,W23_filted,W32_filted)
      deallocate(k1,k2,k3)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7)
      !
    end subroutine SGSPi3Dint
    !
    subroutine SGSstress2D(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n,p,q
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname,mname
      character(len=10):: termname
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhol
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: tau11_term,tau12_term,tau21_term,tau22_term
      real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:,:,:) :: tau_bis ! 1:2,1:2,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
      real(8) :: result,norm2,norm2bis
      real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
      !
      !
      type(C_PTR) :: forward_plan,backward_plan
      type(C_PTR) :: c_w1,c_w2,c_rhol
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
      type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
      type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau21_term,c_tau22_term
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
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
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
          print *,"Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      !
      allocate(tau(1:2,1:2,1:im,1:jm),tau_bis(1:2,1:2,1:im,1:jm))
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      do i=1,num_l
        dl_alpha(i,1) = l_sqrtalpha(i,1)**2 
        !
        do j=2,num_alphas(i)
          dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
        enddo
      enddo
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhol = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw])
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11, [imfftw,jmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12, [imfftw,jmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21, [imfftw,jmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22, [imfftw,jmfftw])
      !
      c_tau11_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw])
      c_tau12_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw])
      c_tau21_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw])
      c_tau22_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! wi,wiwj,rhol in physical space, not filted
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhol(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
        w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
        w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
        w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
        !
        rhol(i,j)=rhol(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      do m=1,num_l
        !
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        ! Method 1: tauij = rho uiuj - rho ui uj
        ! Filter scale: l
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! 
          !
          w1_filted(i,j)    = w1(i,j)   *Gl
          w2_filted(i,j)    = w2(i,j)   *Gl
          !
          w1w1_filted(i,j)  = w1w1(i,j) *Gl
          w1w2_filted(i,j)  = w1w2(i,j) *Gl
          w2w1_filted(i,j)  = w2w1(i,j) *Gl
          w2w2_filted(i,j)  = w2w2(i,j) *Gl
          !
          rho_filted(i,j)   = rhol(i,j) *Gl
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          !
          tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
          tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
          tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
          tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
          !
          do p=1,2
          do q=1,2
            tau_bis(p,q,i,j) = 0.d0
          enddo
          enddo
          !
        enddo
        enddo
        !
        ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
        ! Filter scale: l
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j)  = w1(i,j)  *Galpha
            w2_filted(i,j)  = w2(i,j)  *Galpha
            rho_filted(i,j) = rhol(i,j)*Galpha
          enddo
          enddo
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
          do j=1,jm
          do i=1,im
            w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
            w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          enddo
          enddo
          !
          ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
          ! rho_filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
            w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
            !
            A11(i,j) = imag*w1_filted(i,j)*k1(i,j)
            A21(i,j) = imag*w2_filted(i,j)*k1(i,j)
            A12(i,j) = imag*w1_filted(i,j)*k2(i,j)
            A22(i,j) = imag*w2_filted(i,j)*k2(i,j)
            !
          end do
          end do
          !
          ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
          ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(backward_plan,A11,A11)
          call fftw_mpi_execute_dft(backward_plan,A21,A21)
          call fftw_mpi_execute_dft(backward_plan,A12,A12)
          call fftw_mpi_execute_dft(backward_plan,A22,A22)
          !
          do j=1,jm
          do i=1,im
            !
            tau11_term(i,j) = rho_filted(i,j) * (A11(i,j)*A11(i,j)+A12(i,j)*A12(i,j))
            tau12_term(i,j) = rho_filted(i,j) * (A11(i,j)*A21(i,j)+A12(i,j)*A22(i,j))
            tau21_term(i,j) = rho_filted(i,j) * (A21(i,j)*A11(i,j)+A22(i,j)*A12(i,j))
            tau22_term(i,j) = rho_filted(i,j) * (A21(i,j)*A21(i,j)+A22(i,j)*A22(i,j))
            !
          end do
          end do
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
          !
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            tau11_term(i,j) = tau11_term(i,j)*Gphi/(1.d0*ia*ja)
            tau12_term(i,j) = tau12_term(i,j)*Gphi/(1.d0*ia*ja)
            tau21_term(i,j) = tau21_term(i,j)*Gphi/(1.d0*ia*ja)
            tau22_term(i,j) = tau22_term(i,j)*Gphi/(1.d0*ia*ja)
            !
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
          !
          !
          do j=1,jm
          do i=1,im
            !
            tau_bis(1,1,i,j) = tau_bis(1,1,i,j) + dreal(tau11_term(i,j)) * dl_alpha(m,n)
            tau_bis(1,2,i,j) = tau_bis(1,2,i,j) + dreal(tau12_term(i,j)) * dl_alpha(m,n)
            tau_bis(2,1,i,j) = tau_bis(2,1,i,j) + dreal(tau21_term(i,j)) * dl_alpha(m,n)
            tau_bis(2,2,i,j) = tau_bis(2,2,i,j) + dreal(tau22_term(i,j)) * dl_alpha(m,n)
            !
          enddo
          enddo
          !
        enddo ! loop of integral (alpha)
        !
        do p=1,2
        do q=1,2
          erroravg(p,q)=0.d0
          errormax(p,q)=0.d0
          errorgtr10(p,q)=0.d0
          errorgtr100(p,q)=0.d0
          errornorm2max=0.d0
          errornorm2avg=0.d0
          errornorm2gtr10=0.d0
          errornorm2gtr100=0.d0
        enddo
        enddo
        !
        ! Output comparaison results
        do j=1,jm
        do i=1,im
          !
          norm2 = 0.d0
          norm2bis = 0.d0
          do p=1,2
          do q=1,2
            if(abs(tau(p,q,i,j))>1.d-6)then
              result = abs(tau_bis(p,q,i,j)-tau(p,q,i,j))/(abs(tau_bis(p,q,i,j))+abs(tau(p,q,i,j)))
              norm2 = norm2 + tau(p,q,i,j)**2
              norm2bis = norm2bis + tau_bis(p,q,i,j)**2
              errormax(p,q) = max(errormax(p,q),result)
              erroravg(p,q) = erroravg(p,q) + result
              if(result > 0.1)then
                errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
              endif
              if(result > 1)then
                errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
              endif
            endif
          enddo
          enddo
          result = abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
          errornorm2max = max(errornorm2max,result)
          errornorm2avg = errornorm2avg + result
          if(result > 0.1)then
            errornorm2gtr10 = errornorm2gtr10 + 1.d0
          endif
          if(result > 1)then
            errornorm2gtr100 = errornorm2gtr100 + 1.d0
          endif
          !
        enddo
        enddo
        !
        !
        do p=1,2
        do q=1,2
          errormax(p,q) = pmax(errormax(p,q))
          erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja)
          errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja)
          errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja)
          errornorm2max = psum(errornorm2max)
          errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja)
          errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja)
          errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja)
        enddo
        enddo
        !
        if(mpirank==0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
          endif
          !
          if(mpirank == 0)then
            open(fh,file=outfilename,form='formatted')
            write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
            write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
            write(fh,"(A8,1x,5(A20,1x))")'type','tau11','tau12','tau21','tau22','norm2'
            write(fh,"(A8,1x,5(E20.13E2,1x))")'max',((errormax(p,q),q=1,2),p=1,2),errornorm2max
            write(fh,"(A8,1x,5(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,2),p=1,2),errornorm2avg
            write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,2),p=1,2),errornorm2gtr10
            write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,2),p=1,2),errornorm2gtr100
            close(fh)
            print *, '>>>>', outfilename
          endif
          !
          !
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(loutput)then
          !
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
          else
            outfilename = 'pp/SGS_stress_'//mname//'.h5'
          endif
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          do p=1,2
          do q= 1,2
            write (termname, "(A3,I1,I1)") "tau",p,q
            call h5wa2d_r8(varname=termname,var=tau(p,q,1:im,1:jm),    dir='k')
            write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
            call h5wa2d_r8(varname=termname,var=tau_bis(p,q,1:im,1:jm),dir='k')
          enddo
          enddo
          !
          call h5io_end
          !
        endif
        !
      enddo ! loop of filter point l
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhol)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_tau11_term)
      call fftw_free(c_tau12_term)
      call fftw_free(c_tau21_term)
      call fftw_free(c_tau22_term)
      call mpistop
      deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
      !
    end subroutine SGSstress2D
    !
    subroutine SGSstress3D(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n,p,q
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname,mname
      character(len=10):: termname
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhol
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau11_term,tau12_term,tau13_term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau21_term,tau22_term,tau23_term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau31_term,tau32_term,tau33_term
      real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:,:,:,:) :: tau_bis ! 1:3,1:3,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
      real(8) :: result,norm2,norm2bis
      real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
      !
      !
      type(C_PTR) :: forward_plan,backward_plan
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhol
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
      type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
      type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
      type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
      type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
      type(C_PTR) :: c_A11,c_A12,c_A13
      type(C_PTR) :: c_A21,c_A22,c_A23
      type(C_PTR) :: c_A31,c_A32,c_A33
      type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau13_term,c_tau21_term,c_tau22_term
      type(C_PTR) :: c_tau23_term,c_tau31_term,c_tau32_term,c_tau33_term
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      allocate(tau(1:3,1:3,1:im,1:jm,1:km),tau_bis(1:3,1:3,1:im,1:jm,1:km))
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhol = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw,kmfftw])
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted, [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
      c_w1w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
      c_w2w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
      c_w3w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
      c_w3w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
      c_w3w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
      c_w1w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
      c_w2w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
      c_w3w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
      c_w3w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11, [imfftw,jmfftw,kmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12, [imfftw,jmfftw,kmfftw])
      c_A13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13, A13, [imfftw,jmfftw,kmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21, [imfftw,jmfftw,kmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22, [imfftw,jmfftw,kmfftw])
      c_A23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23, A23, [imfftw,jmfftw,kmfftw])
      c_A31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31, A31, [imfftw,jmfftw,kmfftw])
      c_A32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32, A32, [imfftw,jmfftw,kmfftw])
      c_A33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33, A33, [imfftw,jmfftw,kmfftw])
      !
      c_tau11_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw,kmfftw])
      c_tau12_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw,kmfftw])
      c_tau13_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau13_term, tau13_term, [imfftw,jmfftw,kmfftw])
      c_tau21_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw,kmfftw])
      c_tau22_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw,kmfftw])
      c_tau23_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau23_term, tau23_term, [imfftw,jmfftw,kmfftw])
      c_tau31_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau31_term, tau31_term, [imfftw,jmfftw,kmfftw])
      c_tau32_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau32_term, tau32_term, [imfftw,jmfftw,kmfftw])
      c_tau33_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau33_term, tau33_term, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! wi,wiwj,rhol in physical space, not filted
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhol(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
      call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
      call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
      call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
      call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
        w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
        w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
        w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
        w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
        w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
        w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
        w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
        w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhol(i,j,k)=rhol(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
      !
      do m=1,num_l
        !
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        ! Method 1: tauij = rho uiuj - rho ui uj
        ! Filter scale: l
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0)
          !
          w1_filted(i,j,k)    = w1(i,j,k)   *Gl
          w2_filted(i,j,k)    = w2(i,j,k)   *Gl
          w3_filted(i,j,k)    = w3(i,j,k)   *Gl
          !
          w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
          w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
          w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
          w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
          w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
          w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
          w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
          w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
          w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
          !
          rho_filted(i,j,k)   = rhol(i,j,k) *Gl
        enddo
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w1_filted(i,j,k))
          tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w2_filted(i,j,k))
          tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w3_filted(i,j,k))
          tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w1_filted(i,j,k))
          tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w2_filted(i,j,k))
          tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w3_filted(i,j,k))
          tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w1_filted(i,j,k))
          tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w2_filted(i,j,k))
          tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w3_filted(i,j,k))
          !
          do p=1,3
          do q=1,3
            tau_bis(p,q,i,j,k) = 0.d0
          enddo
          enddo
          !
        enddo
        enddo
        enddo
        !
        ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
        ! Filter scale: l
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
          do k=1,km
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0)
            w1_filted(i,j,k)  = w1(i,j,k)  *Galpha
            w2_filted(i,j,k)  = w2(i,j,k)  *Galpha
            w3_filted(i,j,k)  = w3(i,j,k)  *Galpha
            rho_filted(i,j,k) = rhol(i,j,k)*Galpha
          enddo
          enddo
          enddo
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
          do k=1,km
          do j=1,jm
          do i=1,im
            w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
            w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
            w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
          ! rho_filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
            w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
            w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
            !
            A11(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
            A21(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
            A31(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
            A12(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
            A22(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
            A32(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
            A13(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
            A23(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
            A33(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
            !
          end do
          end do
          end do
          !
          ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
          ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(backward_plan,A11,A11)
          call fftw_mpi_execute_dft(backward_plan,A21,A21)
          call fftw_mpi_execute_dft(backward_plan,A31,A31)
          call fftw_mpi_execute_dft(backward_plan,A12,A12)
          call fftw_mpi_execute_dft(backward_plan,A22,A22)
          call fftw_mpi_execute_dft(backward_plan,A32,A32)
          call fftw_mpi_execute_dft(backward_plan,A13,A13)
          call fftw_mpi_execute_dft(backward_plan,A23,A23)
          call fftw_mpi_execute_dft(backward_plan,A33,A33)
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            tau11_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A11(i,j,k)+A12(i,j,k)*A12(i,j,k)+A13(i,j,k)*A13(i,j,k))
            tau12_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A21(i,j,k)+A12(i,j,k)*A22(i,j,k)+A13(i,j,k)*A23(i,j,k))
            tau13_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A31(i,j,k)+A12(i,j,k)*A32(i,j,k)+A13(i,j,k)*A33(i,j,k))
            tau21_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A11(i,j,k)+A22(i,j,k)*A12(i,j,k)+A23(i,j,k)*A13(i,j,k))
            tau22_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A21(i,j,k)+A22(i,j,k)*A22(i,j,k)+A23(i,j,k)*A23(i,j,k))
            tau23_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A31(i,j,k)+A22(i,j,k)*A32(i,j,k)+A23(i,j,k)*A33(i,j,k))
            tau31_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A11(i,j,k)+A32(i,j,k)*A12(i,j,k)+A33(i,j,k)*A13(i,j,k))
            tau32_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A21(i,j,k)+A32(i,j,k)*A22(i,j,k)+A33(i,j,k)*A23(i,j,k))
            tau33_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A31(i,j,k)+A32(i,j,k)*A32(i,j,k)+A33(i,j,k)*A33(i,j,k))
            !
          end do
          end do
          end do
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(forward_plan,tau13_term,tau13_term)
          call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
          call fftw_mpi_execute_dft(forward_plan,tau23_term,tau23_term)
          call fftw_mpi_execute_dft(forward_plan,tau31_term,tau31_term)
          call fftw_mpi_execute_dft(forward_plan,tau32_term,tau32_term)
          call fftw_mpi_execute_dft(forward_plan,tau33_term,tau33_term)
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            tau11_term(i,j,k) = tau11_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau12_term(i,j,k) = tau12_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau13_term(i,j,k) = tau13_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau21_term(i,j,k) = tau21_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau22_term(i,j,k) = tau22_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau23_term(i,j,k) = tau23_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau31_term(i,j,k) = tau31_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau32_term(i,j,k) = tau32_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau33_term(i,j,k) = tau33_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
          enddo
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(backward_plan,tau13_term,tau13_term)
          call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
          call fftw_mpi_execute_dft(backward_plan,tau23_term,tau23_term)
          call fftw_mpi_execute_dft(backward_plan,tau31_term,tau31_term)
          call fftw_mpi_execute_dft(backward_plan,tau32_term,tau32_term)
          call fftw_mpi_execute_dft(backward_plan,tau33_term,tau33_term)
          !
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            tau_bis(1,1,i,j,k) = tau_bis(1,1,i,j,k) + dreal(tau11_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(1,2,i,j,k) = tau_bis(1,2,i,j,k) + dreal(tau12_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(1,3,i,j,k) = tau_bis(1,3,i,j,k) + dreal(tau13_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,1,i,j,k) = tau_bis(2,1,i,j,k) + dreal(tau21_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,2,i,j,k) = tau_bis(2,2,i,j,k) + dreal(tau22_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,3,i,j,k) = tau_bis(2,3,i,j,k) + dreal(tau23_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,1,i,j,k) = tau_bis(3,1,i,j,k) + dreal(tau31_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,2,i,j,k) = tau_bis(3,2,i,j,k) + dreal(tau32_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,3,i,j,k) = tau_bis(3,3,i,j,k) + dreal(tau33_term(i,j,k)) * dl_alpha(m,n)
            !
          enddo
          enddo
          enddo
          !
        enddo ! loop of integral (alpha)
        !
        do p=1,3
        do q=1,3
          erroravg(p,q)=0.d0
          errormax(p,q)=0.d0
          errorgtr10(p,q)=0.d0
          errorgtr100(p,q)=0.d0
        enddo
        enddo
        errornorm2max=0.d0
        errornorm2avg=0.d0
        errornorm2gtr10=0.d0
        errornorm2gtr100=0.d0
        !
        ! Output comparaison results
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          norm2 = 0.d0
          norm2bis = 0.d0
          do p=1,3
          do q=1,3
            if(abs(tau(p,q,i,j,k))>2.d-4)then
              result = 2*abs(tau_bis(p,q,i,j,k)-tau(p,q,i,j,k))/(abs(tau_bis(p,q,i,j,k))+abs(tau(p,q,i,j,k)))
              norm2 = norm2 + tau(p,q,i,j,k)**2
              norm2bis = norm2bis + tau_bis(p,q,i,j,k)**2
              errormax(p,q) = max(errormax(p,q),result)
              erroravg(p,q) = erroravg(p,q) + result
              if(result > 0.1)then
                errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
              endif
              if(result > 1)then
                errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
                print *, 'p',p,'q',q,tau_bis(p,q,i,j,k),tau(p,q,i,j,k),result, '*'
              endif
            endif
          enddo
          enddo
          result = 2*abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
          errornorm2max = max(errornorm2max,result)
          errornorm2avg = errornorm2avg + result
          if(result > 0.1)then
            errornorm2gtr10 = errornorm2gtr10 + 1.d0
          endif
          if(result > 1)then
            errornorm2gtr100 = errornorm2gtr100 + 1.d0
            print *, norm2, norm2bis, result, '*'
          endif
          !
        enddo
        enddo
        enddo
        !
        !
        do p=1,3
        do q=1,3
          errormax(p,q) = pmax(errormax(p,q))
          erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja*ka)
          errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja*ka)
          errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja*ka)
        enddo
        enddo
        errornorm2max = pmax(errornorm2max)
        errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja*ka)
        errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja*ka)
        errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja*ka)
        !
        if(mpirank==0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
          endif
          !
          if(mpirank == 0)then
            open(fh,file=outfilename,form='formatted')
            write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
            write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
            write(fh,"(A8,1x,10(A20,1x))")'type','tau11','tau12','tau13','tau21','tau22','tau23','tau31','tau32','tau33','norm2'
            write(fh,"(A8,1x,10(E20.13E2,1x))")'max',((errormax(p,q),q=1,3),p=1,3),errornorm2max
            write(fh,"(A8,1x,10(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,3),p=1,3),errornorm2avg
            write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,3),p=1,3),errornorm2gtr10
            write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,3),p=1,3),errornorm2gtr100
            close(fh)
            print *, '>>>>', outfilename
          endif
          !
          !
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(loutput)then
          !
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
          else
            outfilename = 'pp/SGS_stress_'//mname//'.h5'
          endif
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          do p=1,3
            do q= 1,3
              write (termname, "(A3,I1,I1)") "tau",p,q
              call h5write(var=tau(p,q,1:im,1:jm,1:km),      varname=termname,    mode = modeio) 
              write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
              call h5write(var=tau_bis(p,q,1:im,1:jm,1:km),  varname=termname,    mode = modeio)
            enddo
          enddo
          !
          call h5io_end
          !
        endif
        !
      enddo ! loop of filter point l
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhol)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w1w3)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w2w3)
      call fftw_free(c_w3w1)
      call fftw_free(c_w3w2)
      call fftw_free(c_w3w3)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w1w3_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_w2w3_filted)
      call fftw_free(c_w3w1_filted)
      call fftw_free(c_w3w2_filted)
      call fftw_free(c_w3w3_filted)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A13)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_A23)
      call fftw_free(c_A31)
      call fftw_free(c_A32)
      call fftw_free(c_A33)
      call fftw_free(c_tau11_term)
      call fftw_free(c_tau12_term)
      call fftw_free(c_tau13_term)
      call fftw_free(c_tau21_term)
      call fftw_free(c_tau22_term)
      call fftw_free(c_tau23_term)
      call fftw_free(c_tau31_term)
      call fftw_free(c_tau32_term)
      call fftw_free(c_tau33_term)
      call mpistop
      deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
      !
    end subroutine SGSstress3D
    !
    subroutine SGST3D(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho, prs
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom,pcom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: sqrtalphas,dalphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max,ratio_min
      real(8) :: Galpha
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted,p_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p11_filted,p12_filted,p13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p21_filted,p22_filted,p23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p31_filted,p32_filted,p33_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p1_filted,p2_filted,p3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: rho1_filted,rho2_filted,rho3_filted
      complex(8), allocatable, dimension(:,:,:) :: theta_filted,M11_filted,M22_filted,M33_filted
      complex(8), allocatable, dimension(:,:,:) :: M12_filted,M13_filted,M21_filted
      complex(8), allocatable, dimension(:,:,:) :: M23_filted,M31_filted,M32_filted
      real(8) :: Es,Ec,Ts,Tc,Ps1,Pc1,Pc2
      !
      complex(8) :: vxr_D
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,c_pcom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted,c_p_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      type(C_PTR) :: c_p11_filted,c_p12_filted,c_p13_filted
      type(C_PTR) :: c_p21_filted,c_p22_filted,c_p23_filted
      type(C_PTR) :: c_p31_filted,c_p32_filted,c_p33_filted
      type(C_PTR) :: c_p1_filted,c_p2_filted,c_p3_filted
      type(C_PTR) :: c_rho1_filted,c_rho2_filted,c_rho3_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km),prs(0:im,0:jm,0:km))
      !
      if (thefilenumb .ne. 0) then
        write(stepname,'(i4.4)')thefilenumb
        infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
      else
        infilename='outdat/flowfield.'//modeio//'5'
      endif
      !
      call h5io_init(filename=infilename,mode='read')
      !
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='p',  var=prs(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      c_pcom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_pcom, pcom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T)
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T)
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T)
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T)
        pcom(i,j,k)=CMPLX(prs(i,j,k),0.d0,C_INTPTR_T)
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      call fftw_mpi_execute_dft(forward_plan,pcom,pcom)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        pcom(i,j,k)=pcom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      do k = 1,km
      do j = 1,jm
      do i = 1,im
        !
        if(im .ne. ia)then
          stop "error! im /= ia"
        endif
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
        else if(i<=(ia)) then
          k2(i,j,k) = real(j-ja-1,8)
        else
          print *,"Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
          k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=(ka)) then
          k3(i,j,k) = real(k+k0-ka-1,8)
        else
          print *,"Error, no wave number possible, (k+k0) must smaller than ja-1 !"
        end if
        !
      end do
      end do
      end do
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(sqrtalphas(num_alpha),dalphas(num_alpha))
      !
      do i=1,num_alpha
        sqrtalphas(i) = sqrt( exp(log(ratio_max**2) * (i-1) / (num_alpha-1)) ) * l_min
      enddo
      !
      dalphas(1) = sqrtalphas(1)**2 
      !
      do i=2,num_alpha
        dalphas(i) = sqrtalphas(i)**2 - sqrtalphas(i-1)**2 
      enddo
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !!!!
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      c_p_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p_filted, p_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      c_p11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p11_filted, p11_filted,[imfftw,jmfftw,kmfftw])
      c_p12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p12_filted, p12_filted,[imfftw,jmfftw,kmfftw])
      c_p13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p13_filted, p13_filted,[imfftw,jmfftw,kmfftw])
      c_p21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p21_filted, p21_filted,[imfftw,jmfftw,kmfftw])
      c_p22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p22_filted, p22_filted,[imfftw,jmfftw,kmfftw])
      c_p23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p23_filted, p23_filted,[imfftw,jmfftw,kmfftw])
      c_p31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p31_filted, p31_filted,[imfftw,jmfftw,kmfftw])
      c_p32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p32_filted, p32_filted,[imfftw,jmfftw,kmfftw])
      c_p33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p33_filted, p33_filted,[imfftw,jmfftw,kmfftw])
      c_p1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p1_filted, p1_filted,[imfftw,jmfftw,kmfftw])
      c_p2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p2_filted, p2_filted,[imfftw,jmfftw,kmfftw])
      c_p3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p3_filted, p3_filted,[imfftw,jmfftw,kmfftw])
      c_rho1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho1_filted, rho1_filted,[imfftw,jmfftw,kmfftw])
      c_rho2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho2_filted, rho2_filted,[imfftw,jmfftw,kmfftw])
      c_rho3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho3_filted, rho3_filted,[imfftw,jmfftw,kmfftw])
      !
      !
      allocate(theta_filted(1:im,1:jm,1:km),&
      M11_filted(1:im,1:jm,1:km),M22_filted(1:im,1:jm,1:km),M33_filted(1:im,1:jm,1:km),&
      M12_filted(1:im,1:jm,1:km),M13_filted(1:im,1:jm,1:km),M21_filted(1:im,1:jm,1:km),&
      M23_filted(1:im,1:jm,1:km),M31_filted(1:im,1:jm,1:km),M32_filted(1:im,1:jm,1:km))
      !
      !
      !
      Es = 0.d0
      Ec = 0.d0
      Ts = 0.d0
      Tc = 0.d0
      Ps1 = 0.d0
      Pc1 = 0.d0
      Pc2 = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do n=1,num_alpha
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alpha,',now is ',&
                                value(5), ':', value(6),':',value(7)
        !
        !!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do i=1,im
        do j=1,jm
        do k=1,km
          Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*sqrtalphas(n)**2/2.d0) ! Filtre scale :sqrtalpha
          w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
          w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
          w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
          rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
          p_filted(i,j,k)   = pcom(i,j,k)  *Galpha
          !
          rho1_filted(i,j,k) = imag*rho_filted(i,j,k)*k1(i,j,k)
          rho2_filted(i,j,k) = imag*rho_filted(i,j,k)*k2(i,j,k)
          rho3_filted(i,j,k) = imag*rho_filted(i,j,k)*k3(i,j,k)
          !
          p11_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k1(i,j,k)
          p21_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k1(i,j,k)
          p31_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k1(i,j,k)
          p12_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k2(i,j,k)
          p22_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k2(i,j,k)
          p32_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k2(i,j,k)
          p13_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k3(i,j,k)
          p23_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k3(i,j,k)
          p33_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k3(i,j,k)
          p1_filted(i,j,k)   = imag*p_filted(i,j,k)*k1(i,j,k)
          p2_filted(i,j,k)   = imag*p_filted(i,j,k)*k2(i,j,k)
          p3_filted(i,j,k)   = imag*p_filted(i,j,k)*k3(i,j,k)
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        call fftw_mpi_execute_dft(backward_plan,p_filted,p_filted)
        !
        call fftw_mpi_execute_dft(backward_plan,p11_filted,p11_filted)
        call fftw_mpi_execute_dft(backward_plan,p21_filted,p21_filted)
        call fftw_mpi_execute_dft(backward_plan,p31_filted,p31_filted)
        call fftw_mpi_execute_dft(backward_plan,p12_filted,p12_filted)
        call fftw_mpi_execute_dft(backward_plan,p22_filted,p22_filted)
        call fftw_mpi_execute_dft(backward_plan,p32_filted,p32_filted)
        call fftw_mpi_execute_dft(backward_plan,p13_filted,p13_filted)
        call fftw_mpi_execute_dft(backward_plan,p23_filted,p23_filted)
        call fftw_mpi_execute_dft(backward_plan,p33_filted,p33_filted)
        call fftw_mpi_execute_dft(backward_plan,p1_filted,p1_filted)
        call fftw_mpi_execute_dft(backward_plan,p2_filted,p2_filted)
        call fftw_mpi_execute_dft(backward_plan,p3_filted,p3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho1_filted,rho1_filted)
        call fftw_mpi_execute_dft(backward_plan,rho2_filted,rho2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho3_filted,rho3_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do i=1,im
        do j=1,jm
        do k=1,km
          !
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)   = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)   = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)   = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          theta_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
          !
          M11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          M22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          M33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          !
          M12_filted(i,j,k) = A12_filted(i,j,k)
          M13_filted(i,j,k) = A13_filted(i,j,k)
          M21_filted(i,j,k) = A21_filted(i,j,k)
          M23_filted(i,j,k) = A23_filted(i,j,k)
          M31_filted(i,j,k) = A31_filted(i,j,k)
          M32_filted(i,j,k) = A32_filted(i,j,k)
          !
        end do
        end do
        end do
        !
        !!!! T terms
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Es = Es + dreal(0.5d0 * rho_filted(i,j,k)*&
          (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + &
          M13_filted(i,j,k)*M13_filted(i,j,k) + M21_filted(i,j,k)*M21_filted(i,j,k) + &
          M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + &
          M33_filted(i,j,k)*M33_filted(i,j,k))*dalphas(n))
          !
          Ec = Ec + dreal(1.d0/6.d0 * rho_filted(i,j,k)* theta_filted(i,j,k)*theta_filted(i,j,k) *dalphas(n))
          !
          Ts = Ts + dreal(rho_filted(i,j,k)*&
          (M11_filted(i,j,k)*M11_filted(i,j,k)*M11_filted(i,j,k)+&
          M12_filted(i,j,k)*M21_filted(i,j,k)*M11_filted(i,j,k)+&
          M13_filted(i,j,k)*M31_filted(i,j,k)*M11_filted(i,j,k)+&
          M11_filted(i,j,k)*M12_filted(i,j,k)*M12_filted(i,j,k)+&
          M12_filted(i,j,k)*M22_filted(i,j,k)*M12_filted(i,j,k)+&
          M13_filted(i,j,k)*M32_filted(i,j,k)*M12_filted(i,j,k)+&
          M11_filted(i,j,k)*M13_filted(i,j,k)*M13_filted(i,j,k)+&
          M12_filted(i,j,k)*M23_filted(i,j,k)*M13_filted(i,j,k)+&
          M13_filted(i,j,k)*M33_filted(i,j,k)*M13_filted(i,j,k)+&
          M21_filted(i,j,k)*M11_filted(i,j,k)*M21_filted(i,j,k)+&
          M22_filted(i,j,k)*M21_filted(i,j,k)*M21_filted(i,j,k)+&
          M23_filted(i,j,k)*M31_filted(i,j,k)*M21_filted(i,j,k)+&
          M21_filted(i,j,k)*M12_filted(i,j,k)*M22_filted(i,j,k)+&
          M22_filted(i,j,k)*M22_filted(i,j,k)*M22_filted(i,j,k)+&
          M23_filted(i,j,k)*M32_filted(i,j,k)*M22_filted(i,j,k)+&
          M21_filted(i,j,k)*M13_filted(i,j,k)*M23_filted(i,j,k)+&
          M22_filted(i,j,k)*M23_filted(i,j,k)*M23_filted(i,j,k)+&
          M23_filted(i,j,k)*M33_filted(i,j,k)*M23_filted(i,j,k)+&
          M31_filted(i,j,k)*M11_filted(i,j,k)*M31_filted(i,j,k)+&
          M32_filted(i,j,k)*M21_filted(i,j,k)*M31_filted(i,j,k)+&
          M33_filted(i,j,k)*M31_filted(i,j,k)*M31_filted(i,j,k)+&
          M31_filted(i,j,k)*M12_filted(i,j,k)*M32_filted(i,j,k)+&
          M32_filted(i,j,k)*M22_filted(i,j,k)*M32_filted(i,j,k)+&
          M33_filted(i,j,k)*M32_filted(i,j,k)*M32_filted(i,j,k)+&
          M31_filted(i,j,k)*M13_filted(i,j,k)*M33_filted(i,j,k)+&
          M32_filted(i,j,k)*M23_filted(i,j,k)*M33_filted(i,j,k)+&
          M33_filted(i,j,k)*M33_filted(i,j,k)*M33_filted(i,j,k)+&
          2.d0/3.d0 * &
          (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k)+ &
          M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
          *theta_filted(i,j,k))*dalphas(n))
          !
          Tc = Tc + dreal(rho_filted(i,j,k)*&
          ((M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k) + &
          M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
          *theta_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)**3)*dalphas(n))
          !
          Ps1 = Ps1 + dreal((&
          (p11_filted(i,j,k)*A11_filted(i,j,k)+p12_filted(i,j,k)*A12_filted(i,j,k)+p13_filted(i,j,k)*A13_filted(i,j,k) &
          +p21_filted(i,j,k)*A21_filted(i,j,k)+p22_filted(i,j,k)*A22_filted(i,j,k)+p23_filted(i,j,k)*A23_filted(i,j,k) &
          +p31_filted(i,j,k)*A31_filted(i,j,k)+p32_filted(i,j,k)*A32_filted(i,j,k)+p33_filted(i,j,k)*A33_filted(i,j,k))& 
          - 1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * theta_filted(i,j,k)&
          )*dalphas(n))
          !
          Pc1 = Pc1 + dreal(1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * &
                theta_filted(i,j,k)*dalphas(n))
          !
          Pc2 = Pc2 + dreal(1.d0/3.d0 / rho_filted(i,j,k) * (rho1_filted(i,j,k)*p1_filted(i,j,k) + &
          rho2_filted(i,j,k)*p2_filted(i,j,k) + rho3_filted(i,j,k)*p3_filted(i,j,k))* theta_filted(i,j,k) *dalphas(n))
          !
        enddo
        enddo
        enddo
        !
      enddo
      !
      Es = psum(Es)/(1.d0*ia*ja*ka)
      Ec = psum(Ec)/(1.d0*ia*ja*ka)
      Ts = psum(Ts)/(1.d0*ia*ja*ka)
      Tc = psum(Tc)/(1.d0*ia*ja*ka)
      Ps1 = psum(Ps1)/(1.d0*ia*ja*ka)
      Pc1 = psum(Pc1)/(1.d0*ia*ja*ka)
      Pc2 = psum(Pc2)/(1.d0*ia*ja*ka)
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_T_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_T.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time Es Ec Ts Tc Ps1 Pc1 Pc2')
        call listwrite(hand_a,Es,Ec,Ts,Tc,Ps1,Pc1,Pc2)
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call fftw_free(c_p11_filted)
      call fftw_free(c_p12_filted)
      call fftw_free(c_p13_filted)
      call fftw_free(c_p21_filted)
      call fftw_free(c_p22_filted)
      call fftw_free(c_p23_filted)
      call fftw_free(c_p31_filted)
      call fftw_free(c_p32_filted)
      call fftw_free(c_p33_filted)
      call fftw_free(c_p1_filted)
      call fftw_free(c_p2_filted)
      call fftw_free(c_p3_filted)
      call fftw_free(c_rho1_filted)
      call fftw_free(c_rho2_filted)
      call fftw_free(c_rho3_filted)
      call mpistop
      deallocate(theta_filted,M11_filted,M22_filted,M33_filted)
      deallocate(M12_filted,M13_filted,M21_filted)
      deallocate(M23_filted,M31_filted,M32_filted)
      deallocate(k1,k2,k3)
      deallocate(sqrtalphas,dalphas)
      !
    end subroutine SGST3D
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This function is used to read the SGS postprocess 
    ! calculation file at datin/SGSinput
    ! Outputs: 
    !   num_l        : number of ell
    !   num_alpha    : number of maximum alpha for integration
    !   num_alphamin : number of minimum alpha for integration
    !   ratio_max    : maximum ell/lmin 
    !   ratio_min    : minimum ell/lmin
    !   loutput      : logical, output stress hdf5 file or not
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      !
      use parallel,only: bcast,mpirank
      !
      ! local data
      integer, intent(out) :: num_l,num_alpha,num_alphamin
      real(8), intent(out) :: ratio_max,ratio_min
      logical, intent(out) :: loutput
      character(len=64) :: inputfile
      integer :: fh
      !
      inputfile='datin/SGSinput'
      !
      if(mpirank==0) then
        !
        fh=get_unit()
        !
        open(fh,file=trim(inputfile),action='read')
        read(fh,'(//)')
        read(fh,*)num_l,num_alpha,num_alphamin
        read(fh,'(/)')
        read(fh,*)ratio_max,ratio_min
        read(fh,'(/)')
        read(fh,*)loutput
        close(fh)
        print*,' >> ',trim(inputfile),' ... done'
        print*,' >>> Get: Number of l is',num_l,'Number of alpha is',num_alpha,'Minimum number of alpha is',num_alphamin
        print*,' >>> Ratio max:',ratio_max,'Ratio min',ratio_min
        if(loutput)then
          print *, ' >>> No output stress'
        else
          print *, ' >>> Output stress'
        endif
        !
      endif
      !
      call bcast(num_l)
      call bcast(num_alpha)
      call bcast(num_alphamin)
      call bcast(ratio_max)
      call bcast(ratio_min)
      call bcast(loutput)
      !
    end subroutine readSGSinput
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This function is used to generate ell and alpha used for integration 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,&
                                num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      ! 
      integer, intent(in) :: num_l
      real(8), intent(in) :: l_min
      real(8), intent(in) :: ratio_max,ratio_min
      real(8), intent(out) :: l_lim(:)
      integer, intent(in),optional :: num_alpha,num_alphamin
      real(8), intent(out),optional :: l_sqrtalpha(:,:),l_phi(:,:),dl_alpha(:,:)
      integer, intent(out),optional :: num_alphas(:)
      integer :: i,j
      !
      do i=1,num_l
          l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
          if(present(num_alpha)) then
              num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
              !
              do j=1,num_alphas(i)
                  l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
                  (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
                  l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
              enddo
          endif
      enddo
      !
      if(present(num_alpha)) then
          do i=1,num_l
              dl_alpha(i,1) = l_sqrtalpha(i,1)**2 
              !
              do j=2,num_alphas(i)
                dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
              enddo
            enddo
      endif
    end subroutine SGSscale_allocate
    !
end module udf_pp_SGS