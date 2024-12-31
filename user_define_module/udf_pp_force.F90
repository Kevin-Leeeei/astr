!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| initial field generating.                                           |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_force
    !
    use constdef
    use stlaio,  only: get_unit
    use udf_tool
    !
    implicit none
    !
    contains
    !
    subroutine ppForceentrance
        !
        use cmdefne
        use parallel,        only : mpirank,bcast,mpisize,lio
        !
        ! local data
        character(len=64) ::  readmode,inputfile
        integer :: filenumb
        !
        !
        if(mpirank == 0) then
          call readkeyboad(readmode) 
        endif
        call bcast(readmode)
        !
        if(trim(readmode)=='spectral') then
            !
            if(mpirank == 0) then
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
            endif
            call bcast(filenumb)
            call force_spectral(filenumb)
            !
        elseif(trim(readmode)=='selfmade') then
            !
            if(mpirank == 0) then
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
            endif
            call bcast(filenumb)
            call force_selfmade(filenumb)
            !
        elseif(trim(readmode)=='newspectral') then
            !
            if(mpirank == 0) then
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
            endif
            call bcast(filenumb)
            call force_newspectral(filenumb)
            !
        elseif(trim(readmode)=='newspectraltest') then
            !
            if(mpirank == 0) then
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
            endif
            call bcast(filenumb)
            call force_newspectraltest(filenumb)
            !
        else
            print* ,"Readmode is not defined!", readmode
        endif
    end subroutine ppForceentrance
    !
    !
    subroutine force_spectral(thefilenumb)
        !
        use, intrinsic :: iso_c_binding
        use readwrite, only : readinput, readic
        use fftwlink
        use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,time,&
                         nstep
        use commarray, only : vel, rho, forcep
        use hdf5io
        use solver,    only : refcal
        use utility,   only : listinit,listwrite
        use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop, dataswap
        include 'fftw3-mpi.f03'
        !
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=1) :: modeio
        integer :: allkmax
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
        complex(8), allocatable, dimension(:,:) :: usspe,udspe
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1d,u2d,u1s,u2s
        character(len=128) :: outfilename
        type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u1d, c_u2d, c_u1s, c_u2s
        real(8), allocatable, dimension(:,:) :: k1,k2
        real(8) :: Ed, Es, dk, kk, alpha
        integer :: i,j,n
        integer :: values(8)
        integer :: forcek=5
        real(8) :: forcespes=0.0035d0
        ! 
        call readinput
        !
        call readic
        !
        if(mpirank==0) print*, "Warning, this is only a test and fix forcek=", forcek, "forcespes=", forcespes
        call refcal
        if(mpirank==0)  print*, '** refcal done!'
        !
        !
        modeio='h'
        !
        dk = 1.d0
        allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
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
        allocate(vel(0:im,0:jm,0:km,1:2))
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
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *, 'Begin!' ,values(5),":",values(6),":",values(7),":",values(8)
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
            u1spe(i,j)=CMPLX(vel(i,j,0,1),0.d0,C_INTPTR_T);
            u2spe(i,j)=CMPLX(vel(i,j,0,2),0.d0,C_INTPTR_T);
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
        !
        Ed = 0.d0
        Es = 0.d0
        do j=1,jm
        do i=1,im
            kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
            if(kint(kk,dk,2,1.d0)==forcek)then
                usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
                udspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
                u1d(i,j)=  udspe(i,j)*k1(i,j)/kk
                u2d(i,j)=  udspe(i,j)*k2(i,j)/kk
                u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
                u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
                Es = Es + usspe(i,j)*conjg(usspe(i,j))/2
                Ed = Ed + udspe(i,j)*conjg(udspe(i,j))/2
            else
                usspe(i,j) = 0.d0
                udspe(i,j) = 0.d0
                u1d(i,j)= 0.d0
                u2d(i,j)= 0.d0
                u1s(i,j)= 0.d0
                u2s(i,j)= 0.d0
            endif
            !
        end do
        end do
        !
        call fftw_mpi_execute_dft(backward_plan,u1d,u1d)
        call fftw_mpi_execute_dft(backward_plan,u2d,u2d)
        call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
        call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
        !
        Es = psum(Es)
        Ed = psum(Ed)
        !
        !
        alpha = sqrt((forcespes - Ed)/Es)
        !
        allocate(forcep(0:im,0:jm,0:km,1:3))
        do j=1,jm
        do i=1,im
            forcep(i,j,0,1) = alpha * real(u1s(i,j))
            forcep(i,j,0,2) = alpha * real(u2s(i,j))
        enddo
        enddo
        !
        forcep(0:(im-1),1:jm,0,1)=forcep(1:im,1:jm,0,1)
        forcep(0:(im-1),1:jm,0,2)=forcep(1:im,1:jm,0,2)
        !
        forcep(0:im,0:(jm-1),0,1)=forcep(0:im,1:jm,0,1)
        forcep(0:im,0:(jm-1),0,2)=forcep(0:im,1:jm,0,2)
        forcep(0:im,0:jm,0,3) = 0.d0
        !
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *,'Finish!' , values(5),":",values(6),":",values(7),":",values(8)
        !
        if (thefilenumb .ne. 0) then
            outfilename = 'pp/forcespectra'//stepname//'.h5'
        else
            outfilename = 'pp/forcespectra.h5'
        endif
        call h5io_init(trim(outfilename),mode='write')
        call h5wa2d_r8(varname='f1',var=forcep(0:im,0:jm,0,1), dir='k')
        call h5wa2d_r8(varname='f2',var=forcep(0:im,0:jm,0,2),dir='k')
        call h5io_end
        !
        if(mpirank==0) print *, '>>> alpha = ', alpha, forcespes, Ed, Es
        !
        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        call fftw_mpi_cleanup()
        call fftw_free(c_u1spe)
        call fftw_free(c_u2spe)
        call fftw_free(c_u1d)
        call fftw_free(c_u1s)
        call fftw_free(c_u2d)
        call fftw_free(c_u2s)
        call mpistop
        !
        deallocate(k1,k2,usspe,udspe)
        !
    end subroutine force_spectral
    !
    ! subroutine mympitest
    !     use readwrite, only : readinput
    !     use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,time,&
    !                      nstep,forcek,forceamp
    !     use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop, dataswap
    !     ! 
    !     call readinput
    !     !
    !     call mpisizedis
    !     if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !     !
    !     call parapp
    !     if(mpirank==0) print*, '** parapp done!'
    !     !
    !     call parallelini
    !     if(mpirank==0)  print*, '** parallelini done!'
    !     !
    !     call refcal
    !     if(mpirank==0)  print*, '** refcal done!'
    !     !
    !     allocate(vel(0:im,0:jm,0,1:2))
    !     !
    !     do i=0,im
    !     do j=0,jm
    !         vel(i,j,0,1) = mpirank
    !         vel(i,j,0,2) = 100-mpirank
    !     enddo
    !     enddo
    ! end subroutine mympitest
    subroutine force_newspectral(thefilenumb)
        !
        use, intrinsic :: iso_c_binding
        use readwrite, only : readinput, readic
        use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,time,&
                         nstep
        use commarray, only : vel, rho, forcep
        use hdf5io
        use solver,    only : refcal
        use utility,   only : listinit,listwrite
        use parallel!,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop, dataswap, ig0
        use fftwlink
        include 'fftw3-mpi.f03'
        !
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=1) :: modeio
        integer :: allkmax
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
        complex(8), allocatable, dimension(:,:) :: usspe,udspe
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1d,u2d,u1s,u2s
        character(len=128) :: outfilename
        type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u1d, c_u2d, c_u1s, c_u2s
        real(8), allocatable, dimension(:,:) :: k1,k2
        real(8), allocatable, dimension(:,:) :: localvel1t, localvel2t, force1t, force2t
        real(8), allocatable, dimension(:,:) :: fftvel1, fftvel2, fftforce1, fftforce2
        real(8) :: Ed, Es, dk, kk, alpha
        integer :: i,j,n
        integer :: values(8)
        integer :: forcek=5
        real(8) :: forcespes=0.0035d0
        ! 
        call readinput
        !
        call readic
        !
        if(mpirank==0) print*, "Warning, this is only a test and fix forcek=", forcek, "forcespes=", forcespes
        !
        modeio='h'
        !
        dk = 1.d0
        allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
        !
        !
        call mpisizedis
        if(mpirank==0)  print*, '** mpisizedis & parapp done!'
        !
        call parapp
        if(mpirank==0) print*, '** parapp done!'
        !
        call parallelini
        if(mpirank==0)  print*, '** parallelini done!'
        !
        call refcal
        if(mpirank==0)  print*, '** refcal done!'
        !
        allocate(vel(0:im,0:jm,0:km,1:2))
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
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0)  print *, "Field read finish!"
        !
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *, 'Begin!' ,values(5),":",values(6),":",values(7),":",values(8)
        !
        !
        call fftwprepare_forcing
        ! !
        allocate(localvel1t(1:jm,1:im),localvel2t(1:jm,1:im))
        allocate(force1t(1:jm,1:im),force2t(1:jm,1:im))
        allocate(fftvel1(1:ia,1:jmf),fftvel2(1:ia,1:jmf))
        allocate(fftforce1(1:ia,1:jmf),fftforce2(1:ia,1:jmf))
        !
        do j=1,jm
        do i=1,im
            localvel1t(j,i)=vel(i,j,0,1)
            localvel2t(j,i)=vel(i,j,0,2)
            !
        enddo
        enddo
        !
        fftvel1 = fftw_grid_fence(localvel1t)
        fftvel2 = fftw_grid_fence(localvel2t)
        !
        !
        ! Begin FFTW
        !
        c_u1spe = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u1spe, u1spe, [iafftw,jmfftw])
        c_u2spe = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u2spe, u2spe, [iafftw,jmfftw])
        !
        !!!! Do S-C decomposition
        c_u1d = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u1d, u1d, [iafftw,jmfftw])
        c_u2d = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u2d, u2d, [iafftw,jmfftw])
        c_u1s = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u1s, u1s, [iafftw,jmfftw])
        c_u2s = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u2s, u2s, [iafftw,jmfftw])
        ! planning
        forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
        backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        !
        if(mpirank==0)  print *, 'Prepare fft finish'
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
        if(mpirank==0)  print *, 'Before FFT forward'
        !
        !!!! Do 2d FFT
        call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
        call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
        !
        if(mpirank==0)  print *, 'After FFT forward'
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
        Ed = 0.d0
        Es = 0.d0
        !
        do j=1,jmf
        do i=1,ia
            kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
            if(kint(kk,dk,2,1.d0)==forcek)then
                usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
                udspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
                u1d(i,j)=  udspe(i,j)*k1(i,j)/kk
                u2d(i,j)=  udspe(i,j)*k2(i,j)/kk
                u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
                u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
                Es = Es + usspe(i,j)*conjg(usspe(i,j))/2
                Ed = Ed + udspe(i,j)*conjg(udspe(i,j))/2
            else
                usspe(i,j) = 0.d0
                udspe(i,j) = 0.d0
                u1d(i,j)= 0.d0
                u2d(i,j)= 0.d0
                u1s(i,j)= 0.d0
                u2s(i,j)= 0.d0
            endif
            !
        end do
        end do
        !
        if(mpirank==0)  print *, 'Before FFT backward'
        !
        call fftw_mpi_execute_dft(backward_plan,u1d,u1d)
        call fftw_mpi_execute_dft(backward_plan,u2d,u2d)
        call fftw_mpi_execute_dft(backward_plan,u1s,u1s)
        call fftw_mpi_execute_dft(backward_plan,u2s,u2s)
        !
        if(mpirank==0)  print *, 'After FFT backward'
        !
        Es = psum(Es)
        Ed = psum(Ed)
        alpha = sqrt(abs(forcespes - Ed)/Es)
        !
        !
        do j=1,jmf
        do i=1,ia
            fftforce1(i,j) = alpha * real(u1s(i,j))
            fftforce2(i,j) = alpha * real(u2s(i,j))
        enddo
        enddo
        !
        !
        force1t = fftw_fence_grid(fftforce1)
        force2t = fftw_fence_grid(fftforce2)
        !
        if(mpirank==0)  print *, 'MPI scatter finish'
        !
        !
        allocate(forcep(0:im,0:jm,0:km,1:3))
        !
        do j=1,jm
        do i=1,im
            forcep(i,j,0,1) = force1t(j,i)
            forcep(i,j,0,2) = force2t(j,i)
        enddo
        enddo
        !
        forcep(0:(im-1),1:jm,0,1)=forcep(1:im,1:jm,0,1)
        forcep(0:(im-1),1:jm,0,2)=forcep(1:im,1:jm,0,2)
        !
        forcep(0:im,0:(jm-1),0,1)=forcep(0:im,1:jm,0,1)
        forcep(0:im,0:(jm-1),0,2)=forcep(0:im,1:jm,0,2)
        forcep(0:im,0:jm,0,3) = 0.d0
        !
        call date_and_time(values=values)
        if(mpirank==0) print *,'Finish!' , values(5),":",values(6),":",values(7),":",values(8)
        !
        if(mpirank==0) print *,'>>> alpha = ', alpha, forcespes, Ed, Es
        !
        if (thefilenumb .ne. 0) then
            outfilename = 'pp/forcenewspectra'//stepname//'.h5'
        else
            outfilename = 'pp/forcenewspectra.h5'
        endif
        call h5io_init(trim(outfilename),mode='write')
        call h5wa2d_r8(varname='f1',var=forcep(0:im,0:jm,0,1), dir='k')
        call h5wa2d_r8(varname='f2',var=forcep(0:im,0:jm,0,2),dir='k')
        call h5io_end
        !
        !
        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        call fftw_mpi_cleanup()
        call fftw_free(c_u1spe)
        call fftw_free(c_u2spe)
        call fftw_free(c_u1d)
        call fftw_free(c_u1s)
        call fftw_free(c_u2d)
        call fftw_free(c_u2s)
        call mpistop
        !
        deallocate(k1,k2,usspe,udspe)
        !
    end subroutine force_newspectral
    !
    subroutine force_newspectraltest(thefilenumb)
        !
        use readwrite, only : readinput, readic
        use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,time,nstep,forcenum
        use commarray, only : vel, rho, forcep
        use hdf5io
        use solver,    only : refcal
        use utility,   only : listinit,listwrite
        use parallel!,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop, dataswap, ig0
        use userdefine, only: udf_generate_force
        use fftwlink, only: fftwprepare_forcing
        !
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=1) :: modeio
        character(len=128) :: outfilename
        integer :: values(8), ierr
        real(8), allocatable,dimension(:) :: alphas, alphad
        !
        ! 
        call readinput
        !
        call readic
        !
        allocate(alphas(1:forcenum), alphad(1:forcenum))
        !
        modeio='h'
        !
        !
        call mpisizedis
        if(mpirank==0)  print*, '** mpisizedis & parapp done!'
        !
        call parapp
        if(mpirank==0) print*, '** parapp done!'
        !
        call parallelini
        if(mpirank==0)  print*, '** parallelini done!'
        !
        call refcal
        if(mpirank==0)  print*, '** refcal done!'
        !
        allocate(vel(0:im,0:jm,0:km,1:2))
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
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0)  print *, "Field read finish!"
        !
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *, 'Begin!' ,values(5),":",values(6),":",values(7),":",values(8)
        !
        call fftwprepare_forcing
        !
        allocate(forcep(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        !
        call udf_generate_force(alphas, alphad)
        !
        call date_and_time(values=values)
        if(mpirank==0)  print *,'Finish!' , values(5),":",values(6),":",values(7),":",values(8)
        !
        if(mpirank==0) print *, alphas(1), alphad(1)
        !
        !
        if (thefilenumb .ne. 0) then
            outfilename = 'pp/forcenewspectra'//stepname//'.h5'
        else
            outfilename = 'pp/forcenewspectra.h5'
        endif
        call h5io_init(trim(outfilename),mode='write')
        call h5write(varname='f1',var=forcep(0:im,0:jm,0:km,1),mode='h')
        call h5write(varname='f2',var=forcep(0:im,0:jm,0:km,2),mode='h')
        call h5io_end
        !
        call mpistop
        !
    end subroutine force_newspectraltest
    !
    subroutine force_selfmade(thefilenumb)
        !
        use readwrite, only : readinput, readic
        use commvar,  only : im,jm,km,hm,ndims,deltat,ia,ja,ka,time,&
                         nstep
        use commarray, only : vel, rho, x, forcep
        use hdf5io
        use geom,      only : geomcal
        use solver,    only : refcal
        use utility,   only : listinit,listwrite
        use parallel,  only : bcast, pmax, pmin, psum, lio, parallelini,mpistop, dataswap, parapp, mpisizedis
        use gridgeneration
        !
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        integer, dimension(1:2) :: sgns
        integer :: i,j,k1,k2, sgn1,sgn2, jlow, jtop, kkx, kky
        character(len=1) :: modeio
        real(8) :: u1hat_r,u1hat_i,u2hat_r,u2hat_i, ushat_r,ushat_i,udhat_r,udhat_i, &
        us1hat_r,us1hat_i,us2hat_r,us2hat_i,ud1hat_r,ud1hat_i,ud2hat_r,ud2hat_i, &
        Es, Ed, k, dk, alpha
        real(8), allocatable,dimension(:,:) :: us1,us2,ud1,ud2
        character(len=128) :: outfilename
        integer :: ierr
        integer :: values(8)
        integer :: forcek=5
        real(8) :: forcespes=0.0035d0
        !
        call readinput
        !
        call readic
        !
        if(mpirank==0) print*, "Warning, this is only a test and fix forcek=", forcek, "forcespes=", forcespes
        !
        call mpisizedis
        if(mpirank==0) print*, '** mpisizedis done!'
        !
        call parapp
        if(mpirank==0) print*, '** parapp done!'
        !
        !
        call parallelini
        if(mpirank==0)  print*, '** parallelini done!'
        !
        call refcal
        if(mpirank==0)  print*, '** refcal done!'
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        allocate(vel(0:im,0:jm,0:km,1:2))
        !
        call gridsquare(2.d0*pi,2.d0*pi)
        !
        call geomcal
        !
        modeio='h'
        !
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
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        !
        if(mpirank==0)  print *, "Field read finish!"
        !
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *, 'Begin!' ,values(5),":",values(6),":",values(7),":",values(8)
        !
        sgns(1) = -1
        sgns(2) = 1
        !
        dk = 0.5d0
        !
        allocate(us1(0:im,0:jm),us2(0:im,0:jm),ud1(0:im,0:jm),ud2(0:im,0:jm))
        do i = 1,im
        do j = 1,jm 
            us1(i,j) = 0.d0
            us2(i,j) = 0.d0
            ud1(i,j) = 0.d0
            ud2(i,j) = 0.d0
        enddo
        enddo
        !
        Es = 0.d0
        Ed = 0.d0
        !
        do k1 = 0 ,forcek
        jlow = floor(sqrt(max(real((forcek-dk) ** 2 - k1 ** 2),1e-10)))
        jtop = ceiling(sqrt(real((forcek+dk) ** 2 - k1 ** 2)))
        do k2 = jlow,jtop
        k=sqrt(real(k1**2+k2**2))
        if(kint(k,dk,2,0.d0)==forcek)then
        do sgn1 = 1,2
        do sgn2 = 1,2
        if((.not. (k1==0 .and. sgn1 == 2)) .and. (.not. (k2 == 0 .and. sgn2 ==2)))then 
            kkx = sgns(sgn1) * k1
            kky = sgns(sgn2) * k2
            u1hat_r = 0.d0
            u1hat_i = 0.d0
            u2hat_r = 0.d0
            u2hat_i = 0.d0
            do i = 1,im
            do j = 1,jm
                u1hat_r = u1hat_r + cos(-kkx*x(i,j,0,1) - kky*x(i,j,0,2)) * vel(i,j,0,1)
                u1hat_i = u1hat_i + sin(-kkx*x(i,j,0,1) - kky*x(i,j,0,2)) * vel(i,j,0,1)
                u2hat_r = u2hat_r + cos(-kkx*x(i,j,0,1) - kky*x(i,j,0,2)) * vel(i,j,0,2)
                u2hat_i = u2hat_i + sin(-kkx*x(i,j,0,1) - kky*x(i,j,0,2)) * vel(i,j,0,2)
            enddo
            enddo
            !
            u1hat_r = psum(u1hat_r)/(ia*ja)
            u1hat_i = psum(u1hat_i)/(ia*ja)
            u2hat_r = psum(u2hat_r)/(ia*ja)
            u2hat_i = psum(u2hat_i)/(ia*ja)
            !
            ushat_r = u1hat_r*kky/k - u2hat_r*kkx/k
            ushat_i = u1hat_i*kky/k - u2hat_i*kkx/k
            udhat_r = u1hat_r*kkx/k + u2hat_r*kky/k
            udhat_i = u1hat_i*kkx/k + u2hat_i*kky/k
            !
            us1hat_r = ushat_r*kky/k
            us1hat_i = ushat_i*kky/k
            us2hat_r = -ushat_r*kkx/k
            us2hat_i = -ushat_i*kkx/k
            !
            ud1hat_r = udhat_r*kkx/k
            ud1hat_i = udhat_i*kkx/k
            ud2hat_r = udhat_r*kky/k
            ud2hat_i = udhat_i*kky/k
            !
            Es = Es + (ushat_r**2 + ushat_i**2)/2
            Ed = Ed + (udhat_r**2 + udhat_i**2)/2
            !
            do i = 0,im
            do j = 0,jm
                us1(i,j) = us1(i,j) + us1hat_r * cos(kkx*x(i,j,0,1)+kky*x(i,j,0,2)) &
                        - us1hat_i * sin(kkx*x(i,j,0,1)+kky*x(i,j,0,2))
                us2(i,j) = us2(i,j) + us2hat_r * cos(kkx*x(i,j,0,1)+kky*x(i,j,0,2)) &
                        - us2hat_i * sin(kkx*x(i,j,0,1)+kky*x(i,j,0,2))
                ud1(i,j) = ud1(i,j) + ud1hat_r * cos(kkx*x(i,j,0,1)+kky*x(i,j,0,2)) &
                        - ud1hat_i * sin(kkx*x(i,j,0,1)+kky*x(i,j,0,2))
                ud2(i,j) = ud2(i,j) + ud2hat_r * cos(kkx*x(i,j,0,1)+kky*x(i,j,0,2)) &
                        - ud2hat_i * sin(kkx*x(i,j,0,1)+kky*x(i,j,0,2))
            enddo
            enddo
            !
        endif
        enddo
        enddo
        endif
        enddo
        enddo
        !

        !
        alpha = sqrt((forcespes - Ed)/Es)
        do j=0,jm
        do i=0,im
            forcep(i,j,0,1) = alpha * real(us1(i,j))
            forcep(i,j,0,2) = alpha * real(us2(i,j))
        enddo
        enddo
        !
        call mpi_barrier(mpi_comm_world,ierr)
        call date_and_time(values=values)
        if(mpirank==0)  print *, 'Finish!' ,values(5),":",values(6),":",values(7),":",values(8)
        !
        if (thefilenumb .ne. 0) then
            outfilename = 'pp/forceselfmade'//stepname//'.h5'
        else
            outfilename = 'pp/forceselfmade.h5'
        endif
        call h5io_init(trim(outfilename),mode='write')
        !
        call h5wa2d_r8(varname='f1',var=forcep(0:(im-1),0:(jm-1),0,1),dir='k')
        call h5wa2d_r8(varname='f2',var=forcep(0:(im-1),0:(jm-1),0,2),dir='k')
        call h5io_end
        !
        if(mpirank==0) print *, '>>> alpha = ', alpha, forcespes, Ed, Es
        !
        deallocate(us1,us2,ud1,ud2)
        !
    end subroutine force_selfmade
    !
end module udf_pp_force