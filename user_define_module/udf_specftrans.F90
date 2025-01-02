module udf_specftrans
  !
  !
  use constdef
  use stlaio,  only: get_unit
  !
  implicit none
  !
  contains 
  subroutine specftrans
    use cmdefne
    use, intrinsic :: iso_c_binding
    use readwrite, only : readinput
    use fftwlink
    use commvar,only : time,nstep,im,jm,km,ia,ja,ka
    use commarray, only: vel, rho, prs
    use hdf5io
    use solver,    only : refcal
    use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
    include 'fftw3-mpi.f03'
    !
    ! arguments
    integer:: filenumb, targetsize
    character(len=128) :: infilename
    character(len=64) ::numb1,numb2
    character(len=4) :: stepname
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe,u2spe,u3spe,pspe,rhospe
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1spe_new, u2spe_new, u3spe_new, pspe_new, rhospe_new
    real(8), allocatable, dimension(:,:,:) :: u1_new, u2_new, u3_new, rho_new, p_new
    character(len=128) :: outfilename
    character(len=1) :: modeio
    integer :: i,j,k,n
    integer :: ic,jc,kc
    integer ::ii,jj,kk
    integer(C_INTPTR_T) :: icfftw, kcfftw, jcfftw
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_u3spe, c_pspe, c_rhospe, &
                     c_u1spe_new, c_u2spe_new, c_u3spe_new, c_pspe_new, c_rhospe_new
    !
    if(mpirank == 0) then
        call readkeyboad(numb1) 
        read(numb1,'(i4)') filenumb
    endif
    call bcast(filenumb)
    !
    if(mpirank == 0) then
      call readkeyboad(numb2) 
      read(numb2,'(i4)') targetsize
    endif
    call bcast(targetsize)
    !
    call readinput
    !
    modeio='h'
    ! Initialization
    call fftw_mpi_init()
    if(mpirank==0)  print *, "fftw_mpi initialized"
    if(ka==0) stop 'Please add 2d case by yourself, kevin is too lazy'
    !
    if(mpirank==0) then 
      if(ia .NE. ja .OR. ja .NE. ka .OR. ia .NE. ka) then
        stop 'error!! only for cube fields'
      else
        print *, "ia:",ia,",ja:",ja,",ka:", ka
      endif
    endif
    !
    if(mpirank==0)  print *, "size of target field:", targetsize
    ic=targetsize
    jc=targetsize
    kc=targetsize
    icfftw=targetsize
    jcfftw=targetsize
    kcfftw=targetsize
    !
    call mpisizedis_fftw
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0)  print*, '** refcal done!'
    !
    allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km), prs(0:im,0:jm,0:km))
    !
    !
    if (filenumb .ne. 0) then
      write(stepname,'(i4.4)') filenumb
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
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [imfftw,jmfftw,kmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [imfftw,jmfftw,kmfftw])
    c_u3spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe, u3spe, [imfftw,jmfftw,kmfftw])
    c_pspe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe, pspe, [imfftw,jmfftw,kmfftw])
    c_rhospe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhospe, rhospe, [imfftw,jmfftw,kmfftw])
    !
    c_u1spe_new = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe_new, u1spe_new, [icfftw,jcfftw,kcfftw])
    c_u2spe_new = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe_new, u2spe_new, [icfftw,jcfftw,kcfftw])
    c_u3spe_new = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3spe_new, u3spe_new, [icfftw,jcfftw,kcfftw])
    c_pspe_new = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_pspe_new, pspe_new, [icfftw,jcfftw,kcfftw])
    c_rhospe_new = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_rhospe_new, rhospe_new, [icfftw,jcfftw,kcfftw])
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
      u1spe(i,j,k)=CMPLX(vel(i,j,k,1),0.d0,C_INTPTR_T);
      u2spe(i,j,k)=CMPLX(vel(i,j,k,2),0.d0,C_INTPTR_T);
      u3spe(i,j,k)=CMPLX(vel(i,j,k,3),0.d0,C_INTPTR_T);
      pspe(i,j,k)=CMPLX(prs(i,j,k),0.d0,C_INTPTR_T);
      rhospe(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
      !
    enddo
    enddo
    enddo
    !
    !
    do k=1,kc
    do j=1,jc
    do i=1,ic
      !
      u1spe_new(i,j,k)=CMPLX(0.d0,0.d0,C_INTPTR_T);
      u2spe_new(i,j,k)=CMPLX(0.d0,0.d0,C_INTPTR_T);
      u3spe_new(i,j,k)=CMPLX(0.d0,0.d0,C_INTPTR_T);
      pspe_new(i,j,k)=CMPLX(0.d0,0.d0,C_INTPTR_T);
      rhospe_new(i,j,k)=CMPLX(0.d0,0.d0,C_INTPTR_T);
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
    call fftw_mpi_execute_dft(forward_plan,rhospe,rhospe)
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      u1spe(i,j,k)=u1spe(i,j,k)/(1.d0*ia*ja*ka)
      u2spe(i,j,k)=u2spe(i,j,k)/(1.d0*ia*ja*ka)
      u3spe(i,j,k)=u3spe(i,j,k)/(1.d0*ia*ja*ka)
      pspe(i,j,k)=pspe(i,j,k)/(1.d0*ia*ja*ka)
      rhospe(i,j,k)=rhospe(i,j,k)/(1.d0*ia*ja*ka)
      !
    enddo
    enddo
    enddo
    !  
    !
    !field transform 
    do k=1,kc
      do j=1,jc
        do i=1,ic
            if (i <= (ic/2+1)) then
                ii = i
            else
                ii = i + ia - ic
            end if
            if (j <= (jc/2+1)) then
                jj = j
            else
                jj = j + ja - jc
            end if
            if (k <= (kc/2+1)) then
                kk = k
            else
                kk = k + ka - kc
            end if
            u1spe_new(i,j,k) = u1spe(ii,jj,kk)
            u2spe_new(i,j,k) = u2spe(ii,jj,kk)
            u3spe_new(i,j,k) = u3spe(ii,jj,kk)
            pspe_new(i,j,k) = pspe(ii,jj,kk)
            rhospe_new(i,j,k) = rhospe(ii,jj,kk)
        enddo
      enddo
    enddo
 !!!! Do inverse FFT
    !
    call fftw_mpi_execute_dft(backward_plan,u1spe_new,u1spe_new)
    call fftw_mpi_execute_dft(backward_plan,u2spe_new,u2spe_new)
    call fftw_mpi_execute_dft(backward_plan,u3spe_new,u3spe_new)
    call fftw_mpi_execute_dft(backward_plan,rhospe_new,rhospe_new)
    call fftw_mpi_execute_dft(backward_plan,pspe_new,pspe_new)
!
    allocate(u1_new(0:ic,0:jc,0:kc), u2_new(0:ic,0:jc,0:kc), u3_new(0:ic,0:jc,0:kc), p_new(0:ic,0:jc,0:kc), rho_new(0:ic,0:jc,0:kc))
    !
    do k=1,kc
      do j=1,jc
        do i=1,ic
          u1_new(i,j,k) = real(u1spe_new(i,j,k))
          u2_new(i,j,k) = real(u2spe_new(i,j,k))
          u3_new(i,j,k) = real(u3spe_new(i,j,k))
          p_new(i,j,k) = real(pspe_new(i,j,k))
          rho_new(i,j,k) = real(rhospe_new(i,j,k))
        enddo
      enddo
    enddo
    !
    do j=1,jc
      do i=1,ic
        u1_new(i,j,0) = u1_new(i,j,kc)
        u2_new(i,j,0) = u2_new(i,j,kc)
        u3_new(i,j,0) = u3_new(i,j,kc)
        p_new(i,j,0) = p_new(i,j,kc)
        rho_new(i,j,0) = rho_new(i,j,kc)
      enddo
    enddo
    !
    do k=1,kc
      do i=1,ic
        u1_new(i,0,k) = u1_new(i,jc,k)
        u2_new(i,0,k) = u2_new(i,jc,k)
        u3_new(i,0,k) = u3_new(i,jc,k)
        p_new(i,0,k) = p_new(i,jc,k)
        rho_new(i,0,k) = rho_new(i,jc,k)
      enddo
    enddo
    !
    do k=1,kc
      do j=1,jc
        u1_new(0,j,k) = u1_new(ic,j,k)
        u2_new(0,j,k) = u2_new(ic,j,k)
        u3_new(0,j,k) = u3_new(ic,j,k)
        p_new(0,j,k) = p_new(ic,j,k)
        rho_new(0,j,k) = rho_new(ic,j,k)
      enddo
    enddo
    !
    do i=1,ic
      u1_new(i,0,0) = u1_new(i,jc,kc)
      u2_new(i,0,0) = u2_new(i,jc,kc)
      u3_new(i,0,0) = u3_new(i,jc,kc)
      p_new(i,0,0) = p_new(i,jc,kc)
      rho_new(i,0,0) = rho_new(i,jc,kc)
    enddo
    !
    do j=1,jc
      u1_new(0,j,0) = u1_new(ic,j,kc)
      u2_new(0,j,0) = u2_new(ic,j,kc)
      u3_new(0,j,0) = u3_new(ic,j,kc)
      p_new(0,j,0) = p_new(ic,j,kc)
      rho_new(0,j,0) = rho_new(ic,j,kc)
    enddo
    !
    do k=1,kc
      u1_new(0,0,k) = u1_new(ic,jc,k)
      u2_new(0,0,k) = u2_new(ic,jc,k)
      u3_new(0,0,k) = u3_new(ic,jc,k)
      p_new(0,0,k) = p_new(ic,jc,k)
      rho_new(0,0,k) = rho_new(ic,jc,k)
    enddo
    !
    u1_new(0,0,0) = u1_new(ic,jc,kc)
    u2_new(0,0,0) = u2_new(ic,jc,kc)
    u3_new(0,0,0) = u3_new(ic,jc,kc)
    p_new(0,0,0) = p_new(ic,jc,kc)
    rho_new(0,0,0) = rho_new(ic,jc,kc)
    !
    if (filenumb .ne. 0) then
      outfilename = 'pp/transedfield'//stepname//'.'//modeio//'5'
    else
      outfilename = 'pp/transedfield.'//modeio//'5'
    endif
    !
    ia=ic
    ja=jc
    ka=kc
    !
    call h5io_init(trim(outfilename),mode='write')
    call h5write(varname='ro', var=rho_new(0:ic,0:jc,0:kc),  mode = modeio)
    call h5write(varname='u1', var=u1_new(0:ic,0:jc,0:kc),mode = modeio)
    call h5write(varname='u2', var=u2_new(0:ic,0:jc,0:kc),mode = modeio)
    call h5write(varname='u3', var=u3_new(0:ic,0:jc,0:kc),mode = modeio)
    call h5write(varname='p',  var=p_new(0:ic,0:jc,0:kc),mode = modeio)
    call h5io_end
    if(mpirank==0) then
      call h5srite(varname='time',var=time,filename=trim(outfilename))
      call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
    endif
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_u3spe)
    call fftw_free(c_pspe)
    call fftw_free(c_rhospe)
    call fftw_free(c_u1spe_new)
    call fftw_free(c_u2spe_new)
    call fftw_free(c_u3spe_new)
    call fftw_free(c_pspe_new)
    call fftw_free(c_rhospe_new)
    call mpistop
    !
    deallocate(u1_new, u2_new, u3_new, rho_new, p_new)
    !
  end subroutine specftrans
end module udf_specftrans


