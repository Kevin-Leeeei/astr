!+---------------------------------------------------------------------+
!| This module contains subroutines used to help the development of FFT with FFTW      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Mar-2024  | Created by C.S. Luo @ Beihang University             |
!| 02-Dec-2024  | Add for forcing
!+---------------------------------------------------------------------+
module fftwlink
    !
    use, intrinsic :: iso_c_binding
    use commvar,   only : im,jm,km,hm,ia,ja,ka, ndims
    use parallel,  only : isize,jsize,ksize,irkm,jrkm,krkm,mpirank,bcast, mpirankmax
    !
    implicit none
    !
    !
    interface fftw_grid_fence
        module procedure fftw_grid_fence_2D
        module procedure fftw_grid_fence_3D
    end interface
    !
    interface fftw_fence_grid
        module procedure fftw_fence_grid_2D
        module procedure fftw_fence_grid_3D
    end interface
    !
    integer(C_INTPTR_T) :: alloc_local,iafftw,jafftw,kafftw,imfftw,jmfftw,kmfftw
    integer :: nproc, myid, ierr
    integer :: i0f,j0f,k0f
    ! Reserved for forcing part
    integer :: rooti0, imf, jmf, kmf
    integer, allocatable, dimension(:) :: counts_grid,disp_grid,counts_fence,disp_fence
    !
    contains
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is used to assign the distributions size on the ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writen by Chensheng Luo, 2024-03-06.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mpisizedis_fftw
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        ! 
        iafftw = ia
        jafftw = ja
        kafftw = ka
        !
        !
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        isize = 1
        jsize = 1
        ksize = 1
        !
        if(ndims==3) then
            ! 
            alloc_local = fftw_mpi_local_size_3d(kafftw,jafftw,iafftw,MPI_COMM_WORLD, kmfftw,k0fftw)
            jmfftw = jafftw
            imfftw = iafftw
            krk = myid
            ksize = nproc
            !
        elseif(ndims==2) then 
            alloc_local = fftw_mpi_local_size_2d(jafftw,iafftw,MPI_COMM_WORLD, jmfftw, j0fftw)
            imfftw = iafftw
            kmfftw = 0
            jrk = myid
            jsize = nproc
            !
        else
            stop '1D not implemented!'
            !alloc_local = fftw_mpi_local_size(myid,iafftw,MPI_COMM_WORLD, imfftw, i0fftw)
            jmfftw = 0
            kmfftw = 0
            irk = myid
            isize = nproc
            !
        endif
        !
        if(mpirank == 0)then
            write(*,'(3(A,I0))')'  ** mpi size= ',isize,' x ',jsize,' x ',ksize
        endif
        !
        imf = imfftw
        jmf = jmfftw
        kmf = kmfftw
        im = imf
        jm = jmf
        km = kmf
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        allocate(ids(0:(nproc-1)),irks(0:(nproc-1)),jrks(0:(nproc-1)),krks(0:(nproc-1)),&
                ims(0:(nproc-1)),jms(0:(nproc-1)),kms(0:(nproc-1)),i0s(0:(nproc-1)),    &
                j0s(0:(nproc-1)),k0s(0:(nproc-1)))
        !
        do i=0,(nproc-1)
            if(myid == i) then
                ids(i) = myid
                irks(i) = irk
                jrks(i) = jrk
                krks(i) = krk
                ims(i) = im 
                jms(i) = jm
                kms(i) = km 
                i0s(i) = i0fftw
                j0s(i) = j0fftw
                k0s(i) = k0fftw
            endif
            call mpi_bcast(ids(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(irks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jrks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(krks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(ims(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(kms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(i0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(j0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(k0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
        enddo
        !
        if(mpirank == 0)then
            open(fh,file='datin/parallel.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")isize,jsize,ksize
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,(nproc-1)
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print*,' << parallel.info ... done !'
        endif
        !
        irkm=isize-1
        jrkm=jsize-1
        krkm=ksize-1
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(mpirank == 0) print*,' ** parallel processing ... done.'
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine mpisizedis_fftw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The end of the subroutine mpisizedis_fftw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    subroutine mpisizedis_half_fftw
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        ! 
        iafftw = ia
        jafftw = ja
        kafftw = ka
        !
        !
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        isize = 1
        jsize = 1
        ksize = 1
        !
        if(ndims==3) then
            ! 
            alloc_local = fftw_mpi_local_size_3d(kafftw,jafftw,iafftw/2+1,MPI_COMM_WORLD, kmfftw, k0fftw)
            k0fftw = k0fftw
            jmfftw = jafftw
            imfftw = iafftw/2+1
            krk = myid
            ksize = nproc
            !
        elseif(ndims==2) then 
            alloc_local = fftw_mpi_local_size_2d(jafftw,iafftw/2+1,MPI_COMM_WORLD, jmfftw, j0fftw)
            j0fftw = j0fftw
            imfftw = iafftw/2+1
            kmfftw = 0
            jrk = myid
            jsize = nproc
            !
        else
            stop '1D not implemented!'
            !alloc_local = fftw_mpi_local_size(myid,iafftw,MPI_COMM_WORLD, imfftw, i0fftw)
            jmfftw = 0
            kmfftw = 0
            irk = myid
            isize = nproc
            !
        endif
        !
        if(mpirank == 0)then
            write(*,'(3(A,I0))')'  ** mpi size= ',isize,' x ',jsize,' x ',ksize
        endif
        !
        imf = imfftw
        jmf = jmfftw
        kmf = kmfftw
        im = imf
        jm = jmf
        km = kmf
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        allocate(ids(0:(nproc-1)),irks(0:(nproc-1)),jrks(0:(nproc-1)),krks(0:(nproc-1)),&
                ims(0:(nproc-1)),jms(0:(nproc-1)),kms(0:(nproc-1)),i0s(0:(nproc-1)),    &
                j0s(0:(nproc-1)),k0s(0:(nproc-1)))
        !
        do i=0,(nproc-1)
            if(myid == i) then
                ids(i) = myid
                irks(i) = irk
                jrks(i) = jrk
                krks(i) = krk
                ims(i) = im
                jms(i) = jm
                kms(i) = km
                i0s(i) = i0fftw
                j0s(i) = j0fftw
                k0s(i) = k0fftw
            endif
            call mpi_bcast(ids(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(irks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jrks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(krks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(ims(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(kms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(i0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(j0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(k0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
        enddo
        !
        if(mpirank == 0)then
            open(fh,file='datin/parallel.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")isize,jsize,ksize
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,(nproc-1)
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print*,' << parallel.info ... done !'
        endif
        !
        irkm=isize-1
        jrkm=jsize-1
        krkm=ksize-1
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(mpirank == 0) print*,' ** parallel processing ... done.'
        !
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine mpisizedis_half_fftw
    ! !
    subroutine fftwprepare_forcing
        !
        if(ndims == 2) then
            call fftwprepare_forcing2D
        else
            stop "Not implemented error! fftwprepare_forcing3D"
        endif
    end subroutine fftwprepare_forcing
    !
    subroutine fftwprepare_forcing2D
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i,ierr
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        rooti0 = isize * (mpirank/isize)
        jmf = jm/isize
        !
        if((jmf*ia) .ne. (im*jm))then
            print *,'jmf = ', jmf, 'ia=', ia, 'im=', im, 'jm=', jm
            stop 'fftwprepare_forcing2D: Parallel ranking not match!'
        endif
        !
        allocate(counts_grid(0:mpirankmax),disp_grid(0:mpirankmax))
        !
        counts_grid = 0
        disp_grid = 0
        counts_grid(rooti0) = im*jm
        do i = (rooti0+1), (rooti0+isize-1)
            counts_grid(i) = im*jm
            disp_grid(i) =disp_grid(i-1) + counts_grid(i-1)
        enddo
        !
        allocate(counts_fence(0:mpirankmax),disp_fence(0:mpirankmax))
        counts_fence = 0
        disp_fence = 0
        counts_fence(rooti0) = ia*jmf
        do i = (rooti0+1), (rooti0+isize-1)
            counts_fence(i) = ia*jmf
            disp_fence(i) =disp_fence(i-1) + counts_fence(i-1)
        enddo
        !
        iafftw = ia
        jafftw = ja
        kafftw = ka
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        !
        alloc_local = fftw_mpi_local_size_2d(jafftw, iafftw,MPI_COMM_WORLD, jmfftw, j0fftw)
        !
        imfftw = iafftw
        kmfftw = 0
        !
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        imf = ia
        kmf = 0
        if(jmfftw .ne. jmf)then
            print *, 'mpirank = ', mpirank, 'jmfftw=',jmfftw, 'jmf=', jmf
            stop "fftwprepare_forcing2D: Rank not matching!"
        endif
        !
        allocate(ids(0:mpirankmax),irks(0:mpirankmax),jrks(0:mpirankmax),krks(0:mpirankmax),&
                ims(0:mpirankmax),jms(0:mpirankmax),kms(0:mpirankmax),i0s(0:mpirankmax),    &
                j0s(0:mpirankmax),k0s(0:mpirankmax))
        !
        do i=0,mpirankmax
            if(mpirank == i) then
                ids(i) = mpirank
                irks(i) = 0
                jrks(i) = mpirank
                krks(i) = 0
                ims(i) = imf
                jms(i) = jmf
                kms(i) = kmf 
                i0s(i) = i0f
                j0s(i) = j0f
                k0s(i) = k0f
            endif
            call mpi_bcast(ids(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(irks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jrks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(krks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(ims(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(kms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(i0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(j0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(k0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
        enddo
        !
        if(mpirank == 0)then
            open(fh,file='datin/parallelfftw.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")1,mpirank,1
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,mpirankmax
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print *,' << parallelfftw.info ... done !'
        endif
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine fftwprepare_forcing2D
    !
    function fftw_grid_fence_2D(gridarray_t) result(fencearray)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:)  ::  gridarray_t 
        real(8), allocatable, dimension(:,:) ::  fencearray
        real(8), allocatable, dimension(:,:) :: global_t, global
        !
        allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        !
        if (.not. allocated(fencearray)) allocate(fencearray(1:ia,1:jmf))
        !
        call mpi_gatherv(gridarray_t,im*jm, MPI_REAL8, &
        global_t,counts_grid,disp_grid,MPI_REAL8,rooti0, &
        MPI_COMM_WORLD, ierr)
        !
        global = transpose(global_t)
        !
        call mpi_scatterv(global,counts_fence,disp_fence,MPI_REAL8,&
                        fencearray,ia*jmf,MPI_REAL8,rooti0,&
                            MPI_COMM_WORLD,ierr)
        !
        deallocate(global_t,global)
        !
    end function fftw_grid_fence_2D
    !
    function fftw_fence_grid_2D(fencearray) result(gridarray_t)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:)  :: fencearray
        real(8), allocatable, dimension(:,:) :: gridarray_t
        real(8), allocatable, dimension(:,:) :: global_t, global
        !
        allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        if (.not. allocated(gridarray_t)) allocate(gridarray_t(1:jm,1:im))
        !
        call mpi_gatherv(fencearray,ia*jmf, MPI_REAL8, &
        global,counts_fence,disp_fence,MPI_REAL8,rooti0,  &
        MPI_COMM_WORLD, ierr)
        !
        global_t = transpose(global)
        !
        call mpi_scatterv(global_t,counts_grid,disp_grid,MPI_REAL8,&
                            gridarray_t,im*jm,MPI_REAL8,rooti0,  &
                            MPI_COMM_WORLD,ierr)
        !
        deallocate(global_t,global)
        !
    end function fftw_fence_grid_2D
    !
    function fftw_grid_fence_3D(gridarray_t) result(fencearray)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:,:)  ::  gridarray_t 
        real(8), allocatable, dimension(:,:,:) ::  fencearray
        real(8), allocatable, dimension(:,:,:) :: global_t, global
        !
        ! allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        ! !
        ! if (.not. allocated(fencearray)) allocate(fencearray(1:ia,1:jmf))
        ! !
        ! call mpi_gatherv(gridarray_t,im*jm, MPI_REAL8, &
        ! global_t,counts_grid,disp_grid,MPI_REAL8,rooti0, &
        ! MPI_COMM_WORLD, ierr)
        ! !
        ! global = transpose(global_t)
        ! !
        ! call mpi_scatterv(global,counts_fence,disp_fence,MPI_REAL8,&
        !                 fencearray,ia*jmf,MPI_REAL8,rooti0,&
        !                     MPI_COMM_WORLD,ierr)
        ! !
        ! deallocate(global_t,global)
        !
    end function fftw_grid_fence_3D
    !
    function fftw_fence_grid_3D(fencearray) result(gridarray_t)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:,:)  :: fencearray
        real(8), allocatable, dimension(:,:,:) :: gridarray_t
        real(8), allocatable, dimension(:,:,:) :: global_t, global
        ! !
        ! allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        ! if (.not. allocated(gridarray_t)) allocate(gridarray_t(1:jm,1:im))
        ! !
        ! call mpi_gatherv(fencearray,ia*jmf, MPI_REAL8, &
        ! global,counts_fence,disp_fence,MPI_REAL8,rooti0,  &
        ! MPI_COMM_WORLD, ierr)
        ! !
        ! global_t = transpose(global)
        ! !
        ! call mpi_scatterv(global,counts_grid,disp_grid,MPI_REAL8,&
        !                     gridarray_t,im*jm,MPI_REAL8,rooti0,  &
        !                     MPI_COMM_WORLD,ierr)
        ! !
        ! deallocate(global_t,global)
        !
    end function fftw_fence_grid_3D
!     subroutine fftwalloc_find_position(n, dir, proc, pos_in_proc)
!         implicit none
!         !
!         include 'mpif.h' 
!         !
!         integer, intent(in) :: n
!         character(len=1), intent(in) :: dir
!         integer, intent(out):: proc, pos_in_proc
!         integer :: i
!         !
!         proc = -1
!         call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
!         !
!         if(dir == 'i') then
!             do i=0,(nproc-2)
!                 if(i0s(i) <= n .and. i0s(i+1)> n)then
!                     proc = i
!                     pos_in_proc = n - i0s(i)
!                     exit
!                 endif
!             enddo
!             if(proc == -1) then
!                 if(i0s(nproc-1) <= n .and. ia > n)then
!                     proc = nproc - 1
!                     pos_in_proc = n - i0s(i)
!                 else
!                     print *, 'Error at n=',n, 'dir=', dir
!                     stop "error in finding allocation"
!                 endif
!             endif
!         elseif(dir == 'j') then
!             do i=0,(nproc-2)
!                 if(j0s(i) <= n .and. j0s(i+1)> n)then
!                     proc = i
!                     pos_in_proc = n - j0s(i)
!                     exit
!                 endif
!             enddo
!             if(proc == -1) then
!                 if(j0s(nproc-1) <= n .and. ja > n)then
!                     proc = nproc - 1
!                     pos_in_proc = n - j0s(i)
!                 else
!                     print *, 'Error at n=',n, 'dir=', dir
!                     stop "error in finding allocation"
!                 endif
!             endif
!         elseif(dir == 'k' ) then
!             do i=0,(nproc-2)
!                 if(k0s(i) <= n .and. k0s(i+1)> n)then
!                     proc = i
!                     pos_in_proc = n - k0s(i)
!                     exit
!                 endif
!             enddo
!             if(proc == -1) then
!                 if(k0s(nproc-1) <= n .and. ka > n)then
!                     proc = nproc - 1
!                     pos_in_proc = n - k0s(i)
!                 else
!                     print *, 'Error at n=',n, 'dir=', dir
!                     stop "error in finding allocation"
!                 endif
!             endif
!         else
!             stop 'Direction not correct'
!         endif
!         !
!   end subroutine fftwalloc_find_position
  !
end module fftwlink