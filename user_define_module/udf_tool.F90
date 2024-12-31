module udf_tool
    !
    !
    implicit none
    !
    interface ProjectP2
        module procedure ProjectP2_2D
        module procedure ProjectP2_3D
    end interface
    !
    interface ProjectP3
        module procedure ProjectP3_2D
        module procedure ProjectP3_3D
    end interface
    !
    interface ProjectPi2
        module procedure ProjectPi2_2D
        module procedure ProjectPi2_3D
    end interface
    !
    interface ProjectPi3
        module procedure ProjectPi3_2D
        module procedure ProjectPi3_3D
    end interface
    !
    interface GenerateWave
        module procedure GenerateWave_2D
        module procedure GenerateWave_3D
    end interface
    !
    contains
    !
    subroutine GenerateWave_2D(im,jm,ia,ja,j0,k1,k2)
        implicit none
        integer, intent(in) :: im,jm,ia,ja,j0
        real(8), dimension(:,:),intent(out) :: k1,k2
        integer :: i,j
        !
        do j=1,jm
        do i=1,im
        !
        if(im .ne. ia)then
            stop "GenerateWave Error! im /= ia"
        endif
        !
        if(i <= (ia/2+1)) then
            k1(i,j) = real(i-1,8)
        else if(i<=(ia)) then
            k1(i,j) = real(i-ia-1,8)
        else
            stop "GenerateWave Error, no wave number possible, i must smaller than ia-1 !"
        end if
        !
        if((j+j0) <= (ja/2+1)) then
            k2(i,j) = real(j+j0-1,8)
        else if((j+j0)<=(ja)) then
            k2(i,j) = real(j+j0-ja-1,8)
        else
            stop "GenerateWave Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
        end do
        end do
    end subroutine GenerateWave_2D
    !
    subroutine GenerateWave_3D(im,jm,km,ia,ja,ka,k0,k1,k2,k3)
        implicit none
        integer, intent(in) :: im,jm,km,ia,ja,ka,k0
        real(8), dimension(:,:,:),intent(out) :: k1,k2,k3
        integer :: i,j,k
        !
        do k=1,km
        do j=1,jm
        do i=1,im
        !
        if(im .ne. ia)then
            stop "GenerateWave Error! im /= ia"
        endif
        !
        if(i <= (ia/2+1)) then
            k1(i,j,k) = real(i-1,8)
        else if(i<=ia) then
            k1(i,j,k) = real(i-ia-1,8)
        else
            stop "GenerateWave Error, no wave number possible, i must smaller than ia-1 !"
        end if
        !
        if(j <= (ja/2+1)) then
            k2(i,j,k) = real(j-1,8)
        else if(j<=ja) then
            k2(i,j,k) = real(j-ja-1,8)
        else
            stop "GenerateWave Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
            k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=ka) then
            k3(i,j,k) = real(k+k0-ka-1,8)
        else
            stop "GenerateWave Error, no wave number possible, (k+k0) must smaller than ka-1 !"
        end if
        !
        enddo
        enddo
        enddo
    end subroutine GenerateWave_3D
    !
    real(8) function wav(i,im)
    ! This function gives the wave number of index i
    ! with the maximum im
        implicit none
        integer, intent(in) :: i, im
        !
        if(i <= (im/2+1)) then
        wav = real(i-1,8)
        else if(i<=im) then
        wav = real(i-im-1,8)
        else
        print *,"Error, no wave number possible, i must smaller than im !"
        end if
    end function wav
    !
    integer function invwav(i,im)
    ! This function gives the wave number of index i
    ! with the maximum im
        implicit none
        integer, intent(in) :: i,im
        !
        if(i < 0) then
        invwav = i + im + 1;
        else
        invwav = i + 1;
        end if
    end function invwav
    !
    integer function kint(k,dk,method,lambda)
    ! This function gives the nearby k
    !!
        implicit none
        real(8), intent(in) :: k,dk
        integer, intent(in) :: method
        real(8), intent(in), optional :: lambda
        !
        if(method == 1)then
        if( (.not. present(lambda)) .or. (lambda <= 0))then
            stop "Error! kint method 1 with no lambda or lambda is negative!"
        else
            if(k<(dk/2))then
            kint = 0
            else
            kint = floor(log(k/dk)/log(lambda))+1
            endif
        endif
        else
        kint = nint(k/dk)
        endif
    end function kint
    !
    real(8) function ProjectP2_2D(i,j,kx,ky)
        !
        !!
        implicit none
        integer, intent(in) :: i,j
        real(8), intent(in) ::kx,ky
        real(8) :: k
        !
        k = dsqrt(kx**2+ky**2+1.d-15)
        !
        if((i .eq. 1) .and. (j .eq. 1)) then
        ProjectP2_2D = 1.d0 - kx*kx/k/k
        else if((i .eq. 1) .and. (j .eq. 2)) then
        ProjectP2_2D = - kx*ky/k/k
        else if((i .eq. 2) .and. (j .eq. 1)) then
        ProjectP2_2D = - kx*ky/k/k
        else if((i .eq. 2) .and. (j .eq. 2)) then
        ProjectP2_2D = 1.d0 - ky*ky/k/k
        else
        stop "ProjectP2_2D error: i,j"
        end if
        !
    end function ProjectP2_2D
    !
    real(8) function ProjectP2_3D(i,j,kx,ky,kz)
    !
    !!
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: kx,ky,kz
    real(8) :: k
    !
    k = dsqrt(kx**2+ky**2+kz**2+1.d-15)
    !
    if((i .eq. 1) .and. (j .eq. 1)) then
        ProjectP2_3D = 1.d0 - kx*kx/k/k
    else if((i .eq. 1) .and. (j .eq. 2)) then
        ProjectP2_3D = - kx*ky/k/k
    else if((i .eq. 1) .and. (j .eq. 3)) then
        ProjectP2_3D = - kx*kz/k/k
    else if((i .eq. 2) .and. (j .eq. 1)) then
        ProjectP2_3D = - ky*kx/k/k
    else if((i .eq. 2) .and. (j .eq. 2)) then
        ProjectP2_3D = 1.d0 - ky*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 3)) then
        ProjectP2_3D = - ky*kz/k/k
    else if((i .eq. 3) .and. (j .eq. 1)) then
        ProjectP2_3D = - kz*kx/k/k
    else if((i .eq. 3) .and. (j .eq. 2)) then
        ProjectP2_3D = - kz*ky/k/k
    else if((i .eq. 3) .and. (j .eq. 3)) then
        ProjectP2_3D = 1.d0 - kz*kz/k/k
    else
        stop "ProjectP2_3D error: i,j"
    end if
    !
    end function ProjectP2_3D
    !
    real(8) function ProjectP3_2D(i,j,m,kx,ky)
        !
        !!
        implicit none
        integer, intent(in) :: i,j,m
        real(8), intent(in) :: kx,ky
        !
        if((j .eq. 1) .and. (m .eq. 1))then
        ProjectP3_2D = 0.5d0 * kx * ProjectP2(i,j,kx,ky)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky)
        else if((j .eq. 2) .and. (m .eq. 1))then
        ProjectP3_2D = 0.5d0 * kx * ProjectP2(i,j,kx,ky)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky)
        else if((j .eq. 1) .and. (m .eq. 2))then
        ProjectP3_2D = 0.5d0 * ky * ProjectP2(i,j,kx,ky)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky)
        else if((j .eq. 2) .and. (m .eq. 2))then
        ProjectP3_2D = 0.5d0 * ky * ProjectP2(i,j,kx,ky)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky)
        else
        stop "ProjectP3_2D error: j,m"
        endif
    end function ProjectP3_2D
    !
    real(8) function ProjectP3_3D(i,j,m,kx,ky,kz)
    !
    !!
    implicit none
    integer, intent(in) :: i,j,m
    real(8), intent(in) :: kx,ky,kz
    !
    if((j .eq. 1) .and. (m .eq. 1))then
        ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 2) .and. (m .eq. 1))then
        ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 3) .and. (m .eq. 1))then
        ProjectP3_3D = 0.5d0 * kx * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 1) .and. (m .eq. 2))then
        ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 2) .and. (m .eq. 2))then
        ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 3) .and. (m .eq. 2))then
        ProjectP3_3D = 0.5d0 * ky * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 1) .and. (m .eq. 3))then
        ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kx * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 2) .and. (m .eq. 3))then
        ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * ky * ProjectP2(i,m,kx,ky,kz)
    else if((j .eq. 3) .and. (m .eq. 3))then
        ProjectP3_3D = 0.5d0 * kz * ProjectP2(i,j,kx,ky,kz)+ 0.5d0 * kz * ProjectP2(i,m,kx,ky,kz)
    else
        stop "ProjectP3_3D error: j,m"
    endif
    end function ProjectP3_3D
    !
    real(8) function ProjectPi2_2D(i,j,kx,ky)
        !
        !!
        implicit none
        integer, intent(in) :: i,j
        real(8), intent(in) :: kx,ky
        real(8) :: k
        !
        k = dsqrt(kx**2+ky**2+1.d-15)
        !
        if((i .eq. 1) .and. (j .eq. 1)) then
        ProjectPi2_2D = kx*kx/k/k
        else if((i .eq. 1) .and. (j .eq. 2)) then
        ProjectPi2_2D = kx*ky/k/k
        else if((i .eq. 2) .and. (j .eq. 1)) then
        ProjectPi2_2D = kx*ky/k/k
        else if((i .eq. 2) .and. (j .eq. 2)) then
        ProjectPi2_2D = ky*ky/k/k
        else
        stop "ProjectPi2_2D error: i,j"
        end if
        !
    end function ProjectPi2_2D
    !
    real(8) function ProjectPi2_3D(i,j,kx,ky,kz)
    !
    !!
    implicit none
    integer, intent(in) :: i,j
    real(8), intent(in) :: kx,ky,kz
    real(8) :: k
    !
    k = dsqrt(kx**2+ky**2+kz**2+1.d-15)
    !
    if((i .eq. 1) .and. (j .eq. 1)) then
        ProjectPi2_3D = kx*kx/k/k
    else if((i .eq. 1) .and. (j .eq. 2)) then
        ProjectPi2_3D = kx*ky/k/k
    else if((i .eq. 1) .and. (j .eq. 3)) then
        ProjectPi2_3D = kx*kz/k/k
    else if((i .eq. 2) .and. (j .eq. 1)) then
        ProjectPi2_3D = ky*kx/k/k
    else if((i .eq. 2) .and. (j .eq. 2)) then
        ProjectPi2_3D = ky*ky/k/k
    else if((i .eq. 2) .and. (j .eq. 3)) then
        ProjectPi2_3D = ky*kz/k/k
    else if((i .eq. 3) .and. (j .eq. 1)) then
        ProjectPi2_3D = kz*kx/k/k
    else if((i .eq. 3) .and. (j .eq. 2)) then
        ProjectPi2_3D = kz*ky/k/k
    else if((i .eq. 3) .and. (j .eq. 3)) then
        ProjectPi2_3D = kz*kz/k/k
    else
        stop "ProjectPi2_3D error: i,j"
    end if
    !
    end function ProjectPi2_3D
    !
    real(8) function ProjectPi3_2D(i,j,m,kx,ky)
        !
        !!
        implicit none
        integer, intent(in) :: i,j,m
        real(8), intent(in) :: kx,ky
        !
        if(m .eq. 1)then
        ProjectPi3_2D = kx * ProjectP2(i,j,kx,ky)
        else if(m .eq. 2)then
        ProjectPi3_2D = ky * ProjectP2(i,j,kx,ky)
        else
        stop "ProjectPi3_2D error: m"
        endif
    end function ProjectPi3_2D
    !
    real(8) function ProjectPi3_3D(i,j,m,kx,ky,kz)
        !
        !!
        implicit none
        integer, intent(in) :: i,j,m
        real(8), intent(in) :: kx,ky,kz
        !
        if(m .eq. 1)then
        ProjectPi3_3D = kx * ProjectP2(i,j,kx,ky,kz)
        else if(m .eq. 2)then
        ProjectPi3_3D = ky * ProjectP2(i,j,kx,ky,kz)
        else if(m .eq. 3)then
        ProjectPi3_3D = kz * ProjectP2(i,j,kx,ky,kz)
        else
        stop "ProjectPi3_3D error: m"
        endif
    end function ProjectPi3_3D
end module udf_tool