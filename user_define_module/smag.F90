module smag
  !
  use parallel, only : mpirank,mpistop,irk,jrk,krk
  use commarray,only : miut
  use constdef

  implicit none
  !
  type :: smag_coef
   real(8) :: C_s,C_k,C_epsil
   real(8),allocatable,dimension(:,:,:) :: nu_sgs,k_sgs
  end type smag_coef
  !
  type(smag_coef) :: smag_model
  !
  contains
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_smag
    !
    use commvar,only : im,jm,km,hm
    !
    integer :: lallo
    !
    allocate(miut(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating miut'
    !
    allocate( smag_model%k_sgs(0:im,0:jm,0:km),    & 
    smag_model%nu_sgs(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating nu_sgs k_sgs'
    !
    smag_model%C_s = 0.182d0
    smag_model%C_k =0.094d0
    smag_model%C_epsil=1.048d0
    !
  end subroutine init_smag
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine src_smag
    !
    use commvar,  only : im,jm,km
    use commarray,only :rho,dvel,cell
    !
    ! local data
    integer :: i,j,k
    real(8) :: a,b,c
    real(8) :: s11,s12,s13,s22,s23,s33,tr,D,ro
    real(8) :: delta
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ro=rho(i,j,k)
      delta=(cell(1,1,1)%vol)**(1.d0/3.d0)
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      tr=s11+s22+s33
      D=s11*s11+s22*s22+s33*s33+2.d0*(s12*s12+s13*s13*s23*s23)
      !
      smag_model%nu_sgs(i,j,k)=((smag_model%C_s*delta)**2)*sqrt(2*D)
      smag_model%k_sgs(i,j,k)=ro*smag_model%nu_sgs(i,j,k)*tr
      miut(i,j,k)=ro*smag_model%nu_sgs(i,j,k)
      !
      !D=D-(tr*tr)/3.d0 !dev(D):D reference:openFOAM smagorinsky.C
      !
      ! a=smag_model%C_epsil/delta
      ! b=num2d3*tr
      ! c=2.d0*smag_model%C_k*delta*D
      !
      ! smag_model%k_sgs(i,j,k)=((-b+sqrt(b*b+4.d0*a*c))/(2.d0*a))**2.d0
      ! smag_model%nu_sgs(i,j,k)=smag_model%C_k*delta*sqrt(smag_model%k_sgs(i,j,k))
      ! miut(i,j,k)=ro*smag_model%nu_sgs(i,j,k)
    enddo
    enddo
    enddo
      !
    end subroutine src_smag
    !
end module smag

