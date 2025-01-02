module udf_pp_Skewness
    !
    !
    use constdef
    use stlaio,  only: get_unit
    !
    implicit none
    !
    contains 
    !
    subroutine skewnesscal
        !
        !
        use cmdefne
        use singleton
        use readwrite, only : readinput
        use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka
        use commarray, only :x,vel,dvel
        use hdf5io
        use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,psum,&
                                mpirank,bcast,mpisize,lio
        use comsolver, only : solvrinit,grad
        use solver,    only : refcal
        use geom,      only : geomcal
        use utility,  only : listinit,listwrite
        use gridgeneration
        !
        !local data
        character(len=64) :: inputrange,numb1,numb2
        integer :: start_numb,end_numb,filenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=1) :: modeio
        integer :: i,j,k
        !
        real(8)::du11,du12,du21,du22,du33,du23,du13,du32,du31,du3,du2
        real(8)::skewness,rsamples
        logical,save :: linit=.true.
        integer,save :: hand_fs
        !
        if(mpirank == 0) then
            call readkeyboad(numb1) 
            read(numb1,'(i4)') start_numb
        endif
        call bcast(start_numb)
		!
        if(mpirank == 0) then
            call readkeyboad(numb2) 
            read(numb2,'(i4)') end_numb
        endif
        call bcast(end_numb)
        print*,' ** input file range: ',start_numb,' to ',end_numb
        !
        !
        call readinput
        !
        call mpisizedis
        if(mpirank==0)then
            print*, '** mpisizedis done!'
        endif
        !
        call parapp
        if(mpirank==0)then
            print*, '** parapp done!'
        endif
        !
        call parallelini
        if(mpirank==0)then
            print*, '** parallelini done!'
        endif
        !
        call refcal
        if(mpirank==0)then
            print*, '** refcal done!'
        endif
        !
        modeio='h'
        !
        if(mpirank==0)then
            if(ka==0)then
                print *,"2D, ia:",ia,",ja:",ja
            else
                print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
            endif
        endif
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
        allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
        !
        if(ka==0)then
            call gridsquare(2.d0*pi,2.d0*pi)
        else
            call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
        endif
        !
        if(mpirank==0)then
            print*, '** gird generation done!'
        endif
        !
        call geomcal
        !
        if(mpirank==0)then
            print*, '** geomcal done!'
        endif
        !
        call solvrinit
        !
        if(linit) then
            !
            if(lio) then
                call listinit(filename='log/skewness.dat',handle=hand_fs, &
                 firstline='nstep time skewness')
            endif
            !
            linit=.false.
            !
        endif
        !
        do filenumb=start_numb,end_numb
            if (filenumb .ne. 0) then
                write(stepname,'(i4.4)')filenumb
                infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
            else
                infilename='outdat/flowfield.'//modeio//'5'
            endif
            !
            call h5io_init(filename=infilename,mode='read')
            !
            call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
            call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
            call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
            call h5read(varname='time',var=time)
            call h5read(varname='nstep',var=nstep)
            !
            call h5io_end
            !
            if(mpirank==0)then
                print *, "Swap velocity"
            endif
            !
            call dataswap(vel)
            !
            if(mpirank==0)then
                print *, "Calculate skewness"
            endif
            !
            dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
            dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
            if(ka .ne. 0)then
                dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
            endif
            !
            if(ka==0) then
                rsamples=dble(ia*ja)
                du3=0.d0
                du2=0.d0
                k=0
                do j=1,jm
                do i=1,im
                    du11=dvel(i,j,k,1,1)
                    du12=dvel(i,j,k,1,2)
                    du21=dvel(i,j,k,2,1)
                    du22=dvel(i,j,k,2,2)
                    du3=du3+(du11*du11*du11+du22*du22*du22)
                    du2=du2+(du11*du11+du22*du22)
                enddo
                enddo
                skewness=psum(du3)/(2.d0*rsamples)/sqrt((psum(du2)/(2.d0*rsamples))**3)
                if(lio) then 
                    call listwrite(hand_fs,skewness)
                endif
            else
                rsamples=dble(ia*ja*ka)
                du3=0.d0
                du2=0.d0
                do k=1,km
                do j=1,jm
                do i=1,im
                  du11=dvel(i,j,k,1,1)
                  du12=dvel(i,j,k,1,2)
                  du13=dvel(i,j,k,1,3)
                  du21=dvel(i,j,k,2,1)
                  du22=dvel(i,j,k,2,2)
                  du23=dvel(i,j,k,2,3)
                  du31=dvel(i,j,k,3,1)
                  du32=dvel(i,j,k,3,2)
                  du33=dvel(i,j,k,3,3)
                  du3=du3+(du11*du11*du11+du22*du22*du22+du33*du33*du33)
                  du2=du2+(du11*du11+du22*du22+du33*du33)
                enddo
                enddo
                enddo
                skewness=psum(du3)/(3.d0*rsamples)/sqrt((psum(du2)/(3.d0*rsamples))**3)
                if(lio) then 
                    call listwrite(hand_fs,skewness)
                endif
            endif
              !
        end do
        deallocate(x,vel,dvel)
        call mpistop
    end subroutine skewnesscal

end module udf_pp_Skewness

