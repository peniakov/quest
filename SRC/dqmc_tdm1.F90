module DQMC_TDM1
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_STRUCT
  use DQMC_GTAU

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gtau (up and down).
  !
  ! Two measurements are taken in this module: Green's function (G) and Chi
  ! (chi) function.  Both are stored in a 3 dimensional array. 
  ! The first dimension is for space, the second is for time, and the 
  ! last dimension is for bins.
  !
  type TDM1
     integer  :: n
     integer  :: nClass
     integer  :: L
     integer  :: nbin
     integer  :: avg
     integer  :: err
     integer  :: idx

     integer  :: tmp
     integer  :: cnt

     logical  :: compute=.false.

     real(wp) :: dtau
     real(wp), pointer :: sgn(:)
     real(wp), pointer :: gnl(:,:,:)    ! Green's function
     real(wp), pointer :: chinl(:,:,:)  ! S+S- susceptibility

     integer,  pointer :: D(:,:)
     integer,  pointer :: F(:)

  end type TDM1
  
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDM1_Init(n, L, dtau, T1, nBin, S)
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM1. 
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! time dependent measurement
    integer, intent(in)       :: n       ! No. of sites
    integer, intent(in)       :: L       ! No of time slice
    integer, intent(in)       :: nBin    ! No of Bins
    real(wp), intent(in)      :: dtau
    type(Struct), intent(in)  :: S

    ! ... local variables ...
    integer     :: nClass

    ! ... Executable ...

    T1%n      =  n
    nClass    =  S%nClass
    T1%nClass =  nClass
    T1%L      =  L
    T1%D      => S%D
    T1%F      => S%F
    T1%dtau   =  dtau
    T1%nbin   =  nBin

    T1%tmp    =  nBin + 1
    T1%avg    =  nBin + 1
    T1%err    =  nBin + 2
    T1%idx    =  1

    T1%compute = .true.
    
    ! Allocate storages
    allocate(T1%sgn(nBin+2))
    allocate(T1%gnl(nClass, 0:L-1, nBin+2))
    allocate(T1%chinl(nClass, 0:L-1, nBin+2))

    ! initialize values
    T1%sgn   = ZERO
    T1%gnl   = ZERO
    T1%chinl = ZERO
    
  end subroutine DQMC_TDM1_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Free(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! TDM to be freed

    ! ... Executable ...

    if (.not.T1%compute) return
    deallocate(T1%gnl, T1%chinl)

  end subroutine DQMC_TDM1_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Meas(T1, tau)
    !
    ! Purpose
    ! =======
    !    This subroutine fills the bin. It assumes that tau%A_up and,
    !    when necessary, tau%A_dn have been filled.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)   :: t1
    type(Gtau), intent(inout)   :: tau
    
    ! ... Local var ...
    integer  :: i, k, m, L, cnt, dt, i0, it, j0, jt, dtau
    real(wp) :: sgn, factor
    real(wp),pointer :: up0t(:,:), upt0(:,:), dn0t(:,:), dnt0(:,:)

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(tau%north,2)
    m   = (tau%north-k) / 2

    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t

    blocks: do i0 = 1, tau%nb
       do dtau = 0, tau%nb-1
          it = mod(i0+dtau-1,tau%nb) + 1

          ! Stored value
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up
          call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, jt, j0)
          !call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, it, i0)

          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if ( tau%comp_dn ) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif ( .not.tau%neg_u ) then
                ! 05/15/2012, C.C.:
                ! Note that here we have (.not.neg_u) and (.not.comp_dn).
                ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
                ! So we can safely call DQMC_Gfun_CopyUp() which uses particle-hole symmetry.
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, jt, j0)
          enddo

          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif
          ! Increment index tau%it
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if ( tau%comp_dn ) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif ( .not.tau%neg_u ) then
                ! 05/15/2012, C.C.:
                ! Note that here we have (.not.neg_u) and (.not.comp_dn).
                ! This implies that S%P is defined, i.e. S%checklist(STRUCT_PHASE) = 'T'.
                ! So we can safely call DQMC_Gfun_CopyUp() which uses particle-hole symmetry.
                call DQMC_Gtau_CopyUp(tau)
             endif
             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, jt, j0)
                
          enddo

       enddo

       cnt = cnt + 1
    enddo blocks
 
    if (i0 .ne. tau%nb+1) then
       write(*,*) "Up and down time slices are mismatched. Stop"
       stop
    endif

    sgn = tau%sgnup * tau%sgndn
    do it = 0, L-1
       do i = 1, T1%nClass
          factor = sgn/(T1%F(i)*cnt)
          T1%gnl(i,it,T1%idx)   = T1%gnl(i,it,T1%idx)   + factor*T1%gnl(i,it,T1%tmp)
          T1%chinl(i,it,T1%idx) = T1%chinl(i,it,T1%idx) + factor*T1%chinl(i,it,T1%tmp)
       end do
    end do

    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

    ! Clean up
    T1%gnl(:,:,T1%tmp)   = ZERO
    T1%chinl(:,:,T1%tmp) = ZERO
    
  end subroutine DQMC_TDM1_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, it, i0)
    !
    ! Purpose
    ! =======
    !    This subroutine assembles the time dependent properties
    !    starting from the 1-body Green's function
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    integer, intent(in)          :: it, i0

    ! ... Local scalar ...

    integer  :: i, j, n, k, dt, dt1, dt2
    integer,  pointer :: D(:,:)
    real(wp), pointer :: gnl1(:), chi1(:)
    real(wp), pointer :: gnl2(:), chi2(:)
    real(wp) :: factor

    ! ... Executable ...
    if (.not.T1%compute) return

    dt = it - i0
    if (dt .gt. 0) then
       ! it > i0
       dt1  =  dt
       dt2  =  T1%L - dt
       factor  = 0.25d0
    elseif (dt .lt. 0) then
       ! it < i0
       dt1  =  dt + T1%L
       dt2  = -dt
       factor = -0.25d0
    else
       dt1 = 0
       dt2 = 0
       factor = 0.5d0
    endif

    ! Initialization
    n     =  T1%n
    D     => T1%D
    gnl1  => T1%gnl(:, dt1, T1%tmp)
    gnl2  => T1%gnl(:, dt2, T1%tmp)
    chi1  => T1%chinl(:, dt1, T1%tmp)
    chi2  => T1%chinl(:, dt2, T1%tmp)

    ! Compute Green's function and Chi function
    if (dt .ne. 0) then

       do i = 1, n
          do j = 1, n
             ! k is the distance index of site i and site j
             k = D(i,j)
             gnl1(k)  = gnl1(k) + factor*(upt0(i,j) + dnt0(i,j))
             gnl2(k)  = gnl2(k) - factor*(up0t(i,j) + dn0t(i,j))
             chi1(k)  = chi1(k) - (up0t(j,i)*dnt0(i,j) &
                  + dn0t(j,i)*upt0(i,j))/2
             chi2(k)  = chi2(k) - (up0t(i,j)*dnt0(j,i) &
                  + dn0t(i,j)*upt0(j,i))/2
          end do
       end do

    else

       do i = 1, n
          do j = 1, n
             ! k is the distance index of site i and site j
             k = D(i,j)
             gnl1(k)  = gnl1(k) + factor*(upt0(i,j) + dnt0(i,j))
             chi1(k)  = chi1(k) - up0t(j,i)*dnt0(i,j) &
                  - dn0t(j,i)*upt0(i,j)
          end do
       end do

    endif

  end subroutine DQMC_TDM1_Compute

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Avg(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine average properties in a bin and
    !    increment the bin count (idx).
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer  :: nl, idx
    real(wp) :: factor

    ! ... Executable ...
    if (.not.T1%compute) return
    idx    = T1%idx
    factor = ONE/T1%cnt
    nl     = T1%nClass*T1%L

    ! Compute average on Green's function
    call blas_dscal(nl, factor, T1%gnl(:,:,idx), 1)
    ! Compute average on Chi 
    call blas_dscal(nl, factor, T1%chinl(:,:,idx), 1)

    T1%sgn(idx) = T1%sgn(idx)*factor
    T1%cnt = 0
    T1%idx = T1%idx + 1

  end subroutine DQMC_TDM1_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetErr(T1)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine compute the error in tdm using the jackknife
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer   :: i, j
    integer   :: nproc, n, avg, err, mpi_err
    real(wp)  :: sum_sgn, sgn(T1%nBin), y(T1%nBin), data(T1%nBin)

    if (.not.T1%compute) return
    ! ... Executable ...
    nproc  = qmc_sim%size
    n      = T1%nbin
    avg    = T1%avg
    err    = T1%err

    if (nproc .eq. 1) then

       data = T1%sgn(1:n)
       call DQMC_JackKnife(n, T1%sgn(avg), T1%sgn(err), data , &
            y, sgn, sum_sgn)

       do i = 1, T1%nClass
          do j = 0, T1%L-1
             data =  T1%gnl(i, j, 1:n)
             call DQMC_SignJackKnife(n, T1%gnl(i, j, avg), T1%gnl(i, j, err), &
                  data, y, sgn, sum_sgn)
          enddo
       end do

       do i = 1, T1%nClass
          do j = 0, T1%L-1
             data =  T1%chinl(i, j, 1:n)
             call DQMC_SignJackKnife(n, T1%chinl(i, j, avg), T1%chinl(i, j, err), &
                  data, y, sgn, sum_sgn)
          enddo
       end do

    else

       mpi_err = 0

#      ifdef _QMC_MPI
          
          !Average sign
          call mpi_allreduce(T1%sgn(1), T1%sgn(avg), 1, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Average properties
          call mpi_allreduce(T1%gnl(:,:,1), T1%gnl(:,:,avg), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Average properties
          call mpi_allreduce(T1%chinl(:,:,1), T1%chinl(:,:,avg), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Compute average over n-1 processors
          T1%gnl(:,:,1)   = (T1%gnl(:,:,avg) - T1%gnl(:,:,1)) / dble(nproc - 1)
          T1%chinl(:,:,1) = (T1%chinl(:,:,avg) - T1%chinl(:,:,1)) / dble(nproc - 1)
          T1%sgn(1)   = (T1%sgn(avg) - T1%sgn(1)) / dble(nproc - 1)

          !Store average amongst all processors
          T1%gnl(:,:,avg)   = T1%gnl(:,:,avg) / T1%sgn(avg) 
          T1%chinl(:,:,avg) = T1%chinl(:,:,avg) / T1%sgn(avg) 
          T1%sgn(:,avg)     = T1%sgn(:,avg) / dble(nproc)

          !Store jackknife in the processor bin
          T1%gnl(:,:,1)   = T1%(:,:,1) / T1%sgn(1) 

          !Compute error
          call mpi_allreduce(T1%gnl(:,:,1)**2, T1%gnl(:,err), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)
          T1%gnl(:,:,err) = T1%gnl(:,:,err) / dble(nproc) - T1%gnl(:,:,avg)**2 
          T1%gnl(:,:,err) = sqrt(T1%gnl(:,:,err) * dble(nproc-1))

          call mpi_allreduce(T1%chinl(:,:,1)**2, T1%chinl(:,err), n, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)
          T1%chinl(:,:,err) = T1%chinl(:,:,err) / dble(nproc) - T1%chinl(:,:,avg)**2 
          T1%chinl(:,:,err) = sqrt(T1%chinl(:,:,err) * dble(nproc-1))

#      endif

    endif

  end subroutine DQMC_TDM1_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Print(T1, S, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    type(struct), intent(in) :: S
    integer, intent(in)      :: OPT

    integer             :: i, j
    real(wp)            :: tmp(T1%L, 2)
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
    enddo

    do i = 1, T1%nclass
       do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%gnl(i, j, T1%avg:T1%err)
       enddo
       title="GFun "//trim(adjustl(S%clabel(i)))
       call DQMC_Print_RealArray(0, T1%L , title, label, &
          tmp(:, 1:1), tmp(:, 2:2), OPT)
       write(OPT,'(1x)')
    enddo

    do i = 1, T1%nclass
       do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%chinl(i, j, T1%avg:T1%err)
       enddo
       title="S+S- "//trim(adjustl(S%clabel(i)))
       call DQMC_Print_RealArray(0, T1%L , title, label, &
          tmp(:, 1:1), tmp(:, 2:2), OPT)
       write(OPT,'(1x)')
    enddo

  end subroutine DQMC_TDM1_Print

  !--------------------------------------------------------------------!

end module DQMC_TDM1
