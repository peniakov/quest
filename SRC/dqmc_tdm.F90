module DQMC_TDM
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_Struct
  use DQMC_TDM1
  use DQMC_TDM2
  use DQMC_Gtau
  use DQMC_CFG

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
  type TDM
     type(TDM1) :: T1
     type(TDM2) :: T2
     type(gtau) :: tau

     integer    :: nTausk
     integer    :: which

     ! for average
     integer    :: L
     integer    :: nBin
     integer    :: itvl
     
     real(wp)   :: dtau
     real(wp), pointer   :: sgn(:)
     integer    :: avg, err, idx, cnt
     
     ! for Fourier transformation
     real(wp), pointer :: CS(:,:)       ! FT matrix for space
     logical    :: FTspace
     character(30), pointer :: clabel(:)
     character(30), pointer :: wlabel(:)
  end type TDM
  
  integer, parameter :: TDM_NONE = 0
  integer, parameter :: TDM_TDM1 = 1
  integer, parameter :: TDM_TDM2 = 2
  integer, parameter :: TDM_ALL  = 3

contains

  !--------------------------------------------------------------------!
  
  !subroutine DQMC_TDM_Init(which, tm, S, WS, B_up, B_dn, cfg)
  !
  ! 05/08/2012: Fixed by C.C.
  subroutine DQMC_TDM_Init(which, Hub, tm, S, WS, B_up, B_dn, cfg)
  use dqmc_hubbard
    !
    ! Purpose
    ! =======
    !    This subroutine initializes time dependent measurements. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)       :: which   ! which tdm to perform
    type(TDM), intent(inout)  :: tm      ! time dependent measurement
    type(Struct), intent(in)  :: S
    type(MatB), intent(in)    :: B_up, B_dn
    type(WSpace), intent(in)  :: WS
    type(Config), intent(in)  :: cfg
    type(Hubbard), intent(in) :: Hub

    ! ... Local Variables ...
    integer  :: n, L, nOrth, nWrap, nBin, itvl
    real(wp) :: dtau

    ! ... Executable ...
    
    call CFG_Get(cfg,  "n", n)
    call CFG_Get(cfg,  "L", L)
    call CFG_Get(cfg,  "north", nOrth)
    call CFG_Get(cfg,  "nwrap", nWrap)
    call CFG_Get(cfg,  "nbin",  nBin)
    call CFG_Get(cfg,  "dtau",  dtau)
    call CFG_Get(cfg,  "nitvl", itvl)
    tm%which = which
    
    tm%L    = L
    tm%nBin = nBin
    allocate(tm%sgn(nBin+2))
    tm%avg = nBin + 1
    tm%err = nBin + 2
    tm%idx = 1
    tm%cnt = 0
    tm%itvl = itvl
    
    tm%clabel  => S%clabel
    tm%wlabel => S%wlabel
    tm%dtau = dtau

    ! initialize TDM1 and TDM2, and gtau
    if (which .gt. TDM_NONE) then
       !call DQMC_TDM1_Init(n, L, tm%T1, nBin, S, tm%avg)
       call DQMC_TDM1_Init(n, L, tm%T1%dtau, tm%T1, nBin, S)   ! 05/08/2012: Fixed by C.C 
    end if

    if (which .gt. TDM_TDM1) then
       if (S%checklist(STRUCT_WAVE) .and. S%checklist(STRUCT_BOND)) then
          call DQMC_TDM2_Init(n, L, nBin, tm%T2, S, WS%R4)
       else
          tm%which = TDM_TDM1
       end if
    end if
    
    if (which .gt. TDM_NONE) then
       !call DQMC_Gtau_Init(n, L, TAU_BOTH, nOrth, nWrap, tm%tau, B_up, B_dn, WS)
       call DQMC_Gtau_Init(Hub,tm%tau)
    end if
    
    ! allocate for FT
    tm%FTspace = S%checklist(STRUCT_FT)
    if (tm%FTspace) then
       tm%CS => S%FT
    else
       nullify(tm%CS)
    end if

  end subroutine DQMC_TDM_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Free(tm)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: tm      ! TDM to be freed

    ! ... Executable ...

    call DQMC_Gtau_Free(tm%tau)
    if (tm%which .gt. TDM_NONE) then
       call DQMC_TDM1_Free(tm%t1)
       if (tm%which .gt. TDM_TDM1) then
          call DQMC_TDM2_Free(tm%t2)
       end if
    end if

    ! free FT space
    ! deallocate(tm%data, tm%FTM, tm%Wspl)

  end subroutine DQMC_TDM_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_Make_FTM(L, itvl, dtau, FTM, base)
    !
    ! Purpose
    ! =======
    !    This subroutine constructs Fourier transformation matrix (FTM).
    !    It implements numerical itegration on FT.
    !
    ! Arguments
    ! =========
    !
    real(wp), intent(inout)   :: FTM(:,:)
    real(wp), intent(in)      :: base, dtau
    integer, intent(in)       :: L
    integer, intent(in)       :: itvl

    ! ... parameters ...
    real(wp), parameter  :: ONETHIRD  = ONE/3.0_wp    
    real(wp), parameter  :: TWOTHIRD  = TWO/3.0_wp
    real(wp), parameter  :: FOURTHIRD = 4.0_wp/3.0_wp    
    real(wp), parameter  :: TWOPI     = 6.283185307179586
    complex(wp), parameter  :: EI     = cmplx(ZERO, ONE)

    ! ... Local variables ...
    integer  :: i, j, img_idx, longL
    real(wp) :: r
    complex(wp) :: coef, omega

    
    ! ... Executable ...
    
    r       = dtau / itvl
    img_idx = L/2 + 1
    longL   = itvl*L+1

    do j = 1, L/2+1
       ! The first term
       FTM(1, j) = ONETHIRD * r

       ! The middle terms
       omega = EI*(j-base)*TWOPI/L/itvl
       do i = 2, longL-1, 2
          coef = FOURTHIRD * exp(omega*(i-1)) * r
          FTM(i,j)         = real (coef)          ! real part
          FTM(i,j+img_idx) = aimag(coef)          ! imagianry part

          coef =  TWOTHIRD * exp(omega*i) * r
          FTM(i+1,j)         = real (coef)        ! real part
          FTM(i+1,j+img_idx) = aimag(coef)        ! imagianry part
       end do

       ! The last second term
       if (mod(longL-1,2) .ne. 0) then
          coef = FOURTHIRD * exp(omega*(longL-2)) * r
          FTM(longL-1, j)        = real (coef)    ! real part
          FTM(longL-1,j+img_idx) = aimag(coef)    ! imagianry part
       end if
       
       ! The last term
       coef = ONETHIRD * exp(omega*(longL-1)) * r
       FTM(longL, j)        = real (coef)         ! real part
       FTM(longL,j+img_idx) = aimag(coef)         ! imagianry part
    end do

  end subroutine DQMC_Make_FTM

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Meas(tm, G_up, G_dn)
    !
    ! Purpose
    ! =======
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout)    :: tm
    type(G_fun), intent(inout)  :: G_up, G_dn
    

    ! ... Local var ...
    integer :: i, j, k, m, L, it, itplus, itminus
    real(wp):: sgn
    real(wp),pointer :: up0t(:,:), upt0(:,:), dn0t(:,:), dnt0(:,:)

    ! ... executable ...
    L = tm%tau%L
    m = tm%tau%nb
    
    upt0 => tm%tau%upt0
    up0t => tm%tau%up0t
    dnt0 => tm%tau%dnt0
    dn0t => tm%tau%dn0t

    itplus   = tm%tau%north/2
    itminus  = itplus
    if (mod(tm%tau%north,2) .eq.  0) itminus = itminus + 1

    do i = 1, m

       do j = 0, m-1

          k = i + j
          if (k .lt. m) k = k - m

          !Stored value
          !call DQMC_Gtau_LoadA(tm%tau, k, i)
          !call DQMC_TDM1_Meas(tm%T1, upt0, up0t, dnt0, dn0t, tm%tau%ii, tm%tau%ib)
          call DQMC_Gtau_LoadA(tm%tau, k, i, sgn)  ! 05/08/2012: Fixed by C.C.
          call DQMC_TDM1_Meas(tm%T1, tm%tau)

          do it = 1, itminus
             !call DQMC_change_gtau_time(2, tau, G_up, G_dn)
             !call DQMC_TDM1_Meas(tm%T1, upt0, up0t, dnt0, dn0t, tm%tau%ii, tm%tau%ib)
             !
             ! 05/08/2012: Fixed by C.C.
             call DQMC_change_gtau_time(tm%tau, 2, k)
             call DQMC_TDM1_Meas(tm%T1, tm%tau)
          enddo

          call DQMC_Gtau_LoadA(tm%tau, k, i, sgn) ! 05/08/2012: Fixed by C.C.
          do it = 1, itplus
             !call DQMC_change_gtau_time(1, tau, G_up, G_dn)
             !call DQMC_TDM1_Meas(tm%T1, upt0, up0t, dnt0, dn0t, tm%tau%ii, tm%tau%ib)
             !
             ! 05/08/2012: Fixed by C.C.
             call DQMC_change_gtau_time(tm%tau, 1, k)
             call DQMC_TDM1_Meas(tm%T1, tm%tau)
          enddo

       enddo
    end do

    sgn = G_up%sgn * G_dn%sgn

    call DQMC_TDM1_Acc(tm%T1, sgn, tm%idx)
    if (tm%which .gt. TDM_TDM1) then
       call DQMC_TDM2_Acc(tm%T2,sgn)
    end if

    tm%sgn(tm%idx) =  tm%sgn(tm%idx) + sgn
    tm%cnt = tm%cnt + 1
    
  end subroutine DQMC_TDM_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Avg(tm)
    !
    ! Purpose
    ! =======
    !    This subroutine averges the pair measurements.
    !    which is stored in T1%sus2.
    !    The averaging process runs as follows. 
    !    1. Let T = sus2*W
    !    2. Compute b(i) = W(:,i)'*T(:,i)
    !    3. Averaging b  = b/nMeas 
    !    4. Store b into sus3(:,idx)
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: tm
    
    ! ... Local variables ...
    integer  :: idx
    real(wp) :: factor

    ! ... Executable ...

    idx = tm%idx
    factor = ONE/tm%cnt
    
    !call DQMC_TDM1_Avg(tm%T1, idx, factor)
    call DQMC_TDM1_Avg(tm%T1)    ! 05/08/2012: Fixed by C.C.
    if (tm%which .gt. TDM_TDM1) then
       call DQMC_TDM2_Avg(tm%T2, tm%T1%gnl(:,:,idx), idx, factor)
    end if

    tm%sgn(idx) = tm%sgn(idx)*factor
    tm%cnt = 0
    tm%idx = tm%idx + 1

  end subroutine DQMC_TDM_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_Spline(n, L, dtau, input, output, W, itvl)
    !  
    ! Purpose
    ! =======
    !    This subroutine computes spline interpolation of a set of data.
    !    Fx is a nxL matrix where L is the time slice. For the details
    !    of spline interpolation, see [1]
    !
    !    1. Solve a tridiagonal linear system (L-2)x(L-2)
    !       The following fomula only shows the first row of Fx. 
    !       Let F(i) = Fx(1,i)
    !       
    !                 [4 1       ]          [   F(1)-2F(2)+F(3)  ]
    !                 [1 4 1     ]      6   [   F(2)-2F(3)+F(4)  ]
    !               X [  ......  ] = ------ [          :         ]
    !                 [     1 4 1]   dtau^2 [          :         ]
    !                 [       1 4]          [ F(L-2)-2F(L-1)+F(L)]
    !
    !    2. Formulate a_i*dx^3 + b_i*dx^2+c_i*dx+d_i
    !       where a_i = (X_{i+1}-M_i)/(6*dtau)
    !             b_i = X_i/2
    !             c_i = (F(i+1)-F(i))/dtau - (X_{i+1}+2X_{i})*dtau/6
    !             d_i = F(i)
    !
    !   [1] Van Loan, Charles F. Introduction to Scientific Computing. 
    !       1997 New Jersey: Prentice Hall
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: n, L
    real(wp), intent(in)    :: dtau
    real(wp), intent(in)    :: input(:,:)         ! original data
    real(wp), target, intent(inout) &
                      :: output(n,itvl*(L-1)+1)   ! interpolated points
    real(wp), intent(inout) :: W(:)               ! working space
    integer,  intent(in)    :: itvl               ! itvl must >= 2

    ! ... Local ...
    real(wp) :: r, h               ! temp variable
    integer  :: i, j, k
    real(wp) :: a, b, c, d         ! coefficient of spline
    real(wp), pointer :: Pt(:,:)

    ! ... Executable ...
    
    ! We use the Thomas algorithm to solve the tridiagonal linear
    ! system.

    !! Reuse the last block of the output array as a working space
    Pt => output(1:n, (itvl-1)*(L-1)+1:itvl*(L-1)+1)

    r = 6.0_wp/dtau/dtau
    do i = 1,L-2
       do j = 1,n
          Pt(j,i+1) = r*(input(j,i)-TWO*input(j,i+1)+input(j,i+2))
       end do
    end do

    ! forward subsitution
    W(2) = 4.0_wp
    do i = 3,L-1
       r = ONE/W(i-1)
       W(i) = 4.0_wp - r
       do j = 1, n
          Pt(j,i) = Pt(j,i) - r*Pt(j,i-1)
       end do
    end do

    ! back subsitution
    r = ONE/W(L-1)
    do j = 1, n
       Pt(j, 1)   = ZERO
       Pt(j, L)   = ZERO
       Pt(j,L-1)  = r*Pt(j,L-1)
    end do

    do i = L-2,2,-1
       r = ONE/W(i)
       do j = 1, n
          Pt(j,i) = r*(Pt(j,i)-Pt(j,i+1))
       end do
    end do
    
    ! Find interpolated points

    ! There are L-1 curves for L points
    r = dtau/itvl        ! step size
    do i = 1, L-1
       do j = 1, n
          a = (Pt(j,i+1)-Pt(j,i))/6.0_wp/dtau
          b = Pt(j,i)/TWO
          c = (input(j,i+1)-input(j,i))/dtau-&
               (Pt(j,i+1)+2*Pt(j,i))*dtau/6.0_wp
          d = input(j,i)
          
          ! The Ft(j,itvl*(i-1)+1) = input(j,i) 
          output(j,itvl*(i-1)+1) = input(j,i)
          h = r
          do k = itvl*(i-1)+2, itvl*i
             output(j,k) = ((a*h+b)*h+c)*h+d
             h = h + r
          end do
       end do
    end do

    ! The Ft(:,itvl*L) = Pt(:,L) 
    output(:,itvl*(L-1)+1) = input(:,L)

  end subroutine DQMC_Spline

  !--------------------------------------------------------------------!

  subroutine DQMC_DCT_Space(n, L, CS, data, dct)
    !  
    ! Purpose
    ! =======
    !    This subroutine computes DCT for given data.
    !
    !                           
    !  dct(i,k)=sum_{j=1:n}data(j,k)*cos(2*pi*x(i)*xx(j)/nx+2*pi*y(i)*yy(j)/ny)
    !   
    !    where x, y, xx, yy are indices for sites in the lattice.
    !    The element (i,j) of array CS contains the factor 
    !
    !    cos(2*pi*x(i)*xx(j)/nx+2*pi*y(i)*yy(j)/ny) 
    !
    ! Preassumption
    ! =============
    !   1. for 2d rectangular data only. 
    !   2. The domain must be rectangular.
    !   3. cstbl must store "correct" cos data
    !   4. If want to use this fucntion, call dqmc_init_2dperl_adv to 
    !      initialize the data.
    !   
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: n, L
    real(wp), intent(in)    :: data(:, :)
    real(wp), intent(in)    :: CS(:, :)
    real(wp), intent(inout) :: dct(:, :)

    call blas_dgemm('N', 'N', n, L, n, ONE, CS, n, data, n, ZERO, dct, n)

  end subroutine DQMC_DCT_Space

  !--------------------------------------------------------------------!

  subroutine DQMC_FT_Time(n, L, data, FT, W1, W2, W3, itvl, dtau)
    !  
    ! Purpose
    ! =======
    !    This subroutine computes FT on given data.
    !  
    !   
    ! Arguments
    ! =========
    !
    integer,  intent(in)    :: n, L          ! number of sites and time slices
    real(wp), intent(in)    :: data(:, :)    ! nClass*longL, input data
    real(wp), intent(inout) :: FT(:, :)      ! longL*(L+2), output of FT
    real(wp), intent(inout) :: W1(:, :)      ! nClass*longL, for spline output
    real(wp), intent(inout) :: W2(:, :)      ! longL*(L+2), Fourier matrix 
    real(wp), intent(inout) :: W3(:)         ! working space: dim=L+1
    real(wp), intent(in)    :: dtau          ! step size
    integer,  intent(in)    :: itvl          ! spline interval

    ! ... local varaibles ...
    integer  :: longL

    ! ... Executable ...
    longL    = L*itvl+1
    call DQMC_Spline(n, L+1, dtau, data, W1, W3, itvl)

    ! FT
    call blas_dgemm('N', 'N', n, L+2, longL, ONE, W1, n, W2, longL, ZERO, FT, n)
    
  end subroutine DQMC_FT_Time

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Postprocessing(tm, OPT)
    !  
    ! Purpose
    ! =======
    !    This subroutine post process TDM1 and TDM2.
    !    1. Compute Fourier transfomation
    !    2. Compute avg and error of measurements
    !    3. If OPT >=0 then print out the results
    !       Otherwise, store the results
    !  
    ! Arguments
    ! =========
    !  
    type(TDM), intent(inout)    :: tm
    integer, intent(in)         :: OPT

    ! ... Local variables ...
    integer :: longL, L, itvl, nClass, nWave, nBin, L2
    integer :: bin, avg, err
    real(wp) :: dtau
    real(wp), allocatable :: W1(:,:), FTM(:,:), W2(:), DCT(:,:,:), FT(:,:,:)
    logical :: FTS, bPrint
    real(wp):: sgn(tm%nBin), y(tm%nBin), sum_sgn ! for statistics


    ! ... Executable ...

    L      = tm%L
    L2     = L/2+1
    itvl   = tm%itvl
    longL  = L*itvl+1 
    nClass = tm%T1%nClass
    nWave  = tm%T2%nWave
    
    dtau   = tm%dtau
    nBin   = tm%nBin
    avg    = tm%avg
    err    = tm%err

    ! allocate space for FT
    allocate(W1 (nClass, longL))
    allocate(FTM(longL, L+2))
    allocate(W2 (L+1))
    allocate(FT (max(nClass,nWave), L+2, nBin+2))
    FTS = tm%FTSpace
    ! if CS matrix is provided, FT K space
    if (FTS) then
       allocate(DCT(nClass, L+1, nBin+2))
    end if
    bPrint = (OPT .ge. 0)

    ! print header
    if (bPrint) then
       write(OPT,FMT_DBLINE) 
       write(OPT,*) "TIME DEPENDENT MEASUREMENTS"
       write(OPT,FMT_DBLINE) 
    end if

    ! Average sign, using Jackknif
    call DQMC_JackKnife(tm%nBin, tm%sgn(tm%avg), tm%sgn(tm%err), &
         tm%sgn(1:tm%nBin), y, sgn, sum_sgn)

    ! =========================================================
    ! make FTM for Green's function
    call DQMC_Make_FTM(L, itvl, dtau, FTM, HALF)
        
    ! G_NL    
    call DQMC_TDM_Err(nClass, nBin, L+1, tm%T1%gnl, y, &
         sgn, avg, err, sum_sgn)
    if (bPrint) then
       call DQMC_Print_RealArray(nClass, L+1, "G(nx,ny,ti)", &
            tm%clabel, tm%T1%gnl(:,:,avg), tm%T1%gnl(:,:,err) , OPT)
    end if

    ! G_KL
    if (FTS) then ! FT on space
       do bin = 1, nBin
          call DQMC_DCT_Space(nClass, L+1, tm%CS, &
               tm%T1%gnl(:,:,bin), DCT(:,:,bin))
       end do
       call DQMC_TDM_Err(nClass, nBin, L+1, DCT, y, &
            sgn, avg, err, sum_sgn)
       if (bPrint) then
          call DQMC_Print_RealArray(nClass, L+1, "G(qx,qy,ti)", &
               tm%clabel, DCT(:,:,avg), DCT(:,:,err) , OPT)
       end if
    end if
    
    ! G_NW
    do bin = 1, nBin
       call DQMC_FT_Time(nClass, L, tm%T1%gnl(:,:,bin),  &
            FT(:,:,bin), W1, FTM, W2, itvl, dtau)
    end do
    call DQMC_TDM_Err(nClass, nBin, L+2, FT, y, &
         sgn, avg, err, sum_sgn)    
    if (bPrint) then
       call DQMC_Print_ComplexArray(nClass, L2, "G(nx,ny,omega)", &
            tm%clabel, FT(:,:,avg), FT(:,:,err) , OPT)
    end if


    ! G_QW
    if (FTS) then
       do bin = 1, nBin
          call DQMC_FT_Time(nClass, L, DCT(:,:,bin),  &
               FT(:,:,bin), W1, FTM, W2, itvl, dtau)
       end do
       call DQMC_TDM_Err(nClass, nBin, L+2, FT, y, &
            sgn, avg, err, sum_sgn)    
       if (bPrint) then
          call DQMC_Print_ComplexArray(nClass, L2, "G(qx,qy,omega)", &
               tm%clabel, FT(:,:,avg), FT(:,:,err) , OPT)
       end if
    end if

    ! =========================================================
    ! make FTM for Chi function
    call DQMC_Make_FTM(L, itvl, dtau, FTM, ONE)

    ! chi_NL    
    call DQMC_TDM_Err(nClass, nBin, L+1, tm%T1%chinl, y, &
         sgn, avg, err, sum_sgn)
    if (bPrint) then
       call DQMC_Print_RealArray(nClass, L+1, "chi(nx,ny,ti)", &
            tm%clabel, tm%T1%chinl(:,:,avg), tm%T1%chinl(:,:,err) , OPT)
    end if

    ! chi_KL
    if (FTS) then ! FT on space
       do bin = 1, nBin
          call DQMC_DCT_Space(nClass, L+1, tm%CS, &
               tm%T1%chinl(:,:,bin), DCT(:,:,bin))
       end do
       call DQMC_TDM_Err(nClass, nBin, L+1, DCT, y, &
            sgn, avg, err, sum_sgn)
       if (bPrint) then
          call DQMC_Print_RealArray(nClass, L+1, "chi(qx,qy,ti)", &
               tm%clabel, DCT(:,:,avg), DCT(:,:,err) , OPT)
       end if
    end if
    
    ! chi_NW
    do bin = 1, nBin
       call DQMC_FT_Time(nClass, L, tm%T1%chinl(:,:,bin),  &
            FT(:,:,bin), W1, FTM, W2, itvl, dtau)
    end do
    call DQMC_TDM_Err(nClass, nBin, L+2, FT, y, &
         sgn, avg, err, sum_sgn)    
    if (bPrint) then
       call DQMC_Print_ComplexArray(nClass, L2, "chi(nx,ny,omega)", &
            tm%clabel, FT(:,:,avg), FT(:,:,err) , OPT)
    end if

    ! chi_QW
    if (FTS) then
       do bin = 1, nBin
          call DQMC_FT_Time(nClass, L, DCT(:,:,bin),  &
               FT(:,:,bin), W1, FTM, W2, itvl, dtau)
       end do
       call DQMC_TDM_Err(nClass, nBin, L+2, FT, y, &
            sgn, avg, err, sum_sgn)    
       if (bPrint) then
          call DQMC_Print_ComplexArray(nClass, L2, "chi(qx,qy,omega)", &
               tm%clabel, FT(:,:,avg), FT(:,:,err) , OPT)
       end if
    end if


    ! =========================================================
    ! check if to process TDM2
    if (tm%which .gt. TDM_TDM1) then
       ! with vertex sus
       call DQMC_TDM_Err(nWave, nBin, L+1, tm%T2%Bpair, y, sgn, &
            avg, err, sum_sgn)
       call DQMC_TDM_Err(nWave, nBin, L+1, tm%T2%Npair, y, sgn, &
            avg, err, sum_sgn)
       if (bPrint) then
          call DQMC_Print_RealArray(nWave, L, "P_dd'(L)  (vertex)", &
               tm%wlabel, tm%T2%Bpair(:,:,avg), tm%T2%Bpair(:,:,err) , OPT)
          call DQMC_Print_RealArray(nWave, L, "P_dd'(L)  (non-vertex)", &
               tm%wlabel, tm%T2%Npair(:,:,avg), tm%T2%Npair(:,:,err) , OPT)
       end if

       ! vertex one
       do bin = 1, nBin
          call DQMC_FT_Time(nWave, L, tm%T2%Bpair(:,:,bin), &
               FT(:,:,bin), W1, FTM, W2, itvl, dtau)
       end do
       call DQMC_TDM_Err(nWave, nBin, L+2, FT, y, &
            sgn, avg, err, sum_sgn)    
       if (bPrint) then
          call DQMC_Print_ComplexArray(nWave, L2, "P_dd'(w)  (vertex)", &
               tm%wlabel, FT(:,:,avg), FT(:,:,err) , OPT)
       end if

       ! nonvertex one
       do bin = 1, nBin
          call DQMC_FT_Time(nWave, L, tm%T2%Npair(:,:,bin), &
               FT(:,:,bin), W1, FTM, W2, itvl, dtau)
       end do
       call DQMC_TDM_Err(nWave, nBin, L+2, FT, y, &
            sgn, avg, err, sum_sgn)    
       if (bPrint) then
          call DQMC_Print_ComplexArray(nWave, L2, "P_dd'(w)  (nonvertex)", &
               tm%wlabel, FT(:,:,avg), FT(:,:,err) , OPT)
       end if
    end if

    ! =========================================================

    deallocate(W1, W2, FTM, FT)
    if (FTS) then
       deallocate(DCT)
    end if
    
  end subroutine DQMC_TDM_Postprocessing

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Err(n, nBin, L, G, y, ts, avg, err, denom)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the average and errors
    !    of time dependent measurements. 
    !
    ! Argument
    ! ========
    !
    integer, intent(in)     :: n, nBin, L, avg, err
    real(wp), intent(inout), target :: G(n,L,nBin+2)
    real(wp), intent(inout) :: y(nBin)
    real(wp), intent(in)    :: denom, ts(nBin)
    
    ! ... Local scalar
    integer :: i, j
    real(wp), pointer :: gpt(:)

    do j = 1, L
       do i = 1, n
          gpt => G(i,j,1:nBin)
          call DQMC_SignJackKnife(nBin, G(i,j,avg), G(i,j,err), &
               gpt, y, ts, denom)
       end do
    end do

  end subroutine DQMC_TDM_Err

  !--------------------------------------------------------------------!

end module DQMC_TDM
