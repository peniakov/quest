module DQMC_GFun
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE
  use _DQMC_MATB
  use DQMC_SEQB

  implicit none 

  ! 
  ! This module defines the data type and subroutines for 
  ! constructiong and manipulating the Green's function.
  ! The Green's function in Hubbard's model is defined as
  ! 
  !     G = inv((I+B_{i-1}...B_1*B_l...B_{i})),   
  !
  ! where B_{i} = V_{i}*B.  Matrix V_{i} is diagonal and B is
  ! defined as exp(-t*K), where K is an adjacent matrix for
  ! underline lattice. 
  !
  !     G = inv((I+inv(B_{i-1}...B_1*B_l...B_{i}))),   
  !
  ! For more details, see the working note.
  !
  ! [1] Z. Bai, W.Chen, R. Scalettar, I. Yamazaki, "Lecture Notes 
  !     on Advances of Numerical Methods for Hubbard Quantum Monte
  !     Carlo Simulation." 
  !
  !  Subroutine List
  !  ===============
  !    DQMC_GFun_Init(n, l, GF) : initialize the data type.
  !    DQMC_GetV(G_up, G_dn, explook, hub) : construct matrix V
  !    DQMC_GetG(sl, G, R, Q, T, rw, D, 
  !              tau, piv1, piv2, lwork, nOrth) : construct matrix G
  !
  !  Data Type
  !  =========
  !
  type G_fun
     integer  :: n                      ! Number of sites in Hubbard's model
     integer  :: L                      ! Number of slice
     real(wp) :: sgn                    ! Sign of det(G)
     real(wp), pointer :: G(:,:)        ! Matrix of Green's function
     logical  :: owns_G                 ! does this G_fun own its G's storage?
     real(wp), pointer :: V(:,:)        ! Matrix info of V(1)...V(l)
     integer  :: ilb                    ! index of left most B

     ! working space, in addition to SB's
     integer, pointer  :: pvt(:)
     real(wp),pointer  :: tmp(:,:)

     ! For numerical stab
     integer  :: nWrap, wps, lastwr
     real(wp) :: difflim, errrate
     integer  :: redo, noredo

     ! Block size of delayed update
     integer  :: nBlk   
     integer  :: blkSz
     real(wp),pointer  :: U(:,:), W(:,:)
     integer  :: nModify

  end type G_fun

  logical, parameter :: GMAT_UP = .true.
  logical, parameter :: GMAT_DN = .false.

contains

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Init(n, L, G, V, WS, nWrap, difflim, errrate, up)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n              ! Number of sites
    integer, intent(in)         :: L              ! Number of slice
    type(G_fun), intent(inout)  :: G              ! Green's function
    real(wp), intent(in), target:: V(n,L)         ! V matrix
    type(Wspace), intent(in), target :: WS        ! working space
    real(wp), intent(in)        :: difflim        ! limit of mat diff
    real(wp), intent(in)        :: errrate        ! toleratble error
    integer, intent(in)         :: nwrap          ! safe wrap number
    logical, intent(in)         :: up             ! for G_up or G_dn

    ! ... LAPACK subroutine
    integer  :: ILAENV

    ! ... Local variables ...
    integer  :: nb

    ! ... Executable ...

    G%n   = n
    G%L   = L
    G%sgn = ONE
    G%pvt => WS%I2
    G%V   => V

    G%ilb = -L-1
    G%difflim = difflim
    G%nWrap   = nWrap
    G%wps     = 1
    G%tmp => WS%R8
    G%errrate = errrate

    G%redo   = 0
    G%noRedo = 1
    G%lastwr = 0
    G%blkSz  = 0
    G%nModify = 0
    
    nb = ILAENV( 1, "DGETRI", " ", n, -1, -1, -1 )
    G%nblk = min(nb, n)

    allocate(G%G(n,n))
    G%owns_G = .true.
    G%sgn = ONE

    if (up) then
       G%U   => WS%R1
       G%W   => WS%R2
    else
       G%U   => WS%R3
       G%W   => WS%R4
    end if
    
  end subroutine DQMC_GFun_Init

  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Init_original(n, L, G, V, WS, nWrap, difflim, errrate, up)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: n              ! Number of sites
    integer, intent(in)         :: L              ! Number of slice
    type(G_fun), intent(inout)  :: G              ! Green's function
    real(wp), intent(in), target:: V(n,L)         ! V matrix
    type(Wspace), intent(in), target :: WS        ! working space
    real(wp), intent(in)        :: difflim        ! limit of mat diff
    real(wp), intent(in)        :: errrate        ! toleratble error
    integer, intent(in)         :: nwrap          ! safe wrap number
    logical, intent(in)         :: up             ! for G_up or G_dn

    ! ... LAPACK subroutine
    integer  :: ILAENV

    ! ... Local variables ...
    integer  :: nb

    ! ... Executable ...

    G%n   = n
    G%L   = L
    G%sgn = ONE
    G%pvt => WS%I2
    G%V   => V

    G%ilb = -L-1
    G%difflim = difflim
    G%nWrap   = nWrap
    G%wps     = nWrap
    G%tmp => WS%R8
    G%errrate = errrate

    G%redo   = 0
    G%noRedo = 1
    G%lastwr = 0
    G%blkSz  = 0
    G%nModify = 0
    
    nb = ILAENV( 1, "DGETRI", " ", n, -1, -1, -1 )
    G%nblk = min(nb, n)

    allocate(G%G(n,n))
    G%sgn = ONE

    if (up) then
       G%U   => WS%R1
       G%W   => WS%R2
    else
       G%U   => WS%R3
       G%W   => WS%R4
    end if
    
  end subroutine DQMC_GFun_Init_original


  !-------------------------------------------------------------------------!

  subroutine DQMC_GFun_Clone(G1, G2)
    !
    ! Purpose
    ! =======
    !    G1 = G2
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G1              ! Green's function
    type(G_fun), intent(in)     :: G2              ! Green's function 

    ! ... Executable ...

    G1%n   = G2%n
    G1%L   = G2%L

    G1%sgn = G2%sgn
    G1%G   => G2%G
    G1%owns_G = .false.
    G1%V   => G2%V

  end subroutine DQMC_GFun_Clone

  !-------------------------------------------------------------------------!

  subroutine DQMC_Gfun_Free(G)
    !
    ! Purpose
    ! =======
    !    This subroutine frees memory of G.
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout)  :: G  ! Green's function

    ! ... Executable ...

    if (G%owns_G) then
       deallocate(G%G)
    end if

  end subroutine DQMC_Gfun_Free
  
  !--------------------------------------------------------------------------!

  subroutine DQMC_GetG(il, G, SB)
    !
    ! Purpose
    ! =======
    !    This subroutine returns 
    !    
    !         G = inv((I+B_{il}...B_1*B_L...B_{il+1})),
    !
    !    where B_{i} = V_{i}*B. Matrix V_{i} is diagonal,
    !    whose elements are stored in a vector V.
    !    It also returns the sign of det(G).
    !
    ! Arguments
    ! =========
    !
    integer,  intent(in)       :: il            ! starting slice
    type(G_fun), intent(inout) :: G
    type(SeqB), intent(inout)  :: SB

    ! ... local scalars    ...
    logical  :: match
    real(wp) :: diff, det, sgn

    ! ... Executable ...

    match = ((G%ilb .eq. G%L .and. il .eq. 1) .or. (G%ilb .eq. il-1)) 

    if (match) then
       ! Swap the slice of G_up
       call DQMC_MultB_Left  (G%n, G%G, SB%B, G%V(:,il), SB%W1)
       call DQMC_MultBi_Right(G%n, G%G, SB%B, G%V(:,il), SB%W1)

       ! check if the swapping is safe
       G%wps = G%wps - 1

       ! need to recompute
       if (G%wps .eq. 0) then
          ! store current result
          ! G%tmp = G%G
          call blas_dcopy(G%n*G%n, G%G, 1, G%tmp, 1)
          ! recompute G
          call DQMC_ComputeG(il, G%n, sgn, G%G, G%V, SB, G%pvt, .false., det)

          G%wps = G%nWrap
          G%sgn = sgn

          ! evaluate the difference
          diff = DQMC_MatDiff(G%n, G%tmp, G%G)
          
          ! increase the statistics
          if(diff .gt. G%difflim) then
             G%redo = G%redo + 1
          else
             G%noredo = G%noredo + 1
          end if          
       end if
    else
       ! cannot use swap, recompute
       call DQMC_ComputeG(il, G%n, G%sgn, G%G, G%V, SB, G%pvt, .false., det)
       G%wps = G%nWrap
    end if

    G%ilb = il
       
  end subroutine DQMC_GetG

  !--------------------------------------------------------------------------!

  subroutine DQMC_ComputeG(il, n, sgn, G, V, SB, pvt2, compDet, det)
    !
    ! Purpose
    ! =======
    !    This subroutine computes the matrix 
    !    
    !         G = inv((I+B_{il}...B_1*B_l...B_{il+1})),
    !
    !    where B_{i} = V_{i}*B. Matrix V_{i} is diagonal,
    !    whose elements are stored in a vector V.
    !    It also returns the sign of det(G).
    !
    !    If check is true, it will compare the previous G and 
    !    the new computed one, and put the difference to diff.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix G, B, work are square with the same order n.
    !
    ! Steps
    ! =====
    !
    ! 1. Compute A = B_{i-1}...B_1*B_l...B_{i}.
    !    To keep the multiplication stable, A is decomposed as
    !        
    !       A = U*D*T,
    !
    !    where U is an orthogonal matrix, D is a diagonal matrix.
    !    This UDT decompisition is performed in every nOrth steps. 
    !
    ! 2. Compute C = inv(U)*inv(T) + D.
    !
    ! 3. The final G is computed as inv(T)*inv(C)*inv(U), which is what 
    !    we want.
    !
    !       G = inv(T)*inv(C)*inv(U)
    !         = inv(U*C*T)
    !         = inv(U*(inv(U)*inv(T)+D )*T)
    !         = inv(I+U*D*T)
    !         = inv(I+A)
    !         = inv((I+B_{i-1}...B_1*B_l...B_{i}))
    !
    ! 4. The sign of det(G) is determined by the LU decomposition.
    !
    ! Arguments
    ! =========
    !

#ifdef _SXX

    integer,  intent(in)       :: il            ! starting slice
    integer, intent(in)        :: n             ! dim of G
    real(wp), intent(inout)    :: sgn           ! sign(det(G)) 
    real(wp), intent(in)       :: V(:,:)        ! HSF
    real(wp), intent(inout)    :: G(:,:)        ! Green's function
    type(SeqB), intent(inout)  :: SB            
    integer, intent(inout)     :: pvt2(:)       ! working space 
    logical, intent(in)        :: compDet
    real(wp), intent(inout)    :: det           ! det(G)


    ! ... local scalars    ...
    integer :: info           ! parameters for lapack's sub
    integer :: i, j           ! iterator

    ! ... local variables for dgesvxx ...
    !
    ! Control parameters (use default values)
    integer, parameter ::  nparam = 0
    real(wp) :: param(1)
    ! Error estimates
    integer, parameter :: nerrb = 2
    real(wp) :: rcond, rpwgr
    real(wp) :: berr(n), errbn(n, nerrb), errbc(nerrb)
    ! info on equilibration
    character(len=1) :: eqt, eqw
    real(wp) :: Rw(n), Cw(n)
    ! Additional working space
    real(wp) :: work(4*n)
    integer  :: iwork(n)

    !Alias to SB's working space
    real(wp), pointer :: U(:,:), T(:,:), D(:), Rt(:), Ct(:)
    real(wp), pointer :: W1(:,:), W2(:,:), W3(:,:)
    integer,  pointer :: pvt1(:)

    ! ... Executable ...
    U    => SB%U
    D    => SB%D
    T    => SB%T
    W1   => SB%W1
    W2   => SB%W2
    W3   => SB%W3
    pvt1 => SB%piv
    Ct   => SB%tau
    Rt   => SB%rw

    sgn   = ONE
    det   = ZERO
    info  = 0

    !Compute UDT decomposition of B_{il}...B_1*B_L...B_{il+1}
    call DQMC_SeqMultB(il, il+1, SB, V)

    !Equilibrate T to improve conditioning
    call lapack_dggequb(n, n, T, n, Rt, Ct, work(1), work(2), work(3), info)
    call lapack_dlaqge(n, n, T, n, Rt, Ct, work(1), work(2), work(3), eqt)
    call blas_dcopy(n*n, T, 1, W3, 1)

    !Save U in W2 and scale its rows
    call blas_dcopy(n*n, U, 1, W2, 1)
    if(eqt=='C'.or.eqt=='B') call DQMC_ScaleRow(n, U, Ct)
    call blas_dcopy(n*n, U, 1, W1, 1)

    !LU factorize and solve
    call lapack_dgetrf(n, n, W3, n, pvt1, info)
    call lapack_dgetrs('T', n, n, W3, n, pvt1, W1, n, info)
    call lapack_dgerfs('T', n, n, T, n, W3, n, pvt1, U, n, W1, n, errbn(1:n,1), &
     &  berr, work, iwork, info)

    !Rescale so we have the solution to original problem.
    if(eqt=='R'.or.eqt=='B') call DQMC_ScaleRow(n, W1, Rt)

    !!Solve trans(T) W1 = U so that trans(W1) = inv(U) * inv(T)
    !call dgesvxx('E','T', n, n, T, n, W3, n, pvt1, eqt, Rt, Ct, &
    ! & U, n, W1, n, rcond, rpwgr, berr, nerrb, errbn, errbc,    &
    ! & nparam, param, work, iwork, info)

    ! determine the sign of det(T)
    call DQMC_DetSgn(n, W3, Rt, Ct, eqt, pvt1, sgn)
    !and log(abs(det(T))) 
    if (compDet) call DQMC_DetLog(n, W3, Rt, Ct, eqt, det)

    !Sum matrix of singular values to W1: W1 = W1 + D
    do i = 1,n
       W1(i,i) = W1(i,i) + D(i)
    end do

    !Restore U to its original self (may have been modified by equilibration)
    call blas_dcopy(n*n, W2, 1, U, 1)
    !and invert it. Store inverse in G.
    call DQMC_Trans(n, G, U)

    ! Solve trans(W1) Y = inv(U) so that Y = inv(inv(U)*inv(T) + D) inv(U)
    call lapack_dgesvxx('E','T', n, n, W1, n, W2, n, pvt2, eqw, Rw, Cw, &
     & G, n, U, n, rcond, rpwgr, berr, nerrb, errbn, errbc,      &
     & nparam, param, work, iwork, info)

    ! determine the sign of det(W1)
    call DQMC_DetSgn(n, W2, Rw, Cw, eqw, pvt2, sgn)
    !and log(abs(det(W1))) 
    if (compDet) call DQMC_DetLog(n, W2, Rw, Cw, eqw, det)

    !LU factorize and solve
    if(eqt=='R'.or.eqt=='B') call DQMC_ScaleRow(n, U, Rt)
    call blas_dcopy(n*n, U, 1, G, 1)
    call lapack_dgetrs('N', n, n, W3, n, pvt1, G, n, info)
    call lapack_dgerfs('T', n, n, T, n, W3, n, pvt1, U, n, G, n, errbn(1:n,1), &
     &  berr, work, iwork, info)

    !Rescale so we have the solution to original problem.
    if(eqt=='C'.or.eqt=='B') call DQMC_ScaleRow(n, G, Ct)

    ! Solve T G = Y so that G = inv(T) * inv(inv(U)*inv(T) + D) * inv(U)
    !call dgesvxx('F','N', n, n, T, n, W3, n, pvt1, eqt, Rt, Ct, &
    ! & U, n, G, n, rcond, rpwgr, berr, nerrb, errbn, errbc,     &
    ! & nparam, param, work, iwork, info)

    !Get determinant of the inverse
    if (compDet) det = -det

contains

    subroutine DQMC_DetSgn(n, A, R, C, equed, pvt, detsgn)

      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      real(wp), intent(in)    :: R(n), C(n)
      integer, intent(in)     :: pvt(n)
      character(len=1), intent(in) :: equed
      real(wp), intent(inout) :: detsgn
      
      ! ... local scalar 
      integer :: i

      do i = 1, n
         if (pvt(i).ne. i) then
            detsgn = -detsgn
         end if
         if (A(i,i) .lt. ZERO) then
            detsgn = -detsgn
         end if
      end do
      
      if(equed=='B'.or.equed=='C')then
        do i = 1, n
           if (C(i)<0.d0) detsgn = -detsgn
        enddo
      endif
    
      if(equed=='B'.or.equed=='R')then
        do i = 1, n
           if (R(i)<0.d0) sgn = -sgn
        enddo
      endif

    end subroutine DQMC_DetSgn


    subroutine DQMC_DetLog(n, A, R, C, equed, detlog)

      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      real(wp), intent(in)    :: R(n), C(n)
      character(len=1), intent(in) :: equed
      real(wp), intent(inout) :: detlog
      
      ! ... local scalar 
      integer :: i
      
      do i = 1, n
         detlog = detlog + log(abs(A(i,i)))
      end do

      if(equed=='B'.or.equed=='C')then
         do i = 1, n
            detlog = detlog + log(abs(C(i)))
         enddo
      endif
    
      if(equed=='B'.or.equed=='R')then
         do i = 1, n
            detlog = detlog + log(abs(R(i)))
         enddo
      endif

    end subroutine DQMC_DetLog

#else

    integer,  intent(in)       :: il            ! starting slice
    integer, intent(in)        :: n             ! dim of G
    real(wp), intent(inout)    :: sgn           ! sign(det(G)) 
    real(wp), intent(in)       :: V(:,:)        ! HSF
    real(wp), intent(inout)    :: G(:,:)        ! Green's function
    type(SeqB), intent(inout)  :: SB            
    integer, intent(inout)     :: pvt2(:)       ! working space 
    logical, intent(in)        :: compDet
    real(wp), intent(inout)    :: det           ! det(G)


    ! ... local scalars    ...
    integer :: info           ! parameters for lapack's sub
    integer :: i, j           ! iterator

    !Alias to SB's working space
    real(wp), pointer :: U(:,:), T(:,:), D(:), D2(:) 
    real(wp), pointer :: W1(:,:)
    integer,  pointer :: pvt1(:)

    ! ... Executable ...
    U   => SB%U
    D   => SB%D
    T   => SB%T
    W1  => SB%W1
    pvt1 => SB%piv
    D2  => SB%tau

    !
    ! STEP 0. Initialization
    ! ======================
    !! Setup index for B_{i}
    sgn = ONE
    info  = 0
    !
    ! STEP 1. Cmpute A = B_{il}...B_1*B_l...B_{il+1}.
    !         and its UDT decomposition, A = UDT
    ! =============================================
    call DQMC_SeqMultB(il, il+1, SB, V)

    ! STEP 2. Compute C = inv(U)*inv(T) + D.
    ! ======================================

    !! Compute the LU decomposition of T first.
    call lapack_dgetrf(n, n, T, n, pvt1, info)
    if (info .ne. 0) then
       T(info, info) = epsilon(ONE)
       call DQMC_Warning("T is singular, dgetrf in dqmc_ComputeG.", info)
    end if

    ! determine the sign of det(T)
    call DQMC_DetSgn(n, T, pvt1, sgn)
    if (compDet) then
       det   = ZERO
       call DQMC_DetLog(n, T, det)
    end if

    !! Solve T'W = U. W = transpose(inv(U)inv(T))
    ! W1 = U
    call blas_dcopy(n*n, U, 1, W1, 1)
    call lapack_dgetrs('T', n, n, T, n, pvt1, W1, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs in dqmc_getgp.", info)
    end if

    ! STEP 2.5 Here we need to decompose D = D2/D1=D1\D2.
    !       if   abs(D(i,i))>1,
    !            D1(i,i) = 1/abs(D(i,i))
    !            D2(i,i) = sign(D(i,i))
    !       else
    !            D1(i,i) = 1
    !            D2(i,i) = D(i,i) 
    ! ======================================

    !! Compute W*D1+D2  
    !!   = transpose(inv(U)inv(T))*D1+D2
    !!   = transpose(D1*inv(U)inv(T) + D2)
  
    do i = 1,n
       if (abs(D(i)) .gt. ONE) then
          if (D(i) .gt. ZERO) then
             D2(i) = ONE
          else
             D2(i) = -ONE
          end if
          D(i)  = ONE/abs(D(i))
       else
          D2(i) = D(i)
          D(i)  = ONE
       end if

      ! W' = W'*D1 
       do j = 1, n
          W1(j,i) = W1(j,i)*D(i)
       end do
       ! W' = W' + D2 = W'*D1+D2
       W1(i,i) = W1(i,i) + D2(i)
    end do
    
    !   
    ! STEP 3. Compute G = inv(T)*D1*inv(C)*inv(U)
    ! ===========================================   
    
    !! Compute the LU decomposition of W1.
    call lapack_dgetrf(n, n, W1, n, pvt2, info)

    if (info .ne. 0) then
       call DQMC_Warning("matrix D+inv(U)inv(T) is singular. &
            & dgetrf(2) in dqmc_computeG.", info)
       W1(info, info) = epsilon(ONE)
    end if

    ! determine the sign of det(W)
    call DQMC_DetSgn(n, W1, pvt2, sgn)
    if (compDet) then
       call DQMC_DetLog(n, W1, det)
    end if

    !! G = inv(U), which is transpose of U.
    call DQMC_Trans(n, G, U)

    !! G = D1*G = D1*inv(U)
    call DQMC_ScaleRow(n, G, D)
    U = G

    !! Solve W'G = inv(U), G = inv(W')*G
    call lapack_dgetrs('T', n, n, W1, n, pvt2, G, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs(2) in dqmc_getgp.", info)
    end if

    !! Solve TG = inv(W')*inv(U), G = inv(T)inv(W')*D1*inv(U)
    call lapack_dgetrs('N', n, n, T, n, pvt1, G, n, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrs(3) in dqmc_getgp.", info)
    end if

    !
    ! STEP 4. Compute the sign of det(G)
    ! ==================================
    
    !! Compute the LU decomposition of U first.
    call lapack_dgetrf(n, n, U, n, pvt2, info)
    if (info .ne. 0) then
       call DQMC_Error("Error: dgetrf(4) in dqmc_getgp.", info)
    end if
    call DQMC_DetSgn(n, U, pvt2, sgn)

    ! det = det(D1)/(det(W)*det(T))
    if (compDet) then
       det = -det
       do i = 1, n
          det = det + log(abs(D(i)))
       end do
    end if
    
  contains

    !--------------------------------------------------------!
    
    subroutine DQMC_DetSgn(n, A, pvt, sgn)
      !
      ! Purpose
      ! =======
      !    This subroutine computes sign(DET(A))
      !    where A is the content after dgetrf;
      !    pvt is the vipoting vector.
      ! 
      ! Pre-assumption
      ! ==============
      !    A contains the content of dgetrf
      !
      ! Argument
      ! ========
      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      integer, intent(in)     :: pvt(n)
      real(wp), intent(inout) :: sgn
      
      ! ... local scalar 
      integer :: i
      
      !! decide the sgn of det(Q)
      do i = 1, n
         if (pvt(i).ne. i) then
            sgn = -sgn
         end if

         !! ?? stable ??
         if (A(i,i) .lt. ZERO) then
            sgn = -sgn
         end if
      end do
      
    end subroutine DQMC_DetSgn

    !--------------------------------------------------------!
    
    subroutine DQMC_DetLog(n, A, det)
      !
      ! Purpose
      ! =======
      !    This subroutine computes log(abs(DET(A)))
      !    where A is the content after dgetrf;
      !    pvt is the vipoting vector.
      ! 
      ! Pre-assumption
      ! ==============
      !    A contains the content of dgetrf
      !
      ! Argument
      ! ========
      integer, intent(in)     :: n
      real(wp), intent(in)    :: A(n,n)
      real(wp), intent(inout) :: det
      
      ! ... local scalar 
      integer :: i
      
      !! decide the sgn of det(Q)
      do i = 1, n
         det = det + log(abs(A(i,i)))
      end do
      
    end subroutine DQMC_DetLog

#endif

  end subroutine DQMC_ComputeG

  !-----------------------------------------------------------------------!
  subroutine DQMC_UpdateWraps(G)
    !
    ! Purpose
    ! =======
    !    This subroutine updates G%nWraps depends on previous 
    !    execution results
    !
    ! Arguments
    ! =========
    !
    type(G_fun), intent(inout) :: G
        
    ! ... paremeters ...
    integer, parameter  :: DQMC_REDO_F      = 20
    integer, parameter  :: DQMC_NOREDO_F    = 500
    real(wp), parameter :: DQMC_REDO_RATE   = 0.20_wp
  

    ! ... local scalar ...
    real(wp) :: redorat     

    ! ... Executable ...
    if(G%redo .gt. DQMC_REDO_F) then
       redorat = dble(G%redo)/dble(G%redo+G%noredo)

       if(redorat.gt. G%errrate)then
          G%nwrap = G%nwrap - 1
          G%redo  = 0
          G%noredo = 1
       endif
    endif
    
    if(G%NoRedo .gt. DQMC_NOREDO_F) then
       redorat = dble(G%redo)/dble(G%redo+G%noredo)
       if(redorat .lt. DQMC_REDO_RATE*G%errrate) then
          if (G%nwrap .ge. G%lastwr) then
             G%nwrap = G%nwrap + 2
             G%redo = 0
             G%noredo = 1
             G%lastwr = G%nwrap
          end if
       endif
    endif

  end subroutine DQMC_UpdateWraps


  !-----------------------------------------------------------------------!
  ! The following three subroutines are used for delayed update.          !
  !-----------------------------------------------------------------------!
  
  function DQMC_Gfun_Getjj(n, j, blksz, G, U, W) result(gjj)
    !
    ! Purpose
    ! =======
    !    This function returns the (j,j)th element of G to gjj,
    !    in which 
    !
    !       G = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !    
    !    Therefore, 
    !
    !       gjj = G_1(j,j)+U(j,:)*V(j,:)'
    !
    ! Arguments
    ! =========
    !
    real(wp)                :: gjj
    integer , intent(in)    :: n, j, blksz
    real(wp), intent(in)    :: G(n,n)
    real(wp), intent(in)    :: U(n,n)
    real(wp), intent(in)    :: W(n,n)

    ! ... BLAS function ...
    real(wp), external :: blas_ddot

    ! ... Executable ...

    gjj = G(j,j)
    if (blkSz .gt. 0) then
       gjj = gjj + blas_ddot(blkSz, U(j,1), n, W(j,1), n)
    end if

  end function DQMC_Gfun_Getjj
  
 !--------------------------------------------------------------------------!
  
  subroutine DQMC_UpdateG(j, gamma, G)
    !
    ! Purpose
    ! =======
    !    This subroutine updates U, V, and D, which accumulate
    !    rank-1 updates of G.
    !
    !    Matrix G has the form
    !
    !       G_k = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !    
    !    The new matrix is 
    !
    !       G_{k+1} = G_{k} + gamma*xy'
    !    
    !    where 
    !    
    !       x = G_1(:,j) + UV'(:,j)
    !       y' = G_1(j,:) - e_j' + U(j,:)V'
    !
    !    The update of U, V, and D is as follows.
    !
    !       U = [U x], V = [V gamma*y]
    !
    !    If blkSz==nBlk, then apply the update
    ! Arguments
    ! =========
    !

    integer , intent(in)       :: j
    real(wp), intent(in)       :: gamma
    type(g_fun), intent(inout), target :: G

    ! ... local variables ...
    real(wp), pointer :: x(:), y(:)
    real(wp)  :: xx(G%n)
    integer   :: n, blksz

    ! ... Executable ...

    n     = G%n
    blksz = G%blkSz
    x => G%U(1:n, blkSz+1)
    y => G%W(1:n, blkSz+1)

    x = G%G(1:n,j)
    y = G%G(j,1:n)
    y(j) = y(j) - ONE

    if (blkSz .gt. 0) then
       ! if U, V are not empty, add their effects
       xx = G%W(j,1:n)
       call blas_dgemv('N', n, blkSz, ONE, G%U, n, xx, 1, ONE, x, 1)
       xx = G%U(j,1:n)
       call blas_dgemv('N', n, blkSz, ONE, G%W, n, xx, 1, ONE, y, 1)
    end if
    call blas_dscal(n, gamma, y, 1)
    
    G%blkSz = G%blkSz + 1
 
    ! apply the update when necessary
    call DQMC_ApplyUpdate(G, .false.)
    !call DQMC_ApplyUpdate(G, .true.)
        
  end subroutine DQMC_UpdateG

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_ApplyUpdate(G, forced)
    !
    ! Purpose
    ! =======
    !    This subroutine updates G with U, V, and D.
    !    Matrix G has the form
    !
    !       G_k = G_1 + UV'
    !
    !    where U,V are n*BlkSz.
    !
    ! Arguments
    ! =========
    !
    type(g_fun), intent(inout) :: G
    logical, intent(in)        :: forced

    ! ... local variables ...
    integer :: n

    ! ... Executable ...

    if (forced .or. G%blkSz .eq. G%nBlk) then
       n = G%n
       ! apply the update when necessary
       call blas_dgemm('N','T', n, n, G%blkSz, ONE, G%U, n, G%W, n, ONE, G%G, n)
       
       ! reset the block size
       G%blkSz = 0
    end if

  end subroutine DQMC_ApplyUpdate

  !--------------------------------------------------------------------------!

end module DQMC_GFun
