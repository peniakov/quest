module DQMC_CKB
#include "dqmc_include.h"

  use DQMC_UTIL
  use DQMC_WSPACE

  implicit none 

  ! 
  ! This module defines the type and subroutines for propagator B.
  ! B is defined by the checkerboard method. See [1]
  !
  ! For a 2-dimensional model
  !
  ! B = exp(K) = exp(K_{4})exp(K_{3})exp(K_{2})exp(K_{1})
  !            = B_4B_3B_2B_1
  !            
  !    [1] W. Hanke and Yu.V. Kopaev, "Stable Numerical Simulation 
  !        of Models of Interating Electrons in Condensed-Matter Physics."
  !        Chapter 4. Elsevier Science Pub. 
  !
  !  Data Type
  !  =========
  !
  type matB
     integer  :: n                      ! dim of B
     integer  :: m                      ! number of neighbors of lattice
     integer, pointer :: A(:,:,:)            ! Adjacency info
     real(wp), pointer :: factor(:)      ! temp variable in multiplication
     real(wp), pointer :: sinht(:), cosht(:), f1(:), f2(:)            ! parameters for checkerboard method
     character(12) :: name
  end type matB

contains


  !-------------------------------------------------------------------------!

  subroutine DQMC_B_Init(n, B, WS, Adj, ckb, t, mu, dtau)
    !
    ! Purpose
    ! =======
    !    This subroutine initiliazes the data type of Green function.
    !
    ! Pre-assumption
    ! ==============
    !    This module only used for one band Hubbard's model.
    !
    ! Arguments
    ! =========
    !

    use DQMC_STRUCT

    integer, intent(in)       :: n                     ! Number of sites
    type(MatB), intent(out)   :: B                     ! MatB
    type(CCS), intent(in)     :: adj                   ! adjacent info
    type(CCS), intent(in)     :: ckb
    real(wp), intent(in)      :: t(*)              ! model parameter
    real(wp), intent(in)      :: mu(n), dtau  
    type(WSpace), intent(in), target  :: WS                    ! shared working space

    ! ... local scalars    ...
    integer :: k, i, j, n_t
    real(wp), pointer :: dum(:,:)
    
    ! ... Executable ...

    dum   => WS%R1

    B%n    = n
    B%name = "Checkerboard"

    B%m   = maxval(ckb%A)
    n_t   = maxval(Adj%A)

    allocate(B%A(2,n,B%m))
    B%A = 0

    allocate(B%sinht(0:n_t), B%cosht(0:n_t))
    B%sinht(0)=0.d0
    B%cosht(0)=1.d0
    if(n_t .ge. 1)then
       B%sinht(1:n_t)=sinh(dtau*t(1:n_t))
       B%cosht(1:n_t)=cosh(dtau*t(1:n_t))
    endif

    if(sum((ckb%row-adj%row)**2)+sum((ckb%cstart-adj%cstart)**2)>0)then
       write(*,*)'ckb and adj do not conform. Stop.'
       stop
    endif

    !Initialize B%A so that, if nothing else is done, the hopping
    !part of B is the identity matrix
    do i = 1, n
      !The value of the r.h.s is unimportant (must be between 1 and n)
      !although this choice may be optimal memory-wise.
       B%A(1, i, 1:B%m) = i 
    enddo
    B%A(2, 1:n, 1:B%m) = 0

    !k is the column index
    k = 0
    do i = 1, ckb%nnz
       !When a new column start we increase the column index
       if(ckb%cstart(k+1) == i) k = k + 1
       !j is the row index
       j = ckb%row(i)
       B%A(1, j, ckb%A(i)) = k
       B%A(2, j, ckb%A(i)) = adj%A(i)
    enddo

    allocate(B%f1(n), B%f2(n))
    do i=1,n
       B%f1(i) = exp(dtau*mu(i))
       B%f2(i) = 1.d0/B%f1(i)
    enddo

    allocate(B%factor(n))

  end subroutine DQMC_B_Init

  !-------------------------------------------------------------------------!

  subroutine DQMC_B_Free(B)
    !
    ! Purpose
    ! =======
    !    This subroutine frees memory of B.
    !
    ! Arguments
    ! =========
    !
    type(MatB), intent(inout)  :: B  ! MatB

    ! ... Executable ...

    deallocate(B%factor,B%A)

  end subroutine DQMC_B_Free


  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    leftside.
    !
    !      B_i = V_i*B
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         ! 

    ! ... Local variables ...
    integer  :: i, j, k
    integer,  pointer :: A(:,:,:) 
    real(wp), pointer :: P(:,:), Q(:,:), R(:,:)

    ! ... executable ...
    A => B%A

    ! Combin V_i with the coefficient of B
    B%factor(1:n) = B%f1(1:n) * V_i(1:n)
    
    ! Multiply B from left-hand-side

    do j = 1, n ! do it column by column         

       P => C
       Q => M
       R => P
       do k = 1, B%m
          do i = 1, n
             P(i,j) = B%cosht(A(2,i,k)) * Q(i,j) + B%sinht(A(2,i,k)) * Q(A(1,i,k),j)
          end do
          P => Q
          Q => R
          R => P
       enddo
       
       ! Multiply V (rescale rows of M)
       M(1:n,j) = B%factor(1:n) * Q(1:n,j)

    end do

  end subroutine DQMC_MultB_Left

  !----------------------------------------------------------------------!
  
  subroutine DQMC_MultB_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply B from
    !    righthand side.
    !
    !      B_i = V_i*B
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !    3. Flag lr is either 'l' or 'r'.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !    

    ! ... Local variables ...
    integer  :: i, j, k
    integer,  pointer :: A(:,:,:) 
    real(wp), pointer :: P(:,:), Q(:,:), R(:,:)

    ! ... executable ...
    A => B%A

    ! Combin V_i with the coefficient of B
    B%factor(1:n) = B%f1(1:n) * V_i(1:n)

    call DQMC_ScaleCol(n, M, B%factor)

    ! Transpose M first
    call DQMC_Trans(n, C, M)
    
    ! Multiply B from left-hand-side

    do j = 1, n ! do it column by column         
     
       P => M
       Q => C
       R => P
       do k = B%m, 1, -1
          do i = 1, n
             P(i,j) = B%cosht(A(2,i,k)) * Q(i,j) + B%sinht(A(2,i,k)) * Q(A(1,i,k),j)
          end do
          P => Q
          Q => R
          R => P
       enddo

    end do
    
    ! transpose back
    call DQMC_Trans(n, P, Q)

    if(associated(P, target=C)) M = C

  end subroutine DQMC_MultB_Right

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Left(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
    !    leftside.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k
    integer,  pointer :: A(:,:,:) 
    real(wp), pointer :: P(:,:), Q(:,:), R(:,:)

    ! ... executable ...

    A => B%A

    ! Combin V_i with the coefficient of B
    B%factor(1:n) = B%f2(1:n) / V_i(1:n)

    do j = 1, n ! do it column by column         

       M(1:n,j) = B%factor(1:n) * M(1:n,j)
       
       P => C
       Q => M
       R => P
       do k = B%m, 1, -1
          do i = 1, n
             P(i,j) = B%cosht(A(2,i,k)) * Q(i,j) - B%sinht(A(2,i,k)) * Q(A(1,i,k),j)
          end do
          P => Q
          Q => R
          R => P
       enddo
       
    end do

    if(associated(Q, target=C)) M = C
    
  end subroutine DQMC_MultBi_Left

  !--------------------------------------------------------------------------!
  
  subroutine DQMC_MultBi_Right(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine uses checkerboard method to multiply inv(B) from
    !    right hand side.
    !
    ! Pre-assumption
    ! ==============
    !    1. Matrix M, B and C are square with the same order n.
    !    2. Vector V_i is of length n. 
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)              :: n              ! The order of matrix M  
    type(MatB), intent(inout)        :: B              !
    real(wp), intent(in)             :: V_i(n)         !
    real(wp), intent(inout), target  :: M(n,n)         !  
    real(wp), intent(inout), target  :: C(n,n)         !
    
    ! ... Local variables ...
    integer  :: i, j, k
    integer,  pointer :: A(:,:,:) 
    real(wp), pointer :: P(:,:), Q(:,:), R(:,:)

    ! ... executable ...

    A => B%A

    ! Combin V_i with the coefficient of B
    B%factor(1:n) = B%f2(1:n) / V_i(1:n)

    call DQMC_Trans(n, C, M)

    do j = 1, n ! do it column by column         

       P => M
       Q => C
       R => P
       do k = 1, B%m
          do i = 1, n
             P(i,j) = B%cosht(A(2,i,k)) * Q(i,j) - B%sinht(A(2,i,k)) * Q(A(1,i,k),j)
          end do
          P => Q
          Q => R
          R => P
       enddo

    end do

    ! transpose back
    call DQMC_Trans(n, P, Q)

    do i = 1, n
       M(1:n,i) = P(1:n,i) * B%factor(i)
    enddo

  end subroutine DQMC_MultBi_Right

  !----------------------------------------------------------------------!

  subroutine DQMC_MultrtB0_Left(n, M, B, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix sqrt(B0) from
    !    leftside.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(inout)  :: C(n,n)         ! working space   

    integer :: i

    ! Dummy executable to avoid warning
    C = M
    i = B%n

  end subroutine DQMC_MultrtB0_Left

  !----------------------------------------------------------------------!

  subroutine DQMC_MultrtB0i_Right(n, M, B, C)
    !
    ! Purpose
    ! =======
    !    This subroutine mulitplies matrix M by the matrix sqrt(B0_i) from
    !    rightside.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(in)   :: B              !
    real(wp), intent(inout)  :: C(n,n)         ! working space  

    integer :: i

    ! Dummy executable to avoid warning
    C = M
    i = B%n

  end subroutine DQMC_MultrtB0i_Right

  !-------------------------------------------------------------------------!

  subroutine DQMC_GetB(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = V_iB
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(inout):: B              !
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultB_Left(n, M, B, V_i, C)

  end subroutine DQMC_GetB

  !-----------------------------------------------------------------------!

  subroutine DQMC_GetBi(n, M, B, V_i, C)
    !
    ! Purpose
    ! =======
    !    This subroutine returns M = inv(B)inv(V_i)
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)      :: n              ! The order of matrix M  
    real(wp), intent(inout)  :: M(n,n)         !  
    type(MatB), intent(inout):: B              ! MatB 
    real(wp), intent(in)     :: V_i(n)         !
    real(wp), intent(inout)  :: C(n,n)         !
    ! ... Executable ...

    call DQMC_Eye(n, M)
    call DQMC_MultBi_Left(n, M, B, V_i, C)

  end subroutine DQMC_GetBi

  !-----------------------------------------------------------------------!
  
end module DQMC_CKB
