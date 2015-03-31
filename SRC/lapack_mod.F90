module lapack_mod

  use blas_mod

  integer, parameter :: DLARNV_UNI_0_1  = 1
  integer, parameter :: DLARNV_UNI_N1_1 = 2
  integer, parameter :: DLARNV_NORMAL   = 3

#define DB kind(1.0d0)

  interface

     ! INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
     ! choose problem-dependent parameters for the local environment
     integer function ilaenv(ispec, name, opts, n1, n2, n3, n4)
       character*(*) :: name, opts
       integer   :: ispec, n1, n2, n3, n4
     end function ilaenv

     ! random number generator
     subroutine dlarnv(idist, iseed, n, x)
       integer   :: idist, n, iseed(4)
       real(DB)  :: x(n)
     end subroutine dlarnv
  
  end interface

#if DQMC_PROFILE

#define LAPACK_BEGIN if (lapack_profile) call system_clock(lapack_t1)
#define LAPACK_END(c, t) if (lapack_profile) call lapack_end(c, t)

    real(DB) :: lapack_dsyev_time = 0, lapack_dgesv_time = 0, lapack_dgeqp3_time = 0
    real(DB) :: lapack_dorgqr_time = 0, lapack_dormqr_time = 0, lapack_dgetrf_time = 0
    real(DB) :: lapack_dgetri_time = 0, lapack_dgetrs_time = 0
    real(DB) :: lapack_dgejsv_time = 0, lapack_dgerfsx_time = 0
    integer  :: lapack_dsyev_count = 0, lapack_dgesv_count = 0, lapack_dgeqp3_count = 0
    integer  :: lapack_dorgqr_count = 0, lapack_dormqr_count = 0, lapack_dgetrf_count = 0
    integer  :: lapack_dgetri_count = 0, lapack_dgetrs_count = 0
    integer  :: lapack_dgejsv_count = 0, lapack_dgerfsx_count = 0
    integer  :: lapack_t1, lapack_t2, lapack_profile = 0, lapack_rate

  contains

    subroutine lapack_end(c, t)
      integer :: c
      real(DB) :: t
      call system_clock(lapack_t2, lapack_rate)
      c = c + 1
      t = t + (lapack_t2 - lapack_t1) /  REAL(lapack_rate)
    end subroutine lapack_end

     subroutine lapack_print()
       write(*,*) "DSYEV         ", lapack_dsyev_count, lapack_dsyev_time
       write(*,*) "DGESV         ", lapack_dgesv_count, lapack_dgesv_time
       write(*,*) "DGEQP3        ", lapack_dgeqp3_count, lapack_dgeqp3_time
       write(*,*) "DORGQR        ", lapack_dorgqr_count, lapack_dorgqr_time
       write(*,*) "DORMQR        ", lapack_dormqr_count, lapack_dormqr_time
       write(*,*) "DGETRF        ", lapack_dgetrf_count, lapack_dgetrf_time
       write(*,*) "DGETRI        ", lapack_dgetri_count, lapack_dgetri_time
       write(*,*) "DGETRS        ", lapack_dgetrs_count, lapack_dgetrs_time
       write(*,*) "DGEJSV        ", lapack_dgejsv_count, lapack_dgejsv_time
       write(*,*) "DGERFSX       ", lapack_dgerfsx_count, lapack_dgerfsx_time
     end subroutine lapack_print

#else 

#define LAPACK_BEGIN 
#define LAPACK_END(c, t) 

  contains

#endif
       
     ! compute eigenvalues/eigenvectors of symmetric A 
     subroutine lapack_dsyev(jobZ, uplo, n, A, lda, W, Work, lwork, info)
       character :: jobZ, uplo
       integer   :: n, lda, lwork, info
       real(DB)  :: A(lda,*), W(*), Work(*)
       LAPACK_BEGIN
       call dsyev(jobZ, uplo, n, A, lda, W, Work, lwork, info)
       LAPACK_END(lapack_dsyev_count, lapack_dsyev_time) 
     end subroutine lapack_dsyev

     ! Linear system solver for general matrix
     subroutine lapack_dgesv(n, rhs, A, lda, pvt, B, ldb, info)
       integer   :: n, rhs, lda, ldb, info, pvt(*)
       real(DB)  :: A(lda,*), B(ldb, *)
       LAPACK_BEGIN
       call dgesv(n, rhs, A, lda, pvt, B, ldb, info)
       LAPACK_END(lapack_dgesv_count, lapack_dgesv_time) 
     end subroutine lapack_dgesv

     ! QR-factorization with column pivoting
     subroutine lapack_dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
       integer   :: m, n, lda, lwork, info, jpvt(*)
       real(DB)  :: A(lda,*), tau(*), work(*)
       LAPACK_BEGIN
       call dqmc_dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
       LAPACK_END(lapack_dgeqp3_count, lapack_dgeqp3_time) 
     end subroutine lapack_dgeqp3

     ! generate Q-factor
     subroutine lapack_dorgqr(m, n, k, A, lda, tau, work, lwork, info)
       integer   :: m, n, k, lda, lwork, info
       real(DB)  :: A(lda,*), work(*), tau(*)
       LAPACK_BEGIN
       call dorgqr(m, n, k, A, lda, tau, work, lwork, info)
       LAPACK_END(lapack_dorgqr_count, lapack_dorgqr_time) 
     end subroutine lapack_dorgqr
     
     ! multiply Q-factor to a matrix
     subroutine lapack_dormqr(side, trans, m, n, k, A, lda, tau, c, ldc, &
          work, lwork, info)
       character :: side, trans
       integer   :: info, k, lda, ldc, lwork, m, n
       real(DB)  :: A(lda,*), C(ldc,*), tau(*), work(*)
       LAPACK_BEGIN
       call dormqr(side, trans, m, n, k, A, lda, tau, c, ldc, &
          work, lwork, info)
       LAPACK_END(lapack_dormqr_count, lapack_dormqr_time) 
     end subroutine lapack_dormqr

     ! LU-decomposition
     subroutine lapack_dgetrf(m, n, A, lda, ipiv, info)
       integer   :: info, lda, m, n, ipiv(*)
       real(DB)  :: A(lda,*)
       LAPACK_BEGIN
       call dgetrf(m, n, A, lda, ipiv, info)
       LAPACK_END(lapack_dgetrf_count, lapack_dgetrf_time) 
     end subroutine lapack_dgetrf

     ! Compute the inverse
     subroutine lapack_dgetri(n, A, lda, ipiv, work, lwork, info)
       integer   :: info, lda, n, ipiv(*), lwork
       real(DB)  :: A(lda,*), work(*)
       LAPACK_BEGIN
       call dgetri(n, A, lda, ipiv, work, lwork, info)
       LAPACK_END(lapack_dgetri_count, lapack_dgetri_time) 
     end subroutine lapack_dgetri

     ! solves a system of linear equations
     subroutine lapack_dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
       character :: trans
       integer   :: info, lda, ldb, n, nrhs, ipiv(*)
       real(DB)  :: A(lda,*), B(ldb,*)
       LAPACK_BEGIN
       call dgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, info)
       LAPACK_END(lapack_dgetrs_count, lapack_dgetrs_time) 
    end subroutine lapack_dgetrs

     subroutine lapack_dgejsv(joba, jobu, jobv, jobr, jobt, jobp,     &
                       m, n, a, lda, sva, u, ldu, v, ldv,      &
		       work, lwork, iwork, info)
       integer   :: info, lda, ldu, ldv, lwork, m, n
       real(DB)  :: a(lda,*), sva(n), u(ldu,*), v(ldv,*), work(lwork)
       integer   :: iwork(*)
       character :: joba, jobp, jobr, jobt, jobu, jobv
       LAPACK_BEGIN
       call dgejsv(joba, jobu, jobv, jobr, jobt, jobp,     &
                       m, n, a, lda, sva, u, ldu, v, ldv,      &
		       work, lwork, iwork, info)
       LAPACK_END(lapack_dgejsv_count, lapack_dgejsv_time) 
     end subroutine lapack_dgejsv

     subroutine lapack_dgerfsx(trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, &
		     r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, &
		     err_bnds_norm, err_bnds_comp, nparams, params, &
		     work, iwork, info) 
       character :: trans, equed
       integer   :: info, lda, ldaf, ldb, ldx, n, nrhs, nparams, n_err_bnds
       real(DB)  :: rcond
       integer   :: ipiv(*), iwork(* )
       real(DB)  :: a(lda,*), af(ldaf,*), b(ldb,*), x(ldx,*), work(*)
       real(DB)  :: r(*), c(*), params(*), berr(*), &
                    err_bnds_norm(nrhs,*), err_bnds_comp(nrhs,*)
       LAPACK_BEGIN
       call dgerfsx( trans, equed, n, nrhs, a, lda, af, ldaf, ipiv, &
		     r, c, b, ldb, x, ldx, rcond, berr, n_err_bnds, &
		     err_bnds_norm, err_bnds_comp, nparams, params, &
		     work, iwork, info )  
       LAPACK_END(lapack_dgerfsx_count, lapack_dgerfsx_time) 
     end subroutine lapack_dgerfsx
       
end module lapack_mod
