   module FullG
  
      type Gkw_t
         integer :: na   ! number of basis functions
         integer :: nk   ! number of k-points (or translations)
         integer :: nw   ! number of frequencies (or time slices)
         real*8  :: minw

         complex*16, pointer :: wgt_rk(:,:)  ! Real-to-k space Fourier coefficients
         complex*16, pointer :: wgt_tw(:,:)  ! Time-to-frequency Fourier coefficients

         complex*16, pointer :: Gkw_up(:,:)  ! Green's functions in k-space and 
         complex*16, pointer :: Gkw_dn(:,:)  ! frequency for up and down
      end type

      type wdm_t
         integer :: nclass
         integer :: nclass2

         complex*16, pointer :: Gfun(:,:)      ! Average Green's function
         complex*16, pointer :: G2zero(:,:,:,:)! Average two-ptcl Green's function
         complex*16, pointer :: G2pi(:,:,:,:)  ! Average two-ptcl Green's function
      end type

   end module

   program dqmc_test

     use FullG
     use DQMC_Cfg
     use DQMC_Geom_Wrap
     use DQMC_Hubbard
     use DQMC_gtau
     use DQMC_kbonds
     implicit none

     integer, parameter :: nmom = 2
     real(wp) :: ksum(nmom, rdim)
   
     real                :: t1, t2
     type(Config)        :: cfg
     type(Hubbard)       :: Hub
     type(GeomWrap)      :: Gwrap
     type(Gtau)          :: tau
     type(Gkw_t)         :: Gkw
     type(kbonds_t)      :: kbonds
     character(len=slen) :: gfile
     logical             :: tformat
     integer             :: OPT, na, nt, nkt, nkg, ntausk
     character(len=30), allocatable   :: clabelt(:), clabelg(:)

     integer :: iao, iat, iko, io, ikt, it
     integer :: jao, jat, jko, jo, jkt, jt
     integer :: imom, ib, jb, i

     !Hardwire the number of total momentum channel. (0,0) and (pi,pi)
     ksum = 0.d0
     ksum(2,1) = acos(-1.d0)
     ksum(2,2) = ksum(2,1)
   
     call cpu_time(t1)  
     OPT=6
     
     ! Read input file
     call DQMC_Read_Config(cfg, STDIN)
      
     ! Initialize geometry
     call CFG_Get(cfg, "gfile", gfile)
     call DQMC_Geom_Read_Def(Hub%S, gfile, tformat)
     if(.not.tformat)then
        call DQMC_Geom_Fill(Gwrap, gfile, cfg)
        call DQMC_Geom_Init(Gwrap, Hub%S, cfg)
        call DQMC_init_kmate(Gwrap%RecipLattice, nmom, ksum)
        call init_kbonds(Gwrap%symmOp, Gwrap%lattice, Gwrap%RecipLattice, kbonds)
        call construct_kbonds(Gwrap%RecipLattice, kbonds)
        call map_symm_kbond(kbonds)
        call construct_kbond_classes(kbonds)
     endif

     na=Gwrap%lattice%natom
     nkt=Gwrap%RecipLattice%nkpts
     do imom = 1, kbonds%nmomenta
        do i = 1, kBonds%nclass(imom)
           write(*,*)
           write(*,*)
           write(*,'(A,1x,i3,1x,i3)')'Class',imom,i
           do ib = 1, kbonds%nbonds
              do jb = 1, kbonds%nbonds
                 if (kbonds%myclass(ib,jb,imom) == i) then
                    iao = kbonds%bond_origin(ib, imom)
                    iat = kbonds%bond_target(ib, imom)
                    iko = mod(iao -1, nkt) + 1
                    io  = (iao - 1)/ nkt
                    ikt = mod(iat -1, nkt) + 1
                    it  = (iat - 1)/ nkt
                    jao = kbonds%bond_origin(jb, imom)
                    jat = kbonds%bond_target(jb, imom)
                    jko = mod(jao -1, nkt) + 1
                    jo  = (jao - 1) / nkt
                    jkt = mod(jat -1, nkt) + 1
                    jt  = (jat - 1) / nkt
                    write(*,'(2(i3,1x,3(f10.5)))') io, Gwrap%RecipLattice%klist(iko,1:3),&
                               it, Gwrap%RecipLattice%klist(ikt,1:3)
                    write(*,'(2(i3,1x,3(f10.5)))') jo, Gwrap%RecipLattice%klist(jko,1:3),&
                               jt, Gwrap%RecipLattice%klist(jkt,1:3)
                    write(*,*)
                 endif
              enddo
           enddo
        enddo
     enddo
   
     stop
     
     ! Initialize the rest data
     call DQMC_Hub_Config(Hub, cfg)

     ! Initialize Time dependent Green's function
     call DQMC_Gtau_Init(Hub%n, Hub%L, TAU_BOTH, Hub%SB%north, Hub%G_up%nWrap, tau, Hub%B, Hub%WS)

     !Initialize Gkw
     call Gkw_Init(Gkw, Gwrap%RecipLattice, Gwrap%Lattice, Hub%L, Hub%dtau)
   
     ! Execution MC loop
     call CFG_Get(cfg, "tausk", ntausk)
     !call DQMC_Hub_Run(Hub)!, tau, ntausk)
     call Hub_Run(Hub, tau, Gkw, ntausk)
   
     ! Print computed results
     call DQMC_Hub_OutputParam(Hub, OPT)
     write(OPT, FMT_DBLINE)
   
     call DQMC_Phy0_Print(Hub%P0, Hub%S, OPT)
     
     ! Fill the matrix of fourier weights
     call DQMC_Fill_FourierC(Gwrap%RecipLattice, Gwrap%Lattice)
     call DQMC_Fill_FourierC(Gwrap%GammaLattice, Gwrap%Lattice)
   
     write(*,*)'Filled FC'
   
     na=Gwrap%lattice%natom
   
     nkt=Gwrap%RecipLattice%nclass_k
     nkg=Gwrap%GammaLattice%nclass_k
   
     !Compute Fourier transform
     call DQMC_phy0_Get_FT(Hub%P0, Hub%S%D, Hub%S%gf_phase, Gwrap%RecipLattice%FourierC, &
          Gwrap%GammaLattice%FourierC, nkt, nkg, na, nt)
   
     write(*,*)'Computed FT'
   
     allocate(clabelt(na*nkt),clabelg(na*nkg))  
   
     !Print info on k-points and construct clabel
     call DQMC_Print_HeaderFT(Gwrap,OPT,clabelt,.true.)
     call DQMC_Print_HeaderFT(Gwrap,OPT,clabelg,.false.)
   
     !Print computed reciprocal space properties
     call DQMC_Phy0_Print_FT(Hub%P0, na, nkt, nkg, clabelt, clabelg, OPT)
   
     if(Hub%P2%compute)then
   
      if(Hub%P2%diagonalize)then
       
        call DQMC_Phy2_Get_Wave(Hub%P2, Hub%P0%G_fun(:,Hub%P0%avg), Hub%S)
        !call DQMC_Phy2_Get_Wave(Hub%P2,  Hub%S)
   
        !Get error for waves
        call DQMC_Phy2_GetErr_Symm(Hub%P2, Hub%P0%G_fun, Hub%S)
   
        !Analyze symmetry of pairing modes
        call DQMC_Phy2_WaveSymm(Hub%S,Hub%P2,Gwrap%SymmOp)
                               
        !Print Pairing info
        call dqmc_phy2_PrintSymm(Hub%S, Hub%P2, OPT)
      else
        call dqmc_phy2_print(Hub%P2, Hub%S%wlabel ,OPT)
      endif
   
     endif
     
     ! Clean up the used storage
     call DQMC_Hub_Free(Hub)
     call DQMC_Config_Free(cfg)
     
     call cpu_time(t2)
     write(STDOUT,*) "Running time:",  t2-t1, "(second)"
   
   end program dqmc_test

  !---------------------------------------------------------------------!

  subroutine Hub_Run(Hub, tau, Gkw, ntausk)
    use DQMC_Hubbard
    use DQMC_Gtau
    use FullG
    !
    ! Purpose
    ! =======
    !   This subroutine is the main subroutine for DQMC.
    !   There are four major wroks
    !
    !      1. Compute Green function.
    !      2. Perform warmup sweep.
    !      3. Perform actual sweep.
    !      4. Analyze the measurement. (see DQMC_Phy0)
    !
    ! Arguments
    ! =========
    !
    type(Hubbard), intent(inout) :: Hub    ! Hubbard model
    type(Gtau), intent(inout)    :: tau    ! Hubbard model
    type(Gkw_t), intent(inout)   :: Gkw
    integer, intent(in)          :: ntausk

    ! ... local scalar ...
    integer  :: i, j, k, nIter, nBin

    ! ... Executable ...

    ! Warmup sweep
    do i = 1, Hub%nWarm
       ! The second parameter means no measurement should be made.
       call DQMC_Hub_Sweep(Hub, NO_MEAS0)
       call DQMC_Hub_Sweep2(Hub, Hub%nTry)
    end do
 
    ! We divide all the measurement into nBin,
    ! each having nPass/nBin pass.
    nBin   = Hub%P0%nBin 
    nIter  = Hub%nPass/nBin/ntausk
    do i = 1, nBin
       do j = 1, nIter
          do k = 1, ntausk
            call DQMC_Hub_Sweep(Hub, Hub%nMeas)
            call DQMC_Hub_Sweep2(Hub, Hub%nTry)
          enddo
          !Compute full Green's function 
          call Get_Gkw(tau, Hub%G_up, Hub%G_dn, Gkw)
          !Compute frequency dependent properties
          !call Get_wdm(Gkw_up, Gkw_dn, wdm)
       end do

       ! Accumulate results for each bin
       call DQMC_Phy0_Avg(Hub%P0)
       if (Hub%meas2) then
          if(Hub%P2%diagonalize)then
            call DQMC_Phy2_Avg(Hub%P2, Hub%S)
          else
            call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
          endif
       end if
    end do

    ! Get average result
    call DQMC_Phy0_GetErr(Hub%P0)
    if (Hub%meas2) then
       call DQMC_Phy2_GetErr(Hub%P2)
    end if

  end subroutine Hub_Run
  
  !-------------------------------------------------------------------!

  subroutine Get_Gkw(tau, G_up, G_dn, Gkw)
     use DQMC_Gtau
     use DQMC_Gfun
     use FullG
     implicit none

     type(Gtau), intent(inout)        :: tau
     type(G_fun), intent(inout)       :: G_up, G_dn
     type(Gkw_t), intent(inout)       :: Gkw
  
     integer :: ii, jj, i, j, idir, n, m, L, idx, jdx, is, js, north, nt, na

     real*8 :: a, b, c

     complex*16, allocatable :: up_ij(:, :, :, :)   !  temporary G(t_i,t_j) storage
     complex*16, allocatable :: up_ji(:, :, :, :)   ! 
     complex*16, allocatable :: dn_ij(:, :, :, :)   !
     complex*16, allocatable :: dn_ji(:, :, :, :)   !

     !real*8, allocatable :: Aup(:, :)
     !real*8, allocatable :: Adn(:, :)
     complex*16, pointer :: Gptr(:, :)

     complex*16, pointer :: fwgt(:, :)
     complex*16, pointer :: tmp(:, :)

     integer :: it, jt
     complex*16, pointer :: tmp2(:, :, :)

     n = tau%n
     L = tau%L
     north = tau%SB1%north
     m = (north-1)/2
     nt = Gkw%nk
     na = Gkw%na
      
     allocate(up_ij(n, n, -m:m, -m:m)) 
     allocate(dn_ij(n, n, -m:m, -m:m)) 
     allocate(up_ji(n, n, -m:m, -m:m)) 
     allocate(dn_ji(n, n, -m:m, -m:m))

     !allocate(Aup(n*L, n*L))
     !allocate(Adn(n*L, n*L))

     allocate(tmp(nt, nt))
     fwgt => Gkw%wgt_rk

     !DEBUG
     !Compare to straight inversion
     !call DQMC_Gtau_Big(tau, Aup, Adn, G_up, G_dn)

     !loop over blocks north slices apart. keep jj .ge. ii
     do ii = m+1, L, north
        do jj = ii, L, north

           !Get Gtau from scratch
           call DQMC_GetGtau2(ii, jj, tau%upt0, tau%up0t, G_up%V, tau)
           call DQMC_GetGtau2(ii, jj, tau%dnt0, tau%dn0t, G_dn%V, tau)

           !DEBUG
           !do i=1,n
           !   do j = 1, n
           !      idx = (ii-1) * n + i
           !      jdx = (jj-1) * n + j
           !      a = Aup(idx, jdx)
           !      b = tau%upt0(i,j)
           !      if(abs(a-b)>1.d-12)then
           !        write(*,*) ii, jj, i, j
           !        write(*,*) a, b
           !      endif
           !      idx = (jj-1) * n + i
           !      jdx = (ii-1) * n + j
           !      a = Aup(idx, jdx)
           !      b = tau%up0t(i,j)
           !      if(abs(a-b)>1.d-12)then
           !        write(*,*) jj, ii, i, j
           !        write(*,*) a, b
           !      endif
           !   enddo
           !enddo

           !Store the block label in tau
           tau%ii = ii
           tau%ib = jj

           call save_gtau(0,0)
           
           !For G(tau,0) and G(0,tau) fill the surrouding box
           idir = 1
           do i = 1, m
             call change_gtau_time(idir, tau, G_up, G_dn)
             call save_gtau(i,0)
             !DEBUG
             !do is=1,n
             !   do js = 1, n
             !      idx = (ii-1+i) * n + is
             !      jdx = (jj-1) * n + js
             !      a = Aup(idx, jdx)
             !      b = tau%upt0(is,js)
             !      if(abs(a-b)>1.d-12)then
             !        write(*,*) ii+i, jj, is, js
             !        write(*,*) a, b
             !      endif
             !      idx = (jj-1) * n + is
             !      jdx = (ii-1+i) * n + js
             !      a = Aup(idx, jdx)
             !      b = tau%up0t(is,js)
             !      if(abs(a-b)>1.d-12)then
             !        write(*,*) jj, ii+1, is, js
             !        write(*,*) a, b
             !      endif
             !   enddo
             !enddo
           enddo

           idir = 2
           call reset_gtau(0, 0)
           do i = -1, -m, -1
             call change_gtau_time(idir, tau, G_up, G_dn)
             call save_gtau(i,0)
             !DEBUG
             !do is=1,n
             !   do js = 1, n
             !      idx = (ii-1+i) * n + is
             !      jdx = (jj-1) * n + js
             !      a = Aup(idx, jdx)
             !      b = tau%upt0(is,js)
             !      if(abs(a-b)>1.d-12)then
             !        write(*,*) ii+i, jj, is, js
             !        write(*,*) a, b
             !      endif
             !      idx = (jj-1) * n + is
             !      jdx = (ii-1+i) * n + js
             !      a = Aup(idx, jdx)
             !      b = tau%up0t(is,js)
             !      if(abs(a-b)>1.d-12)then
             !        write(*,*) jj, ii+i, is, js
             !        write(*,*) a, b
             !      endif
             !   enddo
             !enddo
           enddo

           idir = 3
           do i = -m, m
              call reset_gtau(i,0)
              do j = 1, m
                call change_gtau_time(idir, tau, G_up, G_dn)
                call save_gtau(i,j)
                !DEBUG
                !do is=1,n
                !   do js = 1, n
                !      idx = (ii-1+i) * n + is
                !      jdx = (jj-1+j) * n + js
                !      a = Aup(idx, jdx)
                !      b = tau%upt0(is,js)
                !      if(abs(a-b)>1.d-12)then
                !        write(*,*) ii+i, jj+j, is, js
                !        write(*,*) a, b
                !      endif
                !      idx = (jj-1+j) * n + is
                !      jdx = (ii-1+i) * n + js
                !      a = Aup(idx, jdx)
                !      b = tau%up0t(is,js)
                !      if(abs(a-b)>1.d-12)then
                !        write(*,*) jj+j, ii+i, is, js
                !        write(*,*) a, b
                !      endif
                !   enddo
                !enddo
              enddo
           enddo

           idir = 4
           do i = -m, m
              call reset_gtau(i,0)
              do j = -1, -m, -1
                call change_gtau_time(idir, tau, G_up, G_dn)
                call save_gtau(i,j)
                !DEBUG
                !do is=1,n
                !   do js = 1, n
                !      idx = (ii-1+i) * n + is
                !      jdx = (jj-1+j) * n + js
                !      a = Aup(idx, jdx)
                !      b = tau%upt0(is,js)
                !      if(abs(a-b)>1.d-12)then
                !        write(*,'(A,4i4)') 't0', ii+i, jj+j, is, js
                !        write(*,*) a, b
                !      endif
                !      idx = (jj-1+j) * n + is
                !      jdx = (ii-1+i) * n + js
                !      a = Aup(idx, jdx)
                !      b = tau%up0t(is,js)
                !      if(abs(a-b)>1.d-12)then
                !        write(*,'(A,4i4)') '0t',jj+j, ii+i, is, js
                !        write(*,*) a, b
                !      endif
                !   enddo
                !enddo
              enddo
           enddo

           ! Compute Fourier transform in space

           do i = -m, m
              do j = -m, m
                 do is = 1, na
                    do js =  1, na

                       ! Fourier transform G_up(i,j)
                       tmp = up_ij(is:n:na, js:n:na ,i, j)
                       call Get_Fourier_Trans(tmp, fwgt, nt)
                       up_ij(is:n:na, js:n:na ,i, j) = tmp

                       ! Fourier transform G_up(j,i)
                       tmp = up_ji(is:n:na, js:n:na ,i, j)
                       call Get_Fourier_Trans(tmp, fwgt, nt)
                       up_ji(is:n:na, js:n:na ,i, j) = tmp

                       ! Fourier transform G_dn(i,j)
                       tmp = dn_ij(is:n:na, js:n:na ,i, j)
                       call Get_Fourier_Trans(tmp, fwgt, nt)
                       dn_ij(is:n:na, js:n:na ,i, j) = tmp

                       ! Fourier transform G_dn(j,i)
                       tmp = dn_ji(is:n:na, js:n:na ,i, j)
                       call Get_Fourier_Trans(tmp, fwgt, nt)
                       dn_ji(is:n:na, js:n:na ,i, j) = tmp

                    enddo
                 enddo
              enddo
           enddo
           
           ! Fill in the full G : now a function of k and time.
           do i = -m , m
              do j = -m, m
                 idx = (ii + i -1) * n
                 jdx = (jj + j -1) * n
                 do is = 1, n
                    do js = 1, n
                       Gkw%Gkw_up(idx+is,jdx+js) = up_ij(is,js,i,j) / nt
                       Gkw%Gkw_up(jdx+is,idx+js) = up_ji(is,js,j,i) / nt
                       Gkw%Gkw_dn(idx+is,jdx+js) = dn_ij(is,js,i,j) / nt
                       Gkw%Gkw_dn(jdx+is,idx+js) = dn_ji(is,js,j,i) / nt
                    enddo
                 enddo
              enddo
           enddo

        enddo
     enddo

     !Write G_ij(k,k,tau)
     open(unit=98, file='gktau.dat')
     jj = 2
     b = acos(-1.d0) / (Gkw%minw * L)
     do idx = 1, nt
        a = 0.d0
        do ii = jj, jj + L -1
           write(98,'(e18.9)', advance='no') a 
           it = mod( ii -1 , L ) + 1
           i = (it - 1) * n + (idx - 1) * na
           j = (jj - 1) * n + (idx - 1) * na
           do is = 1, na
              do js = 1, na
                write(98,'(e18.9)', advance='no')real(Gkw%Gkw_up(i+is,j+js))
              enddo
           enddo
           write(98,*)
           a = a + b
        enddo
        write(98,*)
        write(98,*)
     enddo
     close(98)

     ! Fill in matrix of Fourier coefficients
     deallocate(tmp)
     allocate(tmp(L, L))
     fwgt => Gkw%wgt_tw

     allocate(tmp2(na*na, L, nt))

     open(unit=99, file='gkw.dat')

     ! Transform up-spin matrix first : i = 1
     c = Gkw%minw / acos(-1.d0)
     Gptr => Gkw%Gkw_up
     do i = 1, 2

        do is = 1, n
           do js = 1, n

              ! Load the L-by-L matrix to be Fourier transformed.
              do ii = 1, L
                 do jj = 1, L
                    idx = (ii - 1) * n + is
                    jdx = (jj - 1) * n + js
                    tmp(ii, jj) = Gptr(idx, jdx)
                 enddo
              enddo

              call Get_Fourier_Trans(tmp, fwgt, L)

              ! Transfer the content back in B
              do ii = 1, L
                 do jj = 1, L
                    idx = (ii - 1) * n + is
                    jdx = (jj - 1) * n + js
                    Gptr(idx, jdx) = tmp(ii, jj) * c
                 enddo
              enddo

           enddo
        enddo

        do ii = 1, L
           do jj = 1, L
              do it  = 1, nt
                 do jt = 1, nt
                    m = 0
                    do is = 1, na
                       do js = 1, na
                          m = m + 1
                          idx = (ii - 1) * n + (it-1) * na + is
                          jdx = (jj - 1) * n + (jt-1) * na + js 
                          if(ii/=jj.or.it/=jt)then
                             if(abs(Gptr(idx,jdx))>1.d-10)then
                                write(*,*)ii,jj,it,jt
                                write(*,*)Gptr(idx,jdx)
                             endif
                          else
                             tmp2(m, ii, it) = Gptr(idx, jdx)
                          endif
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo

        if (i==1)then
           do it = 1, nt
              a = Gkw%minw
              b = 2*a
              do ii = 1, L
                  write(99,'(9e18.9)') a, (tmp2(is, ii, it), is=1,na*na)
                  a = a + b
              enddo
              write(99,*)
              write(99,*)
           enddo
        endif

        ! Switch to down spin matrix : i = 2
        Gptr => Gkw%Gkw_dn
     enddo
 
     close(99)
   
     !write(*,*)'Checking Full G matrix'
     !do ii = 1, L
     !   do jj = 1, L
     !      do i = 1, n
     !         do j = 1, n
     !            idx = (ii-1) * n + i
     !            jdx = (jj-1) * n + j
     !            a = Aup(idx, jdx)
     !            b = Gkw%Gkw_up(idx, jdx)
     !            if(abs(a-b)>1.d-12)then
     !               write(*,*)'Block  :',ii, jj
     !               write(*,*)'Sites  :',i, j
     !               write(*,*)'Aup', a
     !               write(*,*)'Bup', b
     !               write(*,*)
     !            endif
     !            a = Adn(idx, jdx)
     !            b = Gkw%Gkw_dn(idx, jdx)
     !            if(abs(a-b)>1.d-12)then
     !               write(*,*)'Block  :',ii, jj
     !               write(*,*)'Sites  :',i, j
     !               write(*,*)'Adn', a
     !               write(*,*)'Bdn', b
     !               write(*,*)
     !            endif
     !         enddo
     !      enddo
     !   enddo
     !enddo
    
     deallocate(up_ij) 
     deallocate(dn_ij) 
     deallocate(up_ji) 
     deallocate(dn_ji)

     !deallocate(Aup)
     !deallocate(Adn)
 
     stop

  contains

    !-----------------------------------!

    subroutine save_gtau(il,jl)

       integer, intent(in) :: il,jl

       up_ij(:,:,il,jl)=dcmplx(tau%upt0(:,:), 0.d0)
       dn_ij(:,:,il,jl)=dcmplx(tau%dnt0(:,:), 0.d0)
       up_ji(:,:,jl,il)=dcmplx(tau%up0t(:,:), 0.d0)
       dn_ji(:,:,jl,il)=dcmplx(tau%dn0t(:,:), 0.d0)

    end subroutine

    !-----------------------------------!

    subroutine reset_gtau(il, jl)

       integer, intent(in) :: il,jl

       tau%upt0 = real(up_ij(:,:,il,jl))
       tau%dnt0 = real(dn_ij(:,:,il,jl))
       tau%up0t = real(up_ji(:,:,jl,il))
       tau%dn0t = real(dn_ji(:,:,jl,il))

       tau%ii = ii + il
       tau%ib = jj + jl

    end subroutine

  end subroutine

  !-------------------------------------------------------------!

  subroutine Get_Fourier_Trans( gt, fwgt, L )
     ! This subroutine computes the Fourier transform of
     ! a function of two variables as two successive
     ! transformation. It does not use any fast-transform
     ! algorithm. Scaling is L^3. Fourier weights 
     ! are computed internally.

     implicit none
     
     ! ... Arguments ...
     integer, intent(in)       :: L 
     complex*16, intent(inout) :: gt(L, L)
     complex*16, intent(inout)    :: fwgt(L, L)

     ! ... Local variables ...
     integer :: iom, it, jom, jt
     complex*16 :: gom(L, L)
     
     ! ... Executable ...

     ! Transform 1st coordinate
     gom = 0.d0
     do iom = 1, L
        do jt = 1, L
           do it = 1, L
              gom(jt,iom) = gom(jt,iom) + gt(it, jt) * fwgt(it, iom)
           enddo
        enddo
     enddo

     ! Make complex conjugate
     fwgt = dconjg(fwgt)

     ! Transform 2nd coordinate
     gt = 0.d0
     do jom = 1, L
        do iom = 1, L
           do jt = 1, L
              gt(iom,jom) = gt(iom, jom) + gom(jt, iom) * fwgt(jt, jom)
           enddo
        enddo
     enddo

     ! Convert back to original
     fwgt = dconjg(fwgt)

  end subroutine

  !-------------------------------------------------------------!

  subroutine Gkw_Init(Gkw, RecipLatt, Lattice, L, dtau)
     use dqmc_Latt
     use dqmc_Reclatt
     use FullG

     implicit none

     type(lattice_t), intent(in) :: Lattice
     type(recip_lattice_t), intent(in) :: RecipLatt
     type(gkw_t), intent(inout) :: Gkw
     integer, intent(in) :: L 
     real*8, intent(in) :: dtau

     integer :: nt, it, ik, natl

     real*8 :: k_dot_r, minw
     real*8, pointer :: tr(:,:), kpts(:,:)

     complex*16 :: wgt, dwgt

     complex*16, pointer :: twgt(:, :)
 
     nt   =  Lattice%ncell
     tr   => Lattice%translation
     kpts => RecipLatt%klist

     allocate(Gkw%wgt_rk(nt,nt))

     do it = 0, nt-1
        do ik = 1, RecipLatt%nkpts
           k_dot_r = sum(tr(:,it)*kpts(ik,:))
           Gkw%wgt_rk(it+1,ik) = dcmplx(cos(k_dot_r),sin(k_dot_r))
        enddo
     enddo

     allocate(Gkw%wgt_tw(L, L))
     twgt => Gkw%wgt_tw
 
     minw = pi / L
     
     wgt = dcmplx ( cos(minw) , sin(minw) )
     dwgt = wgt * wgt
     do ik = 1, L
        twgt(1,ik) = dcmplx(1.d0, 0.d0)
        do it = 2, L
           twgt(it, ik) = twgt(it-1,ik) * wgt
        enddo
        wgt = wgt * dwgt
     enddo
     twgt = twgt * dtau

     Gkw%nk = nt
     Gkw%nw = L
     Gkw%na = Lattice%natom
     Gkw%minw = minw / dtau

     natl = Gkw%na * nt * L

     allocate(Gkw%Gkw_up(natl, natl))
     allocate(Gkw%Gkw_dn(natl, natl))

  end subroutine

  subroutine GetG2(G, nk, nw, Q, W, G2)
    ! Given G(k,k',w,w'), the total momentum Q and total frequency,
    ! it forms G2(k,k',w,w') for a total
    write(*,*)"To be started"
  end subroutine
