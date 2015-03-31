module DQMC_HAMILT

use DQMC_GEOM_PARAM
use DQMC_LATT
use DQMC_RECLATT
use DQMC_Cfg

implicit none

 type :: Hamiltonian_t
 
   !number of sites having a non-zero J,U,t
   integer                :: nJsites, nUsites, ntsites            
 
   !maximum number of neighbors
   integer                :: maxneig                           
 
   !Number of neighbors of each site (0:nsites-1)
   integer, pointer       :: tnneig(:), Unneig(:), Jnneig(:)
 
   !Sites havind a non-zero J,U,t (0:nsites-1)
   integer, pointer       :: tsite(:), Usite(:), Jsite(:)
 
   !Sites neighboring each site   (0:nsites-1,nsites)
   integer, pointer       :: tneig(:,:), Uneig(:,:), Jneig(:,:)  
 
   !chemical potential
   real*8                 :: mu_up, mu_dn
 
   !value of U and J for each pair (0:nsites-1, 0:nsites-1)
   real*8, pointer        :: Uv(:,:), Jv(:,:)                   
 
   !values of U and mu inequivalent by symmetry (nlocclass)
   real*8, pointer        :: Uvalue(:), muupvalue(:), mudnvalue(:)
 
   !value of t for each pair (0:nsites-1, 0:nsites-1)
   complex*16, pointer    :: hopup(:,:), hopdn(:,:)
 
   !values of t inequivalent by symmetry (nhopclass)
   complex*16, pointer    :: tupvalue(:), tdnvalue(:)
 
   !wave function phase(not used in QMC)
   complex*16, pointer    :: phase(:)
 
   !number of different hoppings
   integer                :: nhopclass
 
   !number of sites with different U/mu
   integer                :: nlocclass
 
   !class for each site having U/mu (0:nsites-1)
   integer, pointer       :: mylocclass(:)
 
   !hopping class for each pair of neighboring site (0:nsites-1,maxval(tnneig))
   integer, pointer       :: myhopclass(:,:)
 
   logical                :: constructed
   logical                :: neig_found
   logical                :: analyzed
 
 end type

contains

!-----------------------------------------------------------------!

 subroutine free_hamilt(hamilt)
   !
   ! Free the space taken by pointers in Hamilt.
   !
    type(hamiltonian_t), intent(inout) :: hamilt

    if(associated(hamilt%tnneig))     deallocate(hamilt%tnneig)
    if(associated(hamilt%Unneig))     deallocate(hamilt%Unneig)
    if(associated(hamilt%Jnneig))     deallocate(hamilt%Jnneig) 
    if(associated(hamilt%tsite))      deallocate(hamilt%tsite)
    if(associated(hamilt%Usite))      deallocate(hamilt%Usite)
    if(associated(hamilt%Jsite))      deallocate(hamilt%Jsite)  
    if(associated(hamilt%Uvalue))     deallocate(hamilt%Uvalue) 
    if(associated(hamilt%muupvalue))  deallocate(hamilt%muupvalue) 
    if(associated(hamilt%mudnvalue))  deallocate(hamilt%mudnvalue) 
    if(associated(hamilt%tupvalue))   deallocate(hamilt%tupvalue)    
    if(associated(hamilt%tdnvalue))   deallocate(hamilt%tdnvalue)    
    if(associated(hamilt%phase))      deallocate(hamilt%phase)     
    if(associated(hamilt%mylocclass)) deallocate(hamilt%mylocclass) 
    if(associated(hamilt%tneig))      deallocate(hamilt%tneig)
    if(associated(hamilt%Uneig))      deallocate(hamilt%Uneig)
    if(associated(hamilt%Jneig))      deallocate(hamilt%Jneig) 
    if(associated(hamilt%Uv))         deallocate(hamilt%Uv)
    if(associated(hamilt%Jv))         deallocate(hamilt%Jv)   
    if(associated(hamilt%hopup))      deallocate(hamilt%hopup)  
    if(associated(hamilt%hopdn))      deallocate(hamilt%hopdn)  
    if(associated(hamilt%myhopclass)) deallocate(hamilt%myhopclass) 

 end subroutine free_hamilt

!-----------------------------------------------------------------!

 subroutine construct_hamilt(hamilt, lattice, recip_lattice, cfg)
   !
   ! Fill the hamiltonian : Fill all the variables making up 
   ! hamiltonian_t except those whose name end in "class'. Set 
   ! constructed and neig_found to true
   !
    type(lattice_t), target, intent(in) :: lattice
    type(recip_lattice_t),intent(in)    :: recip_lattice
    type(Hamiltonian_t),intent(out)     :: hamilt
    type(config), intent(inout)         :: cfg

    integer              :: maxhop, ntcfg, natom, nsites
    integer              :: ios, icount, ihop, j, iat, jat
    real*8               :: hop3d(rdim), tijtmpup, tijtmpdn, Utmp 
    real*8               :: ktwist(rdim), kpoint(rdim)
    real*8, pointer      :: Uv(:,:), Jv(:,:), tcfg(:), pos(:,:)
    complex*16, pointer  :: hopup(:,:), hopdn(:,:), phase(:)
    integer, allocatable :: nhop(:), hopto(:,:), hoptarget(:,:)
    real*8, allocatable  :: tijup(:,:), tijdn(:,:), Uat(:,:), Jint(:,:) 
    real*8, allocatable  :: twisthop(:,:), ahop(:,:,:)
    character*50         :: string
    logical              :: ldum

    character(len=*), parameter :: mu(2) = (/'mu_up','mu_dn'/)

    pos => lattice%cartpos

    if(.not.lattice%constructed) &
       stop'Need to construct lattice before building Hamiltonian'
    if(.not.recip_lattice%initialized) &
       stop'Need to initialize recip_lattice before building Hamiltonian'

    !Set local alias
    natom = lattice%natom
    nsites = lattice%nsites
    ktwist(1:rdim) = recip_lattice%ktwist(1:rdim)
    kpoint(1:rdim) = recip_lattice%kpoint(1:rdim)

    allocate(hopup(0:nsites-1, 0:nsites-1))
    allocate(hopdn(0:nsites-1, 0:nsites-1))
    allocate(Uv(0:nsites-1, 0:nsites-1)) 
    allocate(Jv(0:nsites-1, 0:nsites-1))
    allocate(phase(0:nsites-1))

    !Read chemical potential
    do iat = 1, 2
       call CFG_Get(cfg, mu(iat), ntcfg, tcfg)
       if(ntcfg > 1)then
         write(*,'(A)')'WARNING: Only 1st entry for mu in input file is considered.'
       endif
       if (iat == 1) hamilt%mu_up = tcfg(1)
       if (iat == 2) hamilt%mu_dn = tcfg(1)
    enddo
    deallocate(tcfg)
    nullify(tcfg)

    !Find the hamiltonian field
    ldum = move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)
    allocate(nhop(0:natom-1))

    !Count number of entry in Hamilt
    call count_nhop(nhop, natom) 
    maxhop=maxval(nhop)

    allocate(hopto(0:nsites-1, maxhop))
    allocate(hoptarget(0:natom-1, maxhop))
    allocate(tijup(0:natom-1, maxhop))
    allocate(tijdn(0:natom-1, maxhop))
    allocate(Uat(0:natom-1, maxhop))
    allocate(twisthop(0:nsites-1, maxhop))
    allocate(Jint(0:natom-1, maxhop))
    allocate(ahop(rdim,0:natom-1, maxhop))

    nhop(:) = 0
    ahop    = 0.d0 
    tijup   = 0.d0 
    tijdn   = 0.d0 
    Jint    = 0.d0 
    Uat     = 0.d0

    ldum=move_to_record(INPUT_FIELDS(HAMILT_F),inpunit)

    do
       !Read next line in Hamilt
       read(inpunit,'(A)') string
       read(string, *, iostat=ios) iat, jat, hop3d(1), hop3d(2), hop3d(3), &
           tijtmpup, tijtmpdn, Utmp
       if(ios .ne. 0) exit

       !number of hopping from iat
       nhop(iat) = nhop(iat) + 1

       !species of hopping target
       hoptarget(iat, nhop(iat)) = jat

       !hopping vector
       ahop(1, iat, nhop(iat)) = hop3d(1)
       ahop(2, iat, nhop(iat)) = hop3d(2)
       ahop(3, iat, nhop(iat)) = hop3d(3)

       !hopping matrix element
       tijup(iat, nhop(iat)) = tijtmpup
       tijdn(iat, nhop(iat)) = tijtmpdn
       !Interaction matrix element
       Uat(iat, nhop(iat)) = Utmp

       if( jat==iat .and. sum(hop3d**2)<1.d-10 )cycle
       !Define inverse hopping

       nhop(jat) = nhop(jat) + 1

       hoptarget(jat,nhop(jat)) = iat

       ahop(1, jat, nhop(jat)) = -hop3d(1)
       ahop(2, jat, nhop(jat)) = -hop3d(2)
       ahop(3, jat, nhop(jat)) = -hop3d(3)

       tijup(jat,nhop(jat)) = tijtmpup
       tijdn(jat,nhop(jat)) = tijtmpdn
       Uat(jat,nhop(jat)) = Utmp

    enddo

    !Define phase on each atom compatibly with BC
    do iat = 0, nsites - 1
       jat = mod(iat,natom)
       phase(iat) = exp(im*sum((kpoint(:)-ktwist(:))*(pos(:,iat)-pos(:,jat))))
       phase(iat)= phase(iat) * exp(im*(-sum(ktwist(:)*pos(:,iat))))
    enddo

    !construct hopping table
    do iat = 0, nsites-1
       icount = mod(iat,natom)
       do ihop = 1, nhop(icount)
          jat = hoptarget(icount,ihop)
          hopto(iat,ihop) = hoptowho(iat,ahop(1,icount,ihop),jat,lattice)
          jat = hopto(iat,ihop)
          hop3d(:) = pos(:,jat) - pos(:,iat) - ahop(:,icount,ihop)
          !modify hopping to take phase difference into account
          twisthop(iat,ihop) = sum(ktwist(:)*hop3d(:))
       enddo
    enddo

    !construct hamiltonian matrix for t, U and J
    hopup(:,:) = 0.d0 
    hopdn(:,:) = 0.d0 
    Jv(:,:)  = 0.d0 
    Uv(:,:)  = 0.d0
    do iat = 0, nsites-1
       !iat -> neigbors
       icount = mod(iat,natom)
       do ihop = 1, nhop(icount)
        j = hopto(iat,ihop)
        !hop is the hopping matrix element <j|T|iat> for c_j+ c_iat
        hopup(iat,j) = hopup(iat,j) + &
        & tijup(icount,ihop) * exp(im*twisthop(iat,ihop))
        hopdn(iat,j) = hopdn(iat,j) + &
        & tijdn(icount,ihop) * exp(im*twisthop(iat,ihop))
        Jv(iat,j)  = Jv(iat,j)  + Jint(icount,ihop)
        Uv(iat,j)  = Uv(iat,j)  + Uat(icount,ihop)
       enddo
    enddo

    hamilt%Uv => Uv
    hamilt%Jv => Jv
    hamilt%hopup => hopup
    hamilt%hopdn => hopdn
    hamilt%phase => phase

    hamilt%constructed = .true.

    call find_neighbors(hamilt)
    hamilt%neig_found = .true.

    deallocate(nhop)
    deallocate(hopto)
    deallocate(hoptarget)
    deallocate(tijup)
    deallocate(tijdn)
    deallocate(Uat)
    deallocate(twisthop)
    deallocate(Jint)
    deallocate(ahop)

 end subroutine 

!-------------------------------------------------------------------!

 subroutine count_nhop(nhop,natom)
   !
   ! Count the number of hops that are possible 
   ! from each orbitals in the primitive cell.
   ! Fills up nhop
   !
     integer, intent(in)  :: natom
     integer, intent(out) :: nhop(0:natom-1)

     integer              :: iat, jat, ios
     character*50         :: string
     real*8               :: hop(rdim), ttmp, Utmp

     nhop=0
     do
        !Read next line in input
        read(inpunit,'(A)')string
        read(string,*,iostat=ios)iat,jat,hop(1),hop(2),hop(3),ttmp,ttmp,Utmp

        !Line is not part of Hamilt. Quit.
        if (ios.ne.0) exit

        !Atom label is outside range
        if (iat>natom-1.or.jat>natom-1) then
           write(*,*)'One of the atom in hopping is unspecified' 
           stop
        endif

        !Update value of nhop
        nhop(iat) = nhop(iat) + 1
        if (jat.ne.iat .or. sum(hop**2)>1.d-10) nhop(jat) = nhop(jat) + 1

     enddo

     rewind(inpunit)

 end subroutine count_nhop

!-------------------------------------------------------------------!

 subroutine find_neighbors(hamilt)
   !
   ! Given the matrices Uv,Jv and hop it finds which atoms 
   ! are neighbors and classify neighbors as with respect 
   ! to interaction or hopping. Returns tnneig(is): number
   ! of neighbors of is; tneig(is,j): j-th neighbor of "is".
   ! Analogous definition for Unneig,Uneig and Jnneig,Jneig.
   !
     type(Hamiltonian_t),intent(inout) :: hamilt

     integer is, js, ntsites, nUsites, nJsites, n
     integer, pointer :: tnneig(:),Unneig(:),Jnneig(:)
     integer, pointer :: tneig(:,:),Uneig(:,:),Jneig(:,:)
     integer, pointer :: tsite(:),Usite(:),Jsite(:)
     
     if(.not.hamilt%constructed) stop'Hamiltonian needs to be &
        &constructed before neig can be found'
    
     n=size(hamilt%hopup,1)

     if(associated(tnneig)) then
        deallocate(tnneig, Unneig, Jnneig)
        deallocate(tneig, Uneig, Jneig)
        deallocate(tsite, Usite, Jsite)
     endif

     allocate(tnneig(0:n-1),Unneig(0:n-1),Jnneig(0:n-1))
     allocate(tneig(0:n-1,n),Uneig(0:n-1,n),Jneig(0:n-1,n))
     allocate(tsite(n),Usite(n),Jsite(n))
    
     !Count the number of sites connected to "is" by t, U and J and store site label
     tnneig = 0 
     Unneig = 0 
     Jnneig = 0
     do is = 0, n-1
        do js = 0, n-1
           if (abs(hamilt%hopup(is,js)).gt.1d-9 .and. is/=js) then
              tnneig(is) = tnneig(is) + 1
              tneig(is,tnneig(is)) = js
           elseif (abs(hamilt%hopdn(is,js)).gt.1d-9 .and. is/=js) then
              tnneig(is) = tnneig(is) + 1
              tneig(is,tnneig(is)) = js
           endif  
           if (abs(hamilt%Uv(is,js)).gt.1d-9) then
              Unneig(is) = Unneig(is) + 1
              Uneig(is,Unneig(is)) = js
           endif
           if (abs(hamilt%Jv(is,js)).gt.1d-9) then
              Jnneig(is) = Jnneig(is) + 1
              Jneig(is,Jnneig(is)) = js
           endif
        enddo
     enddo
    
     !Count the number of sites effectively having t,U or J on them
     ntsites = 0 
     nUsites = 0 
     nJsites = 0
     do is = 0, n-1
        if (tnneig(is).gt.0) then
           ntsites = ntsites + 1
           tsite(ntsites) = is
        endif 
        if(Unneig(is).gt.0)then
           nUsites = nUsites + 1
           Usite(nUsites) = is
        endif
        if(Jnneig(is).gt.0)then
           nJsites = nJsites + 1
           Jsite(nJsites) = is
        endif
     enddo
    
     !Fill hamiltonian variables
     hamilt%maxneig = max(maxval(tnneig), maxval(Unneig), maxval(Jnneig))
    
     !number of neighbors of each site
     hamilt%tnneig => tnneig
     hamilt%Unneig => Unneig
     hamilt%Jnneig => Jnneig
    
     !neighbors of each site
     hamilt%tneig => tneig
     hamilt%Uneig => Uneig
     hamilt%Jneig => Jneig
    
     !number of sites on which either t or U or J is different from 0
     hamilt%ntsites = ntsites
     hamilt%nUsites = nUsites
     hamilt%nJsites = nJsites
    
     !sites on which either t or U or J is different from 0
     hamilt%tsite => tsite
     hamilt%Usite => Usite
     hamilt%Jsite => Jsite
    
 end subroutine find_neighbors

!-------------------------------------------------------------------!

 subroutine count_hop_class(lattice, hamilt)
   !
   ! Construct a list of only hopping classes. These are a 
   ! subset of all the distance classes for which the hopping 
   ! is non zero. Find which hoppings are equivalent by 
   ! symmetry. nhopsite(iat,iclass) returns how many hoppings 
   ! in class iclass iat has. ordered_neig(iat,ineig) returns 
   ! the ineig-th neighbor of iat. Neighbors are here ordered 
   ! according to classes. So if nhopsite(2,1)=3 the first 
   ! three neigbors of 2 listed in ordered_neig are those
   ! belonging to class 1. tneig has a similar content but 
   ! neighbors are not ordered.
   !
    type(Hamiltonian_t)  :: hamilt
    type(lattice_t)      :: lattice

    integer              :: isite, ineig, jsite, newclass, ihc
    integer              :: maxtneig, nsites 
    integer              :: hopclass(lattice%nclass),nhopclass
    integer, allocatable :: pairhopclass(:,:)
    complex*16           :: tvaluetmpup(lattice%nclass)
    complex*16           :: tvaluetmpdn(lattice%nclass)
    
    nsites      = size(hamilt%tnneig)
    nhopclass   = 0
    hopclass(:) = 0
    
    allocate(pairhopclass(0:nsites-1,nsites))
    !class for the pair (isite, t-neighbour of site). 
    pairhopclass(:,:) = -1
    
    !Select classes for which hopping is different from 0
    !Uses symmetry info found in lattice.
    do isite = 0, nsites-1
       do ineig = 1, hamilt%tnneig(isite)
    
          jsite = hamilt%tneig(isite,ineig) 
          !(isite,jsite) is a pair of sites with non-zero hopping
          if(isite == jsite) cycle
          newclass = lattice%myclass(isite,jsite)
    
          !First check that this class has not been already found
          do ihc = 1, nhopclass
             if(newclass == hopclass(ihc))exit
          enddo
    
          !If hopclass is new, increase number of classes
          if(ihc == nhopclass+1)then 
             hopclass(ihc)  = newclass
             nhopclass      = ihc
             tvaluetmpup(ihc) = hamilt%hopup(isite,jsite)
             tvaluetmpdn(ihc) = hamilt%hopdn(isite,jsite)
          endif
    
          !Assign pair to the class
          pairhopclass(isite,ineig) = ihc
       enddo
    enddo
    hamilt%nhopclass = nhopclass
    maxtneig = maxval(hamilt%tnneig(:))
    
    allocate(hamilt%myhopclass(0:nsites-1,maxtneig))
    allocate(hamilt%tupvalue(nhopclass))
    allocate(hamilt%tdnvalue(nhopclass))

    !Assign t value to each class and class number to each t-pair
    hamilt%tupvalue(1:nhopclass) = tvaluetmpup(1:nhopclass)
    hamilt%tdnvalue(1:nhopclass) = tvaluetmpdn(1:nhopclass)

    do ineig = 1, maxtneig
       hamilt%myhopclass(0:nsites-1,ineig) = pairhopclass(0:nsites-1,ineig)
    enddo

    deallocate(pairhopclass)

 end subroutine count_hop_class

!-------------------------------------------------------------------!

 subroutine count_local_classes(lattice,hamilt)
  !
  ! Construct a list of only local classes i.e. the subset 
  ! of the distance classes defined on the same site. Note
  ! that we use the symmetry instead of the value to group
  ! sites together. This is because local classes refers to 
  ! both U and mu and two sites may have same value of U 
  ! but different mu. If equivalent by symmetry, the two 
  ! sites have, however, necessarily same u and mu.
  !
    type(Hamiltonian_t), intent(inout) :: hamilt
    type(lattice_t), intent(in)        :: lattice

    integer :: isite, iclass, jclass, nsites, nlocclass
    integer :: localtmp(lattice%nclass)
    real*8  :: Utmp(lattice%nclass) 
    real*8  :: muup(lattice%nclass)
    real*8  :: mudn(lattice%nclass)

    nlocclass = 0
    nsites    = size(hamilt%tnneig)
    allocate(hamilt%mylocclass(0:nsites-1))

    !Loop over all sites
    do isite = 0, nsites-1

       iclass = lattice%myclass(isite,isite)
    
       !Check whether class was already found
       do jclass = 1, nlocclass
          if(iclass.eq.localtmp(jclass))exit
       enddo
    
       !if not augment the number of local classes
       if(jclass>nlocclass)then
          nlocclass = nlocclass+1
          localtmp(jclass) = iclass
          !Save value of mu and U for this class
          Utmp(jclass) = hamilt%Uv(isite,isite)
          muup(jclass) = -dble(hamilt%hopup(isite,isite))+hamilt%mu_up
          mudn(jclass) = -dble(hamilt%hopdn(isite,isite))+hamilt%mu_dn
       endif
    
       !assign site to a class
       hamilt%mylocclass(isite) = jclass
    enddo
    
    !Store the number of different classes
    hamilt%nlocclass = nlocclass
    
    !Store the value of U and on-site energy (shifted by mu) for each class
    allocate(hamilt%Uvalue(nlocclass))
    allocate(hamilt%muupvalue(nlocclass))
    allocate(hamilt%mudnvalue(nlocclass))
    hamilt%Uvalue(1:nlocclass)  = Utmp(1:nlocclass)
    hamilt%muupvalue(1:nlocclass) = muup(1:nlocclass)
    hamilt%mudnvalue(1:nlocclass) = mudn(1:nlocclass)

 end subroutine

!-------------------------------------------------------------------!

 subroutine group_hopping(hamilt, n, nt, tmap, tupvalue, tdnvalue)
  !
  !  Assign to every non-zero hopping a class based 
  !  on its value. Equal hopping matrix elements belong 
  !  to the same class : tmap(i,j)=it.
  !  tvalue(it) contains the value for a given class.
  !  Exclude local site energies.
  !
    type(Hamiltonian_t), intent(in) :: hamilt
    integer, intent(in)             :: n
    integer, intent(out)            :: tmap(n,n), nt
!
! temporary fix: FORTRAN 90 does not support "intent" for pointers
!
!    real*8, pointer, intent(out)    :: tupvalue(:)
!    real*8, pointer, intent(out)    :: tdnvalue(:)
    real*8, pointer    :: tupvalue(:)
    real*8, pointer    :: tdnvalue(:)

    integer         :: is, js, it, jn
    real*8          :: tup, tdn

    !Initialize
    tmap = 0
    nt   = 0

    if (associated(tupvalue)) deallocate(tupvalue)
    if (associated(tdnvalue)) deallocate(tdnvalue)

    if(hamilt%maxneig > 0) then
       allocate(tupvalue(hamilt%ntsites*hamilt%maxneig))
       allocate(tdnvalue(hamilt%ntsites*hamilt%maxneig))
       !Assign same class to identical matrix elements
       do is = 0, n-1
          do jn = 1, hamilt%tnneig(is)
             js = hamilt%tneig(is,jn)
             if (is .eq. js) cycle
             tup  = dble(hamilt%hopup(is,js))
             tdn  = dble(hamilt%hopdn(is,js))
             do it = 1, nt
                if (abs(tup-tupvalue(it)) < 1.d-6 .and. &
                 &  abs(tdn-tdnvalue(it)) < 1.d-6) exit
             enddo
             if (it .eq. nt+1) then
               nt = it
               tupvalue(it) = tup
               tdnvalue(it) = tdn
             endif
             tmap(is+1, js+1) = it
          enddo
       enddo
    else
       nt = 1
       allocate(tupvalue(1))
       allocate(tdnvalue(1))
       tupvalue = 0.d0
       tdnvalue = 0.d0
       do is = 1, n 
          tmap(is,is) = 1
       enddo
    endif

 end subroutine group_hopping

!-------------------------------------------------------------------!

end module DQMC_HAMILT
