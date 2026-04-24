! fmm2d_plan.f90 -- Plan-based FMM interface for cfmm2d, lfmm2d, rfmm2d
!
! When source/target geometry is fixed but charges change between calls,
! build a plan once to cache the tree, sorted coordinates, and interaction
! lists. Subsequent execute calls skip those steps and run only the FMM
! traversal (steps 1-8), which is O(N).
!
! Usage (cfmm2d):
!   call cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
!   call cfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, dipstr, iper,
!        ifpgh, pot, grad, hess, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
!   call cfmm2d_destroy_plan(ier)
!
! For lfmm2d and rfmm2d, replace cfmm2d with lfmm2d or rfmm2d.
! All three FMM types share the same underlying plan (the tree is geometry-
! only). cfmm2d_destroy_plan frees the plan for all three.
!
! C callers: use the Fortran-mangled names with a trailing underscore, and
! pass all scalar arguments as pointers. See include/fmm2d/laplace.h.

module cfmm2d_plan_mod
  implicit none

  type cfmm2d_plan_t
    ! Tree dimensions
    integer :: ns      = 0
    integer :: nt      = 0
    integer :: nlevels = 0
    integer :: nboxes  = 0
    integer :: ltree   = 0
    integer :: lmptot_1 = 0  ! lmptot for nd=1; scale by nd at execute time
    integer :: nmax    = 0
    integer :: ndiv    = 0
    integer :: ldc     = 100
    integer :: iptr(8)
    integer :: mnlist1 = 0
    integer :: mnlist2 = 0
    integer :: mnlist3 = 0
    integer :: mnlist4 = 0
    real*8  :: thresh  = 0.0d0
    real*8  :: eps     = 0.0d0

    ! Tree arrays
    integer, allocatable :: itree(:)
    real*8,  allocatable :: tcenters(:,:)
    real*8,  allocatable :: boxsize(:)
    real*8,  allocatable :: rscales(:)
    integer, allocatable :: nterms(:)
    integer, allocatable :: iaddr_1(:,:)  ! iaddr for nd=1

    ! Sort permutations + sorted coordinates
    integer, allocatable :: isrc(:)
    integer, allocatable :: itarg(:)
    integer, allocatable :: isrcse(:,:)
    integer, allocatable :: itargse(:,:)
    integer, allocatable :: iexpcse(:,:)
    real*8,  allocatable :: sourcesort(:,:)
    real*8,  allocatable :: targsort(:,:)

    ! Interaction lists
    integer, allocatable :: nlist1s(:)
    integer, allocatable :: list1(:,:)
    integer, allocatable :: nlist2s(:)
    integer, allocatable :: list2(:,:)
    integer, allocatable :: nlist3s(:)
    integer, allocatable :: list3(:,:)
    integer, allocatable :: nlist4s(:)
    integer, allocatable :: list4(:,:)

    ! Binomial table
    real*8, allocatable :: carray(:,:)

    ! Expansion buffer (nd-dependent; reallocated when nd changes)
    real*8,  allocatable :: rmlexp(:)
    integer :: rmlexp_nd = 0

    ! Pre-allocated M2L scratch arrays for step 4.
    ! Indexed as mploc_hexp1tmp(nd, 0:nmax, nthreads) so each OpenMP
    ! thread has its own private workspace without repeated allocation.
    ! Reallocated at execute time when nd changes.
    complex*16, allocatable :: mploc_hexp1tmp(:,:,:)
    complex*16, allocatable :: mploc_jexp2tmp(:,:,:)
    integer :: mploc_nd = 0
  end type cfmm2d_plan_t

  type(cfmm2d_plan_t), save :: the_plan
  logical,             save :: plan_built = .false.

end module cfmm2d_plan_mod


!-----------------------------------------------------------------------
! cfmm2d_build_plan
!   Build plan for cfmm2d (Cauchy FMM).
!   Computes and stores: quadtree, Morton-sorted source/target coords,
!   interaction lists (list 1-4), binomial table, and expansion address
!   table for nd=1 (scaled at execute time for the actual nd).
!
!   INPUT:
!     eps         FMM precision
!     ns          number of sources
!     sources(2,ns) source locations
!     nt          number of targets (0 if none)
!     targ(2,nt)  target locations
!   OUTPUT:
!     ier         0 on success, non-zero on error
!-----------------------------------------------------------------------
subroutine cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  use cfmm2d_plan_mod
  implicit none

  real*8,  intent(in)  :: eps
  integer, intent(in)  :: ns, nt
  real*8,  intent(in)  :: sources(2,ns), targ(2,max(nt,1))
  integer, intent(out) :: ier

  integer :: nlmin, nlmax, ifunif, iper
  integer :: nlevels, nboxes, ltree, ndiv, idivflag
  integer :: ifcharge_t, ifdipole_t, ifpgh_t, ifpghtarg_t
  integer :: i, ibox, nmax, lmptot_1

  ier = 0
  nlmin = 0
  nlmax = 51
  ifunif = 0
  iper = 0
  ! Use worst-case flags to get the smallest (most conservative) ndiv
  ifcharge_t  = 1
  ifdipole_t  = 1
  ifpgh_t     = 1
  ifpghtarg_t = 1

  call lndiv2d(eps, ns, nt, ifcharge_t, ifdipole_t, &
       ifpgh_t, ifpghtarg_t, ndiv, idivflag)

  call pts_tree_mem(sources, ns, targ, nt, idivflag, ndiv, &
       nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree)

  ! Free any existing plan before building a new one
  if (plan_built) call cfmm2d_destroy_plan(i)

  ! Store scalars
  the_plan%ns      = ns
  the_plan%nt      = nt
  the_plan%nlevels = nlevels
  the_plan%nboxes  = nboxes
  the_plan%ltree   = ltree
  the_plan%ndiv    = ndiv
  the_plan%eps     = eps
  the_plan%ldc     = 100

  ! Build quadtree
  allocate(the_plan%itree(ltree))
  allocate(the_plan%boxsize(0:nlevels))
  allocate(the_plan%tcenters(2,nboxes))

  call pts_tree_build(sources, ns, targ, nt, idivflag, ndiv, &
       nlmin, nlmax, ifunif, iper, nlevels, nboxes, ltree, &
       the_plan%itree, the_plan%iptr, the_plan%tcenters, &
       the_plan%boxsize)

  ! Compute sort permutations
  allocate(the_plan%isrc(ns))
  allocate(the_plan%isrcse(2,nboxes))
  allocate(the_plan%itarg(max(nt,1)))
  allocate(the_plan%itargse(2,nboxes))
  allocate(the_plan%iexpcse(2,nboxes))

  do ibox = 1, nboxes
    the_plan%iexpcse(1,ibox) = 1
    the_plan%iexpcse(2,ibox) = 0
  enddo

  call pts_tree_sort(ns, sources, the_plan%itree, ltree, nboxes, &
       nlevels, the_plan%iptr, the_plan%tcenters, &
       the_plan%isrc, the_plan%isrcse)

  if (nt .gt. 0) then
    call pts_tree_sort(nt, targ, the_plan%itree, ltree, nboxes, &
         nlevels, the_plan%iptr, the_plan%tcenters, &
         the_plan%itarg, the_plan%itargse)
  endif

  ! Store sorted coordinates
  allocate(the_plan%sourcesort(2,ns))
  allocate(the_plan%targsort(2,max(nt,1)))

  call dreorderf(2, ns, sources, the_plan%sourcesort, the_plan%isrc)
  if (nt .gt. 0) &
    call dreorderf(2, nt, targ, the_plan%targsort, the_plan%itarg)

  ! Expansion orders at each level
  allocate(the_plan%rscales(0:nlevels))
  allocate(the_plan%nterms(0:nlevels))

  nmax = 0
  do i = 0, nlevels
    the_plan%rscales(i) = the_plan%boxsize(i)
    call l2dterms(eps, the_plan%nterms(i), ier)
    if (the_plan%nterms(i) .gt. nmax) nmax = the_plan%nterms(i)
  enddo
  the_plan%nmax = nmax

  ! iaddr for nd=1 (at execute time: iaddr_nd(k,j) = 1 + (iaddr_1(k,j)-1)*nd)
  allocate(the_plan%iaddr_1(2,nboxes))
  call l2dmpalloc(1, the_plan%itree, the_plan%iaddr_1, &
       nlevels, lmptot_1, the_plan%nterms)
  the_plan%lmptot_1 = lmptot_1

  ! Interaction lists
  call computemnlists(nlevels, nboxes, the_plan%itree, ltree, &
       the_plan%iptr, the_plan%tcenters, the_plan%boxsize, iper, &
       the_plan%mnlist1, the_plan%mnlist2, &
       the_plan%mnlist3, the_plan%mnlist4)

  allocate(the_plan%nlist1s(nboxes))
  allocate(the_plan%list1(the_plan%mnlist1, nboxes))
  allocate(the_plan%nlist2s(nboxes))
  allocate(the_plan%list2(the_plan%mnlist2, nboxes))
  allocate(the_plan%nlist3s(nboxes))
  allocate(the_plan%list3(the_plan%mnlist3, nboxes))
  allocate(the_plan%nlist4s(nboxes))
  allocate(the_plan%list4(the_plan%mnlist4, nboxes))

  call computelists(nlevels, nboxes, the_plan%itree, ltree, &
       the_plan%iptr, the_plan%tcenters, the_plan%boxsize, iper, &
       the_plan%mnlist1, the_plan%nlist1s, the_plan%list1, &
       the_plan%mnlist2, the_plan%nlist2s, the_plan%list2, &
       the_plan%mnlist3, the_plan%nlist3s, the_plan%list3, &
       the_plan%mnlist4, the_plan%nlist4s, the_plan%list4)

  ! Binomial table
  allocate(the_plan%carray(0:the_plan%ldc, 0:the_plan%ldc))
  call init_carray(the_plan%carray, the_plan%ldc)

  the_plan%thresh = the_plan%boxsize(0) * 2.0d0**(-51)

  plan_built = .true.
  return
end subroutine cfmm2d_build_plan


!-----------------------------------------------------------------------
! cfmm2d_execute_plan
!   Run the Cauchy FMM using a previously built plan.
!   Reorders charges, allocates output buffers, calls cfmm2dmain_pre,
!   then reorders outputs back to original ordering.
!   The expansion buffer (rmlexp) is reallocated when nd changes.
!
!   INPUT:
!     nd          number of charge densities
!     ifcharge    1 to include charge interactions
!     charge(nd,ns) complex charge strengths
!     ifdipole    1 to include dipole interactions
!     dipstr(nd,ns) complex dipole strengths
!     iper        flag for periodic (currently unused, pass 0)
!     ifpgh       1/2/3 = pot / pot+grad / pot+grad+hess at sources
!     ifpghtarg   1/2/3 = pot / pot+grad / pot+grad+hess at targets
!   OUTPUT:
!     pot(nd,ns)       potential at sources  (if ifpgh >= 1)
!     grad(nd,ns)      d/dz at sources       (if ifpgh >= 2)
!     hess(nd,ns)      d^2/dz^2 at sources   (if ifpgh >= 3)
!     pottarg(nd,nt)   potential at targets  (if ifpghtarg >= 1)
!     gradtarg(nd,nt)  d/dz at targets       (if ifpghtarg >= 2)
!     hesstarg(nd,nt)  d^2/dz^2 at targets   (if ifpghtarg >= 3)
!     ier         0 on success, 4 if no plan built
!-----------------------------------------------------------------------
subroutine cfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, dipstr, &
     iper, ifpgh, pot, grad, hess, &
     ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  use cfmm2d_plan_mod
  implicit none

  integer,     intent(in)    :: nd, ifcharge, ifdipole, iper
  integer,     intent(in)    :: ifpgh, ifpghtarg
  complex*16,  intent(in)    :: charge(nd,*), dipstr(nd,*)
  complex*16,  intent(inout) :: pot(nd,*), grad(nd,*), hess(nd,*)
  complex*16,  intent(inout) :: pottarg(nd,*), gradtarg(nd,*), hesstarg(nd,*)
  integer,     intent(out)   :: ier

  integer :: ns, nt, nlevels, nboxes, lmptot, lmptmp, nthreads
  integer :: i, idim, ibox
  integer, allocatable :: iaddr(:,:)
  integer :: omp_get_max_threads
  external omp_get_max_threads

  complex*16, allocatable :: chargesort(:,:), dipstrsort(:,:)
  complex*16, allocatable :: potsort(:,:), gradsort(:,:), hesssort(:,:)
  complex*16, allocatable :: pottargsort(:,:), gradtargsort(:,:)
  complex*16, allocatable :: hesstargsort(:,:)
  complex*16, allocatable :: mptemp(:)

  ! Dummy variables for unused expansion-center feature (nexpc=0)
  real*8     :: expc(2), scj(1)
  complex*16 :: jexps(1,0:0,1)
  integer    :: nexpc, ntj, ifnear
  real*8     :: timeinfo(8)

  ier = 0

  if (.not. plan_built) then
    ier = 4
    return
  endif

  ns      = the_plan%ns
  nt      = the_plan%nt
  nlevels = the_plan%nlevels
  nboxes  = the_plan%nboxes
  nexpc   = 0
  ntj     = 0
  ifnear  = 1
  do i = 1, 8
    timeinfo(i) = 0
  enddo

  ! Scale iaddr from nd=1 to actual nd
  allocate(iaddr(2, nboxes))
  do ibox = 1, nboxes
    iaddr(1,ibox) = 1 + (the_plan%iaddr_1(1,ibox) - 1) * nd
    iaddr(2,ibox) = 1 + (the_plan%iaddr_1(2,ibox) - 1) * nd
  enddo

  ! (Re)allocate expansion buffer if nd changed
  lmptot = 1 + (the_plan%lmptot_1 - 1) * nd
  if (the_plan%rmlexp_nd .ne. nd) then
    if (allocated(the_plan%rmlexp)) deallocate(the_plan%rmlexp)
    allocate(the_plan%rmlexp(lmptot), stat=ier)
    if (ier .ne. 0) return
    the_plan%rmlexp_nd = nd
  endif

  ! (Re)allocate per-thread M2L scratch arrays if nd changed
  nthreads = omp_get_max_threads()
  if (the_plan%mploc_nd .ne. nd) then
    if (allocated(the_plan%mploc_hexp1tmp)) &
      deallocate(the_plan%mploc_hexp1tmp)
    if (allocated(the_plan%mploc_jexp2tmp)) &
      deallocate(the_plan%mploc_jexp2tmp)
    allocate(the_plan%mploc_hexp1tmp(nd, 0:the_plan%nmax, nthreads), &
             stat=ier)
    if (ier .ne. 0) return
    allocate(the_plan%mploc_jexp2tmp(nd, 0:the_plan%nmax, nthreads), &
             stat=ier)
    if (ier .ne. 0) return
    the_plan%mploc_nd = nd
  endif

  ! Allocate and reorder charges/dipoles
  if (ifcharge.eq.1 .and. ifdipole.eq.0) then
    allocate(chargesort(nd,ns), dipstrsort(nd,1))
  else if (ifcharge.eq.0 .and. ifdipole.eq.1) then
    allocate(chargesort(nd,1), dipstrsort(nd,ns))
  else if (ifcharge.eq.1 .and. ifdipole.eq.1) then
    allocate(chargesort(nd,ns), dipstrsort(nd,ns))
  else
    allocate(chargesort(nd,1), dipstrsort(nd,1))
  endif

  if (ifcharge .eq. 1) &
    call dreorderf(2*nd, ns, charge,  chargesort, the_plan%isrc)
  if (ifdipole .eq. 1) &
    call dreorderf(2*nd, ns, dipstr, dipstrsort, the_plan%isrc)

  ! Allocate output sort buffers
  if (ifpgh .eq. 1) then
    allocate(potsort(nd,ns), gradsort(nd,1), hesssort(nd,1))
  else if (ifpgh .eq. 2) then
    allocate(potsort(nd,ns), gradsort(nd,ns), hesssort(nd,1))
  else if (ifpgh .eq. 3) then
    allocate(potsort(nd,ns), gradsort(nd,ns), hesssort(nd,ns))
  else
    allocate(potsort(nd,1), gradsort(nd,1), hesssort(nd,1))
  endif

  if (ifpghtarg .eq. 1) then
    allocate(pottargsort(nd,max(nt,1)), gradtargsort(nd,1), &
         hesstargsort(nd,1))
  else if (ifpghtarg .eq. 2) then
    allocate(pottargsort(nd,max(nt,1)), gradtargsort(nd,max(nt,1)), &
         hesstargsort(nd,1))
  else if (ifpghtarg .eq. 3) then
    allocate(pottargsort(nd,max(nt,1)), gradtargsort(nd,max(nt,1)), &
         hesstargsort(nd,max(nt,1)))
  else
    allocate(pottargsort(nd,1), gradtargsort(nd,1), hesstargsort(nd,1))
  endif

  ! Zero output buffers
  if (ifpgh .ge. 1) then
    do i = 1, ns
      do idim = 1, nd
        potsort(idim,i) = 0
      enddo
    enddo
  endif
  if (ifpgh .ge. 2) then
    do i = 1, ns
      do idim = 1, nd
        gradsort(idim,i) = 0
      enddo
    enddo
  endif
  if (ifpgh .ge. 3) then
    do i = 1, ns
      do idim = 1, nd
        hesssort(idim,i) = 0
      enddo
    enddo
  endif

  if (ifpghtarg .ge. 1 .and. nt .gt. 0) then
    do i = 1, nt
      do idim = 1, nd
        pottargsort(idim,i) = 0
      enddo
    enddo
  endif
  if (ifpghtarg .ge. 2 .and. nt .gt. 0) then
    do i = 1, nt
      do idim = 1, nd
        gradtargsort(idim,i) = 0
      enddo
    enddo
  endif
  if (ifpghtarg .ge. 3 .and. nt .gt. 0) then
    do i = 1, nt
      do idim = 1, nd
        hesstargsort(idim,i) = 0
      enddo
    enddo
  endif

  lmptmp = (the_plan%nmax + 1) * nd
  allocate(mptemp(lmptmp))

  ! Run FMM with pre-computed tree and lists
  call cfmm2dmain_pre(nd, the_plan%eps, &
       ns, the_plan%sourcesort, &
       ifcharge, chargesort, &
       ifdipole, dipstrsort, &
       nt, the_plan%targsort, nexpc, expc, &
       iaddr, the_plan%rmlexp, mptemp, lmptmp, &
       the_plan%itree, the_plan%ltree, the_plan%iptr, &
       the_plan%ndiv, nlevels, nboxes, iper, &
       the_plan%boxsize, the_plan%rscales, the_plan%tcenters, &
       the_plan%itree(the_plan%iptr(1)), &
       the_plan%isrcse, the_plan%itargse, the_plan%iexpcse, &
       the_plan%nterms, ntj, &
       ifpgh, potsort, gradsort, hesssort, &
       ifpghtarg, pottargsort, gradtargsort, hesstargsort, &
       jexps, scj, ifnear, timeinfo, ier, &
       the_plan%carray, the_plan%ldc, &
       the_plan%mnlist1, the_plan%nlist1s, the_plan%list1, &
       the_plan%mnlist2, the_plan%nlist2s, the_plan%list2, &
       the_plan%mnlist3, the_plan%nlist3s, the_plan%list3, &
       the_plan%mnlist4, the_plan%nlist4s, the_plan%list4, &
       the_plan%nmax, the_plan%mploc_hexp1tmp, &
       the_plan%mploc_jexp2tmp)

  ! Reorder outputs back to original index order
  if (ifpgh .eq. 1) then
    call dreorderi(2*nd, ns, potsort, pot, the_plan%isrc)
  else if (ifpgh .eq. 2) then
    call dreorderi(2*nd, ns, potsort,  pot,  the_plan%isrc)
    call dreorderi(2*nd, ns, gradsort, grad, the_plan%isrc)
  else if (ifpgh .eq. 3) then
    call dreorderi(2*nd, ns, potsort,  pot,  the_plan%isrc)
    call dreorderi(2*nd, ns, gradsort, grad, the_plan%isrc)
    call dreorderi(2*nd, ns, hesssort, hess, the_plan%isrc)
  endif

  if (ifpghtarg .eq. 1) then
    call dreorderi(2*nd, nt, pottargsort, pottarg, the_plan%itarg)
  else if (ifpghtarg .eq. 2) then
    call dreorderi(2*nd, nt, pottargsort,  pottarg,  the_plan%itarg)
    call dreorderi(2*nd, nt, gradtargsort, gradtarg, the_plan%itarg)
  else if (ifpghtarg .eq. 3) then
    call dreorderi(2*nd, nt, pottargsort,   pottarg,   the_plan%itarg)
    call dreorderi(2*nd, nt, gradtargsort,  gradtarg,  the_plan%itarg)
    call dreorderi(2*nd, nt, hesstargsort,  hesstarg,  the_plan%itarg)
  endif

  return
end subroutine cfmm2d_execute_plan


!-----------------------------------------------------------------------
! cfmm2d_destroy_plan
!   Deallocate all plan memory. Safe to call even if no plan is built.
!   Also invalidates any lfmm2d or rfmm2d plan (they share this plan).
!-----------------------------------------------------------------------
subroutine cfmm2d_destroy_plan(ier)
  use cfmm2d_plan_mod
  implicit none

  integer, intent(out) :: ier

  ier = 0

  if (.not. plan_built) return

  if (allocated(the_plan%itree))      deallocate(the_plan%itree)
  if (allocated(the_plan%tcenters))   deallocate(the_plan%tcenters)
  if (allocated(the_plan%boxsize))    deallocate(the_plan%boxsize)
  if (allocated(the_plan%rscales))    deallocate(the_plan%rscales)
  if (allocated(the_plan%nterms))     deallocate(the_plan%nterms)
  if (allocated(the_plan%iaddr_1))    deallocate(the_plan%iaddr_1)
  if (allocated(the_plan%isrc))       deallocate(the_plan%isrc)
  if (allocated(the_plan%itarg))      deallocate(the_plan%itarg)
  if (allocated(the_plan%isrcse))     deallocate(the_plan%isrcse)
  if (allocated(the_plan%itargse))    deallocate(the_plan%itargse)
  if (allocated(the_plan%iexpcse))    deallocate(the_plan%iexpcse)
  if (allocated(the_plan%sourcesort)) deallocate(the_plan%sourcesort)
  if (allocated(the_plan%targsort))   deallocate(the_plan%targsort)
  if (allocated(the_plan%nlist1s))    deallocate(the_plan%nlist1s)
  if (allocated(the_plan%list1))      deallocate(the_plan%list1)
  if (allocated(the_plan%nlist2s))    deallocate(the_plan%nlist2s)
  if (allocated(the_plan%list2))      deallocate(the_plan%list2)
  if (allocated(the_plan%nlist3s))    deallocate(the_plan%nlist3s)
  if (allocated(the_plan%list3))      deallocate(the_plan%list3)
  if (allocated(the_plan%nlist4s))    deallocate(the_plan%nlist4s)
  if (allocated(the_plan%list4))      deallocate(the_plan%list4)
  if (allocated(the_plan%carray))           deallocate(the_plan%carray)
  if (allocated(the_plan%rmlexp))           deallocate(the_plan%rmlexp)
  if (allocated(the_plan%mploc_hexp1tmp))   deallocate(the_plan%mploc_hexp1tmp)
  if (allocated(the_plan%mploc_jexp2tmp))   deallocate(the_plan%mploc_jexp2tmp)

  the_plan%ns        = 0
  the_plan%nt        = 0
  the_plan%nlevels   = 0
  the_plan%nboxes    = 0
  the_plan%ltree     = 0
  the_plan%rmlexp_nd = 0
  the_plan%mploc_nd  = 0
  plan_built = .false.

  return
end subroutine cfmm2d_destroy_plan


!-----------------------------------------------------------------------
! lfmm2d_build_plan
!   Build plan for lfmm2d (complex-charge Laplace FMM).
!   Internally calls cfmm2d_build_plan since the tree depends only
!   on source/target geometry.
!-----------------------------------------------------------------------
subroutine lfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  implicit none
  real*8,  intent(in)  :: eps
  integer, intent(in)  :: ns, nt
  real*8,  intent(in)  :: sources(2,ns), targ(2,max(nt,1))
  integer, intent(out) :: ier
  call cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  return
end subroutine lfmm2d_build_plan


!-----------------------------------------------------------------------
! lfmm2d_execute_plan
!   Run the complex-charge Laplace FMM using a previously built plan.
!   Mirrors the charge/output conversion in lfmm2d.f, then calls
!   cfmm2d_execute_plan and unpacks the Cauchy FMM results.
!
!   INPUT:
!     nd          number of densities
!     ifcharge    1 to include charge interactions
!     charge(nd,ns) complex charge strengths
!     ifdipole    1 to include dipole interactions
!     dipstr(nd,ns) complex dipole strengths
!     dipvec(nd,2,ns) real dipole orientation vectors
!     iper        flag for periodic (pass 0)
!     ifpgh       1/2/3 = pot / pot+grad / pot+grad+hess at sources
!     ifpghtarg   1/2/3 = pot / pot+grad / pot+grad+hess at targets
!   OUTPUT:
!     pot(nd,ns)        potential at sources
!     grad(nd,2,ns)     gradient at sources
!     hess(nd,3,ns)     hessian at sources
!     pottarg(nd,nt)    potential at targets
!     gradtarg(nd,2,nt) gradient at targets
!     hesstarg(nd,3,nt) hessian at targets
!     ier         0 on success
!-----------------------------------------------------------------------
subroutine lfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, dipstr, &
     dipvec, iper, ifpgh, pot, grad, hess, &
     ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  use cfmm2d_plan_mod
  implicit none

  integer,    intent(in)    :: nd, ifcharge, ifdipole, iper
  integer,    intent(in)    :: ifpgh, ifpghtarg
  complex*16, intent(in)    :: charge(nd,*), dipstr(nd,*)
  real*8,     intent(in)    :: dipvec(nd,2,*)
  complex*16, intent(inout) :: pot(nd,*), grad(nd,2,*), hess(nd,3,*)
  complex*16, intent(inout) :: pottarg(nd,*), gradtarg(nd,2,*), &
                                hesstarg(nd,3,*)
  integer,    intent(out)   :: ier

  integer :: nd2, ns, nt, i, j
  complex*16, parameter :: eye = (0.0d0, 1.0d0)
  complex*16 :: ztmp

  complex*16, allocatable :: charge1(:,:,:), dipstr1(:,:,:)
  complex*16, allocatable :: pot1(:,:,:), grad1(:,:,:), hess1(:,:,:)
  complex*16, allocatable :: pottarg1(:,:,:), gradtarg1(:,:,:)
  complex*16, allocatable :: hesstarg1(:,:,:)

  ier = 0

  if (.not. plan_built) then
    ier = 4
    return
  endif

  ns  = the_plan%ns
  nt  = the_plan%nt
  nd2 = 2 * nd

  ! Allocate charge arrays split into real and imaginary parts
  if (ifcharge .eq. 1) then
    allocate(charge1(2,nd,ns))
  else
    allocate(charge1(2,nd,1))
  endif

  if (ifdipole .eq. 1) then
    allocate(dipstr1(2,nd,ns))
  else
    allocate(dipstr1(2,nd,1))
  endif

  ! Allocate output arrays
  if (ifpgh .eq. 1) then
    allocate(pot1(2,nd,ns),  grad1(2,nd,1),  hess1(2,nd,1))
  else if (ifpgh .eq. 2) then
    allocate(pot1(2,nd,ns),  grad1(2,nd,ns), hess1(2,nd,1))
  else if (ifpgh .eq. 3) then
    allocate(pot1(2,nd,ns),  grad1(2,nd,ns), hess1(2,nd,ns))
  else
    allocate(pot1(2,nd,1),   grad1(2,nd,1),  hess1(2,nd,1))
  endif

  if (ifpghtarg .eq. 1) then
    allocate(pottarg1(2,nd,max(nt,1)), gradtarg1(2,nd,1), &
         hesstarg1(2,nd,1))
  else if (ifpghtarg .eq. 2) then
    allocate(pottarg1(2,nd,max(nt,1)), gradtarg1(2,nd,max(nt,1)), &
         hesstarg1(2,nd,1))
  else if (ifpghtarg .eq. 3) then
    allocate(pottarg1(2,nd,max(nt,1)), gradtarg1(2,nd,max(nt,1)), &
         hesstarg1(2,nd,max(nt,1)))
  else
    allocate(pottarg1(2,nd,1), gradtarg1(2,nd,1), hesstarg1(2,nd,1))
  endif

  ! Split charges: charge1(1,j,i)=Re(charge), charge1(2,j,i)=Im(charge)
  if (ifcharge .eq. 1) then
    do i = 1, ns
      do j = 1, nd
        charge1(1,j,i) = dble(charge(j,i))
        charge1(2,j,i) = imag(charge(j,i))
      enddo
    enddo
  endif

  ! Convert dipoles: dipstr1 = dipstr * -(dipvec_x + i*dipvec_y)
  if (ifdipole .eq. 1) then
    do i = 1, ns
      do j = 1, nd
        ztmp = -(dipvec(j,1,i) + eye*dipvec(j,2,i))
        dipstr1(1,j,i) = dble(dipstr(j,i))  * ztmp
        dipstr1(2,j,i) = imag(dipstr(j,i)) * ztmp
      enddo
    enddo
  endif

  ! Run Cauchy FMM via the plan (nd2 = 2*nd)
  call cfmm2d_execute_plan(nd2, ifcharge, charge1, ifdipole, dipstr1, &
       iper, ifpgh, pot1, grad1, hess1, &
       ifpghtarg, pottarg1, gradtarg1, hesstarg1, ier)
  if (ier .ne. 0) return

  ! Unpack potential: pot = Re(pot1_re) + i*Re(pot1_im)
  if (ifpgh .ge. 1) then
    do i = 1, ns
      do j = 1, nd
        pot(j,i) = dble(pot1(1,j,i)) + eye*dble(pot1(2,j,i))
      enddo
    enddo
  endif

  ! Unpack gradient: d/dx = Re(d/dz) + i*Re(d/dz_im);  d/dy = -Im(...)
  if (ifpgh .ge. 2) then
    do i = 1, ns
      do j = 1, nd
        grad(j,1,i) =  dble(grad1(1,j,i)) + eye*dble(grad1(2,j,i))
        grad(j,2,i) = -imag(grad1(1,j,i)) - eye*imag(grad1(2,j,i))
      enddo
    enddo
  endif
  if (ifpgh .eq. 3) then
    do i = 1, ns
      do j = 1, nd
        hess(j,1,i) =  dble(hess1(1,j,i)) + eye*dble(hess1(2,j,i))
        hess(j,2,i) = -imag(hess1(1,j,i)) - eye*imag(hess1(2,j,i))
        hess(j,3,i) = -hess(j,1,i)
      enddo
    enddo
  endif

  if (ifpghtarg .ge. 1 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        pottarg(j,i) = dble(pottarg1(1,j,i)) + eye*dble(pottarg1(2,j,i))
      enddo
    enddo
  endif
  if (ifpghtarg .ge. 2 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        gradtarg(j,1,i) =  dble(gradtarg1(1,j,i)) + &
                            eye*dble(gradtarg1(2,j,i))
        gradtarg(j,2,i) = -imag(gradtarg1(1,j,i)) - &
                            eye*imag(gradtarg1(2,j,i))
      enddo
    enddo
  endif
  if (ifpghtarg .eq. 3 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        hesstarg(j,1,i) =  dble(hesstarg1(1,j,i)) + &
                            eye*dble(hesstarg1(2,j,i))
        hesstarg(j,2,i) = -imag(hesstarg1(1,j,i)) - &
                            eye*imag(hesstarg1(2,j,i))
        hesstarg(j,3,i) = -hesstarg(j,1,i)
      enddo
    enddo
  endif

  return
end subroutine lfmm2d_execute_plan


!-----------------------------------------------------------------------
! rfmm2d_build_plan
!   Build plan for rfmm2d (real-charge Laplace FMM).
!   Delegates to cfmm2d_build_plan.
!-----------------------------------------------------------------------
subroutine rfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  implicit none
  real*8,  intent(in)  :: eps
  integer, intent(in)  :: ns, nt
  real*8,  intent(in)  :: sources(2,ns), targ(2,max(nt,1))
  integer, intent(out) :: ier
  call cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  return
end subroutine rfmm2d_build_plan


!-----------------------------------------------------------------------
! rfmm2d_execute_plan
!   Run the real-charge Laplace FMM using a previously built plan.
!   Mirrors the charge/output conversion in rfmm2d.f.
!
!   INPUT:
!     nd          number of densities
!     ifcharge    1 to include charge interactions
!     charge(nd,ns) real charge strengths
!     ifdipole    1 to include dipole interactions
!     dipstr(nd,ns) real dipole strengths
!     dipvec(nd,2,ns) real dipole orientation vectors
!     iper        flag for periodic (pass 0)
!     ifpgh       1/2/3 = pot / pot+grad / pot+grad+hess at sources
!     ifpghtarg   1/2/3 = pot / pot+grad / pot+grad+hess at targets
!   OUTPUT:
!     pot(nd,ns)        potential at sources
!     grad(nd,2,ns)     gradient at sources
!     hess(nd,3,ns)     hessian at sources
!     pottarg(nd,nt)    potential at targets
!     gradtarg(nd,2,nt) gradient at targets
!     hesstarg(nd,3,nt) hessian at targets
!     ier         0 on success
!-----------------------------------------------------------------------
subroutine rfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, dipstr, &
     dipvec, iper, ifpgh, pot, grad, hess, &
     ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  use cfmm2d_plan_mod
  implicit none

  integer, intent(in)  :: nd, ifcharge, ifdipole, iper
  integer, intent(in)  :: ifpgh, ifpghtarg
  real*8,  intent(in)  :: charge(nd,*), dipstr(nd,*)
  real*8,  intent(in)  :: dipvec(nd,2,*)
  real*8,  intent(inout) :: pot(nd,*), grad(nd,2,*), hess(nd,3,*)
  real*8,  intent(inout) :: pottarg(nd,*), gradtarg(nd,2,*), &
                             hesstarg(nd,3,*)
  integer, intent(out) :: ier

  integer :: ns, nt, i, j
  complex*16, parameter :: eye = (0.0d0, 1.0d0)
  complex*16 :: ztmp

  complex*16, allocatable :: charge1(:,:), dipstr1(:,:)
  complex*16, allocatable :: pot1(:,:), grad1(:,:), hess1(:,:)
  complex*16, allocatable :: pottarg1(:,:), gradtarg1(:,:), hesstarg1(:,:)

  ier = 0

  if (.not. plan_built) then
    ier = 4
    return
  endif

  ns = the_plan%ns
  nt = the_plan%nt

  ! Allocate charge arrays as complex
  if (ifcharge .eq. 1) then
    allocate(charge1(nd,ns))
  else
    allocate(charge1(nd,1))
  endif

  if (ifdipole .eq. 1) then
    allocate(dipstr1(nd,ns))
  else
    allocate(dipstr1(nd,1))
  endif

  ! Allocate output arrays
  if (ifpgh .eq. 1) then
    allocate(pot1(nd,ns),  grad1(nd,1),  hess1(nd,1))
  else if (ifpgh .eq. 2) then
    allocate(pot1(nd,ns),  grad1(nd,ns), hess1(nd,1))
  else if (ifpgh .eq. 3) then
    allocate(pot1(nd,ns),  grad1(nd,ns), hess1(nd,ns))
  else
    allocate(pot1(nd,1),   grad1(nd,1),  hess1(nd,1))
  endif

  if (ifpghtarg .eq. 1) then
    allocate(pottarg1(nd,max(nt,1)), gradtarg1(nd,1), hesstarg1(nd,1))
  else if (ifpghtarg .eq. 2) then
    allocate(pottarg1(nd,max(nt,1)), gradtarg1(nd,max(nt,1)), &
         hesstarg1(nd,1))
  else if (ifpghtarg .eq. 3) then
    allocate(pottarg1(nd,max(nt,1)), gradtarg1(nd,max(nt,1)), &
         hesstarg1(nd,max(nt,1)))
  else
    allocate(pottarg1(nd,1), gradtarg1(nd,1), hesstarg1(nd,1))
  endif

  ! Convert real charges to complex
  if (ifcharge .eq. 1) then
    do i = 1, ns
      do j = 1, nd
        charge1(j,i) = charge(j,i)
      enddo
    enddo
  endif

  ! Convert real dipoles: dipstr1 = dipstr * -(dipvec_x + i*dipvec_y)
  if (ifdipole .eq. 1) then
    do i = 1, ns
      do j = 1, nd
        ztmp = -(dipvec(j,1,i) + eye*dipvec(j,2,i))
        dipstr1(j,i) = dipstr(j,i) * ztmp
      enddo
    enddo
  endif

  ! Run Cauchy FMM via the plan
  call cfmm2d_execute_plan(nd, ifcharge, charge1, ifdipole, dipstr1, &
       iper, ifpgh, pot1, grad1, hess1, &
       ifpghtarg, pottarg1, gradtarg1, hesstarg1, ier)
  if (ier .ne. 0) return

  ! Unpack real output
  if (ifpgh .ge. 1) then
    do i = 1, ns
      do j = 1, nd
        pot(j,i) = dble(pot1(j,i))
      enddo
    enddo
  endif
  if (ifpgh .ge. 2) then
    do i = 1, ns
      do j = 1, nd
        grad(j,1,i) =  dble(grad1(j,i))
        grad(j,2,i) = -imag(grad1(j,i))
      enddo
    enddo
  endif
  if (ifpgh .eq. 3) then
    do i = 1, ns
      do j = 1, nd
        hess(j,1,i) =  dble(hess1(j,i))
        hess(j,2,i) = -imag(hess1(j,i))
        hess(j,3,i) = -hess(j,1,i)
      enddo
    enddo
  endif

  if (ifpghtarg .ge. 1 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        pottarg(j,i) = dble(pottarg1(j,i))
      enddo
    enddo
  endif
  if (ifpghtarg .ge. 2 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        gradtarg(j,1,i) =  dble(gradtarg1(j,i))
        gradtarg(j,2,i) = -imag(gradtarg1(j,i))
      enddo
    enddo
  endif
  if (ifpghtarg .eq. 3 .and. nt .gt. 0) then
    do i = 1, nt
      do j = 1, nd
        hesstarg(j,1,i) =  dble(hesstarg1(j,i))
        hesstarg(j,2,i) = -imag(hesstarg1(j,i))
        hesstarg(j,3,i) = -hesstarg(j,1,i)
      enddo
    enddo
  endif

  return
end subroutine rfmm2d_execute_plan
