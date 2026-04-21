! test_plan.f90 -- verify cfmm2d/lfmm2d/rfmm2d plan API
!
! Checks that plan-based execute gives results matching the direct
! (non-plan) FMM call to within a tight tolerance.
program test_plan
  implicit none

  integer, parameter :: ns = 500, nt = 200, nd = 2
  real*8,  parameter :: eps = 1.0d-9
  integer :: iper

  real*8  :: sources(2,ns), targ(2,nt)
  ! cfmm2d
  complex*16 :: charge(nd,ns), dipstr(nd,ns)
  complex*16 :: pot_ref(nd,ns),    pot_plan(nd,ns)
  complex*16 :: grad_ref(nd,ns),   grad_plan(nd,ns)
  complex*16 :: pottarg_ref(nd,nt), pottarg_plan(nd,nt)
  ! lfmm2d
  complex*16 :: lcharge(nd,ns), ldipstr(nd,ns)
  real*8     :: ldipvec(nd,2,ns)
  complex*16 :: lpot_ref(nd,ns),    lpot_plan(nd,ns)
  complex*16 :: lgrad_ref(nd,2,ns), lgrad_plan(nd,2,ns)
  complex*16 :: lpottarg_ref(nd,nt), lpottarg_plan(nd,nt)
  ! rfmm2d
  real*8     :: rcharge(nd,ns), rdipstr(nd,ns)
  real*8     :: rdipvec(nd,2,ns)
  real*8     :: rpot_ref(nd,ns),    rpot_plan(nd,ns)
  real*8     :: rgrad_ref(nd,2,ns), rgrad_plan(nd,2,ns)
  real*8     :: rpottarg_ref(nd,nt), rpottarg_plan(nd,nt)

  real*8    :: err, nrm
  integer   :: ier, i, j, ipass, ifcharge, ifdipole, ifpgh, ifpghtarg
  real*8    :: dlaran
  integer   :: seed(4)

  iper = 0
  seed = (/1,2,3,4/)

  ! Random geometry
  do i = 1, ns
    sources(1,i) = dlaran(seed)
    sources(2,i) = dlaran(seed)
  enddo
  do i = 1, nt
    targ(1,i) = dlaran(seed)
    targ(2,i) = dlaran(seed)
  enddo

  ! Random charges
  do i = 1, ns
    do j = 1, nd
      charge(j,i)  = dcmplx(dlaran(seed), dlaran(seed))
      dipstr(j,i)  = dcmplx(dlaran(seed), dlaran(seed))
      lcharge(j,i) = dcmplx(dlaran(seed), dlaran(seed))
      ldipstr(j,i) = dcmplx(dlaran(seed), dlaran(seed))
      rcharge(j,i) = dlaran(seed)
      rdipstr(j,i) = dlaran(seed)
      ldipvec(j,1,i) = dlaran(seed)
      ldipvec(j,2,i) = dlaran(seed)
      rdipvec(j,1,i) = ldipvec(j,1,i)
      rdipvec(j,2,i) = ldipvec(j,2,i)
    enddo
  enddo

  ipass = 1

  !-----------------------------------------------------------------
  ! cfmm2d: build plan, compare to direct call
  !-----------------------------------------------------------------
  call cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: cfmm2d_build_plan returned ier=', ier
    stop 1
  endif

  ifcharge = 1 ; ifdipole = 1
  ifpgh = 2    ; ifpghtarg = 1

  ! Reference: direct cfmm2d call
  call cfmm2d(nd, eps, ns, sources, ifcharge, charge, &
       ifdipole, dipstr, iper, ifpgh, pot_ref, grad_ref, pot_ref, &
       nt, targ, ifpghtarg, pottarg_ref, grad_ref, pot_ref, ier)

  ! Plan-based call
  call cfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, dipstr, &
       iper, ifpgh, pot_plan, grad_plan, pot_plan, &
       ifpghtarg, pottarg_plan, grad_plan, pot_plan, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: cfmm2d_execute_plan returned ier=', ier
    stop 1
  endif

  err = 0; nrm = 0
  do i = 1, ns
    do j = 1, nd
      err = err + abs(pot_plan(j,i)  - pot_ref(j,i))**2
      nrm = nrm + abs(pot_ref(j,i))**2
    enddo
  enddo
  err = sqrt(err/max(nrm, 1.0d-30))
  if (err .gt. 1.0d-13) then
    print *, 'FAIL cfmm2d plan pot error =', err
    ipass = 0
  else
    print *, 'PASS cfmm2d plan pot error =', err
  endif

  call cfmm2d_destroy_plan(ier)

  !-----------------------------------------------------------------
  ! lfmm2d: build plan, compare to direct call
  !-----------------------------------------------------------------
  call lfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: lfmm2d_build_plan returned ier=', ier
    stop 1
  endif

  ! Reference: direct lfmm2d call
  call lfmm2d(nd, eps, ns, sources, ifcharge, lcharge, &
       ifdipole, ldipstr, ldipvec, iper, ifpgh, lpot_ref, lgrad_ref, &
       lpot_ref, nt, targ, ifpghtarg, lpottarg_ref, lgrad_ref, &
       lpot_ref, ier)

  ! Plan-based
  call lfmm2d_execute_plan(nd, ifcharge, lcharge, ifdipole, ldipstr, &
       ldipvec, iper, ifpgh, lpot_plan, lgrad_plan, lpot_plan, &
       ifpghtarg, lpottarg_plan, lgrad_plan, lpot_plan, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: lfmm2d_execute_plan returned ier=', ier
    stop 1
  endif

  err = 0; nrm = 0
  do i = 1, ns
    do j = 1, nd
      err = err + abs(lpot_plan(j,i) - lpot_ref(j,i))**2
      nrm = nrm + abs(lpot_ref(j,i))**2
    enddo
  enddo
  err = sqrt(err/max(nrm, 1.0d-30))
  if (err .gt. 1.0d-13) then
    print *, 'FAIL lfmm2d plan pot error =', err
    ipass = 0
  else
    print *, 'PASS lfmm2d plan pot error =', err
  endif

  call cfmm2d_destroy_plan(ier)

  !-----------------------------------------------------------------
  ! rfmm2d: build plan, compare to direct call
  !-----------------------------------------------------------------
  call rfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: rfmm2d_build_plan returned ier=', ier
    stop 1
  endif

  ! Reference: direct rfmm2d call
  call rfmm2d(nd, eps, ns, sources, ifcharge, rcharge, &
       ifdipole, rdipstr, rdipvec, iper, ifpgh, rpot_ref, rgrad_ref, &
       rpot_ref, nt, targ, ifpghtarg, rpottarg_ref, rgrad_ref, &
       rpot_ref, ier)

  ! Plan-based
  call rfmm2d_execute_plan(nd, ifcharge, rcharge, ifdipole, rdipstr, &
       rdipvec, iper, ifpgh, rpot_plan, rgrad_plan, rpot_plan, &
       ifpghtarg, rpottarg_plan, rgrad_plan, rpot_plan, ier)
  if (ier .ne. 0) then
    print *, 'FAIL: rfmm2d_execute_plan returned ier=', ier
    stop 1
  endif

  err = 0; nrm = 0
  do i = 1, ns
    do j = 1, nd
      err = err + (rpot_plan(j,i) - rpot_ref(j,i))**2
      nrm = nrm + rpot_ref(j,i)**2
    enddo
  enddo
  err = sqrt(err/max(nrm, 1.0d-30))
  if (err .gt. 1.0d-13) then
    print *, 'FAIL rfmm2d plan pot error =', err
    ipass = 0
  else
    print *, 'PASS rfmm2d plan pot error =', err
  endif

  call cfmm2d_destroy_plan(ier)

  if (ipass .eq. 1) then
    print *, 'All plan tests passed.'
    stop 0
  else
    print *, 'One or more plan tests FAILED.'
    stop 1
  endif

end program test_plan
