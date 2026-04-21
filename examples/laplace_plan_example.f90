! laplace_plan_example.f90
!
! Demonstrates the speedup from the Laplace FMM plan API.
!
! Geometry: 1,000,000 sources and 1,000,000 targets uniformly random in
! [0,1]^2, nd=1, charges only, potentials requested at sources and targets.
!
! The example runs three timed measurements:
!   (1) cfmm2d -- standard call (tree build + FMM traversal)
!   (2) cfmm2d_build_plan -- tree build and precomputation only
!   (3) cfmm2d_execute_plan -- FMM traversal using the cached plan
!
! Speedup is reported as (1) / (3), showing how much time is saved per
! execute call once the plan is built.

program laplace_plan_example
  implicit none

  integer, parameter :: ns = 1000000
  integer, parameter :: nt = 1000000
  integer, parameter :: nd = 1

  real*8,     allocatable :: sources(:,:), targ(:,:)
  complex*16, allocatable :: charge(:,:)
  complex*16, allocatable :: pot(:,:), pottarg(:,:)
  complex*16, allocatable :: pot2(:,:), pottarg2(:,:)

  real*8  :: eps
  integer :: ier, i, ifcharge, ifdipole, ifpgh, ifpghtarg, iper
  real*8  :: t0, t1, t_direct, t_build, t_exec

  real*8  :: omp_get_wtime
  real*8  :: dlaran
  integer :: seed(4)

  seed = (/1, 2, 3, 4/)
  eps  = 1.0d-6
  iper = 0

  allocate(sources(2,ns), targ(2,nt))
  allocate(charge(nd,ns))
  allocate(pot(nd,ns), pottarg(nd,nt))
  allocate(pot2(nd,ns), pottarg2(nd,nt))

  ! Random geometry in [0,1]^2
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
    charge(1,i) = dcmplx(dlaran(seed), dlaran(seed))
  enddo

  ifcharge  = 1
  ifdipole  = 0
  ifpgh     = 1   ! potential at sources
  ifpghtarg = 1   ! potential at targets

  !------------------------------------------------------------------
  ! (1) Standard cfmm2d call (tree build + FMM traversal combined)
  !------------------------------------------------------------------
  t0 = omp_get_wtime()
  call cfmm2d(nd, eps, ns, sources, ifcharge, charge, &
       ifdipole, charge, iper, ifpgh, pot, pot, pot, &
       nt, targ, ifpghtarg, pottarg, pottarg, pottarg, ier)
  t1 = omp_get_wtime()
  t_direct = t1 - t0

  write(*,'(A)') ''
  write(*,'(A)') '  Laplace FMM plan speedup example'
  write(*,'(A)') '  ================================='
  write(*,'(A,I10)') '  Sources : ', ns
  write(*,'(A,I10)') '  Targets : ', nt
  write(*,'(A,ES10.2)') '  Epsilon : ', eps
  write(*,'(A)') ''
  write(*,'(A,F8.3,A)') '  cfmm2d (direct)    : ', t_direct,  ' s'

  !------------------------------------------------------------------
  ! (2) Build plan (tree + interaction lists; one-time cost)
  !------------------------------------------------------------------
  t0 = omp_get_wtime()
  call cfmm2d_build_plan(eps, ns, sources, nt, targ, ier)
  t1 = omp_get_wtime()
  t_build = t1 - t0
  if (ier .ne. 0) then
    write(*,*) 'cfmm2d_build_plan failed, ier=', ier
    stop 1
  endif
  write(*,'(A,F8.3,A)') '  cfmm2d_build_plan  : ', t_build,   ' s'

  !------------------------------------------------------------------
  ! (3) Execute plan (FMM traversal only; repeated call cost)
  !------------------------------------------------------------------
  t0 = omp_get_wtime()
  call cfmm2d_execute_plan(nd, ifcharge, charge, ifdipole, charge, &
       iper, ifpgh, pot2, pot2, pot2, &
       ifpghtarg, pottarg2, pottarg2, pottarg2, ier)
  t1 = omp_get_wtime()
  t_exec = t1 - t0
  if (ier .ne. 0) then
    write(*,*) 'cfmm2d_execute_plan failed, ier=', ier
    stop 1
  endif
  write(*,'(A,F8.3,A)') '  cfmm2d_execute_plan: ', t_exec,    ' s'

  call cfmm2d_destroy_plan(ier)

  write(*,'(A)') ''
  write(*,'(A,F6.2,A)') '  Speedup (direct / execute_plan)   : ', &
       t_direct / t_exec, 'x'
  write(*,'(A,F6.2,A)') '  Speedup (direct / build+execute)  : ', &
       t_direct / (t_build + t_exec), 'x'
  write(*,'(A)') ''

  ! Verify results match
  call check_error(nd, ns, pot, pot2, 'source pot')
  call check_error(nd, nt, pottarg, pottarg2, 'target pot')

end program laplace_plan_example


subroutine check_error(nd, n, a, b, label)
  implicit none
  integer,    intent(in) :: nd, n
  complex*16, intent(in) :: a(nd,n), b(nd,n)
  character(*), intent(in) :: label
  real*8 :: err, nrm
  integer :: i, j
  err = 0.0d0
  nrm = 0.0d0
  do i = 1, n
    do j = 1, nd
      err = err + abs(a(j,i) - b(j,i))**2
      nrm = nrm + abs(a(j,i))**2
    enddo
  enddo
  err = sqrt(err / max(nrm, 1.0d-30))
  write(*,'(A,A,A,ES10.3)') '  Relative l2 error (', label, '): ', err
end subroutine check_error
