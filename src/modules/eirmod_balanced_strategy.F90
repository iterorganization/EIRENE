module eirmod_balanced_strategy
! An optimized workload distribution to reduce MPI waiting time and maximize
! parallel efficiency.
!
! Main ingredients:
!  - Eirmod_calstr_buffered is used to allow non-blocking reductions in calstr.
!  - If we have N processors calculating stratum(k) then we do NOT divide these
!    particles evenly between the processors.
!  - We measure the execution time of the particle loop, calstr, and stratum
!    postprocessing and use this data to optimize how many particles should
!    each PEs process.
!
! There value of the T_EPSILON_FRACTION parameter can influence the work
! distribution. Lower values increase the time for optimization. Higher values
! will give less optimal results.
!
  use eirmod_precision
  use eirmod_comprt &
   , only: iunout
  use eirmod_comsou &
   , only: npts, nstrai
  use eirmod_parmmod &
   , only: nstra

  implicit none

  private

  !> Init_balanced_strategy should be called first, afterwards in every
  !> iteration we optimize
  public init_balanced_strategy
  public opt_balanced_strategy, simple_balanced_strategy

  public allocate_balanced_strategy, deallocate_balanced_strategy

  !> The timings subroutines should be called from mcarlo
  !> The workload distribution will be optimized based on the data
  !> collected by these subroutines
  public time_calstr, time_postproc, time_particles

  !> Throughput for all PEs over all strata (allocated only on PE 0)
  real(kind=dp), allocatable, dimension(:), public :: throughput_pe_all

  !> Overhead time for each PE, used only for printing work distribution table
  !> (allocated only on PE 0)
  real(kind=dp), allocatable, dimension(:), public :: t_overhead

  !> In opt_balanced_strategy we estimate the ideal working time t_ideal,
  !> and then iteratively distribute it among PEs. In every iteration at least
  !> T_EPSILON_FRACTION of the ideal time is distributed. Currently it is 5%.
  real(kind=dp), parameter :: T_EPSILON_FRACTION = 0.05_dp

  ! Arrays to store timing information, all local to the current PE:

  !> The n_particles_loc(k) stores how many particles have been
  !> processed on stratum(k) by the local PE.
  !> It is a commulative quantity, in every iteration we increase
  !> n_particles_loc(k) by the number of particles that were processed.
  integer, allocatable, dimension(:) :: n_particles_loc

  !> t_particles_loc(k) stores the CPU time (in seconds) that was used to
  !> process the particles n_particles_loc(k)
  !> This way we can calculate an average processing time over many iterations
  !> (see calculate_throughput for details).
  real(kind=dp), allocatable, dimension(:) :: t_particles_loc

  !> How many times calstr was called for each stratum
  integer, allocatable, dimension(:) :: n_calstr_loc

  !> If calstr is called for stratum k, then t_calstr(k) is the execution time
  !> of calstr. This is already a weighted average over iterations of eirene
  real(kind=dp), allocatable, dimension(:) :: t_calstr_loc

  !> calstr time for each stratum, averaged over all the iterations
  real(kind=dp), allocatable, dimension(:) :: t_calstr_strat_avg

  !> calstr time for each PE individually
  real(kind=dp), allocatable, dimension(:) :: t_calstr_pe

  !> The number of all the postprocessing steps so far is
  integer :: n_postproc_executed

  !> t_postproc(k) is the postprocessing time (after calstr) for stratum k
  real(kind=dp), allocatable, dimension(:) :: t_postproc_loc

  !> We cannot acces these from CPES (circular reference) so we store them here
  integer :: my_pe, nprs
  integer :: nsteff

  !> Whether to print timing information that was used for the optimization
  logical, parameter :: print_timing_info = .false.

  contains

  subroutine allocate_balanced_strategy
    allocate(n_particles_loc(nstra))
    n_particles_loc = 0
    allocate(t_particles_loc(nstra))
    t_particles_loc = 0
    allocate(n_calstr_loc(nstra))
    n_calstr_loc = 0
    allocate(t_calstr_loc(nstra))
    t_calstr_loc = 0
    allocate(t_postproc_loc(nstra))
    t_postproc_loc = 0
  end subroutine

  subroutine deallocate_balanced_strategy
    use eirmod_calstr_buffered &
     , only: deallocate_calstr_buffer
    deallocate(t_particles_loc)
    deallocate(t_calstr_loc)
    deallocate(t_postproc_loc)
    deallocate(n_particles_loc)
    deallocate(n_calstr_loc)
    if(allocated(t_calstr_strat_avg)) deallocate(t_calstr_strat_avg)
    if(allocated(t_calstr_pe)) deallocate(t_calstr_pe)
    call deallocate_calstr_buffer
    if(allocated(throughput_pe_all)) deallocate(throughput_pe_all)
    if(allocated(t_overhead)) deallocate(t_overhead)
  end subroutine

  subroutine init_balanced_strategy(ierror)
  ! Check if we have correct MPI version, and initializes calstr_buffer.
  ! This subroutine does not define a parallelization strategy. After calling
  ! this subroutine, please use any other strategy (preferably STRATEGY_APCAS)
  ! to initialize the work distribution.
    use eirmod_calstr_buffered
    use eirmod_mpi
    integer, intent(out) :: ierror !< 0 = success, any other values = error
    integer ierr, ierrr
    integer mpi_major, mpi_minor
    external eirene_masage, eirene_exit_own
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external mpi_allreduce
#endif

    ! We cannot access my_pe and nprs from eirmod_cpes, because that would be a
    ! a circular reference. We store locally instead.
    call mpi_comm_rank(MPI_COMM_WORLD, my_pe, ierr)
    if (ierr /= mpi_success) then
      call eirene_masage('ERROR IN init_balanced_strategy at mpi_comm_rank.')
      call eirene_exit_own(1)
    end if
    call mpi_comm_size(MPI_COMM_WORLD, nprs, ierr)
    if (nprs.eq.1) then
      if (my_pe.eq.0) then
        write(iunout,*) &
         'Warning: we are not using balanced strategy in serial mode'
      endif
      ierror = 1
      return
    endif
    nsteff = count(npts(1:nstrai) > 0)
    ierror = 0
    call mpi_get_version(mpi_major, mpi_minor, ierr)
    if (mpi_major < 3) then
      write(iunout,*) 'Error, balanced strategy requires MPI 3'
      ierror = 1
    else
      ! MPI_VERSION >= 3
      if (my_pe.eq.0) write(iunout,*) 'Initializing balanced strategy'
      call allocate_calstr_buffer(ierr)
      call mpi_allreduce(ierr, ierror, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierrr)
      if (ierror.ne.0) then
        if(my_pe.eq.0) then
          write(iunout,*) &
           'Error, not enough memory to allocate calstr buffer.'
        endif
      endif
    endif
  end subroutine

  subroutine opt_balanced_strategy(nparts_loc, npestr, stratum_leader, &
              procforstra)
  ! Based on the measurements by time_particles, time_calstr, time_postproc we
  ! optimize the work distribution table
  !
  ! We calculate the throughput (particles/sec) for each PE and strata. We use
  ! this to estimate how long it will take to process the particle loop.
  !
  ! First we select the stratum leaders, and let them work as long as they can
  ! on their stratum (to minimize the number of PEs per stratum, and this way
  ! minimize the communication). Then we distribute the rest of the work to the
  ! remaining PEs. Within a stratum, we aim to finish all the particle loops
  ! before the stratum leader finishes its particle loop, to minimize waiting
  ! time in calstr (stratum leader has to wait there before continuing with
  ! stratum postprocessing).
  ! We consider the time needed to execute calstr for all PEs, and the
  ! postprocessing for the stratum leaders. This ensures a quasi-optimal
  ! workload distribution.

    use eirmod_mpi
    !> We will define the following parameters:
    integer, intent(inout), dimension(:) :: nparts_loc
    integer, intent(inout), dimension(:) :: npestr
    integer, intent(inout), dimension(:) :: stratum_leader
    logical, intent(inout), dimension(:,0:) :: procforstra

    real(kind=dp), dimension(nstrai) :: t_postproc
    real(kind=dp), dimension(nstrai) :: throughput !< aggregate throughput per stratum
    real(kind=dp), dimension(nstrai) :: throughput_strat_avg !< average per stratum
    real(kind=dp), dimension(nstrai) :: t_ideal
    integer, dimension(nstrai*nprs) :: nparts_loc_all
    real(kind=dp), dimension(0:nprs) :: t_pe
    real(kind=dp) :: t_ideal_tot, t_epsilon, t_calstr_allavg
    real(kind=dp) :: throughput_avg, time, t_target
    integer, dimension(nstrai) :: npts_remaining, n_epsilon
    integer :: i, k, idx, n, ierr
    character*13 hlp_frm
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external :: mpi_bcast, mpi_scatter
#endif

    if (nprs == 1) return ! nothing to optimize for serial mode

    if(my_pe==0) then
      write(iunout,*) 'Optimizing workload distribution'
    endif

    ! For each PE and stratum we estimate the throghput and the time to execute
    ! calstr. Additionally, we calculate the average postprocessing time for
    ! each stratum. All PEs have to call the next three subroutines, because of
    ! the collective MPI calls.
    call calculate_calstr_time(t_calstr_strat_avg, t_calstr_allavg, t_calstr_pe)
    call calculate_postproc_time(t_postproc)
    call calculate_throughput(throughput_pe_all, throughput, throughput_avg, &
                              throughput_strat_avg)
    if (my_pe == 0) then
      ! The rest of the calculation is performed on PE 0, because only PE 0 has
      ! all the data
      call calculate_ideal_time(throughput, t_calstr_strat_avg, t_postproc, &
                                t_ideal, t_ideal_tot)
      t_pe = 0 ! t_pe(i) stores how much time PE i has already worked
      !> we give away at least t_epsilon time in each step
      t_epsilon = t_ideal_tot * T_EPSILON_FRACTION / nprs
      ! Calstr has an overhead, we calculate how many particles could be
      ! processed in the overhead time.
      ! But t_calstr_allavg might be distorted if there is load imbalance
      ! therefore we choose the smaller of t_epsilon and t_calstr_allavg
      n_epsilon = int(min(t_epsilon, t_calstr_allavg) * throughput_strat_avg)
      if (print_timing_info) then
        write(hlp_frm,'(a,i4,a)') '(1x,a,',size(n_epsilon),'i5)'
        write(iunout,hlp_frm) 'n_epsilon', n_epsilon
        write(iunout,*) 't_calstr_epsilon ', t_calstr_allavg * throughput
      endif
      ! Initialize output arrays
      !> the number of particles that needs to be distributed:
      npts_remaining = npts(1:nstrai)
      npestr = 0         !< number of PEs per stratum
      !> The most important output quantity will be calculated in nparts_loc_all
      !> nparts_loc_all(i*nstrai + k) is the number of particles from stratum k
      !> processed by PE i
      nparts_loc_all = 0
      procforstra = .false. !< this is not used, but we can calculate anyways

      if(.not.allocated(t_overhead)) then
        allocate(t_overhead(0:nprs-1))
      endif
      t_overhead = 0.0_dp !< communication and postproc. overhead for each PE
      ! t_overhead is used only for printing overhead information in
      ! print_work_distribution_table

      ! First select the leaders and assign work to them
      do k=1,nstrai
        if (npts(k) > 0) then
          i = get_next_pe(t_pe)
          stratum_leader(k) = i
          idx = i*nstrai + k
          ! we try to assign all the available time
          time = t_ideal_tot - t_pe(i)
          call assign_work(time, t_pe(i), nparts_loc_all(idx), &
                   throughput_pe_all(idx), npts_remaining(k), n_epsilon(k))
        endif
      end do
      ! Now we distribute the rest of the work
      n = 0 !< calculate how many steps were used to distribute the work
      do k = 1, nstrai
        if (npts(k) <= 0) cycle
        ! we try to distribute the work so that all PE who calculate a stratum
        ! finishes by t_target
        t_target = t_pe(stratum_leader(k))
        do while (npts_remaining(k) > 0)
          n = n + 1
          i = get_next_pe(t_pe)
          if (t_target - t_pe(i) < t_epsilon ) then
            time = t_epsilon
          else
            time = t_target - t_pe(i)
          endif
          idx = i*nstrai + k
          ! assign_work is guaranteed to decrease npts_remaining(k), so we will
          ! exit the loop sooner or later
          call assign_work(time, t_pe(i), nparts_loc_all(idx), &
                throughput_pe_all(idx), npts_remaining(k), n_epsilon(k))
        end do
        ! Consider stratum leaders postprocessing overhead
        i = stratum_leader(k)
        t_overhead(i) = t_overhead(i) + t_postproc(k)
        t_pe(i) = t_pe(i) + t_postproc(k)
        ! Now consider calstr overhead, calculate npestr and procforstra
        npestr(k) = 0
        do i = 0, nprs-1
          idx = i*nstrai + k
          if (nparts_loc_all(idx) > 0) then
            npestr(k) = npestr(k) + 1
            procforstra(k,i) = .true.
            ! We use t_calstr_pe instead of t_calstr_strat_avg(k) to have
            ! individual overhead estimate for each PE
            ! (stratum leaders have larger overhead than others, so if the
            ! leaders do not change then it is better estimate than the average)
            t_overhead(i) = t_overhead(i) + t_calstr_pe(idx)
            t_pe(i) = t_pe(i) + t_calstr_pe(idx)
          endif
        end do
      end do
      write(iunout,*) 'Work distribution table optimized in ', n, ' steps'
    endif
    ! Stratum_leader, npestr, nparts_loc are the main output
    call mpi_bcast(stratum_leader, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npestr, size(npestr), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_scatter(nparts_loc_all, nstrai, MPI_INTEGER,  &
             nparts_loc, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! procforstra is not used outside of cpes anymore, we could actually omit it
    CALL MPI_BCAST(PROCFORSTRA, SIZE(PROCFORSTRA), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  end subroutine

  subroutine simple_balanced_strategy(nparts_loc, npestr, stratum_leader)
    use eirmod_mpi
    ! Simple strategy consist of sharing large strata among all proc, while small
    ! strata are reparted between the processes. This strategy needs that strata
    ! are computed in order of decreasing number of histories in main steering routine
    integer, intent(inout), dimension(:) :: nparts_loc
    integer, intent(inout), dimension(:) :: npestr
    integer, intent(inout), dimension(:) :: stratum_leader

    integer, dimension(nstrai*nprs) :: nparts_loc_all
    real(kind=dp), dimension(0:nprs) :: t_pe
    integer, dimension(nstrai) :: npts_remaining
    integer :: i, k, p, idx, n, ierr
    real(kind=dp), dimension(nstrai*nprs) :: t_particles_all !< temporary variables
    integer, dimension(nstrai*nprs) :: nparts_processed_all  !< for MPI communication

    if (nprs == 1) return ! nothing to optimize for serial mode

    if(my_pe==0) then
      write(iunout,*) 'Simple allocation strategy'
    endif

    ! gather the processing time and particle numbers from other PEs
    call MPI_Gather(t_particles_loc, nstrai, MPI_DOUBLE_PRECISION, &
                        t_particles_all, nstrai, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, ierr)
    call MPI_Gather(n_particles_loc, nstrai, MPI_INTEGER, &
                        nparts_processed_all, nstrai, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)

    if (my_pe == 0) then
      t_pe = 0 ! t_pe(i) is used as a cumulated weight to distribute the strata

      ! Initialize output arrays
      !> the number of particles that needs to be distributed:
      npts_remaining = npts(1:nstrai)
      npestr = 0         !< numbre of PEs per stratum
      nparts_loc_all = 0

      n = sum(npts(1:nstrai))
      do k=1,nstrai
        ! If strata is large, share among all procs
        if (npts(k) > n/nprs/4) then
          i = get_next_pe(t_pe)
          stratum_leader(k) = i
          t_pe(i) = t_pe(i) + npts(k)/nprs/4

          do i = 0, nprs - 1
          idx = i*nstrai + k
          npestr(k) = npestr(k) + 1
          call assign_work(real(npts(k),kind=dp)/nprs, t_pe(i), nparts_loc_all(idx), &
                   1d0, npts_remaining(k), npts(k)/nprs/10)
          enddo
        else
          i = get_next_pe(t_pe)
          stratum_leader(k) = i
          t_pe(i) = t_pe(i) + npts(k)/10
          idx = i*nstrai + k
          npestr(k) = 1
          call assign_work(real(npts(k),kind=dp), t_pe(i), nparts_loc_all(idx), &
                   1d0, npts_remaining(k), npts(k))
        endif
      end do
    endif
    ! Stratum_leader, npestr, nparts_loc are the main output
    call mpi_bcast(stratum_leader, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(npestr, size(npestr), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call mpi_scatter(nparts_loc_all, nstrai, MPI_INTEGER,  &
             nparts_loc, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end subroutine

  subroutine assign_work(time, t_pe, nparts_loc, throughput, npts_remaining, &
                          n_epsilon)
  ! The wallclock time (s) to assign is the first argument. We might assign less
  ! time if there are not enough particles in npts_remaining, or slightly more,
  ! if the remaining particles would be less than n_epsilon(k).
    real(kind=dp), intent(in) :: time !< time to assign (sec)
    !> time that PE has worked so far (will be increased by time)
    real(kind=dp), intent(inout) :: t_pe
    !> we will increase nparts loc with the number of particles assgned
    integer, intent(inout) :: nparts_loc
    real(kind=dp), intent(in) :: throughput !< throughput of the PE
    !> the number of particles still remaining in the stratum
    integer, intent(inout) :: npts_remaining
    !> if npts_remaining would be less than this, then we assign these too
    integer, intent(in) :: n_epsilon
    integer n
    if (time < 0) then
      write(iunout,*) 'Warning, trying to assign negative time'
    endif
    n = int(time * throughput) !< number of particles that can be processed in time
    if (n<=0) then
      write(iunout,*) 'Warning, zero particles would be assigned'
      n = max(1, n_epsilon)
    endif
    if (npts_remaining < n) then
      n = npts_remaining
    endif
    t_pe = t_pe + n / throughput
    nparts_loc =  nparts_loc + n
    npts_remaining = npts_remaining - n
    if (npts_remaining.ne.0 .and. npts_remaining < n_epsilon) then
      nparts_loc = nparts_loc + npts_remaining
      t_pe = t_pe + npts_remaining / throughput
      npts_remaining = 0
    endif
  end subroutine

  function get_next_pe(time) result (ipe)
  ! Returns rank of PE who has the smallest working time
  ! ( = has most time available to do work)
    integer :: ipe
    real(kind=dp), intent(in), dimension(0:nprs-1) :: time
    integer, dimension(1) :: tmp
    tmp = minloc(time)
    ipe = tmp(1) - 1 ! because PEs are numbered from 0
  end function

  subroutine calculate_throughput(throughput_pe, throughput, &
             throughput_avg, throughput_strat_avg)
  ! Throughput = particles / sec
  ! It tells how many particles can be calculated in a second.
  ! All arguments are output arguments, but they are only defined on PE 0
  ! threrefore the intent(out) attribute is omitted.
    use eirmod_comsou &
   , only: nstrai
    use eirmod_mpi
    !> throughput_pe(i*nstrai + k) tells the throughput of PE i on stratum k
    real(kind=dp), dimension(:), allocatable :: throughput_pe
    !> throughput(k) is the aggregate throughput while processing stratum k
    !> It gets larger if we have more MPI tasks
    real(kind=dp), dimension(nstrai) :: throughput
    !> throughput_strat_avg(k) is the average throuphput of stratum (k) over all PE
    real(kind=dp), dimension(nstrai) :: throughput_strat_avg
    !> throughput averaged over all PEs and all strata
    real(kind=dp):: throughput_avg
    !> throughput_pe_avg(i) is the average throughput of PE i over all strata
    real(kind=dp), dimension(0:nprs-1) :: throughput_pe_avg
    real(kind=dp), dimension(nstrai*nprs) :: t_particles_all !< temporary variables
    integer, dimension(nstrai*nprs) :: nparts_processed_all  !< for MPI communication
    integer :: k, i, idx, n, ierr, n_strat
    real(kind=dp) :: t, t_strat
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external :: mpi_gather
#endif

    ! gather the processing time and particle numbers from other PEs
    call MPI_Gather(t_particles_loc, nstrai, MPI_DOUBLE_PRECISION, &
                        t_particles_all, nstrai, MPI_DOUBLE_PRECISION, &
                        0, MPI_COMM_WORLD, ierr)
    call MPI_Gather(n_particles_loc, nstrai, MPI_INTEGER, &
                        nparts_processed_all, nstrai, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)
    if (my_pe==0) then
      if(.not.allocated(throughput_pe)) then
        allocate(throughput_pe(nstrai*nprs))
      endif
      ! calculate the local throughput for all PEs
      ! also calculate the average and aggregate throughput
      n = 0      ! total number of particles processed
      t = 0.0D0  ! total processing time
      do k=1, nstrai
        n_strat = 0
        t_strat = 0
        throughput(k) = 0.0D0 ! aggregate throughput of stratum(k)
        do i = 0, nprs-1
          idx = i*nstrai+k
          if (nparts_processed_all(idx) > 0 .and. &
                  t_particles_all(idx) > 0.0D0) then
            n = n + nparts_processed_all(idx)
            t = t + t_particles_all(idx)
            n_strat = n_strat + nparts_processed_all(idx)
            t_strat = t_strat + t_particles_all(idx)
            throughput_pe(idx) = nparts_processed_all(idx) / t_particles_all(idx)
          else
            ! below we will define an estimate for this case too
            throughput_pe(idx) = 0.0D0
          endif
          throughput(k) = throughput(k) + throughput_pe(idx)
        enddo
        if (t_strat > 0) then
          throughput_strat_avg(k) = n_strat / t_strat
        endif
      enddo
      if ( t > 0) then
        throughput_avg = n / t
      else
        write(iunout,*) 'Something is wrong, we do not have any throughput measurements. '
        throughput_avg = 0.0D0
      endif
      if (throughput_avg == 0.0D0) then
        write(iunout,*) &
         'Error, no particles processed so far, cannot optimize'
        call MPI_abort(MPI_COMM_WORLD, -1, ierr)
      endif
      ! Now we can give an estimate for the missing throughput values
      do k = 1, nstrai
        if (throughput(k) == 0.0D0) then
          write(iunout,*) 'No throughput measurement for stratum ',k, &
                     ', assuming average value.'
          throughput_strat_avg(k) = throughput_avg
          throughput(k) = throughput_avg * nprs
          ! we multiply by nprs, because throughput(k) is not an average but
          ! an aggregate value
        endif
      enddo
      do i = 0, nprs-1
        ! We also calculate the average throughput per PE
        n = 0     ! particles processed by PE i
        t = 0.0D0 ! processing time of PE i
        do k = 1, nstrai
          idx = i*nstrai+k
          if (throughput_pe(idx)==0.0D0) then
            ! we estimate this with the average over the stratum
            throughput_pe(idx) = throughput_strat_avg(k) !throughput(k)/nprs
          endif
          n = n + nparts_processed_all(idx)
          t = t + t_particles_all(idx)
        end do
        if (t > 0) then
          throughput_pe_avg(i) = n / t
        else
          ! This value is used only for print_timing_info
          throughput_pe_avg(i) = 0
        endif
      end do
      if (print_timing_info) then
        write(iunout,*) 'throughput (particles/sec)', throughput
        write(iunout,*) 'throughput_strat_avg', throughput_strat_avg
        write(iunout,*) 'throughput_pe_avg', throughput_pe_avg
        write(iunout,*) 'throughput_avg', throughput_avg
        write(iunout,*) 'throughput 0', throughput_pe(1:nstrai)
      endif
    endif
  end subroutine

  subroutine calculate_ideal_time(throughput, t_calstr, &
           t_postproc, t_ideal, t_ideal_tot)
    !> aggregate throughput per stratum
    real(kind=dp), intent(in), dimension(:) :: throughput
    !> time to execute calstr (for every stratum)
    real(kind=dp), intent(in), dimension(:) :: t_calstr
    !> time to execute postprocessing (for every stratum)
    real(kind=dp), intent(in), dimension(:) :: t_postproc
    !> ideal execution time per stratum
    real(kind=dp), intent(out), dimension(:) :: t_ideal
    real(kind=dp), intent(out) :: t_ideal_tot !< total execution time
    integer :: k
    do k = 1, nstrai
     ! Ideally, stratum(k) would be processed in this much of time:
      if (npts(k) > 0) then
        t_ideal(k) = npts(k) / throughput(k)
        ! we could consider the postprocessing and calstr overhead too
        ! from arrays t_calstr and t_postproc, like:
        ! t_ideal(k) = t_ideal(k) + t_postproc(k)
       else
         t_ideal(k) = 0.0D0
       endif
    enddo
    ! In an ideal case, the parallel calculation would take this much time
    t_ideal_tot = sum(t_ideal)
    if (print_timing_info) then
      write(iunout,*) 't_ideal', t_ideal
      write(iunout,*) 't_ideal_tot', t_ideal_tot
    endif
  end subroutine

  subroutine calculate_calstr_time(t_calstr_strat_avg, t_calstr_avg, t_calstr_pe)
    use eirmod_comsou &
     , only: nstrai
    use eirmod_mpi
    ! All arguments are output arguments, but they are only defined at PE 0
    !> t_calstr_strat_avg(k) is the average calstr time of stratum k, over
    !> iterations of Eirene and PEs
    real(kind=dp), dimension(:), allocatable :: t_calstr_strat_avg
    !> average calstr time over all strata, PE, and iterations
    real(kind=dp) :: t_calstr_avg
    !> Estimated calstr time for each PE
    real(kind=dp), dimension(:), allocatable :: t_calstr_pe
    integer, dimension(nstrai) :: n_calstr_sum
    integer :: ierr, k, i, idx
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external :: mpi_gather, mpi_reduce
#endif

    if(.not.allocated(t_calstr_strat_avg)) then
      if (my_pe==0) then
        allocate(t_calstr_strat_avg(nstrai))
        t_calstr_strat_avg = 0.0_dp
      else
        ! t_calstr_strat_avg should be only significant at the root (PE 0)
        ! but intel mpi gives an error if it is not allocated
       allocate(t_calstr_strat_avg(1))
      endif
    endif
    if(.not.allocated(t_calstr_pe)) then
      if(my_pe==0) then
        allocate(t_calstr_pe(nstrai*nprs))
      else
        allocate(t_calstr_pe(1))
      endif
    endif
    call MPI_Gather(t_calstr_loc, nstrai, MPI_DOUBLE_PRECISION, &
                    t_calstr_pe, nstrai, MPI_DOUBLE_PRECISION,  &
                    0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(t_calstr_loc, t_calstr_strat_avg, nstrai, &
                 MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(n_calstr_loc, n_calstr_sum, nstrai, &
                MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! n_calstr_sum(k) is the number of PEs who called calstr for stratum(k)
    if ( my_pe == 0) then
      do k = 1, nstrai
        if (n_calstr_sum(k) > 0) then
          t_calstr_strat_avg(k) = t_calstr_strat_avg(k) / n_calstr_sum(k)
        else
          t_calstr_strat_avg(k) = 0.0_dp
        endif
      end do
      t_calstr_avg = sum(t_calstr_strat_avg) / nsteff
      if (print_timing_info) then
        write(iunout,*) 't_calstr_strat_avg ', t_calstr_strat_avg
        write(iunout,*) 't_calstr_avg ', t_calstr_avg
      endif
      ! now we fill in zero values with estimates
      do k = 1, nstrai
        if (n_calstr_sum(k) == 0) then
          t_calstr_strat_avg(k) = t_calstr_avg
        endif
        do i = 0, nprs-1
          idx = i*nstrai + k
          if (t_calstr_pe(idx) == 0.0_dp) then
            t_calstr_pe(idx) = t_calstr_strat_avg(k)
          endif
        end do
      enddo
    endif
    n_calstr_loc = 0
  end subroutine

  subroutine calculate_postproc_time(t_postproc_avg)
    use eirmod_mpi
    use eirmod_comsou &
     , only: nstrai
    !> The average postprocessing time of stratum k is returned in
    !> t_postproc_avg(k) (average over the iterations of Eirene)
    real(kind=dp), dimension(nstrai) :: t_postproc_avg
    integer :: ierr, k, n_tot
    real(kind=dp) :: n, tmp
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
    external :: mpi_reduce
#endif

    call MPI_REDUCE(t_postproc_loc, t_postproc_avg, nstrai, &
           MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(n_postproc_executed, n_tot, 1, &
           MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (my_pe == 0) then
    ! The number of all the postprocessing steps so far is: n_postproc_executed
    ! Every stratum is postprocessed once in every iteration. The number of
    ! postprocessing steps per stratum is n_postproc_executed/nsteff (which
    ! should be equal to the number of iterations).
    ! We calculate the average over each stratum:
      n = real(n_tot)/real(nsteff)
      t_postproc_avg = t_postproc_avg / n
      if(print_timing_info) then
        write(iunout,*) 'postproc n n', n_postproc_executed, nsteff
        write(iunout,*) 't_postproc ', t_postproc_avg
      endif
      tmp = sum(t_postproc_avg) / nsteff
      do k = 1,nstrai
        if (t_postproc_avg(k)==0.0_dp) then
          t_postproc_avg(k) = tmp
        endif
      end do
    endif
  end subroutine

  subroutine time_calstr(istra, time)
    integer, intent(in) :: istra !< stratum idx
    real(kind=dp), intent(in) :: time !< time (s)
    ! We calculate a weighted average with the previous values
    ! Let t(n) be the average value at iteration n
    ! t(n) = N * (t(n) + t(n-1)/2 + t(n-2)/4 + t(k,n-3)/8 + ... )
    ! Here N is the normalization factor:
    ! N = (1 + 1/2 + 1/4 + 1/8 + ...)^-1 = 1 / 2 (approximately)
    t_calstr_loc(istra) = (t_calstr_loc(istra) + time) / 2
    n_calstr_loc(istra) = n_calstr_loc(istra) + 1
  end subroutine

  subroutine time_particles(istra, time, n)
    integer, intent(in) :: istra !< stratum idx
    real(kind=dp), intent(in) :: time !< time (s)
    integer, intent(in) :: n !< number of particles
    n_particles_loc(istra) = n_particles_loc(istra) + n
    t_particles_loc(istra) = t_particles_loc(istra) + time
  end subroutine

  subroutine time_postproc(istra, time)
    integer, intent(in) :: istra !< stratum idx
    real(kind=dp), intent(in) :: time !< time (s)
    t_postproc_loc(istra) = t_postproc_loc(istra) + time
    n_postproc_executed = n_postproc_executed + 1
  end subroutine
end module
