cdr Nov. 17  commenting started
cdr July 18  remove nsteff, redundant

      MODULE EIRMOD_CPES

      USE EIRMOD_PARMMOD, ONLY: IUNMEM, NSTRA
      USE EIRMOD_MPI, ONLY: MPI_COMM_NULL, MPI_SET_OWN_IO_UNIT
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_COMSOU
      USE EIRMOD_CAI

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CPES, EIRENE_DEALLOC_CPES, EIRENE_INIT_CPES
      public :: I_am_leader
      public :: create_all_communicators
      public :: get_global_particle_idx, set_index_shift
      public :: init_work_distribution
      public :: get_leader_comm, get_stratum_comm
      public :: need_calstr, calc_stratum
      public :: create_leader_comm
      public :: print_work_distribution_table

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
cdr  npesta(istra): no. of master processor for ISTRA
cdr  npestr(istra): total number of processors working on ISTRA
     I         NPESTR(:), NPESTA(:)

      INTEGER, PUBLIC, SAVE ::
     I         NPRS, MY_PE

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     I         LEXIT(:)

      LOGICAL, PUBLIC, SAVE :: NLIDENT

CVKMPI CORRESPONDENCE TABLE "STRATA VERSUS PROCESSOR"
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: PROCFORSTRA(:,:)

      integer, public :: current_particle_idx

      ! Identifiers for work distribution strategies:
      integer, public, parameter :: STRATEGY_UNDEFINED = -1

      !> Assign PEs to strata by proportional allocation.
      !> The particles within a stratum are divided equally
      !> among the PEs that calculate the stratum.
      !> Usually leads to load imbalance.
      integer, public, parameter :: STRATEGY_ORIGINAL = 0

      !> All PEs Calculate All strata.
      !> It is a balanced workload distribution, but
      !> the communication costs can be high.
      integer, public, parameter :: STRATEGY_APCAS = 1

      !> Assign PEs to strata to minimize communication,
      !> and divide the number of
      !> particles between PEs in a way that balances workload.
      integer, public, parameter :: STRATEGY_BALANCED = 2

      !> "Embarrassingly" parallel scheme.
      !> Each PE calculates the same number of histories.
      integer, public, parameter :: STRATEGY_EMBARRASS = 3

      !> Assign PEs to strata to minimize communication,
      !> and divide the number of particles between PEs
      !> evenly by using CPU time of previous steps.
      integer, public, parameter :: STRATEGY_IMPROVED = 4

      !> EDGE2D strategy, for few bigs and many small strata,
      !> all procs shares the bigs at begining, smalls are
      !> distributed between procs
      integer, public, parameter :: STRATEGY_SIMPLE = 5

      !> STRATEGY_DEFAULT and NPRLL_DEFAULT are set in user_defaults
      integer, public, save :: STRATEGY_DEFAULT
      integer, public, save :: NPRLL_DEFAULT

      !> stores which strategy was chosen during initialization
      integer, public :: work_distribution_strategy

      !> strategy read from input file
      integer, public :: input_distribution_strategy

      !> stratum_leader(i) is the rank of the leader for stratum(i)
      integer, dimension(:), allocatable, public :: stratum_leader
      !> number of particles calculated locally
      integer, dimension(:), allocatable, public :: nparts_loc

      !> to help mapping particle idx to a global idx,
      !> we store on each PE what is the starting global idx
      integer, dimension(:), allocatable, public :: index_shift
      integer, dimension(:), allocatable, public :: stratum_comm
      integer, public :: leader_comm = MPI_COMM_NULL

      !> Whether we have new member of stratum leaders
      logical, public :: new_leader = .true.

      CONTAINS

!> creates communicators for PEs working on same stratum and a
!> communicator for all stratum leaders
      subroutine create_all_communicators()
        use eirmod_comsou, only: nstrai, nlsron
        use eirmod_mpi, only: MPI_COMM_NULL
        integer :: k
        call free_communicators(nstrai)
! Initialize communicators
        do k=1, nstrai
          if (nlsron(k)) then
            stratum_comm(k) = create_stratum_comm(k)
          else
            stratum_comm(k) = MPI_COMM_NULL
          endif
        end do
        leader_comm = create_leader_comm()
      end subroutine

!> returns true if the calling PE should do any work on stratum_idx
      logical function calc_stratum(stratum_idx)
        integer, intent(in) :: stratum_idx
        calc_stratum = nparts_loc(stratum_idx) > 0
      end function

!> creates the communicator for all stratum leaders
      function create_leader_comm() result(comm)
        integer :: comm
        if (my_pe==0) then
          write(iunout,*) 'Creating communicator for stratum leaders'
        end if
        comm = create_communicator(I_am_leader())
      end function

!> creates the communicators for PEs working on same stratum
      function create_stratum_comm(stratum_idx) result(comm)
        integer :: comm
        integer, intent(in) :: stratum_idx
        if (my_pe==0) then
          write(iunout,*) 'Creating communicator within stratum ',
     &     stratum_idx
        end if
        comm = create_communicator(nparts_loc(stratum_idx)>0,
     &           stratum_leader(stratum_idx))
      end function

!> creates a communicator for the PEs which call it with .true. argument
!> It is a collective call: every PE in MPI_COMM_WORLD should call it together.
      function create_communicator(include_me, zero_rank) result(comm)
        use eirmod_mpi
!> whether to include the calling PE
        logical, intent(in) :: include_me
!> the rank of the PE which should become rank 0 in the new communicator
        integer, optional :: zero_rank
        integer :: comm
        integer :: n, k, ierr
        integer, dimension(nprs) :: ranks_tmp
        logical, dimension(nprs) :: include_proc
        integer :: group_world
        integer :: group !, gsize, grank
        logical :: found
        character*13 hlp_frm
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
        external :: mpi_allgather
#endif

        call MPI_allGather(include_me, 1, MPI_LOGICAL,
     &          include_proc, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
! now collect all the ranks where include_proc is true
        n = 0
        do k=0, nprs-1
          if (include_proc(k+1)) then
            n = n + 1
            ranks_tmp(n) = k
          end if
        end do

        if (present(zero_rank)) then
          ! put zero_rank into ranks_tmp(1)
          ! but first free the place ranks_tmp(1)
          found = .false.
          do k=1,n
            if (ranks_tmp(k) == zero_rank) then
              ranks_tmp(k) = ranks_tmp(1)
              found = .true.
            endif
          end do
          if (found) then
            ranks_tmp(1) = zero_rank
          else
            write(iunout,*) "Error, rank ", zero_rank,
     &         " is not a member of the new communicator"
          endif
        end if
        if (my_pe == 0) then
          write(hlp_frm,'(a,i4,a)') '(1x,a,',max(n,1),'i5)'
          write(iunout,hlp_frm)
     &     'Creating group from ranks ', ranks_tmp(1:n)
        end if

        call mpi_comm_group(MPI_COMM_WORLD, group_world, ierr)
        call mpi_group_incl(group_world, n, ranks_tmp, group, ierr)
        call mpi_comm_create(MPI_COMM_WORLD, group, comm, ierr)
        !if (include_me) then
        !  call mpi_group_size(group, gsize, ierr)
        !  call mpi_group_rank(group, grank, ierr)
        !  write(iunout,*) 'My rank in the new communicator ', grank,
        !                                                ' (', gsize,')'
        !endif
        call mpi_group_free(group_world, ierr)
        call mpi_group_free(group, ierr)
      end function

      !> frees all communicators
      subroutine free_communicators(nstra)
        use eirmod_mpi
        integer, intent(in) ::nstra
        integer k, ierr
        do k = 1, nstra
          if(stratum_comm(k) .ne. MPI_COMM_NULL) then
            call mpi_comm_free(stratum_comm(k), ierr)
          end if
        end do
        if (leader_comm .ne. MPI_COMM_NULL) then
          call mpi_comm_free(leader_comm, ierr)
        endif
      end subroutine

      function get_global_particle_idx(stratum_idx, idx)
     &  result(global_idx)
        use eirmod_comusr
     ,   , only: iiter
        use eirmod_comsou
     ,   , only: npts, nstrai
        integer, intent(in) :: stratum_idx
        !> local particle index (1..nparts_loc)
        integer, intent(in) :: idx
        integer :: global_idx
        global_idx = index_shift(stratum_idx) + idx
     &    + (iiter-1)*sum(npts(1:nstrai))
      end function

      !> returns the communicator of the stratum leader
      function get_leader_comm() result(comm)
        integer :: comm
        comm = leader_comm
      end function

      !> returns the communicator of the PEs working
      !> on stratum `stratum_idx`
      function get_stratum_comm(stratum_idx) result(comm)
        integer, intent(in) :: stratum_idx
        integer :: comm
        comm = stratum_comm(stratum_idx)
      end function

      !> Checks whether the calling PE is the leader of
      !> stratum with stratum_idx,
      !> or any stratum if stratum_idx is not present
      logical function I_am_leader(stratum_idx)
        use eirmod_comsou, only: nlsron
        integer, optional, intent(in) :: stratum_idx
        if (present(stratum_idx)) then
          I_am_leader = stratum_leader(stratum_idx) == my_pe
        else
          I_am_leader = any(npesta==my_pe .and. nlsron)
        end if
      end function

      subroutine init_apcas_strategy(nstra, npts)
        integer, intent(in) :: nstra !< number of strata
        !> number of points per strata
        integer, intent(in), dimension(NSTRA) :: npts
        integer :: k, diff, proc_id
        if (my_pe==0) then
          write(iunout,*) ''
          write(iunout,*) 'Initializing APCAS workload distribution'
          write(iunout,*) 'All PEs Calculate All Strata'
        endif
        ! every stratum will have nprs processors working on it
        npestr = nprs
        ! divide the particles equally between PEs
        do k=1, nstra
          if (npts(k).gt.0) then
            nparts_loc(k) = npts(k) / nprs
            procforstra(k,:) = .true.
          else
            nparts_loc(k) = 0
            procforstra(k,:) = .false.
            procforstra(k,0) = .true.
          end if
        end do
        ! correct for remainder
        do k=1, nstra
          diff =  mod(npts(k), nprs)
          if (diff > 0) then
            if (nparts_loc(k)==0) then
              ! if there are less particles than PEs then
              ! let proc 0 do the work
              proc_id = 0
            else
              proc_id = nprs - 1
            endif
            if (my_pe == proc_id) then
              nparts_loc(k) = nparts_loc(k) + diff
            endif
          endif
         end do
         ! rank 0 will be the leader for all strata
         npesta = 0
         stratum_leader = 0
      end subroutine

      subroutine init_work_distribution(xtim, xx1, strategy)
        use eirmod_precision
     ,   , only: dp
        use eirmod_parmmod
     ,   , only: nstra
        use eirmod_comsou
     ,   , only: nstrai, npts
        use eirmod_balanced_strategy
        use eirmod_improved_strategy
        use eirmod_mpi
        real(dp), intent(inout) :: xtim(0:nstra)
        real(dp), intent(in) :: xx1
        integer, optional :: strategy
        integer :: ierr
        logical, save :: first_call = .true.
        integer, dimension(NSTRA) :: old_stratum_leader
        integer, dimension(NSTRA) :: ipts
        external :: eirene_pedist, eirene_leer, eirene_masage,
     &              eirene_exit_own
#if( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
        external :: mpi_bcast
#endif

        if (first_call) then
          old_stratum_leader = -1
          if (present(strategy)) then
            if (strategy == STRATEGY_UNDEFINED) then
              NPRLL = NPRLL_DEFAULT
              work_distribution_strategy = STRATEGY_DEFAULT
            else
!pb              NPRLL = 3
              work_distribution_strategy = strategy
            end if
          else
            NPRLL = 1
            work_distribution_strategy = STRATEGY_ORIGINAL
          endif
        else
          old_stratum_leader = stratum_leader
        endif

        call mpi_comm_size(MPI_COMM_WORLD, nprs, ierr)
        IF (NPRS == 1) THEN

! 1 PROCESSOR: ALL STRATA ARE DONE BY PROCESSOR 0
!              XTIM REMAINS UNCHANGED

          CALL EIRENE_LEER(1)
          CALL EIRENE_MASAGE ('Eirene run in serial mode')
          PROCFORSTRA(1:NSTRAI,0) = NLSRON(1:NSTRAI)
          npestr(1:nstrai) = 1
          npesta(1:nstrai) = 0
          nparts_loc(1:nstrai) = npts(1:nstrai)
          stratum_leader = 0

        ELSE

          select case(work_distribution_strategy)
            case (STRATEGY_EMBARRASS, STRATEGY_ORIGINAL)
            ! The correct value of XMCT is found on PE 0
            ! We distribute it, because it might be used on other PEs
              call mpi_bcast(xmct, nstra+1, MPI_DOUBLE_PRECISION, 0,
     &                     MPI_COMM_WORLD, ierr)
              call eirene_pedist(xtim,xx1)
            case (STRATEGY_APCAS)
              ipts(1:nstrai) = npts(1:nstrai)
              call init_apcas_strategy(nstrai, ipts)
            case (STRATEGY_BALANCED, STRATEGY_IMPROVED)
              if (first_call) then
                if (work_distribution_strategy == STRATEGY_BALANCED)
     >           then
                  call init_balanced_strategy(ierr)
                else
                  call init_improved_strategy(ierr)
                endif
                if (ierr.ne.0) then
                  if(my_pe==0) then
                    if (work_distribution_strategy == STRATEGY_BALANCED)
     >               then
                      write(iunout,*)
     &                 'Error in init_balanced_strategy. ',
     &                 'Falling back to STRATEGY_APCAS.'
                    else
                      write(iunout,*)
     &                 'Error in init_improved_strategy. ',
     &                 'Falling back to STRATEGY_APCAS.'
                    endif
                  endif
                  work_distribution_strategy = STRATEGY_APCAS
                endif
              ! we start with APCAS, later we optimize
                ipts(1:nstrai) = npts(1:nstrai)
                call init_apcas_strategy(nstra,ipts)
              else
                if (work_distribution_strategy
     &               == STRATEGY_BALANCED) then
                  call opt_balanced_strategy(nparts_loc, npestr,
     &                 stratum_leader, procforstra)
                else
                  call opt_improved_strategy(nparts_loc, npestr,
     &                 stratum_leader, procforstra)
                endif
                npesta = stratum_leader
              endif
            case (STRATEGY_SIMPLE)
              call init_balanced_strategy(ierr)
              call simple_balanced_strategy(nparts_loc, npestr,
     &                stratum_leader)
              npesta = stratum_leader
            case default
              if (my_pe==0) then
                write(iunout,*)
     &           "Error, wrong selector for parallel strategy"
                call eirene_exit_own(1)
              endif
          end select

        ENDIF

        new_leader = any(old_stratum_leader .ne. stratum_leader)
        call create_all_communicators
        call set_index_shift
        call print_work_distribution_table
        first_call = .false.
      end subroutine

      !> returns true if the calling PE should call eirene_calstr
      logical function need_calstr(stratum_idx)
        integer, intent(in) :: stratum_idx
        need_calstr = npestr(stratum_idx) > 1
      end function

      subroutine print_work_distribution_table
        use eirmod_precision
     ,   , only: dp
        use eirmod_mpi
        use eirmod_comsou
     ,   , only: nstrai
        use eirmod_balanced_strategy
     ,   , only: throughput_pe_all, t_overhead
        use eirmod_improved_strategy
     ,   , only: throughput2, t_overhead2
        integer, dimension(nstrai*nprs) :: nparts_tmp
        integer, parameter :: ncol = 20
        character(len=9), dimension(nstrai+2) :: ch_tmp
        character(len=16) :: hlp_fmt
        character(len=256) :: long_line
        character(len=1) :: leader_mark
        integer j, k, ierr, n, in, it, ie, isa, ise
        real(kind=dp) :: tot, time_avg, time_std
        real(kind=dp), dimension(0:nprs-1) :: time_pe
        external :: eirene_leer
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
        external :: mpi_gather
#endif

        call MPI_Gather(nparts_loc, nstrai, MPI_INTEGER,
     &     nparts_tmp, nstrai, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        if (my_pe==0) then
          call eirene_leer(1)
          write(iunout,*)
     &     'Work distribution table (number of particles, ',
     &     '# marks leader) ========'
          ie = ((nstrai-1)/ncol)+1
          do it = 1, ie
            isa = (it-1)*ncol + 1
            if (it.eq.ie) then
              in = mod(nstrai,ncol)
              if (in.eq.0 .and. nstrai.gt.0) in = ncol
              ise = nstrai
            else
              in = ncol
              ise = it*ncol
            end if
            if (in.eq.0) cycle
            if (it.eq.1) then
              if (in.gt.1) then
                write(hlp_fmt,'(a,i2,a)') '(A15,',in-1,'I9)'
                write(iunout,hlp_fmt) ' PE  \ Strat  1', (k, k=2,in)
              else
                write(iunout,'(a)') ' PE  \ Strat  1'
              end if
            else
              if (in.gt.1) then
                write(hlp_fmt,'(a,i2,a)') '(A13,I2,',in-1,'I9)'
                write(iunout,hlp_fmt) ' PE  \ Strat ', isa,
     .           (k, k=isa+1,ise)
              else
                write(iunout,'(a13,i2)') ' PE  \ Strat ',isa
              end if
            end if
            ch_tmp = '---------'
            write(hlp_fmt,'(a,i2,a)') '(A7,',in,'A9)'
            write(long_line,hlp_fmt) '-------',ch_tmp(isa:ise)
            write(iunout,'(a)') trim(long_line)
            write(hlp_fmt, '(a,i2,a)') '(I4,A3,',in,'A9)'
            do k = 0, nprs-1
              do j = isa, ise
                if (stratum_leader(j)==k) then
                  leader_mark = '#'
                else
                  leader_mark = ' '
                endif
                if (nparts_tmp(k*nstrai + j)>0) then
                  write(ch_tmp(j), fmt='(I8,A1)')
     &              nparts_tmp(k*nstrai+j),leader_mark
                else
                  write(ch_tmp(j), fmt='(A8,A1)') '        ',leader_mark
                end if
              end do
              write(iunout,hlp_fmt) k, ' | ', ch_tmp(isa:ise)
            end do
            ch_tmp = '---------'
            write(hlp_fmt,'(a,i2,a)') '(A7,',in,'A9)'
            write(long_line,hlp_fmt) '-------',ch_tmp(isa:ise)
            write(iunout,'(a)') trim(long_line)
            do j = isa, ise
              n = 0
              do k = 0, nprs - 1
                n = n + nparts_tmp(k*nstrai + j)
              end do
              write(ch_tmp(j), fmt='(I8,A1)') n,' '
            end do
            write(hlp_fmt,'(a,i2,a)') '(A7,',in,'A9)'
            write(long_line,hlp_fmt) ' sum | ', ch_tmp(isa:ise)
            write(iunout,'(a)') trim(long_line)
            call eirene_leer(1)
            if(allocated(throughput_pe_all) .or.
     &         allocated(throughput2)) then
              call eirene_leer(1)
              write(iunout,*)
     .         'Work distribution (estimated working time',
     &         ' in sec, # marks leader) ===='
              do k = isa+1, ise
                write(ch_tmp(k),fmt='(I9)') k
              end do
              if (it.eq.ie) then
                write(ch_tmp(ise+1),fmt='(A)') '  |overhe'
                write(ch_tmp(ise+2),fmt='(A)') 'ad|   sum'
                if (it.eq.1) then
                  write(hlp_fmt,'(a,i2,a)') '(A15,',in+1,'A9)'
                  write(long_line,hlp_fmt) ' PE  \ Strat  1',
     &             ch_tmp(isa+1:ise+2)
                else
                  write(hlp_fmt,'(a,i2,a)') '(A13,I2,',in+1,'A9)'
                  write(long_line,hlp_fmt) ' PE  \ Strat ',isa,
     &             ch_tmp(isa+1:ise+2)
                end if
              else
                if (it.eq.1) then
                  if (in.gt.1) then
                    write(hlp_fmt,'(a,i2,a)') '(A15,',in-1,'A9)'
                    write(long_line,hlp_fmt) ' PE  \ Strat  1',
     &               ch_tmp(isa+1:ise)
                  else
                    write(long_line,'(a)') ' PE \ Strat  1'
                  end if
                else
                  if (in.gt.1) then
                    write(hlp_fmt,'(a,i2,a)') '(A13,I2,',in-1,'A9)'
                    write(long_line,hlp_fmt) ' PE  \ Strat ',isa,
     &               ch_tmp(isa+1:ise)
                  else
                    write(long_line,'(a,i2)') ' PE  \ Strat ',isa
                  end if
                end if
              end if
              write(iunout,'(a)') trim(long_line)
              ch_tmp = '---------'
              if (it.eq.ie) then
                write(hlp_fmt,'(a,i2,a)') '(A7,',in+2,'A9)'
                write(long_line,hlp_fmt) '-------', ch_tmp(isa:ise+2)
              else
                write(hlp_fmt,'(a,i2,a)') '(A7,',in,'A9)'
                write(long_line,hlp_fmt) '-------', ch_tmp(isa:ise)
              end if
              write(iunout,'(a)') trim(long_line)
              if (it.eq.ie) then
                write(hlp_fmt, '(a,i2,a)') '(I4,A3,',in+2,'A9)'
              else
                write(hlp_fmt, '(a,i2,a)') '(I4,A3,',in,'A9)'
              end if
              do k = 0, nprs-1
                tot = 0.0_DP
                do j = isa, ise
                  if (stratum_leader(j)==k) then
                    leader_mark = '#'
                  else
                    leader_mark = ' '
                  endif
                  if (nparts_tmp(k*nstrai + j)>0) then
                    if (allocated(throughput_pe_all)) then
c                     BALANCED strategy
                      write(ch_tmp(j), fmt='(F8.3,A1)')
     &                     real(nparts_tmp(k*nstrai+j))/
     &                     throughput_pe_all(k*nstrai+j),
     &                     leader_mark
                      tot = tot + real(nparts_tmp(k*nstrai+j))/
     &                     throughput_pe_all(k*nstrai+j)
                    else
c                     IMPROVED strategy
                      write(ch_tmp(j), fmt='(F8.3,A1)')
     &                     nparts_tmp(k*nstrai + j) / throughput2(j),
     &                     leader_mark
                      tot = tot + nparts_tmp(k*nstrai+j)/throughput2(j)
                    endif
                  else
                    write(ch_tmp(j), fmt='(A8,A1)') '        ',
     &                 leader_mark
                  end if
                end do
                if (it.eq.ie) then
                  if (allocated(throughput_pe_all)) then
c                   BALANCED strategy
                    if (t_overhead(k).lt.1000.0_DP) then
                      write(ch_tmp(ise+1),fmt='(A2,F7.3)') ' |',
     &                     t_overhead(k)
                    else
                      write(ch_tmp(ise+1),fmt='(A2,F7.2)') ' |',
     &                     t_overhead(k)
                    end if
                    tot = tot + t_overhead(k)
                  else
c                   IMPROVED strategy
                    if (t_overhead2(k).lt.1000.0_DP) then
                      write(ch_tmp(ise+1),fmt='(A2,F7.3)') ' |',
     &                     t_overhead2(k)
                    else
                      write(ch_tmp(ise+1),fmt='(A2,F7.2)') ' |',
     &                     t_overhead2(k)
                    end if
                    tot = tot + t_overhead2(k)
                  endif
                  if (tot.lt.1000.0_DP) then
                    write(ch_tmp(ise+2),fmt='(A2,F7.3)') ' |', tot
                  else
                    write(ch_tmp(ise+2),fmt='(A2,F7.2)') ' |', tot
                  end if
                  write(iunout,hlp_fmt) k, ' | ', ch_tmp(isa:ise+2)
                  time_pe(k) = tot
                else
                  write(iunout,hlp_fmt) k, ' | ', ch_tmp(isa:ise)
                end if
              end do
              call eirene_leer(1)
              if (it.eq.ie) then
                time_avg = sum(time_pe) / nprs
                time_std = 0.0_DP
                do k = 0, nprs - 1
                  time_std = time_std + (time_pe(k) - time_avg)**2
                end do
                time_std = sqrt(time_std / nprs)
                if (time_avg.lt.10000.0_DP) then
                  write(iunout,fmt='(A,F8.3,A4,E9.2,A4)')
     &             'Estimated working time / PE', time_avg, ' +/-',
     &                                            time_std, ' sec'
                else
                  write(iunout,fmt='(A,F8.2,A4,E9.2,A4)')
     &             'Estimated working time / PE', time_avg, ' +/-',
     &                                            time_std, ' sec'
                end if
                call eirene_leer(1)
              end if
            end if
          end do
        endif
      end subroutine

      subroutine set_index_shift !init_global_index
        use eirmod_comsou
     ,   , only: npts, nstrai
        use eirmod_mpi
        integer :: k, ierr
        integer, dimension(nprs) :: tmp
        integer :: grank
        external :: eirene_masage, eirene_exit_own
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
        external :: mpi_allgather
#endif

        index_shift(1) = 0
        do k = 2, nstrai
          index_shift(k) = index_shift(k-1) + npts(k-1)
        end do
        ! now adjust with the starting index within the stratum
        do k = 1, nstrai
          if (nparts_loc(k) > 0) then
            call MPI_allGather(nparts_loc(k), 1, MPI_INTEGER, tmp, 1,
     &             MPI_INTEGER, stratum_comm(k), ierr)
            call mpi_comm_rank(stratum_comm(k), grank, ierr)
            if (ierr /= mpi_success) then
              call eirene_masage
     .         ('ERROR IN SUBROUTINE set_index_shift at mpi_comm_rank.')
              call eirene_exit_own(1)
            end if
            index_shift(k) = index_shift(k) + sum(tmp(1:grank))
            ! Note that tmp is indexed from 1,
            ! while grank is counted from zero.
            ! Therefore, as intended,
            ! the sum does not include the particles
            ! in the local stratum (which is stored in tmp(grank+1))
          endif
        end do
      end subroutine

      SUBROUTINE EIRENE_ALLOC_CPES(ICAL)
      USE EIRMOD_BALANCED_STRATEGY
     , , ONLY: ALLOCATE_BALANCED_STRATEGY
      USE EIRMOD_IMPROVED_STRATEGY
     , , ONLY: ALLOCATE_IMPROVED_STRATEGY

      INTEGER ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(NPESTR)) RETURN

        ALLOCATE (NPESTR(NSTRA))
        ALLOCATE (NPESTA(NSTRA))

        ALLOCATE(PROCFORSTRA(NSTRA,0:NPRS-1))

        ALLOCATE(nparts_loc(NSTRA))
        ALLOCATE(index_shift(NSTRA))
        ALLOCATE(stratum_leader(NSTRA))
        ALLOCATE(STRATUM_COMM(NSTRA))
        STRATUM_COMM = MPI_COMM_NULL

        call allocate_balanced_strategy
        call allocate_improved_strategy
        WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CPES(1) ',(4*NSTRA+2*(NSTRA+1)+NPRS)*4 +
     .                  NSTRA*NPRS*4

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(LEXIT)) RETURN
        ALLOCATE (LEXIT(0:NPRS-1))

        WRITE (IUNMEM,'(A,T25,I15)')
     .      ' CPES(2) ', NPRS*4

      END IF

      CALL EIRENE_INIT_CPES(ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CPES


      SUBROUTINE EIRENE_DEALLOC_CPES
      USE EIRMOD_BALANCED_STRATEGY
     , , ONLY: DEALLOCATE_BALANCED_STRATEGY
      USE EIRMOD_IMPROVED_STRATEGY
     , , ONLY: DEALLOCATE_IMPROVED_STRATEGY

      IF (.NOT.ALLOCATED(NPESTR)) RETURN

      DEALLOCATE (NPESTR)
      DEALLOCATE (NPESTA)

      DEALLOCATE (LEXIT)

      DEALLOCATE(PROCFORSTRA)

      DEALLOCATE(NPARTS_LOC)
      DEALLOCATE(INDEX_SHIFT)
      DEALLOCATE(STRATUM_LEADER)
      CALL FREE_COMMUNICATORS(NSTRA)
      DEALLOCATE(STRATUM_COMM)
      CALL DEALLOCATE_BALANCED_STRATEGY
      CALL DEALLOCATE_IMPROVED_STRATEGY
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CPES

      SUBROUTINE EIRENE_INIT_CPES(ICAL)
      INTEGER ICAL

      IF (ICAL == 1) THEN

        NPESTR = 0
        NPESTA = 0

c  correspondence table: Strata vs. PEs
        PROCFORSTRA=.TRUE. !  Trivial parallelisation:
                           !  All PEs work on all strata

        CALL MPI_SET_OWN_IO_UNIT(IUNOUT)

      ELSE IF (ICAL == 2) THEN

        LEXIT = .FALSE.

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CPES

      END MODULE EIRMOD_CPES
