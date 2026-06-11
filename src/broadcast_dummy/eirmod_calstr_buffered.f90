module eirmod_calstr_buffered
! A buffered version of calstr to allow non-blocking reduction
!
! All the data that is communicated in calstr are packed into a few large
! buffers (one for each data type), and reduced using MPI_IREDUCE.
!
! If you add any variables to calstr, then the same has to be added to this
! module too. Please add them to the following three subroutines:
!  - allocate_calstr_buffer,
!  - calstr_pack_all,
!  - calstr_unpack_all
!
! Note that the reduction operation is not progressing magically in the
! background after we call MPI_IREDUCE. We have to check the status of the
! reductions regularily using MPI_TEST, and this ensures that the reduction is
! progressing. This is implemented in calstr_progress_message.
!
! Alternatively one could use a thread in the background to manage the MPI
! communication (Intel MPI does that automatically if we set
! I_MPI_ASYNC_PROGRESS=1, but that was 40% slower using intel MPI 5.1).
!
! This module has one knob for tuning the communication: n_progress_check.
! This tells how many times we call MPI_TEST during the particle loop.

  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use eirmod_cestim
  use eirmod_cspez
  use eirmod_comprt
  use eirmod_csdvi
  use eirmod_coutau
  use eirmod_mpi
  use eirmod_cstep
  use eirmod_cai

  implicit none

  private
  public :: allocate_calstr_buffer, deallocate_calstr_buffer
  public :: calstr_buffered_finish
  public :: calstr_progress_message
  public :: eirene_calstr_buffered

  !> number of buffers (1 double and 1 logical) defined in eirmod_mpi
  !  integer, parameter :: N_BUFFERS = 2
  !> Identifies the communication operations, one request object for each buffer
  integer, save, dimension(N_BUFFERS) :: calstr_request =  MPI_REQUEST_NULL

  real(kind=dp), allocatable, dimension(:) :: calstr_buffer_d
  integer :: pos_d !< position in the double precision buffer

  logical, allocatable, dimension(:) :: calstr_buffer_l
  integer :: pos_l !< position in the logical buffer

  ! To ensure that the non-blocking reduce operations are progressing, we
  ! periodically call MPI_TEST for the calstr_request array.
  ! We would like to perform n_progress_check test.
  integer, parameter :: n_progress_check = 100

  ! For each stratum we process 1..nparst_loc(istra) particles
  ! so we will test after every nparts_loc(istra)/n_progress_check particles.
  ! Next_check(istra) stores the particle number when the next check is due
  integer, allocatable, dimension(:) :: next_check

  contains

  subroutine calstr_init_progress
    use eirmod_mpi
    integer :: ierr
  end subroutine calstr_init_progress

  subroutine calstr_progress_message(n, istra, nparts_loc)
  ! Performs an MPI_TEST on the non-blocking messages, to ensure that they are
  ! progressing.
  !
  ! The calstr reduction for stratum(k-1) should finish before we start calstr
  ! for stratum(k). Therefore, during the particle loop of stratum(k) we
  ! periodically call MPI_TEST. We aim to call the test subroutine
  ! n_progress_check times. Using the particle loop index n, and the total
  ! number of particles nparts_loc(k) we decide whether to call MPI_TEST or
  ! skip it.
    use eirmod_mpi
    integer, intent(in) :: n !< particle index
    integer, intent(in) :: istra !< stratum number
    !> total number of particles
    integer, intent(in), dimension(:), allocatable :: nparts_loc
    logical :: flag
    integer :: ierr
  end subroutine calstr_progress_message

  subroutine calstr_buffered_finish()
    use eirmod_mpi
    integer :: ierr
  end subroutine calstr_buffered_finish

  subroutine eirene_calstr_buffered(comm, my_pe_gr, leader)
  ! Non-blocking reduction using extra buffers.
  ! The data is first copied to buffers using calstr_pack_all, and the buffers
  ! are reduced using MPI_IREDUCE. The stratum leader waits for the operation
  ! to finish, the others are free to continue with the next stratum.
  !
  ! If you need to add any other reduction operations, then please do so in
  ! allocate_calstr_buffer, calstr_pack_all, calstr_unpack_all.
    use eirmod_mpi
    integer, intent(in) :: comm !< calstr_communicator for the stratum
    integer, intent(in) :: my_pe_gr !< my rank within the calstr communicator
    logical, intent(in) :: leader !< whether the PE is the leader of the stratum
    integer :: ierr, isize

  end subroutine eirene_calstr_buffered

  subroutine calstr_pack_d0(buff)
    real(kind=dp), intent(in) :: buff

  end subroutine calstr_pack_d0

  subroutine calstr_pack_d(n, buff)
    integer :: n
    real(kind=dp), dimension(*), intent(in) :: buff

  end subroutine calstr_pack_d

  subroutine calstr_unpack_d0(buff)
    real(kind=dp), intent(out) :: buff
    buff = 0

  end subroutine calstr_unpack_d0

  subroutine calstr_unpack_d(n, buff)
    integer, intent(in) :: n
    real(kind=dp), dimension(*) :: buff

  end subroutine calstr_unpack_d

  subroutine calstr_pack_l0(buff)
    logical, intent(in) :: buff

  end subroutine calstr_pack_l0

  subroutine calstr_pack_l(n, buff)
    integer, intent(in) :: n
    logical, dimension(*), intent(in) :: buff

  end subroutine calstr_pack_l

  subroutine calstr_unpack_l0(buff)
    logical, intent(out) :: buff
    buff = .FALSE.

  end subroutine calstr_unpack_l0

  subroutine calstr_unpack_l(n, buff)
    integer, intent(in) :: n
    logical, dimension(*), intent(out) :: buff
    buff(1:n) = .FALSE.

  end subroutine

  subroutine allocate_calstr_buffer(ierr)
    integer, intent(out) :: ierr
    integer n, ns, ispc
    integer :: istra
    ierr = 0

    ! variables from call eirene_calstr_usr should be also included, if there are any
  end subroutine allocate_calstr_buffer

  subroutine deallocate_calstr_buffer

  end subroutine deallocate_calstr_buffer

  subroutine calstr_pack_all
    use eirmod_mpi
    integer :: ns, ispc

  end subroutine calstr_pack_all

  subroutine calstr_unpack_all
    use eirmod_mpi
    integer :: ns, ispc
 ! variables from eirene_calstr_usr should be also unpacked
  end subroutine calstr_unpack_all


end module eirmod_calstr_buffered
