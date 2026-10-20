      MODULE EIRMOD_SECOND_OWN
cdr  used for internal run time monitoring

      USE EIRMOD_PRECISION
#ifdef USE_OPENMP
      use omp_lib           !IGNORE
#endif
#ifdef USE_MPI
#ifndef MPI_VERSION
#define MPI_VERSION 3
#endif
#ifndef MPI_MOD
#if MPI_VERSION > 2
#define MPI_MOD 1
#else
#define MPI_MOD 0
#endif
#endif
#if MPI_MOD > 0
      use mpi               !IGNORE
#endif
#endif
      IMPLICIT NONE
#ifdef USE_MPI
#if MPI_MOD < 1
      include 'mpif.h'
#endif
#endif
      PRIVATE

      PUBLIC :: EIRENE_SECOND_OWN, EIRENE_RESET_SECOND

      real(dp), save :: start=0.0_dp

      CONTAINS
c
      FUNCTION EIRENE_SECOND_OWN()
cdr returns wall clock time (in seconds) since the last call to function eirene_reset_second.
      implicit none

      real(dp) :: EIRENE_second_own
      real(dp) :: time
#ifdef USE_OPENMP
#ifdef NAGFOR
      double precision, external :: omp_get_wtime
#endif
      time=real(omp_get_wtime(),DP)
#elif USE_MPI
#if MPI_MOD < 1
      double precision, external :: mpi_wtime
#endif
      time=real(mpi_wtime(),DP)
#else
      call cpu_time(time)
#endif

      EIRENE_second_own=time-start
      return
      END FUNCTION EIRENE_SECOND_OWN
C
      FUNCTION EIRENE_RESET_SECOND()
cdr returns wall clock time (in seconds)
      implicit none
      real(dp) :: EIRENE_reset_second
#ifdef USE_OPENMP
#ifdef NAGFOR
      double precision, external :: omp_get_wtime
#endif
      start=real(omp_get_wtime(),DP)
#elif USE_MPI
#if MPI_MOD < 1
      double precision, external :: mpi_wtime
#endif
      start=real(mpi_wtime(),DP)
#else
      call cpu_time(start)
#endif
      EIRENE_reset_second=start
      return
      END FUNCTION EIRENE_RESET_SECOND

      END MODULE EIRMOD_SECOND_OWN
