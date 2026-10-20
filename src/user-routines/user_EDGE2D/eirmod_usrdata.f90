      MODULE EIRMOD_USRDATA

      USE EIRMOD_PRECISION
      IMPLICIT NONE

      PUBLIC  :: EIRENE_ALLOC_USRDATA, EIRENE_DEALLOC_USRDATA

      !REAL(DP), PUBLIC, ALLOCATABLE :: ARR(:,:)
      INTEGER, PUBLIC, ALLOCATABLE :: SURFTYPE(:)

      integer, public, parameter :: PUFF = 1
      integer, public, parameter :: PUMP = 2
      integer, public, parameter :: CORE = 3
      integer, public, parameter :: TARG = 4
      integer, public, parameter :: MSUR = 5
      integer, public, parameter :: WALL = 6
      integer, public, parameter :: TIME = 0
      
      integer, save :: eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer, eirene_fstoreneutflux, eirene_felemente
      real(dp) :: eirene_phi_offsets(9)
      integer :: nrpla   = 0
      integer :: nbkgusr = 0
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
!                = 0: no wall fluxes are used (old edge2d model)
!                = 1: only ion fluxes are used
!                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
!                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give 
!                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      logical :: eirene_use_e2d_ne = .false.
      logical, public :: lsavetallies = .false.
      integer :: icoarselev = 1

      namelist /eirene_user/eirene_nbirth,eirene_njetto,    &
     &                      eirene_fbirth,eirene_ftransfer, &
     &                      eirene_phi_offsets,             &
     &                      eirene_fstoreneutflux,          &
     &                      eirene_wallFluxModel,           &
     &                      eirene_use_elstepdat_bug,       &
     &                      eirene_use_e2d_ne, icoarselev,  &
     &                      lsavetallies

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_USRDATA(NSURF)
      INTEGER, INTENT(IN) :: NSURF

      IF (ALLOCATED(SURFTYPE)) RETURN
      
      ALLOCATE(SURFTYPE(NSURF))
      SURFTYPE(:) = 0

      END SUBROUTINE EIRENE_ALLOC_USRDATA

      SUBROUTINE EIRENE_DEALLOC_USRDATA
        DEALLOCATE (SURFTYPE)
      END SUBROUTINE EIRENE_DEALLOC_USRDATA

      END MODULE EIRMOD_USRDATA
