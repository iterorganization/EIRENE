c     memory allocation for additional processors

!pb  15.05.20  allocation is now done in the broadcast routines
!pb            directly in the modules
!pb            allocations triggered here belong to modules with no
!pb            broadcasting routine

      SUBROUTINE EIRENE_ALLOCATE_MODULES
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      !USE EIRMOD_BRAEIR
      !USE EIRMOD_EIRBRA
      USE EIRMOD_CLAST
      USE EIRMOD_CPLOT
      USE EIRMOD_CSPEI
      USE EIRMOD_CSPEZ
      USE EIRMOD_CUPD
      IMPLICIT NONE
      INTEGER JTRJ

      CALL EIRENE_ALLOC_COMPRT
      CALL EIRENE_ALLOC_CPES(1)
      CALL EIRENE_ALLOC_CLAST
      CALL EIRENE_ALLOC_CPLOT
      CALL EIRENE_ALLOC_CSPEI
      CALL EIRENE_ALLOC_CSPEZ
      CALL EIRENE_ALLOC_CUPD

!  allocate and initialize storage for trajectories
      IF (.NOT.ALLOCATED(TRAJ)) THEN
        ALLOCATE (TRAJ(NCHOR+NTRJ))

        DO JTRJ = 1, NCHOR+NTRJ
          ALLOCATE(TRAJ(JTRJ)%TRJ)
          TRAJ(JTRJ)%TRJ%NCOU_CELL = 0
          NULLIFY(TRAJ(JTRJ)%TRJ%CELLS)
        END DO
      END IF

      RETURN
      END SUBROUTINE EIRENE_ALLOCATE_MODULES
