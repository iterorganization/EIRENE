cdr  nov. 15:  comments,  irds --> irei
cdr  april 16: added: fail-safe (exit) step in case of more than one (distinct) bulk
cdr            secondaries.
cdr            This is temporarily necessary, as a consequence of making the
cdr            (bulk) ion energy sources eapl, empl, eipl species-dependent
cdr            We are not aware of any application of eirene, in which this new error exit
cdr            would be activated.

      SUBROUTINE EIRENE_SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
CDR  CALLED IN INITIALIZATION PHASE OF RUN
c  ical =   0: allocate storage
c  ical ne. 0: call EIRENE_INIT_CMDTA(2) (contained in eirmod_comxs) cdr: called twice ??
C
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      USE EIRMOD_CGRID, ONLY : NSBOX
      USE EIRMOD_CZT1
      USE EIRMOD_PHOTON

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, J, IRPI, IREI, IERROR
      EXTERNAL :: EIRENE_XSECTA_PARAM, EIRENE_XSECTM_PARAM,
     .            EIRENE_XSECTI_PARAM, EIRENE_XSECTP_PARAM,
     .            EIRENE_XSECTPH_PARAM, EIRENE_CONDENSE,
     .            EIRENE_XSECTA, EIRENE_XSECTM,
     .            EIRENE_XSECTI, EIRENE_XSECTP, EIRENE_XSECTPH,
     .            EIRENE_EXIT_OWN

      IF (ICAL == 0) THEN
        write (iunout,*) 'setamd(0) called '
        NRCX=0
        NREL=0
        NRPI=0
        NREI=0
        NREC=0
        NBGV=0
        NRPH=0
        CALL EIRENE_XSECTA_PARAM
        CALL EIRENE_XSECTM_PARAM
        CALL EIRENE_XSECTI_PARAM
        CALL EIRENE_XSECTP_PARAM
        CALL EIRENE_XSECTPH_PARAM

        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NREI=MAX(1,NREI)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        NRPH=MAX(1,NRPH)

        CALL EIRENE_SET_PARMMOD(2)
        CALL EIRENE_ALLOC_COMXS(2)
        CALL EIRENE_ALLOC_COMSOU(2)
        CALL EIRENE_ALLOC_CZT1(2)

        ! find valid species range via "maxspc(ityp)"
        MAXSPC(0:4) = (/ NPHOTI,NATMI,NMOLI,NIONI,NPLSI /)

        RETURN

      ELSE  ! (ical.ne.0)

cdr Now we have ICAL /=0
        write (iunout,*) 'setamd(ical) ',ical

        CALL EIRENE_INIT_CMDTA(2)

!PB only if storage mode is not switched on
        IF (NSTORDR >= NRAD) THEN

CVK TABLES CHECKING (FOR ELASTIC COLLISIONS)
          DO I=1,NREL
            DO J=1,NSBOX
              IF(TABEL3(I,J,1).GT.23) THEN
                WRITE(iunout,*)
     .             "SETAMD WARNING: REACTION RATE IS TOO BIG",
     .             "IREL,ICELL,TABEL3",I,J,TABEL3(I,J,1)
              END IF
            END DO
          END DO
CVK TABLES CHECKING (FOR CHARGE EXCHANGE)
          DO I=1,NRCX
            DO J=1,NSBOX
              IF(TABCX3(I,J,1).GT.23) THEN
                WRITE(iunout,*)
     .            "SETAMD WARNING: REACTION RATE IS TOO BIG",
     .            "IRCX,ICELL,TABCX",I,J,TABCX3(I,J,1)
              END IF
            END DO
          END DO

        END IF
      END IF

      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRPHI=0             
      NRBGI=0

csw 27jul2011
      NPBGKP=0 !VK
      NPBGKA=0 !VK
      NPBGKM=0 !VK
      NPBGKI=0 !VK
csw

      CALL EIRENE_XSECTA

      CALL EIRENE_XSECTM

      CALL EIRENE_XSECTI

      CALL EIRENE_XSECTP

      CALL EIRENE_XSECTPH

      CALL EIRENE_CONDENSE
c
cdr  set some further assistant arrays, for EI and PI processes:
cdr  Accumulated information from A, M, I, P and PH for particle processes 'EI' and 'PI'.
cdr  These array are stored in comxs and are used for scoring
cdr  tallies in update.f (tracklength) and collide.f (coll. estim) exclusively
cdr  They are for indirect indexing, in loops over secondary species.
cdr  e.g. rather than
cdr                   do iat=1,natmi
cdr         now:
cdr                   do i   =1,ipatei(irei,0)   (<=natmi,  possibly much shorter loop)
cdr                      iat = ipatei(irei,i)    (now we know: iat is a secondary indeed)
cdr                      inum= patei(irei,iat)   (there are inum secondaries of species iat)
cdr                      ...
cdr                   enddo
cdr
      IERROR = 0

      IPATEI = 0
      IPMLEI = 0
      IPIOEI = 0
cdr   IPPHEI = 0   ARRAY IPPHEI IS STILL MISSING, NO PHOTON SECONDARIES IN EI REACTIONS.
      IPPLEI = 0
      DO IREI=1,NREI
cdr amongst all natm species there are ipatei(...,0) (<= natm)
cdr distinct atomic species which appear as secondaries,
cdr with one or more fragments per atomic species iatm
        ipatei(IREI,0)=COUNT(PATEI(IREI,1:) > 0)
        IF (ipatei(IREI,0).GT.0) THEN
             IPATEI(IREI,1:ipatei(IREI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATEI(IREI,1:) > 0)
cdr  IPATEI(IREI,...)=IATM means: one or more secondaries of species IATM
c
cdr  the arrays patei,...,pplei, and p2nd, contain the further information:
cdr  "how many" of secondary species IATM arise after process IREI.
        END IF

cdr  next: same for molecular secondaries
        ipmlei(IREI,0)=COUNT(PMLEI(IREI,1:) > 0)
        IF (ipmlei(IREI,0).GT.0) THEN
             IPMLEI(IREI,1:ipmlei(IREI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLEI(IREI,1:) > 0)
        END IF

cdr  next: same for test ion secondaries
        ipioei(IREI,0)=COUNT(PIOEI(IREI,1:) > 0)
        IF (ipioei(IREI,0).GT.0) THEN
             IPIOEI(IREI,1:ipioei(IREI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOEI(IREI,1:) > 0)
        END IF

cdr  next: same for field particle (bulk) secondaries
        ipplei(IREI,0)=COUNT(PPLEI(IREI,1:) > 0)
        IF (ipplei(IREI,0).GT.0) THEN
             IPPLEI(IREI,1:ipplei(IREI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLEI(IREI,1:) > 0)
        END IF
        if (ipplei(IREI,0) > 1) then
cdr
cdr  This unnecessary constraint arose with extension to species dependent tallies EAPL, EMPL, etc.
cdr  as needed for some EMC3 bgk iterations, around 2016.
cdr  The coding of these species resolved energy source tallies is still incomplete.
          IERROR = IERROR + 1
          write (iunout,*) 'MORE THAN ONE BULK ION SPECIES SPECIFIED ',
     .          'AS SECONDARY PARTICLE OF EI REACTION IREI = ',IREI
        end if
      END DO

cdr   same as above, for PI processes
      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
cdr   IPPHPI = 0   ARRAY IPPHPI IS STILL MISSING, NO PHOTON SECONDARIES IN PI REACTIONS.
      IPPLPI = 0

      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IF (ipatpi(IRPI,0).GT.0) then
          IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                  PATPI(IRPI,1:) > 0)
        endif
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IF (ipmlpi(IRPI,0).GT.0) then
          IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                  PMLPI(IRPI,1:) > 0)
        endif
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IF (ipiopi(IRPI,0).GT.0) then
          IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                  PIOPI(IRPI,1:) > 0)
        endif
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IF (ipplpi(IRPI,0).GT.0) then
          IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                  PPLPI(IRPI,1:) > 0)
        endif
        if (ipplpi(IRPI,0) > 1) then
          IERROR = IERROR + 1
          write (iunout,*) 'MORE THAN ONE BULK ION SPECIES SPECIFIED ',
     .          'AS SECONDARY PARTICLE OF PI REACTION IRPI = ',IRPI
        end if
      END DO

      if (ierror > 0) then
         write (iunout,*) 'only a temporary fail-safe step'
         write (iunout,*) 'contact eirene group at fzj, if this occurs'
         write (iunout,*) 'CALCULATION ABANDONED '
         CALL EIRENE_EXIT_OWN(1)
      end if

      RETURN
      END SUBROUTINE EIRENE_SETAMD
