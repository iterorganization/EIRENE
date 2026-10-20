!pb APR   16:   pplds  -> pplei
!pb APR   16:   eelds  -> eelei
!pb MAY   16:   tabds1 -> tabei1
cdr Nov   16    finalizing notational synchronisation (..DS.. (legacy) --> ..EI..)
cdr Nov.  17:   1) sync with couple_B2 from git repository. done
cdr             2) ESIG array: additional argument IPLS: done.
cdr             3) RTIS% pointer to sploda,.....
cdr             4) rates SEIODA, SEINWA added (were missing, used for ipls total ion energy density)
cdr                now: SEIOD(.., NPLS), SEINW(...,NPLS) added
cdr Dec.  20:   Added LTEST: testing iterations via calling EIRSRT from eirene_main, repeatedly,
cdr             rather than internal eirene iterative mode (for bgk, photons, time-dep)

C  MAIN INTERFACING ROUTINE FOR COUPLED CFD PLASMA - EIRENE APPLICATIONS

C  This routine is called from CFD PLASMA CODE and provides the entry point into EIRENE.

C................................................................................................
cc  ltime = .true.:
C
C   SPECIAL TREATMENT OF FIRST CALL TO EIRENE IN THIS (COUPLED) RUN:
C
C    IFIRST=0  (A RESTART)

C      CALL EIRENE(..)     (main-routines)
C
C      CALL EIRENE_EIRENE(..)     (main-routines), call input,.....,
C                                                  set IITER=ITNR
C
C   LATER CALLS  (IFIRST.GT.1) (INSIDE A B2.5 CYCLING):
C
C      CALL EIRENE_COUPLE  (entry to EIRENE  main-routines, bypassing some
C                           initialization stuff. Currently only in case LTIME=T)

C.................................................................................................


cc   ltime = .false.:   ONLY CALLING: EIRENE_EIRENE(..)     (main-routines), call input,.....,
C                                                           set IITER=ITNR
C                       BOTH FOR IFIRST=0 (RESTART) AND IFIRST.GE.1 (INSIDE A B2 CYCLING)
CDR

C
      SUBROUTINE EIRENE_EIRSRT(LSTOP_in,LTIME_in,DELTAT_in,FLUXES_in,
     .                  B2BRM_in,B2RD_in,B2Q_in,B2VP_in,STEP_CPU_in)

C   INPUT:
C     LSTOP:
C     LTIME: TIME-DEPENDENT MODE. PREPARE TIME-DEPENDENT OPTIONS,
C            AND THEN CALL EIRENE
C     DELTAT: TIME STEP  (IRRELEVANT IN CASE LTIME=.FALSE.)
C     FLUXES(1:NSTRA): Scaling factors for source strength in external code,
cdr                    not to be confused with: FLUX(1: NSTRA), the source strength
cdr                    used in eirene run.
C     STEP_CPU:  if gt. 0.0: alter eirene cpu time (NTCPU), seconds
C
C   ONLY FOR EIRENE ENERGY BALANCE DIAGNOSTICS:
C     B2BRM: TOTAL BREMSSTAHLUNG LOSS IN PREVIOUS B2 STEP
C     B2RD : TOTAL (LINE) RADIATION LOSS IN PREVIOUS B2 STEP
C     B2Q  : VOLUMETRIC ENERGY EXCHANGE (ELECTRONS-IONS) DUE TO COULOMB INTERACTION
C     B2VP : VOLUMETRIC ENERGY EXCHANGE (ELECTRONS-IONS) DUE TO WORK DONE BY ELECTRIC FIELD

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_BRASCL
csw mpi
      use eirmod_eirbra
csw
csw 12apr2011
      use eirmod_ctrig
csw
      use eirmod_extrab25
      use eirmod_infcop
      use eirmod_mpi
      USE EIRMOD_SECOND_OWN
      USE EIRMOD_EIRENE

      IMPLICIT NONE
c
      REAL(DP), INTENT(IN) :: FLUXES_in(NSTRA)
      REAL(DP), INTENT(IN) :: DELTAT_in, B2BRM_in, B2RD_in, B2Q_in,
     .                        B2VP_in,STEP_CPU_in
      LOGICAL, INTENT(IN) :: LSTOP_in, LTIME_in
      integer :: rank_mpi,ierr_mpi,size_mpi
C
      REAL(DP), ALLOCATABLE :: FLUXES(:), FLUXS(:)
      REAL(DP) :: DELTAT, B2BRM, B2RD, B2Q, B2VP, STEP_CPU
      LOGICAL :: LSTOP, LTIME
      logical :: ltrigger

      REAL(DP) :: EIRENE_FTABEI1, EIRENE_FEELEI1, ESIG,
     .            DUMMY, DTIMVO, EN
      INTEGER :: IN, IAEI, IMEI, IIEI, IREI, IFIRST, NDXY,
     .           ITNR, IPLSTI, IST_RATE, IST,
     .           JATM, JMOL, JION, JPLS, IENTRY,
     .           IFRSTR, ISTH, ISTNEW, ISTIN, ISTRAI, ifrst_save, iflg
      LOGICAL :: LSTP, LLST, LPLASM
      LOGICAL :: NLSRON_SAVE(NSTRA), LOGHELP(NSTRA)
      EXTERNAL :: EIRENE_FTABEI1, EIRENE_FEELEI1,
     .            EIRENE_HEADNG, EIRENE_LEER,
     .            EIRENE_MASAGE, EIRENE_MASR1,
     .            EIRENE_PLASMA, EIRENE_PLASMA_DERIV,
     .            EIRENE_PLSTRT, EIRENE_PLEND,
     .            EIRENE_REPLACE_STORED, EIRENE_RPLAM, EIRENE_SETAMD,
     .            EIRENE_WRPLAM, EIRENE_EXIT_OWN
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif
C

cdr pointer :: rtis, FOR "OD" RATES (SHORT CYCLE) PER STRATUM
      TYPE(RATE_STORE), POINTER :: RTIS
C
C
      SAVE
      DATA IFIRST/0/, ifrst_save/0/
csw 28feb2011 mpi
      if(.not.allocated(fluxes)) allocate(fluxes(nstra),fluxs(nstra))
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi == 0) then
        lstop=lstop_in
        ltime=ltime_in
        deltat=deltat_in
        fluxes(1:nstra) = fluxes_in(1:nstra)
        b2brm=b2brm_in
        b2rd=b2rd_in
        b2q=b2q_in
        b2vp=b2vp_in
        step_cpu=step_cpu_in

        ltrigger=.true.
        call mpi_bcast(ltrigger,1,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
      endif
      call mpi_bcast(lstop,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
      call mpi_bcast(ltime,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
      call mpi_bcast(deltat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(fluxes,nstra,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ierr_mpi)
      call mpi_bcast(b2brm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(b2rd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(b2q,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(b2vp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
      call mpi_bcast(step_cpu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     .               ierr_mpi)
csw
C
      IF (LTIME) THEN
csw mpi
        CALL EIRENE_MASAGE ('NO MPI FOR LTIME=.true.')
        CALL EIRENE_EXIT_OWN(1)
csw
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP

        DUMMY=EIRENE_RESET_SECOND()

        IF(IFIRST.EQ.0) THEN
C
          CALL EIRENE_PLSTRT
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION (CENSUS) ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM.
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
          LPLASM=.FALSE.
          LLST=.FALSE.
          ITNR=1
C

c          CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.FALSE.)
csw --> MPI_INIT=.false., NLPLAS=.TRUE.

          CALL EIRENE_EIRENE(DELTAT,.TRUE.,LLST,ITNR,.FALSE.)
C
C  EIRENE RUN DONE. CENSUS ARRAY WRITTEN
C  NOW ITIMV=ITIMV+1, NLPLAS=.TRUE.
C
          IF (.NOT.NLPLAS) THEN
            WRITE (iunout,*) 'INCONSISTENT COUPLING'
            WRITE (iunout,*) 'LTIME=TRUE, BUT NTIME = ', NTIME
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          DO 3 ISTRAI=1,NSTRAI
            FLUXS(ISTRAI)=FLUX(ISTRAI)
    3     CONTINUE
          IFIRST=1

        ELSE  !(IFIRST.GE.1)
C  THIS IS NOT THE FIRST CALL TO EIRENE, and we are in time dep. mode.
C
C  NOW: NLPLAS=.TRUE., I.E., PLASMA DATA EXPECTED ON BRAEIR
C  NOW: ITIMV=ITIMV+1
C  BUT: COMMON BRAEIR REDONE IN EXTERNAL CODE.
C  REACTIVATE INDEX MAPPING, EVEN WITHOUT READING INPUT BLOCK 14 AGAIN
C
          DTIMVO=DTIMV
          IF (DTIMVN.NE.DELTAT.AND.DELTAT.GT.0.0_DP) THEN
            DTIMVN=DELTAT
            WRITE (iunout,'(a,1pe12.5,a,a)')
     .       ' EIRENE TIMESTEP CHANGED TO ',DTIMVN,
     .       ' PLASMA SECONDS BY EIRENE_STEP_DT FROM ',
     .       'b2.neutrals.parameters'
          END IF
          IF (STEP_CPU.GT.0.0_DP) THEN
            NTCPU = nint(STEP_CPU)
            WRITE (iunout,'(a,i6,a,a)')
     .       ' EIRENE RUN TIME CHANGED TO ',NTCPU,
     .       ' CPU SECONDS BY EIRENE_STEP_CPU FROM ',
     .       'b2.neutrals.parameters'
          END IF
C
C-----------------------------------------------------------------------
C
C  STRATA 1 TO NTARGI ARE SCALED IN PLASMA CODE (RECYCLING STRATA)
C
C     RETURN TO PLASMA CODE THE PROFILES PER UNIT SOURCE STRENGTH
C     IE. THE PROFILES ARE SCALED BY 1./FLUX(ISTRA) BEFORE RETURN
C
C  STRATA NTARGI+1 TO NSTRAI-1 ARE SCALED BY EIRENE
C
C     (EG. GAS PUFF, VOLUME RECOMBINATION, ETC.)
C     THEY MAY BE RESCALED BY PLASMA CODE FACTORS: FLUXES(ISTRA)
C     RETURN TO PLASMA CODE THE PROFILES SCALED WITH
C     SOURCE STRENGTH: FLUX(ISTRA) (AMP)
C
C  STRATUM NSTRAI IS RESCALED WITH RATIO OF OLD TO NEW TIMESTEP
C
C     RETURN TO PLASMA CODE THE PROFILES WITH FLUX(ISTRA) (AMP)
C
          DO ISTRAI=NTARGI+1,NSTRAI-1
            ISTRA = ISTRAI
            IF (FLUXES(ISTRA).NE.0.) THEN
              FLUX(ISTRA)=FLUXS(ISTRA)*FLUXES(ISTRA)*ELCHA
            ELSE
              FLUX(ISTRA)=FLUXS(ISTRA)
            ENDIF
          ENDDO
C  ISTRA=NSTRAI: RESCALE CENSUS ARRAY FLUX, DUE TO DIFFERENT TIME STEPS
C                IN PREVIOUS AND CURRENT EIRENE STEP
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
            FLXCEN=FLUX(NSTRAI)
C
            WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW'
            CALL EIRENE_MASR1('FLUX    ',FLUX(NSTRAI))
            CALL EIRENE_LEER(1)
          ENDIF
C
C-----------------------------------------------------------------------
C
          DTIMV=DTIMVN
C
C  RUN EIRENE ON TIMESTEP DTIMV
C  THEN CALL INTERFACING ROUTINE AT ENTRY IF3COP (FROM EIRENE MAIN)
C
          IITER=1
          IPRNLI=0
          NLSRON(1:NSTRAI)=.TRUE.
          ITNR=1
          CALL EIRENE_EIRENE_COUPLE (LSTOP,ITNR,.TRUE.)
          IF (LSTOP) THEN
            CALL EIRENE_PLEND
          ENDIF
        ENDIF
        if (rank_mpi==0) then
          CALL EIRENE_LEER(2)
          WRITE(iunout,*)
     .     'EIRENE USED ',EIRENE_SECOND_OWN(),' CPU SECONDS'
          CALL EIRENE_LEER(2)
        end if
C
        RETURN   ! time dep. option done
c.....................................................................
C
C  MAIN ENTRY POINT FROM B2 INTO EIRENE, IN CASE OF TIME-INDEPENDENT RUNS
C
      ELSEIF (.NOT.LTIME) THEN

!swpb for multiprocessor calculation
        DUMMY=EIRENE_RESET_SECOND()
C
        IF (IFIRST.GE.1) GOTO 10000
C
C  FIRST CALL IN PRESENT RUN.
C  1) INITIALIZE EIRENE
C  2) CALL EIRENE
C  3) PREPARE ARRAYS FOR SEMI-IMPLICIT "SHORT CYCLE" CORRECTION.
C             STORE SOME A&M RATES FROM PRESENT STEP, FOR NEXT STEP

csw mpi
        if(rank_mpi .eq. 0) then
csw
          CALL EIRENE_PLSTRT
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
          B2BREM=B2BRM
          B2RAD=B2RD
          B2QIE=B2Q
          B2VDP=B2VP
          LPLASM=.FALSE.
csw 9sep2011
          LPLASM=.TRUE.
csw
          LLST=LSTOP
          ITNR=1
csw mpi
        endif
csw

   10   CONTINUE  ! from here on:
                  ! both ifirst=0 and ifirst.ge.1 are possible

csw 17feb2011
c        if(itnr > 1 .and. rank_mpi==0) then
c          DO ISTRAI=NTARGI+1,NSTRAI
c            ISTRA = ISTRAI
c            if(.not.nlvol(istra)) then
c              IF (FLUXES(ISTRA).NE.0.d0 ) THEN
c                FLUX(ISTRA)=FLUXES(ISTRA)*ELCHA
c              ELSE
c                FLUX(ISTRA)=1.
c              ENDIF
c            endif
c          ENDDO
c        endif
csw 24oct2011

        if(rank_mpi==0) then
          DO ISTRAI=1,NSTRA
            FLUX_save(ISTRAI)=FLUXES(ISTRAI)
          ENDDO
        endif


!pb     CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.TRUE.)

        CALL EIRENE_EIRENE(DELTAT,LPLASM,LLST,ITNR,.FALSE.)
        if (ifrst_save == 0) then
          if (count(nlsron(1:nstrai-1)) > count(nlvol(1:nstrai-1))) then
            ifrst_save = 1
          else
            if (nprnli == 0) then
              nlsron(1:nstrai) = .true.
            else
            nlsron(1:nstrai-1) = .true.
            end if
          end if
        end if
C
C  IN THIS CALL TO EIRENE ALREADY IF3COP IS CALLED FOR EACH STRATUM
C  THOSE WITH NLSRON(ISTRA) = TRUE  HAVE BEEN RECOMPUTED BY EIRENE
C  THOSE WITH NLSRON(ISTRA) = FALSE HAVE BEEN SHORT-CYCLED
C  AT IFIRST.EQ.0: ALL NLSRON=TRUE
C  AT IFIRST.GE.1: FIRST A SHORT CYCLE TEST IS DONE, AND NLSRON IS FOUND
C
csw mpi
        if(rank_mpi .eq. 0) then

        IF (.NOT.LLST) THEN
cdr  prepare some stuff for next iteration step within this same run,
cdr  such as the ..OD.. parameters for short cycling

        IF (IFIRST.GE.1) NLSRON = NLSRON_SAVE

csw     NDXY=(NDXA-1)*NR1ST+NDYA  ! this version: triang. grids only
        NDXY=NTRII
C
CDR why alloc short cycle data, even if no short cycle is done ??
        CALL EIRENE_ALLOC_BRASCL

! find new calculated stratum with smallest number
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            IFRSTR = IST
            EXIT
          END IF
        END DO

! determine index of rate storage which has been used in the previous iteration
        IST_RATE = ITS(IFRSTR)

! reduce counters of rate storages for all new calculated strata
        DO IST = 1, NSTRAI
          IF (NLSRON(IST)) THEN
            ISTIN = ITS(IST)
            ITS_COUNT(ISTIN) = ITS_COUNT(ISTIN) - 1
          END IF
        END DO

! check how often storage IST_RATE is still used
        ISTH = 0
        IF (IST_RATE > 0) ISTH = ITS_COUNT(IST_RATE)

        IF (ISTH < 1) THEN
! rate storage can be used again
        ELSE
! rate storage still in use, look for an empty slot
          ISTNEW = MINLOC(ITS_COUNT,DIM=1)
          IF (ITS_COUNT(ISTNEW) > 0) THEN
            WRITE (IUNOUT,*) ' PROBLEM IN EIRSRT'
            WRITE (IUNOUT,*) ' ITS_COUNT > 0'
            WRITE (IUNOUT,*) ' ITS_COUNT ',ITS_COUNT
            CALL EIRENE_EXIT_OWN(1)
          END IF
          IST_RATE = ISTNEW
        END IF

        WHERE (NLSRON(1:NSTRAI))
          ITS(1:NSTRAI) = IST_RATE
        END WHERE

        ITS_COUNT(IST_RATE) = COUNT(NLSRON(1:NSTRAI))
C
        CALL EIRENE_ALLOC_RATE_ARRAY(IST_RATE)
        CALL EIRENE_INIT_RATE_ARRAY(IST_RATE)

        RTIS => RTS(IST_RATE)%RTA
C.....................................................................................

cdr
cdr  now start to store rates from present cycle, for future short cycle corrections
cdr
cdr  to be done
cdr  all these "short cycle data" should only be computed
cdr  if "short cycle" option is turned on at all

C
C  CURRENT RUN: STORE ION ENERGY DENSITY: FOR ALL IPLS, BUT TIIN(IPLS) MAY BE THE SAME FOR
C               SOME OR ALL IPLS
C
C
          DO JPLS=1,NPLSI
            IPLSTI= MPLSTI(JPLS)
            DO IN=1,NDXY
!pb           RTIS%SEIOD(IN,JPLS)=DIIN(JPLS,IN)*
!pb     .                    (1.5*TIIN(IPLSTI,IN)+EDRIFT(JPLS,IN))
              EN = 1.5*TIIN(IPLSTI,IN)
              IF (LEDRIFT) EN = EN + EDRIFT(JPLS,IN)
              RTIS%SEIOD(IN,JPLS)=DIIN(JPLS,IN)*EN
            ENDDO
          ENDDO
C
C  NEXT: ATOMS, EI RATES: SPLODA, SEEODA, SEIODA
C
cdr  correct energy exchange with bulk ions: e0* eplei(irei,ipls,1)+ eheavy* eplei(irei,ipls,2)
cdr  sum over ipls:                          e0* eplei(irei,0,1)   + eheavy *eplei(irei,0,2)
cdr  e0 is taken as center of mass (COM)energy (as appropriate in EI processes, but not in PI processes)
cdr  and eheavy is the kinetic energy release (KER) in reaction irei
cdr  the present short cycle correction only accounts for the KER (=0 for atoms), not for the COM part
C
C  CURRENT RUN: PARTICLE RATE: ATOMS, EI PROCESSES, FROM IATM TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
C
        DO JATM=1,NATMI
         DO JPLS=1,NPLSI
          DO 21 IAEI=1,NAEII(JATM)
           IREI=LGAEI(JATM,IAEI)
           IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 21
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODA(IN,JATM,JPLS)=RTIS%SPLODA(IN,JATM,JPLS)+
     .                        TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ELSE
              RTIS%SPLODA(IN,JATM,JPLS)=RTIS%SPLODA(IN,JATM,JPLS)+
     .                        EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            END IF
           ENDDO
   21     CONTINUE
         ENDDO
        ENDDO
C
C
C  CURRENT RUN: ELECTRON COOLING RATE: ATOMS, EI PROCESSES, FROM IATM,
C                                      SUM OVER ALL EI PROCESSES
C
        DO JATM=1,NATMI
         DO IAEI=1,NAEII(JATM)
          IREI=LGAEI(JATM,IAEI)
          DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
                RTIS%SEEODA(IN,JATM)=RTIS%SEEODA(IN,JATM)+
     .                                        EELEI1(IREI,IN)*
     .                                        TABEI1(IREI,IN)
            ELSE
              RTIS%SEEODA(IN,JATM)=RTIS%SEEODA(IN,JATM)+
     .                                        EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
            END IF
          ENDDO
         ENDDO
        ENDDO
C
C
C
C  CURRENT RUN: ION ENERGY EXCHANGE RATE: ATOMS, EI PROCESSES, FROM IATM
C                                         SUM OVER ALL EI PROCESSES
C                                         SUM OVER ALL IPLS
C
        DO JATM=1,NATMI
         DO JPLS=1,NPLSI
          DO IAEI=1,NAEII(JATM)
           IREI=LGAEI(JATM,IAEI)
!pb 09022016            ESIG=EPLEI(IREI,2)  this was incorrect,
cdr                     because it was already summed over ipls
           ESIG=EPLEI(IREI,JPLS,2)  ! only KER part is corrected
                                    ! in short cycle
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
                RTIS%SEIODA(IN,JATM)=RTIS%SEIODA(IN,JATM)+
     .                  TABEI1(IREI,IN)*ESIG
            ELSE
                RTIS%SEIODA(IN,JATM)=RTIS%SEIODA(IN,JATM)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
c....................................................................
C
C  NEXT: TEST IONS, EI RATES: SPLODI, SEEODI, SEIODI

C
C  CURRENT RUN: PARTICLE RATE: TEST IONS, EI PROCESSES, FROM IION TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
C
        DO JION=1,NIONI
         DO JPLS=1,NPLSI
          DO 27 IIEI=1,NIEII(JION)
           IREI=LGIEI(JION,IIEI)
c
           IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 27
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODI(IN,JION,JPLS)=RTIS%SPLODI(IN,JION,JPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ELSE
              RTIS%SPLODI(IN,JION,JPLS)=RTIS%SPLODI(IN,JION,JPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ENDIF
           ENDDO
   27     CONTINUE
         END DO
        END DO
C
C  CURRENT RUN: ELECTRON COOLING RATE: TEST IONS, EI PROCESSES, FROM IION,
C                                      SUM OVER ALL EI PROCESSES
        DO JION=1,NIONI
          DO IIEI=1,NIEII(JION)
            IREI=LGIEI(JION,IIEI)
            DO IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                ESIG=EELEI1(IREI,IN)
                RTIS%SEEODI(IN,JION)=RTIS%SEEODI(IN,JION)+
     .                         TABEI1(IREI,IN)*ESIG
              ELSE
                ESIG=EIRENE_FEELEI1(IREI,IN)
                RTIS%SEEODI(IN,JION)=RTIS%SEEODI(IN,JION)+
     .                         EIRENE_FTABEI1(IREI,IN)*ESIG
              END IF
            ENDDO
          ENDDO
        ENDDO
C
C  CURRENT RUN: ION ENERGY EXCHANGE RATE: TEST IONS, EI PROCESSES, FROM IION
C                                         SUM OVER ALL EI PROCESSES
C                                         SUM OVER ALL IPLS
C
        DO JION=1,NIONI
         DO JPLS=1,NPLSI
          DO IIEI=1,NIEII(JION)
           IREI=LGIEI(JION,IIEI)
!pb 09022016            ESIG=EPLEI(IREI,2)  this was incorrect,
cdr                     because it was already summed over ipls
           ESIG=EPLEI(IREI,JPLS,2)  ! only KER part is corrected
                                    ! in short cycle
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODI(IN,JION)=RTIS%SEIODI(IN,JION)+
     .                        TABEI1(IREI,IN)*ESIG
            ELSE
              RTIS%SEIODI(IN,JION)=RTIS%SEIODI(IN,JION)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
           ENDDO
          END DO
         END DO
        END DO
C
C
C  NEXT: MOLECULES, EI RATES: SPLODM, SEEODM, SEIODM
C
        DO JMOL=1,NMOLI
         DO JPLS=1,NPLSI
          DO 47 IMEI=1,NMEII(JMOL)
           IREI=LGMEI(JMOL,IMEI)
           IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 47
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SPLODM(IN,JMOL,JPLS)=RTIS%SPLODM(IN,JMOL,JPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ELSE
              RTIS%SPLODM(IN,JMOL,JPLS)=RTIS%SPLODM(IN,JMOL,JPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            END IF
           ENDDO
   47     CONTINUE
         END DO
        END DO
C
C
C  CURRENT RUN: ELECTRON COOLING RATE: MOLECULES, EI PROCESSES, FROM IMOL,
C                                      SUM OVER ALL EI PROCESSES
C
        DO JMOL=1,NMOLI
         DO IMEI=1,NMEII(JMOL)
          IREI=LGMEI(JMOL,IMEI)
          DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
                RTIS%SEEODM(IN,JMOL)=RTIS%SEEODM(IN,JMOL)+
     .                                          EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
            ELSE
              RTIS%SEEODM(IN,JMOL)=RTIS%SEEODM(IN,JMOL)+
     .                                        EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
            END IF
          END DO
         END DO
        END DO
C
C  CURRENT RUN: ION ENERGY EXCHANGE RATE: MOLECULES, EI PROCESSES, FROM IMOL
C                                         SUM OVER ALL EI PROCESSES
C                                         SUM OVER ALL IPLS
        DO JMOL=1,NMOLI
         DO JPLS=1,NPLSI
          DO IMEI=1,NMEII(JMOL)
           IREI=LGMEI(JMOL,IMEI)
!pb 09022106         ESIG=EPLEI(IREI,2)  this was incorrect,
cdr                     because it was already summed over ipls
           ESIG=EPLEI(IREI,JPLS,2) ! only KER part is corrected
                                   ! in short cycle
           DO IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              RTIS%SEIODM(IN,JMOL)=RTIS%SEIODM(IN,JMOL)+
     .                        TABEI1(IREI,IN)*ESIG
            ELSE
              RTIS%SEIODM(IN,JMOL)=RTIS%SEIODM(IN,JMOL)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
           END DO
          END DO
         END DO
        END DO

        END IF   ! (not llst)
C
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
C
csw mpi
        endif
csw
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_RATE_ARRAY
          CALL EIRENE_DEALLOC_BRASPOI
          DEALLOCATE(FLUXES,FLUXS)
          CALL EIRENE_PLEND
        END IF
        IF (rank_mpi==0) THEN
          CALL EIRENE_LEER(2)
          WRITE(iunout,*)
     .     'EIRENE USED ',EIRENE_SECOND_OWN(),' CPU SECONDS'
          CALL EIRENE_LEER(2)
        END IF
        RETURN
C
C  NOT THE FIRST CALL IN THIS CYCLE: CHECK: SHORT LOOP CORRECTION
C                                           OR FULL EIRENE, FOR EACH
C                                           STRATUM INDIVIDUALLY
C.......................................................................

10000   CONTINUE ! IFIRST.GE 1 BRANCH
C
        if(rank_mpi .eq. 0) then

        CALL EIRENE_HEADNG ('NEXT EIRENE RUN STARTS HERE',27)

        LSTP = LSTOP
        IF (DTIMVN.NE.DELTAT.AND.DELTAT.GT.0.0_DP) THEN
          DTIMVN=DELTAT
          WRITE (iunout,'(a,1pe12.5,a,a)')
     .     ' EIRENE TIMESTEP CHANGED TO ', DTIMVN,
     .     ' PLASMA SECONDS BY EIRENE_STEP_DT FROM',
     .     ' b2.neutrals.parameters'
          CALL EIRENE_LEER(1)
        END IF
        IF (STEP_CPU.GT.0.0_DP) THEN
          NTCPU = nint(STEP_CPU)
          WRITE (iunout,'(a,i6,a,a)')
     .     ' EIRENE RUN TIME CHANGED TO ', NTCPU,
     .     ' CPU SECONDS BY EIRENE_STEP_CPU FROM',
     .     ' b2.neutrals.parameters'
          CALL EIRENE_LEER(1)
        END IF

        CALL EIRENE_ALLOC_BCKGRND

        IF ((NBGK > 0) .AND. (NFILEL.EQ.3)) then
          IFLG=0
          CALL EIRENE_RPLAM(TRCFLE,IFLG,'EIRSRT')
        END IF

        IENTRY=1
        CALL EIRENE_IF1COP(IENTRY)
C
        CALL EIRENE_PLASMA
C
        CALL EIRENE_PLASMA_DERIV(0,0)
C
        CALL EIRENE_SETAMD(2)

        if (nbgk > 0) call eirene_replace_stored
C
C  IN PLASMA_DERIV THE BACKGROUND PLASMA STATE HAS BEEN
C  WRITTEN TO FORT.13
C  NFILEL HAS BEEN CHANGED TO NFILEL = 3
C  ==> PLASMA AND REACTION DATA ARE READ IN SUBR. INPUT
C  NOW SAVE REACTION DATA AS WELL IN ORDER TO HAVE A
C  CONSISTENT PLASMA STATE ON FORT.13
C
        IF ((NFILEL >=1) .AND. (NFILEL <=5)) THEN
          NFILEL=3
          CALL EIRENE_WRPLAM(TRCFLE,'EIRSRT')
        END IF

csw force no short cycle
csw        IF ( ANY(XMCP_OLD(1:NSTRAI) <= 2._DP)) THEN
c        if(.true.) then
c           IFIRST=0
c           LPLASM=.TRUE.
c           LSTP=LSTOP
c           ITNR=ITNR+1
c           GOTO 10
c        END IF
csw
C
        CALL EIRENE_ALLOC_BRASCL
        CALL EIRENE_INIT_BRASCL

C
cdr
cdr  now start to store rates FOR NEXT cycle, for short cycle corrections
cdr
cdr  to be done
cdr  all these "short cycle data" should only be computed if short cycle is turned on at all
C
C  NEW RUN: ION ENERGY DENSITY: FOR ALL IPLS, BUT TIIN(IPLS) MAY BE THE SAME FOR ALL IPLS
C
        DO JPLS=1,NPLSI
          IPLSTI= MPLSTI(JPLS)
          DO IN=1,NDXY
!pb            SEINW(IN,JPLS)=DIIN(JPLS,IN)*
!pb     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(JPLS,IN))
            EN = 1.5*TIIN(IPLSTI,IN)
            IF (LEDRIFT) EN = EN + EDRIFT(JPLS,IN)
            SEINW(IN,JPLS)=DIIN(JPLS,IN)*EN
          ENDDO
        ENDDO
C
C
C  NEXT: ATOMS, EI RATES: SPLNWA, SEENWA, SEINWA
C
cdr  correct energy exchange with bulk ions: e0*eplei(IREI,ipls,1) + eplei(IREI,ipls,2)
cdr  sum over ipls:                          e0*eplei(IREI,0,1)    + eplei(IREI,0,2)
cdr  The present short cycle correction only accounts for the KER (=0 for atoms)
C
C  NEXT RUN: PARTICLE RATE: ATOMS, EI PROCESSES, FROM IATM TO IPLS,
C                                     SUM OVER ALL EI PROCESSES
        DO JATM=1,NATMI
         DO 101 JPLS=1,NPLSI
          DO IAEI=1,NAEII(JATM)
            IREI=LGAEI(JATM,IAEI)
            IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 101
            DO 102 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SPLNWA(IN,JATM,JPLS)=SPLNWA(IN,JATM,JPLS)+
     .                          TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
              ELSE
                SPLNWA(IN,JATM,JPLS)=SPLNWA(IN,JATM,JPLS)+
     .                          EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
              END IF
  102       CONTINUE
          END DO
  101    END DO
        END DO
C
        DO JATM=1,NATMI
          DO IAEI=1,NAEII(JATM)
            IREI=LGAEI(JATM,IAEI)
            DO 105 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWA(IN,JATM)=SEENWA(IN,JATM)+EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
              ELSE
                SEENWA(IN,JATM)=SEENWA(IN,JATM)+EIRENE_FEELEI1(IREI,IN)*
     .                                          EIRENE_FTABEI1(IREI,IN)
              END IF
  105       CONTINUE
          END DO
        END DO

cdr  There is no seinwa, because currently only KER part is in short cycle correction for EI processes
cdr             and for atoms this is identically 0.0
C
C  NEW: TEST IONS, EI PROCESSES
C
        DO JION=1,NIONI
         DO JPLS=1,NPLSI
          DO 107 IIEI=1,NIEII(JION)
           IREI=LGIEI(JION,IIEI)
           IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 107
           DO 108 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWI(IN,JION,JPLS)=SPLNWI(IN,JION,JPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ELSE
              SPLNWI(IN,JION,JPLS)=SPLNWI(IN,JION,JPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            END IF
  108      CONTINUE
  107     CONTINUE
         END DO
        END DO
C
C
        DO JION=1,NIONI
          DO IIEI=1,NIEII(JION)
            IREI=LGIEI(JION,IIEI)
            DO 106 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWI(IN,JION)=SEENWI(IN,JION)+EELEI1(IREI,IN)*
     .                                          TABEI1(IREI,IN)
              ELSE
                SEENWI(IN,JION)=SEENWI(IN,JION)+EIRENE_FEELEI1(IREI,IN)*
     .                                          EIRENE_FTABEI1(IREI,IN)
              END IF
  106       CONTINUE
          END DO
        END DO
C
        DO JION=1,NIONI
         DO JPLS=1,NPLSI
          DO IIEI=1,NIEII(JION)
           IREI=LGIEI(JION,IIEI)
!pb 09022016          ESIG=EPLEI(IREI,2)
cdr  correct energy exchange with bulk ions: e0* eplei(IREI,ipls,1)+ eplei(IREI,ipls,2)
cdr  sum over ipls:                          e0* eplei(IREI,0,1)   + eplei(IREI,0,2)
cdr  the present short cycle correction only accounts for the KER (=0 for atoms)
           ESIG=EPLEI(IREI,JPLS,2)
           DO 110 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWI(IN,JION)=SEINWI(IN,JION)+TABEI1(IREI,IN)*ESIG
            ELSE
              SEINWI(IN,JION)=SEINWI(IN,JION)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
  110      CONTINUE
          END DO
         END DO
        END DO
C
C  NEW: MOLECULES, EI PROCESSES
C
        DO JMOL=1,NMOLI
         DO JPLS=1,NPLSI
          DO 117 IMEI=1,NMEII(JMOL)
           IREI=LGMEI(JMOL,IMEI)
           IF (PPLEI(IREI,JPLS).EQ.0.) GOTO 117
           DO 118 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWM(IN,JMOL,JPLS)=SPLNWM(IN,JMOL,JPLS)+
     .                             TABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            ELSE
              SPLNWM(IN,JMOL,JPLS)=SPLNWM(IN,JMOL,JPLS)+
     .                         EIRENE_FTABEI1(IREI,IN)*PPLEI(IREI,JPLS)
            END IF
  118      CONTINUE
  117     CONTINUE
         END DO
        END DO
C
C
        DO JMOL=1,NMOLI
         DO IMEI=1,NMEII(JMOL)
          IREI=LGMEI(JMOL,IMEI)
          DO 116 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEENWM(IN,JMOL)=SEENWM(IN,JMOL)+EELEI1(IREI,IN)*
     .                                        TABEI1(IREI,IN)
            ELSE
              SEENWM(IN,JMOL)=SEENWM(IN,JMOL)+EIRENE_FEELEI1(IREI,IN)*
     .                                        EIRENE_FTABEI1(IREI,IN)
            END IF
  116     CONTINUE
         END DO
        END DO

        DO JMOL=1,NMOLI
         DO JPLS=1,NPLSI
          DO IMEI=1,NMEII(JMOL)
           IREI=LGMEI(JMOL,IMEI)
           ESIG=EPLEI(IREI,JPLS,2)
           DO 120 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWM(IN,JMOL)=SEINWM(IN,JMOL)+TABEI1(IREI,IN)*ESIG
            ELSE
              SEINWM(IN,JMOL)=SEINWM(IN,JMOL)+
     .                        EIRENE_FTABEI1(IREI,IN)*ESIG
            END IF
  120      CONTINUE
          END DO
         END DO
        END DO
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        IENTRY=1
        CALL EIRENE_IF3COP(IENTRY,LSTP,IFIRST,1,NSTRAI,0)

        NLSRON_SAVE = NLSRON

        endif ! rank_mpi == 0

!pb copy NLSRON to LOGHELP to avoid warnings from Intel compiler
        LOGHELP(1:NSTRAI) = NLSRON(1:NSTRAI)
        call mpi_bcast(loghelp,nstrai,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
        NLSRON(1:NSTRAI) = LOGHELP(1:NSTRAI)
        call eirene_broadcast_eirbra
csw
        IF (LSTP.or.ANY(NLSRON(1:NSTRAI))) THEN
           IFIRST=0
           LPLASM=.TRUE.
           LSTP=LSTOP
           ITNR=ITNR+1
           GOTO 10
        END IF
C
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_BRASCL
          CALL EIRENE_DEALLOC_RATE_ARRAY
          CALL EIRENE_DEALLOC_BRASPOI
          DEALLOCATE(FLUXES,FLUXS)
          CALL EIRENE_PLEND
        END IF
        IF (rank_mpi==0) THEN
          CALL EIRENE_LEER(2)
          WRITE(iunout,*)
     .     'EIRENE USED ',EIRENE_SECOND_OWN(),' CPU SECONDS'
          CALL EIRENE_LEER(2)
        END IF

        RETURN

CDR     ENDIF   !(IFIRST=1, IFIRST.GE.1 BRANCHING)
C
      ENDIF  !(LTIME)

      RETURN
      END SUBROUTINE EIRENE_EIRSRT

      SUBROUTINE EIRENE_EIRSRT_BROAD
      use eirmod_parmmod
      use eirmod_comsou
      use eirmod_eirbra
      use eirmod_mpi
      implicit none
      integer :: ierr_mpi
      LOGICAL :: LOGHELP(NSTRA)
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

!pb copy NLSRON to LOGHELP to avoid warnings from Intel compiler
      LOGHELP(1:NSTRAI) = NLSRON(1:NSTRAI)
      call mpi_bcast(loghelp,nstrai,MPI_LOGICAL,
     .                 0,MPI_COMM_WORLD,ierr_mpi)
      NLSRON(1:NSTRAI) = LOGHELP(1:NSTRAI)
      call eirene_broadcast_eirbra
      return
      END SUBROUTINE EIRENE_EIRSRT_BROAD
