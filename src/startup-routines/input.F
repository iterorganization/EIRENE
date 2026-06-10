cdr             call eirene_skip_read_comments: not fully implemented, more places
cdr Nov. 19 :   Clearer messages re BGK iterations
cdr             Started to bring back surface-chemistry model
cdr             (unfinished)
cdr             started to add a Planckian density model
cdr             as further option, similar to SAHA,
cdr             Boltzmann,... this one only for background
cdr             radiation field. (unfinished)
cdr Oct.  19:   In block 5 finally: Ti and V-flow multi-species (NPLS) input options completed and corrected
cdr             The treatment is now fully identical for both parameters.
cdr             Backward compatibility is achieved
cdr             by not reading indpro(2), but instead by infering it from the
cdr             number of input cards read. This card counting is done in find_param.
cdr sept. 19:   call alloc bckgrnd: indpro=6 and/or indpro =7 ?,
c               remove one redundant call
cdr  oct 18 :   unify reading of A&M data (reaction decks) from external file:
cdr             Formerly from block 4, block 5 ("density models")
cdr             and block 12 (line emissivities).
cdr             Now: unified interpreter of reaction card, subr. READ_REACLINES.f.
cdr             Reading only from the list in  block 4, IR=1,NREACI.
cdr             ALL reaction data potentially needed in blocks 5 and 12
cdr             must have been read in block 4, as well as their transfer to internal reaction
cdr             reaction data structures (calls to SLREAC.f) are already carried out.
cdr             Block 5 TDMPAR (density model) option: revised
cdr sept. 18:   iopt: ?? further optional input lines at the end of block 5?
cdr             for turning on/off input tally storage, and for gradients of
cdr             input tallies (optional)
cdr  sept.18:   XDR format options for fort.13 stream: removed.
cdr  apr. 18:   fully connected and tested: trchktm option, in block 11.
cdr  july 17 :  GR cleanup: wrmesh option split into writing and plotting
Cdr  april 17:  some cleanup (spelling, trim(character)) adopted from sols_iter version
cdr             added: logical NEXVS (default: F. Unclear meaning, so far...)
c               added: logical TRCRNF (traceback for random seeds for correlated sampling)
cdr  sept. 16:  extend options for extrapolations for A&M data beyond range
cdr             of tables or validity range fit expressions.

!               rename RMN and RMX to R1MN, R1MX, add R2MN, R2MX for range
!               of second variable in fit or data table
!               same with jfexmn,jfexmx (parameters to select extrapolation scheme)

!pb  June  16:  default for NPLSTI changed from 1 to NPLSI
cdr             indpro(2) and indpro(4): try to synchronize the meaning, to be done
!cd  jan   16:  reset census start time to time0, even for time0=0.
!cd  dec.  15:  jj-nlim, rather than jj-nlimi, for non.dev.std. surfaces
!    april 15:  esptcr, esptsr: sputtered particle energy flags introduced
!cd  29.10.14:  reading external file for block 4&5: allow comment lines at the beginning of file
!               (same in find_param)
!cd  22.09.14:  1D case, levgeo=2:  do not call grid(2)
!cd  22.03.14:  option 'include filname' instead of block 4 and 5 tested and verified
!               some minor changes at transition from end of block ***3 and re-entry to block ***6
!pb  01.01.14:  options AMPTS, multiplier for ntcpu, npts(istra): added (input block 7)
!dr  03.03.10:  READING A&M DATA in block 4, commented, sorted,.....
!pb  20.03.08:  allocate and nullify estiml(1) if no input block 10F
!pb             is available
!pb  01.08.07:  save NLSRON in case of coupled run with short cycle
!    20.06.07:  check NVOLPR against NSPEZV_DIM
!pb  22.05.07:  input of NPTSDEL added in block 7
!pb  20.04.07:  check of geometry flags added
!pb  20.04.07:  check for casename of NLFEM or NLTET
!pb  20.04.07:  option for reading triangle grid from Eirene input
!pb             file removed
!pb  22.03.07:  LEVGEO=6 --> LEVGEO=10
!pb  21.03.07:  logicals in block 1 cleaned up
!pb  02.03.07:  fourth secondary group in species specification cards
!pb             introduced
!pb  02.03.07:  remove ESCD2* arrays, use local variables instead
!pb             sum up ESCD1 and ESCD2 contributions
!pb  01.12.06:  bug fix: advance line in input for tetrahedra
!pb  09.10.06:  save NZADD for higher timesteps
!dr  20.04.06:  fort.10 added as density model in block 5.
!dr             Also other density model may now refere to test
!dr             particle tallies, for postprocessing and iteration
!pb  02.03.06:  NLRAY: switch on raytracing method for stratum
!pb  20.01.06:  line of sight for cell-based spectrum introduced
!pb  12.01.06:  flag for cell-based spectrum added in block 10F
C    24.09 05:  CALL TO IF2COP NOW LATER, AFTER ALL PLASMA DATA ARE SET
C               SO THAT IF2COP CAN BE USED WITHOUT HAVING TO USE IF1COP
C
cdr june-05:  spectrum input (10F) extended: SPC_SHIFT,.....
cdr                                SPCPLT_X,SPCPLT_Y,SPCPLT_SAME
cdr           see corresponding changes in CESTIM (ESTIML...)
cdr  28.4.04: nhsts(ispz) introduced, to select species
cdr           for trajectory plot
cdr           default: = 0: "plot trajectory for this species"
cdr           new    : =-1: "do not plot trajectory for this species"
cdr           see modifications in plt2d.f from 28.4.04
cpb sept-05:  specification of filenames and paths for external
cpb           databases in input block 1 added. (CFILE-cards)
cpb
cpb           example:
cpb           CFILE AMJUEL /home/boerner/Database/AMdata/amjuel.tex

      SUBROUTINE EIRENE_INPUT
C
C   READ INPUT DATA AND SET DEFAULT VALUES
c   IN CASE IITER.GT.1 OR ITIMV.GT.1 : SKIP READING NEW INPUT FROM STREAM IUNIN.
C                                      ONLY INPUT DATA PROCESSING (STATEMENT 4000 FF)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CREF
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_CSPEI
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_CPES, ONLY: NPRS
      USE EIRMOD_PHOTON
      USE EIRMOD_TIMEA, ONLY: EIRENE_TIMEA0
      USE EIRMOD_SECOND_OWN, ONLY: EIRENE_SECOND_OWN
      USE EIRMOD_PROFILES, ONLY: EIRENE_PROFR
      USE EIRMOD_JSON
      USE EIRMOD_INFCOP, ONLY: EIRENE_IF1COP, EIRENE_IF2COP
      USE EIRMOD_MPI
      USE EIRMOD_OPENFILE, ONLY : EIRENE_OPENFILE

      IMPLICIT NONE
C
      TYPE(VOLUMEP), POINTER :: VOLCUR
C
      REAL(DP) :: VOLTOT_TAL, RTEST

C  RUN TIME STATISTICS IN INITIALIZATION PHASE, WITHIN INPUT.F
      REAL(DP) :: tpb1, tpb2
cdr   REAL(DP) :: timea

      REAL(DP), ALLOCATABLE :: SAREA_SAVE(:)
      REAL(DP), ALLOCATABLE :: RDUMMY(:,:)
      INTEGER :: IADTYP(0:4)
      INTEGER :: IERROR, IUNIN_SAVE, JSTREAM, NSOPT, ILIMPS,
     .           I, J, I1, I2, I3, JL, IO, IUSR, IFLG,
     .           ITALI, IRAD, IS, ISS, INC,
     .           JPLS, IRE, JSPZ, JTRJ, ISTRAI, NLJ, IENTRY,
     .           IO13

#ifdef USE_MPI
      INTEGER :: IER
#endif
      LOGICAL :: NLSRON_SAVE(NSTRA)
      LOGICAL :: LTEST, EX
      CHARACTER(10) :: CDATE, CTIME
      CHARACTER(80) :: MPI_LINE
      CHARACTER(420) :: ZEILE
      CHARACTER(8) :: FILENAME
      CHARACTER(3) :: FILENUMBER
      EXTERNAL :: EIRENE_CHECK_GEOM_CONSIST, EIRENE_CORRECT_STATS_INPUT,
     .            EIRENE_CUT_ADS_CELL, EIRENE_GEOUSR, EIRENE_GRID,
     .            EIRENE_INIUSR, EIRENE_INTVOL, EIRENE_MULTIG,
     .            EIRENE_MULTIP, EIRENE_NUM_BACKGROUND_SPECTRA,
     .            EIRENE_PLASMA, EIRENE_PLASMA_DERIV,
     .            EIRENE_PLAUSR, EIRENE_PLMESH,
     .            EIRENE_PREP_PLOTTING, EIRENE_PREP_STRATA,
     .            EIRENE_PROUSR, EIRENE_READ_BLK14_FIXFORM,
     .            EIRENE_READ_BLK14_JSON,
     .            EIRENE_READ_JSON, EIRENE_READ_FIXFORM,
     .            EIRENE_READ_MPI_STRATEGY_FIXED,
     .            EIRENE_READ_MPI_STRATEGY_JSON, EIRENE_READTL,
     .            EIRENE_REPLACE_STORED, EIRENE_RGEOM, EIRENE_RPLAM,
     .            EIRENE_SET_CELL_DIAMETER,
     .            EIRENE_SET_DERIVED_INPUT_PARAMETERS,
     .            EIRENE_SET_SPECIFIC_ZONES,
     .            EIRENE_SETAMD, EIRENE_SETEQ, EIRENE_SETFIT,
     .            EIRENE_SETPRM_INTAL,
     .            EIRENE_SETTXT, EIRENE_SETTXT_INTAL,
     .            EIRENE_SETUP_ATOMIC_WEIGHTS,
     .            EIRENE_SETUP_CHORD_SPECTRA,
     .            EIRENE_SETUP_DEFAULT_REACTIONS,
     .            EIRENE_SETUP_INMP, EIRENE_SETUP_ISPEZ,
     .            EIRENE_SETUP_LOC_REF_MODELS,
     .            EIRENE_SETUP_SPEC_IND_DISTRIB,
     .            EIRENE_SETUP_SURFACE_SWITCHES,
     .            EIRENE_SETUP_TIME_SURFACE,
     .            EIRENE_SETUP_VIS_ADD_SURF_RANGES, EIRENE_VOLUME,
     .            EIRENE_WRGEOM, EIRENE_WRITE_JSON_AMDATA,
     .            EIRENE_WRITE_JSON_FILE, EIRENE_WRMESH, EIRENE_WRPLAM,
     .            EIRENE_XSECTPH,
     .            EIRENE_LEER, EIRENE_MASBOX, EIRENE_MASBR2,
     .            EIRENE_MASIR2, EIRENE_MASJ1, EIRENE_MASJ2,
     .            EIRENE_MASJ3, EIRENE_PAGE, EIRENE_EXIT_OWN
      EXTERNAL :: GREND
C
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST ITERATION
C  STEP IN THIS RUN. IITER IS THE ACTUAL ITERATION NUMBER
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST TIMESTEP
C  IN THIS RUN. ITIMV IS THE ACTUAL TIMESTEP NUMBER
C
C  INITIALISE SOME DATA AND SET DEFAULTS
C
      IERROR=0

!pb IUNIN set in COMPRT
      IUNIN_SAVE = IUNIN
      IUSROUT = 0
C
C  UNIT NUMBERS FOR INPUT FILE: MUST BE DIFFERENT FROM: 8,10,11,12
C  13,14,AND 15
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
C  UNIT NUMBERS COME FROM FIND_PARAM
      IF (IUNIN.eq.IUNOUT) THEN
        WRITE(iunout,*)
     .   'IUNIN cannot have the same value as IUNOUT!'
        WRITE(iunout,*) 'IUNIN = ',IUNIN
        WRITE(iunout,*) 'EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
      ELSE IF (IUNIN.EQ.8  .OR. IUNIN.EQ.10 .OR. IUNIN.EQ.11 .OR.
     1         IUNIN.EQ.12 .OR. IUNIN.EQ.13 .OR. IUNIN.EQ.14 .OR.
     2         IUNIN.EQ.15 .OR. IUNIN.EQ.0) THEN
        WRITE(iunout,*)
     .   'IUNIN MUST BE DIFFERENT FROM: 0,8,10,11,12,13,14,15 !'
        WRITE(iunout,*) 'IUNIN = ', IUNIN
        WRITE(iunout,*) 'EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
      ELSE IF (IUNOUT.EQ.8  .OR. IUNOUT.EQ.10 .OR. IUNOUT.EQ.11 .OR.
     1         IUNOUT.EQ.12 .OR. IUNOUT.EQ.13 .OR. IUNOUT.EQ.14 .OR.
     2         IUNOUT.EQ.15 .OR. IUNOUT.EQ.0) THEN
        WRITE(iunout,*)
     .   'IUNOUT MUST BE DIFFERENT FROM: 0,8,10,11,12,13,14,15 !'
        WRITE(iunout,*) 'IUNOUT = ', IUNOUT
        WRITE(iunout,*) 'EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
      END IF
C
      TPB1=EIRENE_SECOND_OWN()
      IF (IITER.GT.1) THEN
        CALL EIRENE_LEER(2)
        CALL EIRENE_MASJ1('IITER   ',IITER)
        CALL EIRENE_MASBOX
     .   ('NEXT ITERATION STARTS, SKIP READING INPUT FILE')
        GOTO 4000
      ENDIF

      IF (ITIMV.GT.1) THEN
        CALL EIRENE_LEER(2)
        CALL EIRENE_MASJ1('ITIMV   ',ITIMV)
        CALL EIRENE_MASBOX
     .   ('NEXT TIME CYCLE STARTS, SKIP READING INPUT FILE ')
        GOTO 4000
      ENDIF
      TPB1=EIRENE_SECOND_OWN()
C
      NAINI=0
      NCPVI=0
      NBGVI=0
      NPLS_FIX = 0    ! set to NFLA in coupled runs, in interfaces.

      MTSURF=0 ! cdr  ????
C
cdr SET MINIMAL REACTION MODELS (H/He plasmas), unless other reaction data are specified
cdr i.e. should be called only in case if any NRCA, NRCM, NRCI = 0 for any IATM, IMOL, IION
C
      CALL EIRENE_SETUP_DEFAULT_REACTIONS

C
C  SET DEFAULT SOURCE MODEL BLOCK 7
C
      NSTRAI=0
C
C  SET DEFAULT 'ADDITIONAL SURFACE' AND 'STANDARD SURFACE' DATA
C
      NBITS=BIT_SIZE(I)
      CALL EIRENE_SET_DEF_SURF_DATA

      ALLOCATE (SAREA_SAVE(NLIMPS))
      SAREA_SAVE = 666.

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

C  AS DEFAULT: SWITCH OFF MOMENTUM DENSITY TALLIES.
C  FOR ACTIVATING THOSE TALLIES THEY NEED TO BE EXPLICITLY
C  SWITCHED ON IN BLOCK 11

C  LV?DEN.. IS AN ALIAS FOR AN ENTRY IN ARRAY LMISTALV
C  THEREFORE .TRUE. MEANS: TALLY IS SWITCHED OFF
      LVXDENA  = .TRUE.
      LVXDENM  = .TRUE.
      LVXDENI  = .TRUE.
      LVXDENPH = .TRUE.
      LVYDENA  = .TRUE.
      LVYDENM  = .TRUE.
      LVYDENI  = .TRUE.
      LVYDENPH = .TRUE.
      LVZDENA  = .TRUE.
      LVZDENM  = .TRUE.
      LVZDENI  = .TRUE.
      LVZDENPH = .TRUE.
C
      CALL EIRENE_LEER(2)

      CALL EIRENE_ALLOC_CTEXT(3)
C
C  READ TEXT DESCRIBING THE RUN, 100--199
C
      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) ' CPU time before reading ',tpb2-tpb1
      tpb1 = tpb2

      CALL DATE_AND_TIME(CDATE,CTIME)
      READ(CDATE(1:4),*) I1
      READ(CDATE(5:6),*) I2
      READ(CDATE(7:8),*) I3
      WRITE (iunout,'(1X,A6,1X,2(I2,1X),I4)') 'DATE: ',I3,I2,I1
      READ(CTIME(1:2),*) I1
      READ(CTIME(3:4),*) I2
      READ(CTIME(5:6),*) I3
      WRITE (iunout,'(1X,A6,1X,3(I2,1X))') 'TIME: ',I1,I2,I3
      CALL EIRENE_LEER(2)

      IF (IUNIN.NE.5) REWIND IUNIN !VK

      READ (IUNIN,'(A80)') ZEILE

      IF (IUNIN.NE.5) REWIND IUNIN !VK
      LRDJSON = .FALSE.
      IF (ZEILE(1:1) == '{') THEN
        LRDJSON = .TRUE.
      ELSEIF (ZEILE(1:1) /= '*') THEN
        WRITE (IUNOUT,*) ' EIRENE INPUT FORT.1 HAS WRONG FORMAT'
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IF (LRDJSON) THEN

        close (unit=1)
        CALL EIRENE_READ_JSON(NLIMPS, SAREA_SAVE, ierror)

      ELSE

        CALL EIRENE_READ_FIXFORM (NLIMPS, SAREA_SAVE, IERROR)

      END IF
C
      CALL EIRENE_SETUP_TIME_SURFACE (IERROR)
C
      IF (IERROR.GT.0) THEN
        WRITE (iunout,*) IERROR,' INPUT OR PARAMETER ERRORS DETECTED'
        WRITE (iunout,*)
     .   ' SEE THE ERROR MESSAGES LISTED ABOVE AND CORRECT'
        WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      CALL EIRENE_PAGE
cdr done with reading formatted input, blocks 1 to 13.
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
cdr This problem-specific part is done further below. First:
cdr Set some derived input data from blocks 1 to 13.

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) ' CPU time after reading ',tpb2-tpb1
      tpb1 = tpb2
C
C
C   MODIFICATION OF INPUT DUE TO EITHER INCONSISTENCIES OR DUE
C   TO COUPLED NEUTRAL-PLASMA (OR NEUTRAL-NEUTRAL) CALCULATIONS
C   SOME FURTHER CONSTANTS ARE SET.      STATEM. NO. 2000 --> 3999

      CALL EIRENE_CHECK_GEOM_CONSIST
C
C  SOURCE PARAMETERS AND (REFLECTING) BOUNDARY CONDITIONS,
C  ON ADDITIONAL AND NON-DEFAULT STANDARD SURFACES
C
      CALL EIRENE_SETUP_LOC_REF_MODELS

      CALL EIRENE_SETUP_SURFACE_SWITCHES
C
C  SET NON-DEFAULT STANDARD SURFACE IDENTIFIERS INMP...
C
      CALL EIRENE_SETUP_INMP
C
      CALL EIRENE_PREP_STRATA
C
C
C
C  SPECIES INDEX DISTRIBUTION OF PRIMARY SOURCE PARTICLES
C  OR FOR THERMAL PARTICLE REFLECTION MODEL

      CALL EIRENE_SETUP_SPEC_IND_DISTRIB

      CALL EIRENE_SETUP_ATOMIC_WEIGHTS
C
C
C  SET SOME ARRAYS TO SPEED UP COMPUTATIONS
C
      CALL EIRENE_SETUP_ISPEZ

!  this should go into setamd(0)
      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(1)
      CALL EIRENE_SETAMD(0)
      CALL EIRENE_ALLOC_CTEXT(2)

      CALL EIRENE_SETTXT_INTAL
      CALL EIRENE_SETPRM_INTAL
      CALL EIRENE_SETTXT
C
C
C  ADDITIONAL INPUT FOR THIS RUN COMES FROM EITHER
C  ANOTHER CODE (DATA FILE) OR FROM AN EARLIER RUN OF EIRENE
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM)
     > write (iunout,*) ' CPU time before block 14 ',tpb2-tpb1
      tpb1 = tpb2

! CALL TO ALLOC_BCKGRND MOVED HERE TO ALLOW SPECIFICATION OF VOL
! IN IF0COP

      IF ((ANY(INDPRO(1:12) == 6)).OR.(ANY(INDPRO(1:12) == 7)))
     .    CALL EIRENE_ALLOC_BCKGRND
cdr  Warning: meanwhile the target bckgrnd is used also for other purposes,
cdr  not just its original purpose of transfer of background fields
cdr  from external codes. E.g.: bgk iterations ? spectra?

      IF (LRDJSON) THEN

        CALL EIRENE_READ_BLK14_JSON(IERROR)

        CALL EIRENE_READ_MPI_STRATEGY_JSON

C  CHECK FOR USER SPECIFIC INPUT
        IUSR = -9999
        CALL EIRENE_OPENFILE(IUSR,FILE='user_data.input',STATUS='OLD',
     .       FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IO)
        IUSROUT = IUSR
        IUNIN_SAVE = IUNIN
        IF (IO == 0) THEN
          IUNIN = IUSROUT
          CALL EIRENE_LEER(1)
          WRITE (IUNOUT,*) 'USER SPECIFIC INPUT READ FROM ',
     .         'user_data.input'
        ELSE
          CALL EIRENE_LEER(1)
          WRITE (IUNOUT,*) 'NO FILE FOR USR SPECIFIC INPUT FOUND'
        END IF

      ELSE

        CALL EIRENE_READ_BLK14_FIXFORM(IERROR)

C  COPY USER SPECIFIC DATA TO FILE user_data.input
        JL = 0
        IO = 0
        IUNIN_SAVE = IUNIN
        MPI_LINE=REPEAT(' ',80)
        DO WHILE (IO == 0)
          READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
          IF (INDEX(ZEILE,'INFORMATION_FOR_MPI') /= 0) THEN
            MPI_LINE = ZEILE(1:80)
            EXIT
          ELSE
            IF (IO == 0) THEN
              JL = JL + 1
              IF (JL == 1) THEN
                IUSR = -9999
                IF (IUSROUT == 0) THEN
                  CALL EIRENE_OPENFILE(IUSR,FILE='user_data.input',
     .               FORM='FORMATTED',ACCESS='SEQUENTIAL')
                ELSE ! trick to avoid a compiler error
                  CLOSE(IUSROUT)
                  CALL EIRENE_OPENFILE(IUSR,FILE='user_data.input',
     .               FORM='FORMATTED',ACCESS='SEQUENTIAL',
     .               STATUS='OLD',POSITION='APPEND')
                END IF
                IUSROUT = IUSR
              END IF
              WRITE (IUSROUT,'(A)') TRIM(ZEILE)
            END IF
          END IF
        END DO

        IF (VERIFY(MPI_LINE,' ') > 0) THEN
          CALL EIRENE_READ_MPI_STRATEGY_fixed(MPI_LINE)
        END IF

        IF (JL > 0) THEN
          REWIND IUSROUT
          IUNIN = IUSROUT
        ELSE
          IUSROUT = 0
        END IF

      END IF

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) ' CPU time after if0cop ',tpb2-tpb1
      tpb1 = tpb2
C
C  INPUT BLOCK 14 DONE
C
      CALL EIRENE_PREP_PLOTTING

      CALL EIRENE_CORRECT_STATS_INPUT

C
C  NO MODIFICATION OF INPUT VARIABLES BEYOND THIS POINT
C  WITHOUT WARNING
C  EXCEPT IN SUBROUTINE MODUSR FOR THE NEXT ITERATION STEP
C         IN SUBROUTINE TMSUSR FOR THE NEXT TIME STEP
C
C
      CALL EIRENE_INIUSR
C
C  SET DERIVED INPUT PARAMETERS, GRIDS AND PROFILES
C
      CALL EIRENE_SET_DERIVED_INPUT_PARAMETERS(IERROR)

      CALL EIRENE_SET_PARMMOD(3)
      CALL EIRENE_ALLOC_CGEOM(2)
      CALL EIRENE_ALLOC_COMUSR(3)

      SAREA(1:NLIMPS) = SAREA_SAVE(1:NLIMPS)
      DEALLOCATE (SAREA_SAVE)

C
C  WRITE JSON FILE
      IF (LRDJSON) THEN
        CALL EIRENE_WRITE_JSON_AMData('eirene_AMData.json')
        CALL EIRENE_WRITE_JSON_FILE('eirene_input_json.out')
      ELSE
        CALL EIRENE_WRITE_JSON_FILE('eirene.input.json')
      END IF

!  REMOVE LIST OF REFLECTION MODELS
!  MOVED HERE AS TO BE AVAILABLE FOR WRITING ON JSON-FILE

      CALL EIRENE_DEALLOC_REFLIST

      CALL EIRENE_PAGE

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) ' CPU time before grid ',tpb2-tpb1
      tpb1 = tpb2

C
      IF (NFILEM.LE.1) THEN
C
C  SET GRIDS AND VOLUMES OF THE CELLS FOR THE STANDARD MESHES
C
C  SET RADIAL OR X GRID
        IF (NLRAD) CALL EIRENE_GRID (1)
C  SET POLOIDAL OR Y GRID
! NO NEED TO SET UP POLYGON GRID IN CASE OF 1D RUN, LEVGEO=2
        IF (NLPOL.OR.(LEVGEO == 3)) CALL EIRENE_GRID (2)
C  SET TOROIDAL OR Z GRID
        IF (NLTOR) CALL EIRENE_GRID (3)
C
        IF (INDPRO(12).LT.4) THEN
C  INITIALISE SUBROUTINE VOLUME FOR LATER CALLS
C  TO BE WRITTEN: CALL VOLUME(0)
C  SET VOLUMES, 1ST DIMENSION (R/X-GRID)
          IF (NLRAD) CALL EIRENE_VOLUME(1)
C  SET VOLUMES, 2ND DIMENSION (THETA/Y-GRID)
          IF (NLPOL) CALL EIRENE_VOLUME(2)
C  SET VOLUMES, 3RD DIMENSION (PHI/Z-GRID)
          IF (NLTOR) CALL EIRENE_VOLUME(3)
C  SET VOLUMES, IN ADDITIONAL CELL REGION
          IF (NLADD) CALL EIRENE_VOLUME(4)
        ELSEIF (INDPRO(12).EQ.4) THEN
C  READ VOLUMES FROM STREAM JSTREAM=VL0
          JSTREAM=NINT(VL0)
          ITALI=NTALO
          CALL EIRENE_READTL(TXTPLS(1,ITALI),TXTPSP(1,ITALI),
     .                TXTPUN(1,ITALI),
     .                VOL,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                3,JSTREAM)
C  TAKE VOLUMES FROM USER-SUPPLIED ROUTINE
        ELSEIF (INDPRO(12).EQ.5) THEN
          CALL EIRENE_PROUSR (VOL,5+5*NPLS,0._DP,0._DP,0._DP,0._DP,
     .                              0._DP,0._DP,0._DP,NSBOX)
C  TAKE VOLUMES FROM EXTERNAL FILE
        ELSEIF (INDPRO(12).EQ.6) THEN
          CALL EIRENE_PROFR (VOL,5+1*NPLS+NPLSTI+3*NPLSV,1,1,NSURF)
          IF (NLADD) CALL EIRENE_VOLUME(4)
        ELSEIF (INDPRO(12).EQ.7) THEN
          CALL EIRENE_PROFR (VOL,5+1*NPLS+NPLSTI+3*NPLSV,1,1,NSBOX)
        ENDIF
C
C  MULTIPLY GEOMETRICAL DATA, IF NBMLT.GT.1
C
        IF (NBMLT.GT.1) CALL EIRENE_MULTIG
C
C  SET CELL DIAMETER
C
        CALL EIRENE_SET_CELL_DIAMETER
C
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES
C
        DO WHILE (ASSOCIATED(VOLLIST))
          VOL(VOLLIST%II) = VOLLIST%VOL
          VOLCUR => VOLLIST
          VOLLIST => VOLLIST%NEXT
          DEALLOCATE(VOLCUR)
        ENDDO
        IF (NPHOTI > 0) CALL EIRENE_PH_INIT(2)
C
C   MODIFY SOME GEOMETRICAL DATA, USER-SUPPLIED ROUTINE
C
        CALL EIRENE_GEOUSR

        IF (LEVGEO == 4) CALL EIRENE_CUT_ADS_CELL
C
C  SET SOME DATA FOR ADDITIONAL SURFACES: INITIALISE SUBR. TIMEA
C
cdr     timea = EIRENE_SECOND_OWN()
        CALL EIRENE_TIMEA0
cdr     WRITE(iunout,*)'CPU time for timea0 ',EIRENE_SECOND_OWN()-timea
cdr     WRITE(iunout,*)
C
C   MODIFY THE BOUNDARIES OF SOME SURFACES TO AVOID ROUND-OFF
C   ERRORS
C
        CALL EIRENE_SETFIT (TRCSUR)
C
C   SET THE COEFFICIENTS OF SOME SURFACES IDENTICAL TO THOSE
C   OF SOME OTHER, TO AVOID ROUND-OFF ERRORS
C
        CALL EIRENE_SETEQ
C
C   WRITE A LIST OF CLOSED POLYGONIAL LINES, DETERMINED
C   FROM EIRENE ADDITIONAL AND NON-DEFAULT STANDARD SURFACES
C   (WITH THEIR ORIENTATION) ONTO STREAM 78+IFOFF
C   FOR FURTHER USE IN TRIANGULARISATION CODES, WHICH MAY THEN
C   PRODUCE GRIDS OF UNSTRUCTURED TRIANGLES INSIDE THESE CLOSED
C   POLYGONS (EXCLUDING THOSE AREAS WHICH ARE DESCRIBED BY
C   CLOSED POLYGONS WITH NEGATIVE ORIENTATION)
C
        IF (NLWRMSH) THEN
          CALL EIRENE_WRMESH
          CALL EIRENE_PLMESH
        ENDIF
C
C
        ALLOCATE(RDUMMY(1,NSBOX))
        RDUMMY(1,1:NSBOX) = VOL(1:NSBOX)
        CALL EIRENE_INTVOL (RDUMMY,1,1,NSBOX,VOLTOT,
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        VOL(1:NSBOX) = RDUMMY(1,1:NSBOX)
        DEALLOCATE(RDUMMY)
        WRITE (iunout,*) 'TOTAL VOLUME, SUM VOL(:)  ',VOLTOT
        CALL EIRENE_LEER(1)
C
C  SET 'VISIBLE ADDITIONAL SURFACES' RANGES nlimii(j),nlimie(j), for each grid cell j
C  FROM INFORMATION ON IGJUM3
        CALL EIRENE_SETUP_VIS_ADD_SURF_RANGES

C
C  ALL GEOMETRICAL DATA (GRIDS, VOLUMES, SWITCHES) ARE DEFINED NOW
C
C   SAVE GEOMETRICAL DATA ON FILE FT12
C
        DO 8006 IRAD=1,NSBOX
          VOLG(IRAD)=VOL(IRAD)
 8006   CONTINUE
C
        DO ILIMPS=1,NLMPGS
          AREAG(ILIMPS)=SAREA(ILIMPS)
        ENDDO
C
        IF (NFILEM.EQ.1) CALL EIRENE_WRGEOM(TRCFLE)

!VK+AK in order to use EIRENE only for mesh generation
        if(nltrimesh) then
          write(iunout,*) 'NLTRIMESH=.true.'
          IF (NPRS > 1) THEN
            CALL EIRENE_EXIT_OWN(1)
          ELSE
            CALL GREND
#ifdef USE_MPI
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, IER)
            CALL MPI_FINALIZE(IER)
#endif
            STOP "EIRENE TRIANGULATION PREPARATION RUN COMPLETE"
          END IF
        END IF
C
      ELSEIF (NFILEM.EQ.2) THEN
C
C   RESTORE GEOMETRICAL DATA FROM FILE FT12
C
        CALL EIRENE_RGEOM(TRCFLE)
C
        DO 8010 IRAD=1,NSBOX
          VOL(IRAD)=VOLG(IRAD)
 8010   CONTINUE
C
        DO ILIMPS=1,NLMPGS
          SAREA(ILIMPS)=AREAG(ILIMPS)
        ENDDO
C
C
      ELSE
C
        WRITE(IUNOUT,*) 'UNRECOGNIZED NFILEM OPTION. NFILEM = ', NFILEM
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      DO 8011 IS=1,NLIMPS
        IF (ILACLL(IS).NE.0.AND..NOT.NLADD) THEN
          WRITE (iunout,*) 'ADDITIONAL CELL SWITCHES DEFINED, BUT NO'
          WRITE (iunout,*) 'ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCK 2E'
          CALL EIRENE_EXIT_OWN(1)
        ELSEIF (ILBLCK(IS).NE.0.AND..NOT.(NLMLT.OR.NLADD)) THEN
          WRITE (iunout,*) 'BLOCK SWITCHES DEFINED, BUT NEITHER BLOCKS'
          WRITE (iunout,*) 'NOR ADDITIONAL CELLS DEFINED'
          ISS=IS
          IF (ISS.GT.NLIM) ISS=-(ISS-NLIM)
          WRITE (iunout,*) 'SURFACE NO. IS= ',ISS
          WRITE (iunout,*) 'CHECK INPUT BLOCKS 2D AND 2E'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
 8011 CONTINUE

C
      IF(TRCREF) THEN
       WRITE(IUNOUT,*) "PROPERTIES OF ALL SURFACES"
       WRITE(IUNOUT,*)
       DO I=1,NLIM+NSTS
        IF(I.GT.NLIMI) THEN
         NLJ=NLIMI-I
        ELSE
         NLJ=I
        END IF
        WRITE(IUNOUT,*) "SURFACE ", NLJ
        WRITE(IUNOUT,*) " ILIIN, ILSPT, ISPUT ",
     w                    ILIIN(I),ILSPT(I),ISPUT(:,I)
        WRITE(IUNOUT,'(A12,10I6/(12x,10I6))') " ISRS ", ISRS(:,I)
        WRITE(IUNOUT,'(A12,10I6/(12x,10I6))') " ISRC ", ISRC(:,I)
        WRITE(IUNOUT,'(A12,10I6/(12x,10I6))')
     w                                    " LCHSPNWL ", LCHSPNWL(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                     " TRANSP1 ", TRANSP(:,1,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                     " TRANSP2 ", TRANSP(:,2,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " RECYCF ", RECYCF(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " RECYCT ", RECYCT(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " RECPRM ", RECPRM(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                       " EXPPL ", EXPPL(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                       " EXPEL ", EXPEL(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                       " EXPIL ", EXPIL(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " RECYCS ", RECYCS(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " RECYCC ", RECYCC(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " SPTPRM ", SPTPRM(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " ESPUTS ", ESPUTS(:,I)
        WRITE(IUNOUT,'(A12,1P,5E12.4/(12x,5E12.4))')
     w                                      " ESPUTC ", ESPUTC(:,I)
        CALL EIRENE_LEER(1)
       END DO
      END IF
CVK END

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM)
     . write (iunout,*) ' CPU time before plasma definition ',tpb2-tpb1
      tpb1 = tpb2
C
 4000 CONTINUE   !  at this point:
cdr                 continue INPUT.f in case iiter > 1 or itimv>1,
cdr                 after having skipped all the rest above.
C
      if ((iiter.gt.1) .or. (itimv.gt.1)) then
        write (iunout,*) 'iterative or time dep.mode. input 4000'
        write (iunout,*) 'iiter,itimv,nfilel,nfilej,nlshrt13,iprnl'
        write (iunout,*)  iiter,itimv,nfilel,nfilej,nlshrt13,iprnl
      endif

!  NOTHING IS DONE IF ARRAYS FOR BACKGROUND FIELDS ARE ALREADY ALLOCATED
      IF (ANY(INDPRO(1:12) == 6) .or. ANY(INDPRO(1:12) == 7))
     .   CALL EIRENE_ALLOC_BCKGRND

      IF ((NMODE.NE.0.AND.IITER.LE.1) .OR.
     .    (ABS(NMODE).EQ.2)) THEN  !dr: for EMC3 coupling?
                                   !    in case IITER .GT. 1
C  READ PLASMA BACKGROUND
c  EITHER: FROM EXTERNAL DATABASE (FT31) (NOT NLPLAS)
C  OR    : FROM COMMON BRAEIR (NLPLAS)
        IENTRY=0
        CALL EIRENE_IF1COP(IENTRY)
      ENDIF

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) ' CPU time for if1cop ',tpb2-tpb1
      tpb1 = tpb2
C
      IF (NSTEP > 0) CALL EIRENE_ALLOC_CSTEP
C
cdr  VOL are cell volumes on the fine grid.
cdr  Coarser grid FOR SCORING may have been set, find cell volumes
cdr  VOLTAL on coarser grid
      if (iiter.le.1 .and. itimv.le.1) then
        VOLTAL = EPS60
        DO I=1,NSBOX
          INC = NCLTAL(I)
          IF (INC > 0) VOLTAL(INC) = VOLTAL(INC) + VOL(I)
        END DO
        ALLOCATE(RDUMMY(1,NSBOX))
        RDUMMY(1,1:NSBOX) = VOLTAL(1:NSBOX)
        CALL EIRENE_INTVOL (RDUMMY,1,1,NSBOX_TAL,VOLTOT_TAL,
     .                      NR1TAL,NP2TAL,NT3TAL,NBMLT)
        VOLTAL(1:NSBOX) = RDUMMY(1,1:NSBOX)
        DEALLOCATE(RDUMMY)
        WRITE (iunout,*) 'TOTAL VOLUME, SUM VOLTAL(:) ',VOLTOT_TAL
        call eirene_leer(1)
      ENDIF

cdr  Jan 2020:
Cdr  Facilitate "cold" initialisation of an iterative run, when no fort.13
c    from an earlier cycle is available.
c    This did work already for NLSHRT13 = T
c   (then: zero virt. densities are set on target PLASMA_BCKGRND).
c    Now try also for long fort.13: NLSHRT13 = F option.
cdr  better:  use initial parameters from block 5, if they are available,
cdr  to initialize the virtual field particles IPLS=NPLS_FIX+1,NPLS
c
C   Test if data stream FORT.13 is prepared.
cdr If not (IO13 = 0), then write (rather than read) fort.13
cdr despite NFILEL .ge. 2

      io13=0
cdr  Do we need to fiddle with fort.13?
      if (nfilel.ne.0) then
cdr  Yes.
        if (13+ifoff.ge.100) then
          WRITE(FILENUMBER,'(I3)') 13+ifoff
        else
          WRITE(FILENUMBER,'(I2)') 13+ifoff
        end if
        FILENAME = FORT_LC//trim(FILENUMBER)
        INQUIRE(file=trim(FILENAME),exist=EX)
        IF (EX) THEN
          OPEN (UNIT=13+ifoff,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
          REWIND 13+ifoff
          if (.not.nlshrt13) then
            READ (13+ifoff,IOSTAT=IO13) LTEST
          else
            READ (13+ifoff,IOSTAT=IO13) RTEST
          endif
        END IF
        if (io13.eq.0 .and. ex) then
          write (iunout,*) 'INPUT.f: '//fort_lc//'13 found'
        else
          write (iunout,*) 'INPUT.f: '//fort_lc//'13 not found'
        end if
        if (ex) CLOSE (UNIT=13+ifoff)
      endif
cdr

cdr set plasma, plasma_deriv, and create fort.13
      IF (NFILEL.LE.1 .OR. NLSHRT13 .OR. IO13.NE.0 .or. .not.ex) THEN
C
C  SET PLASMA PARAMETERS AND PRIMARY SOURCE PARAMETERS
C  FROM INPUT BLOCKS 5 AND 7, RESP.
C
        CALL EIRENE_PLASMA
C
C  MULTIPLY PLASMA PARAMETERS, IF NBLCKS.GT.1
C
        IF (NBLCKS.GT.1) CALL EIRENE_MULTIP
C
C   INCLUDE INFORMATION PROVIDED BY INPUT BLOCK 8: ADDITIONAL
C   DATA FOR SPECIFIC ZONES:  ONLY AT START OF ITERATIVE LOOP
C   OR AT FIRST INTERNAL TIME STEP
C
cdr do this in first iteration step or first time step
cdr of each iteration cycle or time cycle
        IF ((NZADD.GT.0) .AND. (IITER == 1)
     .                   .AND. (ITIMV == 1))
     .    CALL EIRENE_SET_SPECIFIC_ZONES
C
C  MODIFY SOME PLASMA DATA, USER-SUPPLIED ROUTINE
C
        CALL EIRENE_PLAUSR

        TPB2=EIRENE_SECOND_OWN()
        IF (TRCTIM) write (iunout,*) ' CPU time for plausr ',tpb2-tpb1
        tpb1 = tpb2

        IF (NFILEL.EQ.3.OR.NFILEL.EQ.4) then
C  THIS POINT IS ONLY REACHED WITH: (NLSHRT13=.T. .OR. IO13 /= 0)
c  (EITHER SHORT VERSION OF FORT.13 ONLY, OR NO FORT.13 PREPARED AT ALL).

c  in case nlshrt13:
C  READ ONLY PLASMA DATA OF THOSE BACKGROUND SPECIES NPLS_FIX+1:NPLS
C  WHICH ARE NOT CONTAINED IN EXTERNAL PLASMA CODE,
C  I.E. ONLY THOSE WHICH ARE NEEDED FOR INTERNAL EIRENE CYCLING
C      (NONLINEARITIES)
c  in case IO13 /=0:
!dr For nlshrt13=t: No problem if fort.13 is not ready.
cdr Then just do nothing for virt. background.
!dr For nlshrt13=f: return iflg=-1
          IFLG=0
          CALL EIRENE_RPLAM(TRCFLE,IFLG,'INPUT 1')

        ENDIF
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
        CALL EIRENE_PLASMA_DERIV(0)

        TPB2=EIRENE_SECOND_OWN()
        IF (TRCTIM)
     >   write (iunout,*) ' CPU time for plasma_deriv ',tpb2-tpb1
        tpb1 = tpb2
C
C  SET ATOMIC DATA TABLES: only in first internal iteration and in first time step
C
        if ((iiter.le.1) .and. (itimv.le.1)) CALL EIRENE_SETAMD(1)

cdr iteration on virtual species: reset some A&M data,
cdr from previous iteration, set in calls to modbgk.f, modphot.f  etc.

        if (nfla_virt.gt.0 .or. nrea_virt.gt.0)
     .      call eirene_replace_stored

        TPB2=EIRENE_SECOND_OWN()
        IF (TRCTIM) write (iunout,*) ' CPU time for setamd ',tpb2-tpb1
        tpb1 = tpb2
C
C  WRITE PLASMA DATA, ATOMIC DATA, SOURCE DATA FROM FT13
C
        IF (NFILEL.EQ.1) CALL EIRENE_WRPLAM(TRCFLE,'INPUT 1')

        TPB2=EIRENE_SECOND_OWN()
        IF (TRCTIM) write (iunout,*) ' CPU time for wrplam ',tpb2-tpb1
        tpb1 = tpb2

C
      ELSEIF (NFILEL.GE.2.AND.NFILEL.LE.4) THEN   ! =2,3,4,
                                  ! and NLSHRT13=F, necessarily
C
C  READ PLASMA DATA, ATOMIC DATA, SOURCE DATA FROM FT13
C
        NLSRON_SAVE = NLSRON

        CALL EIRENE_ALLOC_CORNERS

        IFLG = 0
        IF ((ITIMV > 1) .AND. NLMOVIE) IFLG = 1
cdr  iflg rather than nfilel  ??
        IF ((NFILEL == 2) .OR. (NFILEL == 3)) THEN
          CALL EIRENE_RPLAM(TRCFLE,IFLG,'INPUT 2')

cdr  from now on: iflg=4,
cdr  skip reading A&M data and primary source data from fort.13
        ELSEIF (NFILEL == 4) THEN
          IFLG=4
          CALL EIRENE_RPLAM(TRCFLE,IFLG,'INPUT 3')
        END IF
        if (iflg.lt.0) then
cdr  we have returned from rplam_long, but could not read fort.13.
cdr  So we try to create a stream fort.13 now (wrplam_long).
          CALL EIRENE_WRPLAM(TRCFLE,'INPUT 4')
        endif

cdr Is xsectph outside setamd?
        CALL EIRENE_XSECTPH

        IF (IITER > 1) NLSRON = NLSRON_SAVE
C
      ENDIF

      TPB2=EIRENE_SECOND_OWN()
      IF (TRCTIM)
     > write (iunout,*) ' CPU time after plasma definition ',tpb2-tpb1
      tpb1 = tpb2

C
C  SET UP TABLE OF CONTRIBUTIONS OF MONTE CARLO PARTICLES TO BACKGROUND SPECIES
CDR: in case of "density models": multiple test particle components for one new background tally,
cdr  e.g. COLRAD models

cdr  independent:
cdr  iadtyp(ityp) = increment in 1D species arrays
cdr  potentially useful in many places
C
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO JPLS = 1, NPLSI
        IF ((LEN_TRIM(CDENMODEL(JPLS)) > 0) .AND.
     .      (INDEX(CDENMODEL(JPLS),'CONSTANT') == 0)) THEN
          DO IRE = 1, TDMPAR(JPLS)%TDM%NRE
            ITYP = TDMPAR(JPLS)%TDM%ITP(IRE)
            ISPZ = IADTYP(ITYP) + TDMPAR(JPLS)%TDM%ISP(IRE)
            ISPZ_BACK(ISPZ,JPLS) = 1
          END DO
        END IF
      END DO

cdr 2020  remove irrelevant/confusing printout
      if (ANY(ISPZ_BACK(:,:) > 0)) then
        CALL EIRENE_LEER(2)
        WRITE (IUNOUT,*) 'LIST OF CONTRIBUTIONS TO BACKGROUND SPECIES'
        DO JSPZ = 1, NSPZ
          DO JPLS = 1, NPLSI
            IF (ISPZ_BACK(JSPZ,JPLS) > 0)
     .        WRITE (IUNOUT,*) TEXTS(JSPZ), ' CONTRIBUTES TO ',
     .                         TEXTS(NSPAMI+JPLS)
          END DO
        END DO
      endif
C
C  AT THIS POINT THE BACKGROUND MEDIUM DATA ARE ALL SET.
C
C  COMPUTE PRIMARY SOURCE DATA (OVERRULE SOME OF INPUT BLOCK 7)

      IF (NMODE.NE.0.AND.((IITER.LE.1) .OR.
     .                    (IITER > NITER))) THEN
        DO ISTRAI=1,NSTRAI
          IF (INDSRC(ISTRAI).GE.0) CALL EIRENE_IF2COP(ISTRAI)
        ENDDO
      ENDIF

      IF (NTIME.GT.0) THEN
cdr time-stepping: reset census flux from previous call to tmstep.f
cdr Here now NSTRAI is the census stratum
        FLUX(NSTRAI)=FLXCEN

        IF (NPTST.EQ.0) THEN
C  SET NUMBER OF PARTICLES FOR RELAUNCH FROM CENSUS EQUAL TO THE NUMBER OF PREVIOUS SCORES ON CENSUS
C  (BUT STILL: SAMPLING WITH REPLACEMENT, BOOTSTRAPPING)
C  N.B.: WE HAVE ALREADY MADE SURE ABOVE, THAT IN CASE NLMOVIE: NPTST = -1
          NPTS(NSTRAI)=NPRNLI
          NMINPTS(NSTRAI)=NPRNLI
        ELSEIF (NPTST.GT.0) THEN
          NPTS(NSTRAI)=NPTST
        ELSEIF (NPTST.LT.0) THEN
C  ONE BY ONE RELAUNCH FROM OLD CENSUS
C  OLD CENSUS CONTAINS IPRNL ENTRIES.
          NPTS(NSTRAI)=IPRNL
          NMINPTS(NSTRAI)=IPRNL
        ENDIF
        CALL EIRENE_MASJ2('NSTRA_CEN, NPTS ',NSTRAI,NPTS(NSTRAI))
C
        IF (NPTS(NSTRAI).GT.0.AND.FLUX(NSTRAI).GT.0) THEN
CDR WE DO HAVE A CENSUS ARRAY FROM A PREVIOUS TIME STEP
          NSRFSI(NSTRAI)=1
          SORWGT(1,NSTRAI)=1.D0
          NLSRON(NSTRAI)=.TRUE.
          IF (NTMSTP.LT.0) THEN
CDR  NTMSTP LT.0 (I.E.: REQUEST TO PRODUCE A STATIONARY CENSUS)
CDR  IS ONLY A RELEVANT OPTION IF THERE IS AN EMPTY CENSUS FROM THE PAST
            NTMSTP=1
            WRITE (IUNOUT,*) 'NTMSTP RESET TO 1: FINITE TIME STEP'
            WRITE (IUNOUT,*) 'BECAUSE A NON-EMPTY CENSUS IS AVAILABLE'
          ENDIF
        ELSE
CDR EMPTY CENSUS. TURN OFF TIME STRATUM
          NLSRON(NSTRAI)=.FALSE.
        ENDIF
        NLSRON_SAVE=NLSRON

        call eirene_leer(1)
        write (iunout,*) 'ntime,nstrai,flux(nstrai),sorwgt(1,nstrai)'
        write (iunout,*)  ntime,nstrai,flux(nstrai),sorwgt(1,nstrai)
        CALL EIRENE_MASJ2('ISTRA_CEN, NPTS=',NSTRAI,NPTS(NSTRAI))
        IF (NLSRON(NSTRAI) .AND. NPTST.GE.0) THEN
          WRITE (IUNOUT,*) 'BOOTSTRAPPING FROM CENSUS'
        ELSEIF (NLSRON(NSTRAI) .AND. NPTST.LT.0) THEN
          WRITE (IUNOUT,*) 'ONE BY ONE RE-LAUNCH FROM CENSUS'
        ELSE
          WRITE (IUNOUT,*) 'NO CENSUS, NO TIME STRATUM'
        ENDIF
        call eirene_leer(1)


cdr   ELSEIF (NTIME.LT.0) THEN
cdr M.Rack option: read census (fort.15) and launch one by one, no time horizon.
cdr All done already above.
      ENDIF
C
C
      IF (TRCSUR) THEN
C
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'COEFFICIENTS FOR ADDITIONAL SURFACES'
        WRITE (iunout,*)
     .    'THIS IS AFTER IF0COP, GEOUSR, SETEQ AND SETFIT ARE CALLED'
        DO 7701 J=1,NLIMI
          CALL EIRENE_LEER(2)
          WRITE (iunout,*) TXTSFL(J)
          CALL EIRENE_LEER(1)
          IF (IGJUM0(J) == 1) THEN
            WRITE (iunout,*) 'THIS SURFACE IS NOT DEFINED'
          ELSE
            WRITE (iunout,*) 'A0       ',A0LM(J)
            WRITE (iunout,*) 'A1,A2,A3 ',A1LM(J),A2LM(J),A3LM(J)
            WRITE (iunout,*) 'A4,A5,A6 ',A4LM(J),A5LM(J),A6LM(J)
            WRITE (iunout,*) 'A7,A8,A9 ',A7LM(J),A8LM(J),A9LM(J)
            WRITE (iunout,*) 'JUMLIM ',JUMLIM(J)
            WRITE (iunout,*) 'ISWICH(1),ISWICH(2),ISWICH(3),ISWICH(4),',
     .                  'ISWICH(5),ISWICH(6)'
            WRITE (iunout,*)  ISWICH(1,J),ISWICH(2,J),ISWICH(3,J),
     .                   ISWICH(4,J),ISWICH(5,J),ISWICH(6,J)
            IF (RLB(J).GT.0) THEN
                WRITE (iunout,*) 'XLIMS1,YLIMS1,ZLIMS1 ',
     .           XLIMS1(1,J),YLIMS1(1,J),ZLIMS1(1,J)
                WRITE (iunout,*) 'XLIMS2,YLIMS2,ZLIMS2 ',
     .           XLIMS2(1,J),YLIMS2(1,J),ZLIMS2(1,J)
            ELSEIF (RLB(J).LT.0) THEN
              DO I=1,ILIN(J)
                WRITE (iunout,*)
     .           'I,ALIMS,XLIMS,YLIMS,ZLIMS'
                WRITE (iunout,*)
     .            I,ALIMS(I,J),XLIMS(I,J),YLIMS(I,J),ZLIMS(I,J)
              ENDDO
              DO I=1,ISCN(J)
                WRITE (iunout,*) 'I,ALIMS0             ',
     .                            I,ALIMS0(I,J)
                WRITE (iunout,*) 'XLIMS1,YLIMS1,ZLIMS1 ',
     .           XLIMS1(I,J),YLIMS1(I,J),ZLIMS1(I,J)
                WRITE (iunout,*) 'XLIMS2,YLIMS2,ZLIMS2 ',
     .           XLIMS2(I,J),YLIMS2(I,J),ZLIMS2(I,J)
                WRITE (iunout,*) 'XLIMS3,YLIMS3,ZLIMS3 ',
     .           XLIMS3(I,J),YLIMS3(I,J),ZLIMS3(I,J)
              ENDDO
            ENDIF
          ENDIF
 7701   CONTINUE
        CALL EIRENE_LEER(2)
        CALL EIRENE_MASIR2('IGJUM0 ',IGJUM0,1,1,1,1,NLIMPS)
        CALL EIRENE_LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL EIRENE_MASBR2('IGJUM1 ',IGJUM1,0,NLIMPS,1,NLIMPS,NLIMPS,
     .                        NBITS)
        END IF
        CALL EIRENE_LEER(1)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS)
        ELSE
          CALL EIRENE_MASBR2('IGJUM2 ',IGJUM2,0,NLIMPS,1,NLIMPS,NLIMPS,
     .                        NBITS)
        END IF
        CALL EIRENE_LEER(1)
        NSOPT=MIN(NSBOX,NOPTIM)
        IF (NLIMPB >= NLIMPS) THEN
          CALL EIRENE_MASIR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS)
        ELSE
          CALL EIRENE_MASBR2('IGJUM3 ',IGJUM3,0,NOPTIM,1,NSOPT,NLIMPS,
     .                        NBITS)
        END IF
        CALL EIRENE_LEER(1)

        DO 7702 J=1,NSOPT
          CALL EIRENE_MASJ3('J,NLIMII,NLIMIE          ',
     .                       J,NLIMII(J),NLIMIE(J))
 7702   CONTINUE
C
      ENDIF

      CALL EIRENE_DEALLOC_BCKGRND
C
      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(3)

!  allocate and initialise storage for trajectories
cdr  cleanup needed: ntrj option is dead.
cdr  only keep traj(jtrj) for line of sight spectra.
      IF (.NOT.ALLOCATED(TRAJ)) THEN
        ALLOCATE (TRAJ(NCHORI+NTRJ))

        DO JTRJ = 1, NCHORI+NTRJ
          ALLOCATE(TRAJ(JTRJ)%TRJ)
          TRAJ(JTRJ)%TRJ%NCOU_CELL = 0
          NULLIFY(TRAJ(JTRJ)%TRJ%CELLS)
        END DO
      END IF

      IF (NCHORI > 0) THEN
        IF (ANY(NLSTCHR)) CALL EIRENE_SETUP_CHORD_SPECTRA
      END IF

!  determine number of background spectra
      CALL EIRENE_NUM_BACKGROUND_SPECTRA(IADTYP)

      IF (.NOT.ALLOCATED(BACK_SPEC) .AND. (NBACK_SPEC > 0))
     .   ALLOCATE(BACK_SPEC(NBACK_SPEC))

      TPB2=EIRENE_SECOND_OWN()
      if (TRCTIM)
     > write (iunout,*) ' CPU time at end of input ',tpb2-tpb1
      tpb1 = tpb2

      IUNIN = IUNIN_SAVE
      IF (IUSROUT /= 0) CLOSE(IUSROUT)
C
      RETURN
      END SUBROUTINE EIRENE_INPUT

c=======================================================================

      subroutine fix_logical_input(a,l)
cxb: * For gfortran: it does not accept empty field for logicals
c      Fill up empty spaces with 'f'
      use eirmod_comprt
      implicit none
      character*(*) a
      integer l
      integer i,j
      j=1
      do i=1,l !{
        if(a(j:j).eq.' ') a(j:j)='f'
        j=j+1
        if(mod(i,5).eq.0) j=j+1
      end do !}
c      write (iunout,'(a,a)') 'zeile :',trim(a)   !###
      return
      end subroutine fix_logical_input

c=======================================================================

      subroutine fix_integer_input(a,l)
cxb: * For pgf90: it does not accept empty field for integers
c      Fill up empty spaces with '0'
      use eirmod_comprt
      implicit none
      character*(*) a
      integer l
      integer i,j,cr
      j=1
      cr=index(a,char(13))
      do i=1,l !{
        if(a(j:j+5).eq.'     '.or.(cr.ge.j.and.cr.lt.j+6)) then !{
          a(j:j+5) = '     0'
        end if !}
        j=j+6
      end do !}
c      write (iunout,'(a,a)') 'zeile :',trim(a)   !###
      return
      end  subroutine fix_integer_input

c========================================================================

      subroutine eirene_replace_stored
cdr For nonlinear mode:
cdr Restore some selected background data and atomic data
cdr for virtual species, which may have been explicitly saved
cdr in iterative routines eirene_modbgk, eirene_tmstep, eirene_modphot

CDR NOTHING IS DONE IN CASE: nfla_virt=0 .and. nrea_virt=0
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comnnl
      use eirmod_comxs
      use eirmod_comusr, only : lgvac, npls_fix
      use eirmod_comprt, only : iunout

      implicit none
      integer :: ifl, irea, jpls, isw, irei, ircx, irpi, irel
      EXTERNAL :: EIRENE_MASJ2, EIRENE_EXIT_OWN

      if (nstordr < nrad) then
        write (iunout,*)
     .    'ERROR encountered in EIRENE_REPLACE_STORED'
        write (iunout,*) 'Reaction data have been read from fort.13'
        write (iunout,*) 'and should be fed into TABEI,CX,EL..'
        write (iunout,*) 'but no space is provided due to'
        write (iunout,*) 'STORAGE SAVE MODE'
        write (iunout,*) 'Calculation abandoned!'
        CALL EIRENE_EXIT_OWN(1)
      end if

      write (iunout,*) 'replace_stored: data for virt. species '
      call eirene_masj2('nfla_vi,nrea_vi ',nfla_virt,nrea_virt)

      if ((nfla_virt > 0) .and. .not.allocated(lg_store)) then
        write (iunout,*) 'NO STORED DATA AVAILABLE FOR LGVAC '//
     .       'IN EIRENE_REPLACE_STORED'
      else
        do ifl = 1, nfla_virt
          jpls = nfla_ipls(ifl)
          if ((jpls < npls_fix) .or. (jpls > npls)) then
            write (iunout,*) 'WRONG IPLS IN EIRENE_REPLACE_STORED'
            write (iunout,'(a,i6)') 'IPLS = ',jpls
          end if
          lgvac(1:nrad,jpls) = lg_store(1:nrad,ifl)
        end do
      end if

      if ((nrea_virt > 0) .and. (.not.allocated(tab_store) .or.
     .                           .not.allocated(e_store))) then
        write (iunout,*) 'NO STORED DATA AVAILABLE FOR '//
     .       'PARTICLE OR ENERGY RATES IN EIRENE_REPLACE_STORED'
      else
        do irea = 1, nrea_virt
          isw = nfla_iswr(irea)
          select case (isw)
c  data for reaction irei, ircx, irpi or irel
          case (1)
            irei=nfla_ir(irea)
            tabei1(irei,1:nrad) = tab_store(irea,1:nrad,1)
            eelei1(irei,1:nrad) = e_store(irea,1:nrad,1)
          case (3)
            ircx=nfla_ir(irea)
            tabcx3(ircx,1:nrad,1:9) = tab_store(irea,1:nrad,1:9)
            eplcx3(ircx,1:nrad,1:9) = e_store(irea,1:nrad,1:9)
          case (4)
            irpi=nfla_ir(irea)
            tabpi3(irei,1:nrad,1:9) = tab_store(irea,1:nrad,1:9)
            eelpi3(irpi,1:nrad,1:9) = e_store(irea,1:nrad,1:9)
          case (5)
            irel=nfla_ir(irea)
            tabel3(irel,1:nrad,1:9) = tab_store(irea,1:nrad,1:9)
            eplel3(irel,1:nrad,1:9) = e_store(irea,1:nrad,1:9)
cdr  case 6,7: to be done
cdr  for photons we may also need "irrc" data
          end select
        end do
      end if
      return

      end subroutine eirene_replace_stored
