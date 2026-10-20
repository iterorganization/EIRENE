Cdr  Purpose: find storage parameters NPARM for a number of allocatable
c             arrays.
c             Set default for NPARM             (for example NATM=0, no atomic species in this run)
c             Read input file, to find NPARMI,  (for example NATMI)
c             Then set NPARM=MAX(NPARM,NPARMI)  (for example NATM=MAX(NATM,NATMI))
c
c    Later these parameters NPARMI may be modified, but NPARM >= NPARMI must be assured always.
C
!pb  11.12.06:  allow letters 'f' or 't' in case name of fem or tetrahedron
!pb             calculation
!pb  27.12.06:  bug fix: increase NSTS in case of time-dependent mode
!pb  15.01.07:  additional line in input block 4 defining HYDKIN model
!pb  02.03.07:  NUMSEC=4 introduced
!pb  20.03.07:  include input block written by HYDKIN model
!pb  22.03.07:  input for NLFEM and NLTET corrected.
!dr  16.01.14:  default NOPTIM changed from 1 to NRAD (automatically), some printout rearranged
!cd  29.10.14:  reading external file for block 4&5: allow comment lines at the beginning of file
!               (same in find_param.f)
cdr  2.2.15:    nflr renamed to nfr (number of TRIM A_on_B files), now same name as in input.f
cdr             because nflr (common CREF) is later used in RDTRIM and REFDAT with a slightly other meaning.
cdr  Jan 2016:  storage for second dimension only if nlpol=true
cdr             to be tested: storage for nplg, if nlpol=false?
cdr             storage for third dimension only if nltor=true
cdr             to be tested:  storage for nltra, if nltor=false?
cdr             to be done: check for comment lines *... synchronized with input.f?
!pb  June 16:   default for NPLSTI changed from 1 to NPLS
!pb  MAY  16:   nrds -> nrei
cdr  March 17:  NPTRGT printed. May have been changed in call to if0parm, block 14.
CDR  May 2017:  try to fix NSTRAI, NSRFSI, consistent with input.f
cdr             same thing: NCPVI, NCPV (and eliminate old parameters NCOP, NCOPI)
cdr  July 17 :  lmulpl: automatic options for multiple ion temperatures,
cdr                     multiple ion velocities in case of BGK nonlinear collisions
cdr  July 17 :  initialize 2D CFD code coupling parameters NDX,....
c               move NRAD=... after call to if0prm, because of 3D CFD (emc3) coupling
cdr  Jun 18  : various corrections, comments in new (generalized) block 12 options.
cdr            nadv=nadv+10: now out, is contained in more general storage settings.
cdr  Oct. 19 : block 5 card counting to infere the setting of INDPRO(2).
cdr            From then on: completely symmetric options for both background parameter
cdr            sets: Ti (temperature) and V.IN (flow field) in multi-species cases.
cdr            Only unsolved case: NPLS=2. This can mean: two different ion temperatures
cdr            or one ion temperature plus one optional VOL (INDPRO(12) card).
cdr            Default: 1 Ti card and 1 Vol card. If two Ti cards are to be read,
cdr            then one VOL card (which normally is optional) must necessarily be present.
cdr  Nov. 19 : Adopted from ITER branch: additional call eirene_init_cinit
cdr            prior to reading optional CFILE cards (paths to external databases).
cdr            Strictly this call should only be after find_param.f is finished, because allocation of
cdr            storage in module CINIT is possible only after NPLS and NSTRA are known.
cdr            But, apparently some of the information in CINIT is needed earlier, e.g. in
cdr            call to plasma code interface-initialization done in block 14.
cdr            Still fiddling with Ti(ipls) input card counting.
cdr            tbd: call eirene_skip_read_comments: not yet implemented here.
cdr  Aug. 20 : remove PART_NAME, BULK_NAME,... fix NSTRA, NSRFS,..
cmg  Feb 1, 23: added ISPOPT and ESPEC%ISPCOPT to turn on scoring of angle distribution
cmg             additional flag in input deck: ISPSRF ... IDIREC ISPOPT in block 10F
chk  March 24: added ISPOPT=2 for Legendre expansion of angular tallies
chk            in each energy bin, and ISPLDEG for the expansion degree

      SUBROUTINE EIRENE_READ_FIXFORM (NL, SAREA_SAVE, IERROR)
C
C   SET DIMENSIONS (PARAMETERS) FOR ALLOCATABLE ARRAYS (TALLIES, GRIDS, ETC..) AND SET SOME FURTHER DEFAULT VALUES
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
      USE EIRMOD_CGRPTL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPL3D
      USE EIRMOD_CPLOT
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CREF
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES, ONLY: NPRS, NLIDENT,
     >    STRATEGY_UNDEFINED, STRATEGY_EMBARRASS,
     >    STRATEGY_ORIGINAL, STRATEGY_APCAS,
     >    STRATEGY_BALANCED, STRATEGY_IMPROVED,
     >    INPUT_DISTRIBUTION_STRATEGY
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMSPL
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON
      USE EIRMOD_JSON
#ifndef NO_JSON
      USE EIRMOD_IOUSR, ONLY: EIRENE_READ_BLOCK_11_USR
      use json_module           !IGNORE
     .   , only  : json_ck
#else
      USE EIRMOD_IOUSR, ONLY: EIRENE_READ_BLOCK_11_USR_FIXED
#endif
      USE EIRMOD_PRESSURELOOP

      IMPLICIT NONE
C
      TYPE(TEMPERATURE), POINTER :: TEMPCUR
      TYPE(DENSITY), POINTER :: DENCUR
      TYPE(VELOCITY), POINTER :: VELCUR
      TYPE(VOLUMEP), POINTER :: VOLCUR
C
      TYPE(TSURFACE), POINTER :: SURFCUR
      TYPE(REFMODEL), POINTER :: REFCUR, REFLAST
      TYPE(SPEC_REF), POINTER :: SPR

      TYPE REFFILE
        CHARACTER(500) :: RFILE
        TYPE(REFFILE), POINTER :: NEXT
      END TYPE REFFILE

      TYPE(REFFILE), POINTER :: REFFILES, CURFILE

      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(TCONTRIB) :: CNT

      INTERFACE
        SUBROUTINE EIRENE_SLREAC (IR, FILNAM, H123, REAC, CRC,
     .             RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .             RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
! for ADAS type data tables
     .             ELNAME, IZ1, BUNDLING,
cdr  optional parameters
     .             M_popesc, M_upper, M_lower,  ! for internal CR models, line emission etc..
     .             IROW_ESC, ICOL_ESC, POP_ESC,
cdr  optional parameters
     .             M_qext,                      ! for internal CR models, line emission etc..
     .             IROW_EXT, ITAL_EXT,
cdr  optional arguments for FILNAM=CONST options
     .             IFTFL, NCOEF, COEF,
cdr  optional arguments for FILNAM=PHOTONS options
     .             IPRFTYPE)

        USE EIRMOD_PRECISION
        INTEGER,      INTENT(IN) :: IR
        CHARACTER(8), INTENT(IN) :: FILNAM
        CHARACTER(4), INTENT(IN) :: H123
        CHARACTER(LEN=*), INTENT(IN) :: REAC
        CHARACTER(3), INTENT(IN) :: CRC
c  optional input for slreac: TAB2D or ADAS model data
        INTEGER,      INTENT(IN) :: IZ1
        CHARACTER(2), INTENT(IN) :: ELNAME
c  optional input for slreac: CR model data
        INTEGER,      INTENT(IN), OPTIONAL :: M_popesc, M_upper,
     .                                           M_lower
        INTEGER,      INTENT(IN), OPTIONAL :: IROW_ESC(:),
     .                                           ICOL_ESC(:)
        REAL(DP),     INTENT(IN), OPTIONAL :: POP_ESC(:)
c  optional input for slreac: CR model data
        INTEGER,      INTENT(IN), OPTIONAL :: M_qext
        INTEGER,      INTENT(IN), OPTIONAL :: IROW_EXT(:),
     .                                           ITAL_EXT(:,:)
c  optional input for slreac: filnam=CONST option
        INTEGER,      INTENT(IN), OPTIONAL :: IFTFL, NCOEF
        REAL(DP),     INTENT(IN), OPTIONAL :: COEF(9)
c  optional input for slreac: filnam=PHOTONS option
        INTEGER,      INTENT(IN), OPTIONAL :: IPRFTYPE
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: BUNDLING
        INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX,JFEX2MN, JFEX2MX
        REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                              RC2MIN, RC2MAX, FP2(6)
        END SUBROUTINE EIRENE_SLREAC


        subroutine eirene_read_reaclines (ZEILE, iunin, iunout, iread,
     .           ir, filnam, h123, reac, crc, mp, mt, dpp,
     .           rc1min, rc1max, jfex1mn, jfex1mx, fp1,
     .           rc2min, rc2max, jfex2mn, jfex2mx, fp2,
     .           iz, elname, bundling,
     .           m_popesc, m_upper, m_lower,
     .           irow_esc, icol_esc, pop_esc,
     .           m_qext,
     .           irow_ext, ital_ext,
     .           iftfl, ncoef, coef,
     .           iprftype)
         USE EIRMOD_PRECISION
         character(420), intent(inout) :: zeile
         integer, intent(inout) :: iread
         integer, intent(in) :: iunin, iunout
c
         integer, intent(out) :: ir, mp, mt,
     .                           jfex1mn, jfex1mx, jfex2mn, jfex2mx
         real(dp), intent(out) :: dpp, rc1min, rc1max, rc2min, rc2max
         real(dp), intent(out) :: fp1(6), fp2(6)
         character(8), intent(out) :: filnam
         character(4), intent(out) :: h123
         character(50), intent(out) :: reac
         character(3), intent(out) :: crc

cdr  additionaly read, in case filnam = TAB2D, ADAS,...
         integer, intent(out) :: iz
         character(2), intent(out) :: elname
         character(60), intent(out) :: bundling
cdr  additionally read, in case filnam = CR...
         integer, intent(out) :: m_popesc, m_upper, m_lower
         integer, allocatable, intent(out) :: irow_esc(:), icol_esc(:)
         real(dp), allocatable, intent(out) :: pop_esc(:)
         integer, intent(out) :: m_qext
         integer, allocatable, intent(out) :: irow_ext(:), ital_ext(:,:)
cdr  additionally read, in case filnam = CONST
         integer, intent(out) :: iftfl, ncoef
         real(dp), intent(out) :: coef(9)
cdr  additionally read, in case filnam = PHOTONS
         integer, intent(out) :: iprftype
        end subroutine eirene_read_reaclines
      END INTERFACE

      TYPE(TRANSFORM), POINTER :: TRF
C
      INTEGER, INTENT(IN) :: NL
      REAL(DP), INTENT(INOUT) :: SAREA_SAVE(NL)
      INTEGER, INTENT(INOUT) :: IERROR
      REAL(DP) :: AFF(3,3), AFFI(3,3), FP1(6), FP2(6), COEF(9)
      REAL(DP) :: ALR, ROTNRM, RPSDL, XSH, YSH, ZSH, REFNRM,
     .          DPP, SPCMN, SPCMX, SPC_SHIFT,
     .          SPCPLT_X, SPCPLT_Y, SPCPLT_SAME, SPCVX, SPCVY, SPCVZ,
     .          VNORM, ESCD2A, ESCD2M, ESCD2I, ESCD2PH, ESCD2P,
     .          RC1MIN, RC1MAX, RC2MIN, RC2MAX, DUMMY
      REAL(DP), ALLOCATABLE :: POP_ESC(:)
c
C  MULTIPLIER FOR BOTH CPU TIME NTCPU AND MAX NUMBER OF MC HISTORIES NPTS, ....
      REAL(DP) :: AMPTS

      INTEGER :: IHELP(NLIMPS), IPRSF(12), NUMTAL(12)
      INTEGER :: II, NFLGS, NSPZS1, NSPZS2,
     .           NTLVF, NSRF, NTLS, NSP, NTL, ICHORI,
     .           IZ, NSMSTRA1,
     .           JATM, JMOL, JION, JPHOT, JPLS, JSTRA,
     .           IRPTA3, IRPTE2, IRPTE3, IH, IDIMP,
     .           JDUMMY, IRPTA1, IRPTA2, IRPTE1, NLJ, I1, I2, I3,
     .           IREAD, I, ISTS, J, IST, IMSG,
     .           IPOS2, IPOS0, NM, K, IPOS1, NRGEN, INUM, NTLSF,
     .           L, INILGJ, INI, ICO, IS, NTLV, ID, NSC,
     .           INELGJ, NPRCSF, MXL, NSPZV1, NSPZV2, NFLGV,
     .           IPRCSF, IR, MT, MP, NDUMM, NUMSEC, NDUMM1, NDUMM2,
     .           NRTAL1, NFR, IPLN,
     .           NDUMM3, NDUMM4, NRE, IFLR,
     .           ISPSRF, ISPTYP, NSPS, NSPSA, IPTYP, IPSPZ,
     .           IANF, IEND, IDEFLT_SPUT, IDEFLT_SPEZ,
     .           IFILE, ISRFCLL,
     .           IDIREC, ISTCHR, ITOK, IER, ILOGS, IO, ISPOPT, ISPLDEG,
     .           IUNIN_SAVE, NLOGIN,
     .           JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX,
     .           IRC, IOPT, ITAL, IRE,
     .           IADV, NUM_COMPO, NUM_CONTRIB, IDMDL, IND,
     .           ILINE, JCOMP, KCONTR,
     .           IFTFL, NCOEF, ITLV, ITLS,
     .           KK, ELONLY,
     .           M_POPESC, M_UPPER, M_LOWER, M_QEXT,
     .           IPRFTYPE
      INTEGER, ALLOCATABLE :: IROW_ESC(:), ICOL_ESC(:)
      INTEGER, ALLOCATABLE :: IROW_EXT(:), ITAL_EXT(:,:)
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LRPS3D, LRPSCN, LINCL45, LMULPL,
     .           LRDMLTI, LRDMLV, NLSPCSCL_DUM, NLSPCSCL_ON_DUM,
     .           NLSPCSPT_DUM
      LOGICAL, ALLOCATABLE :: LOGRDH(:)
      CHARACTER(420) :: ZEILE, ULINE
#ifndef NO_JSON
!cym/cpg avoid type mismatch
      CHARACTER(kind=json_CK,len=420) :: zeileCk
!cym/cpg end
#endif
      CHARACTER(500) :: FILE, FILE45
      CHARACTER(8) :: FILNAM, VARNAME, SPCNAME
      CHARACTER(4) :: H123, CLAB
      CHARACTER(50) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(400) :: PATH, RFILNM
      CHARACTER(6) :: HANDLE
      CHARACTER(2) :: ELNAME
      CHARACTER(60) :: BUNDLING
cdr  intermediate data for those already set in find_param.f
      integer :: noptim_save
      EXTERNAL :: EIRENE_DEKEY, EIRENE_DEKEYB, EIRENE_FILEPATH_USR,
     .            EIRENE_READ_TOKEN,
     .            EIRENE_READ_TRIANG, EIRENE_READ_TETRA,
     .            EIRENE_ROTADD, EIRENE_SETREF, EIRENE_SETROT,
     .            EIRENE_SETUP_DEFAULT_EMISSIVITY,
     .            EIRENE_SKIP_READ_COMMENT, EIRENE_UPPERCASE,
     .            EIRENE_XSHADD, EIRENE_YSHADD, EIRENE_ZSHADD,
     .            FIX_INTEGER_INPUT, FIX_LOGICAL_INPUT,
     .            EIRENE_LEER, EIRENE_MASAGE, EIRENE_MASPRM,
     .            EIRENE_EXIT_OWN

C
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST ITERATION
C  STEP IN THIS RUN. IITER IS THE ACTUAL ITERATION NUMBER
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST TIMESTEP
C  IN THIS RUN. ITIMV IS THE ACTUAL TIMESTEP NUMBER
C
C  INITIALISE SOME DATA AND SET DEFAULTS
C
      IREAD=1

C
!pb IUNIN set in COMPRT
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
cdr IUNOUT has already been set in COMPRT and in subroutine EIRENE
C

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

cdr LMULPL=TRUE: enforce multiple Ti and V..IN, due to BGK reactions (virt. backgrounds)
      LMULPL = .FALSE.
C
      CALL EIRENE_LEER(2)
C
C  READ TEXT DESCRIBING THE RUN, 100--199
C
      READ (IUNIN,'(A72)') TXTRUN
      WRITE (iunout,'(1X,A)') trim(TXTRUN)
  109 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        IF (ZEILE(1:3).NE.'***') THEN
          WRITE (iunout,'(1X,A)') trim(ZEILE)
#ifndef NO_JSON
!  STORE COMMENT LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(cm_stack,trim(zeile))
          zeileCK=trim(zeile)
          call eirene_push_string_stack(cm_stack,trim(zeileCK))
!cym/cpg end
#endif
        END IF
        GOTO 109
      ELSE
        call fix_integer_input(zeile,10)
        READ (ZEILE,6666) NPRLL,NMODE,NTCPU,NFILE,NDUMM1,NITER,
     .                    NDUMM2,NTIME,NEVNTMAX,NVNMXTR
        NEVNTMAX=1000*NEVNTMAX
      ENDIF
      CALL EIRENE_LEER(1)
      IITER=1
      ITIMV=1

      READ (IUNIN,'(A72)') ZEILE
C......................................................................
C  IS THERE AN ADDITIONAL INPUT LINE FOR STORAGE OPTIMIZATION?

C  THESE DEFAULTS HAVE ALREADY BEEN SET IN FIND_PARAM
!     NOPTIM = ...
!     NOPTM1 = 1
!     NGEOM_USR = 0
!     NCOUP_INPUT = 1
!     NSMSTRA = 1    ADDITIONAL STORAGE FOR SUM OVER STRATA IS MADE AVAILABLE.
!     NSTORAM = 9    FULL STORAGE FOR ATOMIC DATA ARRAYS.
!                    SOME OF THAT CAN BE ELIMINATED BY "ON THE FLY A&M COMPUTATIONS"
!     NGSTAL = 0
!
!     NRTAL1 = 0     ONLY FOR FIND_PARAM, NOT USED ANY FURTHER

      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
C  YES
C THIS LINE: NOPTIM, NSTORAM,.... is read in FIND_PARAM and modified in SET_PARMMOD(1)
!pb do not overwrite here
        noptim_save=noptim
        READ (ZEILE,6666) NOPTIM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA1,NSTORAM,NGSTAL,NRTAL1
        noptim=noptim_save
        READ (IUNIN,'(A72)') ZEILE
      END IF

cdr  repeat the same code as in find_param....?
C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE)
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT)
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0

      NOPTM1 = MAX(NOPTM1,1)
C...........................................................................
c   done with this optional "storage save mode card"
C
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,20)
      READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                  NLERG,NLIDENT,NLONE,NLMOVIE,NLDFST,
     .                  NLRANMAR,NLCASCAD,NLOCTREE,NLWRMSH,NEXVS,
     .                  NLTRIMESH,NLSPCSCL_DUM,NLSPCSCL_ON_DUM,
     .                  NLSOLEDGE,NLSPCSPT_DUM

C  OPTIONAL INPUT CARDS, FOR PATHWAYS AND NAME DEFINITIONS
C                        FOR EXTERNAL DATABASES: AMJUEL, HYDHEL,.....

C  IF DATA FILES ARE NOT ENTERED HERE:
c              THEN BY DEFAULT: THE EXTERNAL DATA FILES ARE EXPECTED IN THE SAME
C              DIRECTORY AS THE ONE FROM WHICH THE RUN IS STARTED.
C
cdr scan for optional CFILE lines: path to external database files
cdr Parse input card for PATH = DBFNAME

cdr  probably identical code as just done in find_param.f, repeated here.
cdr  begin

      READ (IUNIN,'(A420)') ZEILE
      IREAD=1
      I1 = INDEX(ZEILE,'CFILE')
      DO WHILE (I1 /= 0)
        I2 = VERIFY(ZEILE(I1+5:),' ') + I1 + 4
        I3 = SCAN(ZEILE(I2+1:),' ')
        IH = MIN(I3,6)
        HANDLE=REPEAT(' ',6)
        HANDLE(1:IH) = ZEILE(I2:I2+IH-1)
c   cfile card found. Is this one of the permitted external files?
        DO IFILE = 1,NDBNAMES
          IF (INDEX(DBHANDLE(IFILE),HANDLE) /= 0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
c   yes, file type no 'ifile' as stored on dbhandle, in eirmod_cinit.
c   currently: 16 types of files are recognized
          IANF = I2+I3+VERIFY(ZEILE(I2+I3:),' ')-1
          IEND = IANF+SCAN(ZEILE(IANF+1:),' ')-1
          DBFNAME(IFILE)(1:IEND-IANF+1) = ZEILE(IANF:IEND)
          DBFNAME_IN(IFILE) = DBFNAME(IFILE)
          CALL EIRENE_FILEPATH_USR(ZEILE,DBFNAME(IFILE),IANF,IEND)
          LDBREAD(IFILE) = .TRUE.
          WRITE (IUNOUT,*) 'PATH SET FOR FILE ',TRIM(HANDLE)
          WRITE (IUNOUT,*) 'PATH = ',ZEILE(IANF:IEND)
        ELSE
          WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED '
          WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',TRIM(HANDLE),
     .                     ' IGNORED '
        END IF
        READ (IUNIN,'(A420)') ZEILE
        I1 = INDEX(ZEILE,'CFILE')
      END DO
cdr end: identical code. Can go out, if just DBFNAME is kept from find_param.f
cdr      Then: just skip CFILE cards here.


C  READING OF OPTIONAL DATABASE NAME AND PATH-CARDS DONE.

C  CONTINUE WITH MANDATORY INPUT.  IREAD=1 AT THIS POINT

C  READING OF INPUT BLOCK 1 DONE

      NFILEN=EIRENE_IDEZ(NFILE,1,5)
      NFILEM=EIRENE_IDEZ(NFILE,2,5)
      NFILEL=EIRENE_IDEZ(NFILE,3,5)
      NFILEK=EIRENE_IDEZ(NFILE,4,5)
      NFILEJ=EIRENE_IDEZ(NFILE,5,5)
      CALL EIRENE_LEER(2)
      CALL EIRENE_MASAGE('*** 1. DATA FOR OPERATING MODE')
      CALL EIRENE_LEER(1)
      IF (NPRS > 1) THEN
        CALL EIRENE_MASAGE('       PARALLELISATION MODE:')
        SELECT CASE( NPRLL )
        CASE( -1 )
          CALL EIRENE_MASAGE('       MPI USER-DEFINED')
          CALL EIRENE_MASAGE('       WILL BE SELECTED FROM END OF FILE')
          input_distribution_strategy = STRATEGY_UNDEFINED
C       CASE( 0 )
C         Reserved for default, see below
        CASE( 1 )
          CALL EIRENE_MASAGE('       MPI PROPORTIONAL ALLOCATION')
          CALL EIRENE_MASAGE('       UNLESS OVERWRITTEN AT END OF FILE')
          input_distribution_strategy = STRATEGY_ORIGINAL
        CASE( 2 )
          CALL EIRENE_MASAGE('       MPI APCAS STRATEGY')
          CALL EIRENE_MASAGE('       ALL PEs CALCULATE ALL STRATA')
          CALL EIRENE_MASAGE('       UNLESS OVERWRITTEN AT END OF FILE')
          input_distribution_strategy = STRATEGY_APCAS
        CASE( 3 )
          CALL EIRENE_MASAGE('       MPI BALANCED ALLOCATION')
          CALL EIRENE_MASAGE('       UNLESS OVERWRITTEN AT END OF FILE')
          input_distribution_strategy = STRATEGY_BALANCED
        CASE( 4 )
          CALL EIRENE_MASAGE('       MPI IMPROVED ALLOCATION')
          CALL EIRENE_MASAGE('       UNLESS OVERWRITTEN AT END OF FILE')
          input_distribution_strategy = STRATEGY_IMPROVED
        CASE DEFAULT
          CALL EIRENE_MASAGE('       MPI "EMBARRASSINGLY PARALLEL"')
          CALL EIRENE_MASAGE('       UNLESS OVERWRITTEN AT END OF FILE')
          input_distribution_strategy = STRATEGY_EMBARRASS
          NPRLL = 0
        END SELECT
      END IF
      WRITE (IUNOUT,*) '       NUMBER OF PROCESSES NPRS= ',NPRS
      CALL EIRENE_LEER(1)
      IF (NMODE.NE.0) THEN
        CALL EIRENE_MASAGE
     .  ('       EIRENE EXPECTS ADDITIONAL DATA. EXTERNAL CODE')
        CALL EIRENE_MASAGE
     .  ('       INTERFACING ROUTINE INFCOP IS CALLED')
        CALL EIRENE_MASAGE
     .  ('       AT ENTRIES IF0COP (GEOMETRY)')
        CALL EIRENE_MASAGE
     .  ('                  IF1COP (BACKGROUND MEDIUM)')
        CALL EIRENE_MASAGE
     .  ('                  IF2COP (BOUNDARY CONDITIONS)')
        IF (NMODE.GT.0) THEN
          CALL EIRENE_MASAGE('       RETURN DATA TO EXTERNAL CODE:')
          CALL EIRENE_MASAGE('       ENTRIES IF3COP AND IF4COP ARE')
          CALL EIRENE_MASAGE('       CALLED AT THE END OF EACH STRATUM')
          CALL EIRENE_MASAGE('       AND AT THE END OF THE RUN, RESP.')
        ELSE
          CALL EIRENE_MASAGE
     .  ('       NO RETURN OF DATA TO EXTERNAL CODE:')
          CALL EIRENE_MASAGE('       ENTRIES IF3COP AND IF4COP')
          CALL EIRENE_MASAGE('       ARE NOT CALLED')
        ENDIF
        WRITE (iunout,*) '       NMODE= ',NMODE
      ELSE
        CALL EIRENE_MASAGE
     .  ('       EIRENE RUN AS STANDALONE CODE')
        CALL EIRENE_MASAGE
     .  ('       INTERFACING ROUTINE INFCOP IS NOT CALLED')
      ENDIF
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       EIRENE ASSUMES A TOTAL CPU TIME'
      WRITE (iunout,*) '       OF ',NTCPU,' SECONDS'
      CALL EIRENE_LEER(1)

      IF (NFILEN.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA'
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*)
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES'
      ELSEIF (NFILEN.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM'
        WRITE (iunout,*)
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11'
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED'
      ELSEIF (NFILEN.EQ.6) THEN
        WRITE (iunout,*) '       EIRENE SAVES OUTPUT DATA'
        WRITE (iunout,*) '       ON FILES FT10 AND FT11 AFTER'
        WRITE (iunout,*)
     .    '       HAVING COMPUTED THE PARTICLE HISTORIES'
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY'
      ELSEIF (NFILEN.EQ.7) THEN
        WRITE (iunout,*) '       EIRENE READS OUTPUT DATA FROM'
        WRITE (iunout,*)
     .    '       AN EARLIER RUN FROM FILES FT10 AND FT11'
        WRITE (iunout,*) '       FOR THE SUM OVER STRATA TALLIES ONLY'
        WRITE (iunout,*) '       NO NEW HISTORIES ARE COMPUTED'
      ENDIF
      IF (NFILEN.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEM.EQ.1) THEN
        WRITE (iunout,*) '       EIRENE SAVES GEOMETRICAL DATA'
        WRITE (iunout,*) '       ON FILE FT12'
      ELSEIF (NFILEM.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS GEOMETRICAL DATA FROM'
        WRITE (iunout,*) '       AN EARLIER RUN FROM FILE FT12'
      ENDIF
      IF (NFILEM.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEL.GE.6.AND.NFILEL.LE.9) THEN
cdr  reset to nfilel=1,2,3,4, and set nlshrt13 switch
cdr  only use short version of fort.13:
cdr  i.e. store only virtual background npls_fix+1,nplsi
cdr  as well as some A&M data for iteration with these species.
        NFILEL = NFILEL - 5
        NLSHRT13 = .TRUE.
        write (iunout,*) 'Switch to short version of FORT.13'
        write (iunout,*) 'Flag NFILEL reset to ',NFILEL
      ELSEIF (NFILEL.GT.0) THEN
        write (iunout,*) 'Use Long version of FORT.13, NFILEL= ',NFILEL
      END IF
      SELECT CASE (NFILEL)
cdr  write field particle parameters, A&M data, and primary source data
        CASE (1)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE SAVES ONLY PLASMA DATA,'
            WRITE (iunout,*) '       STARTING FROM SPECIES NUMBER'
            WRITE (iunout,*) '       NPLS_FIX+1 (COUPLING VARIABLE,'
            WRITE (iunout,*) '       BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
cdr  read field particle parameters, A&M data, and primary source data
        CASE (2)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NPLS_FIX+1 (COUPLING'
            WRITE (iunout,*) '       VARIABLE, BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          ENDIF
          WRITE (iunout,*) '       FROM FILE FT13 '
cdr  read and write field particle parameters, A&M data, and primary source data
        CASE (3)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NPLS_FIX+1 (COUPLING'
            WRITE (iunout,*) '       VARIABLE, BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          ENDIF
          WRITE (iunout,*) '       FROM FILE FT13 AND'
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       SAVES ONLY PLASMA DATA, STARTING'
            WRITE (iunout,*) '       FROM SPECIES NUMBER NPLS_FIX+1,'
          ELSE
            WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
cdr  read and write field particle parameters, A&M data, but no primary source data
        CASE (4)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE READS (AND EXPECTS) ONLY'
            WRITE (iunout,*) '       PLASMA DATA, STARTING FROM SPECIES'
            WRITE (iunout,*) '       NUMBER NPLS_FIX+1 (COUPLING'
            WRITE (iunout,*) '       VARIABLE, BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE READS PLASMA AND A&M DATA'
          END IF
          WRITE (iunout,*) '       FROM FILE FT13 AND'
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       SAVES ONLY PLASMA DATA, STARTING'
            WRITE (iunout,*) '       FROM SPECIES NUMBER NPLS_FIX+1,'
          ELSE
            WRITE (iunout,*) '       SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
          WRITE (iunout,*) '       SOURCE DISTRIBUTION IS NEWLY'
          WRITE (iunout,*) '       DETERMINED'
        CASE DEFAULT
C  ??
      END SELECT
      IF (NFILEL.NE.0) CALL EIRENE_LEER(1)

      SELECT CASE (NFILEK)
        CASE (1)
          WRITE (iunout,*)
     .      '       EIRENE SAVES DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM ON FILE FT14'
        CASE (2)
          WRITE (iunout,*)
     .      '       EIRENE READS DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM FROM FILE'
          WRITE (iunout,*)
     .      '       FT14, AND CARRIES OUT INPUT MODIFICATIONS THIS RUN'
        CASE (3)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM FROM FILE'
          WRITE (iunout,*)
     .      '       FT14, AND CARRIES OUT INPUT MODIFICATIONS THIS RUN'
          WRITE (iunout,*)
     .      '       EIRENE SAVES NEW DATA FOR RECOMMENDED INPUT'
          WRITE (iunout,*)
     .      '       MODIFICATIONS AND RUN TIME PER STRATUM ON FILE FT14'
          WRITE (iunout,*) '       FOR NEXT RUN'
        CASE (4)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD RUN TIME PER STRATUM FROM FILE'
          WRITE (iunout,*) '       FT14'
        CASE (5)
          WRITE (iunout,*)
     .      '       EIRENE READS OLD RUN TIME PER STRATUM FROM FILE'
          WRITE (iunout,*) '       FT14'
          WRITE (iunout,*)
     .      '       EIRENE SAVES NEW RUN TIME PER STRATUM ON FILE FT14'
          WRITE (iunout,*) '       FOR NEXT RUN'
      END SELECT
      IF (NFILEK.NE.0) CALL EIRENE_LEER(1)

      IF (NFILEJ.EQ.1.AND.NTIME.GT.0) THEN
        WRITE (iunout,*) '       EIRENE SAVES SNAPSHOT POPULATION AT'
        WRITE (iunout,*) '       END OF LAST TIMESTEP ON FILE FT15,'
        WRITE (iunout,*) '       STARTING FROM AN EMPTY INITIAL CENSUS'
      ELSEIF (NFILEJ.EQ.2) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FROM'
        WRITE (iunout,*) '       FILE FT15'
      ELSEIF (NFILEJ.EQ.3.AND.NTIME.GT.0) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       STRATUM NSTRAI+1 FOR FIRST TIMESTEP'
        WRITE (iunout,*) '       FROM FILE FT15'
        WRITE (iunout,*) '       EIRENE SAVES NEW SNAPSHOT POPULATION'
        WRITE (iunout,*) '       AT END OF LAST TIMESTEP ON FILE FT15'
      ELSEIF (NFILEJ.EQ.3.AND.NTIME.EQ.0) THEN
        WRITE (iunout,*) '       EIRENE READS SNAPSHOT POPULATION FOR'
        WRITE (iunout,*) '       FIRST TIMESTEP FROM FILE FT15'
        WRITE (iunout,*) '       NO FURTHER SNAPSHOP PRODUCED'
        WRITE (iunout,*) '       DUE TO NTIME=0'
      ENDIF
      IF (NFILEJ.NE.0) CALL EIRENE_LEER(1)

      IF (NITER.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN ITERATIVE MODE.'
        WRITE (iunout,*) '       ITERATIONS: ',IITER,' TO ',NITER
        WRITE (iunout,*)
     .    '       SUBROUTINE "MODUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       ITERATION'
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN NON-ITERATIVE MODE'
      ENDIF
      CALL EIRENE_LEER(1)
      IF (NTIME.GE.1) THEN
        WRITE (iunout,*) '       EIRENE RUN IN TIME DEP. MODE.'
        WRITE (iunout,*) '       TIME CYCLES: ',ITIMV,' TO ',NTIME
        WRITE (iunout,*)
     .    '       SUBROUTINE "TMSUSR" IS CALLED AFTER EACH'
        WRITE (iunout,*) '       TIME CYCLE'
      ELSE
        WRITE (iunout,*) '       EIRENE RUN IN STATIONARY MODE'
      ENDIF
      CALL EIRENE_LEER(1)

      IF (NEVNTMAX.GT.0) THEN
        WRITE(iunout,*) '        EIRENE LIMITS THE NUMBER OF EVENTS ON',
     .    ' A NEUTRAL TRACE to ',NEVNTMAX
        CALL EIRENE_LEER(1)
      ENDIF
C
C  to be done:
c  printout: which default storage settings have we overruled (consequences thereof?)
C
C
C  READ DATA FOR STANDARD MESH, 200---299
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 2. DATA FOR STANDARD MESH')
      CALL EIRENE_LEER(1)
      IREAD=0
C
  201 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        WRITE (iunout,'(1X,A)') trim(ZEILE)
        CALL EIRENE_LEER(1)
        GOTO 201
      ENDIF
      READ (ZEILE,6666) (INDGRD(J),J=1,3)
C
C INPUT SUB-BLOCK 2A
C
C  RADIAL MESH
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLRAD
      IREAD=0
      IF (NLRAD) THEN
C
* For gfortran: it does not accept empty field for logicals
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,8)
        READ (ZEILE,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,NLFEM,NLTET,
     .                    NLGEN
        NLPLG_IN = NLPLG
        NLFEM_IN = NLFEM

! CHECK GEOMETRY SWITCHES
       ILOGS = 0
       IF (NLSLB) ILOGS = ILOGS + 1
       IF (NLCRC) ILOGS = ILOGS + 1
       IF (NLELL) ILOGS = ILOGS + 1
       IF (NLTRI) ILOGS = ILOGS + 1
       IF (NLPLG) ILOGS = ILOGS + 1
       IF (NLFEM) ILOGS = ILOGS + 1
       IF (NLTET) ILOGS = ILOGS + 1
       IF (NLGEN) ILOGS = ILOGS + 1
       IF ((ILOGS == 0) .OR. (ILOGS > 1)) THEN
         WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION'
         WRITE (IUNOUT,*) ' ONE AND ONLY ONE OF THE FOLLOWING',
     .                    ' FLAGS MAY BE .TRUE. '
         WRITE (IUNOUT,*) 'NLSLB = ',NLSLB
         WRITE (IUNOUT,*) 'NLCRC = ',NLCRC
         WRITE (IUNOUT,*) 'NLELL = ',NLELL
         WRITE (IUNOUT,*) 'NLTRI = ',NLTRI
         WRITE (IUNOUT,*) 'NLPLG = ',NLPLG
         WRITE (IUNOUT,*) 'NLFEM = ',NLFEM
         WRITE (IUNOUT,*) 'NLTET = ',NLTET
         WRITE (IUNOUT,*) 'NLGEN = ',NLGEN
         CALL EIRENE_EXIT_OWN(1)
       END IF

       READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
       NR1ST_IN = NR1ST
       IF (NR1ST < 0) THEN
         NR1ST = N1ST
       ENDIF
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,6664) RIA,RGA,RAA,RRA
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,6664) EP1IN,EP1OT,EP1CH,EXEP1
              READ (IUNIN,6664) ELLIN,ELLOT,ELLCH,EXELL
              IF (NLTRI) THEN
                READ (IUNIN,6664) TRIIN,TRIOT,TRICH,EXTRI
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ(IUNIN,6664) XPCOR,YPCOR,ZPCOR,PLREFL
            READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),J=1,NRPLG)
  212       CONTINUE
            IF (PLREFL.GT.0.D0) NR1ST=NR1ST+1
          ENDIF
          IF (NLFEM .OR. NLTET) THEN
            READ (IUNIN,'(A72)') ZEILE
            IREAD = 1
            CLAB = ZEILE(1:4)
            CALL EIRENE_UPPERCASE(CLAB)
            IF (INDEX(ZEILE,'CASE')==0) THEN
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION'
              WRITE (IUNOUT,*)
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON'
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED'
              CALL EIRENE_EXIT_OWN(1)
            END IF

            READ (ZEILE(6:),'(A66)') CASENAME
            IREAD = 0
            CASENAME=ADJUSTL(CASENAME)
            I2=INDEX(CASENAME,' ')

            IF (NLFEM) CALL EIRENE_READ_TRIANG (CASENAME(1:I2))
            IF (NLTET) CALL EIRENE_READ_TETRA (CASENAME(1:I2))

            READ (IUNIN,'(A72)') ZEILE
            IREAD = 1
            IF ( (ZEILE(1:1) .NE. '*') .AND.
C  DOES THE NEXT LINE START WITH A 'LOGICAL' INPUT
     .           (INDEX('FTft',ZEILE(1:1)) == 0) ) THEN
C  NO: READ COORDINATE TRANSLATION CARD
              READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
              IREAD = 0
            END IF
          END IF
          IF (NLGEN) THEN
            NRGEN=NR1ST
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE

CDR  SKIP READING COMMENT INPUT CARDS STARTING WITH *
          do while (zeile(1:1) == '*')
            WRITE (iunout,'(1X,A)') trim(ZEILE)
            READ (IUNIN,'(A72)') ZEILE
          end do

          IREAD = 1
          IPOS1=INDEX(ZEILE,'T')
          IPOS2=INDEX(ZEILE,'F')
          IF (IPOS1.GT.0.OR.IPOS2.GT.0) THEN
            WRITE (iunout,*) 'ONE INPUT LINE MISSING IN BLOCK 2A'
            WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED'
            CALL EIRENE_LEER(1)
            READ (ZEILE,'(L1)') NLPOL
            IREAD = 0
            GOTO 222
          ENDIF
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (ZEILE,6664) RIA,RGA,RAA
            IREAD = 0
          ELSEIF (NLPLG.OR.NLFEM.OR.NLTET) THEN
            READ (ZEILE,6664) XPCOR,YPCOR,ZPCOR
            IREAD = 0
          ENDIF
        ENDIF
      ENDIF
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B

      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)

      READ (ZEILE,6665) NLPOL
      IREAD = 0
C
  222 CONTINUE
      NLPOL_IN = NLPOL
      READ (IUNIN,'(A72)') ZEILE
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLPLY,NLPLA,NLPLP
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      NP2ND_IN = NP2ND
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ELSEIF (INDGRD(2).EQ.6) THEN
      ENDIF
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
      IREAD=0
      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)

      READ (ZEILE,6665) NLTOR
      IREAD=0
C
      READ (IUNIN,'(A80)') ZEILE
* For gfortran: it does not accept empty field for logicals
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      NT3RD_IN = NT3RD
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
      IREAD=0
      CALL EIRENE_SKIP_READ_COMMENT(IREAD,IUNIN,ZEILE)
      READ (ZEILE,6665) NLMLT
      IREAD=0
C
      IF (NLMLT) THEN
        READ (IUNIN,6666) NBMLT
        READ (IUNIN,6664) (VOLCOR(NM),NM=1,NBMLT)
        READ (IUNIN,'(A72)') ZEILE
      ELSE
        NBMLT=1
        VOLCOR(1)=1.D0
C  FIND START OF NEXT INPUT BLOCK: 2E. SEARCH FOR *, T OR F
  241   READ (IUNIN,'(A72)') ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
        IPOS0=INDEX(ZEILE,'*')
        IPOS1=INDEX(ZEILE,'T')
        IPOS2=INDEX(ZEILE,'F')
        IF (IPOS0.EQ.0.AND.IPOS1.EQ.0.AND.IPOS2.EQ.0) GOTO 241
      ENDIF
      IREAD=1
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
C INPUT SUB-BLOCK 2E
C
  250 IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
        IREAD=0
        GOTO 250
      ENDIF

      READ (ZEILE,6665) NLADD
      IREAD=0
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        READ (IUNIN,6664) (VOLADD(NM),NM=1,NRADD)
      ELSE
C  FIND START OF NEXT INPUT BLOCK: 3A
  252   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .NE. '***') GOTO 252
        IREAD=1
      ENDIF
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
      NBMLT=MAX0(1,NBMLT)
      NR1ST=MAX0(1,NR1ST)
      NP2ND=MAX0(1,NP2ND)
      IF (.NOT.NLPOL) NP2ND=1
      NT3RD=MAX0(1,NT3RD)
      IF (.NOT.NLTOR) NT3RD=1
C
      IF (NLPLG.AND.INDGRD(1).LE.4) THEN
        IF (NLPOL.AND.NRPLG.NE.NP2ND) GOTO 994
      ENDIF
C
C  * 3A: READ DATA FOR NON-DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     . ('*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES')
  310 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 310
      IREAD=1
      READ(ZEILE,6666) NSTSI
      IREAD=0
      WRITE (iunout,*) '        NSTSI= ',NSTSI
      CALL EIRENE_LEER(1)
      IMSG=0
      DO 311 ISTS=1,NSTSI
        NLJ=NLIM+ISTS
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSFL(NLJ)
        ELSE
          READ (ZEILE,'(A72)') TXTSFL(NLJ)
          IREAD=0
        ENDIF
        READ (IUNIN,6666) JDUMMY,IDIMP,INUMP(ISTS,IDIMP),IRPTA1,
     .                    IRPTE1,IRPTA2,IRPTE2,IRPTA3,IRPTE3

        IF (.NOT.(NLFEM.OR.NLTET.OR.NLGEN)) THEN

          IF ((IDIMP == 1) .AND. (INUMP(ISTS,IDIMP) > N1ST)) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF RADIAL SURFACE > N1ST'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF ((IDIMP == 2) .AND. (INUMP(ISTS,IDIMP) > N2ND)) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF POLOIDAL SURFACE > N2ND'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          ELSEIF ((IDIMP == 3) .AND.
     .            ((NLTOR.AND.(INUMP(ISTS,IDIMP) > N3RD)) .OR.
     .             (NLTRA.AND.(INUMP(ISTS,IDIMP) > NTTRA)))) THEN
            WRITE (iunout,*) ' ERROR IN SPECIFICATION OF NON-DEFAULT'
            WRITE (iunout,*) ' SURFACE ',ISTS
            WRITE (iunout,*) ' NUMBER OF TOROIDAL SURFACE > N3RD'
            WRITE (iunout,*) ' CHECK INPUT FILE'
            CALL EIRENE_EXIT_OWN(1)
          END IF

        END IF
C
C  OLD INPUT VERSION BEGIN
        IF (IDIMP.EQ.1.AND.IRPTA1.NE.IRPTE1) THEN
          WRITE (iunout,*) 'WARNING FROM INPUT BLOCK 3A, ISTS= ',ISTS
          WRITE (iunout,*) 'NEW INPUT FOR IRPTA,IRPTE....'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
          IMSG=1
          IRPTA2=IRPTA1
          IRPTE2=IRPTE1
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.IRPTE1.AND.
     .          IRPTA2+IRPTE2+IRPTA3+IRPTE3.EQ.0) THEN
          IRPTA2=IRPTA1
          IRPTE2=IRPTE1
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.0.AND.IRPTE1.EQ.0) THEN
          IRPTA1=INUMP(ISTS,1)
          IRPTE1=INUMP(ISTS,1)
        ELSEIF (IDIMP.EQ.2.AND.IRPTA2.EQ.0.AND.IRPTE2.EQ.0) THEN
          IRPTA2=INUMP(ISTS,2)
          IRPTE2=INUMP(ISTS,2)
        ELSEIF (IDIMP.EQ.3.AND.IRPTA3.EQ.0.AND.IRPTE3.EQ.0) THEN
          IRPTA3=INUMP(ISTS,3)
          IRPTE3=INUMP(ISTS,3)
        ENDIF
C  OLD INPUT VERSION DONE
C
C  OVERWRITE DEFAULTS FOR IRPTA, IRPTE ARRAYS
        IF (IDIMP.NE.1) THEN
          IF (IRPTA1.GT.1) IRPTA(ISTS,1)=IRPTA1
          IF (IRPTE1.LE.1.OR.IRPTE1.GT.NR1ST) THEN
            IRPTE(ISTS,1)=MAX(2,NR1ST)
          ELSE
            IRPTE(ISTS,1)=IRPTE1
          ENDIF
          IF ((IDIMP.EQ.2.AND.IRPTE2.GT.NP2ND).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
            IMSG=1
          ENDIF
        ELSE ! IDIMP .EQ. 1
          IRPTA(ISTS,1) = IRPTA1
          IRPTE(ISTS,1) = IRPTE1
        ENDIF
        IF (IDIMP.NE.2) THEN
          IF (IRPTA2.GT.1) IRPTA(ISTS,2)=IRPTA2
          IF (IRPTE2.LE.1.OR.IRPTE2.GT.MAX(NP2ND,NRPLG)) THEN
            IRPTE(ISTS,2)=MAX(2,NP2ND,NRPLG)
          ELSE
            IRPTE(ISTS,2)=IRPTE2
          ENDIF
          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.3.AND.IRPTE3.GT.NT3RD)) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
            IMSG=1
          ENDIF
        ELSE ! IDIMP .EQ. 2
          IRPTA(ISTS,2) = IRPTA2
          IRPTE(ISTS,2) = IRPTE2
        ENDIF
        IF (IDIMP.NE.3) THEN
          IF (IRPTA3.GT.1) IRPTA(ISTS,3)=IRPTA3
          IF (IRPTE3.LE.1.OR.IRPTE3.GT.NT3RD) THEN
            IRPTE(ISTS,3)=MAX(2,NT3RD)
          ELSE
            IRPTE(ISTS,3)=IRPTE3
          ENDIF

          IF ((IDIMP.EQ.1.AND.IRPTE1.GT.NR1ST).OR.
     .        (IDIMP.EQ.2.AND.IRPTE2.GT.MAX(NP2ND,NRPLG))) THEN
            WRITE (iunout,*) 'WARNING:'
            WRITE (iunout,*) 'NON-DEFAULT STANDARD SURFACE NO. ',ISTS
            WRITE (iunout,*) 'OUT OF RANGE. INVISIBLE!'
            IMSG=1
          ENDIF
        ELSE ! IDIMP .EQ. 3
          IRPTA(ISTS,3) = IRPTA3
          IRPTE(ISTS,3) = IRPTE3
        ENDIF

        READ (IUNIN,6666) ILIIN(NLJ),ILSIDE(NLJ),ILSWCH(NLJ),
     .                    ILEQUI(NLJ),ILTOR(NLJ),ILCOL(NLJ),
     .                    ILFIT(NLJ),ILCELL(NLJ),ILBOX(NLJ),
     .                    ILPLG(NLJ)
        ILCOL_IN(NLJ) = ILCOL(NLJ)
        IF (ABS(ILCOL(NLJ)).EQ.7) THEN
          WRITE (iunout,*) 'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',
     .                      NLJ
          WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
          WRITE (iunout,*)
     .      'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
          ILCOL(NLJ)=ILCOL(NLJ)-2
          IMSG=1
        ENDIF
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  NO LOCAL SURFACE INTERACTION MODEL FOUND, USE: DEFAULT
          GOTO 314
C  ASSIGN ONE OF THE SURFACE MODELS (BLOCK 6) TO THIS SURFACE
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = NLJ
cdr increase counter on surflist
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          SMOD_NAME(NLJ)=SURFCUR%MODNAME
          GOTO 314
C  READ LOCAL (FOR THIS SURFACE) SURFACE INTERACTION MODEL, NEXT 3 OR 4 INPUT CARDS
        ELSE
CVK LCHSPNWL
          READ (ZEILE,6666) ILREF(1,NLJ),ILSPT(NLJ),ISRS(1,NLJ),
     .                      ISRC(1,NLJ),LCHSPNWL(1,NLJ)
          IREAD=0
          READ (IUNIN,6664) ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
     .                      TRANSP(1,1,NLJ),TRANSP(1,2,NLJ),FSHEAT(NLJ)
          READ (IUNIN,6664) RECYCF(1,NLJ),RECYCT(1,NLJ),RECPRM(1,NLJ),
     .                      EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*') THEN
            IF (LCHSPNWL(1,NLJ).NE.0) THEN
cxpb Assume old SOLPS4.3 input format
              READ (ZEILE,6664)
     .         RECYCS(1,NLJ),RECYCC(1,NLJ),
     .         ESPUTC(1,NLJ),DUMMY        ,SPTPRM(1,NLJ)
            ELSE
cxpb Assume 3.0.1 version
              READ (ZEILE,6664)
     .         RECYCS(1,NLJ),RECYCC(1,NLJ),SPTPRM(1,NLJ),
     .         ESPUTS(1,NLJ),ESPUTC(1,NLJ)
              LCHSPNWL(1,NLJ)=1
            ENDIF
            IREAD=0
          ELSEIF (ILSPT(NLJ).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING AT NON-DEF. SURFACE ',
     .                        ISTS
            WRITE (iunout,*)
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
            WRITE (iunout,*)
     .        'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(NLJ)=0
            IMSG=1
          ENDIF
        ENDIF
C
  314   CONTINUE
C
  311 CONTINUE
C
C  * 3B: READ DATA FOR ADDITIONAL SURFACES 350--399
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE('*** 3B. DATA FOR ADDITIONAL SURFACES')
  350 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1).EQ.'*') THEN
        WRITE (iunout,'(1X,A)') trim(ZEILE)
        GOTO 350
      ENDIF
      IREAD=0
      READ (ZEILE,6666) NLIMI
      WRITE (iunout,*) '        NLIMI= ',NLIMI
      CALL EIRENE_LEER(1)
C
      IF (NLIMI.GT.0) THEN
        DO 353 I=1,NLIMI
          IHELP(I)=IGJUM0(I)
  353   CONTINUE
  351   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:3).EQ.'CH0') THEN
          IREAD=0
#ifndef NO_JSON
!  STORE CH0 LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(ch0_stack,zeile(1:72))
          zeileCK=zeile
          call eirene_push_string_stack(ch0_stack,zeileCK(1:72))
!cym/cpg end
#endif
          CALL EIRENE_DEKEY (ZEILE(4:72),IHELP,1,1,1,NLIMPS)
          GOTO 351
        ENDIF
        DO 352 I=1,NLIMI
          IGJUM0(I)=IHELP(I)
  352   CONTINUE
      ENDIF
C
      DO 360 I=1,NLIMI
        IF (IREAD.NE.0) THEN
          READ (ZEILE,'(A72)') TXTSFL(I)
        ELSE
          READ (IUNIN,'(A72)') TXTSFL(I)
        ENDIF
        WRITE (iunout,'(1X,A)') trim(TXTSFL(I))
        IREAD=0
!pb read data for switched off surfaces in order to be able to write
!   data onto json files
!        IF (IGJUM0(I).NE.0) THEN
!  361     READ (IUNIN,'(A72)') ZEILE
!          IREAD=1
!          IF (ZEILE(1:1).EQ.'*') GOTO 369
!          IF (ZEILE(1:9).EQ.'TRANSFORM') GOTO 368
!          GOTO 361
!        ENDIF
C
C   GENERAL SURFACE DATA
  362   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3).EQ.'CH1') THEN
#ifndef NO_JSON
!  STORE CH1 LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!           call eirene_push_string_stack(ch1_stack(i),zeile(1:72))
           zeileCK=zeile
           call eirene_push_string_stack(ch1_stack(i),zeileCK(1:72))
!cym/cpg end
#endif
          IF (NLIMPB >= NLIMPS) THEN
            CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL EIRENE_DEKEYB
     .  (ZEILE(4:72),IGJUM1,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362
        ELSEIF (ZEILE(1:3).EQ.'CH2') THEN
#ifndef NO_JSON
!  STORE CH2 LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(ch2_stack(i),zeile(1:72))
          zeileCK=zeile
          call eirene_push_string_stack(ch2_stack(i),zeileCK(1:72))
!cym/cpg end
#endif
          IF (NLIMPB >= NLIMPS) THEN
            CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPS)
          ELSE
            CALL EIRENE_DEKEYB
     .  (ZEILE(4:72),IGJUM2,0,NLIMPS,I,NLIMPB,NBITS)
          END IF
          GOTO 362
        ELSE
          READ (ZEILE,6664) RLB(I),SAREA_SAVE(I),RLWMN(I),RLWMX(I)
          IF (SAREA_SAVE(I).LE.0.D0) SAREA_SAVE(I)=666.
          READ (IUNIN,6666) ILIIN(I),ILSIDE(I),ILSWCH(I),
     .                      ILEQUI(I),ILTOR(I),ILCOL(I),
     .                      ILFIT(I),ILCELL(I),ILBOX(I),
     .                      ILPLG(I)
          IF (ABS(ILCOL(I)).EQ.7) THEN
            WRITE (iunout,*)
     .        'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',I
            WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOG'
            WRITE (iunout,*)
     .        'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
            ILCOL(I)=ILCOL(I)-2
          ENDIF
C  READ SURFACE COEFFICIENTS
          IF (RLB(I).LT.2.) THEN
            READ (IUNIN,6664) A0LM(I),A1LM(I),A2LM(I),A3LM(I),A4LM(I),
     .                        A5LM(I)
            READ (IUNIN,6664) A6LM(I),A7LM(I),A8LM(I),A9LM(I)
          ELSEIF (RLB(I).LT.3.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
          ELSEIF (RLB(I).LT.5.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
          ELSEIF (RLB(I).LT.7.) THEN
            READ (IUNIN,6664) P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),
     .                        P2(3,I)
            READ (IUNIN,6664) P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),
     .                        P4(3,I)
            READ (IUNIN,6664) P5(1,I),P5(2,I),P5(3,I),P6(1,I),P6(2,I),
     .                        P6(3,I)
          ENDIF
        ENDIF
C  READ BOUNDARY DATA
        IF (RLB(I).GT.0..AND.RLB(I).LT.2.) THEN
          READ (IUNIN,6664) XLIMS1(1,I),YLIMS1(1,I),ZLIMS1(1,I),
     .                      XLIMS2(1,I),YLIMS2(1,I),ZLIMS2(1,I)
        ELSEIF (RLB(I).LE.0.D0) THEN
          IH=NINT(-RLB(I))
          ILIN(I)=EIRENE_IDEZ(IH,1,2)
          ISCN(I)=EIRENE_IDEZ(IH,2,2)
          DO 363 J=1,ILIN(I)
            READ (IUNIN,6664) ALIMS(J,I),XLIMS(J,I),YLIMS(J,I),
     .                        ZLIMS(J,I)
  363     CONTINUE
          DO 364 J=1,ISCN(I)
            READ (IUNIN,6664) ALIMS0(J,I),XLIMS1(J,I),YLIMS1(J,I),
     .                        ZLIMS1(J,I),XLIMS2(J,I),YLIMS2(J,I)
            READ (IUNIN,6664) ZLIMS2(J,I),XLIMS3(J,I),YLIMS3(J,I),
     .                        ZLIMS3(J,I)
  364     CONTINUE
        ENDIF

!  STORE VERTICES OF ADDITIONAL SURFACES FOR OUTPUTING THEM ONTO JSON FILE
!  ORIGINAL VALUES MIGHT BE CHANGED LATER ON
        CALL EIRENE_COPY_ADDSRF(I)
C
C READ LOCAL SURFACE INTERACTION MODEL FOR ADDITIONAL SURFACE NO. I
C
  367   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') THEN
C  SKIP READING LOCAL SURFACE INTERACTION MODEL, USE: DEFAULT
          GOTO 369
        ELSEIF (ZEILE(1:9).EQ.'TRANSFORM') THEN
C  TRANSFORM BLOCK
          GOTO 368
        ELSEIF (ILIIN(I).LE.0) THEN
C  TRANSPARENT: SKIP READING SURFACE INTERACTION MODEL
          GOTO 367
        ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
C  ASSIGN LOCAL REFLECTION FROM A LIST READ BELOW
          ALLOCATE(SURFCUR)
          SURFCUR%MODNAME = TRIM(ADJUSTL(ZEILE(9:)))
          SURFCUR%NOSURF = I
          SURFCUR%NEXT => SURFLIST
          SURFLIST => SURFCUR
          IREAD=0
          SMOD_NAME(I)=SURFCUR%MODNAME
          GOTO 368
        ELSE
C  READ LOCAL REFLECTION MODEL, ADDITIONAL SURFACE
CVK LCHSPNWL
          READ (ZEILE,6666) ILREF(1,I),ILSPT(I),ISRS(1,I),ISRC(1,I),
     .                      LCHSPNWL(1,I)
          IREAD=0
          READ (IUNIN,6664) ZNML(I),EWALL(I),EWBIN(I),
     .                      TRANSP(1,1,I),TRANSP(1,2,I),FSHEAT(I)
          READ (IUNIN,6664) RECYCF(1,I),RECYCT(1,I),RECPRM(1,I),
     .                      EXPPL(1,I),EXPEL(1,I),EXPIL(1,I)
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
          IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:9).NE.'TRANSFORM') THEN
            IF (LCHSPNWL(1,I).NE.0) THEN
cxpb Assume old SOLPS4.3 input style
              READ (ZEILE,6664) RECYCS(1,I),RECYCC(1,I),
     .                          ESPUTC(1,I),DUMMY      ,SPTPRM(1,I)
              ESPUTS(1,I)=0.0_DP
            ELSE
cxpb Assume version 3.0.1 input style
              READ (ZEILE,6664) RECYCS(1,I),RECYCC(1,I),SPTPRM(1,I),
     .                          ESPUTS(1,I),ESPUTC(1,I)
              LCHSPNWL(1,I)=1
            ENDIF
            IREAD=0
          ELSEIF (ILSPT(I).NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING FOR ADD. SURFACE ',I
            WRITE (iunout,*)
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
            WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(I)=0
          ENDIF
        ENDIF
C
  368   IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:9).NE.'TRANSFORM') GOTO 369
C
C  TRANSFORM FROM CONVENIENT COORDINATE SYSTEM TO EIRENE COORDINATE
C  SYSTEM. THIS IS POSSIBLE FOR ALL SURFACES READ BY NOW, I.E. FROM
C  ITINI=1 TO ITEND=I
        IREAD=0
        READ (IUNIN,6666) ITINI(I),ITEND(I)
        IF (ITINI(I).LT.1) ITINI(I)=I
        IF (ITEND(I).GT.I) ITEND(I)=I
        READ (IUNIN,6664) XLCOR(I),YLCOR(I),ZLCOR(I)
        READ (IUNIN,6664) XLREF(I),YLREF(I),ZLREF(I)
        READ (IUNIN,6664) XLROT(I),YLROT(I),ZLROT(I),ALROT(I)
C
c  STORE TRANSFORMATIONS FOR OUTPUT TO JSON FILE
        ALLOCATE(TRF)
        TRF%ITINI = ITINI(I)
        TRF%ITEND = ITEND(I)
        TRF%XLCOR = XLCOR(I)
        TRF%YLCOR = YLCOR(I)
        TRF%ZLCOR = ZLCOR(I)
        TRF%XLREF = XLREF(I)
        TRF%YLREF = YLREF(I)
        TRF%ZLREF = ZLREF(I)
        TRF%XLROT = XLROT(I)
        TRF%YLROT = YLROT(I)
        TRF%ZLROT = ZLROT(I)
        TRF%ALROT = ALROT(I)
        NULLIFY(TRF%NEXT)
        IF (.NOT.ASSOCIATED(TRNSFRM(I)%HEAD)) THEN
          TRNSFRM(I)%HEAD => TRF
          TRNSFRM(I)%LAST => TRF
        ELSE
          TRNSFRM(I)%LAST%NEXT => TRF
          TRNSFRM(I)%LAST => TRF
        END IF
C  SHIFT
        XSH=XLCOR(I)
        IF (XSH.NE.0.D0) CALL EIRENE_XSHADD(XSH,ITINI(I),ITEND(I))
        YSH=YLCOR(I)
        IF (YSH.NE.0.D0) CALL EIRENE_YSHADD(YSH,ITINI(I),ITEND(I))
        ZSH=ZLCOR(I)
        IF (ZSH.NE.0.D0) CALL EIRENE_ZSHADD(ZSH,ITINI(I),ITEND(I))
C  REFLECTION
        REFNRM=XLREF(I)*XLREF(I)+YLREF(I)*YLREF(I)+ZLREF(I)*ZLREF(I)
        IF (REFNRM.GT.EPS10) THEN
          CALL EIRENE_SETREF(AFF,AFFI,1,XLREF(I),YLREF(I),ZLREF(I))
          CALL EIRENE_ROTADD(AFF,AFFI,ITINI(I),ITEND(I))
        ENDIF
C  ROTATION
        ROTNRM=XLROT(I)*XLROT(I)+YLROT(I)*YLROT(I)+ZLROT(I)*ZLROT(I)
        ALR=ALROT(I)
        IF (ALR.NE.0..AND.ROTNRM.GT.EPS10) THEN
          CALL EIRENE_SETROT(AFF,AFFI,1,XLROT(I),YLROT(I),ZLROT(I),
     .                                  ALROT(I))
          CALL EIRENE_ROTADD(AFF,AFFI,ITINI(I),ITEND(I))
        ENDIF
        GOTO 368
C
  369   CONTINUE
C
  360 CONTINUE
      CALL EIRENE_LEER(1)

C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
C  AT THIS POINT THE INPUT LINE *** 4.  .... IS EXPECTED
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IMSG=0
      IREAD=0
      IF ((ZEILE(1:3) .NE. '***') .OR. (INDEX(ZEILE,'4.') == 0)) THEN
        WRITE (IUNOUT,*) 'INPUT ERROR, BLOCK *** 4. NOT FOUND'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      CALL EIRENE_MASAGE('*** 4. DATA FOR SPECIES SPECIFICATION AND')
      CALL EIRENE_MASAGE('       ATOMIC PHYSICS MODULE')
      CALL EIRENE_LEER(1)

! CHECK FOR INCLUDE LINE
      READ (IUNIN,'(A420)') ZEILE
!  next input line is already read from iunin, now on 'ZEILE'
      IREAD=1
      ULINE = ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
      I1 = INDEX(ULINE,'INCLUDE')
      LINCL45 = .FALSE.
      IUNIN_SAVE = IUNIN

      IF (I1 > 0) THEN

C  "Include" found. Skip all the rest of block 4 and 5 of input file (fort.iunin),
C   and read this information only from the "include" file instead
C   stream: 2+ifoff
C   Zeile  = INCLUDE 'FILE45'
C
        LINCL45 = .TRUE.
C
        FILE(1:420-I1-6) = ZEILE(I1+7:420)
        FILE(420-I1-5:420) = REPEAT(' ',I1+6)
        CALL EIRENE_READ_TOKEN(FILE,' ',FILE45,ITOK,IER,.FALSE.)
        IREAD = 0

        IUNIN_SAVE = IUNIN
        IUNIN = 2+ifoff !  FILE45 is the include file.
!                          read block 4 and 5 from FILE45 (stream fort.2)
!                          rather than from stream fort.iunin
!                          close fort.2 at end of block 5.
        WRITE (IUNOUT,*) 'EXTERNAL A&M INPUT BLOCK 4 AND 5 FOUND'
        WRITE (IUNOUT,*) 'FILE45 = ',TRIM(FILE45)
        CALL EIRENE_LEER(1)
        OPEN (IUNIN,FILE=FILE45,FORM='FORMATTED',ACCESS='SEQUENTIAL')
c  read comment lines on external A&M data file FILE45, stream fort.2
  401   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 401
        IREAD=1
        GOTO 402
      END IF

C   no "Include" found. Read block 4 and 5 from input
C   stream: IUNIN
C
      IF (IREAD == 0) READ (IUNIN,*)

      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'

      READ (IUNIN,'(A72)') ZEILE  ! THIS IS THE FIRST NON-COMMENT LINE
  402 CALL EIRENE_UPPERCASE(ZEILE)

chk.................................................................
C  special only in case of HYDKIN INTERFACE: find string "DEFAULT"
      IEND=INDEX(ZEILE,'DEFAULT')
      IF (IEND > 0) THEN

CDR  "DEFAULT" is currently only supported in proprietary versions of eirene.
        WRITE (IUNOUT,*) 'DEFAULT: OPTION IS NOT READY, EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
CDR
      END IF
chk..................................................................

C  Normal start of reading database A&M processes

      READ (ZEILE,*) NREACI
      WRITE (iunout,*) '       NREACI= ',NREACI
      CALL EIRENE_LEER(1)

      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(0)
!pbcrm
      IREAD = 0
C
  411 IF (IREAD == 0) READ (IUNIN,'(A420)') ZEILE
      IF (ZEILE(1:1).NE.'*') THEN
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!        call eirene_push_string_stack(crs_stack,trim(zeile))
         zeileCK=zeile
         call eirene_push_string_stack(crs_stack,trim(zeileCK))
!cym/cpg end
#endif
C
C  READ ONE REACTION FROM FILE "FILNAM" AT A TIME. Input card is on "ZEILE"
C
        CALL EIRENE_READ_REACLINES (ZEILE, IUNIN, IUNOUT, IREAD,
     .       IR, FILNAM, H123, REAC, CRC, MP, MT, DPP,
     .       RC1MIN, RC1MAX, JFEX1MN, JFEX1MX, FP1,
     .       RC2MIN, RC2MAX, JFEX2MN, JFEX2MX, FP2,
     .       IZ, ELNAME, BUNDLING,
     .       M_popesc, M_upper, M_lower,  !  read from extra cards: filnam =CR_... format
     .       IROW_ESC, ICOL_ESC, POP_ESC,
     .       M_qext,                      !  read from extra cards: filnam =CR_... format
     .       IROW_EXT, ITAL_EXT,
     .       IFTFL, NCOEF, COEF,
     .       IPRFTYPE)                    !  read from extra cards, filnam =PHOTONS

C  SAVE SOME OF THE INPUT FLAGS FOR LATER
C  PROCESSING (MASS SCALING, POTENTIAL ENERGY INCREMENT) IN XSTCX,XSTEI,...

        MASSP(IR)=MP
        MASST(IR)=MT
        DELPOT(IR)=DPP
C
        CALL EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
! additional input card: asymptotics P1
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
! additional input card: asymptotics P2
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
! additional input card read for TAB2D/ADAS format
     .               ELNAME,IZ,BUNDLING,
! (optional) additional input card read CR internal models
     .               M_POPESC,M_UPPER,M_LOWER,
     .               IROW_ESC,ICOL_ESC,POP_ESC,
     .               M_QEXT,
     .               IROW_EXT,ITAL_EXT,
! (optional) for "CONST models"
     .               IFTFL,NCOEF,COEF,
     .               IPRFTYPE)
        GOTO 411
      ENDIF
      IF (TRCAMD) CALL EIRENE_LEER(1)
C
      WRITE (iunout,*)
     . '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      WRITE (iunout,*) '       NATMI= ',NATMI
      CALL EIRENE_LEER(1)
C
      NSPH=NPHOTI
      DO 421 JATM=1,NATMI
        ISPZ=NSPH+JATM
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSA(JATM),NCHARA(JATM),
     .                       NDUMM1,NDUMM2,  !NPRT=1, NCHRGA=0, DEFAULT
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCA(JATM),NFOLA(JATM),NGENA(JATM),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR ATOMIC SPECIES:
        NPRT(ISPZ)=1
        DO 422 K=1,NRCA(JATM)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCD3A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1A = ',ISCD1A(JATM,K)
            WRITE (iunout,*) ' ISCD2A = ',ISCD2A(JATM,K)
            WRITE (iunout,*) ' ISCD3A = ',ISCD3A(JATM,K)
            WRITE (iunout,*) ' ISCDEA = ',ISCDEA(JATM,K)
            WRITE (iunout,*) ' IESTMA = ',IESTMA(JATM,K)
            WRITE (iunout,*) ' IBGKA  = ',IBGKA(JATM,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACA(JATM,K),IBULKA(JATM,K),
     .                        ISCD1A(JATM,K),ISCD2A(JATM,K),
     .                        ISCD3A(JATM,K),ISCD4A(JATM,K),
     .                        ISCDEA(JATM,K),IESTMA(JATM,K),
     .                        IBGKA(JATM,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1A = ',ISCD1A(JATM,K)
            WRITE (iunout,*) ' ISCD2A = ',ISCD2A(JATM,K)
            WRITE (iunout,*) ' ISCD3A = ',ISCD3A(JATM,K)
            WRITE (iunout,*) ' ISCD4A = ',ISCD4A(JATM,K)
            WRITE (iunout,*) ' ISCDEA = ',ISCDEA(JATM,K)
            WRITE (iunout,*) ' IESTMA = ',IESTMA(JATM,K)
            WRITE (iunout,*) ' IBGKA  = ',IBGKA(JATM,K)
          END IF
          KK=IREACA(JATM,K)
          elonly=iswr(kk)
          IF (IBGKA(JATM,K).NE.0 .AND. elonly.ne.5) then
            WRITE (IUNOUT,*)
     .       'REACT. KK: BGK ITERAT. FOR IATM IGNORED. NOT EL'
            WRITE (IUNOUT,*) 'IATM, KK', JATM,KK
            IBGKA(JATM,K)=0
            IMSG=1
          endif
          LMULPL = LMULPL .OR. (IBGKA(JATM,K) /= 0)

          READ (IUNIN,6664) EELECA(JATM,K),EBULKA(JATM,K),
     .                      ESCD1A(JATM,K),ESCD2A,
     .                      FREACA(JATM,K),EDPOTA(JATM,K)
cdr  .                      FLDLMA(JATM,K)  removed, now controlled by negative ngena
          ESCD1A(JATM,K) = ESCD1A(JATM,K)+ESCD2A
          NSC = 0
          IF (ISCD3A(JATM,K) > 0) NSC = 3
          IF (ISCD4A(JATM,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACA(JATM,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACA(JATM,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACA(JATM,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACA(JATM,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE',
     .              MAX(REACDAT(IREACA(JATM,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACA(JATM,K))%NOSEC =
     .              MAX(REACDAT(IREACA(JATM,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACA(JATM,K))%NOSEC = NSC
            END IF
          END IF
  422   CONTINUE
  421 CONTINUE
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
      IMSG=0
      READ (IUNIN,*)
      WRITE (iunout,*)
     . '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      READ (IUNIN,*) NMOLI
      WRITE (iunout,*) '       NMOLI= ',NMOLI
      CALL EIRENE_LEER(1)
      NSPA=NSPH+NATMI
      DO 431 JMOL=1,NMOLI
        ISPZ=NSPA+JMOL
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSM(JMOL),NCHARM(JMOL),
     .                       NPRT(ISPZ),NDUMM,   !NCHRGM=0, DEFAULT
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCM(JMOL),NFOLM(JMOL),NGENM(JMOL),
     .                       NHSTS(ISPZ),LKINDM(JMOL)
        IF (ISRT(ISPZ,1).LT.0) THEN
          WRITE (iunout,*) 'INPUT ERROR IN BLOCK 4B '
          WRITE (iunout,*) 'MOLECULAR SPECIES ',JMOL,':'
          WRITE (iunout,*) 'ISRT LT 0 OPTION IS NOT AVAILABLE ANYMORE'
          WRITE (iunout,*) 'PROBABLY YOU MEAN: ISRT= ',NMOLI+1
          IERROR=IERROR+1
        ENDIF
        DO 432 K=1,NRCM(JMOL)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCD3M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
            IF (ISCD3M(JMOL,K) .GT. 0) THEN
              WRITE (iunout,*) 'WARNING: REACTION KK= ',IREACM(JMOL,K)
              WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
              WRITE (iunout,*) ' ISCD1M = ',ISCD1M(JMOL,K)
              WRITE (iunout,*) ' ISCD2M = ',ISCD2M(JMOL,K)
              WRITE (iunout,*) ' ISCD3M = ',ISCD3M(JMOL,K)
              WRITE (iunout,*) ' ISCDEM = ',ISCDEM(JMOL,K)
              WRITE (iunout,*) ' IESTMM = ',IESTMM(JMOL,K)
              WRITE (iunout,*) ' IBGKM  = ',IBGKM(JMOL,K)
            ENDIF
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACM(JMOL,K),IBULKM(JMOL,K),
     .                        ISCD1M(JMOL,K),ISCD2M(JMOL,K),
     .                        ISCD3M(JMOL,K),ISCD4M(JMOL,K),
     .                        ISCDEM(JMOL,K),IESTMM(JMOL,K),
     .                        IBGKM(JMOL,K)
            IF (ISCD4M(JMOL,K) .GT. 0) THEN
              WRITE (iunout,*) 'WARNING: REACTION KK= ',IREACM(JMOL,K)
              WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
              WRITE (iunout,*) ' ISCD1M = ',ISCD1M(JMOL,K)
              WRITE (iunout,*) ' ISCD2M = ',ISCD2M(JMOL,K)
              WRITE (iunout,*) ' ISCD3M = ',ISCD3M(JMOL,K)
              WRITE (iunout,*) ' ISCD4M = ',ISCD4M(JMOL,K)
              WRITE (iunout,*) ' ISCDEM = ',ISCDEM(JMOL,K)
              WRITE (iunout,*) ' IESTMM = ',IESTMM(JMOL,K)
              WRITE (iunout,*) ' IBGKM  = ',IBGKM(JMOL,K)
            ENDIF
          END IF
CDR  FOR THE TIME BEING:  NONLINEAR BGK REACTIONS MUST BE ELASTIC
          kk=IREACM(JMOL,K)
          elonly=iswr(kk)
          IF (IBGKM(JMOL,K).NE.0 .AND. elonly.ne.5) then
            WRITE (IUNOUT,*)
     .       'REACTION KK: BGK ITERAT. FOR IMOL IGNORED. NOT EL'
            WRITE (IUNOUT,*) 'IMOL, KK', JMOL,KK, elonly
            IBGKM(JMOL,K)=0
            IMSG=1
          endif
cdr  enforce multiple background field temperature arrays, in case of any BGK reactions.
          LMULPL = LMULPL .OR. (IBGKM(JMOL,K) /= 0)

          READ (IUNIN,6664) EELECM(JMOL,K),EBULKM(JMOL,K),
     .                      ESCD1M(JMOL,K),ESCD2M,
     .                      FREACM(JMOL,K),EDPOTM(JMOL,K)
cdr  for backward compatibility: formerly: two KER values, now one total is used.
          ESCD1M(JMOL,K) = ESCD1M(JMOL,K)+ESCD2M
          NSC = 0
          IF (ISCD3M(JMOL,K) > 0) NSC = 3
          IF (ISCD4M(JMOL,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACM(JMOL,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACM(JMOL,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACM(JMOL,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACM(JMOL,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACM(JMOL,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACM(JMOL,K))%NOSEC =
     .              MAX(REACDAT(IREACM(JMOL,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACM(JMOL,K))%NOSEC = NSC
            END IF
          END IF
  432   CONTINUE
  431 CONTINUE
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
      IMSG=0
      READ (IUNIN,*)
      WRITE (iunout,*) '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'
      READ (IUNIN,*) NIONI
      WRITE (iunout,*) '       NIONI= ',NIONI
      CALL EIRENE_LEER(1)
      NSPAM=NSPH+NATMI+NMOLI
      DO 441 JION=1,NIONI
        ISPZ=NSPAM+JION
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSI(JION),NCHARI(JION),
     .                       NPRT(ISPZ),NCHRGI(JION),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCI(JION),NFOLI(JION),NGENI(JION),
     .                       NHSTS(ISPZ),LKINDI(JION)
        DO 442 K=1,NRCI(JION)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCD3I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
            IF (ISCD3I(JION,K) .GT. 0) THEN
              WRITE (iunout,*) 'WARNING: REACTION KK= ',IREACI(JION,K)
              WRITE (iunout,*) ' THREE SECONDARY GROUPS FOR SPECIES ',
     .                      TEXTS(ISPZ)
              WRITE (iunout,*) ' ISCD1I = ',ISCD1I(JION,K)
              WRITE (iunout,*) ' ISCD2I = ',ISCD2I(JION,K)
              WRITE (iunout,*) ' ISCD3I = ',ISCD3I(JION,K)
              WRITE (iunout,*) ' ISCDEI = ',ISCDEI(JION,K)
              WRITE (iunout,*) ' IESTMI = ',IESTMI(JION,K)
              WRITE (iunout,*) ' IBGKI  = ',IBGKI(JION,K)
            ENDIF
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACI(JION,K),IBULKI(JION,K),
     .                        ISCD1I(JION,K),ISCD2I(JION,K),
     .                        ISCD3I(JION,K),ISCD4I(JION,K),
     .                        ISCDEI(JION,K),IESTMI(JION,K),
     .                        IBGKI(JION,K)
            IF (ISCD4I(JION,K) .GT. 0) THEN
              WRITE (iunout,*) 'WARNING: REACTION KK= ',IREACI(JION,K)
              WRITE (iunout,*) ' FOUR SECONDARY GROUPS FOR SPECIES ',
     .                      TEXTS(ISPZ)
              WRITE (iunout,*) ' ISCD1I = ',ISCD1I(JION,K)
              WRITE (iunout,*) ' ISCD2I = ',ISCD2I(JION,K)
              WRITE (iunout,*) ' ISCD3I = ',ISCD3I(JION,K)
              WRITE (iunout,*) ' ISCD4I = ',ISCD4I(JION,K)
              WRITE (iunout,*) ' ISCDEI = ',ISCDEI(JION,K)
              WRITE (iunout,*) ' IESTMI = ',IESTMI(JION,K)
              WRITE (iunout,*) ' IBGKI  = ',IBGKI(JION,K)
            ENDIF
          END IF
          kk=IREACI(JION,K)
          elonly=iswr(kk)
          IF (IBGKI(JION,K).NE.0 .AND. elonly.ne.5) then
            WRITE (IUNOUT,*)
     .       'REACT. KK: BGK ITERAT. FOR JION IGNORED. NOT EL'
            WRITE (IUNOUT,*) 'IION, KK', JION,KK
            IBGKI(JION,K)=0
            IMSG=1
          endif
          LMULPL = LMULPL .OR. (IBGKI(JION,K) /= 0)

          READ (IUNIN,6664) EELECI(JION,K),EBULKI(JION,K),
     .                      ESCD1I(JION,K),ESCD2I,
     .                      FREACI(JION,K),EDPOTI(JION,K)
          ESCD1I(JION,K) = ESCD1I(JION,K)+ESCD2I
          NSC = 0
          IF (ISCD3I(JION,K) > 0) NSC = 3
          IF (ISCD4I(JION,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACI(JION,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACI(JION,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACI(JION,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACI(JION,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACI(JION,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACI(JION,K))%NOSEC =
     .              MAX(REACDAT(IREACI(JION,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACI(JION,K))%NOSEC = NSC
            END IF
          END IF
  442   CONTINUE
  441 CONTINUE
C
      IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
      IMSG=0
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) == '***') THEN
        IREAD=1
        NPHOTI=0
        GOTO 500
      END IF
      IREAD=0
      WRITE (iunout,*)
     . '*** 4D. PHOTONS SPECIES CARDS, NPHOTI SPECIES'
      READ (IUNIN,*) NPHOTI
      WRITE (iunout,*) '       NPHOTI= ',NPHOTI
      CALL EIRENE_LEER(1)
C
      DO 451 JPHOT=1,NPHOTI
        ISPZ=JPHOT
        READ (IUNIN,66666) I,TEXTS(ISPZ),NDUMM1,NDUMM2,
     .                       NDUMM3,NDUMM4,
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCPH(JPHOT),NFOLPH(JPHOT),NGENPH(JPHOT),
     .                       NHSTS(ISPZ)
C  DEFAULTS FOR PHOTONIC SPECIES:
        NPRT(ISPZ)=1
        DO 452 K=1,NRCPH(JPHOT)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCD3PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1PH = ',ISCD1PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD2PH = ',ISCD2PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD3PH = ',ISCD3PH(JPHOT,K)
            WRITE (iunout,*) ' ISCDEPH = ',ISCDEPH(JPHOT,K)
            WRITE (iunout,*) ' IESTMPH = ',IESTMPH(JPHOT,K)
            WRITE (iunout,*) ' IBGKPH  = ',IBGKPH(JPHOT,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACPH(JPHOT,K),IBULKPH(JPHOT,K),
     .                        ISCD1PH(JPHOT,K),ISCD2PH(JPHOT,K),
     .                        ISCD3PH(JPHOT,K),ISCD4PH(JPHOT,K),
     .                        ISCDEPH(JPHOT,K),IESTMPH(JPHOT,K),
     .                        IBGKPH(JPHOT,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1PH = ',ISCD1PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD2PH = ',ISCD2PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD3PH = ',ISCD3PH(JPHOT,K)
            WRITE (iunout,*) ' ISCD4PH = ',ISCD4PH(JPHOT,K)
            WRITE (iunout,*) ' ISCDEPH = ',ISCDEPH(JPHOT,K)
            WRITE (iunout,*) ' IESTMPH = ',IESTMPH(JPHOT,K)
            WRITE (iunout,*) ' IBGKPH  = ',IBGKPH(JPHOT,K)
          END IF
          kk=IREACPH(JPHOT,K)
          elonly=iswr(kk)
cdr  allow bgk flags also for photon reactions (type: PH)
          IF (IBGKPH(JPHOT,K).NE.0 .AND. elonly .NE.5
     .         .AND. elonly .NE.7) then
            WRITE (IUNOUT,*)
     .       'REACT. KK: BGK ITERAT. FOR JPHOT IGNORED. NOT EL OR PH'
            WRITE (IUNOUT,*) 'IPHOT, KK', JPHOT,KK
            IBGKPH(JPHOT,K)=0
            IMSG=1
          endif
          LMULPL = LMULPL .OR. (IBGKPH(JPHOT,K) /= 0)

          READ (IUNIN,6664) EELECPH(JPHOT,K),EBULKPH(JPHOT,K),
     .                      ESCD1PH(JPHOT,K),ESCD2PH,
     .                      FREACPH(JPHOT,K),EDPOTPH(JPHOT,K)
cdr  .                      FLDLMPH(JPHOT,K)  removed. now controlled by negative ngenph
          ESCD1PH(JPHOT,K) = ESCD1PH(JPHOT,K)+ESCD2PH
          NSC = 0
          IF (ISCD3PH(JPHOT,K) > 0) NSC = 3
          IF (ISCD4PH(JPHOT,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACPH(JPHOT,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACPH(JPHOT,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACPH(JPHOT,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACPH(JPHOT,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACPH(JPHOT,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACPH(JPHOT,K))%NOSEC =
     .              MAX(REACDAT(IREACPH(JPHOT,K))%NOSEC,NSC)
            ELSE
              REACDAT(IREACPH(JPHOT,K))%NOSEC = NSC
            END IF
          END IF
  452   CONTINUE
  451 CONTINUE
C
C  READ DATA FOR PLASMA BACKGROUND, 500--599
C
  500 CONTINUE
C
      IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
      IMSG=0
      IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,'(1X,A)') trim(ZEILE)
      CALL EIRENE_MASAGE('*** 5. DATA FOR PLASMA BACKGROUND')
      CALL EIRENE_LEER(1)
C
C  READ BULK IONS SPECIES CARDS
C
      READ (IUNIN,'(A72)') ZEILE
C     WRITE (iunout,'(1X,A)') trim(ZEILE)
      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'
  510 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 510
      READ(ZEILE,6666) NPLSI
      WRITE (iunout,*) '       NPLSI= ',NPLSI
      CALL EIRENE_LEER(1)
      NSPAMI=NSPAM+NIONI
      NSPTOT=NSPAMI+NPLSI
! counter for additional reaction cards possibly needed for "density models".
      IDMDL = 0
c  loop over background species (except: electrons)
      DO 511 JPLS=1,NPLSI
        ISPZ=NSPAMI+JPLS
        READ (IUNIN,66666) I,TEXTS(ISPZ),NMASSP(JPLS),NCHARP(JPLS),
     .                       NPRT(ISPZ),NCHRGP(JPLS),
     .                       ISRF(ISPZ,1),ISRT(ISPZ,1),NUMSEC,
     .                       NRCP(JPLS),NDUMM1,NDUMM2,
     .                       NHSTS(ISPZ),NDUMM4,
     .                       CDENMODEL(JPLS),NRE
CVK DENSITY LIMIT FOR PLASMA SPECIES
        IF(NDUMM1.EQ.0) THEN
c         DENSLIM(JPLS)=10**(MAXEXPONENT(DENSLIM(JPLS))-1)
         DENSLIM(JPLS)=10e30
        ELSE
         DENSLIM(JPLS)=2._DP**NDUMM1
        END IF
        WRITE(iunout,*) TEXTS(ISPZ), DENSLIM(JPLS) !VK
CVK END
        CALL EIRENE_UPPERCASE (CDENMODEL(JPLS))
        DO 512 K=1,NRCP(JPLS)
          IF (NUMSEC .LT. 3) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCDEP(JPLS,K)
          ELSE IF (NUMSEC .EQ. 3) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCD3P(JPLS,K),ISCDEP(JPLS,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1P = ',ISCD1P(JPLS,K)
            WRITE (iunout,*) ' ISCD2P = ',ISCD2P(JPLS,K)
            WRITE (iunout,*) ' ISCD3P = ',ISCD3P(JPLS,K)
            WRITE (iunout,*) ' ISCDEP = ',ISCDEP(JPLS,K)
          ELSE IF (NUMSEC .EQ. 4) THEN
            READ (IUNIN,6666) IREACP(JPLS,K),IBULKP(JPLS,K),
     .                        ISCD1P(JPLS,K),ISCD2P(JPLS,K),
     .                        ISCD3P(JPLS,K),ISCD4P(JPLS,K),
     .                        ISCDEP(JPLS,K)
            WRITE (iunout,*) ' WARNING !!! '
            WRITE (iunout,*)
     .        ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .          TEXTS(ISPZ)
            WRITE (iunout,*) ' ISCD1P = ',ISCD1P(JPLS,K)
            WRITE (iunout,*) ' ISCD2P = ',ISCD2P(JPLS,K)
            WRITE (iunout,*) ' ISCD3P = ',ISCD3P(JPLS,K)
            WRITE (iunout,*) ' ISCD4P = ',ISCD4P(JPLS,K)
            WRITE (iunout,*) ' ISCDEP = ',ISCDEP(JPLS,K)
          END IF
          READ (IUNIN,6664) EELECP(JPLS,K),EBULKP(JPLS,K),
     .                      ESCD1P(JPLS,K),ESCD2P,
     .                      FREACP(JPLS,K),EDPOTP(JPLS,K)
          ESCD1P(JPLS,K) = ESCD1P(JPLS,K)+ESCD2P
c
cdr  deal with non-default number of secondaries, NSC > 2
          NSC = 0
          IF (ISCD3P(JPLS,K) > 0) NSC = 3
          IF (ISCD4P(JPLS,K) > 0) NSC = 4
          IF (NSC > 0) THEN
            IF ((REACDAT(IREACP(JPLS,K))%NOSEC > 0) .AND.
     .          (REACDAT(IREACP(JPLS,K))%NOSEC /= NSC)) THEN
              WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .              ' REACTION ',IREACP(JPLS,K)
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' PREVIOUSLY WAS ',REACDAT(IREACP(JPLS,K))%NOSEC
              WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .              ' NOW IS         ',NSC
              WRITE (IUNOUT,*) ' USE ',
     .              MAX(REACDAT(IREACP(JPLS,K))%NOSEC,NSC),
     .              ' SECONDARIES'
              REACDAT(IREACP(JPLS,K))%NOSEC =
     .              MAX(REACDAT(IREACP(JPLS,K))%NOSEC,NSC)
            ELSE

              REACDAT(IREACP(JPLS,K))%NOSEC = NSC
            END IF
cdr  NSC: number of secondaries. But here: 0,1,2 secondaries all have NSC=0 ??
          END IF
c
  512   CONTINUE

        IF (LEN_TRIM(CDENMODEL(JPLS)) > 0) THEN
          ALLOCATE (TDMPAR(JPLS)%TDM)
          NRE=MAX(NRE,1)
          TDMPAR(JPLS)%TDM%NRE=NRE
#ifndef LEGACYCOMP
          ALLOCATE (TDMPAR(JPLS)%TDM%ISP(NRE),SOURCE=0)
          ALLOCATE (TDMPAR(JPLS)%TDM%ITP(NRE),SOURCE=0)
          ALLOCATE (TDMPAR(JPLS)%TDM%ISTR(NRE),SOURCE=0)
          ALLOCATE (TDMPAR(JPLS)%TDM%IRC(NRE),SOURCE=0)
#else
          ALLOCATE (TDMPAR(JPLS)%TDM%ISP(NRE))
          ALLOCATE (TDMPAR(JPLS)%TDM%ITP(NRE))
          ALLOCATE (TDMPAR(JPLS)%TDM%ISTR(NRE))
          ALLOCATE (TDMPAR(JPLS)%TDM%IRC(NRE))
          TDMPAR(JPLS)%TDM%ISP(1:NRE) = 0
          TDMPAR(JPLS)%TDM%ITP(1:NRE) = 0
          TDMPAR(JPLS)%TDM%ISTR(1:NRE) = 0
          TDMPAR(JPLS)%TDM%IRC(1:NRE) = 0
#endif

C needs trim, check if possible...
          SELECT CASE (CDENMODEL(JPLS))
          CASE ('FORT.13','FTN13')
            IDMDL = IDMDL + 1
            READ (IUNIN,6666) TDMPAR(JPLS)%TDM%ISP(1)
c  default: only for bulk ions
                              TDMPAR(JPLS)%TDM%ITP(1)=4
          CASE ('FORT.10','FTN10')
            IDMDL = IDMDL + 1
            READ (IUNIN,6666) TDMPAR(JPLS)%TDM%ISP(1),
     .                        TDMPAR(JPLS)%TDM%ITP(1),
     .                        TDMPAR(JPLS)%TDM%ISTR(1)
          CASE ('CONSTANT')
            IDMDL = IDMDL + 1
            READ (IUNIN,6664) TDMPAR(JPLS)%TDM%TVAL,
     .                        TDMPAR(JPLS)%TDM%DVAL,
     .                        TDMPAR(JPLS)%TDM%VXVAL,
     .                        TDMPAR(JPLS)%TDM%VYVAL,
     .                        TDMPAR(JPLS)%TDM%VZVAL
          CASE ('MULTIPLY')
            IDMDL = IDMDL + 1
            READ (IUNIN,'(3I6,6x,3E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%DFACTOR,
     .           TDMPAR(JPLS)%TDM%TFACTOR,
     .           TDMPAR(JPLS)%TDM%VFACTOR
                 TDMPAR(JPLS)%TDM%ITP(1)=4
          CASE ('SAHA')
            IDMDL = IDMDL + 1
!PB   TO BE WRITTEN
          CASE ('PLANCK')
            IDMDL = IDMDL + 1
            READ (IUNIN,'(1E12.4)')
     .           TDMPAR(JPLS)%TDM%G_PLANCK
                 TDMPAR(JPLS)%TDM%ITP(1)=4
!dr   unfinished, TO BE WRITTEN
          CASE ('BOLTZMANN')
            IDMDL = IDMDL + 1
            READ (IUNIN,'(3I6,6x,2E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%G_BOLTZ,
     .           TDMPAR(JPLS)%TDM%DELTAE
          CASE ('CORONA')
            IDMDL = IDMDL + 1  !  ONE MORE H.2 REACTION data set
            READ (IUNIN,'(A72)') ZEILE
            IF (VERIFY(ZEILE,'1234567890eEdD+- ') > 0) THEN
              WRITE (IUNOUT,*) ' OLD VERSION OF INPUT FOR DENSITY',
     .                         ' MODEL CORONA USED'
              WRITE (IUNOUT,*) ' NEW INPUT FORMAT:',
     .                         ' REACTIONS SPECIFIED IN BLOCK 4 ONLY '
              WRITE (IUNOUT,*) ' CHECK IN THE EIRENE MANUAL '
              CALL EIRENE_EXIT_OWN(1)
            END IF
            READ (ZEILE,'(4I6,E12.4)')
     .           TDMPAR(JPLS)%TDM%ISP(1),
     .           TDMPAR(JPLS)%TDM%ITP(1),
     .           TDMPAR(JPLS)%TDM%ISTR(1),
     .           TDMPAR(JPLS)%TDM%IRC(1),
     .           TDMPAR(JPLS)%TDM%A_CORONA
cdr  currently only one single contribution allowed. NRE=1 for CORONA model
cdr         NRE=TDMPAR(JPLS)%TDM%NRE = 1
            IRC = TDMPAR(JPLS)%TDM%IRC(1)
            IF (.NOT.REACDAT(IRC)%LRTC) THEN
              WRITE (iunout,*)
     .          ' WRONG REACTION SPECIFIED FOR CORONA MODEL'
              WRITE (iunout,*) ' ONLY H.2 REACTION DATA ARE PERMITTED'
              WRITE (iunout,*) ' IPLS = ',JPLS
              CALL EIRENE_EXIT_OWN(1)
            END IF
          CASE ('COLRAD')
            IDMDL = IDMDL + 1  ! ONE MORE H.11 or H.12 REACTION data set
cdr  NRE components contributing to field particle species JPLS
            NRE=TDMPAR(JPLS)%TDM%NRE
            DO IRE = 1, NRE
              READ (IUNIN,'(A72)') ZEILE
              IF (VERIFY(ZEILE,'1234567890eEdD+- ') > 0) THEN
                WRITE (IUNOUT,*) ' OLD VERSION OF INPUT FOR DENSITY',
     .                           ' MODEL COLRAD USED'
                WRITE (IUNOUT,*) ' NEW INPUT FORMAT:',
     .                           ' REACTIONS SPECIFIED IN BLOCK 4 ONLY '
                WRITE (IUNOUT,*) ' CHECK IN THE EIRENE MANUAL '
                CALL EIRENE_EXIT_OWN(1)
              END IF
cdr  use "reaction" IRC (population coeff.) and density multiplier ISP,ITP,ISTR
              READ (ZEILE,6666)
     .             TDMPAR(JPLS)%TDM%ISP(IRE),
     .             TDMPAR(JPLS)%TDM%ITP(IRE),
     .             TDMPAR(JPLS)%TDM%ISTR(IRE),
     .             TDMPAR(JPLS)%TDM%IRC(IRE)

              IRC = TDMPAR(JPLS)%TDM%IRC(IRE)
              IF (.NOT.REACDAT(IRC)%LOTH) THEN
                WRITE (iunout,*)
     .            ' WRONG REACTION SPECIFIED FOR COLRAD MODEL'
                WRITE (iunout,*) ' ONLY H.11 OR H.12 (OT) REACTION DATA'
                WRITE (IUNOUT,*) ' ARE PERMITTED'
                WRITE (iunout,*) ' IPLS = ',JPLS
                CALL EIRENE_EXIT_OWN(1)
              END IF
            END DO  ! IRE
          CASE DEFAULT
!PB  NOTHING TO BE DONE
          END SELECT
        END IF
  511 CONTINUE   ! nplsi
cdr  additional reaction cards due to density models: idmdl
cdr   nreaci=nreaci=idmdl  ??
      CALL EIRENE_LEER(1)
C
      DO I=1,NSPZ
        CALL EIRENE_UPPERCASE (TEXTS(I))
      END DO

C
  520 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') THEN
C       WRITE (iunout,'(1x,a)') trim(ZEILE)
        GOTO 520
      ENDIF
      CALL EIRENE_MASAGE('*** 5B. PLASMA BACKGROUND DATA')
      CALL EIRENE_LEER(1)
cdr  indpro(2) has already been determined in find_param.
      READ (ZEILE,6666) (INDPRO(J),J=1,12)
      if (indpro(2).ne.indpro2_save) then
        write (iunout,*) 'indpro(2) reset to ',indpro2_save
        call eirene_leer(1)
        indpro(2)=indpro2_save
      endif

      DO J=1,12
C  Indicate "smoothed" input tallies (interpolation from cell vertices into cells)
C  amongst the 12 input tallies read here.
C  formerly: LSMOPRO(J) flag.
C  Smoothing is currently only for tallies 1 to 7.
        IF (ABS(INDPRO(J)) > 100) THEN
          INDPRO(J) = MOD(INDPRO(J),100)
          SELECT CASE(J)
          CASE (1)
            LTESMO = .TRUE.
          CASE (2)
            LTISMO = .TRUE.
          CASE (3)
            LDESMO = .TRUE.
            LDISMO = .TRUE.
          CASE (4)
            LVXSMO = .TRUE.
            LVYSMO = .TRUE.
            LVZSMO = .TRUE.
          CASE (5)
            LBXSMO = .TRUE.
            LBYSMO = .TRUE.
            LBZSMO = .TRUE.
            LBFSMO = .TRUE.
          CASE (6)
            LADSMO = .TRUE.
          CASE (7)
            LEXSMO = .TRUE.
            LEYSMO = .TRUE.
            LEZSMO = .TRUE.
            LEFSMO = .TRUE.
            LPOTSMO = .TRUE.     ! smoothed electric potential
          END SELECT
        END IF
      ENDDO

C  Te profile
      IF (INDPRO(1).LE.5.AND.NPLSI.GT.0)
     .  READ (IUNIN,6664) TE0,TE1,TE2,TE3,TE4,TE5

C  Ti profile(s)
cdr  default is: eV units for Ti
      NLKELV=INDPRO(2) < 0  ! Kelvin units instead of eV
      INDPRO(2)=IABS(INDPRO(2))

cdr  default is: |indpro(2)| < 10: nplsi Ti fields, one per background species
cdr              |indpro(2)| > 10: only one common Ti field is specified

      NPLSTI = NPLSI
      LRDMLTI=.TRUE.

cdr  only one common Ti profile for all NPLSI species?
      IF (INDPRO(2) > 9) THEN
        NPLSTI = 1
        LRDMLTI=.FALSE.
      ENDIF
      WRITE (iunout,*) '      NPLSTI= ',NPLSTI

      IF ((NPLS > 1) .AND. (NPLSTI == 1) .AND.
     .    (LMULPL .OR. ((NPLS_FIX > 0) .AND. (NPLS_FIX /= NPLSI))
     .            .OR.  ANY(CDENMODEL == FORT//'13'))) THEN

        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'WARNING !'
        WRITE (IUNOUT,*) 'TIIN STORAGE PROVIDED FOR ONE SPECIES ONLY ',
     .                   'DUE TO INDPRO(2) > 10'
        WRITE (IUNOUT,*) 'STORAGE FOR TIIN OVERWRITTEN ',
     .                   'BECAUSE BGK REACTIONS PRESENT ',
     .                   'OR DENSITY MODEL ',FORT//'13',' WAS FOUND'
        LRDMLTI=.FALSE.  ! read only one common Ti card,
                         ! despite storage for NPLS Ti profiles.
        NPLSTI = NPLSI
        WRITE (IUNOUT,*) '      NPLSTI= ',NPLSTI
        CALL EIRENE_LEER(1)
c     ELSE
c  read one Ti card for each of the NPLS Ti profiles
c       LRDMLTI=.TRUE.
      END IF

      IF (INDPRO(2) > 9) INDPRO(2) = MOD(INDPRO(2),10)

      NLMLTI = (NPLSTI > 1)

cdr  indirect addressing: IVL= MPLSTI(IPLS)   ! find Ti field for IPLS on TIIN(IVL)
      IF (NPLSTI == 1) THEN
        MPLSTI = 1                  ! find ALL Ti fields
                                    ! on IPLS=1 storage
      ELSEIF (NPLSTI == NPLS) THEN
        MPLSTI = (/ (I,I=1,NPLS) /) ! find each individual Ti field
                                    ! on its IPLS storage
      ELSE
cdr     more general options: to be written.
cdr     E.g.: only virtual background species with different flow fields
      ENDIF

      MPLSTI=1
      IF (NLMLTI)     MPLSTI = (/ (I,I=1,NPLS) /)

      IF (INDPRO(2).LE.5.AND.NPLSI.GT.0) THEN
        IF (LRDMLTI) THEN
          READ (IUNIN,6664) (TI0(I),TI1(I),TI2(I),TI3(I),TI4(I),TI5(I),
     .                       I=1,NPLSI)
        ELSE
C  READ ONLY ONE COMMON ION TEMPERATURE FOR ALL SPECIES
          READ (IUNIN,6664)  TI0(1),TI1(1),TI2(1),TI3(1),TI4(1),TI5(1)
          DO 530 JPLS=2,NPLSTI
            TI0(JPLS)=TI0(1)
            TI1(JPLS)=TI1(1)
            TI2(JPLS)=TI2(1)
            TI3(JPLS)=TI3(1)
            TI4(JPLS)=TI4(1)
            TI5(JPLS)=TI5(1)
  530     CONTINUE
        ENDIF
      ENDIF

c  di profiles
      IF (INDPRO(3).LE.5)
     .  READ (IUNIN,6664) (DI0(I),DI1(I),DI2(I),DI3(I),DI4(I),DI5(I),
     .                     I=1,NPLSI)
c  V_IN profile(s)
cdr  default is: cm/s units for flow field(s)
      NLMACH=INDPRO(4).LT.0  ! Mach number units instead,
                             ! rather than cm/s
      INDPRO(4)=IABS(INDPRO(4))

cdr  default is: |indpro(4)| < 10: nplsi flow fields, one per background species
cdr              |indpro(4)| > 10: only one common flow field

      NPLSV = NPLSI
      LRDMLV=.TRUE.

cdr  only one common profile for all NPLS species?
      IF (INDPRO(4) > 9) THEN
        NPLSV = 1
        LRDMLV=.FALSE.
      ENDIF

      IF ((NPLS > 1) .AND. (NPLSV == 1) .and. lmulpl) THEN
        WRITE (IUNOUT,*) 'WARNING !'
        WRITE (IUNOUT,*) 'V.IN STORAGE PROVIDED FOR ONE SPECIES ONLY ',
     .                   'DUE TO INDPRO(4) > 10'
        WRITE (IUNOUT,*) 'STORAGE FOR VXIN, VYIN, VZIN OVERWRITTEN ',
     .                   'BECAUSE BGK REACTIONS PRESENT'
        LRDMLV=.FALSE.  !  read only one common V.IN card,
                        !  despite storage for NPLSI V.IN profiles.
        NPLSV = NPLSI
        WRITE (IUNOUT,*) ' NPLSV = ',NPLSV
c     ELSE
c  read one V.IN card for each of the NPLSI V.IN profiles
c       LRDMLV=.TRUE.
      END IF

      IF (INDPRO(4) > 9) INDPRO(4) = MOD(INDPRO(4),10)

      NLMLV = (NPLSV > 1)

cdr  indirect addressing: IVL= MPLSV(IPLS)   ! find flow field for IPLS on V_XYZ(IVL)
      IF (NPLSV == 1) THEN
        MPLSV = 1                     ! find ALL flow fields
                                      ! on IPLS=1 storage
      ELSEIF (NPLSV == NPLSI) THEN
        MPLSV = (/ (I,I=1,NPLSI) /)   ! find each individual flow field
                                      ! on its IPLS storage
      ELSE
cdr     more general options: to be written.
cdr     E.g.: only virtual background species with different flow fields
      ENDIF

      IF (INDPRO(4).LE.5.AND.NPLSI.GT.0) THEN
        IF (LRDMLV) THEN
C  READ NPLSI SETS OF INPUT PARAMETERS, ONE FOR EACH SPECIES IPLS
          READ (IUNIN,6664) (VX0(I),VX1(I),VX2(I),VX3(I),VX4(I),VX5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VY0(I),VY1(I),VY2(I),VY3(I),VY4(I),VY5(I),
     .                       I=1,NPLSI)
          READ (IUNIN,6664) (VZ0(I),VZ1(I),VZ2(I),VZ3(I),VZ4(I),VZ5(I),
     .                       I=1,NPLSI)
        ELSE
C  READ ONLY ONE COMMON FLOW FIELD FOR ALL NPLSI SPECIES
          READ (IUNIN,6664) VX0(1),VX1(1),VX2(1),VX3(1),VX4(1),VX5(1)
          READ (IUNIN,6664) VY0(1),VY1(1),VY2(1),VY3(1),VY4(1),VY5(1)
          READ (IUNIN,6664) VZ0(1),VZ1(1),VZ2(1),VZ3(1),VZ4(1),VZ5(1)
          DO JPLS=2,NPLSTI
            VX0(JPLS)=VX0(1)
            VX1(JPLS)=VX1(1)
            VX2(JPLS)=VX2(1)
            VX3(JPLS)=VX3(1)
            VX4(JPLS)=VX4(1)
            VX5(JPLS)=VX5(1)

            VY0(JPLS)=VY0(1)
            VY1(JPLS)=VY1(1)
            VY2(JPLS)=VY2(1)
            VY3(JPLS)=VY3(1)
            VY4(JPLS)=VY4(1)
            VY5(JPLS)=VY5(1)

            VZ0(JPLS)=VZ0(1)
            VZ1(JPLS)=VZ1(1)
            VZ2(JPLS)=VZ2(1)
            VZ3(JPLS)=VZ3(1)
            VZ4(JPLS)=VZ4(1)
            VZ5(JPLS)=VZ5(1)
          ENDDO
        END IF
      ENDIF

c  pitch or B field profile
c                              !  default:
                               !  B FIELD WITH BX=0
      NLPITCH=INDPRO(5).LT.0   !  for 1D parallel B runs:
                               !  B FIELD WITH BY=0

      INDPRO(5)=IABS(INDPRO(5))
      IF (INDPRO(5).LE.5)
     .  READ (IUNIN,6664) B0,B1,B2,B3,B4,B5

c  zi profiles
      IF (INDPRO(11).NE.0.AND.INDPRO(11).LE.5)
     .  READ (IUNIN,6664) (ZI0(I),ZI1(I),ZI2(I),ZI3(I),ZI4(I),ZI5(I),
     .                     I=1,NPLSI)

c  cell volume profile card, OPTIONAL

      IF (INDPRO(12).LE.5) THEN
        READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
        IREAD=1
        IF ((IO /= 0) .OR. (ZEILE(1:3) .EQ. '***')
     .                .OR. (ZEILE(1:3) .EQ. 'OPT')) THEN
          WRITE (iunout,*) 'ONE INPUT LINE (VOL) MISSING IN BLOCK 5'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED'
          CALL EIRENE_LEER(1)
          VL0=0
        ELSE
          READ (ZEILE,6664) VL0,VL1,VL2,VL3,VL4,VL5
          IREAD=0
        ENDIF
      ENDIF

c  check for "OPTIONAL" input cards.
c  Here: explicitly switch on/off input tallies,
c        or activate gradient tallies: INTLOPTS(ITAL)
c  ITAL: -ntali<ital<0, then: ital=iabs(ital).
c  IOPT: turn off or add gradient tally


      IO = 0
      IF (IREAD == 0) READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
      IREAD = 1
      IF ((IO == 0) .AND. (ZEILE(1:3) .NE. '***')) THEN
        CALL EIRENE_UPPERCASE(ZEILE)
        IF (INDEX(ZEILE,'OPTIONAL') > 0) THEN
cdr  read an arbitrary number of tally-option cards ITAL, IOPT
cdr  and store results on INTLOPTS(ITAL).
cdr  Above some switches LTESMO, LDESMO, etc.. have already been set from INDPRO > 100.
cdr  Possible conflicts, inconsistencies.
          DO
            READ (IUNIN,'(A72)',IOSTAT=IO) ZEILE
            IF ((IO /= 0) .OR. (ZEILE(1:3) .EQ. '***')) EXIT
            READ (ZEILE,6666) ITAL, IOPT
            IF ((ITAL >= 0) .OR. (ITAL < -NTALI)) CYCLE
            INTLOPTS(IABS(ITAL)) = IOPT
          END DO
        END IF
      END IF

      IF (LINCL45) THEN
        CLOSE (IUNIN)
        IUNIN = IUNIN_SAVE
        LINCL45 =.FALSE.

        DO
          READ (IUNIN,'(A72)') ZEILE
          IF ((ZEILE(1:3) == '***') .AND.
     ,        (INDEX(ZEILE,'6.') > 0)) EXIT
        END DO
        IREAD = 1
      END IF

C
C  READ DATA FOR REFLECTION MODEL  600--699
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      NULLIFY(REFFILES)
      CALL EIRENE_MASAGE('*** 6. GENERAL DATA FOR REFLECTION MODEL')
      CALL EIRENE_LEER(1)
C
  610 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 610
      READ (ZEILE,6665) NLTRIM
      IREAD=0
      NFR = 0
      IF (NLTRIM) THEN

c  read TRIM reflection datasets A_on_B
        READ (IUNIN,'(A420)') ZEILE
        IREAD=1
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').EQ.0) THEN
C  NO PATH SPECIFIED FOR REFLECTION DATABASE
          WRITE (iunout,*)
     .      ' NO PATH SPECIFIED FOR TRIM REFLECTION DATABASE'
          WRITE (iunout,*) ' OLD TRIM DATABASE VERSION USED'
          CALL EIRENE_LEER(1)
          LTRIM_OLD=.TRUE.
C  TAKE OLD "TRIM.DAT" FILE WITH NHD6=12 TARGET-PROJECTILE COMBINATIONS
C           "TRIM.DAT" IS EXPECTED ON INPUT STREAM IUN=(21+IFOFF) (SUBR. REFDAT.F)
C  SKIP LINES CONTAINING SPECIFICATIONS FOR DATABASES AND
C  CONTINUE READING with DATD
  615     IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').NE.0) THEN
            READ (IUNIN,'(A72)') ZEILE
            GOTO 615
          ELSE
C  NEXT VALID INPUT CARD FOUND
            IREAD=1
            GOTO 620
          ENDIF
        ELSE
C  PATH SPECIFICATION FOR DATABASE FOUND
          LTRIM_OLD=.FALSE.
          READ (ZEILE(7:LEN(ZEILE)),'(A400)') PATH
!VK TO REMOVE \r
          I=INDEX(PATH,ACHAR(13)) !VK
          IF(I.GT.0.AND.I.LE.LEN_TRIM(PATH)) PATH(I:I)=' '
          PATH=ADJUSTL(PATH)
          IF (LEN_TRIM(PATH) == 0) PATH='./ '
          I2=INDEX(PATH,' ')
          IF ((I2 == 0) .AND. (ZEILE(408:408) /= ' ')) THEN
            WRITE (iunout,*) ' PATH FOR TRIM REFLECTION DATABASE IS',
     .                  ' TOO LONG !'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          FILE=REPEAT(' ',420)
          FILE(1:I2-1)=PATH(1:I2-1)
          WRITE (iunout,*) ' PATH = ',FILE(1:I2-1)
C  PATH FOUND. NEXT: READ ONE OR MORE CARDS FILNAM A_ON_B
          READ (IUNIN,'(A72)') ZEILE
  625     IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').NE.0) THEN
            NFR=NFR+1
            READ (ZEILE,'(A72)') RFILNM
!VK TO REMOVE \r
            I=INDEX(RFILNM,ACHAR(13)) !VK
            IF(I.GT.0.AND.I.LE.LEN_TRIM(RFILNM)) RFILNM(I:I)=' '
            FILE(I2:)=RFILNM
            WRITE (iunout,*) ' NFR =',NFR,' FILNAM = ',FILE(I2:I2+20)
C           WRITE (iunout,'(A,A)') ' FILE = ',FILE
            ALLOCATE (CURFILE)
            CURFILE%RFILE = FILE
            CURFILE%NEXT => REFFILES
            REFFILES => CURFILE
            READ (IUNIN,'(A72)') ZEILE
            GOTO 625
          ELSE
            IREAD=1
            GOTO 620
          ENDIF
        ENDIF
      ENDIF

  620 CONTINUE  !  READING OF REFLECTION DATASETS 'A_ON_B' COMPLETED
      IF (NFR.GE.1) CALL EIRENE_LEER(1)

      IF (ASSOCIATED(REFFILES)) THEN
c  TRIM files for NFR target-projectile combinations are requested.
c  read them one by one in subr. RDTRIM
        NHD6 = NFR
      ELSE
c  old default: all TRIM data in one single file.
c  Read this (single) file "TRIM.DAT" in subr. REFDAT
        NHD6 = 12
      END IF

      NH0=NHD1*NHD2*NHD6
      NH1=NH0*NHD3
      NH2=NH1*NHD4
      NH3=NH2*NHD5
      CALL EIRENE_ALLOC_CREF

      IFLR = 0
      DO WHILE (ASSOCIATED(REFFILES))
        IFLR = IFLR + 1
        REFFIL(IFLR) = REFFILES%RFILE
        CURFILE => REFFILES
        REFFILES => REFFILES%NEXT
        DEALLOCATE(CURFILE)
      END DO
cdr AT THIS POINT: IFLR=NFR OR IFLR=0 ??
      IF (IFLR.NE.NFR.AND.IFLR.NE.0) GOTO 993

c  next: read species index sampling distributions datm, dmol, dion, dpls,
c        and, in case nphot > 0, also read dphot

      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      READ (ZEILE,6664) (DATD(JATM),JATM=1,MIN(NATMI,6))
      IREAD=0
      IF (NATMI.GT.6) THEN
        READ (IUNIN,6664) (DATD(JATM),JATM=7,NATMI)
      ENDIF
cdr next: read mod(nmoli,6) cards for molecular species distribution
      READ (IUNIN,6664) (DMLD(JMOL),JMOL=1,NMOLI)
cdr ditto: test ions
      READ (IUNIN,6664) (DIOD(JION),JION=1,NIONI)
cdr ditto: field particles (bulk ions)
      READ (IUNIN,6664) (DPLD(JPLS),JPLS=1,NPLSI)
      IF (NPHOTI > 0)   !dr try to make this more logical:
                        !dr always read dphd.
cdr                     !dr but then: backward compatible ?
     .  READ (IUNIN,6664) (DPHD(JPHOT),JPHOT=1,NPHOTI)

      DATD_IN = DATD
      DMLD_IN = DMLD
      DIOD_IN = DIOD
      DPLD_IN = DPLD
      DPHD_IN = DPHD

c  next: read universal surface reflection model flags
      READ (IUNIN,6664) ERMIN,ERCUT,RPROB0,RINTEG(1),EINTEG(1),AINTEG(1)

c  read surface models identified by character string 'SURFMOD_...'

      DO
        IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') THEN
          IREAD = 1
cdr look for next input block: no. 7
          EXIT
        END IF
        IF (ZEILE(1:8) /= 'SURFMOD_') CYCLE

        ALLOCATE (REFCUR)
cdr reflection/recycling parameters
        ALLOCATE (REFCUR%JLREF(nspz))
        ALLOCATE (REFCUR%JSRS(nspz))
        ALLOCATE (REFCUR%JSRC(nspz))
        ALLOCATE (REFCUR%JLCHSPNWL(nspz))
        ALLOCATE (REFCUR%TRANSPR(nspz,2))
        ALLOCATE (REFCUR%RCYCFR(nspz))
        ALLOCATE (REFCUR%RCYCTR(nspz))
        ALLOCATE (REFCUR%RCPRMR(nspz))
        ALLOCATE (REFCUR%EXPPLR(nspz))
        ALLOCATE (REFCUR%EXPELR(nspz))
        ALLOCATE (REFCUR%EXPILR(nspz))
cdr sputtering parameters
        ALLOCATE (REFCUR%RCYCSR(nspz))
        ALLOCATE (REFCUR%RCYCCR(nspz))
        ALLOCATE (REFCUR%STPRMR(nspz))
        ALLOCATE (REFCUR%ESPTSR(nspz))
        ALLOCATE (REFCUR%ESPTCR(nspz))
        NULLIFY (REFCUR%SPEC_LINES)

        REFCUR%JLREF = 0
        REFCUR%JSRS = 0
        REFCUR%JSRC = 0
        REFCUR%JLCHSPNWL = 0
        REFCUR%TRANSPR = 0._DP
        REFCUR%RCYCFR  = 0._DP
        REFCUR%RCYCTR  = 0._DP
        REFCUR%RCPRMR  = 0._DP
        REFCUR%EXPPLR  = 0._DP
        REFCUR%EXPELR  = 0._DP
        REFCUR%EXPILR  = 0._DP
cdr  default parameters for sputter models
        REFCUR%RCYCSR  = 0._DP
        REFCUR%RCYCCR  = 0._DP
        REFCUR%STPRMR  = 0._DP
        REFCUR%ESPTSR  = 0._DP
        REFCUR%ESPTCR  = 0._DP

        REFCUR%REFNAME = TRIM(ADJUSTL(ZEILE(9:)))
        IREAD=0
cdr read a full "surfmod_..." model
        READ (IUNIN,'(A72)') ZEILE
        call fix_integer_input(zeile,6)
        READ (ZEILE,6666) REFCUR%JLREF(1),REFCUR%JLSPT,
     .                    REFCUR%JSRS(1),REFCUR%JSRC(1),
     .                    REFCUR%JLCHSPNWL(1),REFCUR%REFCELL
        READ (IUNIN,6664) REFCUR%ZNMLR,REFCUR%EWALLR,REFCUR%EWBINR,
     .                    REFCUR%TRANSPR(1,1),REFCUR%TRANSPR(1,2),
     .                    REFCUR%FSHEATR
        READ (IUNIN,6664) REFCUR%RCYCFR(1),REFCUR%RCYCTR(1),
     .                    REFCUR%RCPRMR(1),REFCUR%EXPPLR(1),
     .                    REFCUR%EXPELR(1),REFCUR%EXPILR(1)

cdr This particular "surface model" SURFMOD_... is now copied NSPZ times,
cdr i.e., made identically the same for all possible incident species ISPZ
cdr At this point we do not know if it includes a sputter card or not.
        DO I=2,NSPZ
          REFCUR%JLREF(I) = REFCUR%JLREF(1)
          REFCUR%JSRS(I) = REFCUR%JSRS(1)
          REFCUR%JSRC(I) = REFCUR%JSRC(1)
          REFCUR%JLCHSPNWL(I) = REFCUR%JLCHSPNWL(1)    !VK
          REFCUR%TRANSPR(I,1:2)=REFCUR%TRANSPR(1,1:2)
          REFCUR%RCYCFR(I) = REFCUR%RCYCFR(1)
          REFCUR%RCYCTR(I) = REFCUR%RCYCTR(1)
          REFCUR%RCPRMR(I) = REFCUR%RCPRMR(1)
          REFCUR%EXPPLR(I) = REFCUR%EXPPLR(1)
          REFCUR%EXPELR(I) = REFCUR%EXPELR(1)
          REFCUR%EXPILR(I) = REFCUR%EXPILR(1)
        END DO
C
C  DEFAULT
        REFCUR%RCYCSR=RECYCS(1,0)
        REFCUR%RCYCCR=RECYCC(1,0)
        REFCUR%STPRMR=SPTPRM(1,0)
        REFCUR%ESPTSR=ESPUTS(1,0)
        REFCUR%ESPTCR=ESPUTC(1,0)
C
        READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        ideflt_sput=-1
        ideflt_spez=-1
C  EITHER READ ONE LINE FOR NON-DEFAULT SPUTTER MODEL
C  AND/OR LINES FOR SPECIES DEPENDENT SURFACE MODEL
        do while (ZEILE(1:1).NE.'*'.AND.ZEILE(1:8).NE.'SURFMOD_')
cdr There is one more line indeed.
cdr Tentatively assume: species dep. surface model
          varname = zeile(1:8)
          call EIRENE_uppercase (varname)
          spcname = zeile(10:17)
          call EIRENE_uppercase (spcname)
          ispz=-1
          ico=0
          do is=1,nspz
            if (texts(is) == spcname) then
              if (ico == 0) then
                ispz = is
                ico = 1
              else
                call EIRENE_leer(2)
                write (iunout,*) ' WARNING !!'
                write (iunout,*)
     .           ' ambiguous species names found in surface model "',
     .                       REFCUR%REFNAME,'"'
                write (iunout,*)
     .           varname,' is modified for species ',ispz
              end if
            end if
          end do     ! end of "is"-loop
          ideflt_spez=1   !tentatively assume: species card
          if (ispz > 0) then
            allocate(spr)
            spr%varname = varname
            spr%spcname = spcname
            if ((varname == 'ISRS') .or. (varname == 'ISRC')) then
              read (zeile(18:),'(I6)') spr%ival
            else
              read (zeile(18:),'(E12.4)') spr%rval
            end if
            spr%next => refcur%spec_lines
            refcur%spec_lines => spr
          end if
          select case (varname)
          case ('ILREF')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JLREF(ispz)
          case ('ISRS')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRS(ispz)
          case ('ISRC')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JSRC(ispz)
          case ('LCHSPNWL')
            if (ispz > 0)
     .          read (zeile(18:),'(I6)') REFCUR%JLCHSPNWL(ispz) !VK
          case ('TRANSP1')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,1)
          case ('TRANSP2')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%TRANSPR(ispz,2)
          case ('RECYCF')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCFR(ispz)
          case ('RECYCT')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCTR(ispz)
          case ('RECPRM')
            if (ispz >0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCPRMR(ispz)
          case ('EXPPL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPPLR(ispz)
          case ('EXPEL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPELR(ispz)
          case ('EXPIL')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%EXPILR(ispz)
          case ('RECYCS')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCSR(ispz)
          case ('RECYCC')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%RCYCCR(ispz)
          case ('SPTPRM')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%STPRMR(ispz)
          case ('ESPUTS')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%ESPTSR(ispz)
          case ('ESPUTC')
            if (ispz > 0)
     .          read (zeile(18:),'(E12.4)') REFCUR%ESPTCR(ispz)

          case default
c  ZEILE is NOT a species card, hence assume: it is a sputter model card
            if (ispz < 0) then
              IF (REFCUR%JLCHSPNWL(1).EQ.0) THEN
cxpb Assume 3.0.1 input format
                READ (ZEILE,6664) REFCUR%RCYCSR(1),REFCUR%RCYCCR(1),
     .                            REFCUR%STPRMR(1),REFCUR%ESPTSR(1),
     .                            REFCUR%ESPTCR(1),REFCUR%REFPRESS
              ELSE
cxpb Assume old SOLPS4.3 format
                READ (ZEILE,6664) REFCUR%RCYCSR(1),REFCUR%RCYCCR(1),
     .                            REFCUR%ESPTCR(1),DUMMY,
     .                            REFCUR%STPRMR(1),REFCUR%REFPRESS
              ENDIF
              DO I=2,NSPZ
                REFCUR%RCYCSR(I) = REFCUR%RCYCSR(1)
                REFCUR%RCYCCR(I) = REFCUR%RCYCCR(1)
                REFCUR%STPRMR(I) = REFCUR%STPRMR(1)
                REFCUR%ESPTSR(I) = REFCUR%ESPTSR(1)
                REFCUR%ESPTCR(I) = REFCUR%ESPTCR(1)
              ENDDO
              ideflt_sput=1
              ideflt_spez=-1
              !Pressure feedback loop model
              IF (REFCUR%JLREF(1) == 4)
     .          WRITE (IUNOUT,*) "PFL with parameters: ",
     .            REFCUR%REFCELL, REFCUR%REFPRESS
            end if
          end select
          if ((ideflt_spez > 0) .and. (ispz < 0)) then
            write (iunout,*)
     .        ' wrong card in species dep. reflection model'
            write (iunout,'(1x,a)') trim(zeile)
          end if
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
C
        END DO    ! end of do-while loop (search for * or for SURFMOD_)
C all lines of this particular SURFMOD_... are read now
C
        IF (ideflt_sput.le.0.and.REFCUR%JLSPT.NE.0) THEN
          WRITE (iunout,*) 'WARNING: SPUTTERING FOR MODEL ',
     .     TRIM(REFCUR%REFNAME)
          WRITE (iunout,*) 'BUT NO PARAMETERS RECYCS, RECYCC ARE READ'
          WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
          WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
          REFCUR%JLSPT=0
        ENDIF
        L=LEN_TRIM(REFCUR%REFNAME)
        WRITE (iunout,*) 'SURFACE MODEL "',REFCUR%REFNAME(1:L),
     .                   '" DEFINED'
!pb     REFCUR%NEXT => REFLIST
!pb     REFLIST => REFCUR
!  add new reflection model at the end of the list
        nullify(refcur%next)
        if (.not.associated(reflist)) then
          reflist => refcur
          reflast => refcur
        else
          reflast%next => refcur
          reflast => reflast%next
        end if
      END DO
      CALL EIRENE_LEER(1)
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      CALL EIRENE_MASAGE
     .  ('*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA')
C
  710 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 710
      READ (ZEILE,6666) NSTRAI
      IREAD=0
      WRITE (iunout,*) '       NSTRAI= ',NSTRAI
      CALL EIRENE_LEER(1)
      READ (IUNIN,6666) (INDSRC(IST),IST=1,NSTRAI)
      READ (IUNIN,6664) ALLOC, AMPTS

C  AMPTS: (option added 2014) common multiplier for max. allowed cpu time NTCPU,
C          and for number of histories NPTS (see below)
      IF(AMPTS.EQ.0.0_DP) AMPTS=1.0_DP
      IF(AMPTS.NE.1.0_DP) THEN
        MPTS_COMSOU=AMPTS
        NTCPU=INT(REAL(NTCPU)*AMPTS)
        WRITE (iunout,*) ' NTCPU, NPTS ENHANCED BY FACTOR AMPTS ',
     .                     AMPTS
      ELSE
        MPTS_COMSOU=1.0_DP
      END IF
C
      DO 712 JSTRA=1,NSTRAI
        ! SKIP READING INPUT FOR THIS STRATUM?
        IF (INDSRC(JSTRA).EQ.6) GOTO 712

        I=JSTRA
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') TXTSOU(JSTRA)
        ELSEIF (IREAD.EQ.1) THEN
          READ (ZEILE,'(A72)') TXTSOU(JSTRA)
          IREAD=0
        ENDIF
  713   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 713
        call fix_logical_input(zeile,4)
        READ (ZEILE,6665) NLAVRP(I),NLAVRT(I),NLSYMP(I),NLSYMT(I)
C    .                   ,NLRAY(I)  ! CDR: UNFINISHED PROPRIETARY OPTION
                          NLRAY(I)=.FALSE.
        IREAD=0
        READ (IUNIN,6666) NPTS(I),NINITL(I),NEMODS(I),NAMODS(I),

c.............................................
c  NMINPTS: ENFORCE MINIMUM NUMBER OF HISTORIES FOR STRATUM.
     .                    NMINPTS(I),
c  NPTSDEL: =0 (VALUES NPTSDEL.GT.0 ONLY FOR INTERNAL TESTING PROCEDURES
C               OF MPI parallelization. Requires 2 runs with special input settings)
cdr See comments in MCARLO.F
     .                    NPTSDEL(I)
C    .                   ,NRAYEN(I)   !CDR  NOT IN USE
c...................................................................

C APPLY MULTIPLIER AMPTS TO THE SELECTED MC PARTICLE NUMBER RANGE
        IF(AMPTS.NE.1.0_DP.AND.AMPTS.GT.0.0_DP) THEN
         NPTS(JSTRA)=INT(REAL(NPTS(JSTRA))*AMPTS)
         NMINPTS(JSTRA)=INT(REAL(NMINPTS(JSTRA))*AMPTS)
        END IF
C

csw if npts < 0 --> npts = infinity
        if(npts(jstra).lt.0) npts(jstra) = huge(1)-1

        READ (IUNIN,66662) FLUX(I),SCALV(I),IVLSF(I),
     .                    ISCLS(I),ISCLT(I),ISCL1(I),
     .                    ISCL2(I),ISCL3(I),ISCLB(I),
     .                    ISCLA(I)
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,5)
        READ (ZEILE,6665) NLATM(I),NLMOL(I),NLION(I),NLPLS(I),
     .                    NLPHOT(I)
        READ (IUNIN,6666) NSPEZ(I)
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,5)
        READ (ZEILE,6665) NLPNT(I),NLLNE(I),
     .                    NLSRF(I),NLVOL(I),NLCNS(I)
C  SAME FOR POINT, LINE, SURFACE AND VOLUME SOURCES
        READ (IUNIN,6666) NSRFSI(JSTRA)
        IF (NSRFSI(JSTRA).GT.NSRFS) CALL EIRENE_MASPRM
     .   ('NSRFS',5,NSRFS,'NSRFSI(I)',9,NSRFSI(I),IERROR)
        DO 715 J=1,NSRFSI(JSTRA)
          READ (IUNIN,6666) INUM,INDIM(J,I),INSOR(J,I),
     .                      INGRDA(J,I,1),INGRDE(J,I,1),
     .                      INGRDA(J,I,2),INGRDE(J,I,2),
     .                      INGRDA(J,I,3),INGRDE(J,I,3)
          READ (IUNIN,6664) SORWGT(J,I),SORLIM(J,I),
     .                      SORIND(J,I),SOREXP(J,I),SORIFL(J,I)
          READ (IUNIN,6666) NRSOR(J,I),NPSOR(J,I),NTSOR(J,I),
     .                      NBSOR(J,I),NASOR(J,I),NISOR(J,I),ISTOR(J,I)
          READ (IUNIN,6664) SORAD1(J,I),SORAD2(J,I),SORAD3(J,I),
     .                      SORAD4(J,I),SORAD5(J,I),SORAD6(J,I)
  715   CONTINUE
C  VELOCITY SPACE DISTRIBUTION
        READ (IUNIN,6664) SORENI(I),SORENE(I),SORVDX(I),SORVDY(I),
     .                    SORVDZ(I)
        READ (IUNIN,6664) SORCOS(I),SORMAX(I),SORCTX(I),SORCTY(I),
     .                    SORCTZ(I),
     .                    RAYFRAC(I)
!  STORE VALUES FOR JSON OUTPUT
        SORCOS_IN(I) = SORCOS(I)
        SORMAX_IN(I) = SORMAX(I)
C
  712 CONTINUE
C
      READ (IUNIN,*)
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
      CALL EIRENE_MASAGE ('*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES')
  810 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 810
      READ (ZEILE,6666) NZADD
      IREAD=0
      WRITE (iunout,*) '       NZADD= ',NZADD
      CALL EIRENE_LEER(1)
      NULLIFY(TEMPLIST)
      NULLIFY(DENLIST)
      NULLIFY(VELLIST)
      NULLIFY(VOLLIST)
      ALLOCATE (INI_ZONE(NZADD))
      ALLOCATE (INE_ZONE(NZADD))
#ifndef NO_JSON
      ALLOCATE (CH3_STACK(NZADD))
#endif
      INI_ZONE = 0
      INE_ZONE = 0
      DO 811 I=1,NZADD
#ifndef NO_JSON
        NULLIFY(CH3_STACK(I)%HEAD)
        NULLIFY(CH3_STACK(I)%LAST)
#endif
CDR  skip reading optional comment lines
  814   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 814
CDR
        READ (ZEILE,6666) INI,INE
        INI_ZONE(I) = INI
        INE_ZONE(I) = INE
        IREAD=0
        IF (INI.GT.NRAD.OR.INI.LE.0) GOTO 998
        IF (INE.GT.NRAD) GOTO 998
        IF (INE.LE.0) INE=INI
  812   READ (IUNIN,'(A72)') ZEILE
C  IGJUM3 FLAG
        IF (ZEILE(1:3).EQ.'CH3') THEN
          INILGJ=MAX(1,MIN(NOPTIM,INI))
          if (inilgj.ne.ini) write (iunout,*)
     .      'CH3 option failed. ini, noptim =',ini,NOPTIM
          INELGJ=MIN(NOPTIM,INE)
          if (inelgj.ne.ine) write (iunout,*)
     .      'CH3 option failed. ine, noptim =',ine,NOPTIM
          DO 821 II=INILGJ,INELGJ
            IF (NLIMPB >= NLIMPS) THEN
              CALL EIRENE_DEKEY (ZEILE(4:72),IGJUM3,0,NOPTIM,II,NLIMPS)
            ELSE
              CALL EIRENE_DEKEYB
     .                    (ZEILE(4:72),IGJUM3,0,NOPTIM,II,NLIMPB,NBITS)
            END IF
  821     CONTINUE
#ifndef NO_JSON
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(ch3_stack(i),zeile(1:72))
          zeileCK=zeile
          call eirene_push_string_stack(ch3_stack(i),zeileCK(1:72))
!cym/cpg end
#endif
          GOTO 812
C  TEMPERATURE
        ELSEIF (ZEILE(1:1).EQ.'T') THEN
          DO 822 II=INI,INE
            ALLOCATE(TEMPCUR)
            READ (ZEILE(2:72),66664) TEMPCUR%IDION,TEMPCUR%TE,
     .                               TEMPCUR%TI
            TEMPCUR%II = II
            TEMPCUR%NEXT => TEMPLIST
            TEMPLIST => TEMPCUR
  822     CONTINUE
          GOTO 812
C  DENSITY
        ELSEIF (ZEILE(1:1).EQ.'D') THEN
          DO 823 II=INI,INE
            ALLOCATE(DENCUR)
            READ (ZEILE(2:72),66664) DENCUR%IDION,DENCUR%DI
            DENCUR%II = II
            DENCUR%NEXT => DENLIST
            DENLIST => DENCUR
  823     CONTINUE
          GOTO 812
C  VELOCITY (CM/SEC OR MACH)
        ELSEIF ((ZEILE(1:1).EQ.'V'.AND.ZEILE(2:2).NE.'L')
     .       .OR.ZEILE(1:1).EQ.'M') THEN
          DO 824 II=INI,INE
            ALLOCATE(VELCUR)
            READ (ZEILE(2:72),66664) VELCUR%IDION,VELCUR%VX,VELCUR%VY,
     .                               VELCUR%VZ
            IF (ZEILE(1:1).EQ.'M') THEN
              VELCUR%IZ=1
            ELSEIF (ZEILE(1:1).EQ.'V') THEN
              VELCUR%IZ=-1
            ENDIF
            VELCUR%II = II
            VELCUR%NEXT => VELLIST
            VELLIST => VELCUR
  824     CONTINUE
          GOTO 812
C  VOLUME
        ELSEIF (ZEILE(1:2).EQ.'VL') THEN
          DO 825 II=INI,INE
            ALLOCATE(VOLCUR)
            READ (ZEILE(3:72),6664) VOLCUR%VOL
            VOLCUR%II = II
            VOLCUR%NEXT => VOLLIST
            VOLLIST => VOLCUR
  825     CONTINUE
          GOTO 812
        ELSEIF (ZEILE(1:1).EQ.'*') THEN
          IREAD=1
          GOTO 811
        ELSE
          GOTO 998
        ENDIF
  811 CONTINUE

C
C  READ DATA FOR STATISTICS AND NON-ANALOG MODEL, 900--999
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 9. DATA FOR STATISTIC AND NON-ANALOG MODEL')
C
  910 IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IREAD = 0
      IF (ZEILE(1:1) .EQ. '*') GOTO 910
C  DATA FOR CONDITIONAL EXPECTATION ESTIMATOR
      NLOGIN = MAX(1,NATMI + NMOLI + NIONI + NPHOTI)
C  read 60 or more logical flags
      ALLOCATE (LOGRDH(NLOGIN))
      LOGRDH = .FALSE.
      DO J=1, NLOGIN, 60
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) LOGRDH(J:MIN(J+59,NLOGIN))
        READ (IUNIN,'(A72)') ZEILE
      END DO
c
      IF (NATMI > 0) NLPRCA(1:NATMI) = LOGRDH(1:NATMI)
      IF (NMOLI > 0)
     .    NLPRCM(1:NMOLI) = LOGRDH(NATMI+1 : NATMI+NMOLI)
      IF (NIONI > 0) NLPRCI(1:NIONI) =
     .       LOGRDH(NATMI+NMOLI+1 : NATMI+NMOLI+NIONI)
      IF (NPHOTI > 0) NLPRCPH(1:NPHOTI) =
     .       LOGRDH(NATMI+NMOLI+NIONI+1 : NLOGIN)
      DEALLOCATE (LOGRDH)
C  READ NON-DEFAULT OR ADDITIONAL SURFACES, TOWARDS WHICH COND. EXP. EST. IS USED
C  DEFAULT: ???
      READ (ZEILE,6666) NPRCSF

      NPRCSF=MIN0(NLIMPS,NPRCSF)
c  read 12 or more integer flags
      IPRCSF=1
  911 CONTINUE
      IF (IPRCSF.LE.NPRCSF) THEN
        READ (IUNIN,6666) (IPRSF(J),J=1,12)
        DO J=1,12
        IF (IPRSF(J).GT.0.AND.IPRSF(J).LE.NLIMPS)
     .     NLPRCS(IPRSF(J))=.TRUE.
        ENDDO
        IPRCSF=IPRCSF+12
        GOTO 911
      ENDIF

C  DATA FOR SPLITTING AND RUSSIAN ROULETTE
      READ (IUNIN,6666) MAXLEV,MAXRAD,MAXPOL,MAXTOR,MAXADD
      MXL=15
      IF (MAXLEV.GT.MXL)
     .    CALL EIRENE_MASPRM('MXL',3,MXL,'MAXLEV',6,MAXLEV,IERROR)
      IF (ABS(MAXRAD).GT.N1ST)
     .    CALL EIRENE_MASPRM('N1ST',4,N1ST,'MAXRAD',6,MAXRAD,IERROR)
      IF (MAXPOL.GT.N2ND)
     .    CALL EIRENE_MASPRM('N2ND',4,N2ND,'MAXPOL',6,MAXPOL,IERROR)
      IF (MAXTOR.GT.N3RD)
     .    CALL EIRENE_MASPRM('N3RD',4,N3RD,'MAXTOR',6,MAXTOR,IERROR)
      IF (MAXADD.GT.NLIM)
     .    CALL EIRENE_MASPRM('NLIM',4,NLIM,'MAXADD',6,MAXADD,IERROR)
      DO I=1,MAXRAD
        READ (IUNIN,66665) ID,NSSPL(I),PRMSPL(I)
      ENDDO
      DO I=1,MAXPOL
        READ (IUNIN,66665) ID,NSSPL(N1ST+I),PRMSPL(N1ST+I)
      ENDDO
      DO I=1,MAXTOR
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+I),PRMSPL(N1ST+N2ND+I)
      ENDDO
      DO I=1,MAXADD
        READ (IUNIN,66665) ID,NSSPL(N1ST+N2ND+N3RD+I),
     .                        PRMSPL(N1ST+N2ND+N3RD+I)
      ENDDO
C  DATA FOR BIAS SAMPLING
      READ (IUNIN,6664) WMINV,WMINS,WMINC,WMINL
      WMINV=MAX(WMINV,EPS60)
      WMINS=MAX(WMINS,EPS60)
      WMINC=MAX(WMINC,EPS60)
      READ (IUNIN,6664) SPLPAR
C  DATA FOR STANDARD DEVIATION
      READ (IUNIN,*)
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_SPC
      WRITE (iunout,*) '       NSIGVI,NSIGSI,NSIGCI = ',
     .                         NSIGVI,NSIGSI,NSIGCI
      WRITE (iunout,*) '       NSIGI_SPC            = ',NSIGI_SPC
      CALL EIRENE_LEER(1)
      DO 913 J=1,NSIGVI
        READ (IUNIN,6666) IGH(J),IIH(J)
  913 CONTINUE
      DO 914 J=1,NSIGSI
        READ (IUNIN,6666) IGHW(J),IIHW(J)
  914 CONTINUE
      DO 915 J=1,NSIGCI
        READ (IUNIN,6666) IGHC(1,J),IIHC(1,J),IGHC(2,J),IIHC(2,J)
  915 CONTINUE
C
C   READ DATA FOR ADDITIONAL AND SURFACE-AVERAGED TALLIES
      READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION')
      CALL EIRENE_MASAGE
     .  ('        ESTIMATORS AND ALGEBRAIC EXPRESSIONS')
C
 1010 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1010
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       NADVI,NCLVI,NALVI    = ',
     .                         NADVI,NCLVI,NALVI
      WRITE (iunout,*) '       NADSI,NALSI,NADSPC   = ',
     .                         NADSI,NALSI,NADSPC
      CALL EIRENE_LEER(1)
C
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10A. DATA FOR ADDITIONAL TALLIES')

      DO 1020 J=1,NADVI
 1021   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1021
        READ (ZEILE,6666) IADVE(J),IADVS(J),IADVT(J),IADRC(J)
        READ (IUNIN,'(A72)') TXTTLA(J)
        READ (IUNIN,'(2A24)') TXTSCA(J),TXTUTA(J)
 1020 CONTINUE

      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10B. DATA FOR COLLISION ESTIMATORS')
      DO 1030 J=1,NCLVI
 1031   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1031
        READ (ZEILE,6666) ICLVE(J),ICLVS(J),ICLVT(J),ICLRC(J)
        READ (IUNIN,'(A72)') TXTTLC(J)
        READ (IUNIN,'(2A24)') TXTSCC(J),TXTUTC(J)
 1030 CONTINUE
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS')
      DO 1040 J=1,NALVI
 1041   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1041
        READ (ZEILE,'(A72)') CHRTAL(J)
        READ (IUNIN,'(A72)') TXTTLR(J)
        READ (IUNIN,'(2A24)') TXTSCR(J),TXTUTR(J)
 1040 CONTINUE
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES')
      DO 1050 J=1,NADSI
 1051   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1051
        READ (ZEILE,6666) IADSE(J),IADSS(J),IADST(J),IADSC(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSA)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSA),TXTUNW(J,NTLSA)
 1050 CONTINUE
C
      READ (IUNIN,*)
      CALL EIRENE_MASAGE('*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES')
      DO 1060 J=1,NALSI
 1061   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1061
        READ (ZEILE,'(A72)') CHRTLS(J)
        READ (IUNIN,'(A72)') TXTTLW(J,NTLSR)
        READ (IUNIN,'(2A24)') TXTSPW(J,NTLSR),TXTUNW(J,NTLSR)
 1060 CONTINUE
C
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) == '** ') THEN
        CALL EIRENE_MASAGE('*** 10F. DATA FOR SPECTRA')
        IF (NADSPC > 0) THEN
          ALLOCATE(ESTIML(NADSPC))
          IF (NSMSTRA > 0) ALLOCATE(SMESTL(NADSPC))
        ELSE
          ALLOCATE(ESTIML(1))
        END IF
C
        DO J=1,NADSPC
          DO
            READ (IUNIN,'(A72)') ZEILE
            IF (ZEILE(1:1) .NE. '*') EXIT
          END DO
cdr  this input and this code segment should be modified:
c  first read isrfcll: "type of spectrum",
c  then next: ispsrf "where is this spectrum scored"
c  then further details (species: iptyp,ipspz, and then: direction,....)
c
          call fix_integer_input(zeile,8)
          READ (ZEILE,'(12I6)') ISPSRF, IPTYP, IPSPZ, ISPTYP, NSPS,
     .                          ISRFCLL, IDIREC, ISPOPT
          READ (IUNIN,'(6E12.4)') SPCMN, SPCMX, SPC_SHIFT,
     .                            SPCPLT_X, SPCPLT_Y, SPCPLT_SAME
          SPCVX = 0._DP
          SPCVY = 0._DP
          SPCVZ = 0._DP
          IF (IDIREC /= 0) THEN
            READ (IUNIN,'(6E12.4)') SPCVX, SPCVY, SPCVZ
            VNORM = SQRT(SPCVX**2+SPCVY**2+SPCVZ**2)+EPS60
            SPCVX = SPCVX / VNORM
            SPCVY = SPCVY / VNORM
            SPCVZ = SPCVZ / VNORM
cdr  surface tallies: no directional spectral tally scoring yet implemented
            IF (ISRFCLL == 0) THEN
              WRITE (IUNOUT,*) ' SPECTRUM NUMBER ',J
              WRITE (IUNOUT,*) ' DEFINITION OF LINE OF SIGHT FOR',
     .                         ' SURFACE SPECTRUM IS NOT FORESEEN'
              WRITE (IUNOUT,*) ' SPECIFICATION OF DIRECTION IS',
     .                         ' IGNORED'
              IDIREC = 0
            END IF
          END IF
          IF (ISPOPT == 2) THEN
            READ (IUNIN,'(12I6)') ISPLDEG
          END IF

c   isrcfll=0 : surface-averaged tally
c   isrfcll=1 : volume-averaged tally, integrated over all directions
c   isrfcll=2 : volume-averaged tally, along a specific direction
c  it seems: surface-averaged directional tallies: not yet foreseen
c  one may try to score coefficients of orthogonal angular expansion, per energy bin.

cdr  better: first discriminate by isrfcll, then, for each value of isrfcll: do the rest
          IF (ISPSRF > 0) THEN
            IF ((ISRFCLL == 0) .AND. (ISPSRF > NLIMI)) THEN
C  SPECTRUM AT AN ADDITIONAL SURFACE
              WRITE (iunout,*)
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER  = ',ISPSRF
              IERROR = IERROR + 1
            END IF
C  SPECTRUM IN CELL, POSSIBLY ALONG A CERTAIN DIRECTION
cdr:  next 2 lines: why different condition for cell-based spectra cell numbers ??
! cell-based spectrum in scoring cell ISPSRF
            IF (((ISRFCLL == 1) .AND. (ISPSRF > NRTAL)) .OR.
! directional cell-based spectrum in scoring cell ISPSRF
     .          ((ISRFCLL == 2) .AND. (ISPSRF > NRAD))) THEN
              WRITE (iunout,*)
     .          ' CELL INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' CELL NUMBER     = ',ISPSRF
              IERROR = IERROR + 1
            END IF
          ELSEIF (ISPSRF < 0) THEN
            ISPSRF = ABS(ISPSRF)
            IF ((ISRFCLL == 0) .AND. (ISPSRF > NSTSI)) THEN
C  SPECTRUM AT A NON-DEFAULT STANDARD SURFACE
              WRITE (iunout,*)
     .          ' SURFACE INDEX FOR SPECTRUM OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SURFACE NUMBER  = ',ISPSRF
              IERROR = IERROR + 1
            ELSEIF (ISRFCLL == 0) THEN
              ISPSRF = NLIM+ISPSRF
            ELSE
C           IF (((ISRFCLL == 1) .AND. (ISPSRF < 0)) .OR.
C    .          ((ISRFCLL == 2) .AND. (ISPSRF < 0))) THEN
C    .      ... WRONG INPUT !
            END IF
          ELSEIF (ISPSRF == 0) THEN
            WRITE (iunout,*) ' SURFACE OR CELL INDEX = 0: NOT FORESEEN'
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          END IF

          IF ((IPTYP < 0) .OR. (IPTYP > 4)) THEN
            WRITE (iunout,*) ' PARTICLE TYPE ',IPTYP,' NOT FORESEEN'
            WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
            IERROR = IERROR + 1
          ELSE
            IF (((IPTYP == 0)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NPHOTI))) .OR.
     .          ((IPTYP == 1)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NATMI))) .OR.
     .          ((IPTYP == 2)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NMOLI))) .OR.
     .          ((IPTYP == 3)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NIONI))) .OR.
     .          ((IPTYP == 4)
     .        .AND.((IPSPZ < 0).OR.(IPSPZ > NPLSI)))) THEN
              WRITE (iunout,*) ' PARTICLE SPECIES INDEX OUT OF BOUNDS'
              WRITE (iunout,*) ' SPECTRUM NUMBER = ',J
              WRITE (iunout,*) ' SPECIES NUMBER  = ',IPSPZ
              IERROR = IERROR + 1
            END IF
          END IF

          IF ((ISRFCLL == 0) .AND. (ISPTYP > 2)) THEN
            WRITE (IUNOUT,*) ' WRONG TYPE OF SPECTRUM SPECIFIED'
            WRITE (IUNOUT,*) ' SPECTRUM NUMBER = ',J
            WRITE (IUNOUT,*) ' SPECTRUM TYPE   = ',ISPTYP
          END IF

          ESPEC => ESTIML(J)

          ESPEC%ISPCSRF = ISPSRF
          ESPEC%IPRTYP = IPTYP
          ESPEC%IPRSP = IPSPZ
          ESPEC%ISPCTYP = ISPTYP
          ESPEC%NSPC = ABS(NSPS)
          ESPEC%ISRFCLL = ISRFCLL
          ESPEC%IDIREC = IDIREC
          ESPEC%ISPCOPT = ISPOPT
          ESPEC%ISPLDEG = ISPLDEG
cdr
cdr       ESPEC%LOG = .FALSE. ! this was too restrictive !
cdr  Option     LOG = .TRUE. WAS ALREADY AVAILABLE IN SCORING/UPDATE_SPECTRUM

cdr   X.B. correction Sept 17, from SOLPS-ITER branch,
          IF (NSPS > 0) THEN
            NSPSA=NSPS
            ESPEC%LOG = .FALSE.
            ESPEC%SPCMIN=SPCMN
            ESPEC%SPCMAX=SPCMX
            ESPEC%SPCDEL=(SPCMX-SPCMN)/REAL(NSPSA,DP)
          ELSEIF (NSPS < 0) THEN
            NSPSA=-NSPS
            ESPEC%LOG = .TRUE.
cdr  missing here: check for SPCMN, SPCMX > 0.
cdr  In particular in case IDIREC > 0 scores may fall into negative bins.
cdr  tbd: Rule out combination of IDIREC > 0 and NSPS < 0
            ESPEC%SPCMIN = log10(SPCMN)
            ESPEC%SPCMAX = log10(SPCMX)
            ESPEC%SPCDEL = log10(SPCMX/SPCMN)/REAL(NSPSA,DP)
          ELSE
            WRITE (iunout,*) ' SPECTRUM TALLY NUMBER = ',J
            WRITE (iunout,*) ' ZERO NO. OF SPECTRAL BINS "NSPS". EXIT '
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          ESPEC%ESP_00=SPC_SHIFT
          ESPEC%SPC_XPLT=SPCPLT_X
          ESPEC%SPC_YPLT=SPCPLT_Y
          ESPEC%SPC_SAME=SPCPLT_SAME
          ESPEC%SPCVX=SPCVX
          ESPEC%SPCVY=SPCVY
          ESPEC%SPCVZ=SPCVZ
C  MIN AND MAX ENERGY SCORE (EV) ON THIS TALLY
          ESPEC%ESP_MIN=1.E30
          ESPEC%ESP_MAX=-1.E30

          ESPEC%SPCDELI=1._DP/(ESPEC%SPCDEL+EPS60)
          ALLOCATE(ESPEC%SPC(0:NSPSA+1))
c  standard deviation of spectrally resolved tallies
          IF (NSIGI_SPC > 0) THEN
            ALLOCATE(ESPEC%SDV(0:NSPSA+1))
            ALLOCATE(ESPEC%SGM(0:NSPSA+1))
            ALLOCATE(ESPEC%STV(0:NSPSA+1))
            ALLOCATE(ESPEC%GG(0:NSPSA+1))
          ELSE
            NULLIFY(ESPEC%SDV, ESPEC%SGM, ESPEC%STV, ESPEC%GG)
          END IF
          ESPEC%SPC(0:NSPSA+1) = 0._DP
          IF (ISPOPT==2) THEN
            ALLOCATE(ESPEC%SPCAN(0:ISPLDEG,0:NSPSA+1))
            ESPEC%SPCAN(0:ISPLDEG,0:NSPSA+1) = 0._DP
          END IF
          ESPEC%IMETSP = 0

          IF (NSMSTRA > 0) THEN
c  sum over strata
            SSPEC => SMESTL(J)
            ALLOCATE(SSPEC%SPC(0:NSPSA+1))
c  standard deviation of spectra tallies, sum over strata intermediate storage
            IF (NSIGI_SPC > 0) THEN
              ALLOCATE(SSPEC%SDV(0:NSPSA+1))
              ALLOCATE(SSPEC%SGM(0:NSPSA+1))
              ALLOCATE(SSPEC%STV(0:NSPSA+1))
              ALLOCATE(SSPEC%GG(0:NSPSA+1))
            ELSE
              NULLIFY(SSPEC%SDV, SSPEC%SGM, SSPEC%STV, SSPEC%GG)
            END IF
            SMESTL(J) = ESTIML(J)
          END IF
C
        END DO
        IREAD=0
      ELSE
Cdr  No block 10F found, i.e. no spectrally resolved tallies at all.
cdr  Why do we allocate estiml in input.f and not in eirmod_cestim ?
        IF (.NOT.ALLOCATED(ESTIML)) THEN
          ALLOCATE(ESTIML(1))
        END IF
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      IF (IREAD == 0) READ (IUNIN,*)
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT')
c
c  search for input block 11a
 1110 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1110
* For gfortran: it does not accept empty field for logical
      call fix_logical_input(zeile,35)
      READ (ZEILE,6665) TRCPLT,TRCHST,TRCNAL,TRCMOD,TRCSIG,
     .                  TRCGRD,TRCSUR,TRCREF,TRCFLE,TRCAMD,
     .                  TRCINT,TRCLST,TRCSOU,TRCREC,TRCTIM,
     .                  TRCBLA,TRCBLM,TRCBLI,TRCBLP,TRCBLE,
     .                  TRCBLPH,TRCTAL,TRCOCT,TRCCEN,TRCRNF,
CVK TRACING FOR DEBUGGING, V.Kotov: not in use in present EIRENE version
cdr  .                  TRCDBG2,TRCDBGE,TRCDBGM,TRCDBGF,TRCDBGL,
cdr  .                  TRCDBGS,TRCDBGG,TRCDBGMPI,TRCDBGC,
CPB  ACTIVATE SPECIES-RESOLVED CPU CONSUMPTION OPTION
     .                  TRCHKTIM,TRCSCL
      do I = 0, NSTRA, 60
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) (TRCSRC(J),J=I,MIN(I+59,NSTRA))
      end do
C
      READ (IUNIN,6666) NVOLPR, NSPCPR
      CALL EIRENE_LEER(1)
      IF(NVOLPR > NVLPR) THEN
         WRITE (iunout,*) 'NVOLPR > NVLPR, EXIT.'
         WRITE (iunout,*) ' NVOLPR, NVLPR ', NVOLPR, NVLPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
      WRITE (iunout,*) '       NVOLPR= ',NVOLPR
      DO 1120 J=1,NVOLPR
        READ (IUNIN,6666) NTLV,NFLGV,NSPZV1,NSPZV2,NTLVF
        IF (NTLV.LT.-NTALI.OR.NTLV.GT.NTALV) GOTO 990
cdr  output stream for particular tallies: allow only ntlv>=70
        IF ((NTLVF.NE.0).AND.(NTLVF.LT.70)) GOTO 995
        NPRTLV(J)  =NTLV
        NFLAGV(J)  =NFLGV
        NSPEZV(J,1)=NSPZV1
        NSPEZV(J,2)=NSPZV2
        NTLVFL(J)  =NTLVF
 1120 CONTINUE
C
      READ (IUNIN,6666) NSURPR
      IF(NSURPR > NSRPR) THEN
         WRITE (iunout,*) 'NSURPR > NSRPR, EXIT.'
         WRITE (iunout,*) ' NSURPR, NSRPR ', NSURPR, NSRPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      WRITE (iunout,*) '       NSURPR= ',NSURPR
      DO 1130 J=1,NSURPR
        READ (IUNIN,6666) NSRF,NTLS,NFLGS,NSPZS1,NSPZS2,NTLSF
c  negative surface numbers may be used for non-default standard surfaces
        IF (NSRF.LT.0) NSRF=NLIM+IABS(NSRF)
c  zero as surface number for time horizon
        IF (NSRF.EQ.0.AND.NLIM+NSTSI.LT.NLIMPS) NSRF=NLIM+NSTSI+1

        IF (NSRF.LE.0.OR.NSRF.GT.NLIMPS) GOTO 991
        NPRSRF(J)=NSRF
        NPRTLS(J)=NTLS
        NFLAGS(J)=NFLGS
        NSPEZS(J,1)=NSPZS1
        NSPEZS(J,2)=NSPZS2
        NTLSFL(J)=NTLSF
 1130 CONTINUE

c  overrule default switching on/off of volume-averaged tallies

 1131 READ (IUNIN,'(A72)') ZEILE
      CALL EIRENE_UPPERCASE(ZEILE)
      IREAD=1
      IF ((ZEILE(1:1) .NE.'*') .AND. (SCAN(ZEILE,'FT') == 0)) THEN

c  there are ntlvout volume tallies to be dealt with
c  and ntlsout surface tallies

        READ (ZEILE,6666) NTLVOUT
        IREAD=0
        ITLV=0
        IF (NTLVOUT > 0) THEN
          ALLOCATE (ITLVOUT(NTLVOUT))
          ITLVOUT = 0
        END IF
        DO WHILE ((NTLVOUT > 0) .AND. (ITLV < NTLVOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALV)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALV(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALV(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
              ITLVOUT(ITLV+J) = NUMTAL(J)
            END IF
          END DO
          ITLV=ITLV+12
        END DO

c  overrule default switching on/off of surface-averaged tallies

        READ (IUNIN,6666) NTLSOUT
        ITLS=0
        IF (NTLSOUT > 0) THEN
          ALLOCATE (ITLSOUT(NTLSOUT))
          ITLSOUT = 0
        END IF
        DO WHILE ((NTLSOUT > 0) .AND. (ITLS < NTLSOUT))
          READ (IUNIN,6666) (NUMTAL(J),J=1,12)
          DO J=1,12
            IF ((NUMTAL(J) /= 0) .AND. (ABS(NUMTAL(J)) <= NTALS)) THEN
              IF (NUMTAL(J) < 0) THEN
                LMISTALS(ABS(NUMTAL(J))) = .TRUE.  ! SWITCHED OFF
              ELSE
                LMISTALS(NUMTAL(J)) = .FALSE.      ! SWITCHED ON
              END IF
              ITLSOUT(ITLS+J) = NUMTAL(J)
            END IF
          END DO
          ITLS=ITLS+12
        END DO

      ELSE
C EITHER A COMMENT LINE, OR THE NEXT INPUT CARD (LOGICALS FOR PLOTTING) IS ON 'zeile'
        IF (ZEILE(1:1) .eq.'*')  goto 1131
c  reading 'switch tallies off' done
      END IF

C  search for input block 11b

 1132 IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:1) .EQ. '*') GOTO 1132
C  2D GEOMETRY PLOT
      call fix_logical_input(zeile,16)
      READ (ZEILE,6665) PL1ST,PL2ND,PL3RD,PLADD,PLHST,
     .                  PLCUT(1),PLCUT(2),PLCUT(3),PLBOX,PLSTOR,
     .                  PLNUMV,PLNUMS,PLARR,LRPSCUT,PLIDL,
     .                  PLVTK
      READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                  NPLINT,NPLOTT,NPLDLT

#ifndef NO_JSON
      CALL EIRENE_READ_BLOCK_11_USR
#else
      CALL EIRENE_READ_BLOCK_11_USR_FIXED
#endif
!PB moved to subroutine EIRENE_READ_BLOCK_11_USR
!pb  .                  NPLINT,NPLOTT,NPLDLT,nsrfcls !csw 22jul2011 added nsrfcls
C  3D GEOMETRY PLOT
      DO 1140 J=1,5
        READ (IUNIN,6662) PL3A(J),TEXTLA(J),IPLTA(J),
     .                (IPLAA(J,I),IPLEA(J,I),I=1,IPLTA(J))
 1140 CONTINUE
      DO 1141 J=1,3
        READ (IUNIN,6662) PL3S(J),TEXTLS(J),IPLTS(J),
     .                (IPLAS(J,I),IPLES(J,I),I=1,IPLTS(J))
 1141 CONTINUE
C
      READ (IUNIN,6664) CH2MX,CH2MY,      CH2X0,CH2Y0,CH2Z0
      READ (IUNIN,6664) CH3MX,CH3MY,CH3MZ,CH3X0,CH3Y0,CH3Z0
      READ (IUNIN,6664) ANGLE1,ANGLE2,ANGLE3
C
C  PARTICLE HISTORY PLOTS IN 2D OR 3D GEOMETRY PLOTS
      READ (IUNIN,6666) I1TRC,I2TRC,(ISYPLT(J),J=1,8),ILINIE
C

c  search for input block 11c

 1151 READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:1) .EQ. '*') GOTO 1151
C  DATA FOR PLOTS OF VOLUME-AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
      WRITE (iunout,*) '       NVOLPL= ',NVOLPL
      CALL EIRENE_LEER(1)
C
      IF (NVOLPL.GT.NPTAL)
     .  CALL EIRENE_MASPRM('NPTAL',5,NPTAL,'NVOLPT',6,NVOLPL,IERROR)
C
      IF (NVOLPL.LE.0) GOTO 1175
      DO I = 0, NSTRA, 60
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,60)
        READ (ZEILE,6665) (PLTSRC(J),J=I,MIN(I+59,NSTRA))
      END DO
C
      IF (LRPSCUT) THEN
        READ (IUNIN,6664) CUTPLANE(1:4)
!  TEST IF CUTPLANE IS AXIS-PARALLEL
        IF ((SUM(ABS(CUTPLANE(2:4)))-1._DP > EPS10) .OR.
     .      (CUTPLANE(2)*CUTPLANE(3) > 0._DP) .OR.
     .      (CUTPLANE(2)*CUTPLANE(4) > 0._DP) .OR.
     .      (CUTPLANE(3)*CUTPLANE(4) > 0._DP)) THEN
          WRITE (iunout,*) ' PLANE IS NOT PARALLEL TO ONE AXIS'
          WRITE (iunout,*) ' CUTTING OF TETRAHEDRA ABANDONED'
          WRITE (iunout,*) ' LRPSCUT SET TO .FALSE.'
          LRPSCUT = .FALSE.
        END IF
      END IF
C
      IPLANE=1
      LRAPS3D=.FALSE.
      LR3DCON=.FALSE.
      RAPSDEL=-HUGE(1._DP)
C
      DO 1150 J=1,NVOLPL
 1152   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1152
        READ (ZEILE,6666) NSP
        IF (NSP.GT.NPLT)
     .    CALL EIRENE_MASPRM('NPLT',4,NPLT,'NSP',3,NSP,IERROR)
        NSPTAL(J)=NSP
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,4)
        READ (ZEILE,6665) PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
        READ (IUNIN,6664) TALZMI(J),TALZMA(J),TALXMI(J),TALXMA(J),
     .                    TALYMI(J),TALYMA(J)
        IF (PLTL2D(J)) THEN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,2)
          READ (ZEILE,6665) LHIST2(J),LSMOT2(J)
          DO 1160 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,
     .                        NPLIN2(J,I),NPLOT2(J,I),NPLDL2(J,I)
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
            NPLIN2(J,I) = MAX(NPLIN2(J,I),1)
            NPLDL2(J,I) = MAX(NPLDL2(J,I),1)
 1160     CONTINUE
        ENDIF
        IF (PLTL3D(J)) THEN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,8)
          READ (ZEILE,6665) LHIST3(J),LCNTR3(J),LSMOT3(J),
     .                      LRAPS3(J),LVECT3(J),LRPVC3(J),
     .                      LRPS3D,LRPSCN
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,3)
          READ (ZEILE,6665) LPRAD3(J),LPPOL3(J),LPTOR3(J)
          DO 1161 I=1,NSPTAL(J)
            READ (IUNIN,6666) ISPTAL(J,I),NTL,IPROJ3(J,I),
     .                        NPLI13(J,I),NPLO13(J,I),
     .                        NPLI23(J,I),NPLO23(J,I),IPLN
            IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
            NPTALI(J,I)=NTL
 1161     CONTINUE
          READ (IUNIN,6664) TALW1(J),TALW2(J),FCABS1(J),FCABS2(J),
     .                      RPSDL
          IF (FCABS1(J).LE.0.0D0) FCABS1(J)=1.0D0
          IF (FCABS2(J).LE.0.0D0) FCABS2(J)=1.0D0
          LRAPS3D=LRAPS3D.OR.LRPS3D
          LR3DCON=LR3DCON.OR.LRPSCN
          IPLANE=MAX(IPLANE,IPLN)
          RAPSDEL=MAX(RAPSDEL,RPSDL)
        ENDIF
 1150 CONTINUE
      IF (NLTRA) RAPSDEL=RAPSDEL*DEGRAD
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
 1175 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:3) .EQ. '***') GOTO 1205
      GOTO 1175
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
 1205 CALL EIRENE_MASAGE('*** 12. DATA FOR DIAGNOSTIC MODULE')
 1210 READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1) .EQ. '*') GOTO 1210

c June 18: new: more general option for definition of emission profiles
cdr        in ADDV tally. Can e.g. be used in line of sight NCHTAL=2 option,
cdr        or for general additional output.
cdr  Based on population coefficients, QSS ratios, etc already read in block 4.

cdr  read further atomic/molecular data: population coefficients, QSS ratios, etc
cdr       needed for setting up volumetric line emissivity profiles
cdr       as further add. tally ADDV (additional to those already defined in block 10a)
cdr       The ADDV tallies may then be used for line of sight integration, along
cdr       the CHORDS defined further below.
cdr  Distinct from the ADDV tallies from block 10a these
cdr  further ADDV tallies (beyond NADVI) have their text
cdr  (tally label, units, species) assigned in subr. emissivity.f.

cdr  check if such emissivity tallies are defined:
cdr  Search for 'DEFINE LINES' in next card
      CALL EIRENE_MASAGE
     .  ('*** 12A. DATA FOR EMISSIVITY PROFILES (ADDV)  ')
      ULINE=ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
cdr  careful: slightly different from find_param
      NLEMIS = NCHOR > 0  !dr  and: NCHTAL=2 ??
      LDEF_LINES = .FALSE.

      IF (INDEX(ULINE,'DEFINE_LINES') > 0) THEN
cdr read volumetric emission profile data
        IADV = NADVI
        NLEMIS = .TRUE.
        READ (IUNIN,6666) NUM_LINES, MOD_ADDV
      END IF

!   Option 1: Non-zero number of lines are defined
      IF ((INDEX(ULINE,'DEFINE_LINES') > 0).and.(NUM_LINES > 0)) THEN
        LDEF_LINES = .TRUE.

        ALLOCATE (EMIS_LINES(NUM_LINES))
        EMIS_LINES%LINE_NAME = REPEAT(' ',80)
        EMIS_LINES%NUM_COMPO = 0
        WRITE (IUNOUT,*) 'FURTHER EMISSIVITY PROFILES SET:'

        DO ILINE=1, NUM_LINES
          READ (IUNIN,'(A80)') ZEILE
          DO WHILE (ZEILE(1:1) == '*')
            READ (IUNIN,'(A80)') ZEILE
            write (IUNOUT,'(A80)') ZEILE
          END DO
          IREAD=1
cdr  further below, a particular emission profile is identified
cdr  (e.g. for a chord ichori) either by its name  (and ch_line%...)
cdr  or, if that fails, by its energy (and EMIN1 flag)
          READ (ZEILE,'(A80)') EMIS_LINES(ILINE)%LINE_NAME
          write (IUNOUT,*) 'EMISSIVITY NAME ',trim(ZEILE)
          IREAD=0
          READ (IUNIN,6666) NUM_COMPO
          READ (IUNIN,6664) EMIS_LINES(ILINE)%EINSTEIN,
     .                      EMIS_LINES(ILINE)%TRANS_EN
C
          EMIS_LINES(ILINE)%NUM_COMPO = NUM_COMPO

c  Deal with ADDV storage for sum over components.
c  The storage on ADDV for individual components is done below (JCOMP loop).
cdr       IADV = IADV + 1
c  if mod_addv=0: reset storage needs for each emissivity profile
c                 back to nadvi+1..nadvi+num_compo+1,
c                 i.e. the ADDV tallies are only saved for one emissivity at a time.
          IF (MOD_ADDV == 0) IADV = NADVI

          EMIS_LINES(ILINE)%IADV_TOTAL = IADV + NUM_COMPO+1

          IF (NUM_COMPO > 0) THEN
            ALLOCATE (EMIS_LINES(ILINE)%COMPO(NUM_COMPO))

            DO JCOMP=1,NUM_COMPO
              READ (IUNIN,'(A72)')
     .              EMIS_LINES(ILINE)%COMPO(JCOMP)%COMPO_NAME

              READ (IUNIN,6666) NUM_CONTRIB, IRC
              EMIS_LINES(ILINE)%COMPO(JCOMP)%NUM_CONTRIB = NUM_CONTRIB
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IRC = IRC
              ALLOCATE
     .         (EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(NUM_CONTRIB))
cdr  for each new component: store emission profile on addv(iadv)
              IADV = IADV + 1
              if (iadv.gt.nadv) then
                write (iunout,*) 'storage conflict on ADDV tally'
                write (iunout,*) 'iadv,nadv ',iadv,nadv
                call eirene_exit_own(1)
              endif
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IADV = IADV

cdr  now we dwell on the contributions to JCOMP:
cdr  e.g. different isotopes,... but same emission rates, same population factors
cdr  within each component JCOMP
              DO KCONTR = 1, NUM_CONTRIB
                CNT%ISP = -1
                CNT%ITP = -1
                CNT%IRATIO = 0
                CNT%IRC_RAT = 0
                CNT%ISP_RAT = -1
                CNT%ITP_RAT = -1
                READ (IUNIN,6666) CNT%ISP, CNT%ITP, CNT%IRATIO

cdr  read QSS ratio between two densities, e.g.: H2+/H2, if nfoli(H2+)=-1.
cdr  If density n_B of parent state (e.g. H2+) for upper level is not amongst the densities
cdr  known in this run,
cdr  but is in an (QSS) equilibrium (nfol<0) with such a density n_A (e.g. H2) instead.
                IF (CNT%IRATIO > 0) THEN
cdr  read QSS density ratio n_B/n_A(Te,ne). Then the
cdr  upper state population is n_B *pop_B(upper)= n_A * ratio * pop_B(upper)
                  READ (IUNIN,6666) CNT%IRC_RAT(1)
cdr   Read a second density ratio.
cdr   It may turn out that, after reading the first density ratio,
cdr   that now the new parent density n_A is still not amongst the density
cdr   known to eirene in this run. Then read a second density ratio.
cdr   Example: n_B=n_H3+, ratio1=nH3+/nH2+. I.e. n_A =nH2+.
cdr   See routine emissivity.f for further explanations.
                  IF (CNT%IRATIO == 2) THEN
                    READ (IUNIN,6666) CNT%IRC_RAT(2),
     .                                CNT%ISP_RAT(1),CNT%ITP_RAT(1),
     .                                CNT%ISP_RAT(2),CNT%ITP_RAT(2)
                    END IF
                  END IF

                EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(KCONTR) = CNT
              END DO  !  kcontr

            END DO    ! jcomp
c  increase counter iadv, for next line, because sum over comp.
c                         is stored on num_compo+1
            IADV=IADV+1
          END IF      ! if jcomp.gt.0
        END DO        ! iline
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
          READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1

!   Option 2: Requesting default lines, there are two ways:
!       - DEFAULT_LINES
!       - DEFINE_LINES && NUM_LINES=0  (backward compatibility)
      ELSEIF ((INDEX(ULINE,'DEFAULT_LINES') > 0)
     .  .or.(INDEX(ULINE,'DEFINE_LINES') > 0).and.(NUM_LINES <= 0)) THEN
        IREAD=0
        NLEMIS=.TRUE.
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
          READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1
      ENDIF

! no definition of emissivity lines was read in
! define OLD default emissivity model for chords, for backward compatibility.
cdr  allocate storage and fill structure EMIS_LINES
cdr  such that old default options are recovered.
cdr This is exclusive: as soon as at least one line emission profile is
cdr read from block "12.0", no default emissivities are set at all.
      IF (NLEMIS.AND..NOT.ALLOCATED(EMIS_LINES))
     .   CALL EIRENE_SETUP_DEFAULT_EMISSIVITY

c  June 18: new: end of new code, further modifications below, for NCHTAL=2 option
      CALL EIRENE_MASAGE
     .  ('*** 12B. DATA FOR LINES OF SIGHT              ')
      IF (IREAD == 0) READ (IUNIN,'(A80)') ZEILE
      READ (ZEILE,6666) NCHORI,NCHENI
      WRITE (iunout,*) '       NCHORI= ',NCHORI
      IREAD=0
      NCHOR = NCHORI
      NCHEN = NCHENI
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(1x,a,2i5)') 'NCHORI,NCHENI= ',NCHORI,NCHENI
      CALL EIRENE_LEER(1)
      IF (IABS(NCHENI).GT.NCHEN)
     .    CALL EIRENE_MASPRM
     .         ('NCHEN',5,NCHEN,'IABS(NCHENI)',12,IABS(NCHENI),
     .                 IERROR)
      IF (NCHORI.LE.0) GOTO 1230
      CALL EIRENE_ALLOC_COMSIG
      CALL EIRENE_INIT_COMSIG
      DO 1220 ICHORI=1,NCHORI
        READ (IUNIN,'(A72)') TXTSIG(ICHORI)
        READ (IUNIN,6666) NCHTAL(ICHORI),NSPSCL(ICHORI),NSPNEW(ICHORI),
     .                    ISTCHR
cdr  new input option: generalized side on line emissivities
cdr                    for nchtal=2 option.
        READ (IUNIN,'(A400)') ZEILE
        IREAD=1

        IF (NCHTAL(ICHORI) == 2) THEN
cdr  old default: neither "define_lines" nor "default_lines" found.
cdr  still: chords with NCHTAL=2 present. Set default now, latest.
          IF (.NOT.ALLOCATED(EMIS_LINES))
     .      CALL EIRENE_SETUP_DEFAULT_EMISSIVITY
cdr  for nchtal=2: one optional extra input card may be read: search for 'USE_LINE'
cdr  and fill CH_LINE_NAME(ICHORI) with the name of that line.
cdr  Alternatively the energy parameters EMIN1 may be used.
          ULINE = ZEILE
          CALL EIRENE_UPPERCASE(ULINE)
          IND = INDEX(ULINE,'USE_LINE')
          IF (IND > 0) THEN
            READ (ZEILE(IND+8:),'(A80)') CH_LINE_NAME(ICHORI)
            READ (IUNIN,'(A400)') ZEILE
           IREAD=1
          END IF
        END IF

! here should come: NSPTP(..), TYPE
        READ (ZEILE,6666) NSPSTR(ICHORI),NSPSPZ(ICHORI),
     .                    NSPINI(ICHORI),NSPEND(ICHORI),
     .                    NSPBLC(ICHORI),NSPADD(ICHORI)
        IREAD=0
        READ (IUNIN,6664) EMIN1(ICHORI),EMAX1(ICHORI),ESHIFT(ICHORI)
        READ (IUNIN,66664) IPIVOT(ICHORI),
     .                     XPIVOT(ICHORI),YPIVOT(ICHORI),ZPIVOT(ICHORI)
        READ (IUNIN,66664) ICHORD(ICHORI),
     .                     XCHORD(ICHORI),YCHORD(ICHORI),ZCHORD(ICHORI)
! automatically add directional cell-based spectra, along line of sight
        NLSTCHR(ICHORI) = ISTCHR > 0
 1220 CONTINUE   ! ichori

      READ (IUNIN,'(A72)') ZEILE
      call fix_logical_input(zeile,5)
      READ (ZEILE,6665) PLCHOR,PLSPEC,PRSPEC,PLARGL,PRARGL

 1230 CONTINUE  ! nchori > 0

C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      IREAD=0
      IF (ZEILE(1:3).NE.'***') GOTO 1230
      IREAD=1
C
C  READ DATA FOR TIME-DEPENDENT AND NONLINEAR MODE  1300--1399
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      IREAD=0
      CALL EIRENE_MASAGE
     .  ('*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION')
C
      NINITL_READ = 0
      READ (IUNIN,6666) NPRNLI_IN, NINITL_READ, NPRMUL
      NPRNLI = NPRNLI_IN
      IF (NPRMUL > 1) NPRNLI=NPRNLI * NPRMUL

      CALL EIRENE_LEER(1)

cdr NLERG: try to use ergodic property to infer cell volumes.
cdr        All absorption processes are turned off, so we need
cdr        a time horizon. Strictly nprnli=1 would be enough,
cdr        but we use a larger default value to estimate variances as well
      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        NPRNLI=100
        IF (NTIME.EQ.0) NTIME=1
        WRITE (iunout,*) '       NPRNLI= ',NPRNLI,
     .                   ' (MODIFIED DUE TO NLERG)'
      ELSE
        WRITE (iunout,*) '       NPRNLI= ',NPRNLI
      ENDIF
      CALL EIRENE_LEER(1)

      IF (NPRNLI.LE.0.OR.NTIME.EQ.0) THEN
C  TURN OFF TIME DEP MODE IF EITHER NTIME=0 OR NPRNLI=0
        IF (NPRNLI.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NTIME=0'
          NPRNLI=0
        ENDIF
        IF (NTIME.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NPRNLI=0'
          NTIME=0
        ENDIF
        GOTO 1350
      ENDIF

cdr  from here: time dep mode is ON: NTIME.NE.0. and also: NPRNLI GT 0
      LDEF_TIME_HORIZON = .FALSE.
      READ (IUNIN,'(A72)') ZEILE
      IREAD=1
      IF (ZEILE(1:1).EQ.'*') THEN
cdr no time-horizon found in block 13.
C  DATA FOR DEFAULT TIME HORIZON
        NTIME=1
        DTIMV=1.D0
        TIME0=0.D0
        NPTST=0
        NTMSTP=1
        NSNVI=0
        IREAD=0
        LDEF_TIME_HORIZON = .TRUE.
        IF (ZEILE(1:3).EQ.'***') THEN
          IREAD=1
          GOTO 1350
        ELSE
          GOTO 1310
        ENDIF
      ENDIF
cdr  skip reading further comment lines...tbd
c
      READ (ZEILE,6666) NPTST,NTMSTP
      IREAD=0

C   ENFORCE ONE-BY-ONE RELAUNCH FROM CENSUS, IN CASE NLMOVIE
      IF (NLMOVIE.AND.NPTST.GE.0) THEN
        WRITE (IUNOUT,*) 'NPTST RESET TO -1, BECAUSE OF NLMOVIE OPTION'
        NPTST=-1
      ENDIF

      READ (IUNIN,6664) DTIMV,TIME0
c
C   READ DATA FOR SNAPSHOT TALLIES (OPTIONAL)
 1310 IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IREAD = 1
      IF (ZEILE(1:3).EQ.'***') GOTO 1350
      CALL EIRENE_MASAGE('*** 13A. DATA FOR SNAPSHOT TALLIES')
      READ (IUNIN,6666) NSNVI
      IREAD=0
      WRITE (iunout,*) '        NSNVI= ',NSNVI
      CALL EIRENE_LEER(1)
      DO 1320 J=1,NSNVI
cdr  scaling and labels of (additional) snapshot tallies, scored in UPNUSR.F
 1321   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 1321
        READ (ZEILE,6666) ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
        IREAD=0
!       READ (IUNIN,'(A72)') TXTTAL(J,NTALT)
!       READ (IUNIN,'(2A24)') TXTSPC(J,NTALT),TXTUNT(J,NTALT)
        READ (IUNIN,'(A72)') TXTTLT(J)
        READ (IUNIN,'(2A24)') TXTSCT(J),TXTUTT(J)
 1320 CONTINUE
C
      IF (NTIME.LE.0) THEN
        WRITE (iunout,*) 'ERROR IN INPUT: TIME DEP. MODE BUT NTIME.LE.0'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

 1350 CONTINUE
C  SKIP READING REST OF THIS BLOCK
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3).NE.'***') GOTO 1350
      IREAD=1
C
      IF (IREAD.EQ.0) READ (IUNIN,*)

 6662 FORMAT (L1,1X,A24,1X,I1,1X,4(2I3,1X))
 6664 FORMAT (6E12.4)
 6665 FORMAT (12(5L1,1X))
 6666 FORMAT (12I6)
66662 FORMAT (2E12.4,10I6)
66664 FORMAT (I6,6X,5E12.4)
66665 FORMAT (2I6,4E12.4)
66666 FORMAT (I2,1X,A8,12(I3),1X,A10,1X,I2)

      RETURN
C
C  ERROR EXITS
C
  990 CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OR PLOT OF'
      WRITE (iunout,*) 'VOLUME-AVERAGED TALLIES OUT OF RANGE'
      CALL EIRENE_EXIT_OWN(1)
  991 CONTINUE
      WRITE (iunout,*) 'TALLY NUMBER FOR PRINTOUT OF SURFACE-AVERAGED'
      WRITE (iunout,*) 'TALLIES OUT OF RANGE'
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (iunout,*)
     .  'INCONSISTENCY WRT. SURFACE REFLECTION DATABASE, BLOCK 6 '
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: NRPLG.NE.NP2ND, BUT NLPOL=TRUE'
      WRITE (iunout,*) 'NRPLG,NP2ND ',NRPLG,NP2ND
      CALL EIRENE_EXIT_OWN(1)
  995 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: ',
     .            'WRONG UNIT NUMBER FOR OUTPUT OF TALLY SPECIFIED'
      WRITE (iunout,*) 'NTLV,NTLVF',NTLV,NTLVF
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT BLOCK FOR ADDITIONAL DATA FOR'
      WRITE (iunout,*) 'SPECIFIC ZONES FOUND AT ZONE NO. ',I
      CALL EIRENE_EXIT_OWN(1)

      END SUBROUTINE EIRENE_READ_FIXFORM


      SUBROUTINE EIRENE_READ_BLK14_FIXFORM(IERROR)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CTEXT
      USE EIRMOD_CSPEI
      USE EIRMOD_INFCOP, ONLY: EIRENE_IF0COP
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: IERROR
      CHARACTER(420) :: ZEILE
      INTEGER :: IDUM, NCOPIE, J
      LOGICAL :: LSHORT
      EXTERNAL :: EIRENE_LEER, EIRENE_MASAGE, EIRENE_MASPRM,
     .            EIRENE_EXIT_OWN
C
C
C  ADDITIONAL INPUT FOR THIS RUN COMES FROM EITHER
C  ANOTHER CODE (DATA FILE) OR FROM AN EARLIER RUN OF EIRENE
C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"')

      IF (NMODE.EQ.0) THEN
C  STANDALONE RUN, READ BLOCK *** 14 HERE
        WRITE (iunout,*) '        SUBR. INFCOP NOT CALLED.'
!pb NCOPII is completely redundant
!       READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        READ (IUNIN,6666) NAINI,IDUM,NCOPIE
        NCPVI=NCOPIE
        WRITE (iunout,*) '        NAINI, NCPVI = ',NAINI,NCPVI
        IF (NAINI.GT.NAIN) THEN
          CALL EIRENE_MASPRM('NAIN',4,NAIN,'NAINI',5,NAINI,IERROR)
!PB          GOTO 1500
          WRITE (iunout,*) IERROR,' INPUT OR PARAMETER ERRORS DETECTED'
          WRITE (iunout,*)
     .      ' SEE THE ERROR MESSAGES LISTED ABOVE AND CORRECT'
          WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

        IF (NCPVI.GT.NCPV) THEN
          CALL EIRENE_MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

        DO 3020 J=1,NAINI
 3021     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') GOTO 3021
          READ (ZEILE,6666) NAINS(J),NAINT(J)
          READ (IUNIN,'(A72)') TXTPLS(J,NTALN)
          READ (IUNIN,'(2A24)') TXTPSP(J,NTALN),TXTPUN(J,NTALN)
 3020   CONTINUE

      ELSEIF (NMODE.NE.0) THEN
C  COUPLED RUN, READ BLOCK *** 14 IN INTERFACING ROUTINE INFCOP (ENTRY IF0COP)
        NCPVI=0
!pb     NCOPII=0
        NCOPIE=0
        NAINI=0
C  READ BLOCK 14 AND GEOMETRY FROM EXTERNAL DATABASE (FT30),
C  also set NAINI, NCOPII, NCOPIE there
cdr  possibly modified as well:  nlshrt13
        LSHORT=.FALSE.
        CALL EIRENE_IF0COP(.TRUE.,LSHORT)
        WRITE (iunout,*) 'SUBR. INFCOP IS CALLED.'
!pb     write (iunout,*) 'NAINI, NCOPII, NCOPIE, NCPVI ',
!pb  .                    NAINI, NCOPII, NCOPIE, NCPVI
        write (iunout,*) 'NAINI, NCOPIE, NCPVI ',
     .                    NAINI, NCOPIE, NCPVI
        CALL EIRENE_LEER(1)
      ENDIF

      RETURN
 6666 FORMAT (12I6)
      END SUBROUTINE EIRENE_READ_BLK14_FIXFORM


      SUBROUTINE EIRENE_READ_MPI_STRATEGY_FIXED(LINE)
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_COMPRT, ONLY: IUNIN
      USE EIRMOD_CPES, ONLY: NPRS,
     >    STRATEGY_UNDEFINED, STRATEGY_EMBARRASS,
     >    STRATEGY_ORIGINAL, STRATEGY_APCAS,
     >    STRATEGY_BALANCED, STRATEGY_IMPROVED,
     >    INPUT_DISTRIBUTION_STRATEGY

      IMPLICIT NONE
      CHARACTER(80), INTENT(IN) :: LINE
      CHARACTER(80) :: ZEILE
      EXTERNAL :: EIRENE_LEER, EIRENE_MASAGE

C
C  READ MPI STRATEGY IF PRESENT
C
      DO
        IF (INDEX(LINE,'INFORMATION_FOR_MPI') /= 0) THEN
          CALL EIRENE_MASAGE('*** INFORMATION_FOR_MPI')
          DO
            READ (IUNIN,'(A80)',END=1596) ZEILE
            IF (INDEX(ZEILE,'ORIGINAL') /= 0) THEN
              NPRLL = 1
              input_distribution_strategy = STRATEGY_ORIGINAL
              CALL EIRENE_MASAGE
     .         ('APPLYING "ORIGINAL" PARALLELIZATION STRATEGY')
              goto 1597
            ELSE IF (INDEX(ZEILE,'APCAS') /= 0) THEN
              NPRLL = 2
              input_distribution_strategy = STRATEGY_APCAS
              CALL EIRENE_MASAGE
     .         ('APPLYING "APCAS" PARALLELIZATION STRATEGY')
              goto 1597
            ELSE IF (INDEX(ZEILE,'BALANCED') /= 0) THEN
              NPRLL = 3
              input_distribution_strategy = STRATEGY_BALANCED
              CALL EIRENE_MASAGE
     .         ('APPLYING "BALANCED" PARALLELIZATION STRATEGY')
              goto 1597
            ELSE IF (INDEX(ZEILE,'IMPROVED') /= 0) THEN
              NPRLL = 4
              input_distribution_strategy = STRATEGY_IMPROVED
              CALL EIRENE_MASAGE
     .         ('APPLYING "IMPROVED" PARALLELIZATION STRATEGY')
              goto 1597
            ELSE IF (INDEX(ZEILE,'EMBARRASS') /= 0) THEN
              NPRLL = 0
              input_distribution_strategy = STRATEGY_EMBARRASS
              CALL EIRENE_MASAGE
     .         ('APPLYING "EMBARRASSINGLY PARALLEL" STRATEGY')
              goto 1597
            ELSE IF (INDEX(ZEILE,'AUTOMATIC') /=0) THEN
              CALL EIRENE_MASAGE
     .         ('APPLYING AUTOMATIC PARALLELIZATION STRATEGY')
#if MPI_VERSION < 3
              NPRLL = 1
              input_distribution_strategy = STRATEGY_ORIGINAL
              CALL EIRENE_MASAGE
     .         ('CURRENT AUTOMATIC ASSIGNMENT IS "ORIGINAL"')
#else
              NPRLL = 3
              input_distribution_strategy = STRATEGY_BALANCED
              CALL EIRENE_MASAGE
     .         ('CURRENT AUTOMATIC ASSIGNMENT IS "BALANCED"')
#endif
              goto 1597
            END IF
          END DO
 1596     CONTINUE
          if (input_distribution_strategy.eq.STRATEGY_UNDEFINED) THEN
            CALL EIRENE_MASAGE
     .       ('NO MPI PARALLELIZATION STRATEGY NAME RECOGNIZED')
#if MPI_VERSION < 3
            NPRLL = 1
            input_distribution_strategy = STRATEGY_ORIGINAL
            CALL EIRENE_MASAGE
     .       ('APPLYING DEFAULT "ORIGINAL" STRATEGY')
#else
            NPRLL = 3
            input_distribution_strategy = STRATEGY_BALANCED
            CALL EIRENE_MASAGE
     .       ('APPLYING DEFAULT "BALANCED" STRATEGY')
#endif
          end if
 1597     CONTINUE
          CALL EIRENE_LEER(1)
          GOTO 1599
        END IF
      END DO
      IF (NPRS > 1 .and.
     &    input_distribution_strategy.eq.STRATEGY_UNDEFINED) THEN
        CALL EIRENE_MASAGE
     .   ('NO MPI PARALLELIZATION STRATEGY PROVIDED')
#if MPI_VERSION < 3
        NPRLL = 1
        input_distribution_strategy = STRATEGY_ORIGINAL
        CALL EIRENE_MASAGE
     .   ('APPLYING DEFAULT "ORIGINAL" STRATEGY')
#else
        NPRLL = 3
        input_distribution_strategy = STRATEGY_BALANCED
        CALL EIRENE_MASAGE
     .   ('APPLYING DEFAULT "BALANCED" STRATEGY')
#endif
        CALL EIRENE_LEER(1)
      END IF
 1599 CONTINUE
      RETURN

      END SUBROUTINE EIRENE_READ_MPI_STRATEGY_FIXED
