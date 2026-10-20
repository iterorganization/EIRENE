!  03.08.06: data structure for reaction data redefined
!  25.04.07: reading of rate coefficients from HYDKIN database added
c  changed in 2011: new atomic/molecular data structure introduced,
c                   REACDAT(IR)% ..., replaces array CREAC(...)
C
c    at the end of this routine, for each reaction card, call:
cdr  SET_REACTION_DATA(IR,..)
cdr  jan.14: started to comment, cleanup
cdr  april 2015: further commenting, cleanup, nov. 15: continued
cdr  jan 16: started to document options for asymptotics
!pb  apr 16: extensions to allow more precise comments in AMJUEL, HYDHEL, METHAN,
cdr                                                       H2VIBR and AMMONX
c            data files,
cdr          such as character strings H.xxx
cdr          taken over from ITER-IO branch
!pb  may 16: bug fix to the extensions (resolving problem reading HYDHEL H.3)

cdr:  possible conflict with file fort.29, which is also used in coupling to B2
cdr:  subr. infcop.f, there to provide extra information regarding grid distortion
cdr:  june 16: added H.5 - H.7 options for H_COL case.
cdr            started to clarify extrapolation options for polynomial fits. Not ready
cdr            some comments corrected
cdr   Aug. 16: reading Tmin, Emin from hydhel disabled.
cdr            May have corrupted extrapolation in some cases
c     Sept.16: two new internal subroutines,
c              a) READ_RANGE:  to read validity range information,
c              b) READ_COEFFS: three parameters for each validity boundary, for extrapolation options
c    June  17: read_colrad (for old H-COL option (now CRM)) moved to separate routine.
cdr  Jan   19: filnam=CRM --> CR... to prepare merge with branch ...emis....,
cdr            H, He internal CR codes, formulation I, II (MS resolved or not)
cdr  Feb   19: remove obsolete (and unfinished) option HYDRTC
cdr  Nov.  19: add parameters for range of poly-data: inep, knep
cdr  Nov.  19: add DE=EARRH prefactor to separate Arrhenius factor exp(-DE/T)
cdr            (Add a minus-1st term to the polynomial series).
cdr            from the polynomial fits of rate coefficients.
cdr            And add KER (kinetic energy release per reaction).
cdr            And add DELP (potential energy gap reactants - products).
cdr            bugfix: isw=7 rather than =8 in one place
cdr  Apr. 22 : Add iprftype, which is now already read in read_reaclines.f
C
C
      SUBROUTINE EIRENE_SLREAC (IR, FILNAM, H123, REAC, CRC,
     .                          RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                          RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .                          ELNAME, IZ1, BUNDLING,
cdr from here on: optional input: CR models, and filnam=CONST option.
     .                          M_popesc, M_upper, M_lower,
     .                          IROW_ESC, ICOL_ESC, POP_ESC,
     .                          M_QEXT,
     .                          IROW_EXT, ITAL_EXT,
     .                          IFTFL, NCOEF, COEF,
     .                          IPRFTYPE)
c
c  open data stream 29 and read atomic dataset no. IR
c          (note: general input-stream/output-stream no. offset ifoff
c           may have been set (for entire eirene run),
c           then stream is "29+ifoff".  default: ifoff=0)
c  and at the end: call to SET_REACTION_DATA.F (in module COMXS)
c                  to fill REACDAT data structure
c
c
C  input
c    IR    : store data on eirene data structure REACDAT(...,...,IR)

c
c
c
C    FILNAM: read A&M data from file filnam,
c            FILNAM=AMJUEL, HYDHEL, METHAN, AMMONX, H2VIBR, CONST: polynomial fits
CC           FILNAM=TAB2D, ADAS: special treatment, see below.
C            FILNAM=CR...: nothing to be done here, use internal CR code xx_colrad.f
c                          currently available: h_colrad.f, he_colrad.f
C
c    H123  : identifier for data type in file FILNAM, e.g. H.1, H.2, H.3, ...


c    REAC  : in case FILNAM = AMJUEL, HYDHEL, METHAN, H2VIBR, AMMONX:
c               number of reaction in data file "filnam", e.g. 2.2.5
c               and parameter fit-flag is found from the datafile (if available)
c    REAC  : in case FILNAM.eq.CONST:
c               reac IS MISUSED AS fit-flag: iftflg.
C               not NICE, VERY CONFUSING.
C               BETTER MAKE AN OWN INPUT PARAMETER IFTFLG IN CASE OPTION FILNAM= "CONST"

cdr what does that mean for CRM? ADAS ?  what about "spectral database"?
cdr where described, where read ?

C            in case FILNAM = TAB2D, ADAS: the file name DSN = REAC_ELNAME.dat is opened
C                                          (stream 29+ifoff)
C                                  and then subroutine read_tab2d.f is called.

c    CRC   : type of process, e.g. EI, CX, EL, PI, RC, OT, PH, etc.
c
c  parameters for extrapolation beyond specified range [RiMN,RiMX, i=1,2] of data (asymptotics),
c  these asymptotics parameters may already have been read from input file, block 4., subroutine input.f
c  or, if not, they will be searched for on the A&M data files read here.
c
c    for i=1,2:
c    RCiMIN: LOG(RiMN), RiMN: lower boundary for indep. dependent variable (energy, temperature, density)
c           Default: RiMN = exp(-20.)
c    RCiMAX: LOG(RiMX), RiMX: upper boundary for indep. dependent variable (energy, temperature, density)
c           Default: RiMX = exp(20.)
c    FPi    Fitting coefficients for extrapolation (three for MIN and three for MAX, each)
c
c    JFEXiMN Flag for selecting extrapolation expression, left end (minimum)
c            =0  : no data yet, try to read extrapolation from atomic data file here
c            else: extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available
c    JFEXiMX Flag for selecting extrapolation expression, right end (maximum)
c            =0  : no data yet, try to read extrapolation from atomic data file here
c            else: extrapolation is set explicitly in input file, block 4a
c                  skip reading extrapolation data from data file, even if they are available

c  specific input, only available in case FILNAM=ADAS
c    ELNAME  : only in case FILNAM=ADAS: the new file name REAC_ELNAME is construced
C    IZ1     : only in case FILNAM=ADAS: ionization stage number within ADAS file
C    BUNDLING: only in case FILNAM=ADAS: bundling scheme name (optional)

C  internal
C    ISW   <-- H123: H.0: ISW=0, H.1: ISW=1, H.2: ISW=2, ..., H.12: ISW=12
C    I0    derived from ISW, initial value of 2nd index in old CREAC arrays

C  output
c    ISWR  : eirene flag for type of process (1,2,...8), coding EI,CX,EL,PI,...,PH,OT

c    CREAC :          (old version) eirene storage array for a&m data CREAC(9,-1:9,IR)
c    REACDAT(IR)%.... (new version) eirene atomic data structure.

c    MODCLF: see below: further information on input a&m data structure
c    DELPOT: ionisation potential difference (for H.10 data),
c            currently handled in input.f. Not nice! also missing still for: H.8, H.9
c
C    IFTFLG=IFTFLG(IR,IH): flag for type of fitting expression ("fit-flag=...")
C    IH  internally derived from ISW, for different types of data:
C          0 for interaction potential, or differential cross-sections (ISW=0)
C          1 for cross-section, (ISW=1)
C          2 for rate coeff, (ISW=2,3,4)
C          3 for mom-weighted rate coeff. (ISW=5,6,7)
C          4 for energy-weighted rate coeff. (ISW=8,9,10)
C          5 for other quantities, population densities, etc.. (ISW=11,12)
c
c
C       IFTFLG(IR,IH) DEFAULTS:

C           CASE IH=0  (H.0):
C       IFTFLG(IR,0)   =2,  FOR INTERACTION POTENTIAL (GEN. MORSE)

C           CASE  IH=1  (H.1):
C       IFTFLG(IR,1)   =0,  CROSS-SECTION (9-POLYNOMIAL)
C                      =3,  cross-section (ionisation/excitation cross-section
C                           formula (METHANE,...)
C           CASE  IH=2,3....,10 (H.2, H.3,....H.10)
C       IFTFLG(IR,...  =0,  FOR RATE COEFFICIENTS (9-POLYNOMIAL, 9X9-DOUBLE POLYNOMIAL)
C                      =10, FOR RATE COEFFICIENTS (CONSTANT)
C                      =100 FOR RATE, not rate coefficient,
c                      =110 FOR RATE, not rate coefficient, (CONSTANT)
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM"
C    DEFINE PARAMETER MODCLF(IR) (5 DIGITS NMLKJ)
C    FIRST DECIMAL  J           =1  POTENTIAL AVAILABLE
C                                   (ON CREAC(..,-1,IR))
C                   J           =0  ELSE
C    SECOND DECIMAL K           =1  CROSS-SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    THIRD  DECIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    FOURTH DECIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FIFTH  DECIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CINIT
c     USE EIRMOD_PHOTON  ! currently not needed
      USE EIRMOD_CTRCEI

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: IR
      CHARACTER(3), INTENT(IN) :: CRC
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REAC

cdr  optional, only needed for CR (COLRAD) format
      INTEGER,      INTENT(IN), OPTIONAL :: M_popesc, M_upper, M_lower
      INTEGER,      INTENT(IN), OPTIONAL :: IROW_ESC(:), ICOL_ESC(:)
      REAL(DP),     INTENT(IN), OPTIONAL :: POP_ESC(:)
      INTEGER,      INTENT(IN), OPTIONAL :: M_QEXT
      INTEGER,      INTENT(IN), OPTIONAL :: IROW_EXT(:), ITAL_EXT(:,:)

!pb  optional, only needed for CONST format
      INTEGER,      INTENT(IN), OPTIONAL :: IFTFL, NCOEF
      REAL(DP),     INTENT(IN), OPTIONAL :: COEF(9)

cdr  optional, only needed for TAB2D (ADAS) format
      INTEGER,      INTENT(IN) :: IZ1
      CHARACTER(2), INTENT(IN) :: ELNAME
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: BUNDLING

cdr  optional, only needed for PHOTONS format
      INTEGER,      INTENT(IN) :: IPRFTYPE

cdr  asymptotics parameters already read from input block 4?
cdr  if not: try to read from external A&M data file
cdr  in either case: store these on data structure REACDAT,
cdr  in call to: set_reaction_data(IR,...)
      INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
      REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                            RC2MIN, RC2MAX, FP2(6)
cdr
      REAL(DP) :: EARRH0, EARRH1, RTMAX, ERTMAX, ETH, KER, DELP
      CHARACTER(50) :: REACSTR
      LOGICAL :: LCONST
!pb   REAL(DP) :: CREACD(9,9)  ! INTERMEDIATE STORAGE FOR FIT PARAMETERS
!PB   AS CREACD SHALL BE ABLE TO HOLD ANY NUMBER OF FIT PARAMETERS
!PB   E.G. (1,1), (12,1) AND (9,9) IT SHOULD BE ALLOCATABLE
      REAL(DP), ALLOCATABLE :: CREACD(:,:)  ! INTERMEDIATE STORAGE FOR FIT PARAMETERS

      INTEGER :: I, J, K, IC, IREAC, ISW,
     .           IFLG, IANF, IFILE, IL, IERR
      CHARACTER(200) :: DSN, DIR
      CHARACTER(1) :: CUT
      CHARACTER(2) :: CC
      LOGICAL :: FOUND

      INTERFACE
        subroutine EIRENE_read_photdbk (ir, reac, isw, iprftype)
        use EIRMOD_precision
        integer, intent(in) :: ir, isw, iprftype
        character(len=*), intent(in) :: reac
        end subroutine EIRENE_read_photdbk

        SUBROUTINE EIRENE_READ_COLRAD (IR,FILNAM,REAC,ISW,
     .                                 M_POPESC, M_UPPER, M_LOWER,
     .                                 IROW_ESC,ICOL_ESC,POP_ESC,
     .                                 M_QEXT,
     .                                 IROW_EXT,ITAL_EXT)
        USE EIRMOD_PRECISION
        INTEGER, INTENT(IN) :: IR, ISW
        CHARACTER(LEN=*), INTENT(IN) :: FILNAM, REAC
c  optional input parameters
        integer, intent(in), optional :: m_popesc, m_upper, m_lower
        integer, intent(in), optional :: irow_esc(:), icol_esc(:)
        real(dp) , intent(in), optional :: pop_esc(:)
        integer, intent(in), optional :: m_qext
        integer, intent(in), optional :: irow_ext(:), ital_ext(:,:)
        END SUBROUTINE EIRENE_READ_COLRAD

        subroutine EIRENE_read_tab2d (ir,reac,isw,iz1)
        use EIRMOD_precision
        integer, intent(in) :: ir, isw, iz1
        character(len=*), intent(in) :: reac
        end subroutine EIRENE_read_tab2d

        SUBROUTINE EIRENE_READ_POLY (IR,FILNAM,H123,REACSTR,ISW,IFLG,
     .                         RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                         RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX)
        USE EIRMOD_PRECISION
        INTEGER,      INTENT(IN) :: IR, ISW, IFLG
        CHARACTER(8), INTENT(IN) :: FILNAM
        CHARACTER(4), INTENT(IN) :: H123
        CHARACTER(LEN=*), INTENT(IN) :: REACSTR
        INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
        REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                              RC2MIN, RC2MAX, FP2(6)
        END SUBROUTINE EIRENE_READ_POLY

        subroutine eirene_lookup_adasdir_usr(DSN, FOUND,
     .                                       REAC_IN, ELNAME, BUNDLING)
        character(*),intent(inout) :: DSN
        character(LEN=*), intent(in) :: REAC_IN
        logical, intent(inout) :: found
        character(LEN=*), intent(in), optional :: ELNAME, BUNDLING
        end subroutine eirene_lookup_adasdir_usr

      END INTERFACE

!rc read photdbk and read_colrad is an interface now
      EXTERNAL :: EIRENE_LEER, EIRENE_MASJ1R, EIRENE_MASRR1,
     .            EIRENE_EXIT_OWN, EIRENE_UPPERCASE
C
      IF(TRCAMD) WRITE(IUNOUT,'(A,1X,I3,1X,A8,1X,A4,1X,A,1X,A2)')
     w                 "IR,FILNAM,H123,REAC,CRC",IR,FILNAM,H123,
     w                  REAC(1:LEN_TRIM(REAC)),CRC !VK

      WRITE (IUNOUT,*) 'SLREAC CALLED'
      WRITE (IUNOUT,*) 'IR ',IR
      WRITE (IUNOUT,*) 'FILNAM ',FILNAM
      WRITE (IUNOUT,*) 'H123 ',H123
      WRITE (IUNOUT,*) 'REAC ',REAC
      CALL EIRENE_LEER(1)
C
C
c  type (class) of reaction process

      ISWR(IR)=0

      CC = ADJUSTL(TRIM(CRC))
      CALL EIRENE_UPPERCASE(CC)
      SELECT CASE(CC)
      CASE('EI','DS')
        ISWR(IR)=1
cdr   CASE('XX')
cdr     ISWR(ir)=2      free, unused
      CASE('CX')
        ISWR(IR)=3
      CASE('II','PI')
        ISWR(IR)=4
      CASE('EL')
        ISWR(IR)=5
      CASE('RC')
        ISWR(IR)=6
      CASE('PH')
         ISWR(IR)=7   !dr fix ph3
      CASE('OT')
         ISWR(IR)=8
      CASE DEFAULT
         WRITE (IUNOUT,*) 'UNRECOGNIZED REACTION TYPE FOR REACTION ',
     .                    'IR = ',IR
         WRITE (IUNOUT,'(A,A)') 'CRC = ',CRC
         CALL EIRENE_EXIT_OWN(1)
      END SELECT
C
      LCONST=.FALSE.
      IF (INDEX(FILNAM,'CONST').NE.0) THEN
cdr  fit data are directly read from input file, not from database
cdr  Read IFTFLG and CREACD(...,1) from IUNIN, further below.
        LCONST=.TRUE.
!  nothing to be done
      ELSEIF (INDEX(FILNAM,'CR').NE.0) THEN
cdr     internal CR-code, H_colrad.f, He_colrad.f, ....
!  nothing to be done
      ELSE   ! in all other cases: open data file, stream 29+ifoff
!  open data file, stream 29+ifoff.
        DO IFILE=1,NDBNAMES
          IF (INDEX(FILNAM,DBHANDLE(IFILE)).NE.0) EXIT
        END DO
        IF (IFILE <= NDBNAMES) THEN
C  proper filnam DBFNAME no. IFILE found
cdr       write (iunout,*) ifile, dbfname(ifile)
          IF (INDEX(FILNAM,'TAB2D') == 0 .AND.
     .        INDEX(FILNAM,'ADAS')  == 0) THEN
! FILNAM=AMJUEL, HYDHEL, METHAN, AMMONX, H2VIBR, PHOTON....: open data file
            DSN=DBFNAME(IFILE)
            inquire (FILE=TRIM(DSN),exist=found)
          ELSEIF (INDEX(FILNAM,'TAB2D').NE.0 .OR.
     .            INDEX(FILNAM,'ADAS') .NE.0) THEN
! FILNAM=TAB2D or FILNAM=ADAS: open data file
! FIRST: FIND NAME OF SPECIFIC TAB2D or ADAS FILE TO BE READ, DSN=abc.dat
! construct 'DSN' from: reac, elname, bundling
            DIR = ' '
            IL = 0
            IF (VERIFY(DBFNAME(IFILE),' ') .NE. 0) THEN
              IC = SCAN(DBFNAME(IFILE),'/\\')
              IF (IC /= 0) THEN
                CUT = DBFNAME(IFILE)(IC:IC)
                IF ((CUT == '').OR.(CUT == ' ')) CUT = '/'
              ELSE
                CUT = '/'
              END IF
              DIR = TRIM(DBFNAME(IFILE)) // CUT //
     .              ADJUSTL(TRIM(REAC)) // CUT
              IL = INDEX(DIR,CUT,.TRUE.)
            END IF
            IF (PRESENT(BUNDLING) .AND. (LEN_TRIM(BUNDLING) > 0)) THEN
              IF (VERIFY(BUNDLING,' ').NE.0) THEN
                IF (IL == 0) THEN
                  DSN = ADJUSTL(TRIM(REAC)) // '_' //
     .                  TRIM(ELNAME) // '_' // TRIM(BUNDLING) // '.dat'
                ELSE
                  DSN = DIR(1:IL) // ADJUSTL(TRIM(REAC)) // '_' //
     .                  TRIM(ELNAME) // '_' // TRIM(BUNDLING) // '.dat'
                END IF
              END IF
            ELSE
              IF (IL == 0) THEN
                DSN = ADJUSTL(TRIM(REAC)) // '_' //
     .                TRIM(ELNAME) // '.dat'
              ELSE
                DSN = DIR(1:IL) // ADJUSTL(TRIM(REAC)) // '_' //
     .                TRIM(ELNAME) // '.dat'
              END IF
            END IF
            inquire (FILE=TRIM(DSN),exist=found)
            if (found)
     >       write (iunout,'(2a)') 'TAB1D OR TAB2D OR ADAS: ',trim(DSN)
          END IF
          if (found) then
            OPEN (UNIT=29+ifoff,FILE=TRIM(DSN),IOSTAT=IERR)
          else
            IERR = 1
          endif
!PB 21.07.2022
!   Here the standard SOLPS tree was searched. The code has been moved into
!   user routine EIRENE_LOOKUP_ADASDIR_USR
          if ( ierr /= 0 .or. .not.found)
     .      call eirene_lookup_adasdir_usr(dsn,found,
     .                                     reac,elname,bundling)

          if(found.and.(INDEX(FILNAM,'TAB2D').NE.0 .OR.
     .                  INDEX(FILNAM,'ADAS') .NE.0)) CALL EIRENE_LEER(1)

C  THE A&M DATA FILE FILNAM IS NOW OPENED, ON STREAM 29 (+ifoff)

        ELSE
          WRITE (iunout,*)
     .      ' NO VALID FILENAME IN REACTION CARD'
          WRITE (iunout,*) ' CHOOSE EITHER'
          WRITE (iunout,*) ' AMJUEL, METHAN, HYDHEL, AMMONX, H2VIBR'
          WRITE (iunout,*) ' OR'
          WRITE (iunout,*) ' TAB1D, TAB2D, ADAS'
          WRITE (iunout,*) ' OR'
          WRITE (iunout,*) ' CR'
          WRITE (iunout,*) ' OR'
          WRITE (iunout,*) ' CONST'
          WRITE (iunout,*) ' OR'
          WRITE (iunout,*) ' PHOTON'
          WRITE (iunout,*) ' FOR ENTERING REACTION DATA VIA'
          WRITE (iunout,*) ' EIRENE INPUT FILE'
          WRITE (iunout,*) ' FILNAM WAS : ', FILNAM
          CALL EIRENE_EXIT_OWN(1)
        END IF
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
c  this must be an integer !
        READ (H123(3:3),'(I1)') ISW
      ELSE
c  this must be an integer !
        READ (H123(3:4),'(I2)') ISW
      ENDIF

C  ENSURE THAT REAC IS NOT EMPTY
C  PREPARE REACSTR FOR READ_POLY
      REACSTR=REPEAT(' ',11)
      IANF=VERIFY(REAC,' ')
      IF (IANF > 0) THEN
        IREAC=INDEX(REAC(IANF:),' ')-1
        IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
        REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
C  ADD ONE MORE BLANK, IF POSSIBLE
        IREAC=IREAC+2
      ELSE
        IF (.NOT.LCONST) THEN
          WRITE (iunout,*) ' NO REACTION SPECIFIED IN REACTION CARD ',IR
          CALL EIRENE_EXIT_OWN (1)
        END IF
      END IF
C
C Set flags IFLG, IFTFLG(IR), MODCLF(IR) indicating the reaction type
crc all string variables moved to read_poly
C  H.0
      SELECT CASE(ISW)
      CASE(0)
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
C  DEFAULT POTENTIAL: GENERALISED MORSE
        IFTFLG(IR,0)=2
C  no asymptotics yet for interaction potentials

C  H.1
      CASE(1)
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
C  DEFAULT CROSS-SECTION: 8TH-ORDER POLYNOMIAL OF LN(SIGMA) VS LN(E)
        IFTFLG(IR,1)=0
C  H.2
      CASE(2)
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH-ORDER POLYNOMIAL OF LN(<SIGMA V>) VS LN(T), FOR E0=0.
        IFTFLG(IR,2)=0
C  H.3
      CASE(3)
        MODCLF(IR)=MODCLF(IR)+200
        IFLG=2
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(E0)
        IFTFLG(IR,2)=0
C  H.4
      CASE(4)
        MODCLF(IR)=MODCLF(IR)+300
        IFLG=2
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(NE)
        IFTFLG(IR,2)=0
C  H.5
      CASE(5)
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
C  H.6
      CASE(6)
        MODCLF(IR)=MODCLF(IR)+2000
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
C  H.7
      CASE(7)
        MODCLF(IR)=MODCLF(IR)+3000
        IFLG=3
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,3)=0
C  H.8
      CASE(8)
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
C  H.9
      CASE(9)
        MODCLF(IR)=MODCLF(IR)+20000
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
C  H.10
      CASE(10)
        MODCLF(IR)=MODCLF(IR)+30000
        IFLG=4
C  ENERGY-WEIGHTED RATE COEFFICIENT
        IFTFLG(IR,4)=0
C  H.11
      CASE(11)
        IFLG=5
C  OTHER COEFFICIENTS (RATIOS, POPULATION COEFFICIENTS, ETC)
        IFTFLG(IR,IFLG)=0
C  H.12
      CASE(12)
        IFLG=5
C   OTHER COEFFICIENTS (RATIOS, POPULATION COEFFICIENTS, ETC)
        IFTFLG(IR,IFLG)=0

      CASE DEFAULT
         WRITE (IUNOUT,*) 'INVALID REACTION TYPE IDENTIFIER FOUND',
     .                    ' IR = ',IR
         WRITE (IUNOUT,*) 'H123 = ', H123
         WRITE (IUNOUT,*) 'REACTION IS IGNORED '
         RETURN
      END SELECT

C
      IF (INDEX(FILNAM,'PHOTON').NE.0) THEN
        CALL EIRENE_READ_PHOTDBK (IR,REAC,ISW,IPRFTYPE)
        RETURN
      END IF

      IF (INDEX(FILNAM,'CR').NE.0) THEN
        CALL EIRENE_READ_COLRAD (IR,FILNAM,REAC,ISW,
     .                           M_POPESC, M_UPPER, M_LOWER,
     .                           IROW_ESC,ICOL_ESC,POP_ESC,
     .                           M_QEXT,
     .                           IROW_EXT,ITAL_EXT)
c  close unit=29+ifoff: done in READ_COLRAD.f
        RETURN
      END IF

      IF (INDEX(FILNAM,'TAB2D').NE.0 .OR.
     .    INDEX(FILNAM,'ADAS') .NE.0) THEN
        CALL EIRENE_READ_TAB2D (IR,REAC,ISW,IZ1)
c  close unit=29+ifoff: done in READ_TAB2D.f
        RETURN
      END IF

C
      IF (INDEX(FILNAM,'CONST').NE.0) THEN
cdr  'CONST' is an A&M data model in which fit parameters are directly
cdr  read from input file block 4, rather than via parsing an external data file
cdr  Oct. 18: rationalization: as with other data, (FILNAM options)
cdr  we now transfer the parameters IFTFL, NCOEF, COEF(NCOEF) as optional
cdr  parameters via the parameter list, rather than reading from input file IUNIN
cdr  here, as it was the case before.
        IF (PRESENT(IFTFL).AND.PRESENT(NCOEF).AND.PRESENT(COEF)) THEN
           IFTFLG(IR,IFLG) = IFTFL
          ALLOCATE(CREACD(1:NCOEF,1))
          CREACD(1:NCOEF,1) = COEF(1:NCOEF)

cdr Put these coefficients onto "poly" data structure of REACDAT(IR)
cdr Probably: only 1 parameter fits allowed here?
          CALL EIRENE_SET_REACTION_DATA    ! this routine sets only
                                           ! "POLY" data
     .          (IR,ISW,IFTFLG(IR,IFLG),CREACD,
     .           NCOEF,1,
     .           IUNOUT,.FALSE.)
c  no optional extrapolation flags here
          IF (TRCAMD) THEN
             WRITE(IUNOUT,*) 'CREACD'
             DO I=1,NCOEF
               WRITE(IUNOUT,'(7X,I1,1P,(1X,ES12.4))')
     .                       I,CREACD(I,1)
             ENDDO
          ENDIF
          DEALLOCATE (CREACD)
          RETURN
        ELSE
          WRITE (iunout,*) 'CONSTANT REACTION REACTION REQUESTED'
          WRITE (iunout,*) 'BUT NO DATA AVAILABLE FOR IR ', IR
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ENDIF
C
C  READ FROM DATA FILE, stream 29
C
C  already ruled out here (done at this point):
C  FILNAM= "CR", "CONST", "ADAS", "TAB2D", "PHOTON"
C  in all these cases: already returned to calling program
C
C......................................................................
C  AT THIS POINT: FILNAM= AMJUEL, HYDHEL, AMMONX, H2VIBR, METHAN,
C     i.e. single or double polynomial fits

      call EIRENE_read_poly (IR, FILNAM, H123, REACSTR, ISW, IFLG,
     .                       RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                       RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX)
C
      RETURN

      END SUBROUTINE EIRENE_SLREAC
