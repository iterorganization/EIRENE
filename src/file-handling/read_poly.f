cdr  cut and paste from slreac.f (pb, summer 2020? Which slreac.f version? ) 
cdr  pb: plus some (Fortran) code modernisation.

cdr Dec. 20:  remove old parts, add comments..., cleanup
cdr           manually merge with own version, which starts from slreac.f, 
cdr           dr-wip-isfn branch Nov. 20)
cdr           i0 removed
cdr           iflg, modflf, iftflg already in calling program slreac.f
cdr Jan. 22:  default density range for H.4, H.7, H.10 added, because
cdr           the hard coded density cut off for all IFIT options was too restrictive, 
cdr           e.g for internal CR codes which cover the full density range
cdr           algebraically correct (ifit=5 option)

      SUBROUTINE EIRENE_READ_POLY (IR,FILNAM,H123,REACSTR,ISW,IFLG,
     .                          RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .                          RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX)
cdr  IFIT = 1,2 option
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
C    FILNAM: read A&M data from file filnam,
c            FILNAM=AMJUEL, HYDHEL, METHAN, AMMONX, H2VIBR, CONST: polynomial fits
c    H123  : identifier for data type in file FILNAM, e.g. H.1, H.2, H.3, ...
c    REACSTR: in case FILNAM = AMJUEL, HYDHEL, METHAN, H2VIBR, AMMONX:
c             number of reaction in data file "filnam", e.g. 2.2.5, etc.
c             and parameter fit-flag is found from the datafile (if available)
c    ISW   :  indicator for data type in file FILNAM, =0,1,2,...,12
c             already derived from H123 (H.0, H.1,...) in calling routine.
c    IFLG  :  indicator for data type, =0,1,2,3,4,5, used in FITFLAG
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

C  in calling program done already:
C    ISW   <-- H123: H.0: ISW=0, H.1: ISW=1, H.2: ISW=2, ..., H.12: ISW=12
C    IFLAG <-- ISW :   

C  output
c    REACDAT(IR)%.... (new version) eirene atomic data structure.

c    DELPOT: ionisation potential (for H.10 data),
c            currently handled in input.f. not nice! also missing still for: H.8, H.9
c
C    IFTFLG=IFTFLG(IR,IFLG): flag for type of fitting expression ("fit-flag=...")
C    IFLG  internally derived from ISW, for different types of data:
C          0 for interaction potential, or differential cross-sections (ISW=0)
C          1 for cross-section, (ISW=1)
C          2 for rate coeff, (ISW=2,3,4)
C          3 for mom-weighted rate coeff. (ISW=5,6,7)
C          4 for energy-weighted rate coeff. (ISW=8,9,10)
C          5 for other quantities, population densities, etc.. (ISW=11,12)
c
c
C       IFTFLG(IR,IFLG) DEFAULTS:

C           CASE IFLG=0  (H.0):
C       IFTFLG(IR,0)   =2,  FOR INTERACTION POTENTIAL (GEN. MORSE)

C           CASE  IFLG=1  (H.1):
C       IFTFLG(IR,1)   =0,  CROSS-SECTION (9-POLYNOMIAL)
C                      =3,  cross-section (ionisation/excitation cross-section
C                           formula (METHANE,...)
C           CASE  IFLG=2,3....,10 (H.2, H.3,....H.10)
C       IFTFLG(IR,IFLG)=0,  FOR RATE COEFFICIENTS (9-POLYNOMIAL, 9X9-DOUBLE POLYNOMIAL)
C                      =10, FOR RATE COEFFICIENTS (CONSTANT)
C                      =100 FOR RATE, not rate coefficient,
c                      =110 FOR RATE, not rate coefficient, (CONSTANT)
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CINIT
c     USE EIRMOD_PHOTON  ! currently not needed
      USE EIRMOD_CTRCEI

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: IR, ISW, IFLG
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REACSTR

cdr  asymptotics parameters already read from input block 4?
cdr  if not: try to read from external A&M data file
cdr  in either case: store these on data structure REACDAT,
cdr  in call to: set_reaction_data(IR,...)
      INTEGER,  INTENT(IN OUT) :: JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
      REAL(DP), INTENT(IN OUT) :: RC1MIN, RC1MAX, FP1(6),
     .                            RC2MIN, RC2MAX, FP2(6)
cdr
      REAL(DP) :: EARRH0, EARRH1, RTMAX, ERTMAX, ETH, KER, DELP

cdr   temporary storage array, CREACD, SHALL BE ABLE TO HOLD ANY NUMBER OF POLYNOMIAL FIT PARAMETERS
!PB   E.G. (1,1), (12,1) AND (9,9)
      REAL(DP), ALLOCATABLE :: CREACD(:,:) ! INTERMEDIATE STORAGE FOR FIT PARAMETERS

cdr  for reading asymptotics parameters from data files
      INTEGER  :: IF1MN, IF1MX, IF2MN, IF2MX
      REAL(DP) :: FP1L(3), FP1R(3), FP2L(3), FP2R(3)
      REAL(DP) :: R1MN, R1MX, R2MN, R2MX

      INTEGER :: I, IND, J, K, IH, IC, IREAC, INDFF,
     .           IANF, IFILE, IL, INDG, INDADD,
     .           INI, INE, KNI, KNE,
     .           INIP, INEP, KNIP, KNEP      ! range of non-zero fit parameters
      CHARACTER(80) :: ZEILE, LAST_TEX, ULINE
      CHARACTER(4) :: CHR, CHRA, CETH, CKER, CDEL, BEND
      CHARACTER(3) :: CH1L, CH1R, CH2L, CH2R
      CHARACTER(1) :: BACK
      CHARACTER(8) :: SECTION
      CHARACTER(9) :: FITFLAG
      CHARACTER(7) :: C1L, C1R, C2L, C2R, CMR, CEMR
      LOGICAL :: LGC1MIN,LGC1MAX,LGC2MIN,LGC2MAX,
     .           LGR1MIN,LGR1MAX,LGR2MIN,LGR2MAX

! defining backslash character
      BACK=ACHAR(92)
      SECTION=BACK // 'section'
      BEND=BACK // 'end'
C
!   set some defaults

      CHR=' l0 '

cdr header of dataset (optional)
      FITFLAG ='fit-flag '
cdr Arrhenius prefactor: EXP(-EARRH/T), for low T limit of rates
cdr this coefficient can appear as "-1st" coeff, e.g b-1
      EARRH0  = 0.0
      EARRH1  = 0.0  ! not yet in use, for density dependent rates with
cdr                    density dependent threshold coefficient.

cdr end of dataset (optional)
c  some additional  (optional) reaction data:  threshold energy,
c                                              max rate coeff sigma*v_rel,
c                                              at E_rel=ERTMAX
      CMR  = 'MAXRATE'
      CEMR = 'ELAB'
      CETH = 'ETH'
      CKER = 'KER'
      CDEL = 'DELP'

      RTMAX = 0._DP
      ERTMAX = -HUGE(1._DP)
      ETH = 0._DP
      KER = 0._DP
      DELP = 0._DP

C DEFAULT FITTING PARAMETER RANGE FOR POLYNOMIAL FITS: 8th order polynoms.
      INIP=1
      INEP=9
      KNIP=1
      KNEP=9

C     Defaults: no asymptotics
      CH1L='ll0'
      CH1R='lr0'
      CH2L='lb0'
      CH2R='lt0'
      C1L = 'XXMIN'
      C1R = 'XXMAX'
      C2L = 'YYMIN'
      C2R = 'YYMAX'
C
C  Set character string identifiers CHR, etc., 
C  to parse numerical coefficients in data files with polynomial fit format. 
C
C Further flags IFLG, IFTFLG(IR), MODCLF(IR) are already set in calling program SLREAC
C as these are not specific to a particular data format. 
 
C  H.0
      SELECT CASE(ISW)
      CASE(0)
C  DEFAULT POTENTIAL: GENERALISED MORSE
        CHR=' p0 '
        CHRA=' p-1'
C  no asymptotics yet for interaction potentials

C  H.1
      CASE(1)
C  DEFAULT CROSS-SECTION: 8TH-ORDER POLYNOMIAL OF LN(SIGMA) VS LN(E)
        CHR=' a0 '
        CHRA=' a-1'
c  (laboratory) energy range, asymptotics
        CH1L='al0'
        CH1R='ar0'
        C1L = 'ELABMIN'
        C1R = 'ELABMAX'

C  H.2
      CASE(2)
C  DEFAULT RATE COEFFICIENT: 8TH-ORDER POLYNOMIAL OF LN(<SIGMA V>) VS LN(T), FOR E0=0.
        CHR=' b0 '
        CHRA=' b-1'
c  temperature range, asymptotics
        CH1L='bl0'
        CH1R='br0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.3
      CASE(3)
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(E0)
        CHR=' c0 '
        CHRA=' c-1'
c  temperature range, asymptotics
c  beam energy range, asymptotics
        CH1L='cl0'
        CH1R='cr0'
        CH2L='cb0'
        CH2R='ct0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.4
      CASE(4)
C  DEFAULT RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(NE)
        CHR=' d0 '
        CHRA=' d-1'
c  temperature range, asymptotics
c  density range, asymptotics
        CH1L='dl0'
        CH1R='dr0'
        CH2L='db0'
        CH2R='dt0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.5
      CASE(5)
C  MOMENTUM-WEIGHTED RATE COEFFICIENT
        CHR=' e0 '
        CHRA=' e-1'
c  temperature range, asymptotics
        CH1L='el0'
        CH1R='er0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.6
      CASE(6)
C  MOMENTUM-WEIGHTED RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(E0)
        CHR=' f0 '
        CHRA=' f-1'
c  temperature range, asymptotics
c  beam energy range, asymptotics
        CH1L='fl0'
        CH1R='fr0'
        CH2L='fb0'
        CH2R='ft0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.7
      CASE(7)
C  MOMENTUM-WEIGHTED RATE COEFFICIENT: DOUBLE POLYNOMIAL OF LN(<SIGMA V>) VS LN(T) AND LN(NE
        CHR=' g0 '
        CHRA=' g-1'
c  temperature range, asymptotics
c  density range, asymptotics
        CH1L='gl0'
        CH1R='gr0'
        CH2L='gb0'
        CH2R='gt0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.8
      CASE(8)
C  ENERGY-WEIGHTED RATE COEFFICIENT
        CHR=' h0 '
        CHRA=' h-1'
c  temperature range, asymptotics
        CH1L='hl0'
        CH1R='hr0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'

C  H.9
      CASE(9)
C  ENERGY-WEIGHTED RATE COEFFICIENT
        CHR=' i0 '
        CHRA=' i-1'
c  temperature range, asymptotics
c  beam energy range, asymptotics
        CH1L='il0'
        CH1R='ir0'
        CH2L='ib0'
        CH2R='it0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'E2MIN'
        C2R = 'E2MAX'
C  H.10
      CASE(10)
C  ENERGY-WEIGHTED RATE COEFFICIENT
        CHR=' j0 '
        CHRA=' j-1'
c  temperature range, asymptotics
c  density range, asymptotics
        CH1L='jl0'
        CH1R='jr0'
        CH2L='jb0'
        CH2R='jt0'
        C1L = 'T1MIN'
        C1R = 'T1MAX'
        C2L = 'N2MIN'
        C2R = 'N2MAX'
C  H.11
      CASE(11)
C  OTHER COEFFICIENTS (RATIOS, POPULATION COEFFICIENTS, ETC)
        CHR=' k0 '
        CHRA=' k-1'
c  asymptotics, one parameter
        CH1L='kl0'
        CH1R='kr0'
        C1L = 'P1MIN'
        C1R = 'P1MAX'
C  H.12
      CASE(12)
C   OTHER COEFFICIENTS (RATIOS, POPULATION COEFFICIENTS, ETC)
        CHR=' l0 '
        CHRA=' l-1'
c  asymptotics, two parameters
        CH1L='ll0'
        CH1R='lr0'
        CH2L='lb0'
        CH2R='lt0'
        C1L = 'P1MIN'
        C1R = 'P1MAX'
        C2L = 'P2MIN'
        C2R = 'P2MAX'

      CASE DEFAULT
         WRITE (IUNOUT,*) 'INVALID REACTION TYPE IDENTIFIER FOUND',
     .                    ' IR = ',IR
         WRITE (IUNOUT,*) 'ISW = ', ISW
         WRITE (IUNOUT,*) 'REACTION IS IGNORED '
         RETURN
      END SELECT

C
C  READ FROM DATA FILE, stream 29+ifoff
C
C......................................................................

C  now identify proper dataset for reaction IR within file FILNAM

      ZEILE=REPEAT(' ',80)
      DO WHILE (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0)
        READ (29+ifoff,'(A80)',END=990) ZEILE
      END DO

      LAST_TEX=REPEAT(' ',80)
    1 READ (29+ifoff,'(A80)',END=990) ZEILE
      IF (INDEX(ZEILE,H123).EQ.0) THEN
        IF (INDEX(ZEILE,BACK).NE.0) LAST_TEX=ZEILE
        GOTO 1
      ELSE
        IF ((INDEX(LAST_TEX,SECTION) .EQ. 0) .AND.
     .      (INDEX(ZEILE,SECTION) .EQ. 0)) GOTO 1
      END IF
C
      IREAC = VERIFY(REACSTR,' ',.TRUE.)+1
      DO
cdr  end of file: goto 990
        READ (29+ifoff,'(A80)',END=990) ZEILE

        IF (INDEX(ZEILE,'H.').NE.0 .and.
     .      INDEX(ZEILE,'section').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').GT.0 .and.
     .      INDEX(ZEILE,REACSTR(1:ireac)).GT.0) EXIT ! infinite loop possible
      END DO

C
C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
C
      SELECT CASE (ISW)

      CASE(0, 1, 2, 5, 8, 11)

    3   READ (29+ifoff,'(A80)',END=990) ZEILE
        INDFF=INDEX(ZEILE,FITFLAG)
c  skip empty lines in header
        IF (INDEX(ZEILE,CHR)+INDEX(ZEILE,CHRA)+INDFF.EQ.0) GOTO 3

c  input line found which either contains fit-flag FITFLAG,
c       or the reaction identifier a0,b0,...k0
c       or the Arrhenius prefactor identifiers a-1,b-1,...k-1
        IF (INDFF > 0) THEN
c  read parameter for type of fitting expression from data file
c  OTHERWISE: use DEFAULT FOR FIT-FLAG: iftflg = 0
          READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
          GOTO 3
        ENDIF
        IF (INDEX(ZEILE,CHRA) > 0) THEN
c  read Arrhenius parameter
          READ (ZEILE((INDFF+8):80),*) EARRH0
          IF (EARRH0.GT.0.0) INIP=0 ! RATHER THAN DEFAULT: INIP=1
          GOTO 3
        ENDIF

C  read only one parameter: (FIT-FLAG = 10, 110, ....)
        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
          IND=INDEX(ZEILE,CHR(2:2)) ! position of a,b,...k
          ALLOCATE (CREACD(1,1),SOURCE=0._DP)
          READ (ZEILE((IND+2):80),'(E20.12)') CREACD(1,1)
          INEP=1    ! RATHER THAN DEFAULT: INEP=9
          KNEP=1

        ELSE
C  READ 9 FIT COEFFICIENTS, SEPARATED BY 'CHR' FIXED FORMAT E20.12
C  THREE LINES WITH THREE DATA PER LINE
c  READ ONLY THE NON-ZERO PARAMETERS
          ALLOCATE (CREACD(9,1),SOURCE=0._DP)
          INEP=0
          KNEP=1
          JLOOP: DO J=0,2
            IND=0
            DO I=1,3
              INDADD=INDEX(ZEILE((IND+1):80),CHR(2:2))
              IF (INDADD.EQ.0) EXIT JLOOP
              INEP=INEP+1
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(2:2))
              READ (ZEILE((IND+2):80),'(E20.12)') CREACD(J*3+I,1)
            ENDDO
            READ (29+ifoff,'(A80)',END=990) ZEILE
          ENDDO JLOOP
        END IF
C

C  SINGLE PARAMETER POLYNOMIAL FITS: DONE

C   AT THIS POINT WE HAVE STORED FOR REACTION IR, DATA TYPE iflg: 0,...,5
C   IFTFLG(IR,IFLG) (DEFAULT: = 0)
C   UP TO 9 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1)
C   AND POSSIBLY (SOME OF) THE EXTRAPOLATION PARAMETERS RCMIN,RCMAX, FP(1:6)
C


C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12

      CASE (3, 4, 6, 7, 9, 10, 12)
cdr  this coding is overly complicated. We need not search for
cdr  for header parts (FITFLAG, ARRHENIUS FACTORS) for each block of
cdr  9 lines. tbd: sync with single parameter fit reading above.

        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
          ALLOCATE (CREACD(1,1),SOURCE=0._DP)
          INEP=1
          KNEP=1
        ELSE
          ALLOCATE (CREACD(9,9),SOURCE=0._DP)
cdr  tentative max size of polynomial fit coef. table.
cdr  possibly reduced values INEP,KNEP will be found below.
          INEP=9
          KNEP=9
        END IF

C READ 3 BLOCKS "J" OF DATA. Each block contains 9 LINES, 3 numbers per line, i.e. 3 sub blocks
C look for FITFLAG or Arrhenius prefactors (tbd) prior to first block
        DO J=0,2
C  SEARCH FOR STRING 'fit-flag'  OR STRING 'Index'
   16     READ (29+ifoff,'(A80)',END=990) ZEILE

cdr tbd: also search for two Arrh. factors, Ethmin and Ethmax wrt to second fit parameter
          INDFF=INDEX(ZEILE,FITFLAG)
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
c  FITFLAG parameter found.
c  read parameter for type of fitting expression from data file
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 16
          ENDIF
cdr  At this point: 'Index'  (mostly: 'E-Index' or 'Ne-Index' was found.

          READ (29+ifoff,'(1X)')  ! skip reading one line (mostly: 'T-Index')
cdr start reading fit parameters for block J here:
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C  IFTFLG = 10, 110, 210,....ETC: READ ONLY ONE CONSTANT PARAMETER
            READ (29+ifoff,*) IH,CREACD(1,1)
            EXIT
          ELSE
            DO I=1,9
C   READ 9 LINES, THREE DATA EACH LINE, UNFORMATTED I.E. READ 3 SUB-BLOCKS K,K+1,K+2
              READ (29+ifoff,*) IH,(CREACD(I,K),K=J*3+1,J*3+3)
            END DO  ! I
c    first  index I: I-th block, vertical, Temp. dependence
c    second index K: from sub block to sub-block (horizontal), ne, eb dependence.
c  d.h. first sub block corresponds to  ln(ne/1e8))=0, ie. ne=1e8, corona rate vs. T
c                                   or  ln(EB)=0, ie. Eb=1 eV
          END IF
        END DO  ! J
        READ (29+ifoff,'(A80)',END=990) ZEILE
C   At this point we have stored for reaction IR: DATA TYPE (IFLG: 0,...,5)
C   IFTFLG(IR,IFLG) (DEFAULT: = 0)
C   81 FIT COEFFICIENTS ON INTERMEDIATE ARRAY CREACD(1...9,1...9)

C
C  DOUBLE PARAMETER POLYNOMIAL FITS: DONE

      END SELECT


C  NEXT: READ ASYMPTOTICS INFORMATION FROM ATOMIC DATA FILE
C        HYDHEL, AMJUEL, AMMONX, H2VIBR, METHANE.

C FOR 1D OR 2D DATASETS. 4 BOUNDARIES, LEFT1, RIGHT1, LEFT2, RIGHT2.
C FOR 1D: ONLY "LEFT1" AND "RIGHT1" ARE USED

      IF (ISW.EQ.0) GOTO 2000    ! NO ASYMPTOTICS FOR POTENTIALS

c  INDICATE, IF EXTRAPOLATION INFORMATION (r1mn,.., if1mn,...) IS FOUND ON DATA FILE
c  DEFAULT: NO DATA FOUND
      LGR1MIN=.FALSE.
      LGR1MAX=.FALSE.
      LGR2MIN=.FALSE.
      LGR2MAX=.FALSE.

c  INDICATE, IF EXTRAPOLATION COEFFICIENTS (fp1l,fp1r,fp2l,fp2r) ARE FOUND ON DATA FILE
c  DEFAULT: NO DATA FOUND
      LGC1MIN=.FALSE.
      LGC1MAX=.FALSE.
      LGC2MIN=.FALSE.
      LGC2MAX=.FALSE.

C  FLAG FOR CHOICE OF EXTRAPOLATION OPTION:
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=0): NO ASYMPTOTICS
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=1): SET TO ZERO BEYOND LAST VALID POINT
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=4): TAKE LAST VALID POINT AT r1mn,r1mx,....
C  DEFAULT ASYMPTOTIC EXPRESSION  (...=5): 2ND-ORDER POLYNOMIAL BEYOND LAST VALID POINT
C                                          (OLD DEFAULT FOR CROSS-SECTIONS)
      IF1MN = 0
      IF1MX = 0
      IF2MN = 0
      IF2MX = 0

C  FURTHER PARAMETERS, NOT RELATED TO ASYMPTOTICS
      RTMAX = 0._DP
      ERTMAX = -HUGE(1._DP)

cdr Try to read asymptotics, unless already read: jfeximn,jfeximx /= 0
cdr Even if jfeximn,jfeximx /= 0, read anyway, but later: do not use

c  data for one "reaction IR" are located between "\begin" and "\end"
!     BEND="\end" stops looking for asymptotics
      DO WHILE(INDEX(ZEILE,BEND) == 0)

        ULINE = ZEILE
        CALL EIRENE_UPPERCASE(ULINE)
c
c  at this point we have found a card ZEILE which contains
c  one of the extrapolation parameter identifiers al0,ar0,....,k0l,k0r
c  that corresponds to the H.1, ....H.12 type of data IR.
c  next: read up to three fit coefficients FP.L OR FP.R.
c  CHR(2:2) is set to either character a,b,c,....,or k
c
        IF (INDEX(ZEILE,CH1L) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP1L)
          LGC1MIN = .TRUE.
        ENDIF
        IF (INDEX(ZEILE,CH1R) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP1R)
          LGC1MAX = .TRUE.
        END IF

        IF (INDEX(ZEILE,CH2L) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP2L)
          LGC2MIN = .TRUE.
        END IF
        IF (INDEX(ZEILE,CH2R) /= 0) THEN
          CALL EIRENE_READ_COEFFS (ZEILE,CHR(2:2),FP2R)
          LGC2MAX = .TRUE.
        END IF
c
c  currently foreseen asymptotic data identifiers in data files:
c  c1l,c2l,c1r,c2r: ELABMIN, ELABMAX,
c                   T1MIN,T1MAX,E2MIN,E2MAX,N2MIN,N2MAX,
c                   P1MIN,P1MAX,P2MIN,P2MAX

cdr  1st parameter, low end
        IF (INDEX(ULINE,TRIM(C1L)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C1L,'EXT-FLG',R1MN,IF1MN)
          LGR1MIN = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL1
          IF (IF1MN == 0 .AND. LGC1MIN) IF1MN = 5
c  default extrapolation from r1mn (by constant continuation) will be: jfex1mn=4
          IF (IF1MN == 0) IF1MN = 4
        END IF
cdr  1st parameter, high end
        IF (INDEX(ULINE,TRIM(C1R)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C1R,'EXT-FLG',R1MX,IF1MX)
          LGR1MAX = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPR1
          IF (IF1MX == 0 .AND. LGC1MAX) IF1MX = 5
c  default extrapolation from r1mx (by constant continuation) will be: jfex1mx=4
          IF (IF1MX == 0) IF1MX = 4
        END IF
cdr  2nd parameter, low end
        IF (INDEX(ULINE,TRIM(C2L)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C2L,'EXT-FLG',R2MN,IF2MN)
          LGR2MIN = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPL2
          IF (IF2MN == 0 .AND. LGC2MIN) IF2MN = 5
c  default extrapolation from r2mn (by constant continuation) will be: jfex2mn=4
          IF (IF2MN == 0) IF2MN = 4
        END IF
cdr  2nd parameter, high end
        IF (INDEX(ULINE,TRIM(C2R)) /= 0) THEN
          CALL EIRENE_READ_RANGE (ULINE,C2R,'EXT-FLG',R2MX,IF2MX)
          LGR2MAX = .TRUE.
c  old default: 2nd-order polynomial beyond valid range, with coefs. FPr2
          IF (IF2MX == 0 .AND. LGC2MAX) IF2MX = 5
c  default extrapolation from r2mx (by constant continuation) will be: jfex2mx=4
          IF (IF2MX == 0) IF2MX = 4
        END IF
C
C  ...AND FURTHER REACTION PARAMETERS, NOT RELATED TO ASYMPTOTICS
C     ETH, KER, DELP,
C     RTMAX
C     ERTMAX
C
        IND = INDEX(ULINE,TRIM(CETH))
        IF (IND /= 0) THEN
          READ (ULINE(IND+3:80),*) ETH
        END IF
        IND = INDEX(ULINE,TRIM(CKER))
        IF (IND /= 0) THEN
          READ (ULINE(IND+3:80),*) KER
        END IF
        IND = INDEX(ULINE,TRIM(CDEL))
        IF (IND /= 0) THEN
          READ (ULINE(IND+3:80),*) DELP
        END IF
c search for cards containing: RTMAX, ERTMAX (MAX OF [SIGMA(E)*SQRT(E)] AT E = ERTMAX)
        IND = INDEX(ULINE,TRIM(CMR))
        IF (IND /= 0) THEN
          INDG = INDEX(ULINE,'=')
          READ (ULINE(INDG+1:80),*) RTMAX
          IND = INDEX(ULINE,TRIM(CEMR))
          IF (IND /= 0) THEN
            INDG = IND + INDEX(ULINE(IND:80),'=')
            READ (ULINE(INDG+1:),*) ERTMAX
          END IF
        END IF

        READ (29+ifoff,'(A80)',END=990) ZEILE
      END DO

cdr special treatment: use a default valid density range from amjuel, h2vibr.
cdr                    unless otherwise specified
cdr for H.4, H.7, H.10: special treatment for nemin
      if (.not. lgr2min) then
        if (isw.eq.4.or.isw.eq.7.or.isw.eq.10) then
cdr Jan.22: take 1e8 as default minimal valid density (this was old default)
          LGR2MIN = .TRUE.
          R2MN=1.0D8
          if2mn=4
          write (iunout,*) 'default lower density limit set to 1.0e8'            
        endif
      endif

      if (.not. lgr2max) then
cdr for H.4, H.7, H.10: special treatment for nemax, if no limit is specified
        if (isw.eq.4.or.isw.eq.7.or.isw.eq.10) then
cdr take 1e16 as default maximal valid density
          LGR2MAX = .TRUE.
          R2MX=1.0D16
          if2mx=4
          write (iunout,*) 'default high density limit set to 1.0e16'            
        endif
      endif

c  unless asymptotics are already explicitly defined in input block 4a
c  put asymptotics information into proper (intermediate) data structure:
c  flags:      jfex1mn,jfex1mx,jfex2mn,jfex2mx
c  boundaries: rc1min,rc1max,rc2min,rc2max
c  parameters: fp1(1:3),fp1(4:6),fp2(1:3),fp2(4:6)

      IF (JFEX1MN == 0) THEN
        IF (LGR1MIN .AND. .NOT. LGC1MIN.and.if1mn.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, 'TYPE ',H123
          WRITE (IUNOUT,*) ' LOWER RANGE FOR 1ST PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF1MN,R1MN      ',if1mn,r1mn)
          IF (IF1MN.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR1MIN) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, 'TYPE ',H123
          WRITE (IUNOUT,*) 'LOWER RANGE FOR 1ST PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF1MN,R1MN      ',if1mn,r1mn)
          if (if1mn.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp1l,3,3)
          CALL EIRENE_LEER(1)
        END IF

        IF (LGR1MIN) RC1MIN = LOG(R1MN)
        IF (LGC1MIN) FP1(1:3) = FP1L
        JFEX1MN = IF1MN
        IF (LGR1MIN .AND. LGC1MIN .AND. (JFEX1MN == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .          JFEX1MN = 5
      END IF

      IF (JFEX1MX == 0) THEN
        IF (LGR1MAX .AND. .NOT. LGC1MAX.and.if1mx.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' UPPER RANGE FOR 1ST PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF1MX,R1MX      ',if1mx,r1mx)
          IF (IF1MX.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR1MAX) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'UPPER RANGE FOR 1ST PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF1MX,R1MX      ',if1mx,r1mx)
          if (if1mx.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp1r,3,3)
          CALL EIRENE_LEER(1)
        END IF

        IF (LGR1MAX) RC1MAX = LOG(R1MX)
        IF (LGC1MAX) FP1(4:6) = FP1R
        JFEX1MX = IF1MX
        IF (LGR1MAX .AND. LGC1MAX .AND. (JFEX1MX == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .          JFEX1MX = 5
      END IF

      IF (JFEX2MN == 0) THEN
        IF (LGR2MIN .AND. .NOT. LGC2MIN.and.if2mn.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' LOWER RANGE FOR 2ND PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF2MN,R2MN      ',if2mn,r2mn)
          IF (IF2MN.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR2MIN) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'LOWER RANGE FOR 2ND PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF2MN,R2MN      ',if2mn,r2mn)
          if (if2mn.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp2l,3,3)
          CALL EIRENE_LEER(1)
        END IF
        IF (LGR2MIN) RC2MIN = LOG(R2MN)
        IF (LGC2MIN) FP2(1:3) = FP2L
        JFEX2MN = IF2MN
        IF (LGR2MIN .AND. LGC2MIN .AND. (JFEX2MN == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .         JFEX2MN = 5
      END IF

      IF (JFEX2MX == 0) THEN
        IF (LGR2MAX .AND. .NOT. LGC2MAX.and.if2mx.ge.3.) THEN
          WRITE (IUNOUT,*) ' WARNING FROM SLREAC'
          WRITE (IUNOUT,*) ' REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) ' UPPER RANGE FOR 2ND PARAMETER OF FIT',
     .          ' SPECIFIED BUT',
     .          ' NO COEFFICIENTS FOR EXTRAPOLATION PROVIDED'
          CALL EIRENE_MASJ1R('IF2MX,R2MX      ',if2mx,r2mx)
          IF (IF2MX.EQ.4)
     .      WRITE (IUNOUT,*) 'I.E.: CONTINUATION AS CONSTANT'
          CALL EIRENE_LEER(1)
        ELSEIF (LGR2MAX) THEN
          WRITE (IUNOUT,*) 'ASYMPTOTICS FROM SLREAC'
          WRITE (IUNOUT,*) 'REACTION ',IR, ' TYPE ',H123
          WRITE (IUNOUT,*) 'UPPER RANGE FOR 2ND PARAMETER OF FIT'
          CALL EIRENE_MASJ1R('IF2MX,R2MX      ',if2mx,r2mx)
          if (if2mx.ge.3)
     .      CALL EIRENE_MASRR1('PARAMETERS ',fp2r,3,3)
          CALL EIRENE_LEER(1)
        END IF
        IF (LGR2MAX) RC2MAX = LOG(R2MX)
        IF (LGC2MAX) FP2(4:6) = FP2R
        JFEX2MX = IF2MX
        IF (LGR2MAX .AND. LGC2MAX .AND. (JFEX2MX == 0))
! DEFAULT EXTRAPOLATION=EXP(FP(1)+FP(2)*PARM+FP(3)*PARM**2), 2ND ORDER ON LOG SCALE
     .        JFEX2MX = 5
      END IF
C
 2000 CONTINUE

      CALL EIRENE_SET_REACTION_DATA   ! this routine sets only "POLY" data
     .            (IR,ISW,IFTFLG(IR,IFLG),CREACD,INEP,KNEP,
     .             IUNOUT,.TRUE.,
c  from here on: optional input to SET_REACTION_DATA
     .             RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .             RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .             RTMAX,ERTMAX,ETH,KER,DELP,
     .             EARRH0,EARRH1)
C
 3000 IF (TRCAMD) THEN
        WRITE(IUNOUT,*) 'CREACD'
c  range for first parameter, DEFAULT:
        ini=inip
        ine=inep
cdr 2nd parameter: density or beam energy
        DO I = ini, ine
          kni=knip
          kne=1
c  double fits, two parameters
          if (isw.eq.3 .or. isw.eq.4 .or.   ! RATE COEFF.
     .        isw.eq.6 .or. isw.eq.7 .or.   ! MOMENTUM RATE COEFF.
     .        isw.eq.9 .or. isw.eq.10 .or.  ! ENERGY RATE COEFF.
     .        isw.eq.12)                    ! DENSITY RATIO, RATE COEFF. RATIO
     .        kne=9
          if (I.eq.0) then
cdr  the Arrhenius factor exp(-EARR/T) is factored out of the fit.
            WRITE(IUNOUT,'(7X,I1,1P,(1X,9E12.4))')
     .                       0,EARRH0    ! = CREACD(0,1)
          else
            WRITE(IUNOUT,'(7X,I1,1P,(1X,9E12.4))')
     .                       I,(CREACD(I,K),K=kni,kne)
          endif
        END DO
        CALL EIRENE_LEER(1)
      END IF

      DEALLOCATE (CREACD)
C
      CLOSE (UNIT=29+ifoff)
C
      RETURN
C
  990 WRITE (iunout,*) ' ERROR IN READ_POLY: '
      WRITE (iunout,*) ' NO DATA FOUND FOR REACTION ',H123,
     .            ' ',REACSTR(1:LEN_TRIM(REACSTR)),
     .            ' IN DATASET ',FILNAM
      WRITE (iunout,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29+ifoff)
      CALL EIRENE_EXIT_OWN(1)

      CONTAINS

cdr next: some stuff for parsing A&M data reading in input block 4
cdr       Perhaps useful in other context as well?
cdr       Comments:  mostly guessing. Not sure.

      SUBROUTINE EIRENE_READ_COEFFS (ZEILE,CH,FP)
c  ZEILE is in upper case, and CH is found.
c  CH is: a,b,c,.....k,l, (depending on H.1, H.2, ...H.12)
c  read up to three parameters FP(i), i=1,3 formatted: '(E20.12)'
c  to be used for extrapolation
      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(1), INTENT(IN) :: CH
      REAL(DP), INTENT(OUT) :: FP(3)
      INTEGER :: I, IND, INC

      IND=0
      DO I=1,3
        INC=INDEX(ZEILE((IND+1):80),CH)
        IF (INC.GT.0) THEN
           IND=IND+INDEX(ZEILE((IND+1):80),CH)
           READ (ZEILE((IND+3):80),'(E20.12)') FP(I)
        ENDIF
      END DO
      RETURN

      END SUBROUTINE EIRENE_READ_COEFFS


      SUBROUTINE EIRENE_READ_RANGE (ZEILE,KEY1,KEY2,RNG,IFX)
cdr comments by DR: not sure !
c  called from slreac, after the original fit coefficients for reaction IR
c  are read.
c  At this point an extrapolation card ZEILE belonging to this reaction IR
c  has already been found.
c
c  This routine:
c  Reads validity range from atomic data file:
c  Search in ZEILE for key1, key2 and return: RNG, IFX

c  key1: 'ELABMIN','ELABMAX','T1MIN','T1MAX','E2MIN','E2MAX',
C        'N2MIN','N2MAX','P1MIN','P1MAX','P2MIN','P2MAX'=,
C                    read RNG (unformatted, real)
c  key2: 'EXT-FLG'= ,read IFX (unformatted, integer)

      CHARACTER(80), INTENT(IN) :: ZEILE
      CHARACTER(7), INTENT(IN) :: KEY1, KEY2
      REAL(DP), INTENT(OUT) :: RNG
      INTEGER, INTENT(OUT) :: IFX
      INTEGER :: IND1, IND2, INDE, INDP, INDX, INDA, INDG
      CHARACTER(20) :: FORM

      IND1 = INDEX(ZEILE,TRIM(KEY1))
      IND2 = INDEX(ZEILE,TRIM(KEY2))

      RNG = 0._DP
      IFX = 0

      IF (IND1 > 0) THEN
        INDG = INDEX(ZEILE,'=')
        INDE = INDG + VERIFY(ZEILE(INDG+1:),'+-0123456789DEed. ') - 1
        INDP = SCAN(ZEILE(INDG+1:),'.')
        INDX = SCAN(ZEILE(INDG+1:),'EDed')
        INDA = SCAN(ZEILE(INDG+1:),'+-0123456789.')
        FORM=REPEAT(' ',20)
        WRITE (FORM,'(A2,I0,A1,I0,A1)')
     .         '(E',INDX+3-INDA+1,'.',INDX-INDP-1,')'
        READ (ZEILE(INDG+INDA:INDE),FORM) RNG
      END IF

      IF (IND2 > 0) READ (ZEILE(IND2+7:),*) IFX
      RETURN

      END SUBROUTINE EIRENE_READ_RANGE

      END SUBROUTINE EIRENE_READ_POLY
