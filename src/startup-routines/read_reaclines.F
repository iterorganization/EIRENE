      subroutine eirene_read_reaclines (ZEILE, iunin, iunout, iread,
     .           ir, filnam, h123, reac, crc, mp, mt, dpp,
     .           rc1min, rc1max, jfex1mn, jfex1mx, fp1,
     .           rc2min, rc2max, jfex2mn, jfex2mx, fp2,
     .           iz, elname, bundling,
     .           m_popesc, m_upper, m_lower,
     .           irow_esc, icol_esc, pop_esc,
     .           m_qext,
     .           irow_ext, ital_ext,
     .           iftflg, ncoef, coef,
     .           iprftype)

c  input: ZEILE:  character(420) (made upper case already)
c  output:  IR
c           FILNAM  (character(8), e.g.: AMJUEL, CR_HE, TAB2D, ADAS,...)
c           H123    (character(4), H.IJ,  IJ one or two integers, e.g. H.3, H.10,..)
c           REAC    (character(50), e.g.: 2.1.8a,...) 
c           CRC     (character(3),  e.g. CX, EI,...)
c           MP, MT, DPP,
cdr  ASYMPTOTICS, if available
cdr  ASYMPTOTICS, if available
c           RC1MIN, RC1MAX, RC2MIN, RC2MAX
c           JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
c           FP1, FP2
cdr  TAB2D (ADAS) additional specifications
c           IZ
c           ELNAME
c           BUNDLING
cdr  INTERNAL CR CODES...   additional specifications: 
cdr                         population escape factors,
cdr                         individual transition line identifiers
c           M_POPESC, M_UPPER, M_LOWER
c           IROW_ESC(:), ICOL_ESC(:), POPESC(:)
cdr  CONST OPTION: COEFFICINTES READ FROM INPUT FILE, NOT FROM DATABASE
c           IFTFLG, NCOEF, COEF(9)
cdr  PHOTONS  additional specifications
c           IPRFTYPE
c
c
c  Interpret (parsing) character(420) ZEILE for the call to SLREAC.f
c  and to fill data structure REACDAT(IR):
C  If necessary, further input cards are read here.

c  read_reaclines.f is called from block 4 input.f

c  In certain special cases, depending on FILNAM, one further input card is read
c   filnam=ADAS/TAB2D:  read IZ, ELNAME         (mandatory)
c   filnam=CR...     :  read irow,icol,pop_esc  (optional)
cdr filnam=PHOTON    :  read IPRFTYPE
c   filnam=CONST     :  read IFTFLG  (coded as FTxxxxxx, IFTFLG=xxxxxx)
c                       read fit coefficients here and hand them back for slreac
c                           if mod(iftflg,100).eq.10
c                           then
c                                only one coefficient is read
c                                ==>  ncoef=1, coef(1)  is read
c                           else
c                                up to 9 coefficients read
c                                ==>  ncoef<=9, coef(1:9)   are read
c                           better: clearly separate iftflg and reac in case of filnam=CONST
c  if R1MN is present     :  read JFEX1MN,(FP1(I),I=1,3)
c  if R1MX is present     :  read JFEX1MX,(FP1(I),I=4,6)
c  if R2MN is present     :  read JFEX2MN,(FP2(I),I=1,3)
c  if R2MX is present     :  read JFEX2MX,(FP2(I),I=4,6)

C
C  READ ONE REACTION FROM FILE "FILNAM" AT A TIME. Input card is on "ZEILE"
C
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comxs, only : REACLINES, IRLINES
#ifndef NO_JSON
      use eirmod_json, only : S_STACK, STRING_STACK, CRS_STACK,
     >                        EIRENE_PUSH_STRING_STACK
!cym/cpg provide access to json_ck type
      use json_module           !IGNORE
     >  , only : json_ck
!end cym/cpg
#endif
      implicit none
      character(420), intent(inout) :: zeile
#ifndef NO_JSON
!cym/cpg avoid type mismatch
      character(kind=json_CK,len=420) :: zeileCK
!cym/cpg end
#endif
      integer, intent(inout) :: iread
      integer, intent(in) :: iunin, iunout
      integer, intent(out) :: ir, mp, mt,
     .                        jfex1mn, jfex1mx, jfex2mn, jfex2mx,
     .                        iz, iftflg, ncoef
      real(dp), intent(out) :: dpp, rc1min, rc1max, rc2min, rc2max
      real(dp), intent(out) :: fp1(6), fp2(6), coef(9)
      character(8), intent(out) :: filnam
      character(4), intent(out) :: h123
      character(50), intent(out) :: reac
      character(3), intent(out) :: crc
      character(2), intent(out) :: elname
      character(60), intent(out) :: bundling
#ifndef NO_JSON
      character(:), allocatable :: cline
!cym/cpg avoid type mismatch
      character(kind=json_CK,len=420) :: clineCK
!cym/cpg end
#endif
cdr  additionally read, in case filnam = CR...
      integer, intent(out) :: m_popesc, m_upper, m_lower
      integer, allocatable, intent(inout) :: irow_esc(:), icol_esc(:) 
      real(dp), allocatable, intent(inout) :: pop_esc(:)
cdr  additionally read, in case filnam = CR...
      integer, intent(out) :: m_qext
      integer, allocatable, intent(inout) :: irow_ext(:), ital_ext(:,:) 
cdr  additionaly read, in case filnam = PHOTON
      integer, intent(out) :: iprftype

      integer :: iend, itok, isw, i1, i, il, ier, ipe
      real(dp) :: r1mn, r1mx, r2mn, r2mx
      real(dp) :: rhelp(12)
      character(72) :: sline
      character(12) :: chr
      character(3), save :: rea_proc(9) =
     .   (/ 'EI ','DS ','CX ','II ','PI ','EL ','RC ','PH ','OT ' /)
      logical :: chk_h123
      external :: eirene_lowercase, eirene_read_token, eirene_leer,
     .            eirene_exit_own

      filnam = repeat(' ',8)
      h123 = repeat(' ',4)
      reac = repeat(' ',50)
      crc = '   '
      elname = '  '
      iz = 0
cdr  additional (optional) data for internal CR-codes
      m_popesc = 0
      m_qext   = 0
cdr  m_upper, m_lower: currently not needed. tbd: To be used to
cdr                    turn population coefficient into an individual line emissivity
cdr  icol and irow are now vectors
      m_upper = 0
      m_lower = 0
cdr  additional data for filnam: PHOTON
      iprftype = 0
      ncoef = 0
      coef = 0._dp                                     
      iftflg = 0
      mp = 0
      mt = 0
      dpp = 0._dp
      fp1 = 0._dp
      fp2 = 0._dp
      rc1min = -huge(1._dp)
      rc1max =  huge(1._dp)
      rc2min = -huge(1._dp)
      rc2max =  huge(1._dp)
      jfex1mn = 0
      jfex1mx = 0
      jfex2mn = 0
      jfex2mx = 0

      iread = 0

c  1) leading integer IR:
      IEND = 1
      CALL EIRENE_READ_TOKEN (ZEILE(IEND:),' ',CHR,ITOK,IER,.FALSE.)
      IEND = IEND + ITOK
      READ (CHR,*) IR
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING REACTION NUMBER ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

      if ((ir < 1) .or. (ir > nreac)) then
        WRITE (iunout,*) ' ERROR REACTION NUMBER IR IS OUT OF RANGE '
        WRITE (iunout,*) ' 1 <= IR <= NREAC IS EXPECTED'
        WRITE (iunout,*) ' IR = ',IR, ' NREAC = ', NREAC
        CALL EIRENE_EXIT_OWN(1)
      end if

c  2) filename
      CALL EIRENE_READ_TOKEN (ZEILE(IEND:),' ',FILNAM,ITOK,IER,.FALSE.)
      IF (IER > 0) THEN
        WRITE (iunout,*)
     .        ' FILENAME STRING FOR REACTION ',IR,' TOO LONG '
        CALL EIRENE_EXIT_OWN(1)
      END IF
      IEND = IEND + ITOK

      IEND = IEND + VERIFY(ZEILE(IEND:),' ') - 1
      READ (ZEILE(IEND:IEND+3),'(A4)') H123
      IEND = IEND + 4

      CHK_H123 = (H123(1:2) == 'H.') .OR. (H123(1:2) == 'P.')
      IF (H123(4:4) == ' ') THEN
        CHK_H123 = CHK_H123 .AND. (VERIFY(H123(3:3),'1234567890') == 0)
      ELSE
        CHK_H123 = CHK_H123 .AND. (VERIFY(H123(3:4),'1234567890') == 0)
        IF (CHK_H123) THEN
          READ (H123(3:4),'(I2)') ISW
        END IF
        CHK_H123 = CHK_H123 .AND. (ISW <= 12)
      END IF
      IF (.NOT.CHK_H123) THEN
        WRITE (iunout,*) ' ERROR IN REACTION NUMBER ', IR
        WRITE (iunout,*) ' WRONG IDENTIFICATION FLAG FOR TYPE OF DATA'
        WRITE (iunout,*) ' H123 = ', H123
        CALL EIRENE_EXIT_OWN(1)
      END IF

c  3) reaction identifier
      IF (INDEX(FILNAM,'CONST') == 0) THEN !  INPUT FROM EXTERNAL
                                           !  A&M DATA FILE

C  THE INPUT FLAG  "FT...." IS NOT AVAILABLE HERE
C  IT MIGHT BE READ LATER FROM A&M DATA FILE AMJUEL, IN SUBR. SLREAC
C  READ INPUT FLAG "REAC", UP TO 50 CHARACTERS ALLOWED.
C  PUT THIS FLAG ON REAC.
        CALL EIRENE_READ_TOKEN (ZEILE(IEND:),' ',REAC,ITOK,IER,.FALSE.)
        IF (IER > 0) THEN
          WRITE (iunout,*)
     .        ' REACTION STRING FOR REACTION ',IR,' TOO LONG'
          CALL EIRENE_EXIT_OWN(1)
        END IF
C  NEXT: FIND POSITION FROM WHICH NEXT INPUT FLAG "CRC" CAN BE READ
        IEND = IEND + ITOK

      ELSE !  FILNAM=CONST OPTION

C  INPUT FOR A CONSTANT CROSS-SECTION, REACTION RATE, ETC.
C      OR FOR A HARD-WIRED POLYNOMIAL FIT. FIT COEFFICIENTS WILL BE READ
C      FROM CURRENT INPUT FILE, STREAM IUNIN,
C      NOT FROM EXTERNAL A&M DATA FILE
        REAC=REPEAT(' ',LEN(REAC))
!  CHECK CURRENT INPUT LINE FOR FT-FLAG (IFTFLG)
        IFTFLG = 0
        ITOK = INDEX(ZEILE(IEND:),'FT')
        IF (ITOK /= 0) THEN
C  STRING "FT" FOUND. READ 9 CHARACTERS INCLUDING THE STARTING "FT"
C                     PUT THESE ON REAC. THE ORIGINAL PARAMETER REAC
C                     IS NOT USED IN CASE OF CONST - OPTION
!  FIND START AND END POSITION OF FITFLAG
          I1 = ITOK + SCAN(ZEILE(IEND+ITOK+2:),'1234567890') - 1
          IF (I1 > 0) READ (ZEILE(IEND+I1:),*) IFTFLG
          IEND = IEND + ITOK + 2
C  NEXT: FIND POSITION FROM WHICH NEXT INPUT FLAG "CRC" CAN BE READ
          CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.FALSE.)
          IEND = IEND + ITOK
        END IF
!  READ FIT COEFFICIENTS
        IF (MOD(IFTFLG,100) == 10) THEN
          NCOEF = 1
          READ (IUNIN,'(6E12.4)') COEF(1)
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=12) :: CLINE)
          WRITE (CLINE,'(6E12.4)') COEF(1)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=trim(cline)
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        ELSE
          NCOEF = 9
          READ (IUNIN,'(6E12.4)') COEF(1:NCOEF)
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=72) :: CLINE)
          WRITE (CLINE,'(6E12.4)') COEF(1:6)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=trim(cline)
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          CLINE=REPEAT(' ',72)
          WRITE (CLINE,'(6E12.4)') COEF(7:9)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=trim(cline)
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        END IF
      END IF
C

c  3) flag for the type of the reaction process
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CRC,ITOK,IER,.FALSE.)
      IEND = IEND + ITOK
      IF (IER > 0) THEN
        WRITE (iunout,*) ' CRC-STRING FOR REACTION ',IR,' TOO LONG'
        CALL EIRENE_EXIT_OWN(1)
      END IF
      IF (ALL(REA_PROC /= CRC)) THEN
        WRITE (iunout,*) ' WRONG CRC-STRING FOR REACTION ',IR,
     .                   ' FOUND '
        CALL EIRENE_EXIT_OWN(1)
      END IF

C  THE REST OF INPUT DATA FROM THIS REACTION CARD IS NOW ON 'CHR'
C  these are flags to modify, scale, (or extrapolate) the input data tables or fits
C  READ FLAGS MP, MT, DPP, R1MN, R1MX, R2MN, R2MX FROM CHR

!  READ MP,   projectile mass
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.FALSE.)
      IEND = IEND + ITOK
      READ (CHR,*) MP
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING MP FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ MT,   target mass
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.FALSE.)
      IEND = IEND + ITOK
      READ (CHR,*) MT
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING MT FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ DPP,  potential energy increment in energy-weighted rate
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.TRUE.)
      IEND = IEND + ITOK
      READ (CHR,'(E12.4)') DPP
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING DPP FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ R1MN, lower boundary for extrapolation, 1st parameter
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.TRUE.)
      IEND = IEND + ITOK
      READ (CHR,'(E12.4)') R1MN
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING R1MN FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ R1MX, upper boundary for extrapolation, 1st parameter
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.TRUE.)
      IEND = IEND + ITOK
      READ (CHR,'(E12.4)') R1MX
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING R1MX FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ R2MN, lower boundary for extrapolation, 2nd parameter
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.TRUE.)
      IEND = IEND + ITOK
      READ (CHR,'(E12.4)') R2MN
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING R2MN FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

!  READ R2MX, upper boundary for extrapolation, 2nd parameter
      CALL EIRENE_READ_TOKEN(ZEILE(IEND:),' ',CHR,ITOK,IER,.TRUE.)
      IEND = IEND + ITOK
      READ (CHR,'(E12.4)') R2MX
      IF (IER > 0) THEN
        WRITE (iunout,*) ' ERROR READING R2MX FOR REACTION ',IR
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IF (INDEX(ZEILE,'ADAS') .NE. 0 .OR.
     .    INDEX(ZEILE,'TAB2D') .NE. 0 ) THEN
!  CHECK FOR ELEMENT NAME AND CHARGE STATE IN TABLE
        READ (IUNIN,'(4X,A2,1X,I3,2X,A60)') ELNAME,IZ,BUNDLING
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
        ALLOCATE (CHARACTER(LEN=72) :: CLINE)
        WRITE (CLINE,'(4X,A2,1X,I3,2X,A60)') ELNAME,IZ,BUNDLING
!cym/cpg avoid type mismatch
!        call eirene_push_string_stack(crs_stack,trim(cline))
        clineCK=trim(cline)
        call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
        DEALLOCATE(CLINE)
#endif
        CALL EIRENE_LOWERCASE(ELNAME)
        CALL EIRENE_LOWERCASE(BUNDLING)
      ELSE
        BUNDLING = REPEAT(' ',60)
        ELNAME = '  '
        IZ = 0
      END IF

      IF (INDEX(ZEILE,'CR') .NE. 0) THEN

cdr     optinal: M_POPESC:  how many pop esc factor to be read....
cdr              M_UPPER, M_LOWER: specify a particular line transition.
cdr              If a population coefficient is obtained from CR code,
cdr              find also Aik, E_upper, E_lower, g_upper, g_lower from CR code. 
cdr              I.e. turn a population coefficient POP into a 
cdr              line emission rate coefficient POP*Aik.
        READ (IUNIN,'(A80)') ZEILE
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
!cym/cpg avoid type mismatch
!        call eirene_push_string_stack(crs_stack,trim(zeile))
         zeileCK=zeile
         call eirene_push_string_stack(crs_stack,trim(zeileCK))
!cym/cpg end
#endif
cdr     CALL EIRENE_UPPERCASE (ZEILE), removed, because leading blank removal wrecks format
! CHECK FOR POPULATION ESCAPE FACTORS
! IF ANY OTHER CHARACTER (NOT A PURE REAL OR INTEGER) IS FOUND,
! VALUE OF VERIFY... GIVES THE FIRST (LEFTMOST) POSITION.
        IF (VERIFY(ZEILE,'+-.edED0123456789 ') > 0) THEN
! "ZEILE" CONTAINS NEXT REACTION LINE
          WRITE (IUNOUT,*)  'No POP_ESC, Q_EXT for internal CR model'
          WRITE (IUNOUT,*)  'REACTION',IR,'TYPE ',H123
          IREAD = 1
          M_POPESC = 0
          M_UPPER = 0
          M_LOWER = 0
        ELSE
! only integer or reals found in string ZEILE
! READ POPULATION ESCAPE FACTOR, FOR TRANSITION "UPPER=IROW --> LOWER=ICOL"

          READ (ZEILE,66666) M_POPESC, M_QEXT, M_UPPER, M_LOWER
          IF (M_POPESC.GT.0) THEN
            IF (ALLOCATED(IROW_ESC)) THEN
              IF (SIZE(IROW_ESC) < M_POPESC) THEN
                DEALLOCATE (IROW_ESC)
                DEALLOCATE (ICOL_ESC)
                DEALLOCATE (POP_ESC)
              END IF
            END IF
            IF (.NOT.ALLOCATED(IROW_ESC)) THEN
              ALLOCATE (IROW_ESC(M_POPESC), SOURCE = 0)
              ALLOCATE (ICOL_ESC(M_POPESC), SOURCE = 0)
              ALLOCATE (POP_ESC(M_POPESC), SOURCE = 1._DP)              
            END IF
          ELSE
cdr  found: m_popesc=0
          ENDIF
          DO IPE=1,M_POPESC
            READ (IUNIN,66665) 
     .            IROW_ESC(IPE), ICOL_ESC(IPE), POP_ESC(IPE)
            IF (IPE .EQ. 1) then
              WRITE (IUNOUT,*)  'POP_ESC found for internal CR model'
              WRITE (IUNOUT,*)  'REACTION',IR,'TYPE ',H123
            ENDIF
            WRITE (IUNOUT,'(A11,I3,A3,I3,A8,1E12.4)')
     .                    ' TRANSITION ',IROW_ESC(IPE),'-->',
     .                      ICOL_ESC(IPE),
     .                    ' POP_ESC',     POP_ESC(IPE)
          ENDDO
          IREAD = 0

C  POP. ESC. FACTORS ARE ONLY AVAILABLE FOR INTERNAL CR MODELS

C  now do the same for external source contributions
          IF (M_QEXT.GT.0) THEN
            IF (ALLOCATED(IROW_EXT)) THEN
              IF (SIZE(IROW_EXT) < M_QEXT) THEN
                DEALLOCATE (IROW_EXT)
                DEALLOCATE (ITAL_EXT)
              END IF
            END IF
            IF (.NOT.ALLOCATED(IROW_EXT)) THEN
              ALLOCATE (IROW_EXT(M_QEXT), SOURCE = 0)
              ALLOCATE (ITAL_EXT(3,M_QEXT), SOURCE = 0)
            END IF
          ELSE
cdr found: m_qext=0
          ENDIF
          DO IPE=1,M_QEXT
            READ (IUNIN,66666)
     .            IROW_EXT(IPE), ITAL_EXT(1,IPE), ITAL_EXT(2,IPE),
     .            ITAL_EXT(3,IPE)
            IF (IPE .EQ. 1) then
              WRITE (IUNOUT,*)  'Q_EXT tally for internal CR model'
              WRITE (IUNOUT,*)  'REACTION',IR,'TYPE ',H123
            ENDIF
            WRITE (IUNOUT,'(A10,I3,A11,2(I3,1X))')
     .           ' to LEVEL ',IROW_EXT(IPE),
     .           ' use TALLY ',ITAL_EXT(1,IPE), ITAL_EXT(2,IPE), 
     .           ITAL_EXT(3,IPE)
          ENDDO
          IREAD = 0

C  Q_EXT  SOURCES ARE ONLY AVAILABLE FOR INTERNAL CR MODELS

        END IF
C  POP. ESC. FACTORS ARE ONLY AVAILABLE FOR INTERNAL CRM MODELS
      END IF

      IF (INDEX(ZEILE,'PHOTON') .NE. 0) THEN
c       WRITE (IUNOUT,*)  'read_reaclines DATA FILE: PHOTON'
c       WRITE (iunout,*)  'READ EXTRA LINE (LINEPROFILE FLAG) HERE'
        read (iunin,'(12i6)') iprftype  ! april 22: was read_photdbk.f 
        WRITE (IUNOUT,*)  'Photon:',IR,'Line shape: iprftype ',iprftype       
      ENDIF

C  ASYMPTOTICS FOR CROSS-SECTIONS OR (WEIGHTED) RATE COEFFICIENTS FOR REACTION IR

C  OVERWRITES ASYMPTOTICS READ FROM EXTERNAL DATA FILES FOR THIS RUN,
C  (IF THERE HAVE BEEN SUCH)

cdr:  tbd: these asymptotics are different form the default asymp. defined
cdr        elsewhere.   There should be no such difference.
cdr   default asymp: means: no asympt. at all. fp1, fp2 never needed !


C  CURRENT DEFAULT ASYMPTOTICS: SET SUCH THAT "NO ASYMPTOTICS" IS USED
      IF (INDEX(H123,'P.').eq.0) then
c  single parametric data, or for first parameter in 2-parametric data
        RC1MIN = -20._DP ! lower ln(E), ln(T) default limit;
                         ! E,T in eV (2E-9 EV)
        RC1MAX =  20._DP ! upper ln(E), ln(T) default limit;
                         ! E,T in eV (5E8  EV)
c  2nd parameter in 2 parametric data
        RC2MIN = -20._DP ! lower ln(E0), ln(N) default limit;
                         ! E0 in eV, N in cm**-3
        RC2MAX = 100._DP ! upper ln(E0), ln(N) default limit;
                         ! E0 in eV, N in cm**-3

C  ARE PARAMETERS FOR ASYMPTOTICS FOR THIS REACTION IR SPECIFIED
C  EXPLICITLY IN input file?
C  IF YES: OVERWRITE DEFAULTS, AND/OR DATA FROM EXTERNAL FILE
C  PARAMETERS ARE E,T,N: ALWAYS POSITIVE
        IF (R1MN.GT.0.D0) THEN
          READ (IUNIN,66664) JFEX1MN,(FP1(I),I=1,3)
          RC1MIN=LOG(R1MN)
          WRITE (IUNOUT,*) 'NON-DEF. R1MN FOR REACTION IR=',IR,R1MN
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=72) :: CLINE)
          WRITE (CLINE,66664) JFEX1MN,(FP1(I),I=1,3)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=cline
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        ENDIF
        IF (R1MX.GT.0.D0) THEN
          READ (IUNIN,66664) JFEX1MX,(FP1(I),I=4,6)
          RC1MAX=LOG(R1MX)
          WRITE (IUNOUT,*) 'NON-DEF. R1MX FOR REACTION IR=',IR,R1MX
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=72) :: CLINE)
          WRITE (CLINE,66664) JFEX1MX,(FP1(I),I=4,6)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=cline
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        ENDIF
        IF (R2MN.GT.0.D0) THEN
          READ (IUNIN,66664) JFEX2MN,(FP2(I),I=1,3)
          RC2MIN=LOG(R2MN)
          WRITE (IUNOUT,*) 'NON-DEF. R2MN FOR REACTION IR=',IR,R2MN
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=72) :: CLINE)
          WRITE (CLINE,66664) JFEX2MN,(FP2(I),I=1,3)
!cym/cpg avoid type mismatch
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=cline
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        ENDIF
        IF (R2MX.GT.0.D0) THEN
          READ (IUNIN,66664) JFEX2MX,(FP2(I),I=4,6)
          RC2MAX=LOG(R2MX)
          WRITE (IUNOUT,*) 'NON-DEF. R2MX FOR REACTION IR=',IR,R2MX
#ifndef NO_JSON
!  STORE REACTION LINES FOR OUTPUTING TO JSON FILE
          ALLOCATE (CHARACTER(LEN=72) :: CLINE)
          WRITE (CLINE,66664) JFEX2MX,(FP2(I),I=4,6)
!cym/cpg end
!          call eirene_push_string_stack(crs_stack,trim(cline))
          clineCK=cline
          call eirene_push_string_stack(crs_stack,trim(clineCK))
!cym/cpg end
          DEALLOCATE(CLINE)
#endif
        ENDIF
        IF (MAX(R1MN,R1MX,R2MN,R2MX).GT.0.D0) CALL EIRENE_LEER(1)

      else
! identifier "P" found in H123. Data for photon processes! No asymptotics available
! defaults already set at beginning of subroutine
      endif

C  ASYMPTOTICS: DONE
C
C  SAVE INPUT LINES ON REACLINES, FOR WRITING ONTO JSON-FILE
C       ONLY NEEDED FOR AUTOMATED INTERFACE TO HYDKIN DATABASE.
Crc  not used at the moment, need revision as IL not used here but in read_json
      IL = IRLINES+1
      REACLINES(IL)%NO = IR
      REACLINES(IL)%FILE = FILNAM
      REACLINES(IL)%H_SELECT = H123
      REACLINES(IL)%REAC_STRING = REAC
      REACLINES(IL)%REACTYP = CRC
      REACLINES(IL)%MP = MP
      REACLINES(IL)%MT = MT
      REACLINES(IL)%DPP = DPP
      REACLINES(IL)%R1MN = R1MN
      REACLINES(IL)%R1MX = R1MX
      REACLINES(IL)%R2MN = R2MN
      REACLINES(IL)%R2MX = R2MX
      REACLINES(IL)%ELEMENT = ELNAME
      REACLINES(IL)%BUNDLING = BUNDLING
      REACLINES(IL)%IZ = IZ
      REACLINES(IL)%JFEX1MN = 0
      REACLINES(IL)%JFEX1MX = 0
      REACLINES(IL)%JFEX2MN = 0
      REACLINES(IL)%JFEX2MX = 0
      REACLINES(IL)%FP1 = 0._DP
      REACLINES(IL)%FP2 = 0._DP
      REACLINES(IL)%M_POPESC = M_POPESC
 	    IF (M_POPESC >0) THEN
        ALLOCATE(REACLINES(IL)%IROW_ESC(M_POPESC))
        ALLOCATE(REACLINES(IL)%ICOL_ESC(M_POPESC))
        ALLOCATE(REACLINES(IL)%POP_ESC(M_POPESC))
         REACLINES(IL)%IROW_ESC(1:M_POPESC) = IROW_ESC(1:M_POPESC)
         REACLINES(IL)%ICOL_ESC(1:M_POPESC) = ICOL_ESC(1:M_POPESC)
         REACLINES(IL)%POP_ESC(1:M_POPESC) = POP_ESC(1:M_POPESC)
      END IF
      REACLINES(IL)%M_QEXT = M_QEXT
      IF (M_QEXT >0) THEN
        ALLOCATE(REACLINES(IL)%IROW_EXT(M_QEXT))
        ALLOCATE(REACLINES(IL)%ITAL_EXT(3,M_QEXT))
         REACLINES(IL)%IROW_EXT(1:M_QEXT) = IROW_EXT(1:M_QEXT)
         REACLINES(IL)%ITAL_EXT(1:3,1:M_QEXT) = ITAL_EXT(1:3,1:M_QEXT)
      END IF
      REACLINES(IL)%IFTFLG = IFTFLG
      REACLINES(IL)%NCOEF = NCOEF
      IF (NCOEF > 0) REACLINES(IL)%COEF(1:NCOEF) = COEF(1:NCOEF)
      IRLINES = IL
C  DONE

66664 FORMAT (I6,6X,5E12.4)
66665 FORMAT (2I6,4E12.4)
66666 FORMAT (4I6)
      return
      end subroutine eirene_read_reaclines
