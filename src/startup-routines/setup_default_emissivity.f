cdr  Oct. 18: further comments:

cdr           This routine implicity makes some assumptions regarding the
cdr           species in input block 4a,b,c,d:
cdr           H2 (type=2)
cdr           H+ (type=4)
cdr           H  (type=1)
cdr   to be checked: contributions? multiple isotopes ?

      subroutine eirene_setup_default_emissivity

cdr originally programmed by PB 2017
cdr Routine is called from subr. INPUT.f, input block 12.
cdr april 18:  the calculation of volumetric line emissitivies
cdr            and their storing on additional tallies ADDV
cdr            has been generalized,
cdr            replacing the former 6 routines:
cdr            ba_alpha.f, ba_beta.f, ba_gamma.f, ba_delta.f,
cdr            ly_alpha.f, ly_beta.f

cdr  This present routine:
cdr  Try to reproduce the content of the old versions of these 6 routines,
cdr  by filling (and later using)
cdr  the new structures EMIS_LINES%....
cdr
cdr  number of lines       6 (BA_AL, BA_BET, ..., LY_BET)
cdr  number of components: 6 (COUPLING TO H,H+,H2,H2+,H-,H3+)
cdr  number of contributions: detected from input file,
cdr                           as in old ba... ly... routines
cdr                          (there sum over contributions only
cdr                           on ADDV tallies),
cdr  hard-coded here: use pop. coeffs from amjuel H.12, and
cdr                   use ratios for short living radicals (H2+, H3+, H-)
cdr                   from amjuel H.11 and H.12
cdr  hard-coded:
cdr              H2+, H- and H3+ QSS states, because of hard-coded density ratios.
cdr              H2+ must be produced from both EI and IC processes, because
cdr              hard-coded ratio H.12 2.0c is used here.
cdr
cdr  tbd: make consistent notation "component vs. contribution": DONE !
cdr
cdr  The ADDV tallies are filled later,
cdr  in calls to emissivity.f from MCARLO, per stratum.

cdr  Apparently we need at least one chord and nchtal=2, (even if unused)
cdr  to fill the ADDV arrays, because of a hidden link between emissivity options
cdr  and line-of-sight options.
cdr

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_comusr
      use eirmod_comxs, only : nreaci

      implicit none

      TYPE(TCONTRIB) :: CNT

      INTERFACE
        SUBROUTINE EIRENE_SLREAC (IR, FILNAM, H123, REAC, CRC,
     .             RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .             RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
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
      END INTERFACE

      integer :: i, NUM_compo, iat, iml, ipl, nat, npl, nml,
     .           nrc, nrc_rat1, nrc_rat2, nrc_rat3
      integer :: mp, mt,
     .           jfex1mn, jfex1mx, jfex2mn, jfex2mx,
     .           iz, irow_esc, icol_esc, iftfl, ncoef
      real(dp) :: dpp, rc1min, rc1max, rc2min, rc2max, pop_esc
      real(dp) :: fp1(6), fp2(6), coef(9)
      character(8) :: filnam
      character(4) :: h123
      character(50) :: reac
      character(3) :: crc
      character(2) :: elname
      character(60) :: bundling

      NUM_lines    = 6  ! Ba_alpha, Ba_beta, Ba_gamma, Ba_delta,
                        ! Ly_alpha, Ly_beta
      NUM_compo    = 6  ! coupling to H, H+, H2, H2+, H- and H3+
c     NUM_contrib  = inferred from input file, species specification block 4a,b,c.
      MOD_ADDV = 0
      NRC = NREACI

!  defaults for reactions which are to be added for default emissivity model
      elname = '  '
      iz = 0
      bundling = repeat(' ',60)
      irow_esc = 0
      icol_esc = 0
      pop_esc = 1._dp
      ncoef = 0
      coef = 0._dp
      iftfl = 0
      mp = 0
      mt = 0
      dpp = 0._dp
c default asymptotics
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

! read data for ratios from AMJUEL,
c ratio H2+/H2, and assuming also an IC contribution to H2+ formation
      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.0c     '
      CRC    = 'OT '
      NRC = NRC + 1
      NRC_RAT1 = NRC
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)

c ratio H-/H2
      FILNAM = 'AMJUEL  '
      H123   = 'H.11'
      REAC   = '7.0a     '
      CRC    = 'OT '
      NRC = NRC + 1
      NRC_RAT2 = NRC
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)

c ratio H3+/H2:  prod rate: [H2+] [H2].  loss rate: [ne] [H3+]
      FILNAM = 'AMJUEL  '
      H123   = 'H.11'
      REAC   = '4.0a     '
      CRC    = 'OT '
      NRC = NRC + 1
      NRC_RAT3 = NRC
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)

      ALLOCATE (EMIS_LINES(NUM_LINES))
      EMIS_LINES%LINE_NAME = REPEAT(' ',80)
      EMIS_LINES%NUM_COMPO = 0

************************************************
* BALMER ALPHA, LINE NO. 1, all 6 components
************************************************

      EMIS_LINES(1)%LINE_NAME = 'BA_ALPHA'
      EMIS_LINES(1)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(1)%EINSTEIN = 4.410E7
C  transition energy
      EMIS_LINES(1)%TRANS_EN = 1.8889_DP
C  identifier of line:
      EMIS_LINES(1)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(1)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=3)/H(n=1)

      EMIS_LINES(1)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(1)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(1)%IRC = NRC

cdr  hard-coded type: atoms. Species: default: =0
      CNT%ISP     = 0
      CNT%ITP     = 1
c
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(1)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(1)%COMPO(1)%NUM_CONTRIB = NAT

cdr all reaction data are the same for all contributions.
cdr only CNT%ISP (species index) may differ for different contributions.
      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+/D+/T+ ION DENSITY
C  ALL BULK ION (ITYP=4) CONTRIBUTIONS WITH
C                        NUCLEAR CHARGE NUMBER=1 AND CHARGE STATE NUMBER=1
C  H(n=3)/H+

      EMIS_LINES(1)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(1)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(2)%IRC = NRC

cdr  hard-coded type: bulk ions. Species: default: =0
      CNT%ISP     = 0
      CNT%ITP     = 4
c
      CNT%IRATIO  = 0
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(1)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(1)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN "H2" MOLEC. DENSITY
C  ALL MOLECULE (ITYP=2) CONTRIBUTIONS WITH
C                        NUCLEAR CHARGE NUMBER=2
C  H(n=3)/H2(g)

      EMIS_LINES(1)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(1)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(3)%IRC = NRC

cdr hard-coded type: molecules. Species: default: =0
      CNT%ISP     = 0
      CNT%ITP     = 2

      CNT%IRATIO  = 0
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN "H2+" MOLEC. ION DENSITY
C  H(n=3)/H2+(g)

      EMIS_LINES(1)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(1)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14a  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(4)%IRC = NRC

cdr  hard-coded type: molecules. Species: default: =0
cdr  test ion (H2+) is in QSS (nfoli < 0) with molecule, ratio H.12 2.0c
cdr  i.e.: H2+ produced from H2, via both channels EI plus IC.
      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%IRC_RAT(1) = NRC_RAT1   !  density ratio H2+/H2

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(2) = -1

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=3)/H-

      EMIS_LINES(1)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(1)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2a     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(5)%IRC = NRC

c  hard-coded: type: molecules. Species: default: =0
cdr  test ion (H-) is in QSS (nfoli < 0) with molecule, ratio H.11 7.0a
cdr  i.e.: H- produced from H2, via dissociative attachment.
      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%IRC_RAT(1) = NRC_RAT2   ! density ratio H-/H2

c  no second density ratio
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(2) = -1

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=3)/H3+

      EMIS_LINES(1)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(1)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15a  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(1)%COMPO(6)%IRC = NRC

c  hard-coded: type: molecules. Species: default: =0
cdr  test ion (H3+) is in QSS (nfoli < 0) with molecule,
cdr  ratio1  H.12 7.0c and  ratio2 H.11 4.0a
cdr  i.e.: H3+ produced from H2 and H2+, via particle exchange channels
      CNT%ISP     = 0
      CNT%ITP     = 2

      CNT%IRATIO          = 2
      CNT%IRC_RAT(1) = NRC_RAT1   ! density ratio H2+/H2

c  two species, [nH2/ne] density ratio to be multiplied
c  H2
      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
c  electron density
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5

      CNT%IRC_RAT(2) = NRC_RAT3   ! ratio loss H3+/ prod H3+

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(1)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(1)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(1)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO


************************************************
* BALMER BETA, LINE NO. 2
************************************************

      EMIS_LINES(2)%LINE_NAME = 'BA_BETA'
      EMIS_LINES(2)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(2)%EINSTEIN = 8.419E6_DP
      EMIS_LINES(2)%TRANS_EN = 2.5500_DP
      EMIS_LINES(2)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(2)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=4)/H(n=1)

      EMIS_LINES(2)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(2)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5c   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(1)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 1
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(2)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(2)%COMPO(1)%NUM_CONTRIB = NAT

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=4)/H+

      EMIS_LINES(2)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(2)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8c   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(2)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 4
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(2)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(2)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=4)/H2(g)

      EMIS_LINES(2)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(2)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5c   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(3)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=4)/H2+(g)

      EMIS_LINES(2)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(2)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14c  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(4)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT1
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=4)/H-

      EMIS_LINES(2)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(2)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2c     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(5)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT2
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=4)/H3+

      EMIS_LINES(2)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(2)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15c  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(2)%COMPO(6)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2

      CNT%IRATIO  = 2
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5
      CNT%IRC_RAT(2) = NRC_RAT3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(2)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(2)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
         IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(2)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO

************************************************
* BALMER GAMMA
************************************************

      EMIS_LINES(3)%LINE_NAME = 'BA_GAMMA'
      EMIS_LINES(3)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(3)%EINSTEIN = 2.530E6_DP
      EMIS_LINES(3)%TRANS_EN = 2.8560_DP
      EMIS_LINES(3)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(3)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=5)/H(n=1)

      EMIS_LINES(3)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(3)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5d   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(1)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 1
      CNT%IRATIO  = 0
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(3)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(3)%COMPO(1)%NUM_CONTRIB = NAT

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=5)/H+

      EMIS_LINES(3)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(3)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8d   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(2)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 4
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(3)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(3)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=5)/H2(g)

      EMIS_LINES(3)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(3)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5d   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(3)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=5)/H2+(g)

      EMIS_LINES(3)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(3)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14d  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(4)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT1
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=5)/H-

      EMIS_LINES(3)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(3)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2d     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(5)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT2
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=5)/H3+

      EMIS_LINES(3)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(3)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15d  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(3)%COMPO(6)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 2
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5

      CNT%IRC_RAT(2) = NRC_RAT3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(3)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(3)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(3)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO



************************************************
* BALMER DELTA
************************************************

      EMIS_LINES(4)%LINE_NAME = 'BA_DELTA'
      EMIS_LINES(4)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(4)%EINSTEIN = 9.732E5_DP
      EMIS_LINES(4)%TRANS_EN = 3.0222_DP
      EMIS_LINES(4)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(4)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=6)/H(n=1)

      EMIS_LINES(4)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(4)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5e   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(1)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 1
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(4)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(4)%COMPO(1)%NUM_CONTRIB = NAT

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=6)/H+

      EMIS_LINES(4)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(4)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8e   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(2)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 4
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(4)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(4)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=6)/H2(g)

      EMIS_LINES(4)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(4)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5e   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(3)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=6)/H2+(g)

      EMIS_LINES(4)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(4)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14e  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(4)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT1
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=6)/H-

      EMIS_LINES(4)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(4)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2e     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(5)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT2
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=6)/H3+

      EMIS_LINES(4)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(4)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15e  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(4)%COMPO(6)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 2
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5
      CNT%IRC_RAT(2) = NRC_RAT3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(4)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(4)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(4)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO

************************************************
* LYMAN ALPHA
************************************************

      EMIS_LINES(5)%LINE_NAME = 'LY_ALPHA'
      EMIS_LINES(5)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(5)%EINSTEIN = 4.699E8_DP
      EMIS_LINES(5)%TRANS_EN = 10.2375_DP
      EMIS_LINES(5)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(5)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=2)/H(n=1)

      EMIS_LINES(5)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(5)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5b   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(1)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 1
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(5)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(5)%COMPO(1)%NUM_CONTRIB = NAT

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=2)/H+

      EMIS_LINES(5)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(5)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8b   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(2)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 4
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(5)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(5)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=2)/H2(g)

      EMIS_LINES(5)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(5)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5b   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(3)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=2)/H2+(g)

      EMIS_LINES(5)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(5)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14b  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(4)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT1
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=2)/H-

      EMIS_LINES(5)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(5)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2b     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(5)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(1) = NRC_RAT2
      CNT%IRC_RAT(2) = -1

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=2)/H3+

      EMIS_LINES(5)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(5)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15b  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(5)%COMPO(6)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 2
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5
      CNT%IRC_RAT(2) = NRC_RAT3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(5)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(5)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(5)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO

************************************************
* LYMAN BETA
************************************************

      EMIS_LINES(6)%LINE_NAME = 'LY_BETA'
      EMIS_LINES(6)%NUM_COMPO = NUM_COMPO
C  RADIATIVE TRANSITION RATE (1/S)
      EMIS_LINES(6)%EINSTEIN = 5.575E7_DP
      EMIS_LINES(6)%TRANS_EN = 12.089_DP
      EMIS_LINES(6)%IADV_TOTAL = NADVI + NUM_COMPO+1

      ALLOCATE (EMIS_LINES(6)%COMPO(NUM_COMPO))

C  COMPONENT 1: LINEAR IN H/D/T ATOM DENSITY
C  ALL TEST ATOM (ITYP=1) CONTRIBUTIONS WITH
C                         NUCLEAR CHARGE NUMBER=1
C  H(n=3)/H(n=1)

      EMIS_LINES(6)%COMPO(1)%COMPO_NAME = 'ATOMIC NEUTRAL HYDR.'
      EMIS_LINES(6)%COMPO(1)%IADV = NADVI + 1

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.5a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(1)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 1
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NAT = COUNT(NCHARA == 1)
      ALLOCATE (EMIS_LINES(6)%COMPO(1)%CONTRIB(NAT))
      EMIS_LINES(6)%COMPO(1)%NUM_CONTRIB = NAT

      IAT = 0
      DO I = 1, NATMI
        IF (NCHARA(I) == 1) THEN
          IAT = IAT + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(1)%CONTRIB(IAT) = CNT
        END IF
      END DO

C  COMPONENT 2: LINEAR IN H+ ION DENSITY
C  H(n=3)/H+

      EMIS_LINES(6)%COMPO(2)%COMPO_NAME = 'ATOMIC HYDR. ION'
      EMIS_LINES(6)%COMPO(2)%IADV = NADVI + 2

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.1.8a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(2)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 4
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

      NPL = COUNT((NCHARP == 1).and.(NCHRGP == 1))
      ALLOCATE (EMIS_LINES(6)%COMPO(2)%CONTRIB(NPL))
      EMIS_LINES(6)%COMPO(2)%NUM_CONTRIB = NPL

      IPL = 0
      DO I = 1, NPLSI
        IF ((NCHARP(I) == 1).and.(NCHRGP(I) == 1)) THEN
          IPL = IPL + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(2)%CONTRIB(IPL) = CNT
        END IF
      END DO

C  COMPONENT 3: LINEAR IN H2 MOLEC. DENSITY
C  H(n=3)/H2(g)

      EMIS_LINES(6)%COMPO(3)%COMPO_NAME = 'DIATOMIC NEUTRAL HYDR. MOL'
      EMIS_LINES(6)%COMPO(3)%IADV = NADVI + 3

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.5a   '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(3)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 0

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT = 0

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(3)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(3)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(3)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 4: LINEAR IN H2+ MOLEC. ION DENSITY
C  H(n=3)/H2+(g)

      EMIS_LINES(6)%COMPO(4)%COMPO_NAME =
     .     'DIATOMIC HYDR. MOL ION'
      EMIS_LINES(6)%COMPO(4)%IADV = NADVI + 4

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.14a  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(4)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(2) = -1

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(4)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(4)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(4)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 5: LINEAR IN H- NEG. ION DENSITY
C  H(n=3)/H-

      EMIS_LINES(6)%COMPO(5)%COMPO_NAME =
     .     'NEGATIVE HYDR. ION'
      EMIS_LINES(6)%COMPO(5)%IADV = NADVI + 5

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '7.2a     '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(5)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 1
      CNT%IRC_RAT(1) = NRC_RAT2

      CNT%ISP_RAT = -1
      CNT%ITP_RAT = -1
      CNT%IRC_RAT(2) = -1

cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(5)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(5)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(5)%CONTRIB(IML) = CNT
        END IF
      END DO

C  COMPONENT 6: LINEAR IN H3+ MOL. ION DENSITY
C  H(n=2)/H3+

      EMIS_LINES(6)%COMPO(6)%COMPO_NAME =
     .     'TRIATOMIC HYDR. ION'
      EMIS_LINES(6)%COMPO(6)%IADV = NADVI + 6

      FILNAM = 'AMJUEL  '
      H123   = 'H.12'
      REAC   = '2.2.15a  '
      CRC    = 'OT '
      NRC = NRC + 1
      CALL EIRENE_SLREAC (NRC,FILNAM,H123,REAC,CRC,
     .               RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX,
     .               RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX,
     .               ELNAME,IZ,BUNDLING)
      EMIS_LINES(6)%COMPO(6)%IRC = NRC

      CNT%ISP     = 0
      CNT%ITP     = 2
      CNT%IRATIO  = 2
      CNT%IRC_RAT(1) = NRC_RAT1

      CNT%ISP_RAT(1) = 1
      CNT%ITP_RAT(1) = 2
      CNT%ISP_RAT(2) = 1
      CNT%ITP_RAT(2) = 5
      CNT%IRC_RAT(2) = NRC_RAT3

      NML = COUNT(NCHARM == 2)
      ALLOCATE (EMIS_LINES(6)%COMPO(6)%CONTRIB(NML))
      EMIS_LINES(6)%COMPO(6)%NUM_CONTRIB = NML

      IML = 0
      DO I = 1, NMOLI
cdr search all neutral molecular species with nucl. charge number=2.
cdr These must be the isotopomers of H2.
        IF (NCHARM(I) == 2) THEN
          IML = IML + 1
          CNT%ISP = I
          EMIS_LINES(6)%COMPO(6)%CONTRIB(IML) = CNT
        END IF
      END DO

      RETURN

      end subroutine eirene_setup_default_emissivity
