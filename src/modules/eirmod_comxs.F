cdr Jan. 23: Remove obsolete FACREA
cdr          Except for REACDAT and IFTFLG, which are (-11:nreac),
cdr          all other input flag parameters are only (1:NREAC),
cdr          ie. KK or IR reaction indices do not cover minimal model reactions.
cdr          write_cmamf and read_cmamf also do not cover minimal model reactions
cdr Feb. 22: some missing photonic tally data added.
cdr          More "unified type pointers" added, to prepare for unification of collide.f
cdr dec. 21: some extra stuff is in this module for IFIT=5 (internal CR codes) option,
cdr          whereas for other IFIT options these corresponding parts
cdr          seem to be elsewhere?  cleanup needed?
cdr dec. 20: clean up nomenclature for reaction energetics flags: nelr.., nplr.., nhvr..
cdr Nov. 20: various data within TYPE.. Data constructs have been made
cdr          ALLOCATABLE, rather than POINTER, due to more recent Fortran capabilities
cdr          e.g.: ADAS_DATA (= TAB2D_DATA), POLY_DATA, TAB1D_DATA
cdr Feb  20: tbd: remove redundant nstor, nstor1 (unused),
cdr          removed: reac_name (unused)

cdr Aug  19: remove redundant data typ: HYDKIN, tbd: TAB1D, started
cdr Oct  18: tbd: separate quantities that vary along trajectories  (i.e. per thread)
cdr               from those that remain fixed after initialization (i.e. per node)
cdr          private per thread: sigv..., i.e.: xstor and xstorv
cdr          public per node: all the rest ?
cdr Sept.18: remove XDR format for stream fort.13
cdr May  18: FLDLM arrays (old fluid limit flags) now replaced by EDPOT arrays,
cdr          for potential energy difference in reactions.
cdr          The old fluid limit critical Knudsen number is now defined
cdr          via negative ngen..(..) flags
cdr Apr. 18: further pointer, targets set for photons, towards code synchronisation
cdr          across particle types, incl. photons
cdr Nov. 17: p2nds --> p2nei (now in full analogy with p2npi)
cdr
cdr  MXCOLLS --> MSTOR0

      MODULE EIRMOD_COMXS

cdr MODULE FOR ALL ATOMIC/MOLECULAR/PHOTONIC DATA STRUCTURES.

!  jan-05: natprc_2,..... introduced
!  07.12.05: bugfix: IFTFLG is now available for default reactions too
!                    via dimensioning IFTFLG(-11:NREAC,0:5)
!  30.08.06: data structure for reaction data redefined
!  12.10.06: modcol revised
!  19.12.06: test functions added which allow to test if a rate coefficient
!            is defined via ADAS database
!  02.03.07: remove ESCD2* arrays
!  02.03.07: fourth secondary group specifier introduced
!
cdr sometime between 2004 and 2007 the atomic data structure was revised.
cdr
cdr  now it is on REACDAT.  Commenting, cleanup started: jan 2016.
!
!  24.03.15: number of default reactions increased from 10 to 11, REACDAT(-11)...
cdr23.04.15: only text, comments.... continued: Nov. 15, still not complete

cdr  Jan  18:  added colrad_data, alloc_fit_form, rp%ifit=5 option: use internal CR code
cdr  sept 18:  prepare reviving "storage save mode" (for large 3D grids):
cdr            first: rationalize naming of integer flags for collision models
cdr            nhvrei, nhvrpi, for KER (heavy particle post-collision kinetics)
cdr            remove redundant flags: JEREARC  (UNUSED)
cdr            remove redundant flags: JEREAEI  (UNUSED)
cdr jun   22:  add splitting post-collision secondaries
cdr            at selected collision processes.
cdr            So far: tested for CX
cnh nov   23:  add macroscopic length scale L_MACRO for hybrid fluid-kinetic model

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMXS, EIRENE_DEALLOC_COMXS,
     .          EIRENE_INIT_CMDTA,  EIRENE_BROADCAST_COMXS,
     .          EIRENE_WRITE_CMDTA, EIRENE_READ_CMDTA,
     .          EIRENE_WRITE_CMAMF, EIRENE_READ_CMAMF,
     .          EIRENE_GET_REACTION, EIRENE_SET_REACTION_DATA,
     .          EIRENE_FREE_REACDAT,
cdr
     .          LINE_DATA, POLY_DATA, ADAS_DATA, TAB1D_DATA,
     .          COLRAD_DATA,
     .          REACTION_DATA,
     .          FIT_FORMS,
     .          REACTION_INPUT_LINE,
cdr
     .          EIRENE_IS_RTC_TAB2D,
     .          EIRENE_IS_RTCEW_TAB2D,
     .          EIRENE_IS_RTCMW_TAB2D,
c
     .          EIRENE_ALLOC_FIT_FORM

cdr Next: various specific data types for the 5 fit forms IFIT=1,2,...,5.

      TYPE LINE_DATA
cdr fix ph4, nothing allocatable in here any more
cdr also: distinct from poly, adas, colrad, tab1d, this is not a fit form
cdr but a bound-bound line shape parametrization.
        REAL(DP) :: E0, E1, AIK, G1, G2, C2, C3, C4, C6, B12, B21
        REAL(DP) :: C6A(12)
        INTEGER :: IGND, IRCART, IPROFILETYPE, IFREMD, NRJPRT, IMESS
        INTEGER :: IPLSC6(12)
        CHARACTER(2) :: KENN(12)
        CHARACTER(50) :: REACNAME
      END TYPE LINE_DATA

      TYPE ADAS_DATA  !dr: currently only H.4, H.10, vs. Temp, Dens
        INTEGER :: NDENS, NTEMP
        REAL(DP), POINTER :: DENS(:), TEMP(:), TAB2D(:,:)
        REAL(DP), POINTER :: DDE(:), DTE(:)
      END TYPE ADAS_DATA

!pb  IFEXMN, IFEXMX, RCMN, RCMX, FPARM removed from POLY_DATA
cdr  The extrapolation options are now made available generally,
cdr  for all types of A&M data input

      TYPE POLY_DATA
        REAL(DP), POINTER :: DBLPOL(:,:)
      END TYPE POLY_DATA

      TYPE TAB1D_DATA  !dr: currently: only H.2, vs. Temp
        INTEGER :: NTEMPS
        REAL(DP), POINTER :: TEMPS(:), TAB1D(:), DIFFQUOT(:)
        CHARACTER(50) :: REAC_STRING, REACNAME
        CHARACTER(100) :: RPRT
      END TYPE TAB1D_DATA

      TYPE COLRAD_DATA  !dr:  internal CR codes, eff. rates & pop.coefs
        INTEGER :: IFLAV, IVARST, IROW_ESC, ICOL_ESC
        REAL(DP) :: POP_ESC
      END TYPE COLRAD_DATA

      TYPE FIT_FORMS
        INTEGER :: IFIT
        TYPE(POLY_DATA),   POINTER :: POLY              ! IFIT = 1,2
        TYPE(ADAS_DATA),   POINTER :: ADAS  ! or: TAB2D ! IFIT = 3
cdr line_data is not a "fit form", but just some parameters that
cdr define a particular bound-bound line (here: "photon test particle").
cdr It should not be part of this data structure
        TYPE(LINE_DATA),   POINTER :: LINE
        TYPE(TAB1D_DATA),  POINTER :: TAB1D             ! IFIT = 4
        TYPE(COLRAD_DATA), POINTER :: CRM               ! IFIT = 5

cdr asymptotics, common to all of these 5 A&M data formats.
        REAL(DP) :: RC1MIN, RC1MAX, RC2MIN, RC2MAX
        REAL(DP) :: FP1L(3), FP1R(3), FP2B(3), FP2T(3)
        INTEGER  :: JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX
      END TYPE FIT_FORMS

      TYPE REACTION_DATA
        TYPE(FIT_FORMS), POINTER :: POT, CRS, RTC, RTCMW, RTCEW,
     T                              OTH, PHR

        LOGICAL :: LPOT, LCRS, LRTC, LRTCMW, LRTCEW, LOTH, LPHR
        REAL(DP) :: ETH, KER, DELP
        REAL(DP) :: RTMAX, ERTMAX, EARRH0, EARRH1
        INTEGER :: NOSEC
      END TYPE REACTION_DATA

      TYPE REACTION_INPUT_LINE
        INTEGER :: NO, MT, MP, IZ, JFEX1MN, JFEX1MX, NCONST,
     .             JFEX2MN, JFEX2MX,
     .             IROW_ESC, ICOL_ESC,
     .             IFTFLG, NCOEF
        REAL(DP) :: R1MN, R1MX, DPP, FP1(6), COEF(9),
     .              R2MN, R2MX,      FP2(6)
        REAL(DP) :: POP_ESC
        CHARACTER(8) :: FILE
        CHARACTER(50) :: REAC_STRING
        CHARACTER(4) :: H_SELECT
        CHARACTER(3) :: REACTYP
        CHARACTER(2) :: ELEMENT
        CHARACTER(60) :: BUNDLING
      END TYPE REACTION_INPUT_LINE

      TYPE(LINE_DATA), POINTER, PUBLIC, SAVE :: REACTION
      INTEGER, PUBLIC, SAVE :: IDREAC,  ! save kk value of previous
                                        ! reacdat(kk) evaluation
     .                         IRLINES

      TYPE(REACTION_DATA), ALLOCATABLE, PUBLIC, SAVE :: REACDAT(:)

      TYPE(REACTION_INPUT_LINE), ALLOCATABLE, PUBLIC, SAVE ::
     .                           REACLINES(:)
cdr....................................................................

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        XSTOR(:,:), XSTORV(:)

cdr  local (on the fly) atomic-molecular reaction data
      REAL(DP), PUBLIC, POINTER, SAVE ::
c  reaction rates, by reaction
     R SIGVCX(:),   SIGVPI(:),   SIGVEI(:),   SIGVEL(:),   SIGVPH(:),
c  energy exchange rates, by reaction
     R ESIGCX(:,:), ESIGPI(:,:), ESIGEI(:,:), ESIGEL(:,:), ESIGPH(:,:),
c  momentum exchange rates, by reaction
     R VSIGCX(:),   VSIGPI(:),   VSIGEL(:),
c  totals
     R SIGCXT,      SIGPIT,      SIGEIT,      SIGELT,      SIGPHT,
     R SIGTOT, SIGBGK,
c  inverse mean free path
     R ZMFPI

      REAL(DP), PUBLIC, SAVE :: ZMFPTHI, TDGTEMX

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(SIGVCX,SIGVPI,SIGVEI,SIGVEL,SIGVPH,XSTOR,XSTORV,
!$OMP& ESIGCX,ESIGPI,ESIGEI,ESIGEL,ESIGPH,VSIGCX,VSIGPI,VSIGEL,SIGCXT,
!$OMP& SIGPIT,SIGEIT,SIGELT,SIGPHT,SIGTOT,SIGBGK,ZMFPI,ZMFPTHI,TDGTEMX)
#endif

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R TABEI1(:,:),   TABRC1(:,:),
     R TABPI3(:,:,:), TABCX3(:,:,:), TABEL3(:,:,:),  ! missing: tabph3
     R FDLMPI(:),     FDLMCX(:),     FDLMEL(:),
     R ADDPI(:,:),    ADDCX(:,:),    ADDEL(:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R FACRRC(:,:), FACRPI(:,:), FACREL(:,:), FACREI(:,:), FACRCX(:,:)

c  secondaries, species distribution, for EI and PI processes
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R PELEI(:),  PATEI(:,:), PMLEI(:,:), PIOEI(:,:), PPLEI(:,:),
     R PELPI(:),  PATPI(:,:), PMLPI(:,:), PIOPI(:,:), PPLPI(:,:),
c  ...and cumulated distributions thereof, for species sampling
     R P2ND(:,:), P2NP(:,:),  P2NEI(:),   P2NPI(:)
c  post-collision energies: to EL (electrons), PI (background) or HV (heavy test particles)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EELEI1(:,:), EHVEI1(:,:),
     R EELRC1(:,:),
     R EPLPI3(:,:,:), EELPI3(:,:,:), EHVPI3(:,:,:),
     R EPLCX3(:,:,:), !  CX processes have no secondary electrons
     R EPLEL3(:,:,:), !  EL processes have no secondary electrons
     R EPLPH3(:,:,:)  !  missing: eelph...,
                      !  PH processes energetics unfinished

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EATEI(:,:,:), EMLEI(:,:,:), EIOEI(:,:,:), EPLEI(:,:,:),
     R EATPI(:,:,:), EMLPI(:,:,:), EIOPI(:,:,:), EPLPI(:,:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I MODCOL(:,:,:),
     I IESTCX(:,:), IESTEL(:,:), IESTPI(:,:), IESTEI(:,:),
     I IESTPH(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I NAEII(:),    NMEII(:),    NIEII(:),  NPHEII(:),
     I NACXI(:),    NMCXI(:),    NICXI(:),  NPHCXI(:),
     I NAELI(:),    NMELI(:),    NIELI(:),  NPHELI(:),
     I NAPII(:),    NMPII(:),    NIPII(:),  NPHPII(:),
     I                                      NPHPHI(:),
     I NPBGKA(:),   NPBGKM(:),   NPBGKI(:), NPBGKPH(:),
     I NPBGKP(:,:)

!  POINTER FOR UNIFIED "A,M,I,PH" SUBROUTINES: UPDATE, ZMFP, (TBD: COLLIDE)
      INTEGER, PUBLIC, POINTER, SAVE ::
cdr  No need for pointer attributes. These are just scalar integers
cdr  Just set these values in SWITCH_PARTINFO.f
     I NXEII, NXCXI, NXELI, NXPII, NXPHI,
     I NXEIIM, NXCXIM, NXELIM, NXPIIM, NXPHIM,
     I NPBGKX
cym p2nei removed from threadprivate list
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(NXEII,NXCXI,NXELI,NXPII,NXPHI,
!$OMP& NXEIIM,NXCXIM,NXELIM,NXPIIM,NXPHIM,NPBGKX)
#endif

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I NAEIIM(:),   NMEIIM(:),   NIEIIM(:),  NPHEIIM(:),
     I NACXIM(:),   NMCXIM(:),   NICXIM(:),  NPHCXIM(:),
     I NAELIM(:),   NMELIM(:),   NIELIM(:),  NPHELIM(:),
     I NAPIIM(:),   NMPIIM(:),   NIPIIM(:),  NPHPIIM(:),
     I                                       NPHPHIM(:),

     I NPRCI(:),    NPRCIM(:),
     I NPPHI(:),    NPPHIM(:)
c  secondaries, species distribution, for RC processes
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NATPRC(:),    NMLPRC(:),   NIOPRC(:),   NPLPRC(:),   NPHPRC(:),
     I NATPRC_2(:),  NMLPRC_2(:), NIOPRC_2(:), NPLPRC_2(:), NPHPRC_2(:),
c  secondaries, species distribution, for CX processes
     I N1STX(:,:), N2NDX(:,:)
c  secondaries, species distribution, for PH processes
cdr  tbd.

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NSEACX(:,:,:), NSEMCX(:,:,:), NSEICX(:,:,:),
     I NSEAEL(:,:,:), NSEMEL(:,:,:), NSEIEL(:,:,:),
     I NSEPRC(:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NREACX(:),NREAPI(:),NREAEL(:),
     I NREAEI(:),NREARC(:),NREAPH(:),
cdr  flags for field particle reaction energetics:
cdr  electrons (net),      ...el...
cdr  incident heavy bulk,  ...hv...
cdr  post coll. heavy bulk,...pl...
     I NELREI(:),JELREI(:),NHVREI(:),
     I                               NELREL(:),   ! --> nplrel, misnomer
     I NELRRC(:),JELRRC(:),
     I NELRPI(:),JELRPI(:),
     I                               NELRCX(:),   ! --> nplrcx, misnomer
     I NELRPH(:),NHVRPI(:),NPLRPI(:),
     I NREACT(:),
     I IPATEI(:,:),IPMLEI(:,:),
     I IPIOEI(:,:),IPPLEI(:,:),
     I IPATPI(:,:),IPMLPI(:,:),
     I IPIOPI(:,:),IPPLPI(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
c  CX type processes
     I LGACX(:,:,:),LGMCX(:,:,:),
     I LGICX(:,:,:),
     I LGPHCX(:,:,:),
c  EI type processes
     I LGAEI(:,:),  LGMEI(:,:),
     I LGIEI(:,:),
     I LGPHEI(:,:),
c  EL type processes
     I LGAEL(:,:,:),LGMEL(:,:,:),
     I LGIEL(:,:,:),
     I LGPHEL(:,:,:),
c  RC type processes
     I LGPRC(:,:),   !  e.g. e + H+ --> H + ph,
                     !  and H is the resulting test particle
c  PI type processes
     I LGAPI(:,:,:),LGMPI(:,:,:),
     I LGIPI(:,:,:),
     I LGPHPI(:,:,:),
c  PH type processes,
c  photon + bulk ph process Bik, Bki, photo-ionisation, etc.
     I LGPHPH(:,:,:),
     I LGPPH(:,:,:)  !  e.g. Aik, or e + H+ --> H + ph,
                     !  and ph (photon) is the resulting test particle

cdr  begin threadprivate here

!  POINTER FOR UNIFIED "A,M,I,PH" SUBROUTINES
cdr values assigned in SWITCH_PARTINFO during particle tracing
cdr April 22: Remove redundant leading dimension in 3d arrays LGX...
cdr           Remove POINTER attribute, use ALLOCATABLE instead
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I LGXCX(:,:), LGXEI(:), LGXEL(:,:), LGXPI(:,:), LGXPH(:,:)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(LGXCX,LGXEI,LGXEL,LGXPI,LGXPH)
#endif

      INTEGER, PUBLIC, SAVE ::
     I NRPII, NREII, NRCXI, NRELI, NRRCI, NRPHI, NROTI, NRBGI

      INTEGER, PUBLIC, SAVE ::
     I NSTOR1, NSTOR,  NSTORV, NTAB, NDAT, NMDTA, MMDTA, NAMF, MAMF,
     I MSTOR0, MSTOR1, MSTOR2

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R DELPOT(:),   FACREA(:,:),
     R FREACA(:,:), FREACM(:,:), FREACI(:,:), FREACP(:,:), FREACPH(:,:),
     R EDPOTA(:,:), EDPOTM(:,:), EDPOTI(:,:), EDPOTP(:,:), EDPOTPH(:,:),
     R EELECA(:,:), EELECM(:,:), EELECI(:,:), EELECP(:,:), EELECPH(:,:),
     R EBULKA(:,:), EBULKM(:,:), EBULKI(:,:), EBULKP(:,:), EBULKPH(:,:),
     R ESCD1A(:,:), ESCD1M(:,:), ESCD1I(:,:), ESCD1P(:,:), ESCD1PH(:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I ISWR(:),     MODCLF(:),   MASSP(:),    MASST(:),
     I IFTFLG(:,:),
     I NRCA(:),     NRCM(:),     NRCI(:),     NRCP(:),     NRCPH(:),
     I IREACA(:,:), IREACM(:,:), IREACI(:,:), IREACP(:,:), IREACPH(:,:),
     I IBULKA(:,:), IBULKM(:,:), IBULKI(:,:), IBULKP(:,:), IBULKPH(:,:),
     I ISCD1A(:,:), ISCD1M(:,:), ISCD1I(:,:), ISCD1P(:,:), ISCD1PH(:,:),
     I ISCD2A(:,:), ISCD2M(:,:), ISCD2I(:,:), ISCD2P(:,:), ISCD2PH(:,:),
     I ISCD3A(:,:), ISCD3M(:,:), ISCD3I(:,:), ISCD3P(:,:), ISCD3PH(:,:),
     I ISCD4A(:,:), ISCD4M(:,:), ISCD4I(:,:), ISCD4P(:,:), ISCD4PH(:,:),
     I ISCDEA(:,:), ISCDEM(:,:), ISCDEI(:,:), ISCDEP(:,:), ISCDEPH(:,:),
     I IESTMA(:,:), IESTMM(:,:), IESTMI(:,:), IESTMPH(:,:),
     I IBGKA (:,:), IBGKM (:,:), IBGKI (:,:), IBGKPH (:,:)

      INTEGER, PUBLIC, SAVE ::
     I NREACI, NHCOL_STORE

      INTEGER, PUBLIC, SAVE :: MAXSPC(0:4)

      INTEGER, ALLOCATABLE, PUBLIC, SAVE :: M_HCOL(:)

cnh  macroscopic length scale
      REAL(DP), PUBLIC, SAVE :: L_MACRO

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COMXS (ICAL)
CDR
C  AUTOMATED ALLOCATION OF STORAGE FOR A&M DATA STRUCTURES AND ARRAYS.
C  CALLED FROM: ALLOCATE_MODULES.F
cdr  ICAL=1: ...?  comments ?
cdr  ICAL=2: ...?  comments ?
cdr  ILONG: predefined type for 8 byte integer variables
cpb         needed as the amount of storage allocated in this module
cpb         exceeds the range of 4 byte integer variables.

      INTEGER, INTENT(IN) :: ICAL
      INTEGER, PARAMETER :: ILONG = SELECTED_INT_KIND(15)
      INTEGER(ILONG) :: MEM

      IF (ICAL == 1) THEN

        IF (ALLOCATED(XSTORV)) RETURN

        NSTORV = 8
C
        NAMF=9*11*11+11+2*(2*11+6*11)+
     P       NREAC*(9*11+18+ 6*NPHOT+ 6*NATM+ 6*NMOL+ 6*NION+ 6*NPLS)
cdr next line: corrected. Was probably incorrect, but unused
        MAMF=NREAC*(      4+10*NPHOT+10*NATM+10*NMOL+10*NION+ 7*NPLS)+
     P       1*NATM+ 1*NMOL+ 1*NION+ 1*NPLS+ 1*NPHOT +
     P       (12+NREAC)*6

        ALLOCATE (XSTORV(NSTORV))

        SIGCXT  => XSTORV(1)
        SIGPIT  => XSTORV(2)
        SIGEIT  => XSTORV(3)
        SIGELT  => XSTORV(4)
        SIGPHT  => XSTORV(5)
        SIGTOT  => XSTORV(6)
        SIGBGK  => XSTORV(7)
        ZMFPI   => XSTORV(8)

cdr  Range for looping over type of processes: DO IXYY=1:NXYYI(NSPZ)
        ALLOCATE (NAEII(NATM))
        ALLOCATE (NMEII(NMOL))
        ALLOCATE (NIEII(NION))
        ALLOCATE (NPHEII(NPHOT))

        ALLOCATE (NACXI(NATM))
        ALLOCATE (NMCXI(NMOL))
        ALLOCATE (NICXI(NION))
        ALLOCATE (NPHCXI(NPHOT))

        ALLOCATE (NAELI(NATM))
        ALLOCATE (NMELI(NMOL))
        ALLOCATE (NIELI(NION))
        ALLOCATE (NPHELI(NPHOT))

        ALLOCATE (NAPII(NATM))
        ALLOCATE (NMPII(NMOL))
        ALLOCATE (NIPII(NION))
        ALLOCATE (NPHPII(NPHOT))

        ALLOCATE (NPHPHI(NPHOT))
        ALLOCATE (NPPHI(NPLS))
c  for background field particles we allow only RC type reactions
        ALLOCATE (NPRCI(NPLS))

cdr  Now the same arrays, all values reduced by one. To avoid unnecessary do loops.
        ALLOCATE (NAEIIM(NATM))
        ALLOCATE (NMEIIM(NMOL))
        ALLOCATE (NIEIIM(NION))
        ALLOCATE (NPHEIIM(NPHOT))

        ALLOCATE (NACXIM(NATM))
        ALLOCATE (NMCXIM(NMOL))
        ALLOCATE (NICXIM(NION))
        ALLOCATE (NPHCXIM(NPHOT))

        ALLOCATE (NAELIM(NATM))
        ALLOCATE (NMELIM(NMOL))
        ALLOCATE (NIELIM(NION))
        ALLOCATE (NPHELIM(NPHOT))

        ALLOCATE (NAPIIM(NATM))
        ALLOCATE (NMPIIM(NMOL))
        ALLOCATE (NIPIIM(NION))
        ALLOCATE (NPHPIIM(NPHOT))
cdr new for photon options
        ALLOCATE (NPHPHIM(NPHOT))
        ALLOCATE (NPPHIM(NPLS))

        ALLOCATE (NPRCIM(NPLS))
c  for BGK reactions
        ALLOCATE (NPBGKA(NATM))
        ALLOCATE (NPBGKM(NMOL))
        ALLOCATE (NPBGKI(NION))
        ALLOCATE (NPBGKPH(NPHOT))
        ALLOCATE (NPBGKP(NPLS,2))

        ALLOCATE (NSEACX(NATM,NPLS,5))
        ALLOCATE (NSEMCX(NMOL,NPLS,5))
        ALLOCATE (NSEICX(NION,NPLS,5))

        ALLOCATE (NSEAEL(NATM,NPLS,5))
        ALLOCATE (NSEMEL(NMOL,NPLS,5))
        ALLOCATE (NSEIEL(NION,NPLS,5))

        ALLOCATE (DELPOT(NREAC))

        ALLOCATE (FACREA(-11:NREAC,2))

        ALLOCATE (FREACA(NATM,NREAC))
        ALLOCATE (FREACM(NMOL,NREAC))
        ALLOCATE (FREACI(NION,NREAC))
        ALLOCATE (FREACP(NPLS,NREAC))
        ALLOCATE (FREACPH(NPHOT,NREAC))

cdr  former fluid limit, now contained in ngen..
cdr  fldlm=10000/[-(1+ngen)](generation limit: negative values)


cdr new: potential difference in a particular reaction.
cdr      allows to derive radiation loss from
cdr                       electron energy loss PELEC (=eelec)
cdr                                            KER (=escd1)
cdr                                            POT (=edpot)
cdr   Warning: strictly PRAD and PELEC refer to a species,
cdr            not to a reaction. In case of multiple excit
cdr            channels (e.g. H2, or MS resolved atoms)
cdr            great care is needed not to double count
cdr            or miss out parts of the electron cooling or radiation.
cdr
cdr
cdr      PRAD= PELEC-KER-POT
cdr
        ALLOCATE (EDPOTA(NATM,NREAC))
        ALLOCATE (EDPOTM(NMOL,NREAC))
        ALLOCATE (EDPOTI(NION,NREAC))
        ALLOCATE (EDPOTP(NPLS,NREAC))
        ALLOCATE (EDPOTPH(NPHOT,NREAC))

        ALLOCATE (EELECA(NATM,NREAC))
        ALLOCATE (EELECM(NMOL,NREAC))
        ALLOCATE (EELECI(NION,NREAC))
        ALLOCATE (EELECP(NPLS,NREAC))
        ALLOCATE (EELECPH(NPHOT,NREAC))

        ALLOCATE (EBULKA(NATM,NREAC))
        ALLOCATE (EBULKM(NMOL,NREAC))
        ALLOCATE (EBULKI(NION,NREAC))
        ALLOCATE (EBULKP(NPLS,NREAC))
        ALLOCATE (EBULKPH(NPHOT,NREAC))

        ALLOCATE (ESCD1A(NATM,NREAC))
        ALLOCATE (ESCD1M(NMOL,NREAC))
        ALLOCATE (ESCD1I(NION,NREAC))
        ALLOCATE (ESCD1P(NPLS,NREAC))
        ALLOCATE (ESCD1PH(NPHOT,NREAC))

        ALLOCATE (ISWR(NREAC))
        ALLOCATE (MODCLF(NREAC))
        ALLOCATE (MASSP(NREAC))
        ALLOCATE (MASST(NREAC))
        ALLOCATE (IFTFLG(-11:NREAC,0:5))

        ALLOCATE (NRCA(NATM))
        ALLOCATE (NRCM(NMOL))
        ALLOCATE (NRCI(NION))
        ALLOCATE (NRCP(NPLS))
        ALLOCATE (NRCPH(NPHOT))

        ALLOCATE (IREACA(NATM,NREAC))
        ALLOCATE (IREACM(NMOL,NREAC))
        ALLOCATE (IREACI(NION,NREAC))
        ALLOCATE (IREACP(NPLS,NREAC))
        ALLOCATE (IREACPH(NPHOT,NREAC))

        ALLOCATE (IBULKA(NATM,NREAC))
        ALLOCATE (IBULKM(NMOL,NREAC))
        ALLOCATE (IBULKI(NION,NREAC))
        ALLOCATE (IBULKP(NPLS,NREAC))
        ALLOCATE (IBULKPH(NPHOT,NREAC))

        ALLOCATE (ISCD1A(NATM,NREAC))
        ALLOCATE (ISCD1M(NMOL,NREAC))
        ALLOCATE (ISCD1I(NION,NREAC))
        ALLOCATE (ISCD1P(NPLS,NREAC))
        ALLOCATE (ISCD1PH(NPHOT,NREAC))

        ALLOCATE (ISCD2A(NATM,NREAC))
        ALLOCATE (ISCD2M(NMOL,NREAC))
        ALLOCATE (ISCD2I(NION,NREAC))
        ALLOCATE (ISCD2P(NPLS,NREAC))
        ALLOCATE (ISCD2PH(NPHOT,NREAC))

        ALLOCATE (ISCD3A(NATM,NREAC))
        ALLOCATE (ISCD3M(NMOL,NREAC))
        ALLOCATE (ISCD3I(NION,NREAC))
        ALLOCATE (ISCD3P(NPLS,NREAC))
        ALLOCATE (ISCD3PH(NPHOT,NREAC))

        ALLOCATE (ISCD4A(NATM,NREAC))
        ALLOCATE (ISCD4M(NMOL,NREAC))
        ALLOCATE (ISCD4I(NION,NREAC))
        ALLOCATE (ISCD4P(NPLS,NREAC))
        ALLOCATE (ISCD4PH(NPHOT,NREAC))

        ALLOCATE (ISCDEA(NATM,NREAC))
        ALLOCATE (ISCDEM(NMOL,NREAC))
        ALLOCATE (ISCDEI(NION,NREAC))
        ALLOCATE (ISCDEP(NPLS,NREAC))
        ALLOCATE (ISCDEPH(NPHOT,NREAC))

        ALLOCATE (IESTMA(NATM,NREAC))
        ALLOCATE (IESTMM(NMOL,NREAC))
        ALLOCATE (IESTMI(NION,NREAC))
        ALLOCATE (IESTMPH(NPHOT,NREAC))

        ALLOCATE (IBGKA (NATM,NREAC))
        ALLOCATE (IBGKM (NMOL,NREAC))
        ALLOCATE (IBGKI (NION,NREAC))
        ALLOCATE (IBGKPH(NPHOT,NREAC))

cdr  -11 ... -1   : internal default atomic-molecular data
cdr    1 ... NREAC: atomic/molecular data read from external data files, input block 4
        ALLOCATE (REACDAT(-11:NREAC))
        ALLOCATE (REACLINES(NREAC_LINES))

        ALLOCATE (M_HCOL(NREAC))

        MEM = (NSTORV+NAMF)*8_ILONG + (MAMF+
     .                      9_ILONG*(NATM+NMOL+NION)+4_ILONG*NPLS+
     .                      10_ILONG*NPLS*(NATM+NMOL+NION))*4_ILONG +
     .                      NREAC*4_ILONG

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMXS(1) ', MEM

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(XSTOR)) RETURN
C  DIMENSION OF FULL REACTION-SPECIFIC ARRAYS: CFLAG, MODCOL,....
        MSTOR0 = MAX(NRPI, NREI, NRCX, NREL, NRPH, NREC)
C  FIRST DIMENSION OF XSTOR ARRAY
        MSTOR1 = MAX(NRPI, NREI, NRCX, NREL, NRPH)
C  SECOND DIMENSION OF XSTOR ARRAY
        MSTOR2 = 24

cdr  ...to be removed
cdr nstor1 and nstor: unused, 16.2.2020 !
        NSTOR1 = NREL+NRCX+NRPI+NREI
        NSTOR  = NSTOR1+
     .           2*(NREL+NRCX+NRPI)+5*NREI+
     .           NREL+NRCX+NRPI
cdr  ...
C
c  collision rates for MFP, free flight sampling
cdr     also: for particle loss/gain rates for tracklength estimators
c  tab..1/3 arrays  (photon processes missing ?)
        NTAB=NSTORDR*(NREI+NREC)+
     P       NSTORDR*NSTORDT*(NRCX+NREL+NRPI)+
     P       (NPLS+1)*(NRPI+NRCX+NREL)+
     P       2*(NREC+NRPI+NREL+NREI+NRCX)
C
cdr rates for energy loss/gain rates for tracklength estimators
c  exx..1/3 arrays  xx=(cx,el,pl,hv),
        NDAT=NSTORDR*(2*NREI+NREC)+
     P       NSTORDR*NSTORDT*(NRCX+NREL+3*NRPI+NRPH)+
c
     P      (NREI+NRPI)*
     P      (NATMP+NMOLP+NIONP+NPLSP+1)+
     P      (NRPI+NREI)*(NSPZP+1)+
     P       NRPI*2*(NATMP+NMOLP+NIONP+1)+
     P       NREI*2*(NATMP+NMOLP+NIONP+1)
C
        NMDTA=NTAB+NDAT
C
        MMDTA=7*5*MSTOR0+3*(NRCX+NREL+NRPI+NREI)+6+
     P        5*NREC+
     P        6*NRCX+
     P        10*NREC+
     P        2*NRCX+4*NRPI+2*NREL+4*NREI+3*NREC+NREAC+2*NRPH+
     P        (NREI+NRPI)*
     P        (NATMP+NMOLP+NIONP+NPLSP)+
C  LG... ARRAYS
     P          (NATMP+NMOLP+NIONP      )*(NREI+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NRCX+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NREL+1)+
     P          (                  NPLSP)*(NREC+1)+
     P        2*(NATMP+NMOLP+NIONP      )*(NRPI+1)

        ALLOCATE (XSTOR(MSTOR1,MSTOR2))

        SIGVCX => XSTOR(:,1)
        SIGVPI => XSTOR(:,2)
        SIGVEI => XSTOR(:,3)
        SIGVEL => XSTOR(:,4)
        SIGVPH => XSTOR(:,22)

        ESIGCX => XSTOR(:,5:6)
        ESIGPI => XSTOR(:,7:11)
        ESIGEI => XSTOR(:,12:16)
        ESIGEL => XSTOR(:,17:18)
        ESIGPH => XSTOR(:,23:24)

        VSIGCX => XSTOR(:,19)
        VSIGPI => XSTOR(:,20)
        VSIGEL => XSTOR(:,21)
cdr     vsigei  : still missing
cdr     vsigph  : still missing

        ALLOCATE (TABEI1(NREI,NSTORDR))
        ALLOCATE (TABRC1(NREC,NSTORDR))
        ALLOCATE (TABPI3(NRPI,NSTORDR,NSTORDT))
        ALLOCATE (TABCX3(NRCX,NSTORDR,NSTORDT))
        ALLOCATE (TABEL3(NREL,NSTORDR,NSTORDT))

        ALLOCATE (FDLMPI(NRPI))
        ALLOCATE (FDLMCX(NRCX))
        ALLOCATE (FDLMEL(NREL))

c  factors for scaling reaction rates to other target (ipls) masses
c  (only for heavy particle impact reaction)
        ALLOCATE (ADDPI(NRPI,NPLS))
        ALLOCATE (ADDCX(NRCX,NPLS))
        ALLOCATE (ADDEL(NREL,NPLS))

c  factors for scaling reaction processes
        ALLOCATE (FACRRC(NREC,2))
        ALLOCATE (FACRPI(NRPI,2))
        ALLOCATE (FACREL(NREL,2))
        ALLOCATE (FACREI(NREI,2))
        ALLOCATE (FACRCX(NRCX,2))

c  secondaries, EI processes
        ALLOCATE (PELEI(NREI))
        ALLOCATE (PATEI(NREI,0:NATM))
        ALLOCATE (PMLEI(NREI,0:NMOL))
        ALLOCATE (PIOEI(NREI,0:NION))
        ALLOCATE (PPLEI(NREI,0:NPLS))
        ALLOCATE (P2ND(NREI,0:NSPZ))
        ALLOCATE (P2NEI(NREI))
c  secondaries, PI processes
        ALLOCATE (PELPI(NRPI))
        ALLOCATE (PATPI(NRPI,0:NATM))
        ALLOCATE (PMLPI(NRPI,0:NMOL))
        ALLOCATE (PIOPI(NRPI,0:NION))
        ALLOCATE (PPLPI(NRPI,0:NPLS))
        ALLOCATE (P2NP(NRPI,0:NSPZ))
        ALLOCATE (P2NPI(NRPI))

        ALLOCATE (EELEI1(NREI,NSTORDR))
        ALLOCATE (EHVEI1(NREI,NSTORDR))

        ALLOCATE (EELRC1(NREC,NSTORDR))

        ALLOCATE (EELPI3(NRPI,NSTORDR,NSTORDT))
        ALLOCATE (EHVPI3(NRPI,NSTORDR,NSTORDT))
        ALLOCATE (EPLPI3(NRPI,NSTORDR,NSTORDT))

        ALLOCATE (EPLCX3(NRCX,NSTORDR,NSTORDT))
        ALLOCATE (EPLEL3(NREL,NSTORDR,NSTORDT))

        ALLOCATE (EPLPH3(NRPH,NSTORDR,NSTORDT))

        ALLOCATE (EATPI(NRPI,0:NATM,2))
        ALLOCATE (EMLPI(NRPI,0:NMOL,2))
        ALLOCATE (EIOPI(NRPI,0:NION,2))
        ALLOCATE (EPLPI(NRPI,0:NPLS,2))

        ALLOCATE (EATEI(NREI,0:NATM,2))
        ALLOCATE (EMLEI(NREI,0:NMOL,2))
        ALLOCATE (EIOEI(NREI,0:NION,2))
        ALLOCATE (EPLEI(NREI,0:NPLS,2))

        ALLOCATE (MODCOL(7,0:4,MSTOR0))

c   flags for collision or tracklength estimators,
c   for particle (1), momentum (2) and energy (3) source rates, resp.
        ALLOCATE (IESTCX(NRCX,3))
        ALLOCATE (IESTEL(NREL,3))
        ALLOCATE (IESTPI(NRPI,3))
        ALLOCATE (IESTEI(NREI,3))
        ALLOCATE (IESTPH(NRPH,3))

        ALLOCATE (NATPRC(NREC))
        ALLOCATE (NMLPRC(NREC))
        ALLOCATE (NIOPRC(NREC))
        ALLOCATE (NPLPRC(NREC))
        ALLOCATE (NPHPRC(NREC))

        ALLOCATE (NATPRC_2(NREC))
        ALLOCATE (NMLPRC_2(NREC))
        ALLOCATE (NIOPRC_2(NREC))
        ALLOCATE (NPLPRC_2(NREC))
        ALLOCATE (NPHPRC_2(NREC))

        ALLOCATE (N1STX(NRCX,3))
        ALLOCATE (N2NDX(NRCX,3))

        ALLOCATE (NSEPRC(NREC,5))

        ALLOCATE (NREACX(NRCX))
        ALLOCATE (NREAPI(NRPI))
        ALLOCATE (NREAEL(NREL))
        ALLOCATE (NREAEI(NREI))
        ALLOCATE (NREARC(NREC))
        ALLOCATE (NREAPH(NRPH))

        ALLOCATE (NELREI(NREI))
        ALLOCATE (JELREI(NREI))
        ALLOCATE (NHVREI(NREI))
        ALLOCATE (NELREL(NREL))
        ALLOCATE (NELRRC(NREC))
        ALLOCATE (JELRRC(NREC))
        ALLOCATE (NELRPI(NRPI))
        ALLOCATE (JELRPI(NRPI))
        ALLOCATE (NHVRPI(NRPI))
        ALLOCATE (NPLRPI(NRPI))
        ALLOCATE (NELRCX(NRCX))
        ALLOCATE (NELRPH(NRPH))

        ALLOCATE (NREACT(NREAC))
c  again: some arrays for species distribution of secondaries
c         derived from P..EI and P..PI, above.
c         for speeding up scoring in update, collide
cdr unclear meaning. Perhaps redundant. But might be useful?
        ALLOCATE (IPATEI(NREI,0:NATM))
        ALLOCATE (IPMLEI(NREI,0:NMOL))
        ALLOCATE (IPIOEI(NREI,0:NION))
        ALLOCATE (IPPLEI(NREI,0:NPLS))
        ALLOCATE (IPATPI(NRPI,0:NATM))
        ALLOCATE (IPMLPI(NRPI,0:NMOL))
        ALLOCATE (IPIOPI(NRPI,0:NION))
        ALLOCATE (IPPLPI(NRPI,0:NPLS))
c
        ALLOCATE (LGACX(0:NATM,0:NRCX,0:1))
        ALLOCATE (LGMCX(0:NMOL,0:NRCX,0:1))
        ALLOCATE (LGICX(0:NION,0:NRCX,0:1))
        ALLOCATE (LGPHCX(0:NPHOT,0:NRCX,0:1))

        ALLOCATE (LGAEI(0:NATM,0:NREI))
        ALLOCATE (LGMEI(0:NMOL,0:NREI))
        ALLOCATE (LGIEI(0:NION,0:NREI))
        ALLOCATE (LGPHEI(0:NPHOT,0:NREI))

        ALLOCATE (LGAEL(0:NATM,0:NREL,0:1))
        ALLOCATE (LGMEL(0:NMOL,0:NREL,0:1))
        ALLOCATE (LGIEL(0:NION,0:NREL,0:1))
        ALLOCATE (LGPHEL(0:NPHOT,0:NREL,0:1))

        ALLOCATE (LGPRC(0:NPLS,0:NREC))

        ALLOCATE (LGAPI(0:NATM,0:NRPI,0:1))
        ALLOCATE (LGMPI(0:NMOL,0:NRPI,0:1))
        ALLOCATE (LGIPI(0:NION,0:NRPI,0:1))
        ALLOCATE (LGPHPI(0:NPHOT,0:NRPI,0:1))

        ALLOCATE (LGPPH(0:NPLS,0:NRPH,0:1))   ! new in
        ALLOCATE (LGPHPH(0:NPHOT,0:NRPH,0:5))  ! new in

cdr  for unified routines UPDATE, FPATH, COLLIDE
cdr  (single code of atoms, molecules and test ions. tbd: photons)
!pb  allocation moved to EIRENE_SWITCH_PARTINFO
!pb  This is necessary for OPENMP as the allocation needs to done
!pb  at a place which all threads visit.
!pb     ALLOCATE (LGXCX(0:NRCX,0:1))
!pb     ALLOCATE (LGXEI(0:NREI))
!pb     ALLOCATE (LGXEL(0:NREL,0:1))
!pb     ALLOCATE (LGXPI(0:NRPI,0:1))
!pb     ALLOCATE (LGXPH(0:NRPH,0:1))

        MEM = (MSTOR1*MSTOR2+NMDTA)*8_ILONG +
     .                       MMDTA*4_ILONG

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMXS(2) ', MEM

      END IF

      CALL EIRENE_INIT_CMDTA (ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMXS


      SUBROUTINE EIRENE_DEALLOC_COMXS

      IF (.NOT.ALLOCATED(XSTOR)) RETURN

      DEALLOCATE (XSTOR)
      DEALLOCATE (XSTORV)

      DEALLOCATE (TABEI1)
      DEALLOCATE (TABRC1)
      DEALLOCATE (TABPI3)
      DEALLOCATE (TABCX3)
      DEALLOCATE (TABEL3)
      DEALLOCATE (FDLMPI)
      DEALLOCATE (FDLMCX)
      DEALLOCATE (FDLMEL)
      DEALLOCATE (ADDPI)
      DEALLOCATE (ADDCX)
      DEALLOCATE (ADDEL)

      DEALLOCATE (FACRRC)
      DEALLOCATE (FACRPI)
      DEALLOCATE (FACREL)
      DEALLOCATE (FACREI)
      DEALLOCATE (FACRCX)

      DEALLOCATE (PELEI)
      DEALLOCATE (PATEI)
      DEALLOCATE (PMLEI)
      DEALLOCATE (PIOEI)
      DEALLOCATE (PPLEI)
      DEALLOCATE (PELPI)
      DEALLOCATE (PATPI)
      DEALLOCATE (PMLPI)
      DEALLOCATE (PIOPI)
      DEALLOCATE (PPLPI)
      DEALLOCATE (P2ND)
      DEALLOCATE (P2NP)
      DEALLOCATE (P2NEI)
      DEALLOCATE (P2NPI)

      DEALLOCATE (EELEI1)
      DEALLOCATE (EHVEI1)

      DEALLOCATE (EELRC1)
      DEALLOCATE (EELPI3)
      DEALLOCATE (EHVPI3)
      DEALLOCATE (EPLPI3)
      DEALLOCATE (EPLCX3)
      DEALLOCATE (EPLEL3)
      DEALLOCATE (EPLPH3)

      DEALLOCATE (EATPI)
      DEALLOCATE (EMLPI)
      DEALLOCATE (EIOPI)
      DEALLOCATE (EPLPI)

      DEALLOCATE (EATEI)
      DEALLOCATE (EMLEI)
      DEALLOCATE (EIOEI)
      DEALLOCATE (EPLEI)

      DEALLOCATE (MODCOL)

      DEALLOCATE (IESTCX)
      DEALLOCATE (IESTEL)
      DEALLOCATE (IESTPI)
      DEALLOCATE (IESTEI)
      DEALLOCATE (IESTPH)

      DEALLOCATE (NAEII)
      DEALLOCATE (NMEII)
      DEALLOCATE (NIEII)
      DEALLOCATE (NPHEII)

      DEALLOCATE (NACXI)
      DEALLOCATE (NMCXI)
      DEALLOCATE (NICXI)
      DEALLOCATE (NPHCXI)

      DEALLOCATE (NAELI)
      DEALLOCATE (NMELI)
      DEALLOCATE (NIELI)
      DEALLOCATE (NPHELI)

      DEALLOCATE (NAPII)
      DEALLOCATE (NMPII)
      DEALLOCATE (NIPII)
      DEALLOCATE (NPHPII)

      DEALLOCATE (NPHPHI)
      DEALLOCATE (NPPHI)

      DEALLOCATE (NPRCI)

      DEALLOCATE (NAEIIM)
      DEALLOCATE (NMEIIM)
      DEALLOCATE (NIEIIM)
      DEALLOCATE (NPHEIIM)

      DEALLOCATE (NACXIM)
      DEALLOCATE (NMCXIM)
      DEALLOCATE (NICXIM)
      DEALLOCATE (NPHCXIM)

      DEALLOCATE (NAELIM)
      DEALLOCATE (NMELIM)
      DEALLOCATE (NIELIM)
      DEALLOCATE (NPHELIM)

      DEALLOCATE (NAPIIM)
      DEALLOCATE (NMPIIM)
      DEALLOCATE (NIPIIM)
      DEALLOCATE (NPHPIIM)

      DEALLOCATE (NPHPHIM)
      DEALLOCATE (NPPHIM)

      DEALLOCATE (NPRCIM)

      DEALLOCATE (NPBGKA)
      DEALLOCATE (NPBGKM)
      DEALLOCATE (NPBGKI)
      DEALLOCATE (NPBGKPH)
      DEALLOCATE (NPBGKP)

      DEALLOCATE (NATPRC_2)
      DEALLOCATE (NMLPRC_2)
      DEALLOCATE (NIOPRC_2)
      DEALLOCATE (NPLPRC_2)
      DEALLOCATE (NPHPRC_2)

      DEALLOCATE (NATPRC)
      DEALLOCATE (NMLPRC)
      DEALLOCATE (NIOPRC)
      DEALLOCATE (NPLPRC)
      DEALLOCATE (NPHPRC)

      DEALLOCATE (N1STX)
      DEALLOCATE (N2NDX)

      DEALLOCATE (NSEACX)
      DEALLOCATE (NSEMCX)
      DEALLOCATE (NSEICX)
      DEALLOCATE (NSEAEL)
      DEALLOCATE (NSEMEL)
      DEALLOCATE (NSEIEL)
      DEALLOCATE (NSEPRC)

      DEALLOCATE (NREACX)
      DEALLOCATE (NREAPI)
      DEALLOCATE (NREAEL)
      DEALLOCATE (NREAEI)
      DEALLOCATE (NREARC)
      DEALLOCATE (NREAPH)

      DEALLOCATE (NELREI)
      DEALLOCATE (JELREI)
      DEALLOCATE (NHVREI)
      DEALLOCATE (NELREL)
      DEALLOCATE (NELRRC)
      DEALLOCATE (JELRRC)
      DEALLOCATE (NELRPI)
      DEALLOCATE (JELRPI)
      DEALLOCATE (NHVRPI)
      DEALLOCATE (NPLRPI)
      DEALLOCATE (NELRCX)
      DEALLOCATE (NELRPH)

      DEALLOCATE (NREACT)

      DEALLOCATE (IPATEI)
      DEALLOCATE (IPMLEI)
      DEALLOCATE (IPIOEI)
      DEALLOCATE (IPPLEI)
      DEALLOCATE (IPATPI)
      DEALLOCATE (IPMLPI)
      DEALLOCATE (IPIOPI)
      DEALLOCATE (IPPLPI)

      DEALLOCATE (LGACX)
      DEALLOCATE (LGMCX)
      DEALLOCATE (LGICX)
      DEALLOCATE (LGPHCX)

      DEALLOCATE (LGAEI)
      DEALLOCATE (LGMEI)
      DEALLOCATE (LGIEI)
      DEALLOCATE (LGPHEI)

      DEALLOCATE (LGAEL)
      DEALLOCATE (LGMEL)
      DEALLOCATE (LGIEL)
      DEALLOCATE (LGPHEL)

      DEALLOCATE (LGPRC)

      DEALLOCATE (LGAPI)
      DEALLOCATE (LGMPI)
      DEALLOCATE (LGIPI)
      DEALLOCATE (LGPHPI)

      DEALLOCATE (LGPPH)   ! new in
      DEALLOCATE (LGPHPH)  ! new in

      IF (ALLOCATED(LGXCX)) THEN
        DEALLOCATE (LGXCX)
        DEALLOCATE (LGXEI)
        DEALLOCATE (LGXEL)
        DEALLOCATE (LGXPI)
        DEALLOCATE (LGXPH)
      END IF

      DEALLOCATE (DELPOT)

      DEALLOCATE (FACREA)

      DEALLOCATE (FREACA)
      DEALLOCATE (FREACM)
      DEALLOCATE (FREACI)
      DEALLOCATE (FREACP)
      DEALLOCATE (FREACPH)

      DEALLOCATE (EDPOTA)
      DEALLOCATE (EDPOTM)
      DEALLOCATE (EDPOTI)
      DEALLOCATE (EDPOTP)
      DEALLOCATE (EDPOTPH)

      DEALLOCATE (EELECA)
      DEALLOCATE (EELECM)
      DEALLOCATE (EELECI)
      DEALLOCATE (EELECP)
      DEALLOCATE (EELECPH)

      DEALLOCATE (EBULKA)
      DEALLOCATE (EBULKM)
      DEALLOCATE (EBULKI)
      DEALLOCATE (EBULKP)
      DEALLOCATE (EBULKPH)

      DEALLOCATE (ESCD1A)
      DEALLOCATE (ESCD1M)
      DEALLOCATE (ESCD1I)
      DEALLOCATE (ESCD1P)
      DEALLOCATE (ESCD1PH)

      DEALLOCATE (ISWR)
      DEALLOCATE (MODCLF)
      DEALLOCATE (MASSP)
      DEALLOCATE (MASST)
      DEALLOCATE (IFTFLG)
      DEALLOCATE (NRCP)
      DEALLOCATE (NRCA)
      DEALLOCATE (NRCM)
      DEALLOCATE (NRCI)
      DEALLOCATE (NRCPH)

      DEALLOCATE (IREACA)
      DEALLOCATE (IREACM)
      DEALLOCATE (IREACI)
      DEALLOCATE (IREACP)
      DEALLOCATE (IREACPH)

      DEALLOCATE (IBULKA)
      DEALLOCATE (IBULKM)
      DEALLOCATE (IBULKI)
      DEALLOCATE (IBULKP)
      DEALLOCATE (IBULKPH)

      DEALLOCATE (ISCD1A)
      DEALLOCATE (ISCD1M)
      DEALLOCATE (ISCD1I)
      DEALLOCATE (ISCD1P)
      DEALLOCATE (ISCD1PH)

      DEALLOCATE (ISCD2A)
      DEALLOCATE (ISCD2M)
      DEALLOCATE (ISCD2I)
      DEALLOCATE (ISCD2P)
      DEALLOCATE (ISCD2PH)

      DEALLOCATE (ISCD3A)
      DEALLOCATE (ISCD3M)
      DEALLOCATE (ISCD3I)
      DEALLOCATE (ISCD3P)
      DEALLOCATE (ISCD3PH)

      DEALLOCATE (ISCD4A)
      DEALLOCATE (ISCD4M)
      DEALLOCATE (ISCD4I)
      DEALLOCATE (ISCD4P)
      DEALLOCATE (ISCD4PH)

      DEALLOCATE (ISCDEA)
      DEALLOCATE (ISCDEM)
      DEALLOCATE (ISCDEI)
      DEALLOCATE (ISCDEP)
      DEALLOCATE (ISCDEPH)

      DEALLOCATE (IESTMA)
      DEALLOCATE (IESTMM)
      DEALLOCATE (IESTMI)
      DEALLOCATE (IESTMPH)

      DEALLOCATE (IBGKA )
      DEALLOCATE (IBGKM )
      DEALLOCATE (IBGKI )
      DEALLOCATE (IBGKPH)

      CALL EIRENE_FREE_REACDAT
      DEALLOCATE (REACLINES)

      DEALLOCATE (M_HCOL)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMXS


      SUBROUTINE EIRENE_INIT_CMDTA (ICAL)
cdr  initialize (nullify) A&M data
cdr  ical=1: ??
cdr  ical=2: ??

      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: IREAC

      IF (ICAL == 1) THEN
        NAEII   = 0
        NMEII   = 0
        NIEII   = 0
        NPHEII  = 0

        NACXI   = 0
        NMCXI   = 0
        NICXI   = 0
        NPHCXI  = 0

        NAELI   = 0
        NMELI   = 0
        NIELI   = 0
        NPHELI  = 0

        NAPII   = 0
        NMPII   = 0
        NIPII   = 0
        NPHPII  = 0

        NPHPHI  = 0
        NPPHI   = 0

        NPRCI   = 0

        NAEIIM  = 0
        NMEIIM  = 0
        NIEIIM  = 0
        NPHEIIM = 0

        NACXIM  = 0
        NMCXIM  = 0
        NICXIM  = 0
        NPHCXIM = 0

        NAELIM  = 0
        NMELIM  = 0
        NIELIM  = 0
        NPHELIM = 0

        NAPIIM  = 0
        NMPIIM  = 0
        NIPIIM  = 0
        NPHPIIM = 0

        NPHPHIM = 0
        NPPHIM  = 0

        NPRCIM  = 0

        NPBGKA  = 0
        NPBGKM  = 0
        NPBGKI  = 0
        NPBGKP  = 0

        NSEACX  = 0
        NSEMCX  = 0
        NSEICX  = 0
        NSEAEL  = 0
        NSEMEL  = 0
        NSEIEL  = 0

        DELPOT  = 0._DP
        FACREA(:,1) = 1._DP
        FACREA(:,2) = 0._DP
        FREACA  = 0._DP
        FREACM  = 0._DP
        FREACI  = 0._DP
        FREACP  = 0._DP
        FREACPH = 0._DP

        EDPOTA  = 0._DP
        EDPOTM  = 0._DP
        EDPOTI  = 0._DP
        EDPOTP  = 0._DP
        EDPOTPH = 0._DP

        EELECA  = 0._DP
        EELECM  = 0._DP
        EELECI  = 0._DP
        EELECP  = 0._DP
        EELECPH = 0._DP

        EBULKA  = 0._DP
        EBULKM  = 0._DP
        EBULKI  = 0._DP
        EBULKP  = 0._DP
        EBULKPH = 0._DP

        ESCD1A  = 0._DP
        ESCD1M  = 0._DP
        ESCD1I  = 0._DP
        ESCD1P  = 0._DP
        ESCD1PH = 0._DP

        ZMFPTHI = 0._DP
        TDGTEMX = 1.E-30_DP

        ISWR    = 0
        MODCLF  = 0
        MASSP   = 0
        MASST   = 0
        IFTFLG  = 0
        NRCP    = 0
        NRCA    = 0
        NRCM    = 0
        NRCI    = 0
        NRCPH   = 0

cdr  input flags for collision kinetics, secondaries, etc.
        IREACA  = 0
        IREACM  = 0
        IREACI  = 0
        IREACP  = 0
        IREACPH = 0
        IBULKA  = 0
        IBULKM  = 0
        IBULKI  = 0
        IBULKP  = 0
        IBULKPH = 0
        ISCD1A  = 0
        ISCD1M  = 0
        ISCD1I  = 0
        ISCD1P  = 0
        ISCD1PH = 0
        ISCD2A  = 0
        ISCD2M  = 0
        ISCD2I  = 0
        ISCD2P  = 0
        ISCD2PH = 0
        ISCD3A  = 0
        ISCD3M  = 0
        ISCD3I  = 0
        ISCD3P  = 0
        ISCD3PH = 0
        ISCD4A  = 0
        ISCD4M  = 0
        ISCD4I  = 0
        ISCD4P  = 0
        ISCD4PH = 0
        ISCDEA  = 0
        ISCDEM  = 0
        ISCDEI  = 0
        ISCDEP  = 0
        ISCDEPH = 0
        IESTMA  = 0
        IESTMM  = 0
        IESTMI  = 0
        IESTMPH = 0
        IBGKA   = 0
        IBGKM   = 0
        IBGKI   = 0
        IBGKPH  = 0

        idreac  = 0

        XSTORV  = 0._DP

cdr  ireac=-11, to ireac=-1   : minimal (hard-coded) set of default reactions
cdr  ireac=1    to ireac=nreac: reaction data sets read from external files
        DO IREAC= -11, NREAC
          REACDAT(IREAC)%LPOT   = .FALSE.
          REACDAT(IREAC)%LCRS   = .FALSE.
          REACDAT(IREAC)%LRTC   = .FALSE.
          REACDAT(IREAC)%LRTCMW = .FALSE.
          REACDAT(IREAC)%LRTCEW = .FALSE.
          REACDAT(IREAC)%LOTH   = .FALSE.
          REACDAT(IREAC)%LPHR   = .FALSE.
          REACDAT(IREAC)%NOSEC  = 0

c  some universal data for reaction no. ireac
c  data needed for rejection sampling in velocx, veloel, velopi
   ! max value of vel times sigma(vel)
          REACDAT(IREAC)%RTMAX  = 0._DP
c  kinetic collision energy at which this maximum is attained.
          REACDAT(IREAC)%ERTMAX = -HUGE(1._DP)
c  reaction threshold (if any)
          REACDAT(IREAC)%ETH    = 0._DP
c  total kinetic energy release (if any)
          REACDAT(IREAC)%KER    = 0._DP
c  potential energy difference reactants -- products
          REACDAT(IREAC)%DELP   = 0._DP
c  Arrhenius prefactor exp(-EARRH/T) separated from rest of fit
          REACDAT(IREAC)%EARRH0 = 0._DP
c  ditto, for upper limit value of 2nd parameter in 2-parameter fits or tables
          REACDAT(IREAC)%EARRH1 = 0._DP

          NULLIFY(REACDAT(IREAC)%POT)
          NULLIFY(REACDAT(IREAC)%CRS)
          NULLIFY(REACDAT(IREAC)%RTC)
          NULLIFY(REACDAT(IREAC)%RTCMW)
          NULLIFY(REACDAT(IREAC)%RTCEW)
          NULLIFY(REACDAT(IREAC)%OTH)
          NULLIFY(REACDAT(IREAC)%PHR)
        END DO

        NHCOL_STORE = 0
        M_HCOL = 0

cnh
        L_MACRO = 0.0_DP

      ELSE IF (ICAL == 2) THEN

        XSTOR  = 0._DP

        TABEI1  = 0._DP
        TABRC1  = 0._DP
        TABPI3  = 0._DP
        TABCX3  = 0._DP
        TABEL3  = 0._DP
        FDLMPI  = 0._DP
        FDLMCX  = 0._DP
        FDLMEL  = 0._DP
        ADDPI   = 0._DP
        ADDCX   = 0._DP
        ADDEL   = 0._DP

        FACRRC(:,1) = 1._DP
        FACRRC(:,2) = 0._DP
        FACRPI(:,1) = 1._DP
        FACRPI(:,2) = 0._DP
        FACREL(:,1) = 1._DP
        FACREL(:,2) = 0._DP
        FACREI(:,1) = 1._DP
        FACREI(:,2) = 0._DP
        FACRCX(:,1) = 1._DP
        FACRCX(:,2) = 0._DP

        PELEI   = 0._DP
        PATEI   = 0._DP
        PMLEI   = 0._DP
        PIOEI   = 0._DP
        PPLEI   = 0._DP
        PELPI   = 0._DP
        PATPI   = 0._DP
        PMLPI   = 0._DP
        PIOPI   = 0._DP
        PPLPI   = 0._DP
        P2ND    = 0._DP
        P2NP    = 0._DP
        P2NEI   = 0._DP
        P2NPI   = 0._DP

        EELEI1  = 0._DP
        EHVEI1  = 0._DP
        EELRC1  = 0._DP
        EELPI3  = 0._DP
        EHVPI3  = 0._DP
        EPLPI3  = 0._DP
        EPLCX3  = 0._DP
        EPLEL3  = 0._DP
        EPLPH3  = 0._DP

        EATPI   = 0._DP
        EMLPI   = 0._DP
        EIOPI   = 0._DP
        EPLPI   = 0._DP
        EATEI   = 0._DP
        EMLEI   = 0._DP
        EIOEI   = 0._DP
        EPLEI   = 0._DP

        MODCOL  = 0
        IESTCX  = 0
        IESTEL  = 0
        IESTPI  = 0
        IESTEI  = 0
        IESTPH  = 0

cdr  RC secondaries
        NATPRC  = 0
        NMLPRC  = 0
        NIOPRC  = 0
        NPLPRC  = 0
        NPHPRC  = 0
        NATPRC_2= 0
        NMLPRC_2= 0
        NIOPRC_2= 0
        NPLPRC_2= 0
        NPHPRC_2= 0

cdr  CX secondaries
        N1STX   = 0
        N2NDX   = 0

cdr  number of reaction of type PI,EI,CX,EL,RC,...
cdr  These will be counted in setamd.f and subprograms
        NRPII   = 0
        NREII   = 0
        NRCXI   = 0
        NRELI   = 0
        NRRCI   = 0
        NRPHI   = 0
        NROTI   = 0
cdr  number of bgk reactions (nonlinear mode), counted in setamd.f and subprograms
        NRBGI   = 0

        NSEPRC  = 0

        NREACX  = 0
        NREAPI  = 0
        NREAEL  = 0
        NREAEI  = 0
        NREARC  = 0
        NREAPH  = 0

        NELREI  = 0
        JELREI  = 0
        NHVREI  = 0
        NELREL  = 0
        NELRRC  = 0
        JELRRC  = 0
        NELRPI  = 0
        JELRPI  = 0
        NHVRPI  = 0
        NPLRPI  = 0
        NELRCX  = 0
        NELRPH  = 0

        NREACT  = 0

        IPATEI  = 0
        IPMLEI  = 0
        IPIOEI  = 0
        IPPLEI  = 0
        IPATPI  = 0
        IPMLPI  = 0
        IPIOPI  = 0
        IPPLPI  = 0

        LGACX   = 0
        LGMCX   = 0
        LGICX   = 0

        LGAEI   = 0
        LGMEI   = 0
        LGIEI   = 0

        LGAEL   = 0
        LGMEL   = 0
        LGIEL   = 0

        LGPRC   = 0

        LGAPI   = 0
        LGMPI   = 0
        LGIPI   = 0

        LGPHPH  = 0
        LGPPH   = 0

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_CMDTA


      SUBROUTINE EIRENE_WRITE_CMDTA
cdr  read and write A&M data onto fort.13, controlled by NFILEL option (input block 1)

      WRITE (13+IFOFF)
     . TABEI1 ,TABRC1 ,TABPI3 ,TABCX3 ,TABEL3 ,
     . FDLMPI ,FDLMCX ,FDLMEL ,
     . ADDPI  ,ADDCX  ,ADDEL  ,
     . FACRRC ,FACRPI ,FACREL ,FACREI ,FACRCX ,

     . PELEI  ,PATEI  ,PMLEI  ,PIOEI  ,PPLEI  ,
     . PELPI  ,PATPI  ,PMLPI  ,PIOPI  ,PPLPI  ,
     . P2ND   ,P2NP   ,P2NEI  ,P2NPI  ,

     . EELEI1 ,EELRC1 ,EELPI3 ,
     . EHVEI1 ,EHVPI3 ,
     . EPLPI3 ,EPLCX3 ,EPLEL3 ,EPLPH3 ,

     . EATPI  ,EMLPI  ,EIOPI  ,EPLPI  ,
     . EATEI  ,EMLEI  ,EIOEI  ,EPLEI

#ifdef CHECKBIN
      write (113,*) 'CMDTA TABEI1 ...'
      WRITE (113,*)
     . TABEI1 ,TABRC1 ,TABPI3 ,TABCX3 ,TABEL3 ,
     . FDLMPI ,FDLMCX ,FDLMEL ,
     . ADDPI  ,ADDCX  ,ADDEL  ,
     . FACRRC ,FACRPI ,FACREL ,FACREI ,FACRCX ,

     . PELEI  ,PATEI  ,PMLEI  ,PIOEI  ,PPLEI  ,
     . PELPI  ,PATPI  ,PMLPI  ,PIOPI  ,PPLPI  ,
     . P2ND   ,P2NP   ,P2NEI  ,P2NPI  ,

     . EELEI1 ,EELRC1 ,EELPI3 ,
     . EHVEI1 ,EHVPI3 ,
     . EPLPI3 ,EPLCX3 ,EPLEL3 ,EPLPH3 ,

     . EATPI  ,EMLPI  ,EIOPI  ,EPLPI  ,
     . EATEI  ,EMLEI  ,EIOEI  ,EPLEI
#endif

      WRITE (13+IFOFF)
     . MODCOL ,IESTCX ,IESTEL ,IESTPI ,IESTEI ,
     . NAEII  ,NMEII  ,NIEII  ,NPHEII  ,
     . NACXI  ,NMCXI  ,NICXI  ,NPHCXI  ,
     . NAELI  ,NMELI  ,NIELI  ,NPHELI  ,
     . NAPII  ,NMPII  ,NIPII  ,NPHPII  ,
     . NPRCI  ,
     . NAEIIM ,NMEIIM ,NIEIIM ,NPHEIIM ,
     . NACXIM ,NMCXIM ,NICXIM ,NPHCXIM ,
     . NAELIM ,NMELIM ,NIELIM ,NPHELIM ,
     . NAPIIM ,NMPIIM ,NIPIIM ,NPHPIIM ,
     . NPRCIM ,
     . NPBGKA ,NPBGKM ,NPBGKI ,NPBGKP ,
     . NATPRC ,NMLPRC ,NIOPRC ,NPLPRC ,NPHPRC ,
     . NATPRC_2 ,NMLPRC_2 ,NIOPRC_2 ,NPLPRC_2 ,NPHPRC_2 ,
     . N1STX  ,N2NDX  ,

     . NRPII  ,NREII  ,NRCXI  ,NRELI  ,NRRCI  ,NRPHI  ,NRBGI  ,

     . NSEACX ,NSEMCX ,NSEICX ,NSEAEL ,NSEMEL ,NSEIEL ,NSEPRC ,
     . NREACX ,NREAPI ,NREAEL ,NREAEI ,NREARC ,NREAPH ,
     . NELREI ,JELREI ,NHVREI ,NELREL ,
     . NELRRC ,JELRRC ,
     . NELRPI ,JELRPI ,
     . NELRCX ,NHVRPI ,NPLRPI ,
     . NELRPH ,NREACT ,
     . IPATEI ,IPMLEI ,IPIOEI ,IPPLEI ,IPATPI ,IPMLPI ,IPIOPI ,IPPLPI ,
     . LGACX  ,LGMCX  ,LGICX  ,
     . LGAEI  ,LGMEI  ,LGIEI  ,
     . LGAEL  ,LGMEL  ,LGIEL  ,
     . LGPRC  ,
     . LGAPI  ,LGMPI  ,LGIPI
#ifdef CHECKBIN
      write (113,*) 'CMDTA MODCOL ...'
      WRITE (113,*)
     . MODCOL ,IESTCX ,IESTEL ,IESTPI ,IESTEI ,
     . NAEII  ,NMEII  ,NIEII  ,NPHEII  ,
     . NACXI  ,NMCXI  ,NICXI  ,NPHCXI  ,
     . NAELI  ,NMELI  ,NIELI  ,NPHELI  ,
     . NAPII  ,NMPII  ,NIPII  ,NPHPII  ,
     . NPRCI  ,
     . NAEIIM ,NMEIIM ,NIEIIM ,NPHEIIM ,
     . NACXIM ,NMCXIM ,NICXIM ,NPHCXIM ,
     . NAELIM ,NMELIM ,NIELIM ,NPHELIM ,
     . NAPIIM ,NMPIIM ,NIPIIM ,NPHPIIM ,
     . NPRCIM ,
     . NPBGKA ,NPBGKM ,NPBGKI ,NPBGKP ,
     . NATPRC ,NMLPRC ,NIOPRC ,NPLPRC ,NPHPRC ,
     . NATPRC_2 ,NMLPRC_2 ,NIOPRC_2 ,NPLPRC_2 ,NPHPRC_2 ,
     . N1STX  ,N2NDX  ,

     . NRPII  ,NREII  ,NRCXI  ,NRELI  ,NRRCI  ,NRPHI  ,NRBGI  ,

     . NSEACX ,NSEMCX ,NSEICX ,NSEAEL ,NSEMEL ,NSEIEL ,NSEPRC ,
     . NREACX ,NREAPI ,NREAEL ,NREAEI ,NREARC ,NREAPH ,
     . NELREI ,JELREI ,NHVREI ,NELREL ,
     . NELRRC ,JELRRC ,
     . NELRPI ,JELRPI ,
     . NELRCX ,NHVRPI ,NPLRPI ,
     . NELRPH ,
     . NREACT ,
     . IPATEI ,IPMLEI ,IPIOEI ,IPPLEI ,IPATPI ,IPMLPI ,IPIOPI ,IPPLPI ,
     . LGACX  ,LGMCX  ,LGICX  ,
     . LGAEI  ,LGMEI  ,LGIEI  ,
     . LGAEL  ,LGMEL  ,LGIEL  ,
     . LGPRC  ,
     . LGAPI  ,LGMPI  ,LGIPI
#endif

      RETURN
      END SUBROUTINE EIRENE_WRITE_CMDTA


      SUBROUTINE EIRENE_READ_CMDTA

      READ (13+IFOFF)
     . TABEI1 ,TABRC1 ,TABPI3 ,TABCX3 ,TABEL3 ,
     . FDLMPI ,FDLMCX ,FDLMEL ,
     . ADDPI  ,ADDCX  ,ADDEL  ,
     . FACRRC ,FACRPI ,FACREL ,FACREI ,FACRCX ,

     . PELEI  ,PATEI  ,PMLEI  ,PIOEI  ,PPLEI  ,
     . PELPI  ,PATPI  ,PMLPI  ,PIOPI  ,PPLPI  ,
     . P2ND   ,P2NP   ,P2NEI  ,P2NPI  ,

     . EELEI1 ,EELRC1 ,EELPI3 ,
     . EHVEI1 ,EHVPI3 ,
     . EPLPI3 ,EPLCX3 ,EPLEL3 ,EPLPH3 ,

     . EATPI  ,EMLPI  ,EIOPI  ,EPLPI  ,
     . EATEI  ,EMLEI  ,EIOEI  ,EPLEI

      READ (13+IFOFF)
     . MODCOL ,IESTCX ,IESTEL ,IESTPI ,IESTEI ,
     . NAEII  ,NMEII  ,NIEII  ,NPHEII  ,
     . NACXI  ,NMCXI  ,NICXI  ,NPHCXI  ,
     . NAELI  ,NMELI  ,NIELI  ,NPHELI  ,
     . NAPII  ,NMPII  ,NIPII  ,NPHPII  ,
     . NPRCI  ,
     . NAEIIM ,NMEIIM ,NIEIIM ,NPHEIIM  ,
     . NACXIM ,NMCXIM ,NICXIM ,NPHCXIM  ,
     . NAELIM ,NMELIM ,NIELIM ,NPHELIM  ,
     . NAPIIM ,NMPIIM ,NIPIIM ,NPHPIIM  ,
     . NPRCIM ,
     . NPBGKA ,NPBGKM ,NPBGKI ,NPBGKP ,
     . NATPRC ,NMLPRC ,NIOPRC ,NPLPRC ,NPHPRC ,
     . NATPRC_2 ,NMLPRC_2 ,NIOPRC_2 ,NPLPRC_2 ,NPHPRC_2 ,
     . N1STX  ,N2NDX  ,

     . NRPII  ,NREII  ,NRCXI  ,NRELI  ,NRRCI  ,NRPHI  ,NRBGI  ,

     . NSEACX ,NSEMCX ,NSEICX ,NSEAEL ,NSEMEL ,NSEIEL ,NSEPRC ,
     . NREACX ,NREAPI ,NREAEL ,NREAEI ,NREARC ,NREAPH ,
     . NELREI ,JELREI ,NHVREI ,NELREL ,
     . NELRRC ,JELRRC ,
     . NELRPI ,JELRPI ,
     . NELRCX ,NHVRPI ,NPLRPI ,
     . NELRPH ,
     . NREACT ,
     . IPATEI ,IPMLEI ,IPIOEI ,IPPLEI ,IPATPI ,IPMLPI ,IPIOPI ,IPPLPI ,
     . LGACX  ,LGMCX  ,LGICX  ,
     . LGAEI  ,LGMEI  ,LGIEI  ,
     . LGAEL  ,LGMEL  ,LGIEL  ,
     . LGPRC  ,
     . LGAPI  ,LGMPI  ,LGIPI

      RETURN
      END SUBROUTINE EIRENE_READ_CMDTA


      SUBROUTINE EIRENE_WRITE_CMAMF

      INTEGER :: IR

      WRITE (13+IFOFF)
     . DELPOT, FACREA,
     . FREACA, FREACM, FREACI, FREACP, FREACPH,
     . EDPOTA, EDPOTM, EDPOTI, EDPOTP, EDPOTPH,
     . EELECA, EELECM, EELECI, EELECP, EELECPH,
     . EBULKA, EBULKM, EBULKI, EBULKP, EBULKPH,
     . ESCD1A, ESCD1M, ESCD1I, ESCD1P, ESCD1PH,

     . NREACI, ISWR,   MODCLF, MASSP,  MASST,  IFTFLG,
     . NRCP,   NRCA,   NRCM,   NRCI,   NRCPH,
     . IREACA, IREACM, IREACI, IREACP, IREACPH,
     . IBULKA, IBULKM, IBULKI, IBULKP, IBULKPH,
     . ISCD1A, ISCD1M, ISCD1I, ISCD1P, ISCD1PH,
     . ISCD2A, ISCD2M, ISCD2I, ISCD2P, ISCD2PH,
     . ISCD3A, ISCD3M, ISCD3I, ISCD3P, ISCD3PH,
     . ISCD4A, ISCD4M, ISCD4I, ISCD4P, ISCD4PH,
     . ISCDEA, ISCDEM, ISCDEI, ISCDEP, ISCDEPH,
     . IESTMA, IESTMM, IESTMI, IESTMPH,
     . IBGKA,  IBGKM,  IBGKI,  IBGKPH

#ifdef CHECKBIN
      WRITE (113,*) 'WRITE_CMAMF DELPOT ...'
      WRITE (113,*)
     . DELPOT, FACREA,
     . FREACA, FREACM, FREACI, FREACP, FREACPH,
     . EDPOTA, EDPOTM, EDPOTI, EDPOTP, EDPOTPH,
     . EELECA, EELECM, EELECI, EELECP, EELECPH,
     . EBULKA, EBULKM, EBULKI, EBULKP, EBULKPH,
     . ESCD1A, ESCD1M, ESCD1I, ESCD1P, ESCD1PH,

     . NREACI, ISWR,   MODCLF, MASSP,  MASST,  IFTFLG,
     . NRCP,   NRCA,   NRCM,   NRCI, NRCPH,
     . IREACA, IREACM, IREACI, IREACP, IREACPH,
     . IBULKA, IBULKM, IBULKI, IBULKP, IBULKPH,
     . ISCD1A, ISCD1M, ISCD1I, ISCD1P, ISCD1PH,
     . ISCD2A, ISCD2M, ISCD2I, ISCD2P, ISCD2PH,
     . ISCD3A, ISCD3M, ISCD3I, ISCD3P, ISCD3PH,
     . ISCD4A, ISCD4M, ISCD4I, ISCD4P, ISCD4PH,
     . ISCDEA, ISCDEM, ISCDEI, ISCDEP, ISCDEPH,
     . IESTMA, IESTMM, IESTMI, IESTMPH, IBGKA , IBGKM , IBGKI, IBGKPH

      WRITE (113,*) 'WRITE_CMAMF REACDAT ...'
#endif
!  reacdat(-11:nreac)%.... minimal model is not written
      DO IR=1,NREACI
        WRITE (13+IFOFF)
     .             REACDAT(IR)%LPOT,
     .             REACDAT(IR)%LCRS,
     .             REACDAT(IR)%LRTC,
     .             REACDAT(IR)%LRTCMW,
     .             REACDAT(IR)%LRTCEW,
c
     .             REACDAT(IR)%LOTH,REACDAT(IR)%LPHR,
c
     .             REACDAT(IR)%ETH,
     .             REACDAT(IR)%KER,
     .             REACDAT(IR)%DELP,
     .             REACDAT(IR)%RTMAX,
     .             REACDAT(IR)%ERTMAX,
     .             REACDAT(IR)%EARRH0,
     .             REACDAT(IR)%EARRH1,
     .             REACDAT(IR)%NOSEC
#ifdef CHECKBIN
        WRITE (113,*) 'IR =',IR
        WRITE (113,*)
     .             REACDAT(IR)%LPOT,
     .             REACDAT(IR)%LCRS,
     .             REACDAT(IR)%LRTC,
     .             REACDAT(IR)%LRTCMW,
     .             REACDAT(IR)%LRTCEW,
c
     .             REACDAT(IR)%LOTH,REACDAT(IR)%LPHR,
c
     .             REACDAT(IR)%ETH,
     .             REACDAT(IR)%KER,
     .             REACDAT(IR)%DELP,
     .             REACDAT(IR)%RTMAX,
     .             REACDAT(IR)%ERTMAX,
     .             REACDAT(IR)%EARRH0,
     .             REACDAT(IR)%EARRH1,
     .             REACDAT(IR)%NOSEC
#endif
c
        IF (REACDAT(IR)%LPOT)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%POT)
        IF (REACDAT(IR)%LCRS)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%CRS)
        IF (REACDAT(IR)%LRTC)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTC)
        IF (REACDAT(IR)%LRTCMW)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTCMW)
        IF (REACDAT(IR)%LRTCEW)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%RTCEW)
        IF (REACDAT(IR)%LOTH)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%OTH)
        IF (REACDAT(IR)%LPHR)
     .    CALL EIRENE_WRITE_FIT_FORM(REACDAT(IR)%PHR)
      END DO

      RETURN

      CONTAINS

        SUBROUTINE EIRENE_WRITE_FIT_FORM (RP)
        TYPE(FIT_FORMS),POINTER :: RP  ! =reacdat(ir)%xxx%

        WRITE (13+IFOFF) RP%IFIT

        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE SHAPE, BROADENING. SPECIAL FORMAT
          WRITE (13+IFOFF) RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .               RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .               RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .               RP%LINE%B12, RP%LINE%B21, RP%LINE%C6A
          WRITE (13+IFOFF) RP%LINE%IGND, RP%LINE%IRCART,
     .               RP%LINE%IPROFILETYPE,
     .               RP%LINE%IFREMD, RP%LINE%NRJPRT, RP%LINE%IPLSC6,
     .               RP%LINE%IMESS
          WRITE (13+IFOFF) RP%LINE%REACNAME, RP%LINE%KENN
#ifdef CHECKBIN
          WRITE (113,*) 'PHOTON LINE'
          WRITE (113,*) RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .               RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .               RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .               RP%LINE%B12, RP%LINE%B21, RP%LINE%C6A
          WRITE (113,*) RP%LINE%IGND, RP%LINE%IRCART,
     .               RP%LINE%IPROFILETYPE,
     .               RP%LINE%IFREMD, RP%LINE%NRJPRT, RP%LINE%IPLSC6,
     .               RP%LINE%IMESS
          WRITE (113,*) RP%LINE%REACNAME, RP%LINE%KENN
#endif

        ELSE IF (1<=RP%IFIT .AND. RP%IFIT <= 2) THEN
! DATA FOR FIT EXPRESSIONS (E.G. POLYNOMIAL, IN CASE OF HYDHEL DATABASE)
          WRITE (13+IFOFF) UBOUND(RP%POLY%DBLPOL)
          WRITE (13+IFOFF)        RP%POLY%DBLPOL
#ifdef CHECKBIN
          WRITE (113,*) 'POLY DATA'
          WRITE (113,*) UBOUND(RP%POLY%DBLPOL)
          WRITE (113,*)        RP%POLY%DBLPOL
#endif

        ELSE IF (RP%IFIT == 3) THEN
! DATA FOR 2D TABULATED A&M ENTRIES (2-parameter tables, e.g. ADAS)
          WRITE (13+IFOFF) RP%ADAS%NDENS,RP%ADAS%NTEMP
          WRITE (13+IFOFF) RP%ADAS%DENS,RP%ADAS%TEMP,RP%ADAS%TAB2D,
     .                     RP%ADAS%DDE,RP%ADAS%DTE
#ifdef CHECKBIN
          WRITE (113,*) 'ADAS DATA'
          WRITE (113,*) RP%ADAS%NDENS,RP%ADAS%NTEMP
          WRITE (113,*) RP%ADAS%DENS,RP%ADAS%TEMP,RP%ADAS%TAB2D,
     .                     RP%ADAS%DDE,RP%ADAS%DTE
#endif
        ELSE IF (RP%IFIT == 4) THEN
cdr former "HYDKIN" interface data type.
cdr HYDKIN is offline since early 2018.
cdr Transform this structure into a more general
cdr 1D tabular data format.
cdr Currently anticipated: only reaction rates vs. Te, Ti (H.2).
! DATA FOR 1D TABULATED A&M ENTRIES (single parameter tables)
          WRITE (13+IFOFF) RP%TAB1D%NTEMPS
          WRITE (13+IFOFF) RP%TAB1D%TEMPS,RP%TAB1D%TAB1D,
cdr precomputed difference-quotient, for interpolation
     .                     RP%TAB1D%DIFFQUOT,
     .                     RP%TAB1D%REAC_STRING,RP%TAB1D%REACNAME,
     .                     RP%TAB1D%RPRT
#ifdef CHECKBIN
          WRITE (113,*) '1D TAB DATA'
          WRITE (113,*) RP%TAB1D%NTEMPS
          WRITE (113,*) RP%TAB1D%TEMPS,RP%TAB1D%TAB1D,
     .                     RP%TAB1D%DIFFQUOT,
     .                     RP%TAB1D%REAC_STRING,RP%TAB1D%REACNAME,
     .                     RP%TAB1D%RPRT
#endif
        ELSE IF (RP%IFIT == 5) THEN
! DATA FOR COLLISIONAL-RADIATIVE MODEL A&M ENTRIES
          WRITE (13+IFOFF) RP%CRM%IFLAV,
     .                     RP%CRM%IVARST,
     .                     RP%CRM%IROW_ESC,  RP%CRM%ICOL_ESC
          WRITE (13+IFOFF) RP%CRM%POP_ESC
#ifdef CHECKBIN
          WRITE (113,*) 'COLRAD DATA'
          WRITE (113,*) RP%CRM%IFLAV,
     .                     RP%CRM%IVARST,
     .                     RP%CRM%IROW_ESC,  RP%CRM%ICOL_ESC
          WRITE (113,*) RP%CRM%POP_ESC
#endif
        ELSE

        END IF

cdr options for extrapolation from data tables or from validity range of fits.
        WRITE (13+IFOFF) RP%RC1MIN, RP%RC1MAX, RP%RC2MIN, RP%RC2MAX,
     .                   RP%FP1L(3), RP%FP1R(3), RP%FP2B(3), RP%FP2T(3),
     .                   RP%JFEX1MN, RP%JFEX1MX, RP%JFEX2MN, RP%JFEX2MX
#ifdef CHECKBIN
        WRITE (113,*) 'RANGES'
        WRITE (113,*) RP%RC1MIN, RP%RC1MAX, RP%RC2MIN, RP%RC2MAX,
     .                   RP%FP1L(3), RP%FP1R(3), RP%FP2B(3), RP%FP2T(3),
     .                   RP%JFEX1MN, RP%JFEX1MX, RP%JFEX2MN, RP%JFEX2MX
#endif

        RETURN
        END SUBROUTINE EIRENE_WRITE_FIT_FORM

      END SUBROUTINE EIRENE_WRITE_CMAMF


      SUBROUTINE EIRENE_READ_CMAMF

      INTEGER :: IR

      READ (13+IFOFF)
     . DELPOT, FACREA,
     . FREACA, FREACM, FREACI, FREACP, FREACPH,
     . EDPOTA, EDPOTM, EDPOTI, EDPOTP, EDPOTPH,
     . EELECA, EELECM, EELECI, EELECP, EELECPH,
     . EBULKA, EBULKM, EBULKI, EBULKP, EBULKPH,
     . ESCD1A, ESCD1M, ESCD1I, ESCD1P, ESCD1PH,

     . NREACI, ISWR,   MODCLF, MASSP,  MASST,  IFTFLG,
     . NRCP,   NRCA,   NRCM,   NRCI,   NRCPH,
     . IREACA, IREACM, IREACI, IREACP, IREACPH,
     . IBULKA, IBULKM, IBULKI, IBULKP, IBULKPH,
     . ISCD1A, ISCD1M, ISCD1I, ISCD1P, ISCD1PH,
     . ISCD2A, ISCD2M, ISCD2I, ISCD2P, ISCD2PH,
     . ISCD3A, ISCD3M, ISCD3I, ISCD3P, ISCD3PH,
     . ISCD4A, ISCD4M, ISCD4I, ISCD4P, ISCD4PH,
     . ISCDEA, ISCDEM, ISCDEI, ISCDEP, ISCDEPH,
     . IESTMA, IESTMM, IESTMI, IESTMPH,
     . IBGKA , IBGKM , IBGKI,  IBGKPH

!  reacdat(-11:nreac)%.... minimal model is not read
      DO IR=1,NREACI
        READ (13+IFOFF)
     .            REACDAT(IR)%LPOT,
     .            REACDAT(IR)%LCRS,
     .            REACDAT(IR)%LRTC,
     .            REACDAT(IR)%LRTCMW,
     .            REACDAT(IR)%LRTCEW,
     .            REACDAT(IR)%LOTH,
     .            REACDAT(IR)%LPHR,

     .            REACDAT(IR)%ETH,
     .            REACDAT(IR)%KER,
     .            REACDAT(IR)%DELP,
     .            REACDAT(IR)%RTMAX,
     .            REACDAT(IR)%ERTMAX,
     .            REACDAT(IR)%EARRH0,
     .            REACDAT(IR)%EARRH1,
     .            REACDAT(IR)%NOSEC

        IF (REACDAT(IR)%LPOT) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%POT)
        IF (REACDAT(IR)%LCRS) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%CRS)
        IF (REACDAT(IR)%LRTC) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTC)
        IF (REACDAT(IR)%LRTCMW)
     .    CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTCMW)
        IF (REACDAT(IR)%LRTCEW)
     .    CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%RTCEW)
        IF (REACDAT(IR)%LOTH) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%OTH)
        IF (REACDAT(IR)%LPHR) CALL EIRENE_READ_FIT_FORM(REACDAT(IR)%PHR)
      END DO

      RETURN

      CONTAINS

        SUBROUTINE EIRENE_READ_FIT_FORM (RP)
!  RP = REACDAT(IR)%A... IN CALLING PROGRAM, A= CRS, RTC, ....
        TYPE(FIT_FORMS),POINTER :: RP
        INTEGER :: ND1, ND2,  ! table range for POLY
     .             ND, NT     ! table range for TAB1D, TAB2D

        READ (13+IFOFF) RP%IFIT

        IF (RP%IFIT < 0) THEN
! DATA FOR PHOTONIC LINE
cdr  most of this stuff: obsolete, from lighting applications,
cdr  high pressure gas discharge lamps, around 2002.. should not be here!
          IF (.NOT.ASSOCIATED(RP%LINE)) ALLOCATE (RP%LINE)
          READ (13+IFOFF) RP%LINE%E0, RP%LINE%E1, RP%LINE%AIK,
     .              RP%LINE%G1, RP%LINE%G2, RP%LINE%C2,
     .              RP%LINE%C3, RP%LINE%C4, RP%LINE%C6,
     .              RP%LINE%B12, RP%LINE%B21, RP%LINE%C6A
          READ (13+IFOFF) RP%LINE%IGND, RP%LINE%IRCART,
     .              RP%LINE%IPROFILETYPE,
     .              RP%LINE%IFREMD, RP%LINE%NRJPRT, RP%LINE%IPLSC6,
     .              RP%LINE%IMESS
          READ (13+IFOFF) RP%LINE%REACNAME, RP%LINE%KENN

        ELSE IF (1<= RP%IFIT .AND. RP%IFIT <= 2) THEN
! DATA FOR FIT EXPRESSIONS (POLYNOMIAL, E.G. IN CASE OF HYDHEL, AMJUEL DATABASE)
          IF (.NOT.ASSOCIATED(RP%POLY)) ALLOCATE (RP%POLY)
          IF (ASSOCIATED(RP%POLY%DBLPOL)) DEALLOCATE (RP%POLY%DBLPOL)

          READ (13+IFOFF) ND1,ND2
          ALLOCATE (RP%POLY%DBLPOL(ND1,ND2))
          READ (13+IFOFF) RP%POLY%DBLPOL

        ELSE IF (RP%IFIT == 3) THEN
! DATA FOR 2D TABULATED A&M ENTRIES (2-parameter tables, e.g. ADAS)
          READ (13+IFOFF) RP%ADAS%NDENS,RP%ADAS%NTEMP
          ND = RP%ADAS%NDENS
          NT = RP%ADAS%NTEMP
          ALLOCATE (RP%ADAS%DENS(ND))
          ALLOCATE (RP%ADAS%TEMP(NT))
          ALLOCATE (RP%ADAS%TAB2D(NT,ND))
          ALLOCATE (RP%ADAS%DDE(ND))
          ALLOCATE (RP%ADAS%DTE(NT))
          READ (13+IFOFF) RP%ADAS%DENS,RP%ADAS%TEMP,RP%ADAS%TAB2D,
     .                    RP%ADAS%DDE, RP%ADAS%DTE

        ELSE IF (RP%IFIT == 4) THEN
! DATA FOR 1D TABULATED A&M ENTRIES (single parameter table)
cdr unfinished and untested.
          READ (13+IFOFF) RP%TAB1D%NTEMPS
          NT = RP%TAB1D%NTEMPS
          ALLOCATE (RP%TAB1D%TEMPS(NT))
          ALLOCATE (RP%TAB1D%TAB1D(NT))
          ALLOCATE (RP%TAB1D%DIFFQUOT(NT))
          READ (13+IFOFF) RP%TAB1D%TEMPS,RP%TAB1D%TAB1D,
     .               RP%TAB1D%DIFFQUOT,
     .               RP%TAB1D%REAC_STRING,RP%TAB1D%REACNAME,
     .               RP%TAB1D%RPRT
cdr: tbd: set difference quotient in table: tab1d vs. temps

        ELSE IF (RP%IFIT == 5) THEN
! DATA FOR COLLISIONAL-RADIATIVE MODEL A&M ENTRIES
cdr may 2020 added: m_popesc,  and pop_esc arrays, rather scalars.
          IF (.NOT.ASSOCIATED(RP%CRM)) ALLOCATE (RP%CRM)
          READ (13+IFOFF) RP%CRM%IFLAV,
     .                    RP%CRM%IVARST,
     .                    RP%CRM%IROW_ESC,  RP%CRM%ICOL_ESC
          READ (13+IFOFF) RP%CRM%POP_ESC
        ELSE

        END IF

        READ (13+IFOFF) RP%RC1MIN, RP%RC1MAX, RP%RC2MIN, RP%RC2MAX,
     .                  RP%FP1L(3), RP%FP1R(3), RP%FP2B(3), RP%FP2T(3),
     .                  RP%JFEX1MN, RP%JFEX1MX, RP%JFEX2MN, RP%JFEX2MX

        RETURN
        END SUBROUTINE EIRENE_READ_FIT_FORM

      END SUBROUTINE EIRENE_READ_CMAMF


      SUBROUTINE EIRENE_GET_REACTION (IR)
cdr ?  still needed? only for "line"?

      INTEGER, INTENT(IN) :: IR

      REACTION => REACDAT(IR)%PHR%LINE

      IDREAC = IR

      RETURN
      END SUBROUTINE EIRENE_GET_REACTION


      SUBROUTINE EIRENE_SET_REACTION_DATA
     .           (IR,ISW,IFTFL,RDATA,INEP,KNEP,
     .            IUNOUT,LTEST,
c  from here on: optional input parameters
     .            RC1MIN, RC1MAX, FP1, JFEX1MN, JFEX1MX,
     .            RC2MIN, RC2MAX, FP2, JFEX2MN, JFEX2MX,
     .            RTMAX, ERTMAX, ETH, KER, DELP, EARRH0, EARRH1)

c  set reaction data structure REACDAT, for reaction no. IR.
c  here only: 1D or 2D polygonial fits for reaction data.
c              RDATA --> REA, and then: REACDAT(IR)%...%POLY => REA
c  and:                         NULLIFY REACDAT(IR)%...%ADAS
c  and:                         NULLIFY REACDAT(IR)%...%LINE
c  and:                         NULLIFY REACDAT(IR)%...%TAB1D
c  and:                         NULLIFY REACDAT(IR)%...%CRM
c
c  1) called from READ_PHOTDBK
c  2) called from SLREAC, option "CONST"
c  3) called from SLREAC, option AMJUEL, HYDHEL, H2VIBR, METHAN, AMMONX
cpb   not called from READ_COLRAD
cpb   not called from READ_TAB2D
cpb   not called from READ_POLY

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IR, ISW, IFTFL, INEP, KNEP, IUNOUT
      INTEGER, OPTIONAL, INTENT(IN) :: JFEX1MN, JFEX1MX,JFEX2MN, JFEX2MX
      REAL(DP), INTENT(IN) :: RDATA(9,*)
      REAL(DP), OPTIONAL, INTENT(IN) :: RC1MIN, RC1MAX, FP1(6),
     .                                  RC2MIN, RC2MAX, FP2(6),
     .                                  RTMAX,  ERTMAX, ETH, KER, DELP,
     .                                  EARRH0, EARRH1
      LOGICAL, INTENT(IN) :: LTEST
      INTEGER :: NDIM1, NDIM2, IFIT
      REAL(DP) :: CTEST
      LOGICAL :: DATA_FOUND
      TYPE(POLY_DATA), POINTER :: REA
      INTEGER, SAVE :: ISW2D(7) = (/ 3, 4, 6, 7, 9, 10, 12 /)
      EXTERNAL :: EIRENE_EXIT_OWN

cdr   NDIM1 = INEP
      NDIM1 = 9
      IF (MOD(IFTFL,100) == 10) NDIM1 = 1

      NDIM2 = 1
      IF (COUNT(ISW2D == ISW) > 0) NDIM2=9

      IF (LTEST) THEN
        CTEST = SUM(ABS(RDATA(1:NDIM1,1:NDIM2)))
        DATA_FOUND = CTEST.GT.1.E-30_DP
        DATA_FOUND = DATA_FOUND .OR.
     .   (PRESENT(RC1MIN) .and. PRESENT(RC1MAX) .and.
     .    RC1MIN.EQ.RC1MAX .and. RC1MIN.ne.0._DP .and. RC1MAX.ne.0._DP)
        DATA_FOUND = DATA_FOUND .OR.
     .   (PRESENT(RC1MAX) .and. RC1MAX.eq.0._DP .and.
     .    PRESENT(JFEX1MX) .and. PRESENT(FP1) .and.
     .  ((JFEX1MX.eq.2 .and. FP1(5).gt.0_DP) .or.
     .   (JFEX1MX.eq.3 .and.(FP1(4).gt.0_DP .or. FP1(5).gt.0_DP)) .or.
     .   (JFEX1MX.eq.5 .and. SUM(ABS(FP1(4:6))).GT.1.E-30_DP)))
        IF (.not.DATA_FOUND) THEN
           WRITE (iunout,*)
     .            'ERROR IN SUBROUTINE EIRENE_SET_REACTION_DATA:',
     .            ' ZERO FIT COEFFICIENTS'
           WRITE (iunout,*) 'IR = ',IR,'  EXIT CALLED!'
           CALL EIRENE_EXIT_OWN(1)
        END IF
      END IF

cdr the reaction parameters ETH, KER,....,ERTMAX
cdr are specific for a given reaction process.
cdr but then, within this physical process,
cdr must be one and the same for all data (pot, crs, rtc, ertc....)

CDR  warning: until now IR can. but need not necessarily be
      IF (PRESENT(ETH) .AND. (REACDAT(IR)%ETH == 0._DP))
     .    REACDAT(IR)%ETH = ETH
      IF (PRESENT(KER) .AND. (REACDAT(IR)%KER == 0._DP))
     .    REACDAT(IR)%KER = KER
      IF (PRESENT(DELP) .AND. (REACDAT(IR)%DELP == 0._DP))
     .    REACDAT(IR)%DELP = DELP
      IF (PRESENT(EARRH0) .AND. (REACDAT(IR)%EARRH0 == 0._DP))
     .    REACDAT(IR)%EARRH0 = EARRH0
      IF (PRESENT(EARRH1) .AND. (REACDAT(IR)%EARRH1 == 0._DP))
     .    REACDAT(IR)%EARRH1 = EARRH1
c  parameters for rejection sampling of bulk collision partner
c  The product sigma(vel) *vel takes ist max. value RTMAX at E=ERTMAX
      IF (PRESENT(RTMAX) .AND. (REACDAT(IR)%RTMAX == 0._DP))
     .    REACDAT(IR)%RTMAX = RTMAX
      IF (PRESENT(ERTMAX) .AND. (REACDAT(IR)%ERTMAX == -HUGE(1._DP)))
     .    REACDAT(IR)%ERTMAX = ERTMAX

      ALLOCATE (REA)
      ALLOCATE (REA%DBLPOL(1:NDIM1,1:NDIM2))

      REA%DBLPOL(1:NDIM1,1:NDIM2) = RDATA(1:NDIM1,1:NDIM2)
      IFIT = 1
      IF (NDIM2 == 9) IFIT = 2

      SELECT CASE(ISW)

      CASE (0)
        IF (REACDAT(IR)%LPOT) THEN
          WRITE (IUNOUT,*) ' POTENTIAL ALREADY SPECIFIED FOR REACTION',
     .                       IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for interaction potentials,
cdr                                            or differential scattering info.
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%POT)

        REACDAT(IR)%LPOT = .TRUE.
        REACDAT(IR)%POT%POLY => REA
        REACDAT(IR)%POT%IFIT = IFIT

cdr  set asymptotics for POT, if specified
cdr        rea => REACDAT(IR)%POT
cdr        call eirene_asympt_fit_form (rea)

        IF (PRESENT(RC1MIN)) REACDAT(IR)%POT%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%POT%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%POT%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%POT%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%POT%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%POT%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%POT%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%POT%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%POT%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%POT%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%POT%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%POT%JFEX2MX = JFEX2MX
c

      CASE (1)
        IF (REACDAT(IR)%LCRS) THEN
          WRITE (IUNOUT,*) ' CROSS-SECTION ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for (total) cross-sections
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%CRS)

        REACDAT(IR)%LCRS = .TRUE.
        REACDAT(IR)%CRS%POLY => REA
        REACDAT(IR)%CRS%IFIT = IFIT

c  extrapolation options
cdr  set asymptotics for CRS, if specified
cdr        rea => REACDAT(IR)%CRS
cdr        call eirene_asympt_fit_form (rea)
c

        IF (PRESENT(RC1MIN)) REACDAT(IR)%CRS%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%CRS%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%CRS%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%CRS%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%CRS%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%CRS%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%CRS%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%CRS%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%CRS%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%CRS%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%CRS%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%CRS%JFEX2MX = JFEX2MX
c

      CASE (2:4)
        IF (REACDAT(IR)%LRTC) THEN
          WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                     ' FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for reaction rate coefficients
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTC)

        REACDAT(IR)%LRTC = .TRUE.
        REACDAT(IR)%RTC%POLY => REA
        REACDAT(IR)%RTC%IFIT = IFIT

c  extrapolation options
c

        IF (PRESENT(RC1MIN)) REACDAT(IR)%RTC%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%RTC%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%RTC%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%RTC%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%RTC%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%RTC%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%RTC%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%RTC%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%RTC%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%RTC%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%RTC%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%RTC%JFEX2MX = JFEX2MX

      CASE (5:7)
        IF (REACDAT(IR)%LRTCMW) THEN
          WRITE (IUNOUT,*) ' MOMENTUM-WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for momentum-weighted rate coefficients
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCMW)

        REACDAT(IR)%LRTCMW = .TRUE.
        REACDAT(IR)%RTCMW%POLY => REA
        REACDAT(IR)%RTCMW%IFIT = IFIT

c  extrapolation options
c

        IF (PRESENT(RC1MIN)) REACDAT(IR)%RTCMW%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%RTCMW%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%RTCMW%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%RTCMW%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%RTCMW%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%RTCMW%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%RTCMW%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%RTCMW%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%RTCMW%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%RTCMW%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%RTCMW%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%RTCMW%JFEX2MX = JFEX2MX

      CASE (8:10)
        IF (REACDAT(IR)%LRTCEW) THEN
          WRITE (IUNOUT,*) ' ENERGY-WEIGHTED RATE COEFFICIENT',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for energy-weighted rate coefficients
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCEW)

        REACDAT(IR)%LRTCEW = .TRUE.
        REACDAT(IR)%RTCEW%POLY => REA
        REACDAT(IR)%RTCEW%IFIT = IFIT

c  extrapolation options
c

        IF (PRESENT(RC1MIN)) REACDAT(IR)%RTCEW%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%RTCEW%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%RTCEW%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%RTCEW%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%RTCEW%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%RTCEW%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%RTCEW%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%RTCEW%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%RTCEW%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%RTCEW%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%RTCEW%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%RTCEW%JFEX2MX = JFEX2MX
c

      CASE (11:12)
        IF (REACDAT(IR)%LOTH) THEN
          WRITE (IUNOUT,*) ' OTHER POLYNOMIAL FIT COEFFICIENTS',
     .                     ' ALREADY SPECIFIED FOR REACTION', IR
          DEALLOCATE (REA)
          IF (IR < 0) RETURN
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
        END IF

cdr  allocate, initialize, default asymptotics for other data (population coefficients,
cdr                        density ratios, etc...)
        CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%OTH)

        REACDAT(IR)%LOTH = .TRUE.
        REACDAT(IR)%OTH%POLY => REA
        REACDAT(IR)%OTH%IFIT = IFIT

c  extrapolation options

        IF (PRESENT(RC1MIN)) REACDAT(IR)%OTH%RC1MIN = RC1MIN
        IF (PRESENT(RC1MAX)) REACDAT(IR)%OTH%RC1MAX = RC1MAX
        IF (PRESENT(RC2MIN)) REACDAT(IR)%OTH%RC2MIN = RC2MIN
        IF (PRESENT(RC2MAX)) REACDAT(IR)%OTH%RC2MAX = RC2MAX
        IF (PRESENT(FP1)) REACDAT(IR)%OTH%FP1L = FP1(1:3)
        IF (PRESENT(FP1)) REACDAT(IR)%OTH%FP1R = FP1(4:6)
        IF (PRESENT(FP2)) REACDAT(IR)%OTH%FP2B = FP2(1:3)
        IF (PRESENT(FP2)) REACDAT(IR)%OTH%FP2T = FP2(4:6)
        IF (PRESENT(JFEX1MN)) REACDAT(IR)%OTH%JFEX1MN = JFEX1MN
        IF (PRESENT(JFEX1MX)) REACDAT(IR)%OTH%JFEX1MX = JFEX1MX
        IF (PRESENT(JFEX2MN)) REACDAT(IR)%OTH%JFEX2MN = JFEX2MN
        IF (PRESENT(JFEX2MX)) REACDAT(IR)%OTH%JFEX2MX = JFEX2MX
c

      CASE DEFAULT
        WRITE (IUNOUT,*) ' WRONG REACTION TYPE SPECIFIED '
        WRITE (IUNOUT,*) ' REACTION NO. ', IR
        WRITE (IUNOUT,*) ' REACTION TYPE H.', ISW
        CALL EIRENE_EXIT_OWN(1)
      END SELECT

      RETURN
      END SUBROUTINE EIRENE_SET_REACTION_DATA

CDR  NEXT THREE ROUTINES: PROBABLY NOT NEEDED ?
      FUNCTION EIRENE_IS_RTC_TAB2D (IREAC) RESULT(RES)
cdr  something special about ADAS ?  unused !

      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES

      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTC) THEN
        RES = REACDAT(IREAC)%RTC%IFIT == 3
      END IF

      RETURN
      END FUNCTION EIRENE_IS_RTC_TAB2D


      FUNCTION EIRENE_IS_RTCEW_TAB2D (IREAC) RESULT(RES)
cdr only used for adding/subtracting bremsstrahlung ? unused ?

      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES

      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTCEW) THEN
        RES = REACDAT(IREAC)%RTCEW%IFIT == 3
      END IF

      RETURN
      END FUNCTION EIRENE_IS_RTCEW_TAB2D


      FUNCTION EIRENE_IS_RTCMW_TAB2D (IREAC) RESULT(RES)
cdr  something special about adas ?  unused !

      INTEGER, INTENT(IN) :: IREAC
      LOGICAL :: RES

      RES = .FALSE.
      IF (REACDAT(IREAC)%LRTCMW) THEN
        RES = REACDAT(IREAC)%RTCMW%IFIT == 3
      END IF

      RETURN
      END FUNCTION EIRENE_IS_RTCMW_TAB2D


      SUBROUTINE EIRENE_FREE_REACDAT
cdr  called from dealloc_comxs: Free data structure REACDAT at the end of a run.

      type(fit_forms), pointer :: rea
      integer :: ir

      DO IR = -11, NREAC
c  interaction potentials, scattering angle information
        IF (REACDAT(IR)%LPOT) THEN
           rea => REACDAT(IR)%POT
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF
c  cross-sections
        IF (REACDAT(IR)%LCRS) THEN
           rea => REACDAT(IR)%CRS
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF
c  rate coefficients
        IF (REACDAT(IR)%LRTC) THEN
           rea => REACDAT(IR)%RTC
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF
c  momentum-weighted rate coefficients
        IF (REACDAT(IR)%LRTCMW) THEN
           rea => REACDAT(IR)%RTCMW
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF
c  energy-weighted rate coefficients
        IF (REACDAT(IR)%LRTCEW) THEN
           rea => REACDAT(IR)%RTCEW
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF
c  "other reaction", e.g. population rate coefficient, density ratio
        IF (REACDAT(IR)%LOTH) THEN
           rea => REACDAT(IR)%OTH
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

cdr probably obsolete:
c  photonic reaction
        IF (REACDAT(IR)%LPHR) THEN
           rea => REACDAT(IR)%PHR
           call eirene_free_fit_form (rea)
           deallocate (rea)
        END IF

      END DO

      deallocate (reacdat)

      RETURN
      END SUBROUTINE EIRENE_FREE_REACDAT


      SUBROUTINE EIRENE_FREE_FIT_FORM(RP)

      TYPE(FIT_FORMS),POINTER :: RP

      IF (ASSOCIATED(RP%POLY)) THEN
         DEALLOCATE (RP%POLY%DBLPOL)
         DEALLOCATE (RP%POLY)
         NULLIFY(RP%POLY)
      END IF

      IF (ASSOCIATED(RP%ADAS)) THEN
         DEALLOCATE (RP%ADAS%DENS)
         DEALLOCATE (RP%ADAS%TEMP)
         DEALLOCATE (RP%ADAS%TAB2D)
         DEALLOCATE (RP%ADAS%DDE)
         DEALLOCATE (RP%ADAS%DTE)
         DEALLOCATE (RP%ADAS)
         NULLIFY(RP%ADAS)
      END IF

      IF (ASSOCIATED(RP%LINE)) THEN
         DEALLOCATE (RP%LINE)
         NULLIFY(RP%LINE)
      END IF

      IF (ASSOCIATED(RP%TAB1D)) THEN
         DEALLOCATE (RP%TAB1D%TEMPS)
         DEALLOCATE (RP%TAB1D%TAB1D)
         DEALLOCATE (RP%TAB1D%DIFFQUOT)
         DEALLOCATE (RP%TAB1D)
         NULLIFY(RP%TAB1D)
      END IF

      IF (ASSOCIATED(RP%CRM)) THEN
         DEALLOCATE (RP%CRM)
         NULLIFY(RP%CRM)
      END IF

      RETURN
      END SUBROUTINE EIRENE_FREE_FIT_FORM


      SUBROUTINE EIRENE_ALLOC_FIT_FORM (RP)
cdr
c  RP = REACDAT(IR)%TYP%...   with TYP= POT, CRS, RTC, ...
C  Initialize pointers (NULLIFY) and
c  set default asymptotics, for reacdat(ir)%typ%...

      TYPE(FIT_FORMS),POINTER :: RP

      IF (.NOT.ASSOCIATED(RP)) THEN
        ALLOCATE (RP)
        NULLIFY (RP%POLY)
        NULLIFY (RP%ADAS)
        NULLIFY (RP%LINE)   ! this should not be here.
                            ! It is not a data fit type
        NULLIFY (RP%TAB1D)
        NULLIFY (RP%CRM)

c  default asymptotics
        RP%RC1MIN = -HUGE(1._DP)
        RP%RC1MAX =  HUGE(1._DP)
        RP%RC2MIN = -HUGE(1._DP)
        RP%RC2MAX =  HUGE(1._DP)
        RP%FP1L = 0._DP
        RP%FP1R = 0._DP
        RP%FP2B = 0._DP
        RP%FP2T = 0._DP
        RP%JFEX1MN = 0
        RP%JFEX1MX = 0
        RP%JFEX2MN = 0
        RP%JFEX2MX = 0
c
      END IF

      RETURN
      END SUBROUTINE EIRENE_ALLOC_FIT_FORM


      SUBROUTINE EIRENE_BROADCAST_COMXS(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER, IR
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_COMXS(1)
        CALL EIRENE_ALLOC_COMXS(2)
      END IF

      CALL MPI_BCAST (XSTOR,MSTOR1*MSTOR2,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (XSTORV,NSTORV,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ZMFPTHI,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TDGTEMX,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (TABEI1,NREI*NSTORDR,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TABRC1,NREC*NSTORDR,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TABPI3,NRPI*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TABCX3,NRCX*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TABEL3,NREL*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (FDLMPI,NRPI,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FDLMCX,NRCX,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FDLMEL,NREL,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (ADDPI,NRPI*NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ADDCX,NRCX*NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ADDEL,NREL*NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (FACRRC,NREC*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FACRPI,NRPI*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FACREL,NREL*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FACREI,NREI*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FACRCX,NRCX*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

c  EI post-collision species distribution
      CALL MPI_BCAST (PELEI,NREI,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PATEI,NREI*NATMP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PMLEI,NREI*NMOLP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PIOEI,NREI*NIONP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PPLEI,NREI*NPLSP,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (P2ND,NREI*NSPZP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (P2NEI,NREI,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

c  PI post-collision species distribution
      CALL MPI_BCAST (PELPI,NRPI,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PATPI,NRPI*NATMP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PMLPI,NRPI*NMOLP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PIOPI,NRPI*NIONP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (PPLPI,NRPI*NPLSP,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (P2NP,NRPI*NSPZP,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (P2NPI,NRPI,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

c  EI post-collision energetics
      CALL MPI_BCAST (EELEI1,NREI*NSTORDR,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EHVEI1,NREI*NSTORDR,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
c  RC post-collision energetics
      CALL MPI_BCAST (EELRC1,NREC*NSTORDR,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
c  PI post-collision energetics
      CALL MPI_BCAST (EELPI3,NRPI*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EHVPI3,NRPI*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EPLPI3,NRPI*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
c  CX post-collision energetics
      CALL MPI_BCAST (EPLCX3,NRCX*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
c  EL post-collision energetics
      CALL MPI_BCAST (EPLEL3,NREL*NSTORDR*NSTORDT,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (EATEI,NREI*NATMP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EMLEI,NREI*NMOLP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EIOEI,NREI*NIONP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EPLEI,NREI*NPLSP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (EATPI,NRPI*NATMP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EMLPI,NRPI*NMOLP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EIOPI,NRPI*NIONP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EPLPI,NRPI*NPLSP*2,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (MODCOL,7*5*MSTOR0,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (IESTCX,3*NRCX,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTEL,3*NREL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTPI,3*NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTEI,3*NREI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NRPII,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NREII,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCXI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRELI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRRCI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRPHI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NROTI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRBGI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NAEII,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMEII,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIEII,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NACXI,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMCXI,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NICXI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NAELI,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMELI,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIELI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NAPII,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMPII,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIPII,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHPHI,NPHOT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPPHI,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPRCI,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NAEIIM,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMEIIM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIEIIM,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NACXIM,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMCXIM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NICXIM,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NAELIM,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMELIM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIELIM,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NAPIIM,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMPIIM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIPIIM,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHPHIM,NPHOT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPPHIM,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPRCIM,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NPBGKA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPBGKM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPBGKI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPBGKP,2*NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NATPRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMLPRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIOPRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPLPRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHPRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NATPRC_2,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMLPRC_2,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIOPRC_2,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPLPRC_2,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHPRC_2,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (N1STX,3*NRCX,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (N2NDX,3*NRCX,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEACX,5*NATM*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEMCX,5*NMOL*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEICX,5*NION*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEAEL,5*NATM*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEMEL,5*NMOL*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEIEL,5*NION*NPLS,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSEPRC,5*NREC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (NREACX,NRCX,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELRCX,NRCX,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  FLAGS FOR PI TYPE REACTIONS:
      CALL MPI_BCAST (NREAPI,NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELRPI,NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (JELRPI,NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NHVRPI,NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPLRPI,NRPI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  FLAGS FOR EL TYPE REACTIONS:
      CALL MPI_BCAST (NREAEL,NREL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELREL,NREL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  FLAGS FOR RC TYPE REACTIONS:
      CALL MPI_BCAST (NREARC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELRRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (JELRRC,NREC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  FLAGS FOR EI TYPE REACTIONS:
      CALL MPI_BCAST (NREAEI,NREI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELREI,NREI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (JELREI,NREI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NHVREI,NREI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
cdr  FLAGS FOR PH TYPE REACTIONS:
      CALL MPI_BCAST (NREAPH,NRPH,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NELRPH,NRPH,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (IPATEI,NREI*NATMP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPMLEI,NREI*NMOLP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPIOEI,NREI*NIONP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPPLEI,NREI*NPLSP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPATPI,NRPI*NATMP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPMLPI,NRPI*NMOLP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPIOPI,NRPI*NIONP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IPPLPI,NRPI*NPLSP,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  CX
      CALL MPI_BCAST (LGACX,NATMP*(NRCX+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGMCX,NMOLP*(NRCX+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGICX,NIONP*(NRCX+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  EI
      CALL MPI_BCAST (LGAEI,NATMP*(NREI+1),MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGMEI,NMOLP*(NREI+1),MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGIEI,NIONP*(NREI+1),MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  EL
      CALL MPI_BCAST (LGAEL,NATMP*(NREL+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGMEL,NMOLP*(NREL+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGIEL,NIONP*(NREL+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  PI
      CALL MPI_BCAST (LGAPI,NATMP*(NRPI+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGMPI,NMOLP*(NRPI+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGIPI,NIONP*(NRPI+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  PH
      CALL MPI_BCAST (LGPPH,NPLSP*(NRPH+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGPHPH,NPHOTP*(NRPH+1)*2,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
c  RC
      CALL MPI_BCAST (LGPRC,NPLSP*(NREC+1),MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

cdr   old data structure CREAC has been replaced by more general data structure REACDAT

      DO IR=-11, NREAC
c interaction potential, differential cross-sections, etc.
        CALL MPI_BCAST (REACDAT(IR)%LPOT,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c total cross-sections
        CALL MPI_BCAST (REACDAT(IR)%LCRS,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c reaction rate coefficients
        CALL MPI_BCAST (REACDAT(IR)%LRTC,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c momentum-weighted rate coefficients
        CALL MPI_BCAST (REACDAT(IR)%LRTCMW,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c energy-weighted rate coefficients
        CALL MPI_BCAST (REACDAT(IR)%LRTCEW,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c other data, such as population coefficients, CR density ratios,....
        CALL MPI_BCAST (REACDAT(IR)%LOTH,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c ??
        CALL MPI_BCAST (REACDAT(IR)%LPHR,1,MPI_LOGICAL,
     .                  0,MPI_COMM_WORLD,ier)
c ??
        CALL MPI_BCAST (REACDAT(IR)%NOSEC,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)

c  reaction threshold, eV
        CALL MPI_BCAST (REACDAT(IR)%ETH,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (REACDAT(IR)%KER,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (REACDAT(IR)%DELP,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
c  max of sigma * v, cm**3/s
        CALL MPI_BCAST (REACDAT(IR)%RTMAX,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
c  collision energy, at which rtmax is taken, eV
        CALL MPI_BCAST (REACDAT(IR)%ERTMAX,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (REACDAT(IR)%EARRH0,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (REACDAT(IR)%EARRH1,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
      END DO

      DO IR = -11, NREAC
c  data for interaction potential
        IF (REACDAT(IR)%LPOT) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%POT)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%POT)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%POT,ME)
        END IF
c  data for cross-sections, cm**2
        IF (REACDAT(IR)%LCRS) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%CRS)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%CRS)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%CRS,ME)
        END IF
c  data for rate coefficients, cm**3/s
        IF (REACDAT(IR)%LRTC) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%RTC)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTC)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%RTC,ME)
        END IF
c  data for momentum-weighted rate coefficients, g cm/s cm**3/s
        IF (REACDAT(IR)%LRTCMW) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%RTCMW)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCMW)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%RTCMW,ME)
        END IF
c  data for energy-weighted rate coefficients, eV cm**3/s
        IF (REACDAT(IR)%LRTCEW) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%RTCEW)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCEW)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%RTCEW,ME)
        END IF
c  other data, such as CR population coefficients
        IF (REACDAT(IR)%LOTH) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%OTH)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%OTH)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%OTH,ME)
        END IF
c  data for photon line transport
        IF (REACDAT(IR)%LPHR) THEN
          IF (ME .NE. 0) THEN
            IF (.NOT.ASSOCIATED(REACDAT(IR)%PHR)) THEN
              CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%PHR)
            END IF
          END IF
          CALL EIRENE_BROAD_FIT_FORM(REACDAT(IR)%PHR,ME)
        END IF
      END DO  ! IR  -11,...,NREAC

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (FACREA,(NREAC+12)*2,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DELPOT,NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FREACA,NATM*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FREACM,NMOL*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FREACI,NION*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FREACP,NPLS*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (FREACPH,NPHOT*NREAC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EDPOTA,NATM*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EDPOTM,NMOL*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EDPOTI,NION*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EDPOTP,NPLS*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EDPOTPH,NPHOT*NREAC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EELECA,NATM*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EELECM,NMOL*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EELECI,NION*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EELECP,NPLS*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EELECPH,NPHOT*NREAC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EBULKA,NATM*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EBULKM,NMOL*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EBULKI,NION*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EBULKP,NPLS*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EBULKPH,NPHOT*NREAC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ESCD1A,NATM*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ESCD1M,NMOL*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ESCD1I,NION*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ESCD1P,NPLS*NREAC,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ESCD1PH,NPHOT*NREAC,MPI_REAL8,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (ISWR,NREAC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (MODCLF,NREAC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (MASSP,NREAC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (MASST,NREAC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IFTFLG,6*(NREAC+12),MPI_INTEGER,0,
     .                MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCP,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NRCPH,NPHOT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IREACA,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IREACM,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IREACI,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IREACP,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IREACPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBULKA,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBULKM,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBULKI,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBULKP,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBULKPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD1A,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD1M,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD1I,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD1P,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD1PH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD2A,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD2M,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD2I,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD2P,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD2PH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD3A,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD3M,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD3I,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD3P,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD3PH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD4A,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD4M,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD4I,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD4P,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCD4PH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCDEA,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCDEM,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCDEI,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCDEP,NPLS*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISCDEPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (IESTMA,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTMM,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTMI,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGKPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IESTMPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (IBGKA,NATM*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGKM,NMOL*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGKI,NION*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGKPH,NPHOT*NREAC,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NREACI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (MAXSPC,5,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

cdr something for the internal CRM options, of blocks 4,12 here: H_Colrad.
      CALL MPI_BCAST (NHCOL_STORE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (M_HCOL,NREAC,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

cnh macroscopic length scale
      CALL MPI_BCAST (L_MACRO,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_COMXS


!++++++ This is a new version of SUBROUTINE EIRENE_BROAD_FIT_FORM,
!++++++ where dynamical data structures are processed with care
!++++++ IYS 27.02.2015

      SUBROUTINE EIRENE_BROAD_FIT_FORM (RP, ME)
cdr  input:  pointer RP = REACDAT(IR)%TYP%.... with TYP = POT, CRS, RTC,...
      USE EIRMOD_MPI
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      IMPLICIT NONE
      TYPE(FIT_FORMS), POINTER :: RP
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER, ND1, ND2
      EXTERNAL :: EIRENE_CHECK_EXIT, EIRENE_EXIT_OWN
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

C.....................................................................
cdr broadcast A&M data, general for a process, independent of data structure RP%IFIT
      CALL MPI_BCAST (RP%IFIT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)


cdr broadcast data for asymptotic extrapolations
      CALL MPI_BCAST (RP%JFEX1MN,1,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%JFEX1MX,1,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%JFEX2MN,1,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%JFEX2MX,1,MPI_INTEGER,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%RC1MIN,1,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%RC1MAX,1,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%RC2MIN,1,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%RC2MAX,1,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%FP1L,3,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%FP1R,3,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%FP2B,3,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RP%FP2T,3,MPI_REAL8,
     .                0,MPI_COMM_WORLD,ier)

C.....................................................................
      IF (RP%IFIT < 0) THEN
C DATA FOR PHOTONIC LINE SHAPE AND LINE TRANSPORT
cdr 2020:  this code is mostly obsolete. This structure stems from high pressure
cdr        gas discharges (lighting) applications 2000-2002.
cdr        Mostly irrelevant for radiation tranfer in fusion applications.
        IF (ME .NE. 0) THEN
          IF (.NOT.ASSOCIATED(RP%LINE)) ALLOCATE (RP%LINE)  ! IYS
        ENDIF

        CALL MPI_BCAST (RP%LINE%E0,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%E1,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%AIK,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%G1,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%G2,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%C2,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%C3,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%C4,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%C6,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%B12,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%B21,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%C6A,12,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IGND,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IRCART,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IPROFILETYPE,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IFREMD,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%NRJPRT,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IMESS,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%IPLSC6,12,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%REACNAME,50,MPI_CHARACTER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%LINE%KENN,24,MPI_CHARACTER,
     .                  0,MPI_COMM_WORLD,ier)

C.....................................................................
      ELSE IF (1<= RP%IFIT.AND.RP%IFIT <= 2) THEN
C  POLYNOMIAL FIT, either 1D  (RP%IFIT=1),
C                  or     2D  (RP%IFIT=2)
        IF (ME == 0) THEN
          ND1 = UBOUND(RP%POLY%DBLPOL,1)
          ND2 = UBOUND(RP%POLY%DBLPOL,2)
        END IF

        CALL MPI_BCAST (ND1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (ND2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

        IF (ME .NE. 0) THEN
          IF (associated(RP%POLY)) THEN ! IYS 27.02.2015
            IF (associated(RP%POLY%DBLPOL)) THEN
              IF (ND1.ne.UBOUND(RP%POLY%DBLPOL,1) .or.
     &            ND2.ne.UBOUND(RP%POLY%DBLPOL,2)) THEN
                DEALLOCATE(RP%POLY%DBLPOL)
                NULLIFY(RP%POLY%DBLPOL)
                ALLOCATE (RP%POLY%DBLPOL(ND1,ND2))
              ENDIF
            ELSE
              ALLOCATE (RP%POLY%DBLPOL(ND1,ND2))
            ENDIF
          ELSE
            ALLOCATE (RP%POLY)
            ALLOCATE (RP%POLY%DBLPOL(ND1,ND2))
          END IF
        END IF

        CALL MPI_BCAST (RP%POLY%DBLPOL,ND1*ND2,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)

C.....................................................................
      ELSE IF (RP%IFIT == 3) THEN
! 2D TABLES, E.G. ADAS DATA
        IF (ME .NE. 0) THEN
!pb          IF(.NOT. ASSOCIATED(RP%ADAS)) ALLOCATE (RP%ADAS)  ! IYS
!pb  NULLIFY NOT YET ALLOCATED POINTER ARRAYS
          IF(.NOT. ASSOCIATED(RP%ADAS)) THEN
             ALLOCATE (RP%ADAS) ! IYS
             NULLIFY(RP%ADAS%DENS)
             NULLIFY(RP%ADAS%TEMP)
             NULLIFY(RP%ADAS%DDE)
             NULLIFY(RP%ADAS%DTE)
             NULLIFY(RP%ADAS%TAB2D)
          END IF
        END IF

        CALL MPI_BCAST (RP%ADAS%NDENS,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%ADAS%NTEMP,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)

        IF (ME .NE. 0) THEN
          IF (ASSOCIATED(RP%ADAS%DENS)) THEN
            IF (RP%ADAS%NDENS.ne.UBOUND(RP%ADAS%DENS,1)) THEN
               DEALLOCATE (RP%ADAS%DENS)
               NULLIFY (RP%ADAS%DENS)
               ALLOCATE (RP%ADAS%DENS(RP%ADAS%NDENS))
            ENDIF
          ELSE
            ALLOCATE (RP%ADAS%DENS(RP%ADAS%NDENS))
          ENDIF
          IF (ASSOCIATED(RP%ADAS%TEMP)) THEN
            IF (RP%ADAS%NTEMP.ne.UBOUND(RP%ADAS%TEMP,1)) THEN
               DEALLOCATE (RP%ADAS%TEMP)
               NULLIFY (RP%ADAS%TEMP)
               ALLOCATE (RP%ADAS%TEMP(RP%ADAS%NTEMP))
            ENDIF
          ELSE
            ALLOCATE (RP%ADAS%TEMP(RP%ADAS%NTEMP))
          ENDIF
          IF (ASSOCIATED(RP%ADAS%DDE)) THEN
            IF (RP%ADAS%NDENS.ne.UBOUND(RP%ADAS%DDE,1)) THEN
               DEALLOCATE (RP%ADAS%DDE)
               NULLIFY (RP%ADAS%DDE)
               ALLOCATE (RP%ADAS%DDE(RP%ADAS%NDENS))
            ENDIF
          ELSE
            ALLOCATE (RP%ADAS%DDE(RP%ADAS%NDENS))
          ENDIF
          IF (ASSOCIATED(RP%ADAS%DTE)) THEN
            IF (RP%ADAS%NTEMP.ne.UBOUND(RP%ADAS%DTE,1)) THEN
               DEALLOCATE (RP%ADAS%DTE)
               NULLIFY (RP%ADAS%DTE)
               ALLOCATE (RP%ADAS%DTE(RP%ADAS%NTEMP))
            ENDIF
          ELSE
            ALLOCATE (RP%ADAS%DTE(RP%ADAS%NTEMP))
          ENDIF
          IF (ASSOCIATED(RP%ADAS%TAB2D)) THEN
            IF ((RP%ADAS%NTEMP.ne.UBOUND(RP%ADAS%TAB2D,1)) .or.
     &          (RP%ADAS%NDENS.ne.UBOUND(RP%ADAS%TAB2D,2))) THEN
              DEALLOCATE (RP%ADAS%TAB2D)
              NULLIFY (RP%ADAS%TAB2D)
              ALLOCATE (RP%ADAS%TAB2D(RP%ADAS%NTEMP,RP%ADAS%NDENS))
            END IF
          ELSE
            ALLOCATE (RP%ADAS%TAB2D(RP%ADAS%NTEMP,RP%ADAS%NDENS))
          ENDIF
        END IF

        CALL MPI_BCAST (RP%ADAS%DENS,RP%ADAS%NDENS,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%ADAS%TEMP,RP%ADAS%NTEMP,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%ADAS%DDE,RP%ADAS%NDENS,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%ADAS%DTE,RP%ADAS%NTEMP,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%ADAS%TAB2D,RP%ADAS%NTEMP*RP%ADAS%NDENS,
     .                  MPI_REAL8,0,MPI_COMM_WORLD,ier)

C.....................................................................
C     ELSE IF (RP%IFIT == 4) THEN
cdr 1D TABLES. E.G. FORMERLY: HYDKIN DATA: out. moved to "snippets_hydkin"

C.....................................................................
      ELSE IF (RP%IFIT == 5 ) THEN
cdr  atomic data from internal CR model codes: H_colrad, He_colrad, etc.
        IF (ME .NE. 0) THEN
          IF (.NOT.ASSOCIATED(RP%CRM)) THEN
            ALLOCATE (RP%CRM)
          END IF
        END IF

        CALL MPI_BCAST (RP%CRM%IFLAV,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%CRM%IVARST,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
C STUFF TO TRANSFER POPULATION ESCAPE FACTORS INTO INTERNAL CRM ROUTINES
        CALL MPI_BCAST (RP%CRM%IROW_ESC,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%CRM%ICOL_ESC,1,MPI_INTEGER,
     .                  0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (RP%CRM%POP_ESC,1,MPI_REAL8,
     .                  0,MPI_COMM_WORLD,ier)

C.....................................................................
      ELSE
cdr     INVALID RP%IFIT

        WRITE (iunout,*) 'ERROR IN BROADCAST:'
        WRITE (iunout,*) 'IFIT DATA FORMAT NOT AVAILABLE'
        WRITE (iunout,*) 'IFIT ',RP%IFIT
        CALL EIRENE_EXIT_OWN(1)
      END IF

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROAD_FIT_FORM

      END MODULE EIRMOD_COMXS
