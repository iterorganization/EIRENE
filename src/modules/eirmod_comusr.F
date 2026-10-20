cdr  may 2017:  preparing for storage reduction by elimination of unnecessary input tallies:
cdr             commenting,
cdr             lusr, musr, nusr, nplpr1, nplpr2, nsfprm made local,
cdr             rather than public
cpb  Dec. 2017: remove type SPECT_ARRAY, not needed in Fortran 2003
cpb  jan 2018:  remove unused arrays TEDTEDX, TEDTEDY, TEDTEDZ
cdr             added: nstpi (formerly: coutou), naddcor
cdr  oct 2018:  BV_VEC moved into LBSMO condition
cdr             POT  moved into LESMO condition
cdr             tbd: BXPERP, BYPERP, PSI: move into LBSMO condition
cdr             missing: dealloc_corners  ??
cdr             remove redundant tally LGDFT (also from LUSR)
cdr  jan 19  :  nains, naint moved here, formerly: ccoupl
cdr             input tally no. 25 added: PSI, poloidal magn. flux.
cdr  jan 20  :  comments, cleanup, and logical connection
cdr             established (started) for physically connected tallies.
cdr             E.g. for vector components of the same field,
cdr             wrt. smoothing options and input tally disabling.
cdr             This is WIP, unfinished.
cdr             See variables LDSMO, LVSMO, LBSMO, LESMO.
cdr             It is still not ensured that LIVTALI and LSMOPRO flags
cdr             lead to physically or mathematically meaningful
cdr             setting always.
cdr  may 20  :  remove unused variables: natmi_in, nmoli_in,...
cdr             add npls_fix:  number of background species kept fixed
cdr             in eirene. Only npls_fix+1,...,npls are stored stream fort.13shrt,
cdr             for iterations (BGK, Photons,...)
cdr jan.22:     generation limit parameters made target, for unified
cdr             routine COLLIDE.f

      MODULE EIRMOD_COMUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COMUSR, EIRENE_DEALLOC_COMUSR,
     P          EIRENE_INIT_COMUSR, EIRENE_BROADCAST_COMUSR,
     P          EIRENE_ALLOC_CORNERS,
     P          EIRENE_ASSOCIATE_COMUSR,
     P          EIRENE_COMUSR_REINIT,
     P          TEMPERATURE, DENSITY, VELOCITY, VOLUMEP
#ifndef NO_JSON
     P         ,AMDOutParsTYPE
#endif

      INTEGER, SAVE :: IFIRST=0
      INTEGER, SAVE ::
     P NPLPR1, NSFPRM, NPLPR2  ! internal, not public.
                               ! former storage tests
                               ! in setprm are abandoned
      INTEGER, PUBLIC, SAVE ::
     P NPLPRM  ! nplprm, is also used in setprm, for a storage test.
c
      INTEGER, SAVE ::
     P MUSR,   LUSR            ! also only local in this module,
                               ! apparently

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R         PLSTLS(:,:)

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R          CEMETERYP(:,:)

      REAL(DP), POINTER, PUBLIC, SAVE ::
C  NPLPRM, REAL.
C  THE FIRST NPLPR1 DATA ARE PRIMARY INPUT PROFILES, SET IN SUBROUTINE PLASMA
     R        TEIN(:),        TIIN(:,:),      DEIN(:),     DIIN(:,:),
     R        VXIN(:,:),      VYIN(:,:),      VZIN(:,:),
     R        BXIN(:),        BYIN(:),        BZIN(:),     BFIN(:),
     R        ADIN(:,:),
cdr Up to this point: also storage is on PLASMA_BCKGRND
cdr June 23:  BXPERP and BYPERP (see below) added to PLASMA_BCKGRND
cdr           as an alternative to "deriving" these from the B field vector.
cdr    (for background data transfer, iteration), except the derived tally DEIN.
cdr So temporarily this storage is allocated twice.
cdr Next: Derived and assistant input tallies:
     R        VOL(:),         WGHT(:,:),
     R        EXIN(:),        EYIN(:),        EZIN(:),     EFIN(:),
     R        POT(:),
c  derived from primary input profils, in subr. PLASMA_DERIV
c  (strictly: tally no. -3: DEIN, is also a derived tally) :
cdr june 2023: BXPERP, BYPERP
cdr            rather than the coding in plasma_deriv, these
cdr            vector components (grad(PSI))
cdr            can now also be directly transferred from external code indpro(6)=9.
cdr            or be read from file (indpro(5)=6,7)
     R        BXPERP(:),      BYPERP(:),
cdr   only set in plasma_deriv.
     R        BV_VEC(:,:),    PMOM_VEC(:,:),    EDRIFT(:,:),
c  B field fluxfunction PSI, corresponds to POT tally for electric field
     R        PSI(:),         ZIIN(:,:),      FREE27(:),
     R        FREE28(:),      FREE29(:),      FREE30(:),

c  Tallies 31 --120: optional: gradients of all scalar input tallies.
c  For the vectorial input tallies (V_IN, B_IN, E_IN) these gradients
c  then provide the full dyad (all nine components).

     R        DTEDX(:),       DTEDY(:),       DTEDZ(:),
     R        DTIDX(:,:),     DTIDY(:,:),     DTIDZ(:,:),
     R        DDEDX(:),       DDEDY(:),       DDEDZ(:),
     R        DDIDX(:,:),     DDIDY(:,:),     DDIDZ(:,:),
     R        DVXDX(:,:),     DVXDY(:,:),     DVXDZ(:,:),
     R        DVYDX(:,:),     DVYDY(:,:),     DVYDZ(:,:),
     R        DVZDX(:,:),     DVZDY(:,:),     DVZDZ(:,:),
     R        DBXDX(:),       DBXDY(:),       DBXDZ(:),
     R        DBYDX(:),       DBYDY(:),       DBYDZ(:),
     R        DBZDX(:),       DBZDY(:),       DBZDZ(:),
     R        DBFDX(:),       DBFDY(:),       DBFDZ(:),
     R        DADINDX(:,:),   DADINDY(:,:),   DADINDZ(:,:),
     R        DVOLDX(:),      DVOLDY(:),      DVOLDZ(:),
     R        DWGHTDX(:,:),   DWGHTDY(:,:),   DWGHTDZ(:,:),
     R        DEXDX(:),       DEXDY(:),       DEXDZ(:),
     R        DEYDX(:),       DEYDY(:),       DEYDZ(:),
     R        DEZDX(:),       DEZDY(:),       DEZDZ(:),
     R        DEFDX(:),       DEFDY(:),       DEFDZ(:),
     R        DPOTDX(:),      DPOTDY(:),      DPOTDZ(:),
c
     R        DBXPERPDX(:),     DBXPERPDY(:),     DBXPERPDZ(:),
     R        DBYPERPDX(:),     DBYPERPDY(:),     DBYPERPDZ(:),
     R        DBV_VECDX(:,:),   DBV_VECDY(:,:),   DBV_VECDZ(:,:),
     R        DPMOM_VECDX(:,:), DPMOM_VECDY(:,:), DPMOM_VECDZ(:,:),
     R        DEDRIFTDX(:,:), DEDRIFTDY(:,:), DEDRIFTDZ(:,:),
     R        DPSIDX(:),      DPSIDY(:),      DPSIDZ(:),
     R        DZIDX(:,:),     DZIDY(:,:),     DZIDZ(:,:),
     R        DFREE27DX(:),   DFREE27DY(:),   DFREE27DZ(:),
     R        DFREE28DX(:),   DFREE28DY(:),   DFREE28DZ(:),
     R        DFREE29DX(:),   DFREE29DY(:),   DFREE29DZ(:),
     R        DFREE30DX(:),   DFREE30DY(:),   DFREE30DZ(:)

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LIVTALI(:)

c  pointer to LIVTALI: active or inactive input tallies
      LOGICAL, PUBLIC, POINTER, SAVE ::
c  background, drifting Maxwellian parameters
     L         LTEIN,      LTIIN,      LDEIN,     LDIIN,
     L         LVXIN,      LVYIN,      LVZIN,
c  magn. field
     L         LBXIN,      LBYIN,      LBZIN,     LBFIN,
c  additional stuff....
     L         LADIN,      LVOL,       LWGHT,
c  electr. field
     L         LEXIN,      LEYIN,      LEZIN,     LEFIN,
     L         LPOT,
c  derived tallies
     L         LBXPERP,    LBYPERP,  !dr will become redundant, then:
                                     !optional gradient tally grad(PSI)
     L         LBV_VEC,    LPMOM_VEC,  LEDRIFT,
c  added later, belongs to magn. field tallies
     L         LPSI,                 !dr not ready,
                                     !   just placeholder by now.
c  average charge
     L         LZIIN,
c  free tallies, unused
     L         LFREE27,
     L         LFREE28,    LFREE29,    LFREE30,

c  gradient tallies
c  background, drifting Maxwellian parameters
     L         LDTEDX,     LDTEDY,     LDTEDZ,
     L         LDTIDX,     LDTIDY,     LDTIDZ,
     L         LDDEDX,     LDDEDY,     LDDEDZ,
     L         LDDIDX,     LDDIDY,     LDDIDZ,
c  flow fields
     L         LDVXDX,     LDVXDY,     LDVXDZ,
     L         LDVYDX,     LDVYDY,     LDVYDZ,
     L         LDVZDX,     LDVZDY,     LDVZDZ,
c  magn. field
     L         LDBXDX,     LDBXDY,     LDBXDZ,
     L         LDBYDX,     LDBYDY,     LDBYDZ,
     L         LDBZDX,     LDBZDY,     LDBZDZ,
     L         LDBFDX,     LDBFDY,     LDBFDZ,
C  gradient of PSI-function tally (for poloidal magn. field)
     L         LDPSIDX,    LDPSIDY,    LDPSIDZ,
c
     L         LDADINDX,   LDADINDY,   LDADINDZ,
     L         LDVOLDX,    LDVOLDY,    LDVOLDZ,
     L         LDWGHTDX,   LDWGHTDY,   LDWGHTDZ,
c  electr. field
     L         LDEXDX,     LDEXDY,     LDEXDZ,
     L         LDEYDX,     LDEYDY,     LDEYDZ,
     L         LDEZDX,     LDEZDY,     LDEZDZ,
     L         LDEFDX,     LDEFDY,     LDEFDZ,
C  gradient of electric potential
     L         LDPOTDX,    LDPOTDY,    LDPOTDZ,
C  gradients of derived tallies
     L         LDBXPERPDX, LDBXPERPDY, LDBXPERPDZ,
     L         LDBYPERPDX, LDBYPERPDY, LDBYPERPDZ,
     L         LDBV_VECDX,   LDBV_VECDY,   LDBV_VECDZ,
     L         LDPMOM_VECDX, LDPMOM_VECDY, LDPMOM_VECDZ,
     L         LDEDRIFTDX, LDEDRIFTDY, LDEDRIFTDZ,
C  gradient of average ion charge
     L         LDZIDX,     LDZIDY,     LDZIDZ,
C  gradient of free tallies
     L         LDFREE27DX, LDFREE27DY, LDFREE27DZ,
     L         LDFREE28DX, LDFREE28DY, LDFREE28DZ,
     L         LDFREE29DX, LDFREE29DY, LDFREE29DZ,
     L         LDFREE30DX, LDFREE30DY, LDFREE30DZ

      LOGICAL, PUBLIC, SAVE ::
cdr Jan. 2020
cdr try to enforce physical consisteny in new option: living "input" tallies,
cdr in case of physically related tallies, e.g. vector components,
cdr flow-field, B-field, E-field  fields.
cdr In case of vector fields, all components (and their modulus) must be
cdr in the same "smoothing category"
     L         LDIN, LVIN,  LBIN,  LEIN

      REAL(DP), ALLOCATABLE, PUBLIC, SAVE ::

C  NPLPR2, REAL.
C  THIS SECOND SET OF DATA ARE INPUT PROFILES, SET IN SUBROUTINE PLASMA_DERIV
C  FROM THE NPLPR1 BLOCK ABOVE, FOR SPEEDING UP A&M EVALUATIONS.
C  SIMILAR (TO BE MOVED HERE (?)) ZTI, ZT1,....

     R        TEINL(:),  TIINL(:,:),  DEINL(:),  DIINL(:,:),

C  NSFPRM: SURFACE-AVERAGED INPUT TALLIES (BY ABUSE OF LANGUAGE).
     R        FLXOUT(:), SAREA(:),
C
     R        DIOD(:),   DATD(:),     DMLD(:),   DPLD(:),    DPHD(:),
     R        DION(:),   DATM(:),     DMOL(:),   DPLS(:),    DPHOT(:)

!  DECLARATION AS TARGET ARRAYS FOR POINTERS USED BY UNIFIED SUBROUTINES
      REAL(DP), TARGET, ALLOCATABLE, PUBLIC, SAVE ::
     R        RMASSI(:), RMASSA(:),   RMASSM(:), RMASSPH(:), RMASSP(:)

!     POINTER FOR UNIFIED SUBROUTINES
      REAL(DP), POINTER, PUBLIC, SAVE :: RMASSX

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(RMASSX)
#endif

C     PLASMA PROFILES ON CELL VERTICES
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        CORNER_PROFILES(:,:)

c  storage for additionally setting input tallies at cell vertices,
c  for FEM-routines (interpolations, gradients,...)
      REAL(DP), POINTER, PUBLIC, SAVE ::
     .        TEINCORNER(:),   TIINCORNER(:,:), DEINCORNER(:),
     .        DIINCORNER(:,:),
     .        VXINCORNER(:,:), VYINCORNER(:,:), VZINCORNER(:,:),
     .        BXINCORNER(:),   BYINCORNER(:),   BZINCORNER(:),
     .        BFINCORNER(:),
     .        EXCORNER(:),     EYCORNER(:),     EZCORNER(:),
     .        EFCORNER(:),     POTCORNER(:),
     .        ADCORNER(:,:),   VOLCORNER(:),    WGHTCORNER(:,:),
     .        BXPERPCORNER(:), BYPERPCORNER(:),
     .        BV_VECCORNER(:,:),
     .        PMOM_VECCORNER(:,:),
     .        EDRIFTCORNER(:,:),
     .        PSICORNER(:),    ZIINCORNER(:,:), FREE27CORNER(:),
     .        FREE28CORNER(:), FREE29CORNER(:), FREE30CORNER(:)

      REAL(DP), PUBLIC, SAVE :: TVAC, DVAC, VVAC, ALLOC, ZVAC

      CHARACTER(8), ALLOCATABLE, PUBLIC, SAVE :: TEXTS(:)

C  MUSR, INTEGER
      INTEGER, PUBLIC, SAVE ::
     I         NSPH  , NPHOTI, NPHOTIM,
cdr NATMA is something specific to SOLPS. Should not be here.
cdr It is part of the old "Kukushkin rescaling magic".
cdr All other parameters for these options are in module extrab25
cdr So move that there too, as part of "eirpbls globals".
     I         NSPA  , NATMI,  NATMIM,  NATMA,
     I         NSPAM , NMOLI,  NMOLIM,
     I         NSPAMI, NIONI,  NIONIM,
     I         NSPTOT, NPLSI,  NPLSIM,  NPLS_FIX,
     I         NSNVI,  NCPVI,  NADVI,   NBGVI,
     I         NALVI,  NCLVI,  NADSI,   NALSI, NAINI, NBITS
      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         NMASSA(:), NCHARA(:), NFOLA(:),
     I         NMASSM(:), NCHARM(:), NFOLM(:),
     I         NMASSI(:), NCHARI(:), NCHRGI(:), NFOLI(:),
     I         NMASSP(:), NCHARP(:), NCHRGP(:),
     I         NFOLPH(:),
     I         NPRT(:),   ISPEZ(:,:,:,:,:,:),   ISPEZI(:,:),
     I         MPLSTI(:), MPLSV(:)

      INTEGER, ALLOCATABLE, TARGET, PUBLIC, SAVE ::
     I         NGENA(:), NGENM(:), NGENI(:), NGENPH(:)

!     POINTER FOR UNIFIED SUBROUTINES
      INTEGER, POINTER, PUBLIC, SAVE :: NGENX
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(NGENX)
#endif

      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         ISPZ_BACK(:,:)

C  LUSR, LOGICAL
      LOGICAL, ALLOCATABLE, PUBLIC, SAVE ::
     L         LGVAC(:,:)

C  FLAGS FOR "SMOOTHED INPUT TALLIES"
C  (for interpolation from cell vertices into cell)
C  (REQUIRES AVAILABILITY OF ...CORNER(:) TALLIES)
      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LSMOPRO(:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L         LTESMO,     LTISMO,     LDESMO,    LDISMO,
     L         LVXSMO,     LVYSMO,     LVZSMO,
     L         LBXSMO,     LBYSMO,     LBZSMO,    LBFSMO,
     L         LADSMO,     LVOLSMO,    LWGHTSMO,
     L         LEXSMO,     LEYSMO,     LEZSMO,    LEFSMO,
     L         LPOTSMO,
C
     L         LBXPSMO,    LBYPSMO,
     L         LBV_VECSMO, LPMOM_VECSMO, LEDRIFTSMO,
     L         LPSISMO,    LZISMO,     LFREE27SMO,
     L         LFREE28SMO, LFREE29SMO, LFREE30SMO

cdr What about smoothed gradient input tallies?

      LOGICAL, PUBLIC, SAVE ::
cdr Jan. 2020
cdr try to enforce physical consistency in new option: smoothed input tallies.
cdr smoothed density, flow-field, B-field, E-field  fields.
cdr In case of vector fieLds, all components (and their modulus) must be
cdr in the same "smoothing category"
     L         LDSMO, LVSMO,  LBSMO,  LESMO

      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         IADVE(:),  IADVS(:), IADVT(:),  IADRC(:),
     I         ICLVE(:),  ICLVS(:), ICLVT(:),  ICLRC(:),
     I         ISNVE(:),  ISNVS(:), ISNVT(:),  ISNRC(:),
     I         ICPVE(:),  ICPVS(:), ICPVT(:),  ICPRC(:),
     I         IBGVE(:),  IBGVS(:), IBGVT(:),  IBGRC(:),
     I         IADSE(:),  IADSS(:), IADST(:),  IADSC(:),
     I         NFRSTP(:), NADDP(:), NFSTPI(:), NADDCOR(:),
     I         NSPAN(:),  NSPEN(:),
     I         NSPANW(:), NSPENW(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I         NAINS(:), NAINT(:)

      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         INTLOPTS(:)

      INTEGER, ALLOCATABLE, PUBLIC, SAVE ::
     I         LKINDP(:), LKINDM(:), LKINDI(:)

      INTEGER, PUBLIC, SAVE ::
     I         NPRLL, NMODE,  NTCPU,
     I         NFILE, NFILEN, NFILEM, NFILEL, NFILEK, NFILEJ,
     I         NITER, IITER,  NTIME,  ITIMV,
     I         NEVNTMAX, NVNMXTR, KCFF

      LOGICAL, PUBLIC, SAVE :: LVNMXTR
      REAL(DP), PUBLIC, SAVE :: XVNMX(3), XVNMX2(3)
      INTEGER, PUBLIC, ALLOCATABLE, SAVE :: KVNMXTR(:)

      TYPE(EIRENE_SPECTRUM), PUBLIC, ALLOCATABLE, SAVE :: BACK_SPEC(:)
      LOGICAL, PUBLIC, ALLOCATABLE, SAVE :: LSPCCLL(:)

      !VK SPECIES DENSITY LIMIT
      REAL(DP), ALLOCATABLE, PUBLIC, SAVE :: DENSLIM(:)

      REAL(DP), ALLOCATABLE, PUBLIC, SAVE :: EION(:)

      LOGICAL, PUBLIC, SAVE :: COMUSR_FIRST_PASS(5)

! TYPE DEFINITIONS MOVED HERE FOR WRITING OF JSON FILE
      TYPE TEMPERATURE
        DOUBLE PRECISION           :: TE, TI
        INTEGER                    :: II, IDION
        TYPE(TEMPERATURE), POINTER :: NEXT
      END TYPE TEMPERATURE
C
      TYPE DENSITY
        DOUBLE PRECISION       :: DI
        INTEGER                :: II, IDION
        TYPE(DENSITY), POINTER :: NEXT
      END TYPE DENSITY
C
      TYPE VELOCITY
        DOUBLE PRECISION        :: VX, VY, VZ
        INTEGER                 :: IZ, II, IDION
        TYPE(VELOCITY), POINTER :: NEXT
      END TYPE VELOCITY
C
      TYPE VOLUMEP
        DOUBLE PRECISION       :: VOL
        INTEGER                :: II
        TYPE(VOLUMEP), POINTER :: NEXT
      END TYPE VOLUMEP
C
      TYPE(TEMPERATURE), POINTER, PUBLIC, SAVE :: TEMPLIST
      TYPE(DENSITY), POINTER, PUBLIC, SAVE :: DENLIST
      TYPE(VELOCITY), POINTER, PUBLIC, SAVE :: VELLIST
      TYPE(VOLUMEP), POINTER, PUBLIC, SAVE :: VOLLIST

      INTEGER, PUBLIC, SAVE :: NZADD
      INTEGER, ALLOCATABLE, PUBLIC, SAVE :: INI_ZONE(:), INE_ZONE(:)

#ifndef NO_JSON
cDB   AMD output
      TYPE AMDOutParsTYPE
#ifdef LEGACYCOMP
        character(LEN=24), allocatable :: AMDoutput
#else
        character(LEN=:), allocatable :: AMDoutput
#endif
        REAL(DP) :: NEmin
        REAL(DP) :: NEmax
        REAL(DP) :: TEmin
        REAL(DP) :: TEmax
        INTEGER  :: NNE
        INTEGER  :: NTE
#ifdef LEGACYCOMP
        character(LEN=24), allocatable :: MeshType
#else
        character(LEN=:), allocatable :: MeshType
#endif
      ENDTYPE AMDOutParsTYPE
      TYPE (AMDOutParsTYPE), public, save :: AMDOutpars
#endif

      CONTAINS

      SUBROUTINE EIRENE_ALLOC_COMUSR (ICAL)

      INTEGER, INTENT(IN) :: ICAL


      IF (.NOT.ALLOCATED(LSMOPRO)) ALLOCATE (LSMOPRO(NTALG))

      IF (ICAL == 1) THEN

        IF (ALLOCATED(RMASSI)) RETURN

! species (test particle and background) related data
        NPLPR2= 3*(NATM+NMOL+NION+NPLS)+4+NSPZ+2*NPHOT
     .            +NPLS

        MUSR=4*NATM+4*NMOL+5*NION+3*NPLS+27+NSPZ+2*NPHOT+
     .       6*(1+NPHOTP)*(1+NATMP)*(1+NMOLP)*(1+NIONP)*(1+NPLSP)+NSPZ*6
     .       +2*NPLS+NSPZ*NPLS
     .       +4*NADV+4*NCLV+4*NSNV+4*NADS+4*NTALI+NTALG+2*NTALV+2*NTALS
     .       +(NMOL+NION+NPLS)*4

        LUSR=NRAD*(NPLS+2)+NRAD+NTALI

        ALLOCATE (RMASSA(MAX(1,NATM)))
        ALLOCATE (RMASSM(MAX(1,NMOL)))
        ALLOCATE (RMASSI(MAX(1,NION)))
        ALLOCATE (RMASSPH(MAX(1,NPHOT)))
        ALLOCATE (RMASSP(MAX(1,NPLS)))

        ALLOCATE (DIOD(MAX(1,NION)))
        ALLOCATE (DATD(MAX(1,NATM)))
        ALLOCATE (DMLD(MAX(1,NMOL)))
        ALLOCATE (DPLD(MAX(1,NPLS)))
        ALLOCATE (DPHD(MAX(1,NPHOT)))
        ALLOCATE (DION(MAX(1,NION)))
        ALLOCATE (DATM(MAX(1,NATM)))
        ALLOCATE (DMOL(MAX(1,NMOL)))
        ALLOCATE (DPLS(MAX(1,NPLS)))
        ALLOCATE (DPHOT(MAX(1,NPHOT)))

        ALLOCATE (TEXTS(NSPZ))
        ALLOCATE (DENSLIM(NPLS)) !VK

c  integer species and background tally data

cym to be evaluated - see calling order
        if (.not.allocated(NMASSA)) THEN
          ALLOCATE (NMASSA(MAX(1,NATM)))
          COMUSR_FIRST_PASS(1) = .TRUE.
        end if
        if (.not.allocated(NCHARA)) THEN
          ALLOCATE (NCHARA(MAX(1,NATM)))
          COMUSR_FIRST_PASS(2) = .TRUE.
        end if
cym to be evaluated
        ALLOCATE (NFOLA(MAX(1,NATM)))
        ALLOCATE (NGENA(MAX(1,NATM)))
cym to be evaluated - see calling order
        if (.not.allocated(NMASSM)) THEN
          ALLOCATE (NMASSM(MAX(1,NMOL)))
          COMUSR_FIRST_PASS(3) = .TRUE.
        end if
cym to be evaluated
        ALLOCATE (NCHARM(MAX(1,NMOL)))
        ALLOCATE (NFOLM(MAX(1,NMOL)))
        ALLOCATE (NGENM(MAX(1,NMOL)))
        ALLOCATE (NMASSP(MAX(1,NPLS)))
        ALLOCATE (NCHARP(MAX(1,NPLS)))
        ALLOCATE (NCHRGP(MAX(1,NPLS)))
cym to be evaluated - see calling order / find_param
        if (.not.allocated(NMASSI)) THEN
          ALLOCATE (NMASSI(MAX(1,NION)))
          ALLOCATE (NCHARI(MAX(1,NION)))
          ALLOCATE (NCHRGI(MAX(1,NION)))
          COMUSR_FIRST_PASS(4) = .TRUE.
        end if
cym to be evaluated
        ALLOCATE (NFOLI(MAX(1,NION)))
        ALLOCATE (NGENI(MAX(1,NION)))
        ALLOCATE (NFOLPH(MAX(1,NPHOT)))
        ALLOCATE (NGENPH(MAX(1,NPHOT)))
        ALLOCATE (NPRT(MAX(1,NSPZ)))
        ALLOCATE (ISPEZ(-1:4,0:NPHOTP,0:NATMP,0:NMOLP,0:NIONP,0:NPLSP))
        ALLOCATE (ISPEZI(NSPZ,-1:4))

        ALLOCATE (MPLSTI(MAX(1,NPLS)))
        ALLOCATE (MPLSV(MAX(1,NPLS)))
        ALLOCATE (ISPZ_BACK(NSPZ,NPLS))
        ALLOCATE (IADVE(MAX(1,NADV)))
        ALLOCATE (IADVS(MAX(1,NADV)))
        ALLOCATE (IADVT(MAX(1,NADV)))
        ALLOCATE (IADRC(MAX(1,NADV)))
        ALLOCATE (ICLVE(MAX(1,NCLV)))
        ALLOCATE (ICLVS(MAX(1,NCLV)))
        ALLOCATE (ICLVT(MAX(1,NCLV)))
        ALLOCATE (ICLRC(MAX(1,NCLV)))
        ALLOCATE (ISNVE(MAX(1,NSNV)))
        ALLOCATE (ISNVS(MAX(1,NSNV)))
        ALLOCATE (ISNVT(MAX(1,NSNV)))
        ALLOCATE (ISNRC(MAX(1,NSNV)))
        ALLOCATE (IADSE(MAX(1,NADS)))
        ALLOCATE (IADSS(MAX(1,NADS)))
        ALLOCATE (IADST(MAX(1,NADS)))
        ALLOCATE (IADSC(MAX(1,NADS)))
        ALLOCATE (NFRSTP(NTALI))
        ALLOCATE (NFSTPI(NTALI))
        ALLOCATE (NADDP(NTALI))
        ALLOCATE (NADDCOR(NTALG))
        ALLOCATE (NSPAN(NTALV))
        ALLOCATE (NSPEN(NTALV))
        ALLOCATE (NSPANW(NTALS))
        ALLOCATE (NSPENW(NTALS))

        ALLOCATE (INTLOPTS(NTALI))

        IF (.NOT.ALLOCATED(LKINDM)) ALLOCATE (LKINDM(NMOL))
        IF (.NOT.ALLOCATED(LKINDI)) ALLOCATE (LKINDI(NION))
        IF (.NOT.ALLOCATED(LKINDP)) THEN
          ALLOCATE (LKINDP(NPLS))
          COMUSR_FIRST_PASS(5) = .TRUE.
        END IF

c  logicals
        ALLOCATE (LGVAC(NRAD,0:NPLS+1))
        ALLOCATE (LSPCCLL(NRAD))
        ALLOCATE (LIVTALI(NTALI))

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMUSR(1) ', NPLPR2*8 + MUSR*4 + LUSR*4

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(PLSTLS)) RETURN

c
!        NPLPR1=(12+1*NPLS+NPLSTI+3*NPLSV)*NRAD  ! background data, set in plasma.f, 17 arrays
!        NPLPRM=NPLPR1+(NAIN+NSPZMC)*NRAD        ! adin, wght,...??? adin is allocated in call with ICAL == 2
cdr BV_VEC: add nplsv to nplpr2 and remove npls from nplprm. tbd: check correct dimension of bv_vec !


        NPLPR1 = NINPTL*NRAD            ! storage only for
                                        ! active input tallies
        ALLOCATE (PLSTLS(NINPTL,NRAD))

        ALLOCATE (CEMETERYP(0:0,NRAD))  ! storage for
                                        ! inactive input tallies
        ALLOCATE (TEINL(NRAD))
        ALLOCATE (TIINL(NPLSTI,NRAD))
        ALLOCATE (DEINL(NRAD))
        ALLOCATE (DIINL(NPLS,NRAD))

        ALLOCATE (NAINS(NAIN))
        ALLOCATE (NAINT(NAIN))

c  NCPV, NBGV are now set
        ALLOCATE (ICPVE(NCPV))
        ALLOCATE (ICPVS(NCPV))
        ALLOCATE (ICPVT(NCPV))
        ALLOCATE (ICPRC(NCPV))
        ALLOCATE (IBGVE(NBGV))
        ALLOCATE (IBGVS(NBGV))
        ALLOCATE (IBGVT(NBGV))
        ALLOCATE (IBGRC(NBGV))

        WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMUSR(2) ',
     .        (NPLPR1+                   ! ACTIVE INPUT TALLIES PLSTLS
     .        (3+NPLSTI+NPLS)*NRAD)*8 +  ! TEINL,TIINL,DEINL,DIINL
     .         4*(NCPV+NBGV)*4 +         ! BGK AND CPV INTEGERS
     .         2*NAIN*4                  ! NAINS, NAINT

      ELSE IF (ICAL == 3) THEN

        IF (ALLOCATED(FLXOUT)) RETURN

        NSFPRM=2*NLMPGS
        ALLOCATE (FLXOUT(NLMPGS))
        ALLOCATE (SAREA(NLMPGS))

        WRITE (IUNMEM,'(A,T25,I15)')
     .         ' COMUSR(3) ',NSFPRM*8

        ALLOCATE(KVNMXTR(NSTRA))
        KVNMXTR=0

      ELSE IF (ICAL == 4) THEN

        IF (ALLOCATED(EION)) RETURN

        ALLOCATE (EION(MAX(1,NATM)))

        WRITE (IUNMEM,'(A,T25,I15)')
     .         ' COMUSR(4) ',NATM*8

      END IF

      CALL EIRENE_INIT_COMUSR(ICAL)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COMUSR


      SUBROUTINE EIRENE_ASSOCIATE_COMUSR
cdr special treatment of Ti: intlopts.....

      IF (LTEIN) THEN
        TEIN => PLSTLS(NADDP(1)+1,:)
      ELSE
        TEIN => CEMETERYP(0,:)
      END IF
      IF (LTIIN) THEN
        IF (INTLOPTS(2) >= 0) THEN
          TIIN => PLSTLS(NADDP(2)+1:NADDP(3),:)
        ELSE
          TIIN => PLSTLS(NADDP(1)+1:NADDP(1)+1,:) ! Ti = Te,
                                    ! no own storage for Ti
        END IF
      END IF
      IF (LDEIN) THEN
        DEIN => PLSTLS(NADDP(3)+1,:)
      ELSE
        DEIN => CEMETERYP(0,:)
      END IF
      IF (LDIIN) THEN
        DIIN => PLSTLS(NADDP(4)+1:NADDP(5),:)
      ELSE
        DIIN => CEMETERYP(0:0,:)
      END IF

      IF (LVXIN) THEN
        VXIN => PLSTLS(NADDP(5)+1:NADDP(6),:)
      ELSE
        VXIN => CEMETERYP(0:0,:)
      END IF
      IF (LVYIN) THEN
        VYIN => PLSTLS(NADDP(6)+1:NADDP(7),:)
      ELSE
        VYIN => CEMETERYP(0:0,:)
      END IF
      IF (LVZIN) THEN
        VZIN => PLSTLS(NADDP(7)+1:NADDP(8),:)
      ELSE
        VZIN => CEMETERYP(0:0,:)
      END IF

      IF (LBXIN) THEN
        BXIN => PLSTLS(NADDP(8)+1,:)
      ELSE
        BXIN => CEMETERYP(0,:)
      END IF
      IF (LBYIN) THEN
        BYIN => PLSTLS(NADDP(9)+1,:)
      ELSE
        BYIN => CEMETERYP(0,:)
      END IF
      IF (LBZIN) THEN
        BZIN => PLSTLS(NADDP(10)+1,:)
      ELSE
        BZIN => CEMETERYP(0,:)
      END IF
      IF (LBFIN) THEN
        BFIN => PLSTLS(NADDP(11)+1,:)
      ELSE
        BFIN => CEMETERYP(0,:)
      END IF

      IF (LADIN) THEN
        ADIN => PLSTLS(NADDP(12)+1:NADDP(13),:)
      ELSE
        ADIN => CEMETERYP(0:0,:)
      END IF
      IF (LEDRIFT) THEN
        EDRIFT => PLSTLS(NADDP(13)+1:NADDP(14),:)
      ELSE
        EDRIFT => CEMETERYP(0:0,:)
      END IF
      IF (LVOL) THEN
        VOL => PLSTLS(NADDP(14)+1,:)
      ELSE
        VOL => CEMETERYP(0,:)
      END IF
      IF (LWGHT) THEN
        WGHT => PLSTLS(NADDP(15)+1:NADDP(16),:)
      ELSE
        WGHT => CEMETERYP(0:0,:)
      END IF
      IF (LBXPERP) THEN
        BXPERP => PLSTLS(NADDP(16)+1,:)
      ELSE
        BXPERP => CEMETERYP(0,:)
      END IF
       IF (LBYPERP) THEN
        BYPERP => PLSTLS(NADDP(17)+1,:)
      ELSE
        BYPERP => CEMETERYP(0,:)
      END IF
cdr  BZPERP == 0, for the time being

      IF (LEXIN) THEN
        EXIN => PLSTLS(NADDP(18)+1,:)
      ELSE
        EXIN => CEMETERYP(0,:)
      END IF
      IF (LEYIN) THEN
        EYIN => PLSTLS(NADDP(19)+1,:)
      ELSE
        EYIN => CEMETERYP(0,:)
      END IF
      IF (LEZIN) THEN
        EZIN => PLSTLS(NADDP(20)+1,:)
      ELSE
        EZIN => CEMETERYP(0,:)
      END IF
      IF (LEFIN) THEN
        EFIN => PLSTLS(NADDP(21)+1,:)
      ELSE
        EFIN => CEMETERYP(0,:)
      END IF
      IF (LPOT) THEN
        POT => PLSTLS(NADDP(22)+1,:)
      ELSE
        POT => CEMETERYP(0,:)
      END IF

      IF (LBV_VEC) THEN
        BV_VEC => PLSTLS(NADDP(23)+1:NADDP(24),:)
      ELSE
        BV_VEC => CEMETERYP(0:0,:)
      END IF
      IF (LPMOM_VEC) THEN
        PMOM_VEC => PLSTLS(NADDP(24)+1:NADDP(25),:)
      ELSE
        PMOM_VEC => CEMETERYP(0:0,:)
      END IF
      IF (LPSI) THEN
        PSI => PLSTLS(NADDP(25)+1,:)
      ELSE
        PSI => CEMETERYP(0,:)
      END IF
      IF (LZIIN) THEN
        ZIIN => PLSTLS(NADDP(26)+1:NADDP(27),:)
      ELSE
        ZIIN => CEMETERYP(0:0,:)
      END IF
      IF (LFREE27) THEN
        FREE27 => PLSTLS(NADDP(27)+1,:)
      ELSE
        FREE27 => CEMETERYP(0,:)
      END IF
      IF (LFREE28) THEN
        FREE28 => PLSTLS(NADDP(28)+1,:)
      ELSE
        FREE28 => CEMETERYP(0,:)
      END IF
      IF (LFREE29) THEN
        FREE29 => PLSTLS(NADDP(29)+1,:)
      ELSE
        FREE29 => CEMETERYP(0,:)
      END IF
      IF (LFREE30) THEN
        FREE30 => PLSTLS(NADDP(30)+1,:)
      ELSE
        FREE30 => CEMETERYP(0,:)
      END IF

C  TALLIES 31--130: DERIVATIVES WRT. X,Y,Z COORDINATES OF TALLIES 1--30
      IF (LDTEDX) THEN
        DTEDX => PLSTLS(NADDP(31)+1,:)
      ELSE
        DTEDX => CEMETERYP(0,:)
      END IF
      IF (LDTEDY) THEN
        DTEDY => PLSTLS(NADDP(32)+1,:)
      ELSE
        DTEDY => CEMETERYP(0,:)
      END IF
      IF (LDTEDZ) THEN
        DTEDZ => PLSTLS(NADDP(33)+1,:)
      ELSE
        DTEDZ => CEMETERYP(0,:)
      END IF
      IF (LDTIDX) THEN
        DTIDX => PLSTLS(NADDP(34)+1:NADDP(35),:)
      ELSE
        DTIDX => CEMETERYP(0:0,:)
      END IF
      IF (LDTIDY) THEN
        DTIDY => PLSTLS(NADDP(35)+1:NADDP(36),:)
      ELSE
        DTIDY => CEMETERYP(0:0,:)
      END IF
      IF (LDTIDZ) THEN
        DTIDZ => PLSTLS(NADDP(36)+1:NADDP(37),:)
      ELSE
        DTIDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDDEDX) THEN
        DDEDX => PLSTLS(NADDP(37)+1,:)
      ELSE
        DDEDX => CEMETERYP(0,:)
      END IF
      IF (LDDEDY) THEN
        DDEDY => PLSTLS(NADDP(38)+1,:)
      ELSE
        DDEDY => CEMETERYP(0,:)
      END IF
      IF (LDDEDZ) THEN
        DDEDZ => PLSTLS(NADDP(39)+1,:)
      ELSE
        DDEDZ => CEMETERYP(0,:)
      END IF
      IF (LDDIDX) THEN
        DDIDX => PLSTLS(NADDP(40)+1:NADDP(41),:)
      ELSE
        DDIDX => CEMETERYP(0:0,:)
      END IF
      IF (LDDIDY) THEN
        DDIDY => PLSTLS(NADDP(41)+1:NADDP(42),:)
      ELSE
        DDIDY => CEMETERYP(0:0,:)
      END IF
      IF (LDDIDZ) THEN
        DDIDZ => PLSTLS(NADDP(42)+1:NADDP(43),:)
      ELSE
        DDIDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDX) THEN
        DVXDX => PLSTLS(NADDP(43)+1:NADDP(44),:)
      ELSE
        DVXDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDY) THEN
        DVXDY => PLSTLS(NADDP(44)+1:NADDP(45),:)
      ELSE
        DVXDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVXDZ) THEN
        DVXDZ => PLSTLS(NADDP(45)+1:NADDP(46),:)
      ELSE
        DVXDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDX) THEN
        DVYDX => PLSTLS(NADDP(46)+1:NADDP(47),:)
      ELSE
        DVYDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDY) THEN
        DVYDY => PLSTLS(NADDP(47)+1:NADDP(48),:)
      ELSE
        DVYDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVYDZ) THEN
        DVYDZ => PLSTLS(NADDP(48)+1:NADDP(49),:)
      ELSE
        DVYDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDX) THEN
        DVZDX => PLSTLS(NADDP(49)+1:NADDP(50),:)
      ELSE
        DVZDX => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDY) THEN
        DVZDY => PLSTLS(NADDP(50)+1:NADDP(51),:)
      ELSE
        DVZDY => CEMETERYP(0:0,:)
      END IF
      IF (LDVZDZ) THEN
        DVZDZ => PLSTLS(NADDP(51)+1:NADDP(52),:)
      ELSE
        DVZDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDBXDX) THEN
        DBXDX => PLSTLS(NADDP(52)+1,:)
      ELSE
        DBXDX => CEMETERYP(0,:)
      END IF
      IF (LDBXDY) THEN
        DBXDY => PLSTLS(NADDP(53)+1,:)
      ELSE
        DBXDY => CEMETERYP(0,:)
      END IF
      IF (LDBXDZ) THEN
        DBXDZ => PLSTLS(NADDP(54)+1,:)
      ELSE
        DBXDZ => CEMETERYP(0,:)
      END IF
      IF (LDBYDX) THEN
        DBYDX => PLSTLS(NADDP(55)+1,:)
      ELSE
        DBYDX => CEMETERYP(0,:)
      END IF
      IF (LDBYDY) THEN
        DBYDY => PLSTLS(NADDP(56)+1,:)
      ELSE
        DBYDY => CEMETERYP(0,:)
      END IF
      IF (LDBYDZ) THEN
        DBYDZ => PLSTLS(NADDP(57)+1,:)
      ELSE
        DBYDZ => CEMETERYP(0,:)
      END IF
      IF (LDBZDX) THEN
        DBZDX => PLSTLS(NADDP(58)+1,:)
      ELSE
        DBZDX => CEMETERYP(0,:)
      END IF
      IF (LDBZDY) THEN
        DBZDY => PLSTLS(NADDP(59)+1,:)
      ELSE
        DBZDY => CEMETERYP(0,:)
      END IF
      IF (LDBZDZ) THEN
        DBZDZ => PLSTLS(NADDP(60)+1,:)
      ELSE
        DBZDZ => CEMETERYP(0,:)
      END IF
      IF (LDBFDX) THEN
        DBFDX => PLSTLS(NADDP(61)+1,:)
      ELSE
        DBFDX => CEMETERYP(0,:)
      END IF
      IF (LDBFDY) THEN
        DBFDY => PLSTLS(NADDP(62)+1,:)
      ELSE
        DBFDY => CEMETERYP(0,:)
      END IF
      IF (LDBFDZ) THEN
        DBFDZ => PLSTLS(NADDP(63)+1,:)
      ELSE
        DBFDZ => CEMETERYP(0,:)
      END IF
      IF (LDADINDX) THEN
        DADINDX => PLSTLS(NADDP(64)+1:NADDP(65),:)
      ELSE
        DADINDX => CEMETERYP(0:0,:)
      END IF
      IF (LDADINDY) THEN
        DADINDY => PLSTLS(NADDP(65)+1:NADDP(66),:)
      ELSE
        DADINDY => CEMETERYP(0:0,:)
      END IF
      IF (LDADINDZ) THEN
        DADINDZ => PLSTLS(NADDP(66)+1:NADDP(67),:)
      ELSE
        DADINDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDX) THEN
        DEDRIFTDX => PLSTLS(NADDP(67)+1:NADDP(68),:)
      ELSE
        DEDRIFTDX => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDY) THEN
        DEDRIFTDY => PLSTLS(NADDP(68)+1:NADDP(69),:)
      ELSE
        DEDRIFTDY => CEMETERYP(0:0,:)
      END IF
      IF (LDEDRIFTDZ) THEN
        DEDRIFTDZ => PLSTLS(NADDP(69)+1:NADDP(70),:)
      ELSE
        DEDRIFTDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDVOLDX) THEN
        DVOLDX => PLSTLS(NADDP(70)+1,:)
      ELSE
        DVOLDX => CEMETERYP(0,:)
      END IF
      IF (LDVOLDY) THEN
        DVOLDY => PLSTLS(NADDP(71)+1,:)
      ELSE
        DVOLDY => CEMETERYP(0,:)
      END IF
      IF (LDVOLDZ) THEN
        DVOLDZ => PLSTLS(NADDP(72)+1,:)
      ELSE
        DVOLDZ => CEMETERYP(0,:)
      END IF
      IF (LDWGHTDX) THEN
        DWGHTDX => PLSTLS(NADDP(73)+1:NADDP(74),:)
      ELSE
        DWGHTDX => CEMETERYP(0:0,:)
      END IF
      IF (LDWGHTDY) THEN
        DWGHTDY => PLSTLS(NADDP(74)+1:NADDP(75),:)
      ELSE
        DWGHTDY => CEMETERYP(0:0,:)
      END IF
      IF (LDWGHTDZ) THEN
        DWGHTDZ => PLSTLS(NADDP(75)+1:NADDP(76),:)
      ELSE
        DWGHTDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDBXPERPDX) THEN
        DBXPERPDX => PLSTLS(NADDP(76)+1,:)
      ELSE
        DBXPERPDX => CEMETERYP(0,:)
      END IF
      IF (LDBXPERPDY) THEN
        DBXPERPDY => PLSTLS(NADDP(77)+1,:)
      ELSE
        DBXPERPDY => CEMETERYP(0,:)
      END IF
      IF (LDBXPERPDZ) THEN
        DBXPERPDZ => PLSTLS(NADDP(78)+1,:)
      ELSE
        DBXPERPDZ => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDX) THEN
        DBYPERPDX => PLSTLS(NADDP(79)+1,:)
      ELSE
        DBYPERPDX => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDY) THEN
        DBYPERPDY => PLSTLS(NADDP(80)+1,:)
      ELSE
        DBYPERPDY => CEMETERYP(0,:)
      END IF
      IF (LDBYPERPDZ) THEN
        DBYPERPDZ => PLSTLS(NADDP(81)+1,:)
      ELSE
        DBYPERPDZ => CEMETERYP(0,:)
      END IF
      IF (LDEXDX) THEN
        DEXDX => PLSTLS(NADDP(82)+1,:)
      ELSE
        DEXDX => CEMETERYP(0,:)
      END IF
      IF (LDEXDY) THEN
        DEXDY => PLSTLS(NADDP(83)+1,:)
      ELSE
        DEXDY => CEMETERYP(0,:)
      END IF
      IF (LDEXDZ) THEN
        DEXDZ => PLSTLS(NADDP(84)+1,:)
      ELSE
        DEXDZ => CEMETERYP(0,:)
      END IF
      IF (LDEYDX) THEN
        DEYDX => PLSTLS(NADDP(85)+1,:)
      ELSE
        DEYDX => CEMETERYP(0,:)
      END IF
      IF (LDEYDY) THEN
        DEYDY => PLSTLS(NADDP(86)+1,:)
      ELSE
        DEYDY => CEMETERYP(0,:)
      END IF
      IF (LDEYDZ) THEN
        DEYDZ => PLSTLS(NADDP(87)+1,:)
      ELSE
        DEYDZ => CEMETERYP(0,:)
      END IF
      IF (LDEZDX) THEN
        DEZDX => PLSTLS(NADDP(88)+1,:)
      ELSE
        DEZDX => CEMETERYP(0,:)
      END IF
      IF (LDEZDY) THEN
        DEZDY => PLSTLS(NADDP(89)+1,:)
      ELSE
        DEZDY => CEMETERYP(0,:)
      END IF
      IF (LDEZDZ) THEN
        DEZDZ => PLSTLS(NADDP(90)+1,:)
      ELSE
        DEZDZ => CEMETERYP(0,:)
      END IF
      IF (LDEFDX) THEN
        DEFDX => PLSTLS(NADDP(91)+1,:)
      ELSE
        DEFDX => CEMETERYP(0,:)
      END IF
      IF (LDEFDY) THEN
        DEFDY => PLSTLS(NADDP(92)+1,:)
      ELSE
        DEFDY => CEMETERYP(0,:)
      END IF
      IF (LDEFDZ) THEN
        DEFDZ => PLSTLS(NADDP(93)+1,:)
      ELSE
        DEFDZ => CEMETERYP(0,:)
      END IF
      IF (LDPOTDX) THEN
        DPOTDX => PLSTLS(NADDP(94)+1,:)
      ELSE
        DPOTDX => CEMETERYP(0,:)
      END IF
      IF (LDPOTDY) THEN
        DPOTDY => PLSTLS(NADDP(95)+1,:)
      ELSE
        DPOTDY => CEMETERYP(0,:)
      END IF
      IF (LDPOTDZ) THEN
        DPOTDZ => PLSTLS(NADDP(96)+1,:)
      ELSE
        DPOTDZ => CEMETERYP(0,:)
      END IF
      IF (LDBV_VECDX) THEN
        DBV_VECDX => PLSTLS(NADDP(97)+1:NADDP(98),:)
      ELSE
        DBV_VECDX => CEMETERYP(0:0,:)
      END IF
      IF (LDBV_VECDY) THEN
        DBV_VECDY => PLSTLS(NADDP(98)+1:NADDP(99),:)
      ELSE
        DBV_VECDY => CEMETERYP(0:0,:)
      END IF
      IF (LDBV_VECDZ) THEN
        DBV_VECDZ => PLSTLS(NADDP(99)+1:NADDP(100),:)
      ELSE
        DBV_VECDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDPMOM_VECDX) THEN
        DPMOM_VECDX => PLSTLS(NADDP(100)+1:NADDP(101),:)
      ELSE
        DPMOM_VECDX => CEMETERYP(0:0,:)
      END IF
      IF (LDPMOM_VECDY) THEN
        DPMOM_VECDY => PLSTLS(NADDP(101)+1:NADDP(102),:)
      ELSE
        DPMOM_VECDY => CEMETERYP(0:0,:)
      END IF
      IF (LDPMOM_VECDZ) THEN
        DPMOM_VECDZ => PLSTLS(NADDP(102)+1:NADDP(103),:)
      ELSE
        DPMOM_VECDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDPSIDX) THEN
        DPSIDX => PLSTLS(NADDP(103)+1,:)
      ELSE
        DPSIDX => CEMETERYP(0,:)
      END IF
      IF (LDPSIDY) THEN
        DPSIDY => PLSTLS(NADDP(104)+1,:)
      ELSE
        DPSIDY => CEMETERYP(0,:)
      END IF
      IF (LDPSIDZ) THEN
        DPSIDZ => PLSTLS(NADDP(105)+1,:)
      ELSE
        DPSIDZ => CEMETERYP(0,:)
      END IF
      IF (LDZIDX) THEN
        DZIDX => PLSTLS(NADDP(106)+1:NADDP(107),:)
      ELSE
        DZIDX => CEMETERYP(0:0,:)
      END IF
      IF (LDZIDY) THEN
        DZIDY => PLSTLS(NADDP(107)+1:NADDP(108),:)
      ELSE
        DZIDY => CEMETERYP(0:0,:)
      END IF
      IF (LDZIDZ) THEN
        DZIDZ => PLSTLS(NADDP(108)+1:NADDP(109),:)
      ELSE
        DZIDZ => CEMETERYP(0:0,:)
      END IF
      IF (LDFREE27DX) THEN
        DFREE27DX => PLSTLS(NADDP(109)+1,:)
      ELSE
        DFREE27DX => CEMETERYP(0,:)
      END IF
      IF (LDFREE27DY) THEN
        DFREE27DY => PLSTLS(NADDP(110)+1,:)
      ELSE
        DFREE27DY => CEMETERYP(0,:)
      END IF
      IF (LDFREE27DZ) THEN
        DFREE27DZ => PLSTLS(NADDP(111)+1,:)
      ELSE
        DFREE27DZ => CEMETERYP(0,:)
      END IF
      IF (LDFREE28DX) THEN
        DFREE28DX => PLSTLS(NADDP(112)+1,:)
      ELSE
        DFREE28DX => CEMETERYP(0,:)
      END IF
      IF (LDFREE28DY) THEN
        DFREE28DY => PLSTLS(NADDP(113)+1,:)
      ELSE
        DFREE28DY => CEMETERYP(0,:)
      END IF
      IF (LDFREE28DZ) THEN
        DFREE28DZ => PLSTLS(NADDP(114)+1,:)
      ELSE
        DFREE28DZ => CEMETERYP(0,:)
      END IF
      IF (LDFREE29DX) THEN
        DFREE29DX => PLSTLS(NADDP(115)+1,:)
      ELSE
        DFREE29DX => CEMETERYP(0,:)
      END IF
      IF (LDFREE29DY) THEN
        DFREE29DY => PLSTLS(NADDP(116)+1,:)
      ELSE
        DFREE29DY => CEMETERYP(0,:)
      END IF
      IF (LDFREE29DZ) THEN
        DFREE29DZ => PLSTLS(NADDP(117)+1,:)
      ELSE
        DFREE29DZ => CEMETERYP(0,:)
      END IF
      IF (LDFREE30DX) THEN
        DFREE30DX => PLSTLS(NADDP(118)+1,:)
      ELSE
        DFREE30DX => CEMETERYP(0,:)
      END IF
      IF (LDFREE30DY) THEN
        DFREE30DY => PLSTLS(NADDP(119)+1,:)
      ELSE
        DFREE30DY => CEMETERYP(0,:)
      END IF
      IF (LDFREE30DZ) THEN
        DFREE30DZ => PLSTLS(NADDP(120)+1,:)
      ELSE
        DFREE30DZ => CEMETERYP(0,:)
      END IF

      RETURN
      END SUBROUTINE EIRENE_ASSOCIATE_COMUSR



      SUBROUTINE EIRENE_ALLOC_CORNERS
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      INTEGER :: NTOT, I, J, NLSTTL, NTOT2
      INTEGER :: IDSMO, IVSMO, IBSMO, IESMO

      IF (ALLOCATED(CORNER_PROFILES)) RETURN

cdr  are there any FEM interpolated background tallies in this run?

!  interpolation to vertices can only be done if input tally is available (active)
      LSMOPRO(1:NTALG) = LSMOPRO(1:NTALG) .AND. LIVTALI(1:NTALG)

!pb  check for inconsistencies in smoothing
      IDSMO = 0
      IF (LDESMO) IDSMO = IDSMO + 1
      IF (LDISMO) IDSMO = IDSMO + 1

      IVSMO = 0
cdr  strictly: a vectorial quantity.
cdr  Smoothing is probably not correctly implemented
      IF (LVXSMO) IVSMO = IVSMO + 1
      IF (LVYSMO) IVSMO = IVSMO + 1
      IF (LVZSMO) IVSMO = IVSMO + 1
      IF (LBV_VECSMO) IVSMO = IVSMO + 1
cdr   PMOM_VEC should also be here

      IBSMO = 0
cdr  strictly: a vectorial quantity.
cdr  Smoothing is probably not correctly implemented
      IF (LBXSMO) IBSMO = IBSMO + 1
      IF (LBYSMO) IBSMO = IBSMO + 1
      IF (LBZSMO) IBSMO = IBSMO + 1
      IF (LBFSMO) IBSMO = IBSMO + 1
cdr    LBXPSMO, LBYPSMO added
      IF (LBXPSMO) IBSMO = IBSMO + 1
      IF (LBYPSMO) IBSMO = IBSMO + 1
cdr   both BV_VEC and PMOM_VEC should also be here,
cdr   as well as PSI

      IESMO = 0
cdr  strictly: a vectorial quantity.
cdr  Smoothing is probably not correctly implemented
      IF (LEXSMO) IESMO = IESMO + 1
      IF (LEYSMO) IESMO = IESMO + 1
      IF (LEZSMO) IESMO = IESMO + 1
      IF (LEFSMO) IESMO = IESMO + 1
cdr   POT should also be here

      IF ((IDSMO > 0) .AND. (IDSMO < 2)) THEN
        WRITE (IUNOUT,*)
     .    'INCONSISTENCY IN SMOOTHING OF DENSITY ENCOUNTERED'
        IF (.NOT.LDESMO) THEN
          LDESMO =.TRUE.
          WRITE (IUNOUT,*)
     .         'SMOOTHING OF ELECTRON DENSITY SWITCHED ON'
        END IF
        IF (.NOT.LDISMO) THEN
          LDISMO =.TRUE.
          WRITE (IUNOUT,*)
     .         'SMOOTHING OF ION DENSITY SWITCHED ON'
        END IF
      END IF

cdr June 2023:
cdr the code below (and several places elsewhere)
cdr handles components of a vecter
cdr as if they were independent scalars.
cdr This can lead to spectacular unphysical results.
cdr A cleanup of input tally options for vectorial quantities
cdr V, B, E is needed, in particular for div-free B field.
cdr 2022: already done for BV_VEC and PMOM_VEC,

      IF ((IVSMO > 0) .AND. (IVSMO < 4)) THEN
        WRITE (IUNOUT,*)
     .    'INCONSISTENCY IN SMOOTHING OF DRIFT VELOCITIES ENCOUNTERED'
        IF (.NOT.LVXSMO) THEN
          LVXSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF DRIFT VELOCITY IN X-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LVYSMO) THEN
          LVYSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF DRIFT VELOCITY IN Y-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LVZSMO) THEN
          LVZSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF DRIFT VELOCITY IN Z-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LBV_VECSMO) THEN
          LBV_VECSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF FLOW VELOCITY WRT B-FIELD SWITCHED ON'
        END IF
      END IF

cdr legal values IBSMO: =0 or =6
      IF ((IBSMO > 0) .AND. (IBSMO < 6)) THEN
        WRITE (IUNOUT,*)
     .       'INCONSISTENCY IN SMOOTHING OF '//
     .       'MAGN. FIELD VECTORS ENCOUNTERED'
        IF (.NOT.LBXSMO) THEN
          LBXSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF MAGN. FIELD VECTOR, X-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LBYSMO) THEN
          LBYSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF MAGN. FIELD VECTOR, Y-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LBZSMO) THEN
          LBZSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF MAGN. FIELD VECTOR, Z-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LBFSMO) THEN
          LBFSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF MAGN. FIELD STRENGTH SWITCHED ON'
        END IF
        IF (.NOT.LBXPSMO) THEN
          LBXPSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF PERP MAGN FIELD DIRECTION X SWITCHED ON'
        END IF
        IF (.NOT.LBYPSMO) THEN
          LBYPSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF PERP MAGN. FIELD DIRECTION Y SWITCHED ON'
        END IF
      END IF

cdr legal values IESMO: =0 or =4
      IF ((IESMO > 0) .AND. (IESMO < 4)) THEN
        WRITE (IUNOUT,*)
     .       'INCONSISTENCY IN SMOOTHING OF '//
     .       'ELEC. FIELD VECTORS ENCOUNTERED'
        IF (.NOT.LEXSMO) THEN
          LEXSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF ELEC. FIELD VECTOR, X-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LEYSMO) THEN
          LEYSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF ELEC. FIELD VECTOR, Y-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LEZSMO) THEN
          LBZSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF ELEC. FIELD VECTOR, Z-DIRECTION SWITCHED ON'
        END IF
        IF (.NOT.LEFSMO) THEN
          LEFSMO =.TRUE.
          WRITE (IUNOUT,*)
     .      'SMOOTHING OF ELEC. FIELD STRENGTH SWITCHED ON'
cdr  electr. potential POT missing here. iesmo=0 or =5 are legal
        END IF
      END IF

cdr changed Jan 2020: to avoid physical inconsistency, enforce links
cdr                   between related input tallies (e.g. vector components)
cdr not finished. WIP.
cdr but at least now: if one component, then automatically all.
      LDSMO = LDESMO .AND. LDISMO
      LVSMO = LVXSMO .AND. LVYSMO .AND. LVZSMO .AND. LBV_VECSMO
      LBSMO = LBXSMO .AND. LBYSMO .AND. LBZSMO .AND. LBFSMO
CDR  .       .OR.LPSISMO
      LESMO = LEXSMO .AND. LEYSMO .AND. LEZSMO .AND. LEFSMO
CDR  .       .OR.LPOTSMO

c  NTOT2: total number of smoothed tallies, counting also with species index
      NTOT2 = 0
      DO I= 1, NTALG
        IF (LSMOPRO(I)) THEN
          NTOT2 = NTOT2 + NFRSTP(I)
        END IF
      END DO

c  NADDCOR: cumulated index of position of smoothed tally J within all smoothed tallies
C  NLSTLL : highest tally index J amongst all smoothed tallies
      NADDCOR(1)=0
      DO 6 J=2,NTALG
        IF (LSMOPRO(J-1)) THEN
          NADDCOR(J)=NADDCOR(J-1)+NFRSTP(J-1)
          NLSTTL=J-1
        ELSE
          NADDCOR(J)=NADDCOR(J-1)
        END IF
    6 CONTINUE

      IF (LSMOPRO(NTALG)) NLSTTL = NTALG
C
c  NTOT: total number of smoothed tallies, counting also with species index
      NTOT = 0
      IF (ANY(LSMOPRO)) THEN
        NTOT = NADDCOR(NTALG)
        IF (NFRSTP(NLSTTL) > 1) NTOT = NTOT+NFRSTP(NLSTTL)
      END IF

      IF (NTOT > 0) THEN
cdr  allocate storage for background tallies on cell vertices
cdr  ncorner is set in GRID.f (levgeo=4,5) or in SNEIGH.f (levgeo=1,2,3)
        ALLOCATE (CORNER_PROFILES(NCORNER,NTOT))
      ELSE
        ALLOCATE (CORNER_PROFILES(1,1))
      END IF

      WRITE (IUNMEM,'(A,T25,I15)')
     .        ' COMUSR(CORNERS) ',SIZE(CORNER_PROFILES)*8

      IF (LTESMO) THEN
        TEINCORNER => CORNER_PROFILES(:,NADDCOR(1)+1)
      ELSE
        NULLIFY(TEINCORNER)
      END IF

      IF (LTISMO) THEN
        TIINCORNER => CORNER_PROFILES(:,NADDCOR(2)+1 : NADDCOR(3))
      ELSE
        NULLIFY(TIINCORNER)
      END IF

      IF (LDSMO) THEN
cdr Enforce smoothing interpolations for all charged particles
cdr incl. electrons, in order to not wreck quasineutrality.
        IF (LDESMO) THEN
          DEINCORNER => CORNER_PROFILES(:,NADDCOR(3)+1)
        ELSE
          NULLIFY(DEINCORNER)
        END IF
        IF (LDISMO) THEN
          DIINCORNER => CORNER_PROFILES(:,NADDCOR(4)+1 : NADDCOR(5))
        ELSE
          NULLIFY(DIINCORNER)
        END IF
      ELSE
        NULLIFY(DEINCORNER)
        NULLIFY(DIINCORNER)
      END IF

      IF (LVSMO) THEN
        IF (LVXSMO) THEN
          VXINCORNER => CORNER_PROFILES(:,NADDCOR(5)+1 : NADDCOR(6))
        ELSE
          NULLIFY(VXINCORNER)
        END IF
        IF (LVYSMO) THEN
          VYINCORNER => CORNER_PROFILES(:,NADDCOR(6)+1 : NADDCOR(7))
        ELSE
          NULLIFY(VYINCORNER)
        END IF
        IF (LVZSMO) THEN
          VZINCORNER => CORNER_PROFILES(:,NADDCOR(7)+1 : NADDCOR(8))
        ELSE
          NULLIFY(VZINCORNER)
        END IF
        IF (LBSMO) THEN
          BV_VECCORNER => CORNER_PROFILES(:,NADDCOR(23)+1 : NADDCOR(24))
        ELSE
          NULLIFY(BV_VECCORNER)
        END IF
      ELSE
        NULLIFY(VXINCORNER)
        NULLIFY(VYINCORNER)
        NULLIFY(VZINCORNER)
        NULLIFY(BV_VECCORNER)
      END IF
cdr  same simplification needed for other vectors...
      IF (LBSMO) THEN
        IF (LBXSMO) THEN
          BXINCORNER => CORNER_PROFILES(:,NADDCOR(8)+1)
        ELSE
          NULLIFY(BXINCORNER)
        END IF
        IF (LBYSMO) THEN
          BYINCORNER => CORNER_PROFILES(:,NADDCOR(9)+1)
        ELSE
          NULLIFY(BYINCORNER)
        END IF
        IF (LBZSMO) THEN
          BZINCORNER => CORNER_PROFILES(:,NADDCOR(10)+1)
        ELSE
          NULLIFY(BZINCORNER)
        END IF
        IF (LBFSMO) THEN
          BFINCORNER => CORNER_PROFILES(:,NADDCOR(11)+1)
        ELSE
          NULLIFY(BFINCORNER)
        END IF
        IF (LBXPSMO) THEN
          BXPERPCORNER => CORNER_PROFILES(:,NADDCOR(16)+1 )
        ELSE
          NULLIFY(BXPERPCORNER)
        END IF
        IF (LBYPSMO) THEN
          BYPERPCORNER => CORNER_PROFILES(:,NADDCOR(17)+1 )
        ELSE
          NULLIFY(BYPERPCORNER)
        END IF
      ELSE
        NULLIFY(BXINCORNER)
        NULLIFY(BYINCORNER)
        NULLIFY(BZINCORNER)
        NULLIFY(BFINCORNER)
        NULLIFY(BXPERPCORNER)
        NULLIFY(BYPERPCORNER)
      END IF

      IF (LADSMO) THEN
        ADCORNER => CORNER_PROFILES(:,NADDCOR(12)+1 : NADDCOR(13))
      ELSE
        NULLIFY(ADCORNER)
      END IF

      IF (LEDRIFTSMO) THEN
        EDRIFTCORNER => CORNER_PROFILES(:,NADDCOR(13)+1 : NADDCOR(14))
      ELSE
        NULLIFY(EDRIFTCORNER)
      END IF

      IF (LVOLSMO) THEN
        VOLCORNER => CORNER_PROFILES(:,NADDCOR(14)+1 )
      ELSE
        NULLIFY(VOLCORNER)
      END IF

      IF (LWGHTSMO) THEN
        WGHTCORNER => CORNER_PROFILES(:,NADDCOR(15)+1 : NADDCOR(16))
      ELSE
        NULLIFY(WGHTCORNER)
      END IF

      IF (LESMO) THEN
        IF (LEXSMO) THEN
          EXCORNER => CORNER_PROFILES(:,NADDCOR(18)+1)
        ELSE
          NULLIFY(EXCORNER)
        END IF
        IF (LEYSMO) THEN
          EYCORNER => CORNER_PROFILES(:,NADDCOR(19)+1)
        ELSE
          NULLIFY(EYCORNER)
        END IF
        IF (LEZSMO) THEN
          EZCORNER => CORNER_PROFILES(:,NADDCOR(20)+1)
        ELSE
          NULLIFY(EZCORNER)
        END IF
        IF (LEFSMO) THEN
          EFCORNER => CORNER_PROFILES(:,NADDCOR(21)+1)
        ELSE
          NULLIFY(EFCORNER)
        END IF
        IF (LPOTSMO) THEN
          POTCORNER => CORNER_PROFILES(:,NADDCOR(22)+1)
        ELSE
          NULLIFY(POTCORNER)
        END IF

      ELSE  !dr no smoothed E field
        NULLIFY(EXCORNER)
        NULLIFY(EYCORNER)
        NULLIFY(EZCORNER)
        NULLIFY(EFCORNER)
        NULLIFY(POTCORNER)
      END IF
cdr  hidden link warning:
cdr  this should go into LVSMO
      IF (LPMOM_VECSMO) THEN
        PMOM_VECCORNER => CORNER_PROFILES(:,NADDCOR(24)+1 : NADDCOR(25))
      ELSE
        NULLIFY(PMOM_VECCORNER)
      END IF
cdr  hidden link warning:
cdr  this should go into LBSMO
      IF (LPSISMO) THEN
        PSICORNER => CORNER_PROFILES(:,NADDCOR(25)+1)
      ELSE
        NULLIFY(PSICORNER)
      END IF

      IF (LZISMO) THEN
        ZIINCORNER => CORNER_PROFILES(:,NADDCOR(26)+1 : NADDCOR(27))
      ELSE
        NULLIFY(ZIINCORNER)
      END IF

      IF (LFREE27SMO) THEN
        FREE27CORNER => CORNER_PROFILES(:,NADDCOR(27)+1)
      ELSE
        NULLIFY(FREE27CORNER)
      END IF

      IF (LFREE28SMO) THEN
        FREE28CORNER => CORNER_PROFILES(:,NADDCOR(28)+1)
      ELSE
        NULLIFY(FREE28CORNER)
      END IF

      IF (LFREE29SMO) THEN
        FREE29CORNER => CORNER_PROFILES(:,NADDCOR(29)+1)
      ELSE
        NULLIFY(FREE29CORNER)
      END IF

      IF (LFREE30SMO) THEN
!       FREE30CORNER => CORNER_PROFILES(:,NADDCOR(30)+1)
        FREE30CORNER => CORNER_PROFILES(:,NTOT)
      ELSE
        NULLIFY(FREE30CORNER)
      END IF

      CORNER_PROFILES = 0._DP

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CORNERS



      SUBROUTINE EIRENE_DEALLOC_COMUSR
C
      IF (ALLOCATED(PLSTLS)) THEN

        DEALLOCATE (PLSTLS)
        DEALLOCATE (CEMETERYP)
c
        DEALLOCATE (TEINL)
        DEALLOCATE (TIINL)
        DEALLOCATE (DEINL)
        DEALLOCATE (DIINL)

        DEALLOCATE (NAINS)
        DEALLOCATE (NAINT)

        DEALLOCATE (ICPVE)
        DEALLOCATE (ICPVS)
        DEALLOCATE (ICPVT)
        DEALLOCATE (ICPRC)
        DEALLOCATE (IBGVE)
        DEALLOCATE (IBGVS)
        DEALLOCATE (IBGVT)
        DEALLOCATE (IBGRC)

      END IF

      IF (ALLOCATED(RMASSA)) THEN

        DEALLOCATE (RMASSA)
        DEALLOCATE (RMASSM)
        DEALLOCATE (RMASSI)
        DEALLOCATE (RMASSPH)
        DEALLOCATE (RMASSP)

        DEALLOCATE (DIOD)
        DEALLOCATE (DATD)
        DEALLOCATE (DMLD)
        DEALLOCATE (DPLD)
        DEALLOCATE (DPHD)
        DEALLOCATE (DION)
        DEALLOCATE (DATM)
        DEALLOCATE (DMOL)
        DEALLOCATE (DPLS)
        DEALLOCATE (DPHOT)

        DEALLOCATE (TEXTS)
        DEALLOCATE (DENSLIM) !VK

        DEALLOCATE (NMASSA)
        DEALLOCATE (NCHARA)
        DEALLOCATE (NFOLA)
        DEALLOCATE (NGENA)
        DEALLOCATE (NMASSM)
        DEALLOCATE (NCHARM)
        DEALLOCATE (NFOLM)
        DEALLOCATE (NGENM)
        DEALLOCATE (NMASSP)
        DEALLOCATE (NCHARP)
        DEALLOCATE (NCHRGP)
        DEALLOCATE (NMASSI)
        DEALLOCATE (NCHARI)
        DEALLOCATE (NCHRGI)
        DEALLOCATE (NFOLI)
        DEALLOCATE (NGENI)
        DEALLOCATE (NFOLPH)
        DEALLOCATE (NGENPH)
        DEALLOCATE (NPRT)
        DEALLOCATE (ISPEZ)
        DEALLOCATE (ISPEZI)

        DEALLOCATE (MPLSTI)
        DEALLOCATE (MPLSV)
        DEALLOCATE (ISPZ_BACK)
        DEALLOCATE (IADVE)
        DEALLOCATE (IADVS)
        DEALLOCATE (IADVT)
        DEALLOCATE (IADRC)
        DEALLOCATE (ICLVE)
        DEALLOCATE (ICLVS)
        DEALLOCATE (ICLVT)
        DEALLOCATE (ICLRC)
        DEALLOCATE (ISNVE)
        DEALLOCATE (ISNVS)
        DEALLOCATE (ISNVT)
        DEALLOCATE (ISNRC)
        DEALLOCATE (IADSE)
        DEALLOCATE (IADSS)
        DEALLOCATE (IADST)
        DEALLOCATE (IADSC)
        DEALLOCATE (NFRSTP)
        DEALLOCATE (NFSTPI)
        DEALLOCATE (NADDP)
        DEALLOCATE (NSPAN)
        DEALLOCATE (NSPEN)
        DEALLOCATE (NSPANW)
        DEALLOCATE (NSPENW)

        DEALLOCATE (INTLOPTS)

        DEALLOCATE (LKINDM)
        DEALLOCATE (LKINDI)
        DEALLOCATE (LKINDP)
C
        DEALLOCATE (LGVAC)
        DEALLOCATE (LSPCCLL)
        DEALLOCATE (LIVTALI)

      END IF

      IF (ALLOCATED(FLXOUT)) THEN

        DEALLOCATE (FLXOUT)
        DEALLOCATE (SAREA)

      END IF

      IF (ALLOCATED(EION)) DEALLOCATE (EION)

      IF (ALLOCATED(LSMOPRO)) DEALLOCATE (LSMOPRO)

      IF (ALLOCATED(BACK_SPEC)) DEALLOCATE (BACK_SPEC)

      IF (ALLOCATED(CORNER_PROFILES)) DEALLOCATE (CORNER_PROFILES)

      IF (ALLOCATED(KVNMXTR)) DEALLOCATE (KVNMXTR)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COMUSR


      SUBROUTINE EIRENE_INIT_COMUSR(ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (IFIRST == 0) THEN
        LSMOPRO = .FALSE.

        LTESMO      => LSMOPRO(1)
        LTISMO      => LSMOPRO(2)
        LDESMO      => LSMOPRO(3)
        LDISMO      => LSMOPRO(4)
c  plasma flow field
        LVXSMO      => LSMOPRO(5)
        LVYSMO      => LSMOPRO(6)
        LVZSMO      => LSMOPRO(7)
c  vectorial plasma flow: parallel, perp and diamagn.
        LBV_VECSMO    => LSMOPRO(23)
        LPMOM_VECSMO  => LSMOPRO(24)

c  B field
        LBXSMO      => LSMOPRO(8)
        LBYSMO      => LSMOPRO(9)
        LBZSMO      => LSMOPRO(10)
        LBFSMO      => LSMOPRO(11)
c  B  perp  unit vector
        LBXPSMO     => LSMOPRO(16)
        LBYPSMO     => LSMOPRO(17)

        LADSMO      => LSMOPRO(12)
        LEDRIFTSMO  => LSMOPRO(13)
        LVOLSMO     => LSMOPRO(14)
        LWGHTSMO    => LSMOPRO(15)

c  E field
        LEXSMO      => LSMOPRO(18)
        LEYSMO      => LSMOPRO(19)
        LEZSMO      => LSMOPRO(20)
        LEFSMO      => LSMOPRO(21)
        LPOTSMO     => LSMOPRO(22)

c  poloidal B-flux function
        LPSISMO     => LSMOPRO(25)

c  average bundle charge
        LZISMO      => LSMOPRO(26)

c  free slots
        LFREE27SMO  => LSMOPRO(27)
        LFREE28SMO  => LSMOPRO(28)
        LFREE29SMO  => LSMOPRO(29)
        LFREE30SMO  => LSMOPRO(30)

        IFIRST = 1
      ENDIF

      IF (ICAL == 1) THEN
cdr oct 18: initialization of input volumetric tallies moved to ICAL==2

        RMASSA  = 0._DP
        RMASSM  = 0._DP
        RMASSI  = 0._DP
        RMASSPH = 0._DP
        RMASSP  = 0._DP

        DIOD   = 0._DP
        DATD   = 0._DP
        DMLD   = 0._DP
        DPLD   = 0._DP
        DPHD   = 0._DP
        DION   = 0._DP
        DATM   = 0._DP
        DMOL   = 0._DP
        DPLS   = 0._DP
        DPHOT  = 0._DP

        TEXTS  = ' '
        IF (COMUSR_FIRST_PASS(1)) NMASSA = 0
        IF (COMUSR_FIRST_PASS(2)) NCHARA = 0
        NFOLA  = 0
        NGENA  = 0
        IF (COMUSR_FIRST_PASS(3)) NMASSM = 0
        NCHARM = 0
        NFOLM  = 0
        NGENM  = 0
        NMASSP = 0
        NCHARP = 0
        NCHRGP = 0
        IF (COMUSR_FIRST_PASS(4)) NMASSI = 0
        IF (COMUSR_FIRST_PASS(4)) NCHARI = 0
        IF (COMUSR_FIRST_PASS(4)) NCHRGI = 0
        NFOLI  = 0
        NGENI  = 0
        NFOLPH = 0
        NGENPH = 0
        NPRT   = 0
        ISPEZ  = 0
        ISPEZI = 0
        MPLSTI = 0
        MPLSV  = 0
        ISPZ_BACK = 0
        IADVE  = 0
        IADVS  = 0
        IADVT  = 0
        IADRC  = 0
        ICLVE  = 0
        ICLVS  = 0
        ICLVT  = 0
        ICLRC  = 0
        ISNVE  = 0
        ISNVS  = 0
        ISNVT  = 0
        ISNRC  = 0
        IADSE  = 0
        IADSS  = 0
        IADST  = 0
        IADSC  = 0
        NFRSTP = 0
        NFSTPI = 0
        NADDP  = 0
        NSPAN  = 0
        NSPEN  = 0
        NSPANW = 0
        NSPENW = 0
        INTLOPTS = 0

        LKINDM = 0
        LKINDI = 0
        IF (COMUSR_FIRST_PASS(5)) LKINDP = 0

        LGVAC   = .FALSE.
        LSPCCLL = .FALSE.
        LIVTALI = .TRUE.

        LTEIN      => LIVTALI(1)
        LTIIN      => LIVTALI(2)

        LDEIN      => LIVTALI(3)
        LDIIN      => LIVTALI(4)
cdr
        LDIN = ldein.and.ldiin
c
        LVXIN      => LIVTALI(5)
        LVYIN      => LIVTALI(6)
        LVZIN      => LIVTALI(7)
cdr
        LVIN = lvxin.and.lvyin.and.lvzin
c
        LBXIN      => LIVTALI(8)
        LBYIN      => LIVTALI(9)
        LBZIN      => LIVTALI(10)
        LBFIN      => LIVTALI(11)
        LPSI       => LIVTALI(25)
cdr
        LBIN = lbxin.and.lbyin.and.lbzin.and.lbfin  !  psi ?
c
        LADIN      => LIVTALI(12)
        LEDRIFT    => LIVTALI(13)
        LVOL       => LIVTALI(14)
        LWGHT      => LIVTALI(15)

cdr June 23: now included in LBIN  (B field components)
        LBXPERP    => LIVTALI(16)
        LBYPERP    => LIVTALI(17)
cdr     LBZPERP

        LEXIN      => LIVTALI(18)
        LEYIN      => LIVTALI(19)
        LEZIN      => LIVTALI(20)
        LEFIN      => LIVTALI(21)
        LPOT       => LIVTALI(22)
cdr
        LEIN = lexin.and.leyin.and.lezin.and.lefin  !  pot ?
c

        LBV_VEC    => LIVTALI(23)
        LPMOM_VEC  => LIVTALI(24)

        LZIIN      => LIVTALI(26)

        LFREE27    => LIVTALI(27)
        LFREE28    => LIVTALI(28)
        LFREE29    => LIVTALI(29)
        LFREE30    => LIVTALI(30)

        LDTEDX     => LIVTALI(31)
        LDTEDY     => LIVTALI(32)
        LDTEDZ     => LIVTALI(33)
        LDTIDX     => LIVTALI(34)
        LDTIDY     => LIVTALI(35)
        LDTIDZ     => LIVTALI(36)

        LDDEDX     => LIVTALI(37)
        LDDEDY     => LIVTALI(38)
        LDDEDZ     => LIVTALI(39)
        LDDIDX     => LIVTALI(40)
        LDDIDY     => LIVTALI(41)
        LDDIDZ     => LIVTALI(42)

        LDVXDX     => LIVTALI(43)
        LDVXDY     => LIVTALI(44)
        LDVXDZ     => LIVTALI(45)
        LDVYDX     => LIVTALI(46)
        LDVYDY     => LIVTALI(47)
        LDVYDZ     => LIVTALI(48)
        LDVZDX     => LIVTALI(49)
        LDVZDY     => LIVTALI(50)
        LDVZDZ     => LIVTALI(51)

        LDBXDX     => LIVTALI(52)
        LDBXDY     => LIVTALI(53)
        LDBXDZ     => LIVTALI(54)
        LDBYDX     => LIVTALI(55)
        LDBYDY     => LIVTALI(56)
        LDBYDZ     => LIVTALI(57)
        LDBZDX     => LIVTALI(58)
        LDBZDY     => LIVTALI(59)
        LDBZDZ     => LIVTALI(60)
        LDBFDX     => LIVTALI(61)
        LDBFDY     => LIVTALI(62)
        LDBFDZ     => LIVTALI(63)
        LDADINDX   => LIVTALI(64)
        LDADINDY   => LIVTALI(65)
        LDADINDZ   => LIVTALI(66)
        LDEDRIFTDX => LIVTALI(67)
        LDEDRIFTDY => LIVTALI(68)
        LDEDRIFTDZ => LIVTALI(69)
        LDVOLDX    => LIVTALI(70)
        LDVOLDY    => LIVTALI(71)
        LDVOLDZ    => LIVTALI(72)
        LDWGHTDX   => LIVTALI(73)
        LDWGHTDY   => LIVTALI(74)
        LDWGHTDZ   => LIVTALI(75)
        LDBXPERPDX => LIVTALI(76)
        LDBXPERPDY => LIVTALI(77)
        LDBXPERPDZ => LIVTALI(78)
        LDBYPERPDX => LIVTALI(79)
        LDBYPERPDY => LIVTALI(80)
        LDBYPERPDZ => LIVTALI(81)
        LDEXDX     => LIVTALI(82)
        LDEXDY     => LIVTALI(83)
        LDEXDZ     => LIVTALI(84)
        LDEYDX     => LIVTALI(85)
        LDEYDY     => LIVTALI(86)
        LDEYDZ     => LIVTALI(87)
        LDEZDX     => LIVTALI(88)
        LDEZDY     => LIVTALI(89)
        LDEZDZ     => LIVTALI(90)
        LDEFDX     => LIVTALI(91)
        LDEFDY     => LIVTALI(92)
        LDEFDZ     => LIVTALI(93)
        LDPOTDX    => LIVTALI(94)
        LDPOTDY    => LIVTALI(95)
        LDPOTDZ    => LIVTALI(96)
cdr  vectorial input tallies:  par, perp, dia
        LDBV_VECDX   => LIVTALI(97)
        LDBV_VECDY   => LIVTALI(98)
        LDBV_VECDZ   => LIVTALI(99)
        LDPMOM_VECDX => LIVTALI(100)
        LDPMOM_VECDY => LIVTALI(101)
        LDPMOM_VECDZ => LIVTALI(102)

        LDPSIDX    => LIVTALI(103)
        LDPSIDY    => LIVTALI(104)
        LDPSIDZ    => LIVTALI(105)

        LDZIDX     => LIVTALI(106)
        LDZIDY     => LIVTALI(107)
        LDZIDZ     => LIVTALI(108)

        LDFREE27DX => LIVTALI(109)
        LDFREE27DY => LIVTALI(110)
        LDFREE27DZ => LIVTALI(111)
        LDFREE28DX => LIVTALI(112)
        LDFREE28DY => LIVTALI(113)
        LDFREE28DZ => LIVTALI(114)
        LDFREE29DX => LIVTALI(115)
        LDFREE29DY => LIVTALI(116)
        LDFREE29DZ => LIVTALI(117)
        LDFREE30DX => LIVTALI(118)
        LDFREE30DY => LIVTALI(119)
        LDFREE30DZ => LIVTALI(120)

      ELSE IF (ICAL == 2) THEN
c  Active volumetric input tallies
        PLSTLS = 0._DP
c  Cemetery for inactive input tallies (no storage)
        CEMETERYP = 0._DP

        TEINL  = 0._DP
        TIINL  = 0._DP
        DEINL  = 0._DP
        DIINL  = 0._DP

        NAINS  = 0
        NAINT  = 0

        ICPVE  = 0
        ICPVS  = 0
        ICPVT  = 0
        ICPRC  = 0
        IBGVE  = 0
        IBGVS  = 0
        IBGVT  = 0
        IBGRC  = 0

      ELSE IF (ICAL == 3) THEN

        FLXOUT = 0._DP
        SAREA  = 666._DP

      ELSE IF (ICAL == 4) THEN

        EION   = 0._DP

      END IF

      RETURN
      END SUBROUTINE EIRENE_INIT_COMUSR

      SUBROUTINE EIRENE_COMUSR_REINIT
      IMPLICIT NONE
      IFIRST = 0
      RETURN

      END SUBROUTINE EIRENE_COMUSR_REINIT


      SUBROUTINE EIRENE_BROADCAST_COMUSR(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) THEN
        CALL EIRENE_ALLOC_COMUSR(0)
        CALL EIRENE_ALLOC_COMUSR(1)
        CALL EIRENE_ALLOC_COMUSR(2)
        CALL EIRENE_ALLOC_COMUSR(3)
        CALL EIRENE_ALLOC_COMUSR(4)
      END IF

! LSMOPRO needs to be broadcasted before corner arrays are allocated
      CALL MPI_BCAST (LSMOPRO,NTALG,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
      IF (ME /= 0) CALL EIRENE_ALLOC_CORNERS

c  active and inactive tallies:
c  INPUT:
      CALL MPI_BCAST (LIVTALI,NTALI,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
cdr   intlopts is only needed on processor 0

      CALL MPI_BCAST (NADDP,NTALI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFSTPI,NTALI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL EIRENE_ASSOCIATE_COMUSR
      CALL MPI_BCAST (PLSTLS,NINPTL*NRAD,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (TEINL,NRAD,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TIINL,NPLSTI*NRAD,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DEINL,NRAD,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DIINL,NPLS*NRAD,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (FLXOUT,NLMPGS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (SAREA,NLMPGS,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (EION,NATM,MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (RMASSI,NION,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RMASSA,NATM,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RMASSM,NMOL,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RMASSP,NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DIOD,NION,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DATD,NATM,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DMLD,NMOL,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DPLD,NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DION,NION,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DATM,NATM,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DMOL,NMOL,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DPLS,NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DPHOT,NPHOT,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (TVAC,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (DVAC,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (VVAC,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ZVAC,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ALLOC,1,MPI_REAL8,0,MPI_COMM_WORLD,ier)

! smoothed input tallies
      CALL MPI_BCAST (CORNER_PROFILES,SIZE(CORNER_PROFILES),
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (TEXTS,8*NSPZ,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NATMI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NATMIM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPAM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMOLI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMOLIM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPAMI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIONI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NIONIM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPTOT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPLSI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPLSIM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPH,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHOTI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPHOTIM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSNVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCPVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NADVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NBGVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NALVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCLVI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NADSI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NALSI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NAINI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NBITS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMASSA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHARA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFOLA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGENA,NATM,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMASSM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHARM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFOLM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGENM,NMOL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMASSP,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHARP,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHRGP,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMASSI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHARI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NCHRGI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFOLI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGENI,NION,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFOLPH,NPHOT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NGENPH,NPHOT,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NPRT,NSPZ,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISPEZ,6*(1+NPHOTP)*(NATMP+1)*(NMOLP+1)*(NIONP+1)*
     .                      (NPLSP+1),MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISPEZI,6*NSPZ,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (MPLSTI,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (MPLSV,NPLS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADVE,NADV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADVS,NADV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADVT,NADV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADRC,NADV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICLVE,NCLV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICLVS,NCLV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICLVT,NCLV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICLRC,NCLV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISNVE,NSNV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISNVS,NSNV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISNVT,NSNV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ISNRC,NSNV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICPVE,NCPV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICPVS,NCPV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICPVT,NCPV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICPRC,NCPV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGVE,NBGV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGVS,NBGV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGVT,NBGV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IBGRC,NBGV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADSE,NADS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADSS,NADS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADST,NADS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IADSC,NADS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFRSTP,NTALI,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPAN,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NSPEN,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LGVAC,NRAD*(NPLS+2),MPI_LOGICAL,0,MPI_COMM_WORLD,
     .                ier)
      CALL MPI_BCAST (DENSLIM,NPLS,MPI_REAL8,0,MPI_COMM_WORLD,ier) !VK
      CALL MPI_BCAST (NPRLL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NMODE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NTCPU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILEN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILEM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILEL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILEK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFILEJ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NITER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (IITER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NTIME,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ITIMV,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      IF (ALLOCATED(NAINS)) THEN
        CALL MPI_BCAST (NAINS,NAIN,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        CALL MPI_BCAST (NAINT,NAIN,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_COMUSR

      END MODULE EIRMOD_COMUSR
