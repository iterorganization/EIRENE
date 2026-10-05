!PB   17.11.05  NEMODS=-2  =>  NEMODS=8
!PB             NEMODS=-3  =>  NEMODS=9
!PB   14.07.11  OVERWRITE NEMODS SPECIFIED BY THE EIRENE INPUT BY
!PB             SETTINGS CONSISTENT WITH SHEATH SETTINGS IN B2
!PB             IF SHEATH PARAMETERS ARE GIVEN BY B2 NEMODS=9
!PB             IF NO SHEATH PARAMETERS WERE TRANSFERRED BY B2 NEMODS=8
CDR  09. 2014 GENERAL RELATIONS BETWEEN FINE (UNSTRUCTURED) AND COARSE (STRUCTURED)
C             GRID MADE MORE EXPLICIT. NEW NSBOX FOR FINE GRID SET.
C             NP2NDQ REMOVED (REDUNDANT: =NP2TAL)
C             NR1STQ REMOVED (REDUNDANT: =NR1TAL),
C             AND ERROR: NR1STQ WAS USED BEFORE DEFINITION --> PROBLEMS WITH NSTGRD ARRAYS
C             VIA FILES FROM FORT.29?

c             plus minor notational cleanup, comments added

cdr 150407:  orientation of B field made optional, additional input
cdr           flags ibrad,ibpol,ibtor in block 14.
CDR 150419    THIS ROUTINE WAS OBTAINED BY MERGING COUPLE_B2.5 AND COUPLE_TRIA
cdr           from FZJ repositories at 2011.
cdr           now: synchronize again with couple_Tria from 2015 master branch started:
cdr           comments, cleanup, nomenclature --> ITRI  (loops 1111,....)

cdr 150409          magnitude of B field (T) transferred.
cdr                 to be checked: orientation of uudiag for reconstruction of cartesian
cdr                 flow velocity components.
cdr                 apparently not used: vvdiag.


cpb  15.09.15: added: default B field =1 (Tesla), if bfield=0, cell-wise.

cdr  start to use species-resolved energy tallies.
cdr  nov.15:  eapl,empl,eipl: now species-resolved
cdr           --> eppl_cop, eploda get an additional species index ipls
C
cdr  dec. 15: not ready. Started to comment, and to extend copv tallies
cdr           for total and internal energy with species index. Not ready...

cdr  tbd:  check storage on copv, whenever used.

cdr  jan.17
cdr  epel(nrad)      --> epel_cop(nrad)
cdr  cppv(ncpv,nrad) --> mppl_cop(npls,nrad), for notational consistency
c    new: write file of triang. data for idl plotting tool.

cdr  to be done: check handling of paat and pppl_cop for wrneutral.
cdr              which radiation estimate?
cdr              should involve only the net source?

c......................................................................
CDR JAN 17:   re-synchronize with couple_tria, identify differences:

c                    additional here:  mpi stuff from S. Wiesen
c                    chpm,...allocatable
c                    lcut in common, and broadcast.

cdr  NFL dependence:  SFNIT(0:NSTEP,NFL) in global particle balance, already implemented

c                    includes 2014 attempt to invert B field by block 14 input
c                    ncpvi tallies different, probably redundant anyway
c                    double grid structure (tria and polyg.) options not cleaned up
cdr Feb 17:  remove duplicated code re call to geousr.
c            remove unused arrays eapl0,empl0,eipl0
c            double grid structure: triangles and polygons, for tallies.
c            parameters nr1st,  np2nd,  ..... for geometry
c            parameters nr1tal, np2tal, ..... for scoring
c            parameters nr1tal_save, ....     for interfacing tallies between b2 and eirene
c                                             is always the b2 (structured) coarse grid
cdr Jan 18 : bug fix re vol.rec., only one ipls per stratum is supported
c            code was correct in solps4.3, and garching versions of couple_b2/b2.5

cdr March 18: new variable LCOARSE: maintain underlying coarse structured grid, scoring
cdr           on coarse grid (NCLTAL array). Otherwise: only fine (triangular) grid structure
cdr Mar 18:  ELTEST from couple_Tria
cdr Jan 22:  vectorial momentum exchange tallies MAPL_VEC....
cdr          received from EIRENE.
cdr          ... plus some notational cleanup
c......................................................................................


C   EIRENE CODE SEGMENT COUPLE_$, $ MAY CURRENTLY STAND FOR B2,
C                                                           B2.5 (SOLPS-ITER),
C                                                           DIVIMP,
C                                                           TRIA,
C                                                           TETRA,
C                                                           TRANSP,
C                                                           DUMMY
C
C   THIS VERSION: $COUPLE_SOLPS-ITER/$COUPLE_TRIA
c                 combined by s.wiesen@fz-juelich.de, from fzj repositories at 2011
c                 additionally modified for MPI use, s.wiesen@fz-juelich.de, apr2011
C
c  geometry data not any longer via work array into eirene
c                due to module structure
c  eliminate cut cells from balances (lcut(..))
c  removed: ncopib, ncopeb

C   UPDATES:
C   OPTION TO EVALUATE B FIELD VECTORS FROM GRIDADAP FILE FT29
C   FOR NON-ORTHOGONAL GRIDS
C
C   THIS CODE SEGMENT CONTAINS VARIOUS SUBROUTINES NEEDED FOR
C   INTERFACING THE EIRENE CODE TO PLASMA FLUID CODES.
C   IT READS GEOMETRICAL DATA (MESHES) FROM FILE FT30
C   AND PRODUCES THE EIRENE INPUT DATA (BLOCK 2).
C   IT READS PLASMA BACKGROUND DATA FROM FILE FT31 OR COMMON BLOCKS,
C   IT MAY (OPTIONAL) ALSO READ PLASMA DATA FROM FILE FT13
C   WRITTEN IN A PREVIOUS EIRENE RUN (E.G. IN ORDER TO ITERATE
C   IN SOME BACKGROUND SPECIES)
C   IT THEN PRODUCES INPUT DATA FOR EIRENE
C   INPUT BLOCK 5 (PLASMA DATA) AND BLOCK 7 (SURFACE RECYCLING SOURCES)
C
C
C   THIS PARTICULAR VERSION LINKS EIRENE TO THE B2.5 2D MULTIFLUID EDGE
C   PLASMA TRANSPORT CODE.
C
C   IT WAS WRITTEN BY D.REITER AND P.BOERNER, FZ-JUELICH
C   E-MAIL: D.REITER @ FZ-JUELICH.DE, AND: www.eirene.de
C
C
C
C   MOST OF THE FORTRAN IN THIS CODE SEGMENT HAS BEEN DEVELOPED
C   UNDER KFA-NET CONTRACT NO. 428/90-8/FU-D
C
C   FINAL REPORT BY: D.REITER(1), P.BOERNER(1), B.KUEPPERS(1),
C                    M.BAELMANS(2) AND G.P.MADDISON(3)
C                    (1992)
C   1: KFA-JUELICH GMBH
C   2: UNIV. LEUVEN, ERM, KFA-JUELICH
C   3: AEA TECHNOLOGY, FUSION, CULHAM, UKAEA FUSION ASSOCIATION
C
*DK COUPLE
C
      MODULE EIRMOD_INFCOP
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_BRASPOI
      USE EIRMOD_CADGEO
      USE EIRMOD_CCOUPL
      USE EIRMOD_CESTIM
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COUTAU
      USE EIRMOD_CPOLYG
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRID
      USE EIRMOD_CINIT
      USE EIRMOD_CLGIN
      USE EIRMOD_COMXS
      USE EIRMOD_CPLOT
      USE EIRMOD_CSDVI
      USE EIRMOD_CSPEI
      USE EIRMOD_CSPEZ
      USE EIRMOD_CSTEP
      USE EIRMOD_COMNNL
      USE EIRMOD_CTEXT
      USE EIRMOD_CTRIG
C  PLASMA DATA: NI,TE,TI,VV,UU,PR,UP,RR,FNIX,FNIY.. (BRAAMS ---> EIRENE)
      USE EIRMOD_BRAEIR
C  NEUTRAL SOURCE TERMS: SNI,SMO,SEE,SEI (EIRENE ---> BRAAMS)
      USE EIRMOD_EIRBRA
      USE EIRMOD_BRASCL
      USE EIRMOD_CAI
      USE EIRMOD_WNEUTRALS
      USE EIRMOD_EXTRAB25
      USE EIRMOD_MPI
      USE EIRMOD_CPES
      USE EIRMOD_SHEATH
      USE EIRMOD_OPENFILE, ONLY: EIRENE_OPENFILE
#ifndef NO_JSON
      USE EIRMOD_JSON

      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck
#endif
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: EIRENE_IF0COP, EIRENE_IF1COP, EIRENE_IF2COP,
     .          EIRENE_IF3COP, EIRENE_IF4COP, EIRENE_IF3COP_SUM,
     .          EIRENE_BROADCAST_INFCOP, EIRENE_BROADCAST_IF3COP,
     .          EIRENE_DEALLOC_INFCOP, EIRENE_INFCOP_PRE_MCARLO,
     .          EIRENE_INFCOP_POST_STRATUM, EIRENE_INFCOP_PRE_STRATA
C
      TYPE :: CELL
        INTEGER :: TRIANGLE
        TYPE(CELL), POINTER :: NEXT
      END TYPE CELL

      TYPE :: POIFELD
        TYPE (CELL), POINTER :: P
      END TYPE POIFELD

      TYPE :: FACE
        INTEGER :: TRIANGLE, SIDE
        TYPE(FACE), POINTER :: NEXT
      END TYPE FACE

      TYPE :: FACE_POINTER
        TYPE(FACE), POINTER :: FP
      END TYPE FACE_POINTER

      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPS(:), SNIS(:), CHMOS(:), SMOS(:), SCALN(:),
     . SEES0(:), SEIS0(:), SNIS0(:,:), SMOS0(:,:),
     . RESSNI(:,:), RESSMO(:,:),
     . RESSEE(:), RESSEI(:),
     . FLXEIR(:)
      REAL(DP), ALLOCATABLE, SAVE ::
     . CHPM(:,:), CHEEM(:), CHEIM(:), CHMOM(:,:)

C pppl_cop, mppl_cop, eppl_cop and epel_cop are the exact
c volumetric primary source tallies,
c while default eirene tallies pppl, mppl, eppl and epel would be
c the corresponding tallies scored from random sampling in eirene
c Aug 22: ...both: mppl and epel tallies do not exist so far.
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPPL_COP(:,:), MPPL_COP(:,:),
     .            EPPL_COP(:,:), EPEL_COP(:)
C
c  for short cycle correction terms, in vol. rec. strata.
      REAL(DP), ALLOCATABLE, SAVE ::
     .            PPLODA(:,:), MPLODA(:,:),
     .            EPLODA(:,:), EPEODA(:)

      real(dp), allocatable :: save_estimv(:,:), save_estims(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: eapl0(:),empl0(:),eipl0(:)
      REAL(DP), SAVE :: V, VL, BRAD, BPOL, BTOR

      INTEGER, SAVE :: NPBS, NPES
      INTEGER, ALLOCATABLE, SAVE :: NPTS_SAVE(:)
      LOGICAL, ALLOCATABLE, SAVE :: NPTS_SAVED(:)
      LOGICAL, SAVE, PUBLIC :: LSHORT
      LOGICAL, SAVE, PRIVATE :: LSTOP, LSTP3, LTEST
#ifdef WG_TODO
      LOGICAL :: LXSRF
#endif
      TYPE(FACE_POINTER), ALLOCATABLE, SAVE :: NCLFC(:)
      integer, allocatable, save :: invncltal(:)

      CONTAINS

C
C     THIS SUBROUTINE DEFINES THE PLASMA MODEL IN CASE OF A COUPLED
C     NEUTRAL-PLASMA CALCULATION
C
C     THE SUBROUTINE "IF0COP" RECEIVES GEOMETRICAL INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE SUBROUTINE "IF1COP" RECEIVES PLASMA INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE SUBROUTINE "IF2COP" PREPARES THE SOURCE SAMPLING DISTRIBUTION
C     FROM THE EXTERNAL DATA, AND MAY OVERWRITE OTHER INPUT
C     DATA FROM BLOCKS 1 TO 13 AS WELL
C
C     THE SUBROUTINES "IF3COP, IF4COP" RETURN RESULTS TO AN EXTERNAL CODE
C

      SUBROUTINE EIRENE_IF0COP(LFIXED, LSHORT)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: LFIXED, LSHORT
C
#ifdef WG_TODO
      INTEGER :: IP, IXNI
      INTEGER, ALLOCATABLE :: XNMTI_BUF(:), YNMTI_BUF(:), TNMTI_BUF(:)
#endif
      INTEGER, SAVE :: IO,
     I           dum7,dum8,dum9,dum10
      INTEGER :: N, I, J, K, L, JPLS, JUN,
     I           IS, IS1, ISC, ISC1, ISC2, ISCS, ISTS, ISTR, IT,
     I           ICOADD, ICOG, IDUMMY, IERROR,
     I           ITRI, NBAR, NBARSI, ICOFC, IC,
     I           IFC1, IFC2, IFC3, ICC,
     I           IUNIN_SAVE, ID3b, ID3bc
#ifndef NO_JSON
      INTEGER :: JS
#endif
      REAL(DP) :: TX, TY, VS, VT,
     R            VSX, VSY, VTX, VTY, VPRO, tx1, ty1, tx2, ty2
      INTEGER, ALLOCATABLE :: JDUMMY(:)
C
      CHARACTER(72) :: ZEILE
      character(100) :: ioerrmsg
      LOGICAL :: L1, L2, FOUND, IFBOUND
      EXTERNAL :: EIRENE_HEADNG, EIRENE_LEER, EIRENE_EXIT_OWN,
     .            EIRENE_GEOMD, EIRENE_GEOUSR,
     .            EIRENE_INDMAP, EIRENE_PLASM, EIRENE_SNEIGH,
     .            EIRENE_MASJ1, EIRENE_MASJ2, EIRENE_MASJ5
#ifdef NO_JSON
      EXTERNAL :: EIRENE_MASAGE
#endif
C
      IERROR=0
      IUSROUT = 0
      IUNIN_SAVE = IUNIN
C
      IF (.NOT.ALLOCATED(CHPM)) THEN
cdr  for short cycle only. Obsolete
        ALLOCATE (CHPM(NPLS,NRAD))
        ALLOCATE (CHEEM(NRAD))
        ALLOCATE (CHEIM(NRAD))
        ALLOCATE (CHMOM(NPLS,NRAD))  !  should be made vectorial tally

cdr  for direct (primary) contributions from vol.rec sources
cdr  Always needed for vol rec strata. Not only for short cycle
        ALLOCATE (PPPL_COP(NPLS,NRAD))
        ALLOCATE (MPPL_COP(3*NPLS,NRAD))
        ALLOCATE (EPPL_COP(NPLS,NRAD))
        ALLOCATE (EPEL_COP(NRAD))

cdr  for short cycle only. Obsolete
        ALLOCATE (PPLODA(NPLS,NRAD))
        ALLOCATE (MPLODA(3*NPLS,NRAD))
        ALLOCATE (EPLODA(NPLS,NRAD))
        ALLOCATE (EPEODA(NRAD))
      END IF
cdr
      NLSHRT13 = .TRUE.  !  only short version of fort13 is used:
                         !  calls WRPLAM_SHRT, RPLAM_SHRT
      ntrfrm = 0

      if (.not.allocated(eapl0)) then
        allocate (eapl0(nrad))
        allocate (empl0(nrad))
        allocate (eipl0(nrad))
      end if
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        IF (LFIXED) THEN
          CALL EIRENE_READ14_FIXED
        ELSE
#ifndef NO_JSON
          js = itree_num(14)
          CALL EIRENE_READ14_JSON(jtrees(js),blks(14)%p)
#else
          CALL EIRENE_MASAGE
     .     ('JSON FORMAT NOT SUPPORTED IF COMPILED WITH NO_JSON OPTION')
          CALL EIRENE_EXIT_OWN(1)
#endif
        ENDIF
      ENDIF
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C
C DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                             AND IN SUBR. COLLIDE)

cdr  already done in if0prm. Hidden link, must be removed....
      NCPVI = NPLSI
      NCPV  = MAX(NCPV,NCPVI)
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      NFL = NFLA
C
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTR=1,NSTRAI
        LNLVOL(ISTR)=NLVOL(ISTR)
   60 CONTINUE
      NMODEI=NMODE
      NFILNN=NFILEN
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      IF (NCPVI.EQ.0) GOTO 70

CDR  SET THE NCPVI= NPLSI COUPLE TALLIES

      DO JPLS=1,NPLSI
        ICPVE(JPLS)=1
        ICPRC(JPLS)=1
        TXTTAL(JPLS,NTALM)=
     .  'ENERGY-WEIGHTED CX RATE OF ATOMS WITH IPLS                  '
        TXTSPC(JPLS,NTALM)=TEXTS(NSPAMI+JPLS)
        TXTUNT(JPLS,NTALM)='AMP                     '
      ENDDO
C
   70 CONTINUE
C
C
C  TRANSFER GEOMETRY
C
      IF (.NOT.(INDGRD(1).EQ.6.OR.INDGRD(2).EQ.6.OR.INDGRD(3).EQ.6))THEN
        IUNIN = IUNIN_SAVE
        IF (IUSROUT /= 0) CLOSE(IUSROUT)
        RETURN
      END IF
C
C  READ IN DATA TO SET UP GEOMETRY FOR NEUTRAL GAS TRANSPORT CODE
C  STATEMENT NUMBER 1000 ---> 1999
C
C 1000 CONTINUE
C
      CALL EIRENE_LEER(1)
C
C  TRANSFER FLAGS
C
      NAINI=NAINB
C
      nrknot=0
      ntrii=0
      if(.not.nltrimesh) then !{
CTRIG A
C  READ DATA FOR TRIANGULAR MESH
C
CVK I ADDED A FORMAT STRING FOR READING FROM FORT.33-35 THE SAME LIKE IN TRIAGEOM
        OPEN (UNIT=33,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1983)
        OPEN (UNIT=34,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1984)
        OPEN (UNIT=35,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1985)
        REWIND 33
        REWIND 34
        REWIND 35
C
        READ(33,*,IOSTAT=IO,ERR=1983) NRKNOT  !VK IO
        IF(IO.NE.0) NRKNOT=0 !VK
        WRITE(iunout,'(a14,i8)') 'NRKNOT      = ',NRKNOT
C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
        READ(34,*,IOSTAT=IO,ERR=1984) NTRII   !VK IO
        IF(IO.NE.0) NTRII=0 !VK
        WRITE(iunout,'(a14,i8)') 'NTRII       = ',NTRII

        READ (35,*,IOSTAT=IO,ERR=1985) IDUMMY
        IF(IO.NE.0) IDUMMY=0 !VK
        IF (IDUMMY /= NTRII) THEN
          WRITE (IUNOUT,*) ' NUMBERS OF TRIANGLES DO NOT MATCH'
          WRITE (IUNOUT,*) ' IN FILES ...ELEMENTE AND ...NEIGHBORS'
          WRITE (IUNOUT,*) ' CHECK THE GEOMETRY, EXIT CALLED'
          CALL EIRENE_EXIT_OWN(1)
        END IF

C
C  EACH ELEMENT (TRIANGLE) IS GIVEN BY 3 POINTS
C
C                   3
C                 /   \
C            3  /       \  2
C             /           \
C           /               \
C          1.................2
C                  1
C
C
        IF (NTRFRM == 0) THEN
          READ(33,*,IOSTAT=IO) (XTRIAN(I),I=1,NRKNOT)
          READ(33,*,IOSTAT=IO) (YTRIAN(I),I=1,NRKNOT)
          IF(IO.GT.0) THEN
C TO READ THE FILE PRODUCED BY DOS
            WRITE(iunout,*) "Cannot read "//fort_lc//"33 in regular way"
            WRITE(iunout,*) "Trying to read line by line..."
            REWIND 33
            READ(33,*) NRKNOT
            N=INT((NRKNOT+3)/4)
            DO I=1,N
              J=(I-1)*4+1
              READ(33,*,IOSTAT=IO) XTRIAN(J:J+4)
            END DO
            DO I=1,N
              J=(I-1)*4+1
              READ(33,*,IOSTAT=IO) YTRIAN(J:J+4)
            END DO
            WRITE(iunout,*) "... ready!"
          END IF
          if (plidl) then
c write file 'triang_new.npco_char' for triang-grid, for idl tool, in appropriate format.
c first: fetch a free file unit number
            jun=-9999
            call eirene_openfile(jun,file='triang_new.npco_char',
     .           access='SEQUENTIAL',form='FORMATTED')
c next: write file 'triang_new.npco_char'
            write (jun,'(i9)') NRKNOT
            DO I=1,NRKNOT
              WRITE(jun,'(i9,2es24.16)') I,XTRIAN(I),YTRIAN(I)
            ENDDO
            close (unit=jun)
          end if
        ELSE
          DO I=1,NRKNOT
            READ(33,*) J,XTRIAN(I),YTRIAN(I)
          ENDDO
        END IF
C
        IF (NTRII.GT.NRAD.OR.NTRII.GT.NTRI) THEN
          WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFUSR '
          WRITE (iunout,*) ' NTRII MUST BE < NRAD AND <= NTRI'
          WRITE (iunout,*) ' NTRII,NRAD,NTRI = ',NTRII,NRAD,NTRI
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

        IXTRI = 0
        IFCTRI = 0

        DO I=1,NTRII
          READ(34,*) J,NECKE(1,I),NECKE(2,I),NECKE(3,I)
        ENDDO

        DO I=1,NTRII
cdr  june 17:
cdr: careful: I ne J possible. Unless triangles are sorted as J= 1,2,3... on fort.35
cdr           This is implicitly assumed here ?
cpb           J is of no consequence here. It could have been called IDUM as well.
          READ(35,*) J,NCHBAR(1,I),NSEITE(1,I),INMTI(1,I),
     >                 NCHBAR(2,I),NSEITE(2,I),INMTI(2,I),
     >                 NCHBAR(3,I),NSEITE(3,I),INMTI(3,I),
     >                 IXTRI(I),IFCTRI(1,I),IFCTRI(2,I),IFCTRI(3,I)

C         WRITE(iunout,*) J,NECKE(1,J),NECKE(2,J),NECKE(3,J),
C    >                 NCHBAR(1,I),NSEITE(1,I),
C    >                 NCHBAR(2,I),NSEITE(2,I),NCHBAR(3,I),NSEITE(3,I)

C THE SPECIAL SURFACE PROPERTY (IF ANY) IS ON INMTI ARRAY, AND TRANSFERRED INTO
C EIRENE VIA COMMON.
        ENDDO

        CLOSE (UNIT=33)
        CLOSE (UNIT=34)
        CLOSE (UNIT=35)
        CALL EIRENE_LEER(1)
      end if !}

      IF (MAXVAL(IXTRI).GT.NCIP .OR. MAXVAL(IFCTRI).GT.NFCP) THEN
        WRITE(IUNOUT,*)
     .   'Grid size implied in '//fort_lc//'35 is too large!'
        CALL EIRENE_MASJ2 ('max(IXTRI),NCIP=',MAXVAL(IXTRI),NCIP)
        CALL EIRENE_MASJ2 ('mx(IFCTRI),NFCP=',MAXVAL(IFCTRI),NFCP)
        CALL EIRENE_EXIT_OWN(1)
      END IF

C
C
C  DETERMINE THE ARRAYS INMTI AND INSPAT FOR ALL NON-DEFAULT STD. SURFACES
C  ISTS=INMTI(ISIDE,NRCELL), ISIDE=1, 2, OR 3
C
#ifdef WG_TODO
      ALLOCATE(XNMTI_BUF(NTRII*3),YNMTI_BUF(NTRII*3),TNMTI_BUF(NTRII*3))
#endif

      inmtinss(1:3,1:ntrii) = sign(1,ifctri(1:3,1:ntrii))
      ifctri(1:3,1:ntrii) = abs(ifctri(1:3,1:ntrii))

      ICOG = 0
      DO ISTS = 1, NSTSI
        DO ITRI = 1, NTRII
          DO IS = 1, 3
            IF (INMTI(IS,ITRI) == -ISTS) THEN
              ICOG = ICOG + 1
              INSPAT(IS,ITRI) = ICOG
              INMTI(IS,ITRI) = NLIM + ISTS
              NBAR = NCHBAR(IS,ITRI)
#ifdef WG_TODO
              XNMTI_BUF(ICOG)=IP !VK
              YNMTI_BUF(ICOG)=IXNI-1 !VK
              TNMTI_BUF(ICOG)=1 !VK
#endif
              IF (NBAR > 0) THEN
                NBARSI = NSEITE(IS,ITRI)
                ICOG = ICOG + 1
                INSPAT(NBARSI,NBAR) = ICOG
                INMTI(NBARSI,NBAR) = NLIM + ISTS
#ifdef WG_TODO
                XNMTI_BUF(ICOG)=IP !VK
                YNMTI_BUF(ICOG)=IXNI-1 !VK
                TNMTI_BUF(ICOG)=1 !VK
                LXSRF=.FALSE.
                CALL CORRECTNSS !VK
#endif
              END IF
            END IF
          END DO
        END DO
      END DO
#ifdef WG_TODO
      NNNMTI=ICOG
#endif
C
C  NOW THE ADJUSTMENTS, WHICH ARE AUTOMATICALLY DONE IN GEOUSR OTHERWISE
C  (INPUT BLOCK 15, B2.5-CODE SPECIFIC)
C

      CALL EIRENE_GEOUSR

      ! Do further adjustments based on fort.36 file to deal with 
      ! precision effects of 3b surfaces due to void polygon 
      ! vertex coordinate precision
      open(unit=36, file='fort.36', status='old', action='read', 
     .    iostat=io, iomsg=ioerrmsg)
      if (io == 0) then 
        write(iunout, *) 'EIRENE_IF0COP: Found fort.36 file, reading...'
        ! File found, read
        ID3bc = 0
        do while (.true.)
          ! Read the next line
          read(36, '(72a)', iostat=io, iomsg=ioerrmsg) ZEILE 
          if (io /= 0) then 
            print *, 'iostat: ', io
            print *, 'ioerrmsg:', ioerrmsg
            print *, ZEILE 
            exit 
          end if 
          backspace(36)
          read(36, *, iostat=io) ID3b, tx1, ty1, tx2, ty2 

          ! Check if read was successful
          if (io /= 0) then 
            write(iunout, *) 'EIRENE_IF0COP: stopped reading after '
            write(iunout, *) ID3bc, 'surfaces read'
            write(iunout, *) 'io status: ', io
            write(iunout, *) 'io msg: ', ioerrmsg
            exit 
          end if 

          ! Change this surface
          write(iunout, *) 'EIRENE_IF0COP: modifying 3b surface ', ID3b 
          write(iunout, *) 'before (x1, y1, x2, y2): ', p1(1, ID3b), 
     .      p1(2, ID3b), p2(1, ID3b), p2(2, ID3b)

          p1(1, ID3b) = tx1 
          p1(2, ID3b) = ty1 
          p2(1, ID3b) = tx2 
          p2(2, ID3b) = ty2

          write(iunout, *) 'after (x1, y1, x2, y2)', p1(1, ID3b), 
     .      p1(2, ID3b), p2(1, ID3b), p2(2, ID3b)
        end do 
      else
        ! File not found, skip
        write(iunout, *) 'EIRENE_IF0COP: Could not open fort.36 file' //
     .   ', continuing without mapping void polygon vertices...'
        write(iunout, *) 'iomsg: ', ioerrmsg
      end if 
      close(unit=36)

C
C  DETERMINE THE ARRAY INMTI FOR ALL ADDITIONAL SURFACES
C  ISTS=INMTI(ISIDE,NRCELL), ISIDE=1, 2, OR 3
C
      DO I=1,NLIMI
        IF (IGJUM0(I)==0) THEN
C  SURFACE I IS ACTIVE
          IF (ILPLG(I).NE.0) THEN
C  SURFACE I IS PART OF A CONTOUR USED FOR THE MESH GENERATOR
            VSX=P2(1,I)-P1(1,I)
            VSY=P2(2,I)-P1(2,I)
            VS=SQRT(VSX**2+VSY**2)+EPS60
C
            DO 1111 ITRI=1,NTRII
              DO IS=1,3
                IF (IS.EQ.1) THEN
                  VTX=XTRIAN(NECKE(2,ITRI))-XTRIAN(NECKE(1,ITRI))
                  VTY=YTRIAN(NECKE(2,ITRI))-YTRIAN(NECKE(1,ITRI))
                  ISCS=1
                  ISC1=1
                  ISC2=2
                ELSEIF (IS.EQ.2) THEN
                  VTX=XTRIAN(NECKE(3,ITRI))-XTRIAN(NECKE(2,ITRI))
                  VTY=YTRIAN(NECKE(3,ITRI))-YTRIAN(NECKE(2,ITRI))
                  ISCS=2
                  ISC1=2
                  ISC2=3
                ELSEIF (IS.EQ.3) THEN
                  VTX=XTRIAN(NECKE(1,ITRI))-XTRIAN(NECKE(3,ITRI))
                  VTY=YTRIAN(NECKE(1,ITRI))-YTRIAN(NECKE(3,ITRI))
                  ISCS=3
                  ISC1=3
                  ISC2=1
                ENDIF
                VT=SQRT(VTX**2+VTY**2)+EPS60
C
                VPRO=(VSX*VTY-VTX*VSY)/(VT*VS)
                IF (ABS(VPRO).LT.1.E-2) THEN
C  SURFACES ARE PARALLEL,
C  TEST IF ONE POINT OF THE APPROPRIATE TRIANGLE SIDE BELONGS TO
C  THE SURFACE
                  ISC=ISC1
                  TX=XTRIAN(NECKE(ISC,ITRI))
                  TY=YTRIAN(NECKE(ISC,ITRI))

                  l1=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                              p2(1,i),p2(2,i))
                  TX=XTRIAN(NECKE(ISC2,ITRI))
                  TY=YTRIAN(NECKE(ISC2,ITRI))
                  l2=eirene_point_on_interval(tx,ty,p1(1,i),p1(2,i),
     .                                              p2(1,i),p2(2,i))
                  if(l1 .and. l2) then
                    IGJUM0(I)=1
                    ICOG=ICOG+1
                    INSPAT(ISCS,ITRI)=ICOG
                    INMTI(ISCS,ITRI)=I
                  endif

                ENDIF
              ENDDO
 1111       CONTINUE  ! END OF ITRI LOOP
          ENDIF
        ENDIF
      ENDDO

C  SET UP MAPPING ARRAY USED FOR DIRECTLY SCORING ON PLASMA CELLS

      icoadd = maxval(ixtri)

      icofc = maxval(ifctri)
      allocate (nclfc(icofc))
      do i = 1, icofc
        nullify(nclfc(i)%fp)
      end do

!  STORE FOR EACH FACE THE NUMBER OF THE TRIANGLE IT BELONGS TO

!  NCLTAL indicates for each triangle the plasma cell it belongs to
!  NCLFC is a list of triangles associated with the face IFC
      do itri = 1, ntrii
        ic = ixtri(itri)
        if (ic > 0) then
          ncltal(itri) = ic
          ifc1 = ifctri(1,itri)
          ifc2 = ifctri(2,itri)
          ifc3 = ifctri(3,itri)
          if (ifc1 > 0) call eirene_pushf(nclfc,icofc,ifc1,itri,1)
          if (ifc2 > 0) call eirene_pushf(nclfc,icofc,ifc2,itri,2)
          if (ifc3 > 0) call eirene_pushf(nclfc,icofc,ifc3,itri,3)
        else
          icoadd = icoadd + 1
          ncltal(itri) = icoadd
        end if
      end do

      ! store total number of (plasma + void triangle) cells
      ncvtot = icoadd

      ncltal(ntrii+1) = 0  ! cell ntrii+1 is used for integrals
      do ic = 1, nradd
        icoadd = icoadd + 1
        icc = ntrii+1 + ic
        ncltal(icc) = icoadd
      end do

!  INVNCLTAL indicates for each plasma cell the index of one triangle
!  belonging to this cell
      allocate (invncltal(ncvp))
      invncltal = 0

      do itri = 1,ntrii
        if (ixtri(itri) > 0) then
          if (ncltal(itri) > 0) invncltal(ncltal(itri)) = itri
        end if
      end do

C  CARRY OUT SOME CONSISTENCY CHECKS ON NEW TRIANGULAR GRID
      DO ITRI=1,NTRII
        DO IS=1,3
          IF (NCHBAR(IS,ITRI).EQ.0.AND.INMTI(IS,ITRI).EQ.0) THEN
            WRITE (iunout,*) ' ERROR IN INFCOP'
            WRITE (iunout,*) ' NCHBAR,INMTI',
     .                         NCHBAR(IS,ITRI),INMTI(IS,ITRI)
            WRITE (iunout,*) ' OPEN SIDE OF TRIANGLE ',ITRI,' SIDE ',IS
            IF (TRCGRD) THEN
              if (ixtri(itri) == 0) then
                write (iunout,*)
     .               'triangle outside original structured grid'
                call eirene_masj1('itri   ',itri)
              else
                write (iunout,*)
     .               'triangle inside original structured grid'
                call eirene_masj2('nr,np           ',
     .                           iytri(itri),ixtri(itri))
              endif
            ENDIF
            write (iunout,*) ' necke ',necke(1:3,itri)
            write (iunout,*) ' xtrian,ytrian(1) ',xtrian(necke(1,itri)),
     .                                            ytrian(necke(1,itri))
            write (iunout,*) ' xtrian,ytrian(2) ',xtrian(necke(2,itri)),
     .                                            ytrian(necke(2,itri))
            write (iunout,*) ' xtrian,ytrian(3) ',xtrian(necke(3,itri)),
     .                                            ytrian(necke(3,itri))
            IS1=IS+1
            IF (IS.EQ.3) IS1=1
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS,ITRI)),
     .                                         YTRIAN(NECKE(IS,ITRI))
            WRITE (iunout,*) ' XTRIAN,YTRIAN ',XTRIAN(NECKE(IS1,ITRI)),
     .                                         YTRIAN(NECKE(IS1,ITRI))
          ENDIF
        ENDDO
      ENDDO
C

csw 09jan2012 missing IGJUM3 correction!
CVK CORRECT IGJUM1 AND IGJUM3 FOR TRIANGLES OUTSIDE STANDARD (B2) MESH

      WRITE(iunout,*) 'IF0COP: CHECKING IGJUM3 FROM SETTINGS IN GEOUSR'
      DO ITRI=1,NTRII

C SET FLAG FOR BOUNDARY CELLS
c ifbound=t means: TRIANGULAR CELL NO. ITRI has AT LEAST ONE neighbor IT, which
c                  is not part of the STANDARD (polygonal) grid

        IFBOUND=.FALSE.
        DO IS=1,3
          IT=NCHBAR(IS,ITRI)
          IF(IT.GT.0) THEN
            IF(IXTRI(IT).LE.0) THEN
              IFBOUND=.TRUE.
              EXIT
            END IF
          END IF
        END DO

        DO I=1,NLIM
C SWITCHING ON ADDITIONAL SURFACE CHECKING FOR TRIANGULAR CELL ITRI (OUTSIDE B2 MESH)
          IF(IXTRI(ITRI).LE.0.AND.IGJUM0(I).EQ.0) THEN
            if (IGJUM3(ITRI,I).ne.0) then
              if (trcint.and.trcsur)
     .         write (iunout,*) 'activated itri,i, 1', itri,i
              IGJUM3(ITRI,I)=0
            endif
          END IF
C SWITCHING ON ADDITIONAL SURFACE CHECKING FOR BOUNDARY CELL ITRI OF B2 GRID
          IF(IFBOUND.AND.IGJUM0(I).EQ.0) THEN
            if (IGJUM3(ITRI,I).ne.0) then
              if (trcint.and.trcsur)
     .         write (iunout,*) 'activated itri,i, 2', itri,i
              IGJUM3(ITRI,I)=0
            endif
          END IF
        END DO
      END DO
CVK END
csw

      IF (NTRII.GT.0) THEN
        NLPLG=.FALSE.
        NLFEM=.TRUE.
        LEVGEO=4
        NR1ST=NTRII+1
        NLPOL=.FALSE.
        NP2ND=1
        NR1TAL=NTRII
        NP2TAL=NP2ND
        NT3TAL=NT3RD
        NSBOX_TAL=NR1TAL*NP2TAL*NT3TAL*NBMLT+NRADD
        NGITT = COUNT(INMTI(1:3,1:NTRII) .NE. 0)
      ELSE IF (.NOT.NLTRIMESH) THEN
        CALL EIRENE_LEER(2)
        WRITE (iunout,*) 'CASE STOPPED !'
        WRITE (iunout,*) 'NO TRIANGULAR MESH INFORMATION FOUND'
        WRITE (iunout,*) 'PLEASE PROVIDE '//FORT//'3[3-5] FILES'
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

#ifdef WG_TODO
      ALLOCATE(XNMTI(NNNMTI),YNMTI(NNNMTI),TNMTI(NNNMTI))
      XNMTI(1:NNNMTI)=XNMTI_BUF(1:NNNMTI)
      YNMTI(1:NNNMTI)=YNMTI_BUF(1:NNNMTI)
      TNMTI(1:NNNMTI)=TNMTI_BUF(1:NNNMTI)
      DEALLOCATE(XNMTI_BUF,YNMTI_BUF,TNMTI_BUF)
#endif

      CALL IMPORT_B2_CHEM_SPUT

      CALL EIRENE_LEER(2)
      CALL EIRENE_HEADNG(' CASE REDEFINED IN COUPLE_SOLPS-ITER:',37)
      WRITE (iunout,*) 'NLPLG,NLFEM ',NLPLG,NLFEM
      WRITE (iunout,*) 'NLPOL       ',NLPOL

      WRITE (IUNOUT,*) 'NEW (FINE, UNSTRUCTURED) GRID'
      CALL EIRENE_MASJ5('NR1ST,   NP2ND,  NSURF,  NRADD,  NSBOX  ',
     .                   NR1ST,   NP2ND,  NSURF,  NRADD,  NSBOX)

      WRITE (IUNOUT,*) 'NR1TAL,NP2TAL,NSBOX_TAL',
     .                  NR1TAL,NP2TAL,NSBOX_TAL
      CALL EIRENE_LEER(2)
CTRIG E

      IUNIN = IUNIN_SAVE
      IF (IUSROUT /= 0) CLOSE(IUSROUT)
csw
      RETURN
 1983 WRITE(IUNOUT,*)
     w     "ERROR IN INFCOP: CANNOT OPEN "//FORT//"33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      CALL EIRENE_EXIT_OWN (1)
      RETURN
 1984 WRITE(IUNOUT,*)
     w     "ERROR IN INFCOP: CANNOT OPEN "//FORT//"34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      CALL EIRENE_EXIT_OWN (1)
      RETURN
 1985 WRITE(IUNOUT,*) "ERROR IN INFCOP: CANNOT OPEN "//FORT//"35 ",
     .                "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      CALL EIRENE_EXIT_OWN (1)
      RETURN

      CONTAINS


       SUBROUTINE EIRENE_READ14_FIXED
       USE EIRMOD_CCOUPL
       USE EIRMOD_COMUSR
       USE EIRMOD_CTEXT
       IMPLICIT NONE
       INTEGER :: J, IAIN, IAOT, IPL, IT
       INTEGER :: ICOND
       LOGICAL :: LDUMMY
       EXTERNAL :: EIRENE_MASPRM, EIRENE_MASR4,
     .             EIRENE_EXIT_OWN, EIRENE_SKIP_READ_COMMENT,
     .             FIX_INTEGER_INPUT, FIX_LOGICAL_INPUT

       WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED'
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
       CALL EIRENE_LEER(1)
       CALL EIRENE_ALLOC_CCOUPL(1)
       READ (IUNIN,'(A72)') ZEILE
       call fix_logical_input(zeile,5)
       READ (ZEILE,'(5L1)') LSYMET,LBALAN,LCOARSE,LDUMMY,LSPRCL
       IF (.NOT.NLTRIMESH) NLTRIMESH=LDUMMY
       IF (TRCINT) THEN
        WRITE (iunout,*) ' LSYMET,LBALAN,LCOARSE,NLTRIMESH,LSPRCL = ',
     .                     LSYMET,LBALAN,LCOARSE,NLTRIMESH,LSPRCL
       END IF

       READ (IUNIN,'(A72)') ZEILE
       call fix_integer_input(zeile,11)
       READ (ZEILE,'(11I6)') NFLA,
     .                      ntrfrm,nfull,ibrad,ibpol,ibtor,
     .                      dum7,dum8,dum9,dum10,NCOND
       NPLS_FIX = NFLA

cdr added in april 2015:
c  flags for orientation of radial, poloidal and toroidal magnetic field components
c  currently: only btor is used, and only in the vacuum region outside the B2.5 grid
cdr everywhere else we now have the full BRADB,BPOLB and BTORB tallies (and the pitch angle
cdr tally RRB should be redundant?
       brad = 1._dp
       if (ibrad < 0) brad = -brad
       bpol = 1._dp
       if (ibpol < 0) bpol = -bpol
       btor = 1._dp
       if (ibtor < 0) btor = -btor

       IF (TRCINT) THEN
         WRITE (iunout,*) ' NFLA = ', NFLA
         WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
       ENDIF
       NFLB=0
       NATMA=0
       DO IPL=1,NPLSI
         READ (IUNIN,'(2I6,2E12.4,I6)')
     -     I,IFLB(IPL),FCTE(IPL),BMASS(IPL),L
         if(iflb(ipl) .gt. 0) lkindp(ipl) = l
         IF (TRCINT) WRITE (iunout,'(2I6,1P,2ES12.4,I6)')
     .     IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL),L
         if(iflb(ipl) .gt. 0) then
           NFLB=NFLB+1
           if(lkindp(ipl).le.0) then
             k=NINT(BMASS(IPL))
             l=0
             do j=1,natmi
               if(k.eq.nmassa(j)) then
                 if(l.eq.0) then
                   l=j
                 else if(l.gt.0) then
                   l=-1
                 end if
               end if
             end do
             if(l.gt.0) then
               lkindp(ipl)=l
             else
               write(iunout,*) 'Plasma species ',ipl,
     -                         ': composition unknown'
             end if
           end if
           found = .false.
           do j=1,ipl-1
             found=found.or.lkindp(j).eq.lkindp(ipl)
           end do
           if (.not.found) NATMA=NATMA+1
         endif
       END DO  ! IPL

       IF (.NOT.ALLOCATED(CONDENSE_FROM)) THEN
         ALLOCATE (CONDENSE_FROM(NCOND))
         ALLOCATE (CONDENSE_TO(NCOND))
       ENDIF

       IF (NCOND.GT.0) THEN
         DO ICOND =1,NCOND
           READ (IUNIN,'(2I6)') CONDENSE_FROM(ICOND),CONDENSE_TO(ICOND)
         ENDDO
       ENDIF

       READ (IUNIN,'(3I6)') NCVP, NCIP, NFCP
       IF (TRCINT) THEN
         WRITE (iunout,*) '  NCVP = ', NCVP
         WRITE (iunout,*) '  NCIP = ', NCIP
         WRITE (iunout,*) '  NFCP = ', NFCP
       ENDIF
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
       READ (IUNIN,'(I6)') NTARGI
       IF (TRCINT) WRITE (iunout,*) 'NTARGI=    ',NTARGI
       CALL EIRENE_LEER(1)
       IF (NTARGI.GT.NSTEP) THEN
         CALL EIRENE_MASPRM ('NSTEP',5,NSTEP,'NTARGI',6,NTARGI,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF
C  NUMBER OF PARTS PER TARGET SOURCE
C  NO LONGER NEEDED, INPUT KEPT FOR THE CONVERSION ROUTINE
       IF (NTARGI.GT.0) THEN
         ALLOCATE (JDUMMY(NTARGI))
         JDUMMY = 0
         READ (IUNIN,'(12I6)') (JDUMMY(IT),IT=1,NTARGI)
         DEALLOCATE (JDUMMY)
       ENDIF
C  ALL INDICES: AFTER INDEX MAPPING
C  NDT: INDEX OF X-CELL (EAST OR NORTH SURFACE OF BRAAMS CELL) OF TARGET
C  NINCT: DIRECTION OF OUTER TARGET NORMAL WITH RESPECT TO POSITIVE DIR.
C  NIXY: SOURCE ON Y SURFACE: NIXY=1; SOURCE ON X SURFACE: NIXY=2
C  NTIN,NTEN: SOURCE RANGE FROM GRID POINT NTIN TO GRID POINT NTEN
C  ITYPE: TYPE OF THE STRATUM (0: BULK IONS, 1: ATOMS, 2: MOLECULES,
C        3: TRACE IONS, 4: BULK IONS, 5: PHOTONS)
       IF (TRCINT)
     .  WRITE (iunout,*) '    IT,  NDT, NPTC, NPTCM',
     .              ',NSPZI,NSPZE,NEMOD,ITYPE'
       DO IT=1,NTARGI
  331    READ (IUNIN,'(A72)') ZEILE
         IF (ZEILE(1:1).EQ.'*') THEN
C          WRITE (iunout,......)
           GOTO 331
         ENDIF
         READ (ZEILE,'(12I6)') I,NDT(IT),
     .                           NPTC(IT),NPTCM(IT),
     .                           NSPZI(IT),NSPZE(IT),
     .                           NEMOD(IT), !VK NPTCM
     .                           ITYPE(IT)  !NH
         NSPZI(IT)=MAX0(1,NSPZI(IT))
         NSPZE(IT)=MIN0(NFLA,NSPZE(IT))
         IF (NSPZE(IT).LT.NSPZI(IT)) THEN
           WRITE (iunout,*) 'WARNING FROM INFCOP: '
           WRITE (iunout,*) 'ITARG : ',IT
           WRITE (iunout,*) 'NSPZI,NSPZE MODIFIED TO 1,NFLA, RESP.'
           NSPZI(IT)=1
           NSPZE(IT)=NFLA
         ENDIF
         IF (TRCINT)
     .     WRITE (iunout,'(1X,12I6)') IT,NDT(IT),
     .                                   NPTC(IT),NPTCM(IT),
     .                                   NSPZI(IT),NSPZE(IT),
     .                                   NEMOD(IT), !VK NPTCM
     .                                   ITYPE(IT)  !NH
       END DO  ! IT
       IF (TRCINT) CALL EIRENE_LEER(1)
       READ (IUNIN,'(6E12.4)')  CHGP,CHGEE,CHGEI,CHGMOM
       IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2.5 INTO EIRENE
C  HERE: B2.5 VOLUME TALLIES
       READ (IUNIN,'(I6)') NAINB
C  ADDITIONAL INPUT TALLY ADIN:  ITAL=12
       NAIN = MAX(NAIN,NAINB)
       CALL EIRENE_ALLOC_CCOUPL(2)
       WRITE (iunout,*) '       NAINI = ',NAINB
       IF (NAINB.GT.NAIN) THEN
         CALL EIRENE_MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF
       IF (TRCINT.AND.NAINB.GT.0)
     .      WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
       DO IAIN=1,NAINB
         READ (IUNIN,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
         READ (IUNIN,'(A72)') TXTPLS(IAIN,12)
         READ (IUNIN,'(2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
         IF (TRCINT) THEN
           WRITE (iunout,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
           WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
           WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
         ENDIF
       END DO  ! IAIN
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
       READ (IUNIN,'(I6)') NAOTB
       WRITE (iunout,*) '       NAOTI = ',NAOTB
       IF (NAOTB.GT.NLIMPS) THEN
         CALL EIRENE_MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF
       IF (TRCINT.AND.NAOTB.GT.0)
     .      WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
       DO IAOT=1,NAOTB
         READ (IUNIN,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
         IF (TRCINT) THEN
           WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
         ENDIF
       END DO  ! IAOT
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
       RETURN
       END SUBROUTINE EIRENE_READ14_FIXED

#ifndef NO_JSON
       SUBROUTINE EIRENE_READ14_JSON(json,me)
       USE EIRMOD_JSON
       use json_module           !IGNORE
!pgf     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

       IMPLICIT NONE

       class(json_core), intent(inout) :: json
       type(json_value), pointer, intent(in) :: me
       type(json_value), pointer :: pflds, pfld, ptrgs, prt,
     .                              pconds, pcond, padds, padd
!pgf       character(kind=CK,len=:), allocatable :: txt
       character(kind=json_CK,len=:), allocatable :: txt
       integer :: j, iain, iaot, ipl, npl, ntrg, icond, nconds
       integer, allocatable :: ihelp(:)
       logical :: found, foundi, foundo
       external :: eirene_masprm, eirene_masr4

       WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED'
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
       CALL EIRENE_LEER(1)
       CALL EIRENE_ALLOC_CCOUPL(1)

       call json%get(me,'LSYMET',lsymet,found)
       call json%get(me,'LBALAN',lbalan,found)
       call json%get(me,'LCOARSE',lcoarse,found)
       call json%get(me,'LSPRCL',lsprcl,found)
       IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN,LCOARSE,LSPRCL = ',
     .                     LSYMET,LBALAN,LCOARSE,LSPRCL

       call json%get(me,'NFLA',nfla,found)
       call json%get(me,'NCOND',ncond,found)
       call json%get(me,'NTRFRM',ntrfrm,found)
       call json%get(me,'NFULL',nfull,found)
       call json%get(me,'IBRAD',ibrad,found)
       call json%get(me,'IBPOL',ibpol,found)
       call json%get(me,'IBTOR',ibtor,found)
       NPLS_FIX = NFLA

cdr added in april 2015:
c  flags for orientation of radial (not in use), poloidal and toroidal magnetic field components
       brad = 1._dp
       if (ibrad < 0) brad = -brad
       bpol = 1._dp
       if (ibpol < 0) bpol = -bpol
       btor = 1._dp
       if (ibtor < 0) btor = -btor

       IF (TRCINT) THEN
         WRITE (iunout,*) ' NFLA = ', NFLA
         WRITE (iunout,*) ' IPLS,IFLB(IPLS),FCTE(IPLS),BMASS(IPLS)'
       ENDIF

       call json%get_child(me,'B2FLUIDS',pflds)
       call json%info(pflds,n_children=npl)
       if (npl /= NPLSI) then
         write (iunout,*) ' NUMBER OF BULK IONS DOES',
     .          ' NOT MATCH NUMBER OF FLUIDS FOUND IN FILE '
         write (iunout,*) 'NPLSI = ',nplsi
         write (iunout,*) 'NPL =   ',npl
         call eirene_exit_own(1)
       end if

       NFLB=0
       NATMA=0
       DO IPL=1,NPLSI
         call json%get_child(pflds,ipl,pfld)

         call json%get(pfld,'IPL',i,found)
         call json%get(pfld,'IFLB',iflb(ipl),found)
         call json%get(pfld,'FCTE',fcte(ipl),found)
         call json%get(pfld,'BMASS',bmass(ipl),found)
         call json%get(pfld,'LKINDP',l,found)
         if (found) then
           if(iflb(ipl).gt.0) then
             lkindp(ipl)=l
             NFLB=NFLB+1
             if(lkindp(ipl).le.0) then
               k=NINT(BMASS(IPL))
               l=0
               do j=1,natmi
                 if(k.eq.nmassa(j)) then
                   if(l.eq.0) then
                     l=j
                   else if(l.gt.0) then
                     l=-1
                   end if
                 end if
               end do
               if(l.gt.0) then
                 lkindp(ipl)=l
               else
                 write(iunout,*) 'Plasma species ',ipl,
     -                           ': composition unknown'
               end if
             end if
             found = .false.
             do j=1,ipl-1
               found=found.or.lkindp(j).eq.lkindp(ipl)
             end do
             if (.not.found) NATMA=NATMA+1
           endif
         end if
         nullify(pfld)
         IF (TRCINT) WRITE (iunout,'(2I6,2ES12.4,I6)')
     .               IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL),L
       END DO
       nullify(pflds)

       IF (NCOND.GT.0) THEN
         IF (.NOT.ALLOCATED(CONDENSE_FROM)) THEN
           ALLOCATE (CONDENSE_FROM(NCOND))
           ALLOCATE (CONDENSE_TO(NCOND))
         ENDIF
         call json%get_child(me,'CONDENSE',pconds)
         call json%info(pconds,n_children=nconds)
         if (nconds /= NCOND) then
           write (iunout,*) ' NUMBER OF CONDENSATION SPECIES DOES',
     .          ' NOT MATCH NUMBER FOUND IN FILE '
           write (iunout,*) 'NCOND  = ',ncond
           write (iunout,*) 'NCONDS = ',nconds
           call eirene_exit_own(1)
         end if
         DO ICOND =1,NCOND
           call json%get_child(pconds,icond,pcond)
           call json%get(pcond,'CONDENSE_FROM',condense_from(icond))
           call json%get(pcond,'CONDENSE_TO',condense_to(icond))
         ENDDO
       ENDIF

       call json%get(me,'NCVP',ncvp,found)
       call json%get(me,'NCIP',ncip,found)
       call json%get(me,'NFCP',nfcp,found)
       IF (TRCINT) THEN
         WRITE (iunout,*) '  NCVP = ', NCVP
         WRITE (iunout,*) '  NCIP = ', NCIP
         WRITE (iunout,*) '  NFCP = ', NFCP
       ENDIF

C     NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
       call json%get(me,'NTARGI',ntargi,found)
       IF (TRCINT) WRITE (iunout,*) 'NTARGI=    ',NTARGI
       CALL EIRENE_LEER(1)
       IF (NTARGI.GT.NSTEP) THEN
         CALL EIRENE_MASPRM ('NSTEP',5,NSTEP,'NTARGI',6,NTARGI,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF

C     NUMBER OF PARTS PER TARGET SOURCE
C     NO LONGER NEEDED, INPUT KEPT FOR THE CONVERSION ROUTINE
       IF (NTARGI.GT.0) THEN
         call json%get(me,'NTGPRT',ihelp,found)
         if (found) deallocate(ihelp)

C  ALL INDICES: AFTER INDEX MAPPING
C  NDT: INDEX OF X-CELL (EAST OR NORTH SURFACE OF BRAAMS CELL) OF TARGET
C  NINCT: DIRECTION OF OUTER TARGET NORMAL WITH RESPECT TO POSITIVE DIR.
C  NIXY: SOURCE ON Y SURFACE: NIXY=1; SOURCE ON X SURFACE: NIXY=2
C  NTIN,NTEN: SOURCE RANGE FROM GRID POINT NTIN TO GRID POINT NTEN
         IF (TRCINT)
     .     WRITE (iunout,*) '    IT,  NDT, NPTC, NPTCM',
     .                       ',NSPZI,NSPZE,NEMOD,ITYPE'

         call json%get_child(me,'TARGETS',ptrgs,found)
         call json%info(ptrgs,n_children=ntrg)
         if (ntrg /= NTARGI) then
           write (iunout,*) ' NUMBER OF TARGETS DOES',
     .          ' NOT MATCH NUMBER OF TARGETS FOUND IN FILE '
           write (iunout,*) 'NTARGI = ',ntargi
           write (iunout,*) 'NTRG =   ',ntrg
           call eirene_exit_own(1)
         end if
         DO IT=1,NTARGI
           call json%get_child(ptrgs,it,prt,found)
           call json%get(prt,'NDT',ndt(it),found)
           call json%get(prt,'NPTC',nptc(it),found)
           call json%get(prt,'NPTCM',nptcm(it),found)
           call json%get(prt,'NSPZI',nspzi(it),found)
           call json%get(prt,'NSPZE',nspze(it),found)
           call json%get(prt,'NEMOD',nemod(it),found)
           call json%get(prt,'ITYPE',itype(it),found)
           nullify(prt)
           NSPZI(IT)=MAX0(1,NSPZI(IT))
           NSPZE(IT)=MIN0(NFLA,NSPZE(IT))
           IF (NSPZE(IT).LT.NSPZI(IT)) THEN
             WRITE (iunout,*) 'WARNING FROM INFCOP: '
             WRITE (iunout,*) 'ITARG : ',IT
             WRITE (iunout,*) 'NSPZI,NSPZE MODIFIED TO 1,NFLA, RESP.'
             NSPZI(IT)=1
             NSPZE(IT)=NFLA
           ENDIF
           IF (TRCINT)
     .      WRITE (iunout,'(1X,12I6)') IT,NDT(IT),
     .                                    NPTC(IT),NPTCM(IT),
     .                                    NSPZI(IT),NSPZE(IT),
     .                                    NEMOD(IT), !VK NPTCM
     .                                    ITYPE(IT)  !NH
         END DO
         IF (TRCINT) CALL EIRENE_LEER(1)
         nullify(ptrgs)
       END IF

       call json%get(me,'CHGP',chgp,found)
       call json%get(me,'CHGEE',chgee,found)
       call json%get(me,'CHGEI',chgei,found)
       call json%get(me,'CHGMOM',chgmom,found)
       IF (TRCINT) CALL EIRENE_MASR4
     .                         ('CHGP,CHGEE,CHGEI,CHGMOM         ',
     .                           CHGP,CHGEE,CHGEI,CHGMOM)

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
       call json%get(me,'NAINB',nainb,found)
C  ADDITIONAL INPUT TALLY ADIN:  ITAL=12
       NAIN = MAX(NAIN,NAINB)
       CALL EIRENE_ALLOC_CCOUPL(2)
       WRITE (iunout,*) '       NAINI = ',NAINB
       IF (NAINB.GT.NAIN) THEN
         CALL EIRENE_MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF
       IF (NAINB > 0) THEN
         IF (TRCINT)
     .     WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
         call json%get(me,'ADD_IN_TAL',padds,foundi)
         IF (FOUNDI) THEN
           DO IAIN=1,NAINB
             call json%get_child(padds,iain,padd,found)
             call json%get(padd,'NAINS',nains(iain),found)
             call json%get(padd,'NAINT',naint(iain),found)
             call json%get(padd,'TXTPLS',txt,found)
             txtpls(iain,12) = txt
             deallocate(txt)
             call json%get(padd,'TXTPSP',txt,found)
             txtpsp(iain,12) = txt
             deallocate(txt)
             call json%get(padd,'TXTPUN',txt,found)
             txtpun(iain,12) = txt
             deallocate(txt)
             nullify(padd)
             IF (TRCINT) THEN
               WRITE(iunout,'(6I6)') IAIN,NAINS(IAIN),NAINT(IAIN)
               WRITE(iunout,'(1X,A72)') TXTPLS(IAIN,12)
               WRITE(iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
             ENDIF
           END DO
           nullify(padds)
         END IF
       END IF

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
       call json%get(me,'NAOTB',naotb,found)
       WRITE (iunout,*) '       NAOTI = ',NAOTB
       IF (NAOTB.GT.NLIMPS) THEN
         CALL EIRENE_MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
         WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
         CALL EIRENE_EXIT_OWN(1)
       ENDIF
       IF (NAOTB > 0) THEN
         IF (TRCINT)
     .     WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
         call json%get(me,'ADD_OUT_TAL',padds,foundo)
         IF (FOUNDO) THEN
           DO IAOT=1,NAOTB
             call json%get_child(padds,iaot,padd,found)
             call json%get(padd,'NAOTS',naots(iaot),found)
             call json%get(padd,'NAOTT',naott(iaot),found)
             nullify(padd)
             IF (TRCINT) THEN
               WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
             ENDIF
           END DO
           nullify(padds)
         END IF
       END IF
C
C  INPUT BLOCK 14 DONE
C
       IUSROUT = -9999
       CALL EIRENE_OPENFILE(IUSROUT,FILE='user_data.input',STATUS='OLD',
     .      FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IO)
       IUNIN_SAVE = IUNIN
       IF (IO == 0) THEN
         IUNIN = IUSROUT
         CALL EIRENE_LEER(1)
         WRITE (IUNOUT,*) 'USER SPECIFIC INPUT READ FROM ',
     .        'user_data.input'
       ELSE
         CALL EIRENE_LEER(1)
         WRITE (IUNOUT,*) 'NO FILE FOR USR SPECIFIC INPUT FOUND'
       END IF

       RETURN
       END SUBROUTINE EIRENE_READ14_JSON
#endif

      END SUBROUTINE EIRENE_IF0COP
C
C   GEOMETRY DEFINITION PART FINISHED
C
      SUBROUTINE EIRENE_IF1COP(IENTRY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IENTRY
      INTEGER :: IN, IT, IFL, IPLSV, IPLSTI, JPLS,
     I           IAIN, IREAD, ISTRAI, ITRI, IC
      REAL(DP) :: BX, BY, BZ, BN
      EXTERNAL :: EIRENE_PLASM, EIRENE_LEER,
     .            EIRENE_EXIT_OWN, EIRENE_READ_TARGET_DATA

C
C   NOW READ THE PLASMA STATE GIVEN BY BRAAMS
C   AT PRESENT THE DATA COME FROM THE FILE FT31
C   FURTHERMORE: SCALING TO EIRENE UNITS AND INDEX MAPPING
C   STATEMENT NO. 2000 ---> 2999
C
C  IN CASE OF "SHORT CYCLE" THE PLASMA STATE IS TRANSFERRED VIA COMMON BRAEIR
C
      IF (IENTRY.EQ.0) THEN
        LSHORT = .FALSE.
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'IF1COP CALLED'
C  SKIP READING PLASMA, IF NLPLAS
        IF (NLPLAS) THEN
          WRITE (iunout,*) 'PLASMA DATA EXPECTED ON BRAEIR'
        ELSE
C
C  TRANSFER PLASMA PROFILES VIA FORT 31
C
          IF (.NOT.(INDPRO(1).EQ.6.OR.INDPRO(2).EQ.6.OR.INDPRO(3).EQ.6
     .          .OR.INDPRO(4).EQ.6.OR.INDPRO(5).EQ.6)) RETURN
C
          WRITE (iunout,*) 'PLASMA DATA EXPECTED ON '//FORT//'31'
C
          OPEN (UNIT=31,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1986)
          REWIND 31
C
          IF (NFLA.GT.NFL) THEN
            WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFCOP'
            WRITE (iunout,*) ' NFLA MUST BE <= NFL'
            WRITE (iunout,*) ' NFLA,NFL = ',NFLA,NFL
            CALL EIRENE_EXIT_OWN(1)
          ENDIF

!         CALL EIRENE_ALLOC_BRAEIR(NCVP,NFCP,NFL)
C
C  B2-BRAAMS CODE SPECIFIC BEGIN

          CALL EIRENE_PLASM (31, IUNOUT, NCVP, NCIP, NFCP, NFLA)

          CALL EIRENE_READ_TARGET_DATA (31, IUNOUT, NTARGI, NFLA)
        END IF ! NLPLAS
C
C  IN CASE OF "SHORT CYCLE" OR TIME DEP. MODE
C  THE PLASMA STATE IS TRANSFERRED VIA COMMON
C  ONLY SCALING TO EIRENE UNITS AND INDEX MAPPING NEEDS TO BE DONE HERE
C
      END IF

!pb  for the time being set default values
csw 26sep2011
c      delta_sheathxb=3.1
c      delta_sheathyb=3.1
csw 06jan2012 NO! not for direct comparison SOLSP4.3!!!

cwd do not set zero here; are initialized to zero, but could be passed through braeir for
c   consistent potential drop on both sides (including effect of currents, secondary electron
c   emission)
c      delta_sheathxb=0.0
c      delta_sheathyb=0.0
csw
C
C  RESET 2D ARRAYS ONTO 1D EIRENE ARRAYS, RESCALE TO EIRENE UNITS
C  AND CONVERT BRAAMS VECTORS INTO CARTESIAN EIRENE VECTORS
C
C  UNITS CONVERSION FACTORS
      V=1.D2                             !pb 1.e2 -> 1.d2
      VL=1.D6                            !pb 1.e6 -> 1.d6
CTRIG A
C  VACUUM DATA NEEDED FOR REGION OUTSIDE B2-MESH
      TVAC=0.02_DP                       !pb 0.02 -> 0.02_dp
      DVAC=1.D2
      VVAC=0.
      ZVAC=0.
CTRIG E
      DO 2105 JPLS=1,NPLSI
        D(JPLS)=1.D-6*FCTE(JPLS)         !pb 1.e-6 -> 1.d-6
        FL(JPLS)=ELCHA*FCTE(JPLS)        !SCALING FOR FLUX: 1/S --> AMP
 2105 CONTINUE
C
C  SET PLASMA BACKGROUND ON TRIANGULAR GRID
C  A TRIANGULAR CELL ITRI RECEIVES THE PLASMA DATA FROM ITS LARGER HOST CELL (IX,IY),
C  WITHOUT ANY WEIGHTING/INTERPOLATION ETC...
C
      BZINTF = 1._DP
      BFINTF = 1._DP
      DO ITRI=1,NTRII
        IC=IXTRI(ITRI)
        IF (IC .GT. 0) THEN
          TEINTF(ITRI)=TEB(IC)
C
C  ION OR NEUTRAL TEMPERATURE AVAILABLE FROM PLASMA FLUID CODE
          DO JPLS=1,NPLSI
            IF (NCHRGP(MPLSTI(JPLS)).EQ.0) THEN !! neutral species
              TIINTF(MPLSTI(JPLS),ITRI)=TNB(IC)
            ELSE                       !! ion species
              TIINTF(MPLSTI(JPLS),ITRI)=TIB(IC)
            ENDIF
          END DO
C
CDR,WD Construct magnetic field from:
C  poloidal field
          BX=BBXB(IC)
          BY=BBYB(IC)   ! + PVY(IN)*BRADB(IN)
c  toroidal field
          BZ=BBZB(IC)
c  normalize B field vector to length 1 (one)
c  avoid using bfeldb, because limited accuracy if coming through fort.31 (standalone cases)
          BN=SQRT(BX**2 + BY**2 + BZ**2)
          IF (BN <= EPS30) THEN
            BZ = 1._DP
            BN = 1._DP
          END IF
          BXINTF(ITRI)=BX/BN
          BYINTF(ITRI)=BY/BN
          BZINTF(ITRI)=BZ/BN
          BFINTF(ITRI)=BN

cdr  try to also bring grad(PSI) direction unit vector into eirene
          BXPERF(ITRI)=PVX(IC)
          BYPERF(ITRI)=PVY(IC)
c
          VLINTF(ITRI)=VOLB(IC)*VL
        ELSE
c  outside the original b2.5 grid:
C    set default vacuum temperatures TVAC
C    set default B field: (0,0,1) or (0,0,-1)
          TEINTF(ITRI)=TVAC
          TIINTF(1,ITRI)=TVAC
c  default unit vector along B
          BXINTF(ITRI)=0.
          BYINTF(ITRI)=0.
cdr  BTOR is a sign factor (=1 or = -1), read from block 14.
cdr  This may perhaps be inconsistent with the use of BTORB above
          BZINTF(ITRI)=1.*BTOR  ! B_INTF is interpreted as
                                ! unit vector inside eirene
c  default modulus of B field
          BFINTF(ITRI)=1.  !dr here:  = BTORB (magnitude),
                           !   e.g. for Zeeman splitting
c
C         VLINTF(ITRI)=1.
        ENDIF
      ENDDO  !ITRI LOOP

C  ntrii+1: for average over ntrii grid

C  additional cells from input block 2e
      DO itri=ntrii+2, nsbox
        TEINTF(ITRI)=TVAC
        TIINTF(1,ITRI)=TVAC
c  default unit vector along B
        BXINTF(ITRI)=0.
        BYINTF(ITRI)=0.
cdr  apply same sign factor BTOR as above
        BZINTF(ITRI)=1.*BTOR
c  default modulus of B field
        BFINTF(ITRI)=1.  !dr here:   = BTORB (magnitude),
                         !   e.g. for Zeeman splitting
c
C       VLINTF(ITRI)=1.
      ENDDO
C
C  SET SAME ION TEMPERATURE FOR ALL EIRENE BACKGROUND SPECIES
C

      DO IPLSTI=1,NPLSTI
       DO ITRI=1,NSBOX
        TIINTF(IPLSTI,ITRI)=TIINTF(1,ITRI)
       END DO
      END DO
C
CDR  set density from B2 array DNIB, for each fluid
CDR  set plasma flow velocity field from B2 arrays UPB (parallel velocity)
c  without drifts:
c  upb * pitch: poloidal velocity (i.e. in cartesian x,y direction).
c  poloidal field direction is given by that of the poloidal cell face PU..(in),
C  i.e. along a flux surface. (PU(...) is cell-centered)
c  and upb*sqrt(1-pitch^2): toroidal velocity  (i.e. cartesian z direction (nltrz) or
c                                                toroidal phi direction (nltra)
c  sign of flow field follows the sign of poloidal grid in B2.
c
c  with drifts:
c  uudia and vvdia are additional flow velocities in B2.5 only.
c
c
c
c
      IREAD=0
      DO 2200 JPLS=1,NPLSI
        IF (IFLB(JPLS).GT.0) THEN  ! DEAL WITH B2.5 ION SPECIES ONLY,
                                   ! EXCLUDE VIRTUAL EIRENE BACKGROUND
          IPLSV=MPLSV(JPLS)
          DO 2201 IFL=1,NFLA
            IF (IFLB(JPLS).NE.IFL) GOTO 2201
            DO ITRI=1,NTRII
              IC=IXTRI(ITRI)
              IF (IC .GT. 0) THEN

C  SET PLASMA DENSITY, ION AVERAGE CHARGE, AND PLASMA VELOCITIES
                DIINTF(JPLS,ITRI)=DNIB(IC,IFL)*D(JPLS)
                ZIIN(JPLS,ITRI)=ZIB(IC,IFL)
cwd projection onto poloidal, radial, toroidal directions now done on the B2.5 side,
c   including correct signs, and effect of drifts
                VXINTF(IPLSV,ITRI)=VXB(IC,IFL)*V
                VYINTF(IPLSV,ITRI)=VYB(IC,IFL)*V
                VZINTF(IPLSV,ITRI)=VZB(IC,IFL)*V
              ELSE
c  region outside B2.5 grid, but inside triangular grid
                DIINTF(JPLS,ITRI)=DVAC
                VXINTF(IPLSV,ITRI)=VVAC
                VYINTF(IPLSV,ITRI)=VVAC
                VZINTF(IPLSV,ITRI)=VVAC
                ZIIN(JPLS,ITRI)=ZVAC
              ENDIF
            ENDDO  !itri loop
c  further additional cells from input block 2e ?
            DO itri=ntrii+2, nsbox
              DIINTF(JPLS,ITRI)=DVAC
              VXINTF(IPLSV,ITRI)=VVAC
              VYINTF(IPLSV,ITRI)=VVAC
              VZINTF(IPLSV,ITRI)=VVAC
              ZIIN(JPLS,ITRI)=ZVAC
            ENDDO

C  EIRENE BACKGROUND SPECIES "JPLS" IS NOW FILLED WITH B2.5 DATA "IFL"
 2201     CONTINUE   !IFL loop

C  NO DATA FOR "JPLS" IN B2 FILES
        ELSEIF (IFLB(JPLS).EQ.-13) THEN
C  READ DATA FOR "JPLS" FROM EIRENE DUMP FILE FT13
csw 09jan2012 NO!!!
csw       IF (IREAD.EQ.0) THEN
csw         OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
csw         REWIND 13
csw         READ (13,IOSTAT=IO) TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN
csw         IREAD=1
csw         IF (TRCFLE) WRITE (iunout,*) 'READ 13: RCMUSR, IO= ',IO
csw         CLOSE (UNIT=13)
csw       ENDIF
csw       IF (IO.EQ.0) THEN
cdr         IF (TRCINT) THEN
cdr            WRITE(IUNOUT,*)
cdr              'PLASMA DATA FOR JPLS READ FROM '//FORT//'13'
cdr            WRITE(IUNOUT,*) 'JPLS,IPLSV,IPLSTI ',JPLS,IPLSV,IPLSTI
cdr         ENDIF

            IPLSTI = MPLSTI(JPLS)
            IPLSV = MPLSV(JPLS)
            DO ITRI=1,NTRII
              DIINTF(JPLS,ITRI)=DIIN(JPLS,ITRI)
              VXINTF(IPLSV,ITRI)=VXIN(IPLSV,ITRI)
              VYINTF(IPLSV,ITRI)=VYIN(IPLSV,ITRI)
              VZINTF(IPLSV,ITRI)=VZIN(IPLSV,ITRI)
              TIINTF(IPLSTI,ITRI)=TIIN(IPLSTI,ITRI)
              ZIIN(JPLS,ITRI)=0.0
            ENDDO
csw       ENDIF
        ELSE
C  SET PARAMETERS FOR SPECIES JPLS TO ZERO
C  NOTHING TO BE DONE HERE
        ENDIF
 2200 CONTINUE
C  B2.5-BRAAMS CODE-SPECIFIC END
C
C
C  READ OTHER B2.5 ARRAYS INTO EIRENE, FOR PRINTOUT AND PLOTTING
C
c  density, species index as in B2 code, cell-centered
      DO 2300 IAIN=1,NAINB
        IF (NAINT(IAIN).EQ.1.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2321 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=DNIB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2321     CONTINUE
#ifdef WG_TODO
c  poloidal (projection) flow velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.2.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2322 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UUB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2322     CONTINUE
c  radial drift velocity, species index as in B2 code, east surface-centered
        ELSEIF (NAINT(IAIN).EQ.3.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2323 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VVB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2323     CONTINUE
c  plasma pressure, cell-centered, no species index
        ELSEIF (NAINT(IAIN).EQ.6) THEN
          DO 2326 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=PRB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2326     CONTINUE
c  parallel velocity, species index as in B2 code, north surface-centered
        ELSEIF (NAINT(IAIN).EQ.7.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2327 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UPB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2327     CONTINUE
c  pitch angle, no species index
        ELSEIF (NAINT(IAIN).EQ.8) THEN
          DO 2328 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=RRB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2328     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.9.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO 2329 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FNIB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2329     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.11) THEN
          DO 2331 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEIB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2331     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.13) THEN
          DO 2333 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=FEEB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2333     CONTINUE
        ELSEIF (NAINT(IAIN).EQ.15.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=UUDIAB(IXTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
          ENDDO
        ELSEIF (NAINT(IAIN).EQ.16.AND.NAINS(IAIN).GT.0.AND.
     .      NAINS(IAIN).LE.NFLA) THEN
          DO IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VVDIAB(IXTRI(IN),IYTRI(IN),NAINS(IAIN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
          ENDDO
c   cell volume as in b2.5 code, no species index (cell-centered)
        ELSEIF (NAINT(IAIN).EQ.17) THEN
          DO 2335 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=VOLB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2335     CONTINUE
#endif
cdr  magnetic field strength, Tesla
        ELSEIF (NAINT(IAIN).EQ.18) THEN
          DO 2336 IN=1,NTRII
            IF (IXTRI(IN).GT.0) THEN
              ADINTF(IAIN,IN)=BFELDB(IXTRI(IN))
            ELSE
              ADINTF(IAIN,IN)=0.
            ENDIF
 2336     CONTINUE

cdr  free: NAINT=20 --39: reserved for AMDIAG: scaled atomic/molecular rate coefficients
cdr                       evaluated on computational grid. See Manual.

        ENDIF
 2300 CONTINUE
C
c.......................................................................
C
!pb 22.01.2014 change of source strength for gas puff as required
!pb            from SOLPS moved here from subroutine EIRENE_MCARLO
csw
csw 24oct2011
      if(my_pe==0) then
        if (.not.allocated(flux_save)) then
           allocate(flux_save(nstra))
           flux_save = 0._dp
        end if
        if (.not.allocated(npts_save)) then
           allocate(npts_save(nstra))
           allocate(npts_saved(nstra))
           npts_saved = .false.
           npts_save = 0
        endif

        DO ISTRAI=NTARGI+1,NSTRAI
          ISTRA = ISTRAI
csw 20mar2013
csw       if(.not.nlvol(istra)) then
          if(.not.nlvol(istra) .and. .not. nlcns(istra)) then
            IF (FLUX_save(ISTRA).NE.0.d0 ) THEN
              IF (NPTS_SAVED(ISTRA)) NPTS(ISTRA) = NPTS_SAVE(ISTRA)
              FLUX(ISTRA)=FLUX_save(ISTRA)*ELCHA
            ELSE
              IF (.NOT.NPTS_SAVED(ISTRA)) THEN
                 NPTS_SAVE(ISTRA) = NPTS(ISTRA)
                 NPTS_SAVED(ISTRA) = .TRUE.
              ENDIF
              NPTS(ISTRA)=0
!pb 22.01.2014 Do NOT overwrite sources which are not mentioned in
!pb            SOLPS input

!pb            ELSE
!pb              FLUX(ISTRA)=1.
            ENDIF
          endif
        ENDDO
      endif
csw

      if(allocated(ixtri).and..not.nltrimesh) then
csw 04dec2014 collecting normals for B2.5/B2 cells per triangle
        write(iunout,*) 'IF1COP: collecting normals'
        if(allocated(plnxtri)) then
          deallocate(plnxtri, plnytri, pplnxtri, pplnytri)
        endif
        allocate(plnxtri(ntrii))
        allocate(plnytri(ntrii))
        allocate(pplnxtri(ntrii))
        allocate(pplnytri(ntrii))
        plnxtri = 0._dp
        plnytri = 0._dp
        pplnxtri= 0._dp
        pplnytri= 0._dp
        do it=1,ntrii
          ic=ixtri(it)
          if (ic<=0) cycle
          plnxtri(it) = pvx(ic)
          plnytri(it) = pvy(ic)
          pplnxtri(it)= pux(ic)
          pplnytri(it)= puy(ic)
        enddo !it
      endif

cnh transfer macroscopic length scale to EIRENE internal variable
      L_MACRO = L_MACRO_B25

      RETURN
C
C 2999 CONTINUE
 1986 WRITE(IUNOUT,*)
     w      "ERROR IN INFCOP: CANNOT OPEN "//FORT//"31 ",
     w      "(PLASMA BACKGROUND)"
      CALL EIRENE_EXIT_OWN (1)
      RETURN
      END SUBROUTINE EIRENE_IF1COP
C
C  PLASMA PROFILES ARE NOW READ IN
C
      SUBROUTINE EIRENE_IF2COP(ITRG)
      USE EIRMOD_COMSOU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ITRG
      INTEGER :: ISP, IT, ITARG, IS, IS1, I34,
     .           NEM, JPLS, IPL, IPLV, IPLSV, IFL,
     .           IG, IGITT, IIPLS, IEPLS, IPLSTI,
     .           ITRI, NDTARG, IFC
      INTEGER, SAVE :: LTARG
      INTEGER :: NRWL(NSTRA)

      REAL(DP) :: EFLX(NSTRA)
      REAL(DP) :: DI(NPLS), VP(NPLS), ZI(NPLS)
      REAL(DP) :: VPX, VPY, VPZ, VR, VT, TI, VTEST, VTEST2, RP1,
     .            DELTE_PARA, DELTI_PARA,
     .            CS, GAMMA, DELY, DRR,
     .            TE, TES, TIS, THMAX, EESHT, EEMAX, EMAXW, CUR,
     .            PARW, PERW, PARWI, PERWI, PM1, PN1,
     .            EADD, ESHEATH, ESUM
      REAL(DP), ALLOCATABLE, SAVE ::
     .            TORL(:,:), ESHT(:,:), ORI(:,:)
      TYPE(FACE), POINTER :: CURFC

      CHARACTER(6) :: CITARG
C
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      REAL(DP), EXTERNAL :: EIRENE_STEP, EIRENE_EMAXW
      EXTERNAL :: EIRENE_LEER, EIRENE_MASR1, EIRENE_MASJ6,
     .            EIRENE_FTCRI, EIRENE_EXIT_OWN
C

      ITARG=ITRG
      IF (ITARG.GT.NTARGI) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'SOURCE DATA FOR STRATUM ISTRA= ',ITARG
        WRITE (iunout,*)
     .    'CANNOT BE DEFINED IN IF2COP. CHANGE INDSRC(ISTRA)'
        CALL EIRENE_LEER(1)
        RETURN
      ENDIF
C
C  NEXT DEFINE FLUXES, TEMPERATURES AND VELOCITIES AT THE TARGETS
C  (FLUXES IN AMP/(CM ALONG TARGET), TEMPERATURES IN EV, VELOCITIES IN CM/SEC)
C   FNIXB*FL (FNIYB*FL) ARE GIVEN IN AMP
C  STATEMENT NO 3000 ---> 3999
C
C 3000 CONTINUE
C
      IF (TRCINT.AND.LTARG.EQ.0) THEN
        LTARG=1
        WRITE (iunout,*) 'ITARG: TARGET NUMBER'
        WRITE (iunout,*)
     .        'NPBS : BRAAMS (SURFACE) X-CELL INDEX OF TARGET'
        WRITE (iunout,*) 'NPES : POLOIDAL SURFACE INDEX OF TARGET'
        WRITE (iunout,*) '       IN EIRENE MESH'
        CALL EIRENE_LEER(1)
      ENDIF
C
      DO JPLS=1,NPLSI
       DO IGITT=1,NGITT
        FLSTEP(JPLS,ITARG,IGITT)=0.
       END DO
      END DO
C
      ALLOCATE (TORL(NSTRA,NGITT))
      ALLOCATE (ESHT(NSTEP,NGITT))
      ALLOCATE (ORI(NSTEP,NGITT))

      RRSTEP(ITARG,1)=0.
      IG=0
      IIPLS=NPLSI
      IEPLS=1
      NDTARG = NDT(ITARG)
      IF (NDTARG < 0) NDTARG = ABS(NDTARG) + NLIM
!
!     filnam(1:9)='vel_targ_'
!     write (filnam(10:11),'(i2.2)') itarg
!     open (newunit=iov,file=filnam)
C
      DO IT=1,TRGT(ITARG)%NTRGDAT
        IFC = TRGT(ITARG)%FACES(IT)
        CURFC => NCLFC(IFC)%FP
        DO WHILE (ASSOCIATED(CURFC))
          ITRI = CURFC%TRIANGLE
          IS = CURFC%SIDE
          IF (INMTI(IS,ITRI) /= NDTARG) THEN
            WRITE (IUNOUT,*) ' TARGET FACE IS NOT PART OF NONDEF.'//
     .                       ' STANDARD SURFACE'
            WRITE (IUNOUT,*) ' ITARG=  ', ITARG
            WRITE (IUNOUT,*) ' IFC=    ', IFC
            WRITE (IUNOUT,*) ' ITRI=   ', ITRI
            WRITE (IUNOUT,*) ' IS=     ', IS
            WRITE (IUNOUT,*) ' NDT=    ', NDT(ITARG)
          END IF
          IS1=IS+1
          IF (IS1.GT.3) IS1=1
          IG=IG+1
          IF (IG.GT.NGITT) GOTO 999
C  TESTEP, TISTEP: ZONE-CENTERED TEMPERATURE IN BOUNDARY ZONE (EV)
          TESTEP(ITARG,IG) = TRGT(ITARG)%TET(IT)
C  RRSTEP,IRSTEP,IPSTEP: GEOMETRICAL INFORMATION ALONG TARGET
C  RRSTEP IS THE ARC LENGTH ALONG THE TARGET (CM)
          RRSTEP(ITARG,IG+1)=RRSTEP(ITARG,IG) +
     .        SQRT((XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .             (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)
C  EIRENE CELL NUMBER INFORMATION ALONG TARGET
          IRSTEP(ITARG,IG)=ITRI
          IPSTEP(ITARG,IG)=IS
          ITSTEP(ITARG,IG)=1
          IASTEP(ITARG,IG)=0
          IBSTEP(ITARG,IG)=1
          IGSTEP(ITARG,IG)=IFC
          IF (INMTI(IS,ITRI).EQ.0) THEN
            WRITE (iunout,*) 'ERROR IN INFCOP '
            WRITE (iunout,*) 'SOURCE NOT ON A KNOWN SURFACE'
            WRITE (iunout,*) 'ITARG,IG ',ITARG,IG
          ENDIF
C FACTOR FOR SHEATH POTENTIAL: Vsh=Te*SHSTEP(ITARG,IG)
          SHSTEP(ITARG,IG)=TRGT(ITARG)%FSH(IT)
C  TORL: TOROIDAL LENGTH (CM) AT TARGET SEGMENT IY: CENTER OF GRAVITY
          TORL(ITARG,IG)=2.*PIA*0.5*(XTRIAN(NECKE(IS,ITRI))+
     .                               XTRIAN(NECKE(IS1,ITRI)))
          DO JPLS=1,NPLSI
            IPLSTI=MPLSTI(JPLS)
            IPLSV=MPLSV(JPLS)
            ELSTEP(JPLS,ITARG,IG)=0.
            IF (NCHRGP(IPLSTI).eq.0) THEN !! neutral background species
              TISTEP(IPLSTI,ITARG,IG) = TRGT(ITARG)%TNT(IT)
            ELSE                          !! ion background species
              TISTEP(IPLSTI,ITARG,IG) = TRGT(ITARG)%TIT(IT)
            ENDIF

C  DISTEP: ZONE-CENTERED DENSITY IN BOUNDARY ZONE
            IFL=IFLB(JPLS)
            IF (IFL.LE.0.OR.IFL.GT.NFLA) CYCLE
            DISTEP(JPLS,ITARG,IG)=TRGT(ITARG)%DNIT(IT,IFL)*D(JPLS)
C  ZISTEP: ZONE-CENTERED CHARGE IN BOUNDARY ZONE
            IF (TRGT(ITARG)%ZIT(IT,IFL).NE.ZVAC) THEN
              ZISTEP(JPLS,ITARG,IG)=TRGT(ITARG)%ZIT(IT,IFL)
            ELSE
              ZISTEP(JPLS,ITARG,IG)=DBLE(NCHRGP(JPLS))
            END IF
C  FLSTEP: SURFACE-CENTERED FLUX (AMP/CM ALONG TARGET)
            IF (NSPZI(ITARG).LE.IFL.AND.
     .                          IFL.LE.NSPZE(ITARG)) THEN
              IIPLS=MIN0(IIPLS,JPLS)
              IEPLS=MAX0(IEPLS,JPLS)
!  LENGTH OF CELL FACE OF B2 CELL
              DELY = TRGT(ITARG)%FLENGTH(IT)
              FLSTEP(JPLS,ITARG,IG)=0.
              IF (DELY.GT.0.) THEN
                FLSTEP(JPLS,ITARG,IG)=MAX(0._DP,
     .                                    TRGT(ITARG)%FLUX(IT,IFL))
     .                                    *FL(JPLS)/DELY

C  SET DEFAULT ION ENERGY FLUXES FROM B2 BOUNDARY CONDITIONS
!pb              delti_para=3
!pb              delte_para=0.5
                delti_para = trgt(itarg)%fi(it,ifl)
                delte_para = trgt(itarg)%fe(it)
!pb 20.09.2011 moved up
!pb csw 26sep2011
!pb                delta_sheathxb=3.1
!pb                delta_sheathyb=3.1
!pb csw
                tis=TISTEP(IPLSTI,ITARG,IG)
                tes=TESTEP(ITARG,IG)
!               ELSTEP(JPLS,ITARG,IG) = ELSTEP(JPLS,ITARG,IG)+
!     .                                FL(JPLS)/DELY*
!     .                                TRGT(ITARG)%FE(IT)*
!     .                                TRGT(ITARG)%FLUX(IT,IFL)
                ELSTEP(JPLS,ITARG,IG) = ELSTEP(JPLS,ITARG,IG)+
     .                                FL(JPLS)/DELY*
     .                                (TIS*delti_para +
     .                                 TES*delte_para) *
     .                                 TRGT(ITARG)%FLUX(IT,IFL)
              ENDIF
            ENDIF
C
C  VXSTEP,VYSTEP,VZSTEP: SURFACE-CENTERED FLOW VELOCITY (CM/S)
            VXSTEP(IPLSV,ITARG,IG)=TRGT(ITARG)%VXT(IT,IFL)*V
            VYSTEP(IPLSV,ITARG,IG)=TRGT(ITARG)%VYT(IT,IFL)*V
            VZSTEP(IPLSV,ITARG,IG)=TRGT(ITARG)%VZT(IT,IFL)*V
!           write (iov,'(i5,6es14.6)') ig,
!    .         XTRIAN(NECKE(IS,ITRI)), YTRIAN(NECKE(IS,ITRI)),
!    .         XTRIAN(NECKE(IS1,ITRI)), YTRIAN(NECKE(IS1,ITRI)),
!    .         VXSTEP(IPLSV,ITARG,IG), VYSTEP(IPLSV,ITARG,IG)
          END DO
          CURFC => CURFC%NEXT
        ENDDO
C
      END DO  ! IT
      NRWL(ITARG)=IG+1
C
      IF (TRCINT) CALL EIRENE_LEER(2)
C
C  INITIALIZE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
C  SET SOME SOURCE PARAMETERS EXPLICITLY TO ENFORCE INPUT CONSISTENCY
C  also: sum over species: flstep(0,...), elstep(0,...) will be set.
C
      FLUX(ITARG)=EIRENE_STEP(IIPLS,IEPLS,NRWL(ITARG),ITARG,4)
C
CNH DEFINING THE SPECIES TYPE FOR THE STRATUM
      IF (ITYPE(ITARG).EQ.1) THEN     ! ATOM SOURCE
        NLATM(ITARG) =.TRUE.
        NLMOL(ITARG) =.FALSE.
        NLION(ITARG) =.FALSE.
        NLPLS(ITARG) =.FALSE.
        NLPHOT(ITARG)=.FALSE.
      ELSEIF (ITYPE(ITARG).EQ.2) THEN ! MOLECULE SOURCE
        NLATM(ITARG) =.FALSE.
        NLMOL(ITARG) =.TRUE.
        NLION(ITARG) =.FALSE.
        NLPLS(ITARG) =.FALSE.
        NLPHOT(ITARG)=.FALSE.
      ELSEIF (ITYPE(ITARG).EQ.3) THEN ! TRACE ION SOURCE
        NLATM(ITARG) =.FALSE.
        NLMOL(ITARG) =.FALSE.
        NLION(ITARG) =.TRUE.
        NLPLS(ITARG) =.FALSE.
        NLPHOT(ITARG)=.FALSE.
      ELSEIF (ITYPE(ITARG).EQ.4) THEN ! BULK ION SOURCE
        NLATM(ITARG) =.FALSE.
        NLMOL(ITARG) =.FALSE.
        NLION(ITARG) =.FALSE.
        NLPLS(ITARG) =.TRUE.
        NLPHOT(ITARG)=.FALSE.
      ELSEIF (ITYPE(ITARG).EQ.5) THEN ! PHOTON SOURCE
        NLATM(ITARG) =.FALSE.
        NLMOL(ITARG) =.FALSE.
        NLION(ITARG) =.FALSE.
        NLPLS(ITARG) =.FALSE.
        NLPHOT(ITARG)=.TRUE.
      ELSE                     ! BULK ION SOURCE (DEFAULT)
        NLATM(ITARG) =.FALSE.
        NLMOL(ITARG) =.FALSE.
        NLION(ITARG) =.FALSE.
        NLPLS(ITARG) =.TRUE.
        NLPHOT(ITARG)=.FALSE.
      ENDIF
C
      NLSRF(ITARG)=.TRUE.
      NLPNT(ITARG)=.FALSE.
      NLLNE(ITARG)=.FALSE.
      NLVOL(ITARG)=.FALSE.
      NLCNS(ITARG)=.FALSE.
C
      NSRFSI(ITARG)=1
CTRIG A
C  IN TRIA OPTION INDIM=4 (for LEVGEO=3) CORRESPONDS TO INDIM=1 (for LEVGEO=4)
      INDIM(1,ITARG)=1
CTRIG E
      IF (INDSRC(ITARG).NE.6) THEN
        I34=EIRENE_IDEZ(INT(SORLIM(1,ITARG)),3,3)
        SORLIM(1,ITARG)=I34*100+40   !sample with step fct.
      ELSEIF (INDSRC(ITARG).EQ.6) THEN
C  SORLIM DEFAULT WAS 0.D0
        SORLIM(1,ITARG)=0240
        SORIND(1,ITARG)=ITARG
CNH ADAPTATION IN CASE OF AN ATOMIC/MOLECULAR SOURCE ==>
C   TO BE SAMPLED FROM FLUID NEUTRAL BACKGROUND
C   DISTRIBUTION INSTEAD OF THE ION DISTRIBUTION
        IF (NLATM(ITARG).OR.NLMOL(ITARG)) THEN
          IF (NSPZI(ITARG).NE.NSPZE(ITARG)) THEN
            WRITE (iunout,*) 'ERROR IN IF2COP: NSPZI HAS TO BE EQUAL ',
     .        'TO NSPZE FOR ATOMIC/MOLECULAR SOURCE'
            CALL EIRENE_EXIT_OWN(1)
          ELSE
            SORIND(1,ITARG) = 100*NSPZI(ITARG) + ITARG
          ENDIF
        ENDIF
      ENDIF
c      SORIND(1,ITARG)=ITARG

C  IN CASE INDIM=4: INSOR,INDGRD,... ARE REDUNDANT
      NRSOR(1,ITARG)=-1
      NPSOR(1,ITARG)=-1
      IF (INDSRC(ITARG).LT.6) THEN
        WRITE (iunout,*) 'MESSAGE FROM IF2COP:'
        WRITE (iunout,'(a,i2,a)')
     .   'SOURCE STRENGTH AND SPATIAL DISTRIBUTION FOR STRATUM ',
     .    ITARG,' MODIFIED.'
        CALL EIRENE_MASR1('FLUX=   ',FLUX(ITARG))
        WRITE (iunout,*) 'USE STEP FUNCTION ISTEP= ',ITARG,
     .                   ' FROM BLOCK 14'
        WRITE (iunout,*) ' FROM BLOCK 14: SORLIM= ',SORLIM(1,ITARG)
        CALL EIRENE_LEER(1)
      ENDIF
C
      IF (INDSRC(ITARG).EQ.6) THEN
C  DEFINE SOURCE FOR TARGET RECYCLING STRATUM ITARG
C  ASSUME NOW: ITARG=ISTRA
C  DEFAULTS ARE ALREADY SET IN SUBR. INPUT.
C
        CALL EIRENE_FTCRI(ITARG,CITARG)
        TXTSOU(ITARG)= 'SURFACE RECYCLING SOURCE NO.'//CITARG
        NPTS(ITARG)=INT(NPTC(ITARG)*MPTS_COMSOU)
        IF (NPTS(ITARG) < 0) NPTS(ITARG) = HUGE(1)
!VK MINIMUM NUMBER OF HISTORIES FOR THE STRATUM
        NMINPTS(ITARG)=INT(NPTCM(ITARG)*MPTS_COMSOU)
        IF (NTIME > 0) THEN
          IF (NINITL(NSTRAI).LT.0.AND.(.NOT.NLCRR)) THEN
            NINITL(ITARG)=NINITL(NSTRAI)
          ELSE
            NINITL(ITARG)=ITARG*1001
          END IF
        ELSE
          IF (NINITL_READ.LT.0.AND.(.NOT.NLCRR)) THEN
            NINITL(ITARG)=NINITL_READ
          ELSE
            NINITL(ITARG)=ITARG*1001
          END IF
        END IF
        NSPEZ(ITARG)=0
!pb        SORIFL(1,ITARG)=NIFLG(ITARG,1)  ??
        SORIFL(1,ITARG)=0
        SORWGT(1,ITARG)=1.
C  USE ENERGY FLUXES SPECIFIED HERE, IE., SORENE, SORENI ARE REDUNDANT
        NEMODS(ITARG)=7  !DEFAULT; INCL. SHEATH ACCELERATION
CNH ADAPTATION FOR AN ATOMIC/MOLECULAR SOURCE ==> NEEDED FOR THE HYBRID APPROACH
C   IT IS ASSUMED THAT THE FIRST ATOM/MOLECULE CORRESPONDS TO THE HYDROGEN
C   ATOM/MOLECULE (NSEPZ=1) ==> SHOULD BE ADAPTED FOR A SIMULATION WITH A
C   COMBINATION OF ISOTOPES
C   MOLECULES: DEFAULT ASSUMPTION OF NON-DRIFTING MAXWELLIAN WITH A TEMPERATURE
C   OF 0.1 EV
        IF (NLATM(ITARG)) THEN
          NSPEZ(ITARG)   = LKINDP(NSPZI(ITARG))
          NEMODS(ITARG)  = 1000*NSPZI(ITARG) + 336
          SORIFL(1,ITARG)= 2000
        ELSE IF (NLMOL(ITARG)) THEN
          NSPEZ(ITARG) =1
          NEMODS(ITARG)=116
          SORENI(ITARG)=0.1_DP
          SORENE(ITARG)=0.1_DP
          SORVDX(ITARG)=0.0_DP
          SORVDY(ITARG)=0.0_DP
          SORVDZ(ITARG)=0.0_DP
        END IF
        IF (NEMOD(ITARG).NE.0) then
          NEMODS(ITARG)=NEMOD(ITARG)
        END IF
        NAMODS(ITARG)=1
C
C  USE POLYGON MESH, IE., SORAD1,...,SORAD4 ARE REDUNDANT.
        SORAD5(1,ITARG)=ZIA
        SORAD6(1,ITARG)=ZAA
C
C  VELOCITY SPACE DISTRIBUTION
        SORCOS(ITARG)=1.
        SORMAX(ITARG)=0.
C
C  DO 2028 LOOP FROM SUBR. INPUT
        THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ITARG)*DEGRAD))
        IF (NAMODS(ITARG).EQ.1) THEN
          RP1=SORCOS(ITARG)+1.
          SORCOS(ITARG)=1./RP1
!pb          IF (ABS(COS(THMAX)).LE.EPS10) THEN
          IF (ABS(COS(THMAX)).LE.EPS5) THEN
            SORMAX(ITARG)=1.
          ELSE
            SORMAX(ITARG)=1.-COS(THMAX)**RP1
          ENDIF
        ELSEIF (NAMODS(ITARG).EQ.2) THEN
          SORCOS(ITARG)=SORCOS(ITARG)*DEGRAD
          SORMAX(ITARG)=THMAX
        ENDIF
        NLSYMT(0)=NLSYMT(0).AND.NLSYMT(ITARG)
        NLSYMP(0)=NLSYMP(0).AND.NLSYMP(ITARG)
C
      ENDIF
C
C  SOURCE DEFINITION FOR TARGET RECYCLING STRATUM ITARG COMPLETED
C
C 3999 CONTINUE
C
C  TARGET DATA ARE DEFINED NOW
C
C
C  COMPUTE MACH NUMBERS, SHEATH POTENTIAL, AND EXACT SURFACE ENERGY FLUXES
C  FOR COMPARISON WITH SAMPLED ENERGY FLUXES
C  E-FLUX "ETOTP". THIS IS ONLY FOR DIAGNOSTICS PURPOSES
C  E.G. TO CHECK CONSISTENCY OF BOUNDARY CONDITIONS
C  STATEMENT NO. 6000 ---> 6500
C
      IF (.NOT.TRCINT) GOTO 6500
C
      EEMAX=0.
      EESHT=0.
C
      NEM=IABS(NEMODS(ITARG))
      DO 6011 IG=1,NRWL(ITARG)-1
C
C  COMPUTE EIRENE SHEATH POTENTIAL ESHT(ITARG,IG)
C  USE ALL NPLSI SPECIES, NOT JUST IFL=NSPZI,NSPZE
C  IN CASE NEMOD=9: USE THE SHEATH PARAMETERS TRANSFERED FROM B2.5
C
        ESHT(ITARG,IG)=0.D0
        IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7) THEN
          DO 6005 IPL=1,NPLSI
            IPLV=MPLSV(IPL)
            VPX=VXSTEP(IPLV,ITARG,IG)
            VPY=VYSTEP(IPLV,ITARG,IG)
            VPZ=VZSTEP(IPLV,ITARG,IG)
            VP(IPL)=SQRT(VPX**2+VPY**2+VPZ**2)
            DI(IPL)=DISTEP(IPL,ITARG,IG)
            ZI(IPL)=ZISTEP(IPL,ITARG,IG)
 6005     CONTINUE
          TE=TESTEP(ITARG,IG)
          CUR=0.
          GAMMA=0.
          IF (SHSTEP(ITARG,IG).NE.0.0_DP) THEN
            ESHT(ITARG,IG)=SHSTEP(ITARG,IG)*TE
          ELSE
            ESHT(ITARG,IG)=EIRENE_SHEATH(TE,DI,VP,ZI,GAMMA,CUR,
     .                         NPLSI,-ITARG)
          ENDIF
        ELSE IF (NEM == 9) THEN
          NPBS=IGSTEP(ITARG,IG)
          TE=TESTEP(ITARG,IG)
          ESHT(ITARG,IG)=DELTA_SHEATHB(NPBS)*TE
        ENDIF
C
        ISP = 0
        DO 6009 JPLS=1,NPLSI
          IF (FLSTEP(JPLS,ITARG,IG).EQ.0.D0) GOTO 6009
C
          IPLSTI=MPLSTI(JPLS)
          IPLSV=MPLSV(JPLS)
C  CHECK BOHM CRITERION AT "POLOIDAL" TARGET SURFACE COMPONENTS
          ITRI=IRSTEP(ITARG,IG)
          NPES=IPSTEP(ITARG,IG)
          VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(JPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO POLOIDAL TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
C  ASSUMING ORTHOGONAL TARGET
          PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLSV,ITARG,IG)+
     .         PTRIY(NPES,ITRI)*VYSTEP(IPLSV,ITARG,IG))
C  VELOCITY COMPONENT PARALLEL TO POLOIDAL TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
C  AGAIN: ASSUMING ORTHOGONAL TARGET
          VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PTRIX(NPES,ITRI)
          VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PTRIY(NPES,ITRI)
          VPZ=VZSTEP(IPLSV,ITARG,IG)-0.
          PN1=SQRT(VPX**2+VPY**2+VPZ**2)
          PERW=0.
          PARW=0.
          IF (VT.GT.0.) THEN
            PERW=PM1/VT
            PARW=PN1/VT
          ENDIF
C
          CS=SQRT((1.*TISTEP(IPLSTI,ITARG,IG)+
     .                TESTEP(ITARG,IG))/BMASS(JPLS))*CVEL2A
C THE MACH NUMBER BOUNDARY CONDITION ONLY AFFECTS THE PARALLEL TO B
C MOMENTUM, I.E., NOT THE RADIAL VELOCITY
          VTEST=SQRT(PM1**2+VPZ**2)
          VTEST=VTEST/(CS+EPS60)
          VR=SQRT(VPX**2+VPY**2)
          VTEST2=VPZ/(CS+EPS60)
          IF (TRCINT.AND.ABS(VTEST2/VTEST).GT.EPS12) THEN
            ISP = ISP + 1
            WRITE (iunout,'(A,I4,2I5,1P,2E12.4)')
     .       'IPLS,ITARG,IG,MACH_PAR,MACH_Z ',
     .        JPLS,ITARG,IG,VTEST,VTEST2
C           WRITE (iunout,'(A,1P,3E12.4)') 'POL., TOR., RAD. (CM/S)',
C                                           PM1,  VPZ,  VR
          ENDIF
C
C  BOHM CRITERION CHECK DONE
C
C  NEXT: TARGET MAXW. ENERGY FLUXES
C  EADD=  IN EV, SUCH THAT EADD*PARTICLE FLUX = ENERGY FLUX
          DRR=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
C  ENERGY FLUX DEFINED WITH PARAMETERS IN INPUT BLOCK 7
          IF (NEM.EQ.1) THEN
            EADD=SORENI(ITARG)
          ELSEIF (NEM.EQ.2.OR.NEM.EQ.3) THEN
            EADD=SORENI(ITARG)*TISTEP(IPLSTI,ITARG,IG)+SORENE(ITARG)*
     .           TESTEP(ITARG,IG)
          ELSEIF (NEM.GE.4 .AND. NEM.LE.7) THEN
            PERWI=PERW/SQRT(BMASS(JPLS)/RMASSP(JPLS))
            PARWI=PARW/SQRT(BMASS(JPLS)/RMASSP(JPLS))
            TI=TISTEP(IPLSTI,ITARG,IG)
            EADD=EIRENE_EMAXW(TI,PERWI,PARWI)
          ELSEIF (NEM.EQ.8 .OR. NEM.EQ.9) THEN
C  ENERGY FLUX ELSTEP IS ALREADY DEFINED BY B2 BOUNDARY CONDITIONS (SUM: JPLS=0?)
            EADD=ELSTEP(JPLS,ITARG,IG)/FLSTEP(JPLS,ITARG,IG)
          ENDIF
          EMAXW=EADD
          ESUM=EMAXW*FLSTEP(JPLS,ITARG,IG)
          EEMAX=EEMAX+ESUM*DRR

C  ADD ENERGY GAIN BY SHEATH ACCELERATION TO TOTAL
          IF (NEM.EQ.3.OR.NEM.EQ.5.OR.NEM.EQ.7.OR.NEM.EQ.9) THEN
            ESHEATH=ZISTEP(JPLS,ITARG,IG)*ESHT(ITARG,IG)
            ESUM=ESHEATH*FLSTEP(JPLS,ITARG,IG)
            EESHT=EESHT+ESUM*DRR
          ENDIF

 6009   CONTINUE  ! JPLS loop
        IF (TRCINT.AND.ISP.GT.1) CALL EIRENE_LEER(1)
 6011 CONTINUE    ! IG, CELL ALONG TARGET
C
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ITARG
      WRITE (iunout,*) TXTSOU(ITARG)
      WRITE (iunout,'(1X,A4,10A11,A3)')
     . 'IG','ARC','P-FLUX','E-FLUX','TE','TI','SHEATH/TE','SHSTEP',
     . 'VXSTEP','VYSTEP','VZSTEP'
      DO 6100 IG=1,NRWL(ITARG)-1

        WRITE (iunout,'(1X,I4,1P,10E11.3)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELSTEP(0,ITARG,IG),
     .             TESTEP(ITARG,IG),TISTEP(1,ITARG,IG),
     .             ESHT(ITARG,IG)/(TESTEP(ITARG,IG)+EPS60),
     .             SHSTEP(ITARG,IG),
     .             VXSTEP(1,ITARG,IG),VYSTEP(1,ITARG,IG),
     .             VZSTEP(1,ITARG,IG)
 6100 CONTINUE
      WRITE (iunout,'(1X,I4,1P,1E11.3)') NRWL(ITARG),
     .                                 RRSTEP(ITARG,NRWL(ITARG))
      CALL EIRENE_MASR1 ('EEMAX   ',EEMAX)
      CALL EIRENE_MASR1 ('EESHT   ',EESHT)
      flush(iunout)
C
      EFLX(ITARG)=EEMAX+EESHT
      WRITE (iunout,*) 'PARTICLE FLUX(IPLS), IPLS=1,NPLSI '
      WRITE (iunout,'(1X,1P,6E12.4)') (FLTOT(ISP,ITARG),ISP=1,NPLSI)
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'ENERGY FLUX '
      WRITE (iunout,'(1X,1P,1E12.4)') EFLX(ITARG)
      CALL EIRENE_LEER(2)
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
 6500 CONTINUE
      flush(iunout)

      DEALLOCATE (TORL)
      DEALLOCATE (ESHT)
      DEALLOCATE (ORI)
C
C
      RETURN
  999 CONTINUE
      WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL'
      CALL EIRENE_EXIT_OWN(1)
      RETURN
      END SUBROUTINE EIRENE_IF2COP
C
C
      SUBROUTINE EIRENE_IF3COP(IENTRY,LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
      IMPLICIT NONE
      INTEGER :: IFRST
      INTEGER, INTENT(IN) :: IENTRY, ISTRAA, ISTRAE, NEW_ITER
      LOGICAL, INTENT(INOUT) :: LSTP
      REAL(DP), SAVE :: SEES, SEIS,
     .          SMOCL_VEC(0:2), SNIRES, SMORES, SEERES, SEIRES,
     .          SUMM_VEC(0:2), SUMN, SUMEI, SUMEE, SUMN_OLD,
     .          TEST, FLX, FLXI,
     .          RECADD, RECTOT, PIADD_VEC(0:2),
     .          PRESS00, SNICL, FLX_EIR,
     .          EAMISUM, EEADD, EIADD, CHE, CHI, CHP, CHEES, CHEIS
C
      logical :: lhit(nrad)
      INTEGER :: I, JATM, JMOL, JION, JPLS, IPLSV, IPLSTI, ISTRAI,
     .           IN, IT, IC, IFL, INC, IIRC, IRRC, IST_RATE, IER,
     .           ISTAT_COP, ICPV, ICOND,
     .           KK, KKK, ICO
      INTEGER, SAVE :: IFIRST

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
      TYPE(RATE_STORE), POINTER :: RTIS
C
      REAL(DP), EXTERNAL :: EIRENE_FTABRC1, EIRENE_FEELRC1
      EXTERNAL :: EIRENE_MASR6, EIRENE_INDMPI, EIRENE_LEER,
     .            EIRENE_SAVE_TALLIES
C
      IF (IENTRY.EQ.0) THEN
        WRITE (iunout,*) 'IF3COP IS CALLED, ISTRAA,ISTRAE'
        WRITE (iunout,*) ISTRAA,ISTRAE
        LSHORT=.FALSE.
        LSTP3=.TRUE.
        LSTOP=LSTP3
        IFIRST=0
C
      ELSE IF (IENTRY.EQ.1) THEN
C
C  ENTRY FOR SHORT CYCLE FROM SUBR. EIRSRT
C
C  IFIRST=0: RESTORE DATA FROM A PREVIOUS EIRENE RUN, SET REFERENCE
C            DATA FOR "STOP CRITERION" SNIS,SEES,SEIS
C  IFIRST>0: MODIFY SOURCE TERMS ACCORDING TO NEW PLASMA CONDITIONS,
C            COMPARE INTEGRALS WITH SNIS,...., AND DECIDE TO STOP OR
C            CONTINUE SHORT CYCLE (LSTOP)
C
        LSHORT=.TRUE.
        LSTP3=LSTP
        LSTOP=LSTP
        IFIRST=IFRST
C
      END IF

csw 14jul2011
      IF (ISTRAA.ne.ISTRAE.and..not.LSHORT) THEN
          WRITE (IUNOUT,*) 'Error, this version of INFCOP '//
     &                     'should process one stratum at a time'
          CALL MPI_ABORT(MPI_COMM_WORLD, -1, IER)
         ! note that calling exit_own would lead to a deadlock
         ! because that has collective mpi calls
      ENDIF

      if (.not. I_am_leader(istraa)) return

      IF (.NOT.ALLOCATED(CHPS)) THEN
        ALLOCATE (CHPS(NFL))
        ALLOCATE (SNIS(0:NFL))
        ALLOCATE (CHMOS(NFL))
        ALLOCATE (SMOS(3*NFL))
        ALLOCATE (SCALN(0:NFL))
        ALLOCATE (SEES0(NSTRA))
        ALLOCATE (SEIS0(NSTRA))
        ALLOCATE (SNIS0(NSTRA,0:NFL))
        ALLOCATE (SMOS0(NSTRA,-3:3*NFL))

        ALLOCATE (RESSNI(0:NSTRA,NFL))
        ALLOCATE (RESSMO(0:NSTRA,NFL))
        ALLOCATE (RESSEE(0:NSTRA))
        ALLOCATE (RESSEI(0:NSTRA))

        ALLOCATE (FLXEIR(NSTRA))

        CALL EIRENE_ALLOC_BRASPOI
        CALL EIRENE_ALLOC_EIRBRA(NCVP,NFL,NATM,NMOL,NION,NSTRA)
C
        RESSNI = 0._DP
        RESSMO = 0._DP
        RESSEE = 0._DP
        RESSEI = 0._DP
      END IF
      if (.not.allocated(chpm)) then
cdr  next 4: for short cycle. Obsolete
        ALLOCATE (CHPM(NPLS,NRAD))
        ALLOCATE (CHEEM(NRAD))
        ALLOCATE (CHEIM(NRAD))
        ALLOCATE (CHMOM(NPLS,NRAD))
cdr  primary volumetric source contributions, from volumetric sources
cdr  Here: not sampled, but exact
cdr  Always needed for vol. rec. strata. Not only for short cycle
        ALLOCATE (PPPL_COP(NPLS,NRAD))
        ALLOCATE (MPPL_COP(3*NPLS,NRAD))
        ALLOCATE (EPPL_COP(NPLS,NRAD))
        ALLOCATE (EPEL_COP(NRAD))
cdr  next 4: for short cycle. Obsolete
        ALLOCATE (PPLODA(NPLS,NRAD))
        ALLOCATE (MPLODA(3*NPLS,NRAD))
        ALLOCATE (EPLODA(NPLS,NRAD))
        ALLOCATE (EPEODA(NRAD))
      END IF
      if (.not.allocated(eapl0)) then
         allocate (eapl0(nrad))
         allocate (empl0(nrad))
         allocate (eipl0(nrad))
      end if
C
      IF (NEW_ITER == 0) THEN
        RESSNI = 0.D0
        RESSMO = 0.D0
        RESSEE = 0.D0
        RESSEI = 0.D0
      ENDIF
C
      IF (.NOT.LSHORT) THEN
        RESSNI(ISTRAA:ISTRAE,:) = 0._DP
        RESSMO(ISTRAA:ISTRAE,:) = 0._DP
        RESSEE(ISTRAA:ISTRAE) = 0._DP
        RESSEI(ISTRAA:ISTRAE) = 0._DP
      ENDIF

      volSUMN(ISTRAA:ISTRAE)=0.0           !dpc
      volSUMM(0:2,ISTRAA:ISTRAE)=0.0       !dr
      volSUMEI(ISTRAA:ISTRAE)=0.0          !dpc
      volSUMEE(ISTRAA:ISTRAE)=0.0          !dpc

      DO 10000 ISTRAI=ISTRAA,ISTRAE
C
csw 20mar2013
        IF (XMCP(ISTRAI).LE.1.) then
         sni(:,:,istrai) = 0.d0
         sne_pael(:,istrai) = 0.d0 !djm Jan2017
         sne_pmel(:,istrai) = 0.d0
         sne_piel(:,istrai) = 0.d0
         sna_paat(:,:,istrai) = 0.d0 !djm Jan2017
         sna_pmat(:,:,istrai) = 0.d0
         sna_piat(:,:,istrai) = 0.d0
         sna_paml(:,:,istrai) = 0.d0 !djm Jan2017
         sna_pmml(:,:,istrai) = 0.d0
         sna_piml(:,:,istrai) = 0.d0
         sna_paio(:,:,istrai) = 0.d0 !djm Jan2017
         sna_pmio(:,:,istrai) = 0.d0
         sna_piio(:,:,istrai) = 0.d0
         sni_papl(:,:,istrai) = 0.d0 !djm Jan2017
         sni_pmpl(:,:,istrai) = 0.d0
         sni_pipl(:,:,istrai) = 0.d0
         sni_pppl(:,:,istrai) = 0.d0
         smo(:,:,istrai) = 0.d0
         smo_mapl(:,:,istrai) = 0.d0 !djm Jan2017
         smo_mmpl(:,:,istrai) = 0.d0
         smo_mipl(:,:,istrai) = 0.d0
         smo_mppl(:,:,istrai) = 0.d0
         see(:,istrai) = 0.d0
         see_eael(:,istrai) = 0.d0 !djm Jan2017
         see_emel(:,istrai) = 0.d0
         see_eiel(:,istrai) = 0.d0
         see_epel(:,istrai) = 0.d0
         sei(:,:,istrai) = 0.d0
         sei_eapl(:,:,istrai) = 0.d0 !djm Jan2017
         sei_empl(:,:,istrai) = 0.d0
         sei_eipl(:,:,istrai) = 0.d0
         sei_eppl(:,:,istrai) = 0.d0
         GOTO 10000
        endif
C
C  DATA TRANSFER BACK FROM EIRENE TO EXTERNAL CODE
C  STATEMENT NO 7000 ---> 7999
C
C 7000   CONTINUE
C
        IF (.NOT.LSHORT) CALL EIRENE_SAVE_TALLIES(ISTRAI)
C
C  SCALE SURFACE SOURCES PER UNIT FLUX, FOR OTHER SOURCES USE
C  EIRENE SCALINGS
        IF (ISTRAI.LE.NTARGI.AND.NLPLS(ISTRAI)) THEN
          IF (WTOTP(0,ISTRAI).NE.0.) THEN
C  FLUX FROM EIRENE TO PLASMA CODE: NEGATIVE
            FLX=-WTOTP(0,ISTRAI)
            FLXI=1./FLX
            FLXEIR(ISTRAI)=1._DP
          ELSE
            WRITE (iunout,*)
     .       'NO PLASMA FLUX FROM STRATUM NO. ISTRAI= ',ISTRAI
            WRITE (iunout,*)
     .       'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
            GOTO 7999
          ENDIF
        ELSEIF (ISTRAI.LE.NTARGI.AND.NLATM(ISTRAI)) THEN
          IF (WTOTA(0,ISTRAI).NE.0.) THEN
CNH FOR HYBRID NEUTRAL MODEL ALSO ATOM SOURCE
            FLX=WTOTA(0,ISTRAI)
            FLXI=1./FLX
            FLXEIR(ISTRAI)=1._DP
          ELSE
            WRITE (iunout,*)
     .       'NO ATOM FLUX FROM STRATUM NO. ISTRAI= ',ISTRAI
            WRITE (iunout,*)
     .       'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
            GOTO 7999
          ENDIF
        ELSEIF (ISTRAI.LE.NTARGI.AND.NLMOL(ISTRAI)) THEN
          IF (WTOTM(0,ISTRAI).NE.0.) THEN
CNH FOR KINETIC MOLECULE SIMULATION
            FLX=WTOTM(0,ISTRAI)
            FLXI=1./FLX
            FLXEIR(ISTRAI)=1._DP
          ELSE
            WRITE (iunout,*)
     .       'NO MOLECULE FLUX FROM STRATUM NO. ISTRAI= ',ISTRAI
            WRITE (iunout,*)
     .       'NO DATA RETURNED TO PLASMA CODE FOR THIS STRATUM'
            GOTO 7999
          ENDIF
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.

C  IF THE SOURCE STRENGTH IS TO BE CHANGED DURING THE SHORT CYCLE (E.G.: VOL-REC)
C  THEN FLXEIR HAS TO BE RESET TO SCALE TO NEW SOURCE STRENGTH DURING SHORT CYCLE
          FLXEIR(ISTRAI)=1._DP
        ENDIF
C
C  FIRSTLY INITIALIZE SOURCE TERM ARRAYS

CDR: already done above !

C
        SEE(:,ISTRAI)=0.
        SNE_PAEL(:,ISTRAI)=0. !djm Jan2017
        SNE_PMEL(:,ISTRAI)=0.
        SNE_PIEL(:,ISTRAI)=0.
        SEE_EAEL(:,ISTRAI)=0. !djm Jan2017
        SEE_EMEL(:,ISTRAI)=0.
        SEE_EIEL(:,ISTRAI)=0.
        SEE_EPEL(:,ISTRAI)=0.
        SEI(:,:,ISTRAI)=0.
        SEI_EAPL(:,:,ISTRAI)=0. !djm Jan2017
        SEI_EMPL(:,:,ISTRAI)=0.
        SEI_EIPL(:,:,ISTRAI)=0.
        SEI_EPPL(:,:,ISTRAI)=0.

        SNI(:,:,ISTRAI)=0.
        SNI_PAPL(:,:,ISTRAI)=0. !djm Jan2017
        SNI_PMPL(:,:,ISTRAI)=0.
        SNI_PIPL(:,:,ISTRAI)=0.
        SNI_PPPL(:,:,ISTRAI)=0.
        SMO(:,:,ISTRAI)=0.
        SMO_MAPL(:,:,ISTRAI)=0. !djm Jan2017
        SMO_MMPL(:,:,ISTRAI)=0.
        SMO_MIPL(:,:,ISTRAI)=0.
        SMO_MPPL(:,:,ISTRAI)=0.
C
        SNA_PAAT(:,:,ISTRAI)=0.
        SNA_PMAT(:,:,ISTRAI)=0.
        SNA_PIAT(:,:,ISTRAI)=0.

        SNA_PAML(:,:,ISTRAI)=0.
        SNA_PMML(:,:,ISTRAI)=0.
        SNA_PIML(:,:,ISTRAI)=0.

        SNA_PAIO(:,:,ISTRAI)=0.
        SNA_PMIO(:,:,ISTRAI)=0.
        SNA_PIIO(:,:,ISTRAI)=0.
C
        CHPM  = 0._DP
        CHMOM = 0._DP
        CHEEM = 0._DP
        CHEIM = 0._DP
C
        IF (.NOT.LSHORT) GOTO 7400

        IST_RATE = ITS(ISTRAI)
        RTIS => RTS(IST_RATE)%RTA

cdr begin
cdr this next code snippet is likely something unfinished for primary vol.rec.
cdr momentum sinks for bulk ions?
        COPV=0.D0
        IF (LCOPV) THEN
          CPMUL => COPVS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            ICPV=CPMUL%IART
            IN=CPMUL%ICM
            COPV(ICPV,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF
cdr end

cdr  mapl, mmpl, mipl are now vectorial tallies.
cdr  The code below only deals with the parallel momentum source component.
cdr  This code (from former "short cycle") was obsolete anyway.
cdr  So I have not bothered to update it.
cdr  If it is used accidentally, then we are in big trouble anyway.
cdr  I think we urgently need to separate the (entirely obsolete)
cdr  short cycle snippets from the active code.

        MAPL_VEC=0.D0
        IF (LMAPL_VEC) THEN
          CPMUL => MAPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            MAPL_VEC(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

        MMPL_VEC=0.D0
        IF (LMMPL_VEC) THEN
          CPMUL => MMPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            MMPL_VEC(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

        MIPL_VEC=0.D0
        IF (LMIPL_VEC) THEN
          CPMUL => MIPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            MIPL_VEC(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

        MPHPL_VEC=0.D0
        IF (LMPHPL_VEC) THEN
          CPMUL => MPHPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            MPHPL_VEC(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT IONISATION OF ATOMS
C                        AND BULK ION CHARGE EXCHANGE WITH ATOMS:
C                        PARTICLE AND ENERGY SOURCES
C
C  PARTICLE SOURCE: SPLIT FOR MULTIPLE IPLS SPECIES
        PAPL=0.D0
        IF (LPAPL) THEN
          CPMUL => PAPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            PAPL(IPLS,IN)=CPMUL%VALUEM
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PAPL2(1:NPLS,0:NATM) => PAPL(:,IN)
              PAPL2(IPLS,1:NATMI) = CPMUL%VALUAM(1:NATMI)
#else
              DO JATM = 1, NATMI
                PAPL(IPLS+JATM*NPLS,IN) = CPMUL%VALUAM(JATM)
              END DO
#endif
            END IF
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

C  ELECTRON ENERGY: SINGLE (ELECTRON) SPECIES ARRAY
        EAEL=0.D0
        IF (LEAEL) THEN
          CPSIM => EAELS(ISTRAI)%PSIM
          DO WHILE (ASSOCIATED(CPSIM))
            IN=CPSIM%ICS
            EAEL(IN)=CPSIM%VALUES
            CPSIM => CPSIM%NXTSIM
          END DO
        END IF

C  ION ENERGY: SPLIT FOR MULTIPLE IPLS SPECIES
        EAPL=0.D0
        IF (LEAPL) THEN
          CPMUL => EAPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=         CPMUL%IART
            IN=           CPMUL%ICM
            EAPL(IPLS,IN)=CPMUL%VALUEM
            CPMUL =>      CPMUL%NXTMUL
          END DO
        END IF

        IF (IFIRST.EQ.0) GOTO 7310

        CPMUL => PDENAS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IATM=CPMUL%IART
          IN=CPMUL%ICM
cdr for particle source, from atoms
          DO JPLS=1,NPLSI
            CHP=CPMUL%VALUEM*
     .          (SPLNWA(IN,IATM,JPLS)-RTIS%SPLODA(IN,IATM,JPLS))*ELCHA
            IF (LPAPL) THEN
              PAPL(JPLS,IN)=PAPL(JPLS,IN)+CHP
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAPL2(1:NPLS,0:NATM) => PAPL(:,IN)
                PAPL2(JPLS,IATM)=PAPL2(JPLS,IATM)+CHP
#else
                PAPL(JPLS+IATM*NPLS,IN)=PAPL(JPLS+IATM*NPLS,IN)+CHP
#endif
              END IF
            END IF
            CHPM(JPLS,IN)=CHPM(JPLS,IN)+CHP
          ENDDO
cdr for electron energy source, from atoms
          CHE=CPMUL%VALUEM*
     .        (SEENWA(IN,IATM)-RTIS%SEEODA(IN,IATM))*ELCHA
          IF (LEAEL) EAEL(IN)=EAEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO
cdr  for ion energy source, from atoms:
        CPMUL => COPVS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          ICPV=CPMUL%IART
          IF (ICPV.LE.NPLSI) THEN
            IPLS=ICPV
            IN=CPMUL%ICM
            CHI=CPMUL%VALUEM*
     .          (SEINW(IN,IPLS)-RTIS%SEIOD(IN,IPLS))*ELCHA
            IF (LEAPL) EAPL(IPLS,IN)=EAPL(IPLS,IN)+CHI
            CHEIM(IN)=CHEIM(IN)+CHI
          ENDIF
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7310   CONTINUE

C
C  CORRECTION FOR ELECTRON IMPACT IONISATION AND CX OF ATOMS FINISHED
C
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT DISSOCIATION OF TEST IONS
C

cdr for particle source, from test ions
        PIPL=0.D0
        IF (LPIPL) THEN
          CPMUL => PIPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            PIPL(IPLS,IN)=CPMUL%VALUEM
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PIPL2(1:NPLS,0:NION) => PIPL(:,IN)
              PIPL2(IPLS,1:NIONI)=CPMUL%VALUIM(1:NIONI)
#else
              DO JION = 1, NIONI
                PIPL(IPLS+JION*NPLS,IN) = CPMUL%VALUIM(JION)
              END DO
#endif
            END IF
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

cdr for electron energy source, from test ions
        EIEL=0.D0
        IF (LEIEL) THEN
          CPSIM => EIELS(ISTRAI)%PSIM
          DO WHILE (ASSOCIATED(CPSIM))
            IN=CPSIM%ICS
            EIEL(IN)=CPSIM%VALUES
            CPSIM => CPSIM%NXTSIM
          END DO
        END IF

cdr for ion energy source, from test ions
        EIPL=0.D0
        IF (LEIPL) THEN
          CPMUL => EIPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            EIPL(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          END DO
        END IF

        IF (IFIRST.EQ.0) GOTO 7330

        CPMUL => PDENIS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IION=CPMUL%IART
          IN=CPMUL%ICM
cdr for particle source, from test ions
          DO JPLS=1,NPLSI
            CHP=CPMUL%VALUEM *
     .          (SPLNWI(IN,IION,JPLS)-RTIS%SPLODI(IN,IION,JPLS))*ELCHA
            IF (LPIPL) THEN
              PIPL(JPLS,IN)=PIPL(JPLS,IN)+CHP
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIPL2(1:NPLS,0:NION) => PIPL(:,IN)
                PIPL2(JPLS,IION)=PIPL2(JPLS,IION)+CHP
#else
                PIPL(JPLS+IION*NPLS,IN)=PIPL(JPLS+IION*NPLS,IN)+CHP
#endif
              END IF
            END IF
            CHPM(JPLS,IN)=CHPM(JPLS,IN)+CHP
cdr for ion energy source, from test ions
            CHI=CPMUL%VALUEM *
     .          (SEINWI(IN,IION)-RTIS%SEIODI(IN,IION))*ELCHA
            IF (LEIPL) EIPL(JPLS,IN)=EIPL(JPLS,IN)+CHI
            CHEIM(IN)=CHEIM(IN)+CHI
          END DO
cdr for electron energy source, from test ions
          CHE=CPMUL%VALUEM *
     .        (SEENWI(IN,IION)-RTIS%SEEODI(IN,IION))*ELCHA
          IF (LEIEL) EIEL(IN)=EIEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7330   CONTINUE
C
C
C  CORRECTION FOR TEST IONS FINISHED
C
C
C  SHORT LOOP CORRECTION FOR ELECTRON IMPACT COLLISIONS
C             OF MOLECULES
C
        PMPL=0.D0
        IF (LPMPL) THEN
          CPMUL => PMPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            PMPL(IPLS,IN)=CPMUL%VALUEM
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PMPL2(1:NPLS,0:NMOL) => PMPL(:,IN)
              PMPL2(IPLS,1:NMOLI)=CPMUL%VALUMM(1:NMOLI)
#else
              DO JMOL = 1, NMOLI
                PMPL(IPLS+JMOL*NPLS,IN) = CPMUL%VALUMM(JMOL)
              END DO
#endif
            END IF
            CPMUL => CPMUL%NXTMUL
          ENDDO
        END IF

        EMEL=0.D0
        IF (LEMEL) THEN
          CPSIM => EMELS(ISTRAI)%PSIM
          DO WHILE (ASSOCIATED(CPSIM))
            IN=CPSIM%ICS
            EMEL(IN)=CPSIM%VALUES
            CPSIM => CPSIM%NXTSIM
          END DO
        END IF

        EMPL=0.D0
        IF (LEMPL) THEN
          CPMUL => EMPLS(ISTRAI)%PMUL
          DO WHILE (ASSOCIATED(CPMUL))
            IPLS=CPMUL%IART
            IN=CPMUL%ICM
            EMPL(IPLS,IN)=CPMUL%VALUEM
            CPMUL => CPMUL%NXTMUL
          END DO
        END IF

        IF (IFIRST.EQ.0) GOTO 7350

        CPMUL => PDENMS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IMOL=CPMUL%IART
          IN=CPMUL%ICM
cdr for particle source, from molecules
          DO JPLS=1,NPLSI
            CHP=CPMUL%VALUEM*
     .          (SPLNWM(IN,IMOL,JPLS)-RTIS%SPLODM(IN,IMOL,JPLS))*ELCHA
            IF (LPMPL) THEN
              PMPL(JPLS,IN)=PMPL(JPLS,IN)+CHP
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMPL2(1:NPLS,0:NMOL) => PMPL(:,IN)
                PMPL2(JPLS,IMOL)=PMPL2(JPLS,IMOL)+CHP
#else
                PMPL(JPLS+IMOL*NPLS,IN)=PMPL(JPLS+IMOL*NPLS,IN)+CHP
#endif
              END IF
            END IF
            CHPM(JPLS,IN)=CHPM(JPLS,IN)+CHP
          END DO
cdr for ion energy source, from molecules
cdr apparently missing, although we do set SEINWM in infcop. Bug?
c
cdr for electron energy source, from molecules
          CHE=CPMUL%VALUEM*
     .        (SEENWM(IN,IMOL)-RTIS%SEEODM(IN,IMOL))*ELCHA
          IF (LEMEL) EMEL(IN)=EMEL(IN)+CHE
          CHEEM(IN)=CHEEM(IN)+CHE
          CPMUL => CPMUL%NXTMUL
        ENDDO

 7350   CONTINUE
C  CORRECTION FOR ELECTRON IMPACT DISSOCIATION OF MOLECULES FINISHED

C
C  SHORT LOOP CORRECTION FOR VOLUME RECOMBINATION PROCESSES (UNFINISHED)
C
cdr begin
cdr  probably obsolete: no short cycle for vol. rec strata ever existed?
        PPLODA=0.D0
        CPMUL => PPPL_COPS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          PPLODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        MPLODA=0.D0
        CPMUL => MPPL_COPS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          MPLODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        ENDDO

        EPLODA=0.D0
        CPMUL => EPPL_COPS(ISTRAI)%PMUL
        DO WHILE (ASSOCIATED(CPMUL))
          IPLS=CPMUL%IART
          IN=CPMUL%ICM
          EPLODA(IPLS,IN)=CPMUL%VALUEM
          CPMUL => CPMUL%NXTMUL
        END DO

        EPEODA=0.D0
        CPSIM => EPELS(ISTRAI)%PSIM
        DO WHILE (ASSOCIATED(CPSIM))
          IN=CPSIM%ICS
          EPEODA(IN)=CPSIM%VALUES
          CPSIM => CPSIM%NXTSIM
        END DO
cdr end
C
C
C  SHORT LOOP CORRECTION FINISHED
C
 7400   CONTINUE
C
C
C  ADD CONTRIBUTIONS TO SOURCE RATES, FROM PRIMARY VOLUME RECOMBINATION SOURCE
C
        PPPL_COP = 0.D0
        MPPL_COP = 0.D0
        EPPL_COP = 0.D0
        EPEL_COP = 0.D0

        IF (NLVOL(ISTRAI).AND.NLPLS(ISTRAI)) THEN
C
          RECTOT = 0._DP
cdr  only one bulk ion species per volume source stratum supported
          IPLS=NSPEZ(ISTRAI)  ! RANGE CHECK FOR IPLS
                              ! ALREADY DONE IN SAMVOL
          WRITE (IUNOUT,*) 'DIRECT VOL.REC. CONTRIBUTION, IPLS= ',IPLS

          IPLSTI = MPLSTI(IPLS)
          DO 7472 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            SUMN=0.0
            SUMM_VEC(0:2)=0.0
            SUMEI=0.0
            SUMEE=0.0
            lhit = .false.
            DO 7471 IN=1,NTRII
              INC=NCLTAL(IN)
              if ((nstgrd(in) /= 0) .or. lgvac(in,ipls)) cycle
              IF (NSTORDR >= NRAD) THEN
                RECADD=-TABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                EEADD=  EELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
              ELSE
                RECADD=-EIRENE_FTABRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
                EEADD=  EIRENE_FEELRC1(IRRC,IN)*DIIN(IPLS,IN)*ELCHA
              END IF
              PIADD_VEC(0:2)=0._DP
              IF (LPMOM_VEC) THEN
cdr vectorial direct momentum sink rates via vol. recomb.
                DO ICO=0,2
                  KK=ICO*NPLS
                  PIADD_VEC(ICO)=PMOM_VEC(KK+IPLS,IN)*RECADD
                ENDDO
              ENDIF
              EIADD=1.5*TIIN(IPLSTI,IN)*RECADD
              IF (LEDRIFT) EIADD=EIADD+EDRIFT(IPLS,IN)*RECADD

!pb 21012013
!  if lcoarse: add contribution to tallies only once per fine grid cell
!  tallies are to be scaled with voltal
cdr to be checked again. Not sure that this is correct
              if (.not.lhit(inc)) then
                PPPL_COP(IPLS,INC)=PPPL_COP(IPLS,INC)+RECADD
cdr  mppl_cop is a vectorial tally: par., perp., and diamagn. component
                DO ICO=0,2
                  KK=ICO*NPLS
                  MPPL_COP(KK+IPLS,INC)=MPPL_COP(KK+IPLS,INC)+
     .                                  PIADD_VEC(ICO)
                ENDDO
                EPPL_COP(IPLS,INC)=EPPL_COP(IPLS,INC)+EIADD
                EPEL_COP(INC)=EPEL_COP(INC)+EEADD
                lhit(inc) = .true.
              end if

! ALL INPUT TALLIES, TAB.., FTAB.., VOL, ARE ON UNDERLYING FINE GRID.
! and hence: RECADD,PIADD_VEC,EIADD,EEADD also on fine grid.
!  here we scale with volume of triangle cell
!  this needs to be done per triangle
              SUMN=SUMN+RECADD*VOL(IN)
              SUMM_VEC(0:2) =SUMM_VEC(0:2) +PIADD_VEC(0:2) *VOL(IN)
              SUMEI=SUMEI+EIADD*VOL(IN)
              SUMEE=SUMEE+EEADD*VOL(IN)

 7471       CONTINUE  ! loop over grid

            RECTOT = RECTOT + SUMN

            WRITE (iunout,*) 'PARTIAL: IRRC ',IRRC
            CALL EIRENE_MASR6
     .           ('SUMN, SUMM_VEC, SUMEI, SUMEE                    ',
     .             SUMN, SUMM_VEC(0), SUMM_VEC(1), SUMM_VEC(2),
     .             SUMEI, SUMEE)
            volSUMN(ISTRAI)=volSUMN(ISTRAI)+SUMN             ! dpc
            volSUMM(0:2,ISTRAI)=volSUMM(0:2,ISTRAI)+SUMM_VEC(0:2)! dr
            volSUMEI(ISTRAI)=volSUMEI(ISTRAI)+SUMEI          ! dpc
            volSUMEE(ISTRAI)=volSUMEE(ISTRAI)+SUMEE          ! dpc
 7472     CONTINUE  !IRRC

          WRITE (iunout,*) 'TOTAL'
          CALL EIRENE_MASR6
     .         ('SUMN, SUMM_VEC, SUMEI, SUMEE                    ',
     .           volSUMN(ISTRAI),   volSUMM(0,ISTRAI),
     .           volSUMM(1,ISTRAI), volSUMM(2,ISTRAI),
     .           volSUMEI(ISTRAI),  volSUMEE(ISTRAI))
cdr
CC SUMN_OLD=WTOTP ???
          IF (.NOT.LSHORT) SUMN_OLD=RECTOT
csw 08jun2010
          sumn_old=rectot

C  RESCALE PPPL_COP,.... TO SOURCE STRENGTH FROM LAST FULL EIRENE RUN
C  BECAUSE ALSO PAPL,.... ARE SCALED LIKE THIS
          PPPL_COP=PPPL_COP*SUMN_OLD/RECTOT
          MPPL_COP=MPPL_COP*SUMN_OLD/RECTOT
          EPPL_COP=EPPL_COP*SUMN_OLD/RECTOT
          EPEL_COP=EPEL_COP*SUMN_OLD/RECTOT

C  RESCALE SOURCES FOR B2 (PPPL_COP,PAPL,....) TO NEW SOURCE STRENGTH
          FLXEIR(ISTRAI)=RECTOT/SUMN_OLD
cdr
          IF (LSHORT) THEN
            DO JPLS=1,NPLSI
              CHPM(JPLS,1:NSBOX_TAL) = CHPM(JPLS,1:NSBOX_TAL) +
     .             PPPL_COP(JPLS,1:NSBOX_TAL) - PPLODA(JPLS,1:NSBOX_TAL)
              CHEIM(1:NSBOX_TAL) = CHEIM(1:NSBOX_TAL) +
     .             EPPL_COP(JPLS,1:NSBOX_TAL) - EPLODA(JPLS,1:NSBOX_TAL)
            END DO

            CHEEM(1:NSBOX_TAL) = CHEEM(1:NSBOX_TAL) +
     .            EPEL_COP(1:NSBOX_TAL) - EPEODA (1:NSBOX_TAL)
          END IF

C
C  SAVE RECOMBINATION SOURCES FOR USE BY SHORT CYCLE
C
          IF (.NOT.LSHORT) THEN
            DO JPLS=1,NPLSI
              DO IN=1,NSBOX_TAL
                IF (PPPL_COP(JPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = JPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = PPPL_COP(JPLS,IN)
                  CPMUL%NXTMUL => PPPL_COPS(ISTRAI)%PMUL
                  PPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
                IF (MPPL_COP(JPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = JPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = MPPL_COP(JPLS,IN)
                  CPMUL%NXTMUL => MPPL_COPS(ISTRAI)%PMUL
                  MPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
                IF (EPPL_COP(JPLS,IN) .NE. 0.D0) THEN
!pb               ALLOCATE(CPMUL)
                  CPMUL => EIRENE_NEW_MULARR()
                  CPMUL%IART = JPLS
                  CPMUL%ICM = IN
                  CPMUL%VALUEM = EPPL_COP(JPLS,IN)
                  CPMUL%NXTMUL => EPPL_COPS(ISTRAI)%PMUL
                  EPPL_COPS(ISTRAI)%PMUL => CPMUL
                ENDIF
              ENDDO
            ENDDO
            DO IN=1,NSBOX_TAL
              IF (EPEL_COP(IN) .NE. 0.D0) THEN
!pb             ALLOCATE(CPSIM)
                CPSIM => EIRENE_NEW_SIMARR()
                CPSIM%ICS = IN
                CPSIM%VALUES = EPEL_COP(IN)
                CPSIM%NXTSIM => EPELS(ISTRAI)%PSIM
                EPELS(ISTRAI)%PSIM => CPSIM
              ENDIF
            ENDDO
          END IF

        ENDIF
C
csw 26jan2011 extra B25
csw 08mar2013 moved to here, i.e. after correction from volume recombination
        if(.not. lshort) then
          call eirene_wneutrals_fill(istrai)
        endif
csw
        IF (.NOT.LSYMET) GOTO 7500
C
C  SECONDLY SYMMETRISE EIRENE ARRAYS IN CASE OF UP-DOWN SYMMETRY IN MODEL
C
C
C   THIRDLY WRITE EIRENE ARRAYS (1D) ONTO BRAAMS ARRAYS (2D)
C   AND RESCALE TO PROPER UNITS: #/CELL/STRATUM FLUX
C   # STANDS FOR PARTICLES (SNI), MOMENTUM (SMO)
C   AND ENERGY (SEE,SEI)
C
 7500   CONTINUE

        IF (ISTRAI <= NTARGI) THEN
          FLX_EIR = 1._DP
        ELSE
          FLX_EIR = FLXEIR(ISTRAI)
        END IF

cc Add contributions from condensation (kinetic -> fluid)
cc CONDENSE_FROM relates to the atom species (IATM, block 4) from which kinetic particles are removed.
cc CONDENSE_TO relates to the background species (block 5) to which the particles, momentum and energy must be transferred.
cc This is later added to the correct B2.5 species using the "fluids_list".

        IF (NCOND.GT.0) THEN
          DO ICOND=1,NCOND
            DO IN=1,NCIP
!! particles
              SNI(IN,CONDENSE_TO(ICOND),ISTRAI)=
     .        SNI(IN,CONDENSE_TO(ICOND),ISTRAI) -
     .        PGENA(CONDENSE_FROM(ICOND),IN)*VOLTAL(IN)*FLX_EIR
!! momentum
              DO ICO = 0, 2
                KK = ICO*NFL
                SMO(IN,KK+CONDENSE_TO(ICOND),ISTRAI)=
     .          SMO(IN,KK+CONDENSE_TO(ICOND),ISTRAI) -
     .          VGENA(CONDENSE_FROM(ICOND),IN)*
     .                          VOLTAL(IN)*1.D-5*FLX_EIR
              END DO
!! energy
              SEI(IN,CONDENSE_TO(ICOND),ISTRAI) =
     .        SEI(IN,CONDENSE_TO(ICOND),ISTRAI) -
     .        EGENA(CONDENSE_FROM(ICOND),IN)*
     .                       VOLTAL(IN)*ELCHA
            ENDDO
          ENDDO
        ENDIF

cdr  fill bulk particle source rate sni(...ifl) from all contributing
cdr  test particle sources papl,pmpl,pipl,pppl (...,ipls)
cdr  test particle may have scored on a finer mesh (lcoarse=.false)

        DO IFL=1,NFLA

          CHPS(IFL)=0.
          SNIS(IFL)=0.
          CHMOS(IFL)=0.
          DO ICO=0,2
            KKK=ICO*NFL+IFL
            SMOS(KKK)=0.
          END DO

          DO JPLS=1,NPLSI
            IF (IFLB(JPLS).NE.IFL) CYCLE
            IPLSV=MPLSV(JPLS)
c  jpls contributes to plasma code species ifl
            DO IN=1,NCIP
              SNICL=0._DP
              IF (LPAPL)
     .         SNICL=SNICL+PAPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              IF (LPMPL)
     .         SNICL=SNICL+PMPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              IF (LPIPL)
     .         SNICL=SNICL+PIPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              SNICL=SNICL+PPPL_COP(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              SNI(IN,IFL,ISTRAI)=SNI(IN,IFL,ISTRAI)+SNICL
cdjm Jan2017
              IF (LPAPL) SNI_PAPL(IN,IFL,ISTRAI)=
     .                   SNI_PAPL(IN,IFL,ISTRAI)+
     .                                PAPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              IF (LPMPL) SNI_PMPL(IN,IFL,ISTRAI)=
     .                   SNI_PMPL(IN,IFL,ISTRAI)+
     .                                PMPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              IF (LPIPL) SNI_PIPL(IN,IFL,ISTRAI)=
     .                   SNI_PIPL(IN,IFL,ISTRAI)+
     .                                PIPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
              SNI_PPPL(IN,IFL,ISTRAI)=SNI_PPPL(IN,IFL,ISTRAI)+
     .                                PPPL_COP(JPLS,IN)*VOLTAL(IN)*
     .                                FLX_EIR

cdr:  totals, these now include pppl contribution (from pppl_cop)
              SNIS(IFL)=SNIS(IFL)+ SNICL
              CHPS(IFL)=CHPS(IFL)+CHPM(JPLS,IN)*VOLTAL(IN)
            ENDDO  ! IN LOOP

            IF (.NOT.LSHORT) THEN

cdr:  try to find stat. variance for particle balance sources
cdr   ntalm is a copv tally.
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == NPLSI+JPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IN=1,NCIP
                  SNIRES=0._DP
                  IF (LPAPL)
     .             SNIRES=SNIRES+PAPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
                  IF (LPMPL)
     .             SNIRES=SNIRES+PMPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
                  IF (LPIPL)
     .             SNIRES=SNIRES+PIPL(JPLS,IN)*VOLTAL(IN)*FLX_EIR
                  RESSNI(ISTRAI,IFL)=RESSNI(ISTRAI,IFL)+
     .                               ABS(SIGMA(ISTAT_COP,IN)*
     .                               SNIRES/100.D0)
                END DO
              end if

            END IF

cdr   particle sources done.

cdr   next: dwell on momentum sources. still inside IFL and JPLS loop
cdr   JPLS contributes to plasma code species IFL

            DO INC=1,NCIP
              IT = INVNCLTAL(INC)
              IF (IT == 0) CYCLE
              SMOCL_VEC(0:2)=0._DP
              DO ICO=0,2
cdr  EIRENE species indexing, for vectorial tallies
                KK=ICO*NPLS
cdr  B2.5 species indexing
                KKK=ICO*NFL+IFL

                IF (LMAPL_VEC)
     .           SMOCL_VEC(ico)=SMOCL_VEC(ico)+
     .            MAPL_VEC(kk+JPLS,INC)*VOLTAL(INC)*1.D-5*FLX_EIR
                IF (LMMPL_VEC)
     .           SMOCL_VEC(ico)=SMOCL_VEC(ico)+
     .            MMPL_VEC(kk+JPLS,INC)*VOLTAL(INC)*1.D-5*FLX_EIR
                IF (LMIPL_VEC)
     .           SMOCL_VEC(ico)=SMOCL_VEC(ico)+
     .            MIPL_VEC(kk+JPLS,INC)*VOLTAL(INC)*1.D-5*FLX_EIR
cdr  sampled tally MPPL does not exist.
cdr  MPPL_COP is analytic tally, and is for vol.rec. only.
                SMOCL_VEC(ico)=SMOCL_VEC(ico)+
     .           MPPL_COP(kk+JPLS,INC)*VOLTAL(INC)*1.D-5*FLX_EIR
cdr  now use B2.5 species indexing: 1:3*NFLA
                SMO(INC,kkk,ISTRAI)=SMO(INC,kkk,ISTRAI)+
     .           SMOCL_VEC(ico)
cdr:  totals, these now include mppl contribution (from mppl_cop)
                SMOS(kkk)=SMOS(kkk)+SMOCL_VEC(ico)
              ENDDO        ! ICO loop over three vector components
                           ! of momentum sources
              CHMOS(IFL)=CHMOS(IFL)+CHMOM(JPLS,INC)*VOLTAL(INC)
cdjm Jan2017
cdr  careful: IFL is B2.5 species index, JPLS is EIRENE species index
              IF (LMAPL_VEC) SMO_MAPL(INC,IFL,ISTRAI)=
     .               SMO_MAPL(INC,IFL,ISTRAI)+MAPL_VEC(JPLS,INC)*
     .                VOLTAL(INC)*1.D-5*FLX_EIR
              IF (LMMPL_VEC) SMO_MMPL(INC,IFL,ISTRAI)=
     .               SMO_MMPL(INC,IFL,ISTRAI)+MMPL_VEC(JPLS,INC)*
     .                VOLTAL(INC)*1.D-5*FLX_EIR
              IF (LMIPL_VEC) SMO_MIPL(INC,IFL,ISTRAI)=
     .               SMO_MIPL(INC,IFL,ISTRAI)+MIPL_VEC(JPLS,INC)*
     .                VOLTAL(INC)*1.D-5*FLX_EIR
cdr  sampled tally MPPL does not exist.
cdr  MPPL_COP is analytic tally, and is for vol.rec. only.
              SMO_MPPL(INC,IFL,ISTRAI)=
     .               SMO_MPPL(INC,IFL,ISTRAI)+MPPL_COP(JPLS,INC)*
     .                VOLTAL(INC)*1.D-5*FLX_EIR
cdr:  totals, these now include MPPL contribution (from mppl_cop)
            ENDDO  ! INC LOOP

            IF (.NOT.LSHORT) THEN

cdr:  try to find stat. variance for momentum balance sources
cdr   ntalm is a copv tally.
              istat_cop = 0
              do i = 1, nsigvi
                if ((iih(i) == ntalm).and.(igh(i) == 2*NPLSI+JPLS)) then
                  istat_cop = i
                  exit
                end if
              end do

              if (istat_cop > 0) then
                DO IN=1,NCIP
                  IT = INVNCLTAL(IN)
                  IF (IT == 0) CYCLE
                  SMORES=0._DP
                  IF (LMAPL_VEC) SMORES=SMORES+MAPL_VEC(JPLS,IN)*
     .                  VOLTAL(INC)*1.D-5*FLX_EIR
                  IF (LMMPL_VEC) SMORES=SMORES+MMPL_VEC(JPLS,IN)*
     .                  VOLTAL(INC)*1.D-5*FLX_EIR
                  IF (LMIPL_VEC) SMORES=SMORES+MIPL_VEC(JPLS,IN)*
     .                  VOLTAL(INC)*1.D-5*FLX_EIR

cdr  no need to add mppl_cop contribution here, because this is exact (MC noise free)
                  RESSMO(ISTRAI,IFL)=RESSMO(ISTRAI,IFL)+
     .                               ABS(SIGMA(ISTAT_COP,IN)*
     .                               SMORES/100.D0*1.D5)
                END DO
              end if

            END IF
          END DO ! JPLS LOOP
        END DO ! IFL LOOP
C
        DO JATM=1,NATM
          DO INC = 1, NCIP
            IT = INVNCLTAL(INC)
            IF (IT == 0) CYCLE
            IF (LPAAT) SNA_PAAT(INC,JATM,ISTRAI)=
     .               SNA_PAAT(INC,JATM,ISTRAI)+PAAT(JATM,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPMAT) SNA_PMAT(INC,JATM,ISTRAI)=
     .               SNA_PMAT(INC,JATM,ISTRAI)+PMAT(JATM,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPIAT) SNA_PIAT(INC,JATM,ISTRAI)=
     .               SNA_PIAT(INC,JATM,ISTRAI)+PIAT(JATM,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
          END DO
        END DO
C
        DO JMOL=1,NMOL
          DO INC = 1, NCIP
            IT = INVNCLTAL(INC)
            IF (IT == 0) CYCLE
            IF (LPAML) SNA_PAML(INC,JMOL,ISTRAI)=
     .               SNA_PAML(INC,JMOL,ISTRAI)+PAML(JMOL,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPMML) SNA_PMML(INC,JMOL,ISTRAI)=
     .               SNA_PMML(INC,JMOL,ISTRAI)+PMML(JMOL,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPIML) SNA_PIML(INC,JMOL,ISTRAI)=
     .               SNA_PIML(INC,JMOL,ISTRAI)+PIML(JMOL,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
          END DO
        END DO
C
        DO JION=1,NION
          DO INC = 1, NCIP
            IT = INVNCLTAL(INC)
            IF (IT == 0) CYCLE
            IF (LPAIO) SNA_PAIO(INC,JION,ISTRAI)=
     .               SNA_PAIO(INC,JION,ISTRAI)+PAIO(JION,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPMIO) SNA_PMIO(INC,JION,ISTRAI)=
     .               SNA_PMIO(INC,JION,ISTRAI)+PMIO(JION,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
            IF (LPIIO) SNA_PIIO(INC,JION,ISTRAI)=
     .               SNA_PIIO(INC,JION,ISTRAI)+PIIO(JION,INC)*
     .                                         VOLTAL(INC)*FLX_EIR
          END DO
        END DO
C
        CHEES=0.
        SEES=0.
        DO IN=1,NCIP
          IF (LEAEL) SEE(IN,ISTRAI)=SEE(IN,ISTRAI)+
     .            EAEL(IN)*VOLTAL(IN)*ELCHA
          IF (LEMEL) SEE(IN,ISTRAI)=SEE(IN,ISTRAI)+
     .            EMEL(IN)*VOLTAL(IN)*ELCHA
          IF (LEIEL) SEE(IN,ISTRAI)=SEE(IN,ISTRAI)+
     .            EIEL(IN)*VOLTAL(IN)*ELCHA
cdr  account for vol.rec electron energy loss
          SEE(IN,ISTRAI)=SEE(IN,ISTRAI)+
     .            EPEL_COP(IN)*VOLTAL(IN)*ELCHA
c
          IF (LPAEL) SNE_PAEL(IN,ISTRAI)=SNE_PAEL(IN,ISTRAI)+
     .                        PAEL(IN)*VOLTAL(IN)*FLX_EIR
          IF (LPMEL) SNE_PMEL(IN,ISTRAI)=SNE_PMEL(IN,ISTRAI)+
     .                        PMEL(IN)*VOLTAL(IN)*FLX_EIR
          IF (LPIEL) SNE_PIEL(IN,ISTRAI)=SNE_PIEL(IN,ISTRAI)+
     .                        PIEL(IN)*VOLTAL(IN)*FLX_EIR
          IF (LEAEL) SEE_EAEL(IN,ISTRAI)=SEE_EAEL(IN,ISTRAI)+
     .                        EAEL(IN)*VOLTAL(IN)*ELCHA
          IF (LEMEL) SEE_EMEL(IN,ISTRAI)=SEE_EMEL(IN,ISTRAI)+
     .                        EMEL(IN)*VOLTAL(IN)*ELCHA
          IF (LEIEL) SEE_EIEL(IN,ISTRAI)=SEE_EIEL(IN,ISTRAI)+
     .                        EIEL(IN)*VOLTAL(IN)*ELCHA
cdr  account for vol.rec electron energy loss
          SEE_EPEL(IN,ISTRAI)=SEE_EPEL(IN,ISTRAI)+
     .                        EPEL_COP(IN)*VOLTAL(IN)*ELCHA
cdr:  totals, these now include EPEL contribution (from epel_cop)
          CHEES=CHEES+CHEEM(IN)*VOLTAL(IN)
          IF (LEAEL) SEES=SEES+EAEL(IN)*VOLTAL(IN)
          IF (LEMEL) SEES=SEES+EMEL(IN)*VOLTAL(IN)
          IF (LEIEL) SEES=SEES+EIEL(IN)*VOLTAL(IN)
          SEES=SEES+EPEL_COP(IN)*VOLTAL(IN)
        ENDDO  ! IN LOOP
C
        SEIS=0.
        CHEIS=0.

        if (leapl) eapl0(1:ntrii) = sum(eapl(1:nplsi,1:ntrii),1)
        if (lempl) empl0(1:ntrii) = sum(empl(1:nplsi,1:ntrii),1)
        if (leipl) eipl0(1:ntrii) = sum(eipl(1:nplsi,1:ntrii),1)

        DO IFL=1,NFLA
          DO JPLS=1,NPLSI
            IF (IFLB(JPLS).NE.IFL) CYCLE
            DO IN=1,NCIP
              IF (LEAPL) SEI(IN,IFL,ISTRAI)=SEI(IN,IFL,ISTRAI)+
     .                  EAPL(JPLS,IN)*VOLTAL(IN)*ELCHA
              IF (LEMPL) SEI(IN,IFL,ISTRAI)=SEI(IN,IFL,ISTRAI)+
     .                  EMPL(JPLS,IN)*VOLTAL(IN)*ELCHA
              IF (LEIPL) SEI(IN,IFL,ISTRAI)=SEI(IN,IFL,ISTRAI)+
     .                  EIPL(JPLS,IN)*VOLTAL(IN)*ELCHA
cdr  account for vol.rec bulk ion energy loss
              SEI(IN,IFL,ISTRAI)=SEI(IN,IFL,ISTRAI)+
     .                  EPPL_COP(JPLS,IN)*VOLTAL(IN)*ELCHA
cdjm Jan2017
              IF (LEAPL) SEI_EAPL(IN,IFL,ISTRAI)=
     .            SEI_EAPL(IN,IFL,ISTRAI)+
     .            EAPL(JPLS,IN)*VOLTAL(IN)*ELCHA
              IF (LEMPL) SEI_EMPL(IN,IFL,ISTRAI)=
     .            SEI_EMPL(IN,IFL,ISTRAI)+
     .            EMPL(JPLS,IN)*VOLTAL(IN)*ELCHA
              IF (LEIPL) SEI_EIPL(IN,IFL,ISTRAI)=
     .            SEI_EIPL(IN,IFL,ISTRAI)+
     .            EIPL(JPLS,IN)*VOLTAL(IN)*ELCHA

              SEI_EPPL(IN,IFL,ISTRAI)=SEI_EPPL(IN,IFL,ISTRAI)+
     .            EPPL_COP(JPLS,IN)*VOLTAL(IN)*ELCHA

cdr:  totals, these now include EPPL contribution (from eppl_cop)
              CHEIS=CHEIS+CHEIM(INC)*VOLTAL(INC)
              IF (LEAPL) SEIS=SEIS+EAPL(JPLS,IN)*VOLTAL(IN)
              IF (LEMPL) SEIS=SEIS+EMPL(JPLS,IN)*VOLTAL(IN)
              IF (LEIPL) SEIS=SEIS+EIPL(JPLS,IN)*VOLTAL(IN)
cdr  account for vol.rec bulk ion energy loss
              SEIS=SEIS+EPPL_COP(JPLS,IN)*VOLTAL(IN)
            ENDDO ! IN LOOP

          ENDDO  ! JPLS LOOP
        ENDDO  ! IFL LOOP

        IF (.NOT.LSHORT) THEN

cdr:  try to find stat. variance for electron energy balance sources
cdr   ntalm is a copv tally.
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+1)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IN=1,NCIP
              SEERES=0._DP
              IF (LEAEL) SEERES=SEERES+EAEL(IN)*VOLTAL(IN)*FLX_EIR
              IF (LEMEL) SEERES=SEERES+EMEL(IN)*VOLTAL(IN)*FLX_EIR
              IF (LEIEL) SEERES=SEERES+EIEL(IN)*VOLTAL(IN)*FLX_EIR
              RESSEE(ISTRAI)=RESSEE(ISTRAI)+
     .                       ABS(SIGMA(ISTAT_COP,IN)*
     .                       SEERES/100.D0)
            END DO
          end if

cdr:  try to find stat. variance for ion energy balance sources
cdr   ntalm is a copv tally.
          istat_cop = 0
          do i = 1, nsigvi
            if ((iih(i) == ntalm).and.(igh(i) == 3*NPLSI+2)) then
              istat_cop = i
              exit
            end if
          end do

          if (istat_cop > 0) then
            DO IN=1,NCIP
              eamisum = 0._dp
              DO IFL=1,NFLA
                DO JPLS=1,NPLSI
                  IF (IFLB(JPLS).EQ.IFL) THEN
                    IF (LEAPL) eamisum = eamisum + EAPL(JPLS,IN)
                    IF (LEMPL) eamisum = eamisum + EMPL(JPLS,IN)
                    IF (LEIPL) eamisum = eamisum + EIPL(JPLS,IN)
                  END IF
                END DO  !NPLSI
              END DO  !nfla
              SEIRES=EAMISUM*VOLTAL(IN)*FLX_EIR
              RESSEI(ISTRAI)=RESSEI(ISTRAI)+
     .                       ABS(SIGMA(ISTAT_COP,IN)*
     .                       SEIRES/100.D0)
            END DO ! IN
          end if   !ISTAT_COP
        END IF  !lshort

CVK CALCULATES TOTAL POWER FOR SCALING
        IF (.NOT.LSHORT) THEN
         PRESS00=0._DP
         DO IT=1,NTRII                  !VK
          IC=IXTRI(IT)                !VK
          IF(IC.GT.0) THEN  !VK
            IN=NCLTAL(IT)
            IF (LEDENA) PRESS00=PRESS00+
     +              SUM(EDENA(:,IN))*VOLTAL(IN)*FLX_EIR
            IF (LEDENM) PRESS00=PRESS00+
     +              SUM(EDENM(:,IN))*VOLTAL(IN)*FLX_EIR
            IF (LEDENI) PRESS00=PRESS00+
     +              SUM(EDENI(:,IN))*VOLTAL(IN)*FLX_EIR
          END IF
         END DO
         PRESS00=PRESS00*ELCHA

CVK SCALING TO MAKE POSSIBLE COMPARISON WITH B2
         IF(PRESS00.GT.EPS12) THEN
          RESSEI(ISTRAI)=RESSEI(ISTRAI)/PRESS00
          RESSEE(ISTRAI)=RESSEE(ISTRAI)/PRESS00
         END IF

        END IF

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'IF3COP: TOTAL SOURCES P,EI,EE,MOM_VEC: ISTRA='
     .                   ,ISTRAI
        WRITE (iunout,'(6E14.6)') sum(SNIS(1:nfla)), SEIS, SEES,
     .                            sum(SMOS(1:nfla)),
     .                            sum(SMOS(1*nfl+1:1*nfl+nfla)),
     .                            sum(SMOS(2*nfl+1:2*nfl+nfla))
C
C   NEXT:
C   IF LSHORT OR IFIRST.GT.0     : CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT OR IFIRST.EQ.0 : ONLY RESCALE SURFACE SOURCE STRATA
C                  UNITS: # PER UNIT TARGET PLATE FLUX
C
        IF (IFIRST.EQ.0) THEN
C
          SNIS0(ISTRAI,0)=0.
          SMOS0(ISTRAI,-3:0)=0.
          DO 7550 IFL=1,NFLA
            SNIS0(ISTRAI,0)=SNIS0(ISTRAI,0)+SNIS(IFL)*FLXI
            SNIS0(ISTRAI,IFL)=SNIS(IFL)*FLXI

            SMOS0(ISTRAI,-1)=SMOS0(ISTRAI,-1)+SMOS(IFL)*FLXI
            SMOS0(ISTRAI,IFL)=SMOS(IFL)*FLXI
            SMOS0(ISTRAI,-2)=SMOS0(ISTRAI,-2)+SMOS(NFL+IFL)*FLXI
            SMOS0(ISTRAI,NFL+IFL)=SMOS(NFL+IFL)*FLXI
            SMOS0(ISTRAI,-3)=SMOS0(ISTRAI,-3)+SMOS(2*NFL+IFL)*FLXI
            SMOS0(ISTRAI,2*NFL+IFL)=SMOS(2*NFL+IFL)*FLXI
 7550     CONTINUE
          SEES0(ISTRAI)=SEES*FLXI
          SEIS0(ISTRAI)=SEIS*FLXI
C
        ELSEIF (LSHORT.AND.IFIRST.GT.0) THEN
C
          SNIS(0)=0.
          DO 7551 IFL=1,NFLA
            SNIS(0)=SNIS(0)+SNIS(IFL)
            SCALN(IFL)=SNIS0(ISTRAI,IFL)/(SNIS(IFL)+EPS60)
 7551     CONTINUE

          SCALN(0)=SNIS0(ISTRAI,0)/(SNIS(0)+EPS60)
C
          DO IC=1,NCIP
            SEE(IC,ISTRAI)=SEE(IC,ISTRAI)*SCALN(0)
cdjm Jan2017
            SNE_PAEL(IC,ISTRAI)=SNE_PAEL(IC,ISTRAI)*SCALN(0)
            SNE_PMEL(IC,ISTRAI)=SNE_PMEL(IC,ISTRAI)*SCALN(0)
            SNE_PIEL(IC,ISTRAI)=SNE_PIEL(IC,ISTRAI)*SCALN(0)
cdjm Jan2017
            SEE_EAEL(IC,ISTRAI)=SEE_EAEL(IC,ISTRAI)*SCALN(0)
            SEE_EMEL(IC,ISTRAI)=SEE_EMEL(IC,ISTRAI)*SCALN(0)
            SEE_EIEL(IC,ISTRAI)=SEE_EIEL(IC,ISTRAI)*SCALN(0)
            SEE_EPEL(IC,ISTRAI)=SEE_EPEL(IC,ISTRAI)*SCALN(0)
          END DO

          DO IFL=1,NFLA
            DO IC=1,NCIP
              SEI(IC,IFL,ISTRAI)=SEI(IC,IFL,ISTRAI)*SCALN(0)
cdjm Jan2017
              SEI_EAPL(IC,IFL,ISTRAI)=SEI_EAPL(IC,IFL,ISTRAI)*SCALN(0)
              SEI_EMPL(IC,IFL,ISTRAI)=SEI_EMPL(IC,IFL,ISTRAI)*SCALN(0)
              SEI_EIPL(IC,IFL,ISTRAI)=SEI_EIPL(IC,IFL,ISTRAI)*SCALN(0)
              SEI_EPPL(IC,IFL,ISTRAI)=SEI_EPPL(IC,IFL,ISTRAI)*SCALN(0)
            END DO  ! IC
          END DO

          DO IFL=1,NFLA
            DO IC=1,NCIP
              SNI(IC,IFL,ISTRAI)=SNI(IC,IFL,ISTRAI)*SCALN(IFL)
cdjm Jan2017
              SNI_PAPL(IC,IFL,ISTRAI)=
     .                 SNI_PAPL(IC,IFL,ISTRAI)*SCALN(IFL)
              SNI_PMPL(IC,IFL,ISTRAI)=
     .                 SNI_PMPL(IC,IFL,ISTRAI)*SCALN(IFL)
              SNI_PIPL(IC,IFL,ISTRAI)=
     .                 SNI_PIPL(IC,IFL,ISTRAI)*SCALN(IFL)
              SNI_PPPL(IC,IFL,ISTRAI)=
     .                 SNI_PPPL(IC,IFL,ISTRAI)*SCALN(IFL)

              KKK=IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*SCALN(IFL)
              KKK=NFL+IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*SCALN(IFL)
              KKK=2*NFL+IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*SCALN(IFL)
cdjm Jan2017
              SMO_MAPL(IC,IFL,ISTRAI)=
     .                 SMO_MAPL(IC,IFL,ISTRAI)*SCALN(IFL)
              SMO_MMPL(IC,IFL,ISTRAI)=
     .                 SMO_MMPL(IC,IFL,ISTRAI)*SCALN(IFL)
              SMO_MIPL(IC,IFL,ISTRAI)=
     .                 SMO_MIPL(IC,IFL,ISTRAI)*SCALN(IFL)
!dr  smo_cppv --> smo_MPPL
              SMO_MPPL(IC,IFL,ISTRAI)=
     .                 SMO_MPPL(IC,IFL,ISTRAI)*SCALN(IFL)
            END DO  ! IC
          END DO
          DO JATM=1,NATM
            IFL=1
            DO WHILE (LKINDP(IFL).NE.JATM.AND.IFL.LT.NFLA)
              IFL=IFL+1
            END DO
            DO IC=1,NCIP
              SNA_PAAT(IC,JATM,ISTRAI)=
     .         SNA_PAAT(IC,JATM,ISTRAI)*SCALN(IFL)
              SNA_PMAT(IC,JATM,ISTRAI)=
     .         SNA_PMAT(IC,JATM,ISTRAI)*SCALN(IFL)
              SNA_PIAT(IC,JATM,ISTRAI)=
     .         SNA_PIAT(IC,JATM,ISTRAI)*SCALN(IFL)
            END DO
          END DO
          DO JMOL=1,NMOL
            IFL=1
            DO WHILE (LKINDP(IFL).NE.LKINDM(JMOL).AND.IFL.LT.NFLA)
              IFL=IFL+1
            END DO
            DO IC=1,NCIP
              SNA_PAML(IC,JMOL,ISTRAI)=
     .         SNA_PAML(IC,JMOL,ISTRAI)*SCALN(IFL)
              SNA_PMML(IC,JMOL,ISTRAI)=
     .         SNA_PMML(IC,JMOL,ISTRAI)*SCALN(IFL)
              SNA_PIML(IC,JMOL,ISTRAI)=
     .         SNA_PIML(IC,JMOL,ISTRAI)*SCALN(IFL)
            END DO
          END DO
          DO JION=1,NION
            IFL=1
            DO WHILE (LKINDP(IFL).NE.LKINDI(JION).AND.IFL.LT.NFLA)
              IFL=IFL+1
            END DO
            DO IC=1,NCIP
              SNA_PAIO(IC,JION,ISTRAI)=
     .         SNA_PAIO(IC,JION,ISTRAI)*SCALN(IFL)
              SNA_PMIO(IC,JION,ISTRAI)=
     .         SNA_PMIO(IC,JION,ISTRAI)*SCALN(IFL)
              SNA_PIIO(IC,JION,ISTRAI)=
     .         SNA_PIIO(IC,JION,ISTRAI)*SCALN(IFL)
            END DO
          END DO
C
          LTEST=.TRUE.  ! tactically assume:
                        ! this stratum will continue in short cycle mode

          IF (LSTOP)
     .      WRITE (iunout,*) 'STOP SHORT CYCLE: ALL B2 TIMESTEPS DONE'

!  DO AT LEAST ONE MORE TIME STEP
CDR  DECIDE FOR THIS CURRENT STRATUM ISTRAI:
CDR  SHORT CYCLE (IMPLICIT CORRECTION) ONLY, OR FULL MONTE CARLO

          DO 7558 IFL=1,NFLA
            TEST=CHPS(IFL)/(SNIS(IFL)+1.D-60)*100.
            write (iunout,*) ' global change in sni,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGP) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE.  ! stop short cycle mode.
                             ! Full new set of trajectories.
              WRITE (iunout,'(A,1P,3E12.4)')
     .                    'STOP SHORT CYCLE: PART. SOURCES: ',
     .                     SNIS(IFL),CHPS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                          ISTRAI,IFL
            ENDIF
            TEST=CHMOS(IFL)/(SMOS(IFL)+1.D-60)*100.
            write (iunout,*) ' global change in smo,ifl ',test,ifl
            IF (ABS(TEST).GT.CHGMOM) THEN
              LSTP3=.TRUE.
              LTEST=.FALSE. ! stop short cycle mode.
                            ! Full new set of trajectories.
              WRITE (iunout,'(A,1P,3E12.4)')
     .                    'STOP SHORT CYCLE: MOMENTUM SOURCE: ',
     .                     SMOS(IFL),CHMOS(IFL),TEST
              WRITE (iunout,*) 'STRATUM ISTRAI, SPECIES IFL ',
     .                          ISTRAI,IFL
            ENDIF
 7558     CONTINUE

          TEST=CHEES/(SEES+1.D-60)*100.
          write (iunout,*) ' global change in see ',test
          IF (ABS(TEST).GT.CHGEE) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE. ! stop short cycle mode.
                          ! Full new set of trajectories.
            WRITE (iunout,'(A,1P,3E12.4)')
     .                  'STOP SHORT CYCLE: EL. EN. SOURCE: ',
     .                   SEES,CHEES,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          TEST=CHEIS/(SEIS+1.D-60)*100.
          write (iunout,*) ' global change in sei ',test
          IF (ABS(TEST).GT.CHGEI) THEN
            LSTP3=.TRUE.
            LTEST=.FALSE. ! stop short cycle mode.
                          ! Full new set of trajectories.
            WRITE (iunout,'(A,1P,3E12.4)')
     .                  'STOP SHORT CYCLE: ION EN. SOURCE: ',
     .                   SEIS,CHEIS,TEST
            WRITE (iunout,*) 'STRATUM ISTRAI ',ISTRAI
          ENDIF
          IF (LSHORT) LSTP=LSTP3

        ELSE
cdr  here: .not.lshort .and. ifirst.gt.0
        ENDIF

csw ignore short cycle attempt
        LSTP3=.true.
        LSTOP=LSTP3
        LTEST=.false.
        if(lshort) then
          LSTP=LSTP3
        end if
csw
        NLSRON(ISTRAI) = .NOT.LTEST
C
        IF (.NOT.LSHORT) THEN
C
          DO IC=1,NCIP
            SEE(IC,ISTRAI)=SEE(IC,ISTRAI)*FLXI
cdjm Jan2017
            SNE_PAEL(IC,ISTRAI)=SNE_PAEL(IC,ISTRAI)*FLXI
            SNE_PMEL(IC,ISTRAI)=SNE_PMEL(IC,ISTRAI)*FLXI
            SNE_PIEL(IC,ISTRAI)=SNE_PIEL(IC,ISTRAI)*FLXI
cdjm Jan2017
            SEE_EAEL(IC,ISTRAI)=SEE_EAEL(IC,ISTRAI)*FLXI
            SEE_EMEL(IC,ISTRAI)=SEE_EMEL(IC,ISTRAI)*FLXI
            SEE_EIEL(IC,ISTRAI)=SEE_EIEL(IC,ISTRAI)*FLXI
            SEE_EPEL(IC,ISTRAI)=SEE_EPEL(IC,ISTRAI)*FLXI
          END DO
          DO IFL=1,NFLA
            DO IC=1,NCIP
              SEI(IC,IFL,ISTRAI)=SEI(IC,IFL,ISTRAI)*FLXI
cdjm Jan2017
              SEI_EAPL(IC,IFL,ISTRAI)=SEI_EAPL(IC,IFL,ISTRAI)*FLXI
              SEI_EMPL(IC,IFL,ISTRAI)=SEI_EMPL(IC,IFL,ISTRAI)*FLXI
              SEI_EIPL(IC,IFL,ISTRAI)=SEI_EIPL(IC,IFL,ISTRAI)*FLXI
              SEI_EPPL(IC,IFL,ISTRAI)=SEI_EPPL(IC,IFL,ISTRAI)*FLXI
            END DO  ! IC
          END DO  ! IFL
          DO IFL=1,NFLA
            DO IC=1,NCIP
              SNI(IC,IFL,ISTRAI)=SNI(IC,IFL,ISTRAI)*FLXI
cdjm Jan2017
              SNI_PAPL(IC,IFL,ISTRAI)=
     .               SNI_PAPL(IC,IFL,ISTRAI)*FLXI
              SNI_PMPL(IC,IFL,ISTRAI)=
     .               SNI_PMPL(IC,IFL,ISTRAI)*FLXI
              SNI_PIPL(IC,IFL,ISTRAI)=
     .               SNI_PIPL(IC,IFL,ISTRAI)*FLXI
              SNI_PPPL(IC,IFL,ISTRAI)=
     .               SNI_PPPL(IC,IFL,ISTRAI)*FLXI

              KKK=IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*FLXI
              KKK=NFL+IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*FLXI
              KKK=2*NFL+IFL
              SMO(IC,KKK,ISTRAI)=SMO(IC,KKK,ISTRAI)*FLXI
cdjm Jan2017
              SMO_MAPL(IC,IFL,ISTRAI)=
     .               SMO_MAPL(IC,IFL,ISTRAI)*FLXI
              SMO_MMPL(IC,IFL,ISTRAI)=
     .               SMO_MMPL(IC,IFL,ISTRAI)*FLXI
              SMO_MIPL(IC,IFL,ISTRAI)=
     .               SMO_MIPL(IC,IFL,ISTRAI)*FLXI
              SMO_MPPL(IC,IFL,ISTRAI)=
     .               SMO_MPPL(IC,IFL,ISTRAI)*FLXI
            END DO  ! IC
          END DO  ! IFL
          DO JATM=1,NATM
            DO IC=1,NCIP
              SNA_PAAT(IC,JATM,ISTRAI)=
     .             SNA_PAAT(IC,JATM,ISTRAI)*FLXI
              SNA_PMAT(IC,JATM,ISTRAI)=
     .             SNA_PMAT(IC,JATM,ISTRAI)*FLXI
              SNA_PIAT(IC,JATM,ISTRAI)=
     .             SNA_PIAT(IC,JATM,ISTRAI)*FLXI
            END DO
          END DO
          DO JMOL=1,NMOL
            DO IC=1,NCIP
              SNA_PAML(IC,JMOL,ISTRAI)=
     .             SNA_PAML(IC,JMOL,ISTRAI)*FLXI
              SNA_PMML(IC,JMOL,ISTRAI)=
     .             SNA_PMML(IC,JMOL,ISTRAI)*FLXI
              SNA_PIML(IC,JMOL,ISTRAI)=
     .             SNA_PIML(IC,JMOL,ISTRAI)*FLXI
            END DO
          END DO
          DO JION=1,NION
            DO IC=1,NCIP
              SNA_PAIO(IC,JION,ISTRAI)=
     .             SNA_PAIO(IC,JION,ISTRAI)*FLXI
              SNA_PMIO(IC,JION,ISTRAI)=
     .             SNA_PMIO(IC,JION,ISTRAI)*FLXI
              SNA_PIIO(IC,JION,ISTRAI)=
     .             SNA_PIIO(IC,JION,ISTRAI)*FLXI
            END DO
          END DO
C
        ENDIF
csw 26jan2011 extra B25, correct particle sources for recycling strata if necessary
        if(my_pe == 0) then
          call eirene_extrab25_eirpbls(istrai)
          srcstrn(istrai)=flux(istrai)
        endif
csw
C
 7999 CONTINUE
C
C  DATA TRANSFER BACK TO PLASMA CODE FINISHED FOR STRATUM NO. ISTRAI
C
csw 14jul2011
        if(allocated(save_estimv)) then
          estimv = save_estimv
          deallocate(save_estimv)
        endif
        if(allocated(save_estims)) then
          estims = save_estims
          deallocate(save_estims)
        endif
csw
10000 CONTINUE
C
      RETURN
C
      END SUBROUTINE EIRENE_IF3COP
C
      SUBROUTINE EIRENE_IF4COP
      USE EIRMOD_REFUSR, ONLY : EIRENE_WRITE_CONBE
      IMPLICIT NONE
      INTEGER, SAVE :: MINSPEZ, MAXSPEZ, NREC11
      REAL(DP) :: SFNIT(0:NSTEP,NFL), SFEIT(0:NSTEP),
     R            SFEET(0:NSTEP), SHEAE(0:NSTEP), SHEAI(0:NSTEP)
      REAL(DP) :: SSNI(NFL)
      REAL(DP), SAVE :: TOT, TOTI, TOTE, FLX,
     R            SSE, SSI, SSEI, SSEE, RE, RI, RNT,
     R            BALAN, BALANI, BALANE
      REAL(DP), ALLOCATABLE :: SFETMP(:)
cdr for species-dependent global particle balance
      REAL(DP) :: SSN(NFL), BALANN(NFL), TOTN(NFL), RN(NFL), FN
cdr  sputter fluxes
      REAL(DP) :: SPAT(0:NATM,0:NSTRA), SPML(0:NMOL,0:NSTRA),
     .            SPIO(0:NION,0:NSTRA), SPPL(0:NPLS,0:NSTRA)
      INTEGER :: I, J, JC, IRC, IFL, ITARG, IC,
     .           IT, IFC, ITRI, IS, JTRI, JS, JFC, IL, ISTR
#ifdef WG_TODO
      INTEGER :: IPRT, IX, IY
      REAL(DP) :: SFNISY(NFL), SFNINY(NFL), SFNIWX(NFL), SFNIEX(NFL)
      REAL(DP), SAVE :: TIFLX,
     R            SFEIWX, SFEEWX, SFEIEX, SFEEEX,
     R            SFEENY, SFEINY, SFEESY, SFEISY
#endif
      INTEGER, ALLOCATABLE, SAVE :: IHELP(:)
      real(dp), allocatable :: helpw(:)
      integer, allocatable, save :: nonrec(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: SFEIN(:), SFEEN(:), SFNIN(:,:)
      CHARACTER(24) :: CFORM
      TYPE(FACE), POINTER :: CURFC
      EXTERNAL :: EIRENE_LEER, EIRENE_MASRR1, EIRENE_MASR1,
     .            EIRENE_MASR2, EIRENE_MASR3, EIRENE_MASR4,
     .            EIRENE_MASJ2, EIRENE_MASJ4, EIRENE_MASJ5,
     .            EIRENE_MASYR1, EIRENE_EXIT_OWN
C
      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11,
     .      ERR=1987)
      IRC=3
C  WRITE RCCPL
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  RCCPL,   IRC= ',IRC
C     IRC=3   STILL
C  WRITE ICCPL1
      ALLOCATE (IHELP(NOUTAU))
      JC=0
      DO J=1,6*NSTEP
        JC=JC+1
        IHELP(JC)=ICCPL1(J)
        IF (JC == NOUTAU) THEN
          IRC=IRC+1
          WRITE (11,REC=IRC) IHELP
          IF (TRCINT.OR.TRCFLE)
     .        WRITE (iunout,*) 'WRITE 11  ICCPL1,  IRC= ',IRC
          JC=0
        END IF
      END DO
c  write last (incomplete) record of ICCPL1
      IF (JC > 0) THEN
        IHELP(JC+1:NOUTAU) = 0   ! fill up last record,
                                 ! up to full length NOUTAU
        IRC=IRC+1
        WRITE (11,REC=IRC) IHELP
        IF (TRCINT.OR.TRCFLE)
     .      WRITE (iunout,*) 'WRITE 11  ICCPL1,  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
C  WRITE ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  ICCPL2,  IRC= ',IRC
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
      IF (TRCINT.OR.TRCFLE)
     .    WRITE (iunout,*) 'WRITE 11  LCCPL,   IRC= ',IRC
C
      IF (LSHORT) LSTOP=LSTP3
C
      IF (.NOT.LSTOP) RETURN
C
      IF (.NOT.(LBALAN)) GOTO 11000
C
C  BALANCES, SHOULD BE DONE ONLY AT THE END OF B2.5 RUN
C  AT THE END OF AN EIRENE RUN THE BALANCES MAY BE OFF AT LEAST AT
C  THE BEGINNING OF THE CYCLING PROCEDURE, BECAUSE THE PLASMA STILL
C  HAS TO ADJUST TO THE NEW SOURCES
C
C
C...............................................................................
C
C  CHECK FLUXES AT THE GRID BOUNDARY, SIGN, RECYCLING OR NON-RECYCLING BOUNDARY
C  COUNT FLUXES FROM OUTSIDE INTO GRID AS POSITIVE
C
C...............................................................................
      IF (.NOT.ALLOCATED(NONREC))
     > ALLOCATE (NONREC(size(INMTI,1),size(INMTI,2)))
      NONREC = INMTI

C  FIRST: FLUXES TO THOSE SURFACES, AT WHICH RECYCLING BOUNDARY
C         CONDITIONS ARE SPECIFIED
      SFEIT(0)=0.
      SFEET(0)=0.
      SFNIT=0.
      SHEAE(0)=0.
      SHEAI(0)=0.

      DO I=1,NTARGI
        SFEIT(I)=0.
        SFEET(I)=0.
        SHEAE(I)=0.
        SHEAI(I)=0.
C  BALANCE CONTRIB. X-GRID REC. SOURCE
        DO IT=1,TRGT(I)%NTRGDAT
          IFC = TRGT(I)%FACES(IT)
          CURFC => NCLFC(IFC)%FP
          IF (.NOT.ASSOCIATED(CURFC)) CYCLE
          DO WHILE (ASSOCIATED(CURFC))
            ITRI = CURFC%TRIANGLE
            IS = CURFC%SIDE
            IC = NCLTAL(ITRI) ! all triangles should point
                              ! to the same plasma cell
            NONREC(IS,ITRI) = 0
            CURFC => CURFC%NEXT
          END DO
          DO IFL=NSPZI(I),NSPZE(I)
            FN = TRGT(I)%FCORI(IT)*FNIB(IFC,IFL)
            IF (FN.GE.0.d0) THEN
              SFNIT(I,IFL)=SFNIT(I,IFL)-FN
!              SHEAE(I)=SHEAE(I)+TEB(IC)*FN*(-DELTA_SHEATHB(IFC))
!              SHEAI(I)=SHEAI(I)+TEB(IC)*FN*  DELTA_SHEATHB(IFC)*elcha
              SHEAE(I)=SHEAE(I)+
     .                 TRGT(I)%TET(IT)*FN*(-DELTA_SHEATHB(IFC))
              SHEAI(I)=SHEAI(I)+
     .                 TRGT(I)%TET(IT)*FN*  DELTA_SHEATHB(IFC)*elcha
            ELSE
              WRITE (iunout,*)
     .              'WRONG ORIENTATION OF TARGET RECYCLING FLUX '
              CALL EIRENE_MASJ4 (
     .             'ITARG, IPLS, NDT, IFC             ',
     .              I    , IFL,  NDT(I), IFC)
              WRITE (iunout,*) 'FNI(IFC) ',FNIB(IFC,IFL)
              FN = 0._DP
              FNIB(IFC,IFL) = 0._DP
            ENDIF
          END DO ! IFL
          SFEIT(I)=SFEIT(I)-TRGT(I)%FCORI(IT)*FEIB(IFC)
          SFEET(I)=SFEET(I)-TRGT(I)%FCORI(IT)*FEEB(IFC)
        END DO ! IT
      END DO ! I
      SHEAE(1:NTARGI) = SHEAE(1:NTARGI) * ELCHA
!      SHEAI(1:NTARGI) = SHEAI(1:NTARGI) * ELCHA
      SFNIT = SFNIT * ELCHA

      IF (.NOT.ALLOCATED(SFEIN)) THEN
#ifndef LEGACYCOMP
        ALLOCATE(SFEIN(0:NLIM+NSTS),SOURCE=0._DP)
        ALLOCATE(SFEEN(0:NLIM+NSTS),SOURCE=0._DP)
        ALLOCATE(SFNIN(0:NLIM+NSTS,NFL),SOURCE=0._DP)
#else
        ALLOCATE(SFEIN(0:NLIM+NSTS))
        ALLOCATE(SFEEN(0:NLIM+NSTS))
        ALLOCATE(SFNIN(0:NLIM+NSTS,NFL))
        SFEIN = 0._DP
        SFEEN = 0._DP
        SFNIN = 0._DP
#endif
      END IF

      DO ITRI = 1, NTRII
        DO IS = 1, 3
          IT = NONREC(IS,ITRI)
          IF (IT /= 0) THEN
! LOCAL REFLECTION MODEL ENCOUNTERED ON TRIANGLE FACE
            IF (NCHBAR(IS,ITRI) == 0) THEN
! NO NEIGHBOR ==> DOMAIN BOUNDARY
              IFC = IFCTRI(IS,ITRI)
              IF (IFC /= 0) THEN
                SFEIN(IT)=SFEIN(IT)+FEIB(IFC)
                SFEEN(IT)=SFEEN(IT)+FEEB(IFC)
                DO IFL=1,NFLA
                  SFNIN(IT,IFL)=SFNIN(IT,IFL)+FNIB(IFC,IFL)
                END DO
              END IF
            ELSE
! CHECK IF TRIANGLE FACE IS ON EDGE OF THE PLASMA REGION
              JTRI = NCHBAR(IS,ITRI)
              JS = NSEITE(IS,ITRI)
!pb           IF (IXTRI(JTRI) == 0) THEN
              IF (IXTRI(JTRI) <= 0) THEN
                JFC = IFCTRI(JS,JTRI)
                IF (JFC /= 0) THEN
                  SFEIN(IT)=SFEIN(IT)+FEIB(JFC)
                  SFEEN(IT)=SFEEN(IT)+FEEB(JFC)
                  DO IFL=1,NFLA
                    SFNIN(IT,IFL)=SFNIN(IT,IFL)+FNIB(JFC,IFL)
                  END DO
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
      SFNIN = SFNIN * ELCHA

      WRITE (37,*) 'NON RECYCLING FLUXES '
      WRITE (37,'(a6,3a14)') 'SRF.NO','SFEI','SFEE','SFNI'
      WRITE (CFORM,'(a,i3,a)') '(I6,2ES14.6,',NFLA,'ES14.6)'
      DO I = 1, NLIM+NSTS
        IF (ABS(SFEIN(I)) < EPS30) CYCLE
        IL = I
        IF (IL > NLIM) IL = -IL + NLIM
        WRITE (37,CFORM) IL,SFEIN(I),SFEEN(I),(SFNIN(I,IFL),IFL=1,NFLA)
        SFEIN(0) = SFEIN(0) + SFEIN(I)
        SFEEN(0) = SFEEN(0) + SFEEN(I)
        SFNIN(0,:) = SFNIN(0,:) + SFNIN(I,:)
      END DO

      WRITE (CFORM,'(a,i3,a)') '(I6,4ES14.6,',NFLA,'ES14.6)'
      WRITE (37,'(a6,5a14)') 'SRF.NO','SFEI','SFEE',
     .                       'SHEAI','SHEAE','SFNI'
      WRITE (37,*) 'RECYCLING FLUXES '
      WRITE (37,'(a6,3a14)') 'SRF.NO','SFEI','SFEE','SFNI'
      DO I = 1, NTARGI
        IF (ABS(SFEIT(I)) < EPS30) CYCLE
        WRITE (37,CFORM) I,SFEIT(I),SFEET(I),SHEAI(I),SHEAE(I),
     .                    (SFNIT(I,IFL),IFL=1,NFLA)
        SFEIT(0)=SFEIT(0)+SFEIT(I)
        SFEET(0)=SFEET(0)+SFEET(I)
        SFNIT(0,:)=SFNIT(0,:)+SFNIT(I,:)
        SHEAE(0)=SHEAE(0)+SHEAE(I)
        SHEAI(0)=SHEAI(0)+SHEAI(I)
      END DO

#ifdef WG_TODO
C
C  FIRST: SOUTH EDGE: IY=0
C
C NON-RECYCLING FLUXES AT SOUTH EDGE: SFEISY,SFEESY,SFNISY
      SFEISY=0.
      SFEESY=0.
      SFNISY=0.
      NPBS=0
      NPES=NPBS+1
C
      DO 10113 IX=1,NDXA
C
C IS (IX,0) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.0) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10110
            ENDIF
          ENDDO
        ENDDO

        ITARG=0
        IF (LCUT(IX)) GOTO 10113
C
C SURFACE NORMAL IS INWARD. HENCE: TAKE ALL FLUXES F...YB POSITIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEISY=SFEISY+FEIYB(IX,0)
        SFEESY=SFEESY+FEEYB(IX,0)
        DO 10111 IFL=1,NFLA
          SI=SIGN(1._DP,FNIYB(IX,0,IFL))
          SFNISY(IFL)=SFNISY(IFL)+
     .                FNIYB(IX,0,IFL)+SI*ABS(FNIY_XB(IX,0,IFL))
10111   CONTINUE
        GOTO 10113
10110   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10112 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            SI=SIGN(1._DP,FNIYB(IX,0,IFL))
            DUMVAL = FNIYB(IX,0,IFL)+SI*ABS(FNIY_XB(IX,0,IFL))
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*)
     .          'SOUTH,IX,ITARG,IPRT,IFL ',IX,ITARG,IPRT,IFL
              SFNISY(IFL)=SFNISY(IFL)+DUMVAL
              TIFLX=TIFLX+FEIYB(IX,0)*DUMVAL
            ENDIF
10112     CONTINUE

          TIFLX=TIFLX/(FLX+EPS60)
          SFEISY=SFEISY+TIFLX
        ENDIF

10113 CONTINUE
C
      SFNISY=SFNISY*ELCHA
C
      WRITE (IUNOUT,*) 'NON-RECYCLING FLUXES FROM SOUTH EDGE'
      WRITE (IUNOUT,8888) sum(SFNISY(1:nfla)),SFEISY,SFEESY
C
C
C  SECOND: NORTH EDGE: IY=NDYA
C
C NON-RECYCLING FLUXES AT NORTH EDGE: SFEINY,SFEENY,SFNINY
      SFEINY=0.
      SFEENY=0.
      SFNINY=0.
      NPBS=NDYA
      NPES=NPBS+1
      DO 10118 IX=1,NDXA
C
C IS (IX,NDYA) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.2.AND.NDT(ITARG,IPRT).EQ.NDYA) THEN
              IF (IX.GE.NTIN(ITARG,IPRT).AND.
     .            IX.LT.NTEN(ITARG,IPRT)) GOTO 10115
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
        IF (LCUT(IX)) GOTO 10118
C
C SURFACE NORMAL IS OUTWARD. HENCE: TAKE ALL FLUXES F...YB NEGATIVE
C SIGN OF ADDITIONAL COMPONENT DUE TO INCLINED GRID AS SIGN OF F...YB
C
        SFEINY=SFEINY-FEIYB(IX,NDYA)
        SFEENY=SFEENY-FEEYB(IX,NDYA)
        DO 10116 IFL=1,NFLA
          SI=SIGN(1._DP,FNIYB(IX,NDYA,IFL))
          SFNINY(IFL)=SFNINY(IFL)-FNIYB(IX,NDYA,IFL)-
     .                SI*ABS(FNIY_XB(IX,NDYA,IFL))
10116   CONTINUE
        GOTO 10118
10115   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10117 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            SI=SIGN(1._DP,FNIYB(IX,NDYA,IFL))
            DUMVAL=FNIYB(IX,NDYA,IFL)+SI*ABS(FNIY_XB(IX,NDYA,IFL))
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*)
     .          'NORTH,IX,ITARG,IPRT,IFL ',IX,ITARG,IPRT,IFL
              SFNINY(IFL)=SFNINY(IFL)-DUMVAL
              TIFLX=TIFLX+FEIYB(IX,NDYA)*(-DUMVAL)
            ENDIF
10117     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEINY=SFEINY-TIFLX
        ENDIF

10118 CONTINUE
C
      SFNINY=SFNINY*ELCHA
C
      WRITE (IUNOUT,*) 'NON-RECYCLING FLUXES TO NORTH EDGE'
      WRITE (IUNOUT,8888) sum(SFNINY(1:nfla)),SFEINY,SFEENY
C
C
C  THIRD: WEST EDGE: IX=0
C
C NON-RECYCLING FLUXES AT WEST EDGE: SFEIWX,SFEEWX,SFNIWX
      SFEIWX=0.
      SFEEWX=0.
      SFNIWX=0.
      DO 10123 IY=1,NDYA
C
C IS (0,IY) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
             IF (NIXY(ITARG,IPRT).EQ.1.AND.NDT(ITARG,IPRT).EQ.0) THEN
               IF (IY.GE.NTIN(ITARG,IPRT).AND.
     .             IY.LT.NTEN(ITARG,IPRT)) GOTO 10120
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
        SFEIWX=SFEIWX+FEIXB(0,IY)
        SFEEWX=SFEEWX+FEEXB(0,IY)
        DO 10121 IFL=1,NFLA
          SI=SIGN(1._DP,FNIXB(0,IY,IFL))
          SFNIWX(IFL)=SFNIWX(IFL)+FNIXB(0,IY,IFL)+
     .                SI*ABS(FNIX_YB(0,IY,IFL))
10121   CONTINUE
10120   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10122 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            SI=SIGN(1._DP,FNIXB(0,IY,IFL))
            DUMVAL=FNIXB(0,IY,IFL)+SI*ABS(FNIX_YB(0,IY,IFL))
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'WEST,IY,ITARG,IPRT,IFL ',
     .                               IY,ITARG,IPRT,IFL
              SFNIWX(IFL)=SFNIWX(IFL)+DUMVAL
              TIFLX=TIFLX+FEIXB(0,IY)*DUMVAL
            ENDIF
10122     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEIWX=SFEIWX+TIFLX
        ENDIF
10123 CONTINUE
C
      SFNIWX=SFNIWX*ELCHA
C
      WRITE (IUNOUT,*) 'NON-RECYCLING FLUXES FROM WEST EDGE'
      WRITE (IUNOUT,8888) sum(SFNIWX(1:nfla)),SFEIWX,SFEEWX
C
C
C  FOURTH: EAST EDGE: IX=NDXA
C
C NON-RECYCLING FLUXES AT EAST EDGE: SFEIEX,SFEEEX,SFNIEX
      SFEIEX=0.
      SFEEEX=0.
      SFNIEX=0.
      DO 10128 IY=1,NDYA
C
C IS (NDXA,IY) A RECYCLING SOURCE? IF YES, DO NOT COUNT HERE
C
        DO ITARG=1,NTARGI
          DO IPRT=1,NTGPRT(ITARG)
            IF (NIXY(ITARG,IPRT).EQ.1.AND.NDT(ITARG,IPRT).EQ.NDXA) THEN
              IF (IY.GE.NTIN(ITARG,IPRT).AND.
     .            IY.LT.NTEN(ITARG,IPRT)) GOTO 10125
            ENDIF
          ENDDO
        ENDDO
        ITARG=0
C
        SFEIEX=SFEIEX-FEIXB(NDXA,IY)
        SFEEEX=SFEEEX-FEEXB(NDXA,IY)
        DO 10126 IFL=1,NFLA
          SI=SIGN(1._DP,FNIXB(NDXA,IY,IFL))
          SFNIEX(IFL)=SFNIEX(IFL)-FNIXB(NDXA,IY,IFL)-
     .                SI*ABS(FNIX_YB(NDXA,IY,IFL))
10126   CONTINUE

10125   CONTINUE

C  DO NOT RECYCLE TARGET FLUXES WITH FALSE ORIENTATION
        IF (ITARG.GT.0) THEN
C  ITARG, IPRT KNOWN FROM ABOVE
          FLX=0.
          TIFLX=0.
          DO 10127 IFL=NSPZI(ITARG,IPRT),NSPZE(ITARG,IPRT)
            SI=SIGN(1._DP,FNIXB(NDXA,IY,IFL))
            DUMVAL=FNIXB(NDXA,IY,IFL)+SI*ABS(FNIX_YB(NDXA,IY,IFL))
            FLX=FLX+DUMVAL
            IF (NINCT(ITARG,IPRT)*DUMVAL.LT.0) THEN
              WRITE (iunout,*)
     .          'RECYCLING TARGET, BUT WRONG FLOW DIRECTION: '
              WRITE (iunout,*) 'EAST,IY,ITARG,IPRT,IFL ',
     .                               IY,ITARG,IPRT,IFL
              SFNIEX(IFL)=SFNIEX(IFL)-DUMVAL
              TIFLX=TIFLX+FEIXB(NDXA,IY)*(-DUMVAL)
            ENDIF
10127     CONTINUE
          TIFLX=TIFLX/(FLX+EPS60)
          SFEIEX=SFEIEX-TIFLX
        ENDIF

10128 CONTINUE
C
      SFNIEX=SFNIEX*ELCHA
C
      WRITE (IUNOUT,*) 'NON-RECYCLING FLUXES TO EAST EDGE'
      WRITE (IUNOUT,8888) sum(SFNIEX(1:nfla)),SFEIEX,SFEEEX
C
C  NEXT: FLUXES TO THOSE SURFACES, AT WHICH RECYCLING BOUNDARY
C        CONDITIONS ARE SPECIFIED
      CALL EIRENE_LEER(2)
C
      SFEIT(0)=0.
      SFEET(0)=0.
      SFNIT(0,:)=0.
      SHEAE(0)=0.
      SHEAI(0)=0.
      DO 10139 I=1,NTARGI
        SFEIT(I)=0.
        SFEET(I)=0.
        SFNIT(I,:)=0.
        SHEAE(I)=0.
        SHEAI(I)=0.
        DO IPRT=1,NTGPRT(I)

          IF (NIXY(I,IPRT).EQ.1) THEN
C  BALANCE CONTRIB. X-GRID REC. SOURCE
            DO 10132 IY=NTIN(I,IPRT),NTEN(I,IPRT)-1
              DO 10131 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL).GE.0.D0)
     >           THEN
                  SFNIT(I,IFL)=SFNIT(I,IFL)-
     .              NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)-
     .              ABS(FNIX_YB(NDT(I,IPRT),IY,IFL))

cdr sheath contributions: count negative for electrons, positive for ions
cdr unfinished: need to account for charge state of ion species IFL
                  SHEAE(I)=SHEAE(I)+TEB(NDT(I,IPRT),IY)*
     .             (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)+
     .              ABS(FNIX_YB(NDT(I,IPRT),IY,IFL)))*
     .             (-DELTA_SHEATHXB(NDT(I,IPRT),IY))
                  SHEAI(I)=SHEAI(I)+TEB(NDT(I,IPRT),IY)*
     .             (NINCT(I,IPRT)*FNIXB(NDT(I,IPRT),IY,IFL)+
     .              ABS(FNIX_YB(NDT(I,IPRT),IY,IFL)))*
     .              DELTA_SHEATHXB(NDT(I,IPRT),IY)
cdr  sheath done
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF W/E-TARGET RECYCLING FLUX '
                  CALL EIRENE_MASJ5 (
     .             'ITARG, IPRT, IPLS, NDT, IY              ',
     .              I    , IPRT, IFL,  NDT(I,IPRT), IY)
                  WRITE (iunout,*) 'FNIX(NDT,IY) ',
     .                              FNIXB(NDT(I,IPRT),IY,IFL)
                ENDIF
10131         CONTINUE

              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIXB(NDT(I,IPRT),IY)
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEXB(NDT(I,IPRT),IY)
10132       CONTINUE

C  BALANCE CONTRIB. FROM Y-GRID RECYCLING SOURCE
          ELSEIF (NIXY(I,IPRT).EQ.2) THEN
            DO 10135 IX=NTIN(I,IPRT),NTEN(I,IPRT)-1
              IF (LCUT(IX)) GOTO 10135
              DO 10136 IFL=NSPZI(I,IPRT),NSPZE(I,IPRT)
                IF (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL).GE.0.D0)
     >           THEN
                  SFNIT(I,IFL)=SFNIT(I,IFL)-
     .             NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)-
     .             ABS(FNIY_XB(IX,NDT(I,IPRT),IFL))

cdr sheath contributions: count negative for electrons, positive for ions
cdr unfinished: need to account for charge state of ion species IFL
                  SHEAE(I)=SHEAE(I)+TEB(IX,NDT(I,IPRT))*
     .             (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)+
     .              ABS(FNIY_XB(IX,NDT(I,IPRT),IFL)))*
     .             (-DELTA_SHEATHYB(IX,NDT(I,IPRT)))
                  SHEAI(I)=SHEAI(I)+TEB(IX,NDT(I,IPRT))*
     .             (NINCT(I,IPRT)*FNIYB(IX,NDT(I,IPRT),IFL)+
     .              ABS(FNIY_XB(IX,NDT(I,IPRT),IFL)))*
     .             DELTA_SHEATHYB(IX,NDT(I,IPRT))
cdr  sheath done
                ELSE
                  WRITE (iunout,*)
     .              'WRONG ORIENTATION OF S/N-TARGET RECYCLING FLUX '
                  CALL EIRENE_MASJ5 (
     .             'ITARG, IPRT, IPLS, IX, NDT              ',
     .              I    , IPRT, IFL,  IX, NDT(I,IPRT))
                  WRITE (iunout,*) 'FNIY(IX,NDT) ',
     .                              FNIYB(IX,NDT(I,IPRT),IFL)
                ENDIF
10136         CONTINUE
              SFEIT(I)=SFEIT(I)-NINCT(I,IPRT)*FEIYB(IX,NDT(I,IPRT))
              SFEET(I)=SFEET(I)-NINCT(I,IPRT)*FEEYB(IX,NDT(I,IPRT))
10135       CONTINUE
          ENDIF
        ENDDO
C
        SFNIT(I,:)=SFNIT(I,:)*ELCHA
C
C
        WRITE (IUNOUT,*) 'FLUXES TO TARGET NO. ',I
        WRITE (IUNOUT,8888) (SFNIT(I,IFL),IFL=1,NFL),SFEIT(I),SFEET(I)
C
        SFEIT(0)=SFEIT(0)+SFEIT(I)
        SFEET(0)=SFEET(0)+SFEET(I)
        SFNIT(0,:)=SFNIT(0,:)+SFNIT(I,:)
        SHEAE(0)=SHEAE(0)+SHEAE(I)
        SHEAI(0)=SHEAI(0)+SHEAI(I)
10139 CONTINUE
#endif
C
      SSNI=0.
      SSEI=0.
      SSEE=0.
      DO ISTR=1,NSTRAI
        ISTRA = ISTR
        IF (XMCP(ISTRA).LE.1) CYCLE
        FLX=0.
        IF (ISTRA.LE.NTARGI) THEN
          FLX=SUM(ABS(SFNIT(ISTRA,1:NFLA)))
        ELSE
          FLX=1.
        ENDIF
        DO IFL=1,NFLA
          SSN(IFL)=SUM(SNI(1:NCIP,IFL,ISTRA))
        END DO
        SSI=SUM(SEI(1:NCIP,1:NFLA,ISTRA))
        SSE=SUM(SEE(1:NCIP,ISTRA))
C
        WRITE (IUNOUT,*)
     .   'RECYCLING SOURCE RATES, POTENTIAL+RAD. EN. ',ISTRA
        WRITE (IUNOUT,8888) SSN*FLX,SSI*FLX/ELCHA,SSE*FLX/ELCHA
C
C  TRENNEN VON RAD. UND POTENTIELLER ENERGY IM ELECTRONENKANAL.
C  DAZU ABER TEILCHENQUELLE SPEZIESAUFGELOEST NOETIG.
C
C
C
C     WRITE (IUNOUT,*) 'RADIATION LOSSES VIA NEUTRAL CHANNEL ',ISTRA
C     WRITE (IUNOUT,8888) 0.,0.,0.
C
        SSNI(1:NFLA)=SSNI(1:NFLA)+SSN(1:NFLA)*FLX
        SSEI=SSEI+SSI*FLX/ELCHA
        SSEE=SSEE+SSE*FLX/ELCHA
      END DO ! ISTR
C
      WRITE (IUNOUT,*) 'EQUILIBRATION '
      WRITE (IUNOUT,8888) 0.,B2QIE,-B2QIE
C
C
      WRITE (IUNOUT,*) 'BREMSSTRAHLUNG '
      WRITE (IUNOUT,8888) 0.,0.,B2BREM
C
      WRITE (IUNOUT,*) 'CHARGED IMPURITY RAD.,IONIZ. AND RECOMB. '
      WRITE (IUNOUT,8888) 0.,0.,B2RAD
C
      WRITE (IUNOUT,*) 'ELECTRIC FIELD TERMS (PRESSURE GRADIENTS)'
      WRITE (IUNOUT,8888) 0.,B2VDP,-B2VDP
C
      BALANI=SFEIN(0)+SFEIT(0)+SSEI+B2QIE+B2VDP
      BALANE=SFEEN(0)+SFEET(0)+SSEE+B2BREM+B2RAD-B2QIE-B2VDP
      BALANN(1:NFLA)=SFNIN(0,1:NFLA)+
     .               SFNIT(0,1:NFLA)+SSNI(1:NFLA)
C
      TOTI=ABS(SFEIN(0))+ABS(SFEIT(0))+
     .     ABS(SHEAI(0))+ABS(SSEI)
      TOTE=ABS(SFEEN(0))+ABS(SFEET(0))+
     .     ABS(SHEAE(0))+ABS(SSEE)
      TOTN(1:NFLA)=ABS(SFNIN(0,1:NFLA))+
     .             ABS(SFNIT(0,1:NFLA))+ABS(SSNI(1:NFLA))
      RE=BALANE/(TOTE+EPS60)*100.
      RI=BALANI/(TOTI+EPS60)*100.
      RN(1:NFLA)=BALANN(1:NFLA)/(TOTN(1:NFLA)+EPS60)*100.

c  residuals, contributions from noise in source terms.
c  sum over strata
      DO IFL=1,NFLA
        RESSNI(0,IFL) = SUM(RESSNI(1:NSTRAI,IFL))
        RESSMO(0,IFL) = SUM(RESSMO(1:NSTRAI,IFL))
      END DO
      RESSEE(0) = SUM(RESSEE(1:NSTRAI))
      RESSEI(0) = SUM(RESSEI(1:NSTRAI))
C
      CALL EIRENE_LEER (1)
      IF (LBALAN) THEN
        WRITE (iunout,*) 'B2-EIRENE GLOBAL BALANCES'
        WRITE (iunout,*) 'PARTICLE FLUXES (SFNI..) IN AMP'
        WRITE (iunout,*) 'ENERGY FLUXES (SFEI..,SFEE..,) IN WATT'
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) ' NON-RECYCLING FLUXES '
        CALL EIRENE_MASR2(' SFEIN,SFEEN    ',SFEIN(0),SFEEN(0))
        DO IFL=1,NFLA
          WRITE(iunout,'(A,I3,A,ES12.4)') 'SFNIN(IFL = ',IFL,') ',
     .                                     SFNIN(0,IFL)
        ENDDO
        ALLOCATE(SFETMP(1:NTARGI))
        SFETMP(1:NTARGI)=SFEIT(1:NTARGI)
        CALL EIRENE_MASRR1 (' TARGETS,EI',SFETMP,NTARGI,5)
        SFETMP(1:NTARGI)=SFEET(1:NTARGI)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SFETMP,NTARGI,5)
        DO ITARG=1,NTARGI
          IF (ANY(SFNIT(ITARG,1:NFLA).NE.0.0)) THEN
            DO IFL=1,NFLA
              WRITE(iunout,'(1X,A,I3,A,I3,A,ES15.8)')
     .            'TARGET ', ITARG, ', NI(IFL =',IFL,') ',
     .             SFNIT(ITARG,IFL)
            ENDDO
          ENDIF
        ENDDO
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASR2(' TOTALS, EI,EE  ',SFEIT(0),SFEET(0))
        DO IFL=1,NFLA
          WRITE(iunout,'(1X,A,I3,A,ES15.8)')
     .          'TOTALS, NI(IFL =',IFL,') ',SFNIT(0,IFL)
        ENDDO

        CALL EIRENE_LEER(2)

        WRITE (iunout,*) ' NEUTRAL PLASMA INTERACTION:'
        CALL EIRENE_MASR2(' SSEI,SSEE      ',SSEI,SSEE)
        DO IFL=1,NFLA
          WRITE(iunout,'(1X,A,I3,A,ES15.8)')
     .           'SSNI(IFL =',IFL,') ',SSNI(IFL)
        ENDDO
        CALL EIRENE_LEER(2)

        WRITE (iunout,*)
     .    'VOLUMETRIC ENERGY SINKS FOR ELECTRONS, FROM B2'
        CALL EIRENE_MASR4('B2BREM,B2RAD,-B2QIE,-B2VDP      ',
     .               B2BREM,B2RAD,-B2QIE,-B2VDP)
        WRITE (iunout,*)
     .    'TARGET SHEATH CONTRIBUTIONS, ELECTRONS AND IONS'
        SFETMP(1:NTARGI)=SHEAI(1:NTARGI)
        CALL EIRENE_MASRR1 (' TARGETS,EI',SFETMP,NTARGI,5)
        SFETMP(1:NTARGI)=SHEAE(1:NTARGI)
        CALL EIRENE_MASRR1 (' TARGETS,EE',SFETMP,NTARGI,5)
        DEALLOCATE(SFETMP)
        CALL EIRENE_MASR2('TOTALS,EI,EE    ',SHEAI(0),SHEAE(0))
        CALL EIRENE_LEER(2)

        CALL EIRENE_MASR2('BALANI,BALANE   ',BALANI,BALANE)
        DO IFL=1,NFLA
           WRITE(iunout,'(1X,A,I3,A,ES15.8)')
     .          'BALANN(IFL =',IFL,') ',BALANN(IFL)
        ENDDO
        CALL EIRENE_LEER(1)

        CALL EIRENE_MASR2('REL.ERR.(%)RI,RE',RI,RE)
        DO IFL=1,NFLA
           WRITE(iunout,'(1X,A,I3,A,ES15.8)')
     .           'RN(IFL =',IFL,') ',RN(IFL)
        ENDDO
        CALL EIRENE_LEER(1)
        MINSPEZ=99
        MAXSPEZ=-1
        MINSPEZ=MIN(MINSPEZ,MINVAL(NSPZI(1:NTARGI)))
        MAXSPEZ=MAX(MAXSPEZ,MAXVAL(NSPZE(1:NTARGI)))
        BALAN=0.
        TOT=0.
        DO IFL=MINSPEZ,MAXSPEZ
          BALAN=BALAN+SFNIN(0,IFL)+SFNIT(0,IFL)+SSNI(IFL)
          TOT=TOT+ABS(SFNIN(0,IFL))+ABS(SFNIT(0,IFL))+
     .            ABS(SSNI(IFL))
        ENDDO
        RNT=BALAN/(TOT+EPS60)*100.
        CALL EIRENE_MASJ2('SUMMED OVER     ',MINSPEZ,MAXSPEZ)
        CALL EIRENE_MASR3('BALAN,TOT,RNT           ',BALAN,TOT,RNT)

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NOISE FROM SOURCE TERMS'

        RESSNI(0,1:NFLA) = RESSNI(0,1:NFLA)/ELCHA
        RESSMO(0,1:NFLA) = RESSMO(0,1:NFLA)/ELCHA
        CALL EIRENE_MASR4('RESSEE,RESSEI,RESSNI,RESSMO     ',
     .        RESSEE(0),RESSEI(0),SUM(RESSNI(0,1:NFLA)),
     .        SUM(RESSMO(0,1:NFLA)))
        CALL EIRENE_LEER(1)

        WRITE (iunout,*) 'RESSNI-CONTRIBUTIONS BY DIFFERENT SPECIES'
        if (.not.allocated(helpw)) allocate (helpw(nfla))
        helpw(1:nfla) = RESSNI(0,1:NFLA)
        CALL EIRENE_MASRR1 ('RESSNI     ',HELPW,NFLA,5)
        CALL EIRENE_LEER(1)

        WRITE (iunout,*) 'RESSMO-CONTRIBUTIONS BY DIFFERENT SPECIES'
        helpw(1:nfla) = RESSMO(0,1:NFLA)
        CALL EIRENE_MASRR1 ('RESSMO     ',HELPW,NFLA,5)
        if (allocated(helpw)) deallocate (helpw)
      ENDIF  !LBALAN
C
c sputtering
!      WRITE(iunout,*) 'Sputtering: Total'
!      CALL EIRENE_MASYR1('ATOMS    ',SPTATI,LOGATM,0,0,NATM,NATM,0,NSTRA,
!     .   TEXTS(NSPH+1))
!      CALL EIRENE_MASYR1('MOLECULES',SPTMLI,LOGMOL,0,0,NMOL,NMOL,0,NSTRA,
!     .   TEXTS(NSPA+1))
!      CALL EIRENE_MASYR1('TEST IONS',SPTIOI,LOGION,0,0,NION,NION,0,NSTRA,
!     .   TEXTS(NSPAM+1))
!      CALL EIRENE_MASYR1('BULK IONS',SPTPLI,LOGPLS,0,0,NPLS,NPLS,0,NSTRA,
!     .   TEXTS(NSPAMI+1))
!      CALL EIRENE_MASR1('TOT. FLX',
!     .   SPTATI(0,0)+SPTMLI(0,0)+SPTIOI(0,0)+SPTPLI(0,0))

      spat = 0._dp
      spml = 0._dp
      spio = 0._dp
      sppl = 0._dp
      spat(0:natmi,0:nstrai) = sptaati(0:natmi,0:nstrai) +
     .                         sptmati(0:natmi,0:nstrai) +
     .                         sptiati(0:natmi,0:nstrai) +
     .                         sptphati(0:natmi,0:nstrai) +
     .                         sptpati(0:natmi,0:nstrai)
      spml(0:nmoli,0:nstrai) = sptamli(0:nmoli,0:nstrai) +
     .                         sptmmli(0:nmoli,0:nstrai) +
     .                         sptimli(0:nmoli,0:nstrai) +
     .                         sptphmli(0:nmoli,0:nstrai) +
     .                         sptpmli(0:nmoli,0:nstrai)
      spio(0:nioni,0:nstrai) = sptaioi(0:nioni,0:nstrai) +
     .                         sptmioi(0:nioni,0:nstrai) +
     .                         sptiioi(0:nioni,0:nstrai) +
     .                         sptphioi(0:nioni,0:nstrai) +
     .                         sptpioi(0:nioni,0:nstrai)
      sppl(0:nplsi,0:nstrai) = sptapli(0:nplsi,0:nstrai) +
     .                         sptmpli(0:nplsi,0:nstrai) +
     .                         sptipli(0:nplsi,0:nstrai) +
     .                         sptphpli(0:nplsi,0:nstrai) +
     .                         sptppli(0:nplsi,0:nstrai)

      CALL EIRENE_LEER (1)
      WRITE(iunout,*) 'Sputtering: Total'
      CALL EIRENE_MASYR1('ATOMS    ',SPAT,LOGATM,0,0,NATM,NATM,0,NSTRA,
     .   TEXTS(NSPH+1))
      CALL EIRENE_MASYR1('MOLECULES',SPML,LOGMOL,0,0,NMOL,NMOL,0,NSTRA,
     .   TEXTS(NSPA+1))
      CALL EIRENE_MASYR1('TEST IONS',SPIO,LOGION,0,0,NION,NION,0,NSTRA,
     .   TEXTS(NSPAM+1))
      CALL EIRENE_MASYR1('BULK IONS',SPPL,LOGPLS,0,0,NPLS,NPLS,0,NSTRA,
     .   TEXTS(NSPAMI+1))
      CALL EIRENE_MASR1('TOT. FLX',
     .   SPAT(0,0)+SPML(0,0)+SPIO(0,0)+SPPL(0,0))
C
      CALL EIRENE_LEER (1)
C
11000 CONTINUE
C
!pb      IF (.not.LSHORT) call EIRENE_wneutrals        ! DPC-ADD
csw 07feb2011 extra B2.5
        if(my_pe == 0) then
          if (nlemis) call eirene_extrab25_emissivity
          call eirene_wneutrals_save
          call write_f44('    ')
          call eirene_write_conbe('    ')
        endif
csw
      IF (ALLOCATED(NONREC)) DEALLOCATE(NONREC)
      IF (ALLOCATED(SFEIN)) THEN
        DEALLOCATE(SFEIN)
        DEALLOCATE(SFEEN)
        DEALLOCATE(SFNIN)
      END IF

      RETURN
C
 8888 FORMAT (3E14.6)
C
 1987 WRITE(IUNOUT,*) "ERROR IN INFCOP: CANNOT OPEN "//FORT//"11 ",
     w                "(EIRENE OUTPUT)"
      CALL EIRENE_EXIT_OWN (1)
      RETURN
      END SUBROUTINE EIRENE_IF4COP

csw mpi 07apr2010
      SUBROUTINE EIRENE_BROADCAST_INFCOP
      IMPLICIT NONE
      integer :: ier
#if( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_bcast
#endif

      if(my_pe > 0) then
#ifdef WG_RM
        if(.not.allocated(lcut)) allocate(lcut(0:NDXP))
#endif

        if(.not.allocated(pux)) then
          ALLOCATE (PUX(NCVP))
          ALLOCATE (PUY(NCVP))
          ALLOCATE (PVX(NCVP))
          ALLOCATE (PVY(NCVP))
        endif

        if(.not.allocated(invncltal)) allocate(invncltal(ncvp))

          if(.not.allocated(condense_to)) then
            ALLOCATE (condense_to(ncond))
            ALLOCATE (condense_from(ncond))
          endif

        endif

#ifdef WG_RM
      call mpi_bcast(lcut(0),ndxp+1,MPI_LOGICAL,
     .               0,MPI_COMM_WORLD,ier)
#endif

      call mpi_bcast(pux,ncvp,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(puy,ncvp,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
#ifdef WG_RM
      call mpi_bcast(puxe,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(puye,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(puxn,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(puyn,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvx,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvy,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvxe,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvye,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvxn,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(pvyn,nrad,MPI_DOUBLE_PRECISION,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(NR1STQ,1,MPI_INTEGER,
     .               0,MPI_COMM_WORLD,ier)
#endif

      call mpi_bcast(invncltal,ncvp,MPI_INTEGER,
     .               0,MPI_COMM_WORLD,ier)
      call mpi_bcast(NCOND, 1,MPI_INTEGER,
     .               0, MPI_COMM_WORLD,ier)
      call mpi_bcast(condense_to(1),NCOND,
     .               MPI_INTEGER,0, MPI_COMM_WORLD,ier)
      call mpi_bcast(condense_from(1),NCOND,
     .               MPI_INTEGER,0, MPI_COMM_WORLD,ier)

      call mpi_barrier(MPI_COMM_WORLD,ier)
      return
      end subroutine eirene_broadcast_infcop
csw
csw mpi 09jun2010
      subroutine eirene_broadcast_if3cop
      implicit none
      integer :: ifl, istrx, ier
      real(dp), allocatable :: dumvec(:)
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_bcast
#endif

      if(.not.allocated(snis0)) then
        ALLOCATE (SEES0(NSTRAi))
        ALLOCATE (SEIS0(NSTRAi))
        ALLOCATE (SNIS0(NSTRAi,0:NFL))
        ALLOCATE (SMOS0(NSTRAi,-3:3*NFL))
        ALLOCATE (RESSNI(0:NSTRAi,NFL))
        ALLOCATE (RESSMO(0:NSTRAi,NFL))
        ALLOCATE (RESSEE(0:NSTRAi))
        ALLOCATE (RESSEI(0:NSTRAi))
        ALLOCATE (FLXEIR(NSTRAi))
      endif

      do istrx=1,nstrai
        if(procforstra(istrx,my_pe) .and. my_pe /= 0) then
          SNIS0(istrx,0) = sum(snis0(istrx,1:nfl))
          SMOS0(istrx,0) = sum(smos0(istrx,1:3*nfl))
          DO IFL=1,NFL
            SMOS0(istrx,-1) = sum(smos0(istrx,      1:  nfl))
            SMOS0(istrx,-2) = sum(smos0(istrx,  nfl+1:2*nfl))
            SMOS0(istrx,-3) = sum(smos0(istrx,2*nfl+1:3*nfl))
            RESSNI(0,IFL) = ressni(0,ifl) + RESSNI(istrx,IFL)
            RESSMO(0,IFL) = ressni(0,ifl) + RESSMO(istrx,IFL)
          END DO
          RESSEE(0) = ressee(0) + RESSEE(istrx)
          RESSEI(0) = ressei(0) + RESSEI(istrx)
        endif
      enddo

      allocate(dumvec(0:nfl))
      do istrx=1,nstrai
        dumvec(0:nfl) = snis0(istrx,0:nfl)
        call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
        snis0(istrx,0:nfl) = dumvec(0:nfl)
      end do
      deallocate(dumvec)

      allocate(dumvec(-3:3*nfl))
      do istrx=1,nstrai
        dumvec(-3:3*nfl) = smos0(istrx,-3:3*nfl)
        call mpi_bcast(dumvec,3*nfl+4,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
        smos0(istrx,0:nfl) = dumvec(0:nfl)
      end do
      deallocate(dumvec)

      allocate(dumvec(0:nfl))
      do istrx=1,nstrai
        call mpi_bcast(sees0(istrx),    1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

        call mpi_bcast(seis0(istrx),    1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

        call mpi_bcast(flxeir(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

        dumvec(1:nfl) = ressni(istrx,1:nfl)
        call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
        ressni(istrx,1:nfl) = dumvec(1:nfl)

        dumvec(1:nfl) = ressmo(istrx,1:nfl)
        call mpi_bcast(dumvec, nfl+1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)
        ressmo(istrx,1:nfl) = dumvec(1:nfl)

        call mpi_bcast(ressee(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

        call mpi_bcast(ressei(istrx),   1,MPI_DOUBLE_PRECISION,
     .                 0, MPI_COMM_WORLD,ier)

      enddo
      deallocate(dumvec)
      return
      end subroutine eirene_broadcast_if3cop

      subroutine eirene_mpirecv_if3cop(istrr,irank)
      implicit none
      integer, intent(in) :: istrr, irank
      integer ier
      real(dp), allocatable :: dumvec(:)
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_recv
#endif

      allocate(dumvec(0:nfl))
      call mpi_recv(xmct(istrr),     1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      call mpi_recv(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)
      snis0(istrr,0:nfl) = dumvec(0:nfl)
      deallocate(dumvec)

      allocate(dumvec(-3:3*nfl))
      call mpi_recv(dumvec, (3*nfl+4), MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)
      smos0(istrr,-3:3*nfl) = dumvec(-3:3*nfl)
      deallocate(dumvec)

      call mpi_recv(sees0(istrr),    1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      call mpi_recv(seis0(istrr),    1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      call mpi_recv(flxeir(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      allocate(dumvec(0:nfl))
      call mpi_recv(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)
      ressni(istrr,1:nfl) = dumvec(1:nfl)

      call mpi_recv(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)
      ressmo(istrr,1:nfl) = dumvec(1:nfl)
      deallocate(dumvec)

      call mpi_recv(ressee(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      call mpi_recv(ressei(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .              ier)

      call mpi_recv(see(1,istrr), size(see(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sne_pael(1,istrr), size(sne_pael(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sne_pmel(1,istrr), size(sne_pmel(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sne_piel(1,istrr), size(sne_piel(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(see_eael(1,istrr), size(see_eael(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(see_emel(1,istrr), size(see_emel(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(see_eiel(1,istrr), size(see_eiel(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(see_epel(1,istrr), size(see_epel(:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sei(1,1,istrr), size(sei(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sei_eapl(1,1,istrr), size(sei_eapl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sei_empl(1,1,istrr), size(sei_empl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sei_eipl(1,1,istrr), size(sei_eipl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sei_eppl(1,1,istrr), size(sei_eppl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sni(1,1,istrr), size(sni(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sni_papl(1,1,istrr),
     &              size(sni_papl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sni_pmpl(1,1,istrr),
     &              size(sni_pmpl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sni_pipl(1,1,istrr),
     &              size(sni_pipl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sni_pppl(1,1,istrr),
     &              size(sni_pppl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(smo(1,1,istrr), size(smo(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(smo_mapl(1,1,istrr),
     &              size(smo_mapl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(smo_mmpl(1,1,istrr),
     &              size(smo_mmpl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(smo_mipl(1,1,istrr),
     &              size(smo_mipl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(smo_mppl(1,1,istrr),
     &              size(smo_mppl(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sna_paat(1,1,istrr),
     &              size(sna_paat(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_pmat(1,1,istrr),
     &              size(sna_pmat(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_piat(1,1,istrr),
     &              size(sna_piat(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sna_paml(1,1,istrr),
     &              size(sna_paml(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_pmml(1,1,istrr),
     &              size(sna_pmml(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_piml(1,1,istrr),
     &              size(sna_piml(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
cdjm Jan2017
      call mpi_recv(sna_paio(1,1,istrr),
     &              size(sna_paio(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_pmio(1,1,istrr),
     &              size(sna_pmio(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call mpi_recv(sna_piio(1,1,istrr),
     &              size(sna_piio(:,:,istrr)),
     &              MPI_DOUBLE_PRECISION, irank, istrr,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      return
      end subroutine eirene_mpirecv_if3cop

      subroutine eirene_mpisend_if3cop(istrr,irank)
      implicit none
      integer, intent(in) :: istrr, irank
      integer :: ier
      real(dp), allocatable :: dumvec(:)
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_send
#endif

      allocate(dumvec(0:nfl))
      call mpi_send(xmct(istrr),     1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      dumvec(0:nfl) = snis0(istrr,0:nfl)
      call mpi_send(dumvec, (nfl+1), MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)
      deallocate(dumvec)

      allocate(dumvec(-3:3*nfl))
      dumvec(-3:3*nfl) = smos0(istrr,-3:3*nfl)
      call mpi_send(dumvec, (3*nfl+4), MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)
      deallocate(dumvec)

      call mpi_send(sees0(istrr),    1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      call mpi_send(seis0(istrr),    1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      call mpi_send(flxeir(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      allocate(dumvec(0:nfl))
      dumvec(1:nfl) = ressni(istrr,1:nfl)
      call mpi_send(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      dumvec(1:nfl) = ressmo(istrr,1:nfl)
      call mpi_send(dumvec, nfl+1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)
      deallocate(dumvec)

      call mpi_send(ressee(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      call mpi_send(ressei(istrr),   1, MPI_DOUBLE_PRECISION,
     .              irank, istrr, MPI_COMM_WORLD, ier)

      call mpi_send(see(1,istrr), size(see(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sne_pael(1,istrr), size(sne_pael(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sne_pmel(1,istrr), size(sne_pmel(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sne_piel(1,istrr), size(sne_piel(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(see_eael(1,istrr), size(see_eael(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(see_emel(1,istrr), size(see_emel(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(see_eiel(1,istrr), size(see_eiel(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(see_epel(1,istrr), size(see_epel(:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sei(1,1,istrr), size(sei(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sei_eapl(1,1,istrr), size(sei_eapl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sei_empl(1,1,istrr), size(sei_empl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sei_eipl(1,1,istrr), size(sei_eipl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sei_eppl(1,1,istrr), size(sei_eppl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sni(1,1,istrr), size(sni(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sni_papl(1,1,istrr),
     &         size(sni_papl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sni_pmpl(1,1,istrr),
     &         size(sni_pmpl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sni_pipl(1,1,istrr),
     &         size(sni_pipl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sni_pppl(1,1,istrr),
     &         size(sni_pppl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(smo(1,1,istrr), size(smo(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(smo_mapl(1,1,istrr),
     &         size(smo_mapl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(smo_mmpl(1,1,istrr),
     &         size(smo_mmpl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(smo_mipl(1,1,istrr),
     &         size(smo_mipl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(smo_mppl(1,1,istrr),
     &         size(smo_mppl(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sna_paat(1,1,istrr),
     &         size(sna_paat(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_pmat(1,1,istrr),
     &         size(sna_pmat(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_piat(1,1,istrr),
     &         size(sna_piat(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sna_paml(1,1,istrr),
     &         size(sna_paml(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_pmml(1,1,istrr),
     &         size(sna_pmml(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_piml(1,1,istrr),
     &         size(sna_piml(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
cdjm Jan2017
      call mpi_send(sna_paio(1,1,istrr),
     &         size(sna_paio(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_pmio(1,1,istrr),
     &         size(sna_pmio(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
      call mpi_send(sna_piio(1,1,istrr),
     &         size(sna_piio(:,:,istrr)),
     &         MPI_DOUBLE_PRECISION, irank, istrr, MPI_COMM_WORLD, ier)
csw
      return
      end subroutine eirene_mpisend_if3cop

      FUNCTION EIRENE_POINT_ON_INTERVAL(X,Y,X1,Y1,X2,Y2)
C
C FINDS IF THE POINT (X,Y) BELONGS TO INTERVAL (X1,Y1)..(X2,Y2)
C
      IMPLICIT NONE
      LOGICAL :: EIRENE_POINT_ON_INTERVAL
      REAL(DP),INTENT(IN) :: X,Y,X1,Y1,X2,Y2
      REAL(DP) :: A,B,C,D,L, EPS4

      EIRENE_POINT_ON_INTERVAL=.FALSE.
      EPS4 = 10._DP * EPS5

      IF(X1.LT.X2.AND.
     f   (X.GT.X2+EPS4.OR.X.LT.X1-EPS4)) RETURN
      IF(X1.GT.X2.AND.
     f   (X.GT.X1+EPS4.OR.X.LT.X2-EPS4)) RETURN
      IF(Y1.LT.Y2.AND.
     f   (Y.GT.Y2+EPS4.OR.Y.LT.Y1-EPS4)) RETURN
      IF(Y1.GT.Y2.AND.
     f   (Y.GT.Y1+EPS4.OR.Y.LT.Y2-EPS4)) RETURN

      A=Y2-Y1
      B=X1-X2
      C=Y1*X2-Y2*X1
      D=A*X+B*Y+C
      D=D*D
      L=A*A+B*B
      IF(L.LT.EPS6) THEN
        WRITE(iunout,*) "WARNING FROM POINT_ON_INTERVAL"
        WRITE(iunout,*) "THE LENGTH OF THE INTERVAL IS TOO SMALL"
        WRITE(iunout,'(1x,a,1p,5(1e14.7,1x))') "L,X1,Y1,X2,Y2 ",
     .                                          L,X1,Y1,X2,Y2
        if ((abs(x1-x)/(x1+eps60) <= eps10) .and.
     .      (abs(y1-y)/(y1+eps60) <= eps10)) then
          write (iunout,*) 'vertex 1 met'
          EIRENE_POINT_ON_INTERVAL=.TRUE.
        endif
        if ((abs(x2-x)/(x2+eps60) <= eps10) .and.
     .      (abs(y2-y)/(y2+eps60) <= eps10)) then
          write (iunout,*) 'vertex 2 met'
          EIRENE_POINT_ON_INTERVAL=.TRUE.
        endif
        RETURN
      END IF
      IF(D.LT.1.D-8*L) EIRENE_POINT_ON_INTERVAL=.TRUE.

      RETURN
      END FUNCTION EIRENE_POINT_ON_INTERVAL

#ifdef WG_TODO
C DEFINE NORMAL DIRECTION FOR SURFACE-AVERAGED TALLIES (SEE FOLNEUT.F)
      SUBROUTINE CORRECTNSS
      IMPLICIT NONE
      INTEGER :: ITRI, NBAR, NBARSI
      integer :: ist, ier, it, ndtarg, itarg, ifc, is
      integer, allocatable :: inmsav(:,:)
      external :: eirene_exit_own
      TYPE(FACE), POINTER :: CURFC

      ier = 0
      inmtinss = 1
      ALLOCATE (INMSAV(size(INMTI,1),size(INMTI,2)))
      INMSAV = INMTI

      do ist = 1, nstsi
        ndtarg = 0
        do itarg = 1, ntargi
          if (-ist == ndt(itarg)) ndtarg = itarg
        end do
        if (ndtarg == 0) then
          write (iunout,*)
          write (iunout,*) ' NONDEFAULT STANDARD SURFACE ISTS=',ist
          write (iunout,*) ' IS NO RECYCLING SURFACE '
          write (iunout,*) ' NO ORIENTATION FOR NET FLUXES AVAILABLE'
          write (iunout,*)
          cycle
        end if

!  surface IST is recycling surface
!  introduce orientation in INMTINSS

        DO IT=1,TRGT(NDTARG)%NTRGDAT
          IFC = TRGT(NDTARG)%FACES(IT)
          CURFC => NCLFC(IFC)%FP
          IF (.NOT.ASSOCIATED(CURFC)) CYCLE
          DO WHILE (ASSOCIATED(CURFC))
            ITRI = CURFC%TRIANGLE
            IS = CURFC%SIDE
            NBAR = NCHBAR(IS,ITRI)
            NBARSI = NSEITE(IS,ITRI)
            INMTINSS(IS,ITRI) = NINT(TRGT(NDTARG)%FCORI(IT))
            INMSAV(IS,ITRI) = 0
            IF (NBAR > 0) THEN
              INMTINSS(NBARSI,NBAR) = -INMTINSS(IS,ITRI)
              INMSAV(NBARSI,NBAR) = 0
            END IF
            CURFC => CURFC%NEXT
          END DO
        END DO

        if (any(inmsav == -ist)) then
          write (iunout,*)
          write (iunout,*) ' NONDEFAULT STANDARD SURFACE ISTS=',ist
          write (iunout,*) ' IS NOT COMPLETELY COVERED BY TARGET ',
     .          NDTARG
          write (iunout,*)
          ier = 1
        end if
      end do

      deallocate (inmsav)

      if (ier == 1) then
        write (iunout,*)
        write (iunout,*) ' CALCULATION ABANDONNED DUE TO PROBLEMS',
     .       ' ENCOUNTERED'
        call eirene_exit_own (1)
      end if

      RETURN
      END SUBROUTINE CORRECTNSS
CVK END
#endif

      SUBROUTINE IMPORT_B2_CHEM_SPUT
      implicit none
      logical :: found
#ifndef NO_B2_CHEM_SPUT
      integer i, l, n, ilim
      logical, save:: lfirst_b2
      data lfirst_b2 /.true./
#endif

      found = .false.
#ifndef NO_B2_CHEM_SPUT
      if(fchar_chemical.gt.0.) then
        if (lfirst_b2)
     .   write(iunout,*) '%% fchar_chemical,igass_chemical =',
     ,                       fchar_chemical,igass_chemical
c
c*** The chemical sputtering data were input in B2 -
c*** let them supercede the original ones
c*** (added for compatibility with dpc)
c
cVK 180204    lchsptgt and lchsptgp definition is omitted
c

        if(igass_chemical.gt.natmi) then
          if (lfirst_b2) write(iunout,*)
     ,                 'INFUSR: igass_chemical > natmi. No use of ',
     ,                 'B2-read data for chemical sputtering!'
        else
          n=nint(fchar_chemical)
          i=1
          do while (i.le.natmi)
            l=0
            if(nchara(i).eq.n) then
              l=i
              found=.true.
            end if
            if (l.ne.0) then
              lchspnwl(l,:)=1
              if(itsput_chemical.eq.0) then
                esputc(l,:)=1 ! typical dis. energy for CD4 about 1eV
              else
                esputc(l,:)=0
              end if
              do ilim=0,nlim+nsts
                recycc(l,ilim)=bchemical_sputter_yield(ilim)
              enddo
            end if
            i = i + 1
          end do
          if(.not.found) then
            if (lfirst_b2)
     .         write(iunout,*) 'INFUSR: no match found for ',
     ,                 'specification by fchar_chemical. ',
     ,          'No use of B2-read data for chemical sputtering!'
          else
            if (lfirst_b2) write(iunout,*)
     ,       'INFUSR: B2-read data for chemical sputtering ',
     ,                                                       'are used'
            if(issput_chemical.eq.0) then
              if (lfirst_b2)
     ,         write(iunout,*) 'INFUSR: issput_chemical = 0. ',
     ,                'Chemical sputtering off carbon assumed!'
              issput_chemical=1206
            end if
            if(itsput_chemical.eq.0) then
              if (lfirst_b2)
     ,         write(iunout,*) 'INFUSR: itsput_chemical = 0. ',
     ,                '1 eV atoms assumed!'
            end if
            do i=0,nlim+nsts
              if(bchemical_sputter_yield(i).gt.0.) then
                isput(2,i)=1
              else
                isput(2,i)=0
              end if
              ilspt(i)=10*isput(2,i)+isput(1,i)
            end do
          end if
        end if
      end if
#endif
      lfirst_b2 = .false.
      return
      END SUBROUTINE IMPORT_B2_CHEM_SPUT


      subroutine eirene_if3cop_sum
      use eirmod_mpi
      use eirmod_comsou
     & , only: nstrai, nlsron
      use eirmod_eirbra
     & , only: eirene_mpisend_eirbra, eirene_mpirecv_eirbra
      use eirmod_cpes
     & , only: I_am_leader, my_pe, stratum_leader, nprs
      use eirmod_wneutrals
     & , only: eirene_wneutrals_send, eirene_wneutrals_recv,
     &         eirene_wneutrals_reduce
      implicit none
      integer :: k

      if (I_am_leader() .and. nprs > 1) then
        call eirene_wneutrals_reduce
        do k = 1, nstrai
          if ((stratum_leader(k) .ne. 0) .and. nlsron(k)) then
            if (I_am_leader(k)) then
              call eirene_mpisend_if3cop(k, 0)
              call eirene_wneutrals_send(k, 0)
              call eirene_mpisend_eirbra(k,0)
            endif
            if (my_pe==0) then
              call eirene_mpirecv_if3cop(k, stratum_leader(k))
              call eirene_wneutrals_recv(k, stratum_leader(k))
              call eirene_mpirecv_eirbra(k, stratum_leader(k))
            endif
          endif
        enddo
      endif
      return
      end subroutine eirene_if3cop_sum

      SUBROUTINE EIRENE_DEALLOC_INFCOP
      IMPLICIT NONE

      IF (ALLOCATED(CHPS)) THEN
        DEALLOCATE (CHPS)
        DEALLOCATE (SNIS)
        DEALLOCATE (CHMOS)
        DEALLOCATE (SMOS)
        DEALLOCATE (SCALN)
        DEALLOCATE (SEES0)
        DEALLOCATE (SEIS0)
        DEALLOCATE (SNIS0)
        DEALLOCATE (SMOS0)

        DEALLOCATE (RESSNI)
        DEALLOCATE (RESSMO)
        DEALLOCATE (RESSEE)
        DEALLOCATE (RESSEI)

        DEALLOCATE (FLXEIR)
      END IF

      IF (ALLOCATED(CHPM)) THEN
        DEALLOCATE (CHPM)
        DEALLOCATE (CHEEM)
        DEALLOCATE (CHEIM)
        DEALLOCATE (CHMOM)

        DEALLOCATE (PPPL_COP)
        DEALLOCATE (MPPL_COP)
        DEALLOCATE (EPPL_COP)
        DEALLOCATE (EPEL_COP)

        DEALLOCATE (PPLODA)
        DEALLOCATE (MPLODA)
        DEALLOCATE (EPLODA)
        DEALLOCATE (EPEODA)
      END IF

      IF (ALLOCATED(EAPL0)) THEN
        DEALLOCATE (EAPL0)
        DEALLOCATE (EMPL0)
        DEALLOCATE (EIPL0)
      END IF

      IF (ALLOCATED(PUX)) THEN
        DEALLOCATE (PUX)
        DEALLOCATE (PUY)
        DEALLOCATE (PVX)
        DEALLOCATE (PVY)
      END IF

      IF (ALLOCATED(CONDENSE_FROM)) THEN
        DEALLOCATE (CONDENSE_FROM)
        DEALLOCATE (CONDENSE_TO)
      END IF

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_INFCOP


      subroutine eirene_pushf (array,ndim,ifc,itri,iside)
      implicit none
      integer, intent(in) :: ndim, ifc, itri, iside
      type(face_pointer), intent(inout) :: array(ndim)
      type(face), pointer :: new

      allocate(new)
      nullify(new%next)
      new%triangle = itri
      new%side = iside

      new%next => array(ifc)%fp
      array(ifc)%fp => new

      end subroutine eirene_pushf

C> \brief Any property requiring hand-over in parallel part.
C>
C> This interfacing routine is called in the parallel part of EIRENE
C> after the broadcast of any other quantity and before MCARLO.
      SUBROUTINE EIRENE_INFCOP_PRE_MCARLO
      USE EIRMOD_WNEUTRALS
      USE EIRMOD_EXTRAB25

csw 28jan2011 extra B25, produce the files for surface properties visualisation if required
      call eirene_extrab25_srfprvsl
csw
cxpb zero out the output arrays to B2.5
      call eirene_wneutrals_clean(.false.)
cxpb
      RETURN
      END SUBROUTINE EIRENE_INFCOP_PRE_MCARLO

C> \brief Return data at the end of an EIRENE stratum
C>
C> In this subroutine data are transferred back from EIRENE to the
C> interfacing module (and from there, after possibly further
C> processing, to the external code). This subroutine is called from the
C> "strata-loop" in subr. MCARLO, after all trajectories for a
C> particular stratum ISTRA have been sampled and after all volume and
C> surface tallies have been scaled and processed to their final form.
C>
C> The call to IF3COP is controlled by the flag NMODE (input block 1).
C> At call data are expected for stratum ISTRA. There they may be
C> further prepared (e.g. normalized, or scaled to other units) for
C> transfer to the external code
      SUBROUTINE EIRENE_INFCOP_POST_STRATUM(ISTRA)
      integer, intent(in) :: istra
      RETURN
      END SUBROUTINE EIRENE_INFCOP_POST_STRATUM

C> \brief Prepare some data prior to calculation of strata but after
C> the distribution of processors has been updated
C>
      SUBROUTINE EIRENE_INFCOP_PRE_STRATA
      USE EIRMOD_WNEUTRALS
      USE EIRMOD_CPES

      if (I_am_leader()) then
         ! The stratum leaders will use extrab25 when they call if3cop.
         ! We need the stratum leader communicator, therefore we can
         ! only call wneutrals_init after pedist.
         call eirene_wneutrals_init(broadcast=.true.)
      endif

      RETURN
      END SUBROUTINE EIRENE_INFCOP_PRE_STRATA

      END MODULE EIRMOD_INFCOP
