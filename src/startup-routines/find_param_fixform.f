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
C
      SUBROUTINE EIRENE_FIND_PARAM_FIXFORM
C
C   SET DIMENSIONS (PARAMETERS) FOR ALLOCATABLE ARRAYS (TALLIES, GRIDS, ETC..) AND SET SOME FURTHER DEFAULT VALUES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMSOU, ONLY: NSTRAI
      USE EIRMOD_COMNNL, ONLY: NPRNLI
      USE EIRMOD_COMPRT, ONLY: IUNIN, IUNOUT
      USE EIRMOD_CLOGAU, ONLY: NLWRMSH, NLSPCSCL, NLSPCSCL_ON,
     .                         NLSOLEDGE, NLSPCSPT
      USE EIRMOD_CTRCEI, ONLY: TRCAMD, TRCINT, NVOLPR, NSURPR,
     .                         EIRENE_ALLOC_CTRCEI
      USE EIRMOD_CINIT, ONLY: CASENAME, DBFNAME, DBHANDLE, NDBNAMES,
     .                        INDPRO2_SAVE, INDSRC,
     .                        EIRENE_INIT_CINIT
      USE EIRMOD_JSON, ONLY : NOPTIM_IN, NRTAL_IN, NSMSTRA_IN,
     .                        INDPRO_IN, NSTRAI_IN, NTIME_IN,
     .                        DBFNAME_IN
      USE EIRMOD_COMXS, ONLY : MAXCRM, NCR_REAC, CR_KEY
      IMPLICIT NONE

      INTEGER :: INDGRD(3), INDPRO(12), IDUM(12)
      INTEGER :: IDUMMY(2), NDUMM1, NDUMM2
      INTEGER :: NFR, ISOR, NSRFSI, NRADD,
     .           NREACI,
     .           NSTSI, NLIMI, NVOLPL, NSP, ICO,
     .           NCHORI, NCHENI, NSIGSI, ID, INDIM, NSIGVI,
     .           NR1ST, NRSEP,
     .           NP1, NP2, NRKNOT, NRPLG, NPPLG,
     .           NTPER, NTTRA, NCOOR, NTET, NBMLT,
     .           NT3RD, NTSEP, NTRII, NP2ND, NPPER, NPSEP, NPPLA,
     .           NSIGCI, IREAD, NCOPIE,
     .           NRC, NRE, NLINES, LL,
     .           INMDL, IEND, ITOK, IER,
     .           IUNIN_SAVE, I1, NPRMUL,
     .           IATM, IMOL, IION, IPHOT, IPLS,
     .           ISTRA, ISPZ,
     .           NUMSEC, NINITL_READ,
     .           MOD_ADDV, NUM_COMPO,
     .           NUM_CONTRIB, ISP, ITP, IRATIO,
     .           I, J, K,
     .           I2, I3, IH, IANF, IFILE,
     .           ILINE, JCOMP, KCONTR, IREAC_ADD, IDUM1, IDUM2
      REAL(DP) :: SORIND, SORLIM, DUMM1, ROA, ZAA, ZZA, ZGA, YAA, YYA,
     .            ZIA, YP, XP, YIA, YGA,
     .            AMPTS, RDUM1, RDUM2
      LOGICAL :: NLSCL, NLTEST, NLANA, NLDRFT, NLCRR, NLERG, NLIDENT,
     .           NLONE, NLMOVIE, LINCL45, NLCASCAD, NLDFST,
     .           NLRANMAR, NLOCTREE, NEXVS, NLTRIMESH
      LOGICAL :: NLSLB, NLCRC,  NLELL, NLTRI,  NLPLG, NLFEM, NLTET,
     .           NLGEN
      LOGICAL :: NLRAD, NLPOL,  NLTOR, NLADD,  NLMLT, NLTRIM
      LOGICAL :: NLTRA, NLTRT, NLTRZ
      LOGICAL :: PLTL2D, PLTL3D, LRPSCUT
      LOGICAL :: LDEFSTOR
      LOGICAL :: UEX, NLEMIS
      LOGICAL :: LMULPL   ! multiple Ti and V..IN (per species)
                          ! due to virt. background iterations
      LOGICAL :: ldum(35)
      CHARACTER(420) :: ULINE
      CHARACTER(420) :: ZEILE, FILE45
      CHARACTER(4) :: CLAB
      CHARACTER(6) :: HANDLE
cym
      character*8, allocatable :: textal(:)
cym
      EXTERNAL :: EIRENE_COUPLE_PARAM_CONSISTENCY,
     .            EIRENE_FILEPATH_USR, EIRENE_FIND_TRIANG_DIM,
     .            EIRENE_FIND_TET_DIM, EIRENE_IF0PRM,
     .            EIRENE_INIT_PARAMS,
     .            EIRENE_READ_TOKEN, EIRENE_UPPERCASE,
     .            EIRENE_LEER, EIRENE_EXIT_OWN,
     .            FIX_INTEGER_INPUT, FIX_LOGICAL_INPUT
C
C  SET DEFAULT VALUES FOR STORAGE PARAMETERS
C
      CALL EIRENE_INIT_PARAMS
      LMULPL = .FALSE.
      NTRII=0
C
c  NEXT: BROWSE INPUT FILE AND IDENTIFY THE REAL STORAGE NEEDS.
c   e.g. NPARMI, then set the storage (for allocatable arrays): NPARM = MAX(NPARM,NPARMI)
c   in most cases then: NPARM=NPARMI

C
C  UNIT NUMBER FOR INPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15
C
      IF (IUNIN.EQ.5.OR.IUNIN.EQ.8.OR.IUNIN.EQ.10.OR.
     .    IUNIN.EQ.11.OR.IUNIN.EQ.12.OR.IUNIN.EQ.13.OR.
     .    IUNIN.EQ.14.OR.IUNIN.EQ.15) THEN
        WRITE (IUNOUT,*) 'INVALID INPUT STREAM IUNIN: ',IUNIN
        WRITE (IUNOUT,*) 'ERROR EXIT FROM FIND_PARAM.F      '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF

      IF (IUNIN.NE.5) THEN
        INQUIRE(UNIT=IUNIN,OPENED=UEX)
        IF (UEX) THEN
          REWIND IUNIN
        ELSE
          OPEN(UNIT=IUNIN,ACCESS='SEQUENTIAL',FORM='FORMATTED',
     .         ERR=7999)
        ENDIF
      ENDIF
C
      CALL EIRENE_LEER(3)

      WRITE (IUNOUT,*) 'PRINTOUT FROM EIRENE PRE-PROCESSING:'
      WRITE (IUNOUT,*) 'BROWSE INPUT FOR STORAGE NEEDS (FIND_PARAM.F)'
      WRITE (IUNOUT,*) 'IUNIN = ', IUNIN
      CALL EIRENE_LEER(1)
C
C  start browsing the header
      READ (IUNIN,'(A80)',END=6999) ZEILE
      WRITE (IUNOUT,'(A)') TRIM(ZEILE)
      CALL EIRENE_LEER(1)

c  skip further comments in header
      DO WHILE (ZEILE(1:1).EQ.'*')
        READ (IUNIN,'(A72)') ZEILE
      END DO

c start browsing block 1
      WRITE (iunout,*) '*** 1. DATA FOR OPERATING MODE'

      READ (ZEILE,6666) NPRLL,NMODE,NTCPU,NFILE,NDUMM1,NITER,
     .                  NDUMM2,NTIME
      NTIME_IN=NTIME

      READ (IUNIN,'(A72)') ZEILE
      LDEFSTOR = .FALSE.   ! INDICATES:
                           ! NO STORAGE OPTIMIZATION INPUT CARD
      IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .     INDEX(ZEILE,'t')) == 0) THEN
        LDEFSTOR = .TRUE.  ! INDICATES:
                           ! STORAGE OPTIMIZATION INPUT CARD IS READ
C   READ OPTIONAL INPUT CARD FOR STORAGE HANDLING.
C   OTHERWISE: USE DEFAULTS DEFINED ABOVE.
        READ (ZEILE,6666) NOPTIM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
     .                    NSMSTRA,NSTORAM,NGSTAL,NRTAL
        READ (IUNIN,'(A72)') ZEILE
      ENDIF
      NOPTIM_IN = NOPTIM
      NSMSTRA_IN = NSMSTRA
      NRTAL_IN = NRTAL

C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE)
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT)
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0
      NOPTM1 = MAX(NOPTM1,1)
C
      call fix_logical_input(zeile,20)
      READ (ZEILE,6665) NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
     .                  NLERG,NLIDENT,NLONE,NLMOVIE,NLDFST,
     .                  NLRANMAR,NLCASCAD,NLOCTREE,NLWRMSH,NEXVS,
     .                  NLTRIMESH,NLSPCSCL,NLSPCSCL_ON,NLSOLEDGE,
     .                  NLSPCSPT

      CALL EIRENE_INIT_CINIT

cdr scan for optional CFILE lines: path to external database files
cdr Parse input card for PATH = DBFNAME
      READ (IUNIN,'(A420)') ZEILE    !dr  probably out
      DO WHILE (ZEILE(1:1) .NE. '*') !dr  probably out
cdr  I think this outer loop is unnecessary. Identical code in input.f, without it.
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

            WRITE (IUNOUT,*) 'PATH SET FOR FILE ',TRIM(HANDLE)
            WRITE (IUNOUT,*) 'PATH = ',ZEILE(IANF:IEND)
          ELSE
            WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED'
            WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',TRIM(HANDLE),
     .                       ' IGNORED'
          END IF
          READ (IUNIN,'(A420)') ZEILE
          I1 = INDEX(ZEILE,'CFILE')
        END DO

      END DO  !dr  unnecessary outer loop: probably out
C
C
C  READ DATA FOR STANDARD MESH, 200---299


      WRITE (iunout,*) '*** 2. DATA FOR VOXEL GRID GENERATION'

C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
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
        READ (IUNIN,'(A72)') ZEILE
        call fix_logical_input(zeile,8)
        READ (ZEILE,6665) NLSLB,NLCRC,NLELL,NLTRI,NLPLG,
     .                    NLFEM,NLTET,NLGEN
        READ (IUNIN,6666) NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT,NCOOR
        N1ST = MAX(N1ST,NR1ST)
        IF (NLPLG) NPLG = MAX(NPLG,NRPLG)
        IF (NLPLG) NPPART = MAX(NPPART,NPPLG)
        NKNOT = MAX(NKNOT,NRKNOT)
        NCOORD = MAX(NCOORD,NCOOR)
        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            READ (IUNIN,*)
            IF (NLELL.OR.NLTRI) THEN
              READ (IUNIN,*)
              READ (IUNIN,*)
              IF (NLTRI) THEN
                READ (IUNIN,*)
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            READ (IUNIN,*)
            READ (IUNIN,6666) (NP1,NP2,K=1,NPPLG)
            DO 212 I=1,NR1ST
              READ (IUNIN,6664) (XP,YP,    J=1,NRPLG)
  212       CONTINUE
          ENDIF

          IF (NLFEM .OR. NLTET) THEN
            READ (IUNIN,'(A72)') ZEILE
            CLAB = ZEILE(1:4)
            CALL EIRENE_UPPERCASE(CLAB)
            IF (INDEX(CLAB,'CASE')==0) THEN
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
              WRITE (IUNOUT,*)
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
              CALL EIRENE_EXIT_OWN(1)
            END IF

            READ (ZEILE(6:),'(A66)') CASENAME
            CASENAME=ADJUSTL(CASENAME)
            LL=LEN_TRIM(CASENAME)
          END IF

          IF (NLFEM) THEN
            CALL EIRENE_FIND_TRIANG_DIM(NR1ST, NTRI, NTRII, NKNOT,
     .                                  NGITT, CASENAME)
          ENDIF

          IF (NLTET) THEN
            CALL EIRENE_FIND_TET_DIM(NR1ST, NTETRA, NTET, NCOORD,
     .                                  NGITT, CASENAME, IFOFF)
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
C  IS THERE ONE MORE LINE, OR IS NLPOL THE NEXT VARIABLE
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .         ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .           INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
             READ (IUNIN,'(A72)') ZEILE
          END DO
          IREAD=1
        ENDIF
      ENDIF
      NTRIS=NTRII
      NKNOTS=NKNOT
      N1ST = MAX(N1ST,NR1ST)
C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B
C
      IF (IREAD == 0) READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLPOL
C
      READ (IUNIN,*)
      READ (IUNIN,6666) NP2ND,NPSEP,NPPLA,NPPER
      IF (INDGRD(2).LE.5) THEN
        READ (IUNIN,6664) YIA,YGA,YAA,YYA
      ENDIF
cdr  storage for 2nd coordinate grid only, if NLPOL=T
      IF (NLPOL) N2ND = MAX(N2ND,NP2ND)
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTOR
C
      READ (IUNIN,'(A72)') ZEILE
      call fix_logical_input(zeile,3)
      READ (ZEILE,6665) NLTRZ,NLTRA,NLTRT
      READ (IUNIN,6666) NT3RD,NTSEP,NTTRA,NTPER
      IF (INDGRD(3).LE.5) THEN
        READ (IUNIN,6664) ZIA,ZGA,ZAA,ZZA,ROA
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
      IF (NLTOR.AND.NLTRA) NTTRA=NT3RD
cdr  storage for 3nd coordinate grid only, if NLTOR=T
      IF (NLTOR) N3RD = MAX(N3RD,NT3RD)

CDR STORAGE FOR TOROIDAL EFFECTS, EVEN IF
CDR NO TOROIDAL RESOLUTION IS USED
      NTOR = MAX(NTOR,NTTRA)
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLMLT

      NBMLT=0
      IF (NLMLT) THEN
        READ (IUNIN,6666) NBMLT
      ENDIF
      NBMLT=MAX0(1,NBMLT)
      NBMAX=NBMLT
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE ((ZEILE(1:1) .NE. '*') .AND.
     .          ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') +
     .            INDEX(ZEILE,'T') + INDEX(ZEILE,'t')) == 0))
         READ (IUNIN,'(A72)') ZEILE
      END DO
      IF (ZEILE(1:1) .EQ. '*') READ (IUNIN,'(A72)') ZEILE
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
      READ (ZEILE,6665) NLADD
C
      IF (NLADD) THEN
        READ (IUNIN,6666) NRADD
        NADD = MAX(NADD,NRADD)
      ENDIF

C  FIND START OF NEXT INPUT BLOCK: 3A

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C
C  READING FOR INPUT BLOCK 2 DONE
C
C
      WRITE (iunout,*) '*** 3. DATA FOR BREP SURFACES'
C  READ DATA FOR NON-DEFAULT SURFACE MODELS ON STANDARD SURFACES
C  300--349
C
      WRITE (iunout,*) '*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES'
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NSTSI
      NSTS = MAX(NSTS,NSTSI)

C  FIND START OF NEXT INPUT BLOCK: 3B

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR ADDITIONAL SURFACES 350--399
C
      WRITE (iunout,*) '*** 3B. DATA FOR ADDITIONAL SURFACES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NLIMI
      NLIM = MAX(NLIM,NLIMI)

C  FIND START OF NEXT INPUT BLOCK: 4

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
C  400--499
C
      WRITE (iunout,*) '*** 4. DATA FOR SPECIES SPECIFICATION AND'
      WRITE (iunout,*) '       ATOMIC PHYSICS MODULE'

! CHECK FOR INCLUDE LINE
      READ (IUNIN,'(A420)') ZEILE
      ULINE = ZEILE
      IREAD=1
      CALL EIRENE_UPPERCASE(ULINE)
      I1 = INDEX(ULINE,'INCLUDE')
      LINCL45 = .FALSE.
      IF (I1 > 0) THEN
c   read block 4 and block 5 from external include file, stream fort.2
        IREAD = 0
        CALL EIRENE_READ_TOKEN(ZEILE(I1+7:),' ',FILE45,ITOK,IER,.FALSE.)
        LINCL45 = .TRUE.
        IUNIN_SAVE = IUNIN
        IUNIN = 2+ifoff
        OPEN (IUNIN,FILE=FILE45,FORM='FORMATTED',ACCESS='SEQUENTIAL')
c  read comment lines on external A&M data file FILE45, stream fort.2
  401   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 401
        IREAD=1  ! now ZEILE contains the first non-comment line
                 ! from fort.2
        GOTO 402
      END IF
C
C
      IF (IREAD == 0) READ (IUNIN,*)
      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      READ (IUNIN,'(A420)') ZEILE

  402 CALL EIRENE_UPPERCASE(ZEILE)
      IEND=INDEX(ZEILE,'DEFAULT')
cdr ............................................
      IF (IEND > 0) THEN

cdr  here error exit: unfinished option, proprietary version only...
        WRITE (IUNOUT,*) 'INVALID OPTION LHYDDEF IN INPUT BLOCK 4 '
        WRITE (IUNOUT,*) 'USE LHYDDEF ONLY IN PROPRIETARY VERSIONS'
        WRITE (IUNOUT,*) 'ERROR EXIT FROM FIND_PARAM.F      '
        CALL EIRENE_EXIT_OWN(1)
      END IF
cdr ....................................
      READ (ZEILE,*) NREACI
!PB   increase number of reactions by 1 as there are still
!PB   calls to SLREAC which use reaction number NREACI+1 (SGNAL and HE_EMISS)
      NREAC = MAX(NREAC,NREACI+1)

c  Count number of reaction data sets using internal coll. rad. models, currently H, He
cdr   NRC_REAC(1): no. of H_colrad reaction rates or pop-coef data,
cdr   NRC_REAC(4): no. of He_colrad reaction rates or pop-coef data
      NCR_REAC = 0
C
cdr  count the number of reaction cards read here.
      NREAC_LINES=0
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        NREAC_LINES=NREAC_LINES+1
        IF (INDEX(ZEILE,'CR') > 0) THEN
!  coll. rad. model found
          IF (INDEX(ZEILE,'CR_') > 0) THEN
!  coll. rad. model specified in more detail
            DO I = MAXCRM, 1, -1
              IF (LEN_TRIM(CR_KEY(I)) > 0) THEN
                IF (INDEX(ZEILE,TRIM(CR_KEY(I))) > 0)
     .            NCR_REAC(I) = NCR_REAC(I) + 1
              END IF
            END DO
          ELSE
!  use default coll. rad. model for hydrogen H, formulation II
            NCR_REAC(1) = NCR_REAC(1) + 1
          END IF
        ENDIF
        READ (IUNIN,'(A72)') ZEILE
      END DO

cdr  start reading species specification block 4a,4b,4c,4d
      WRITE (iunout,*)
     .  '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      READ (IUNIN,*) NATMI
      NATM = MAX(NATM,NATMI)
cxpb New code: we are doing this for the converter that needs to know early the
cxpb  correspondence between the species in the old and new runs
cym moved as above // introduce a local variable to pass as argument to couple_param_...
cym      if(.not.allocated(TEXTA)) allocate(TEXTA(NATM))
      if(.not.allocated(TEXTAL)) allocate(TEXTAL(NATM))
cym to be evaluated - see calling order / find_param
      if(.not.allocated(NMASSA)) then
        allocate(NMASSA(NATM))
        NMASSA = 0
        COMUSR_FIRST_PASS(1) = .FALSE.
      end if
      if(.not.allocated(NCHARA)) then
        allocate(NCHARA(NATM))
        NCHARA = 0
        COMUSR_FIRST_PASS(2) = .FALSE.
      end if
cym to be evaluated

      ISPZ = 0
      DO IATM=1,NATMI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
cym
cym     TEXTA(IATM) = ZEILE(4:11)
        TEXTAL(IATM) = ZEILE(4:11)
cym
        READ (ZEILE(12:17),'(2I3)') NMASSA(IATM),NCHARA(IATM)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cdr  read 2 cards per reaction assigned to IATM, i.e.: NRC*NATMI*2 cards
cpb......................................
cdr: try to identify if there are so-called NONLINEAR BGK collisions, input flag IBGK:
cdr: to be generalized: there may be other reactions, which require multiple Ti, Vi profiles
          READ (IUNIN,'(A72)') ZEILE
          call fix_integer_input(zeile,12)
          READ (ZEILE,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULPL = LMULPL .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULPL = LMULPL .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULPL = LMULPL .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      WRITE (iunout,*)
     .  '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,*) NMOLI
      NMOL = MAX(NMOL,NMOLI)
      if(.not.allocated(NMASSM)) then
        allocate(NMASSM(NMOL))
        NMASSM = 0
        COMUSR_FIRST_PASS(3) = .FALSE.
      end if

      DO IMOL=1,NMOLI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        READ (ZEILE(12:14),'(I3)') NMASSM(IMOL)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cpb......................................
cdr: try to identify if there are so-called BGK collisions, input flag IBGK:
cdr: to be generalized: there may be other reactions, which require multiple (IPLS) profiles
          READ (IUNIN,'(A72)') ZEILE
          call fix_integer_input(zeile,12)
          READ (ZEILE,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULPL = LMULPL .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULPL = LMULPL .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULPL = LMULPL .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      WRITE (iunout,*)
     . '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,*) NIONI
      NION = MAX(NION,NIONI)
      if(.not.allocated(NMASSI)) then
        allocate(NMASSI(NION))
        NMASSI = 0
        allocate(NCHARI(NION))
        NCHARI = 0
        allocate(NCHRGI(NION))
        NCHRGI = 0
cym this variable has to be moved from extraB25 to EIRENE as an extra optional parameter (comusr)
!pb variables LKIND? have been removed from extraB25 and put into EIRMOD_COMUSR
        allocate(LKINDI(NION))
        LKINDI = 0
        COMUSR_FIRST_PASS(4) = .FALSE.
      end if

      DO IION=1,NIONI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        READ (ZEILE(12:14),'(I3)') NMASSI(IION)
        READ (ZEILE(15:17),'(I3)') NCHARI(IION)
        READ (ZEILE(21:23),'(I3)') NCHRGI(IION)
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        READ (ZEILE(45:47),'(I3)') LKINDI(IION)
        DO K=1,NRC
cpb......................................
cdr: try to identify if there are so-called BGK collisions, input flag IBGK:
cdr: to be generalized: there may be other reactions, which require multiple (IPLS) profiles
          READ (IUNIN,'(A72)') ZEILE
          call fix_integer_input(zeile,12)
          READ (ZEILE,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULPL = LMULPL .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULPL = LMULPL .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULPL = LMULPL .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO

C  FIND START OF NEXT INPUT BLOCK: 4D
      WRITE (iunout,*)
     .  '*** 4D. PHOTONS SPECIES CARDS, NPHOTI SPECIES '
      READ (IUNIN,'(A72)') ZEILE
      NPHOTI=0
      IF (ZEILE(1:3) == '***') GOTO 500
      READ (IUNIN,*) NPHOTI
      NPHOT = MAX(NPHOT,NPHOTI)

      DO IPHOT=1,NPHOTI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
cdr
        READ (ZEILE(30:35),'(2I3)') NUMSEC,NRC
        DO K=1,NRC
cpb......................................
cdr: try to identify if there are so-called BGK collisions, input flag IBGK:
cdr: to be generalized: there may be other reactions, which require multiple Ti profiles
          READ (IUNIN,'(A72)') ZEILE
          call fix_integer_input(zeile,12)
          READ (ZEILE,'(12I6)') IDUM(1:12)
          IF (NUMSEC < 3) THEN
            LMULPL = LMULPL .OR. (IDUM(7) /= 0)
          ELSEIF (NUMSEC == 3) THEN
            LMULPL = LMULPL .OR. (IDUM(8) /= 0)
          ELSEIF (NUMSEC == 4) THEN
            LMULPL = LMULPL .OR. (IDUM(9) /= 0)
          END IF
cpb.......................................
          READ (IUNIN,*)
        END DO
      END DO
C
C  READ DATA FOR PLASMA BACKGROUND, 500--599
C
  500 WRITE (iunout,*) '*** 5. DATA FOR PLASMA BACKGROUND'
C
C  READ BULK IONS SPECIES CARDS
C
      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ(ZEILE,6666) NPLSI
      NPLS = MAX(NPLS,NPLSI)

c  count special plasma background models:
      ICO = 0
      DO IPLS=1,NPLSI
        READ (IUNIN,'(A72)') ZEILE
        ISPZ = ISPZ + 1
        ULINE = ZEILE
        CALL EIRENE_UPPERCASE(ULINE)
        INMDL=INDEX(ULINE,'FORT')+
     .        INDEX(ULINE,'FTN')+
     .        INDEX(ULINE,'SAHA')+
     .        INDEX(ULINE,'CORONA')+
     .        INDEX(ULINE,'BOLTZMANN')+
     .        INDEX(ULINE,'PLANCK')+
     .        INDEX(ULINE,'COLRAD')+
     .        INDEX(ULINE,'CONSTANT')
        IF (INMDL > 0) ICO = ICO + 1

        READ (ZEILE(33:35),'(I3)') NRC
        DO K=1,NRC
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO

        IF (INMDL > 0) THEN
          NRE=0
          IF (VERIFY(ULINE(INMDL+11:),' ') > 0)
     .      READ(ZEILE(INMDL+11:),*) NRE
          NRE=MAX(NRE,1)
          NLINES=1
          IF (INDEX(ULINE(INMDL:),'COLRAD') > 0) NLINES=NRE
          DO K=1,NLINES
             READ (IUNIN,*)
          END DO
        END IF
      END DO

cdr  this next line is probably not needed.
cdr  ico > 0 indicates: at least one bulk species has a special
cdr  background data model,  fort.., saha, corona, ...etc...
      IF (ICO > 0) NREAC=NREAC+1

      READ (IUNIN,'(A72)') ZEILE
      WRITE (IUNOUT,*) '*** 5B. PLASMA BACKGROUND DATA'
      DO WHILE (ZEILE(1:1) == '*')
         READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) (INDPRO(J),J=1,12)

!     SAVE VALUES READ FROM INPUT FILE
      INDPRO_IN = INDPRO

      call eirene_Ti_input(indpro2_save,lincl45)
      indpro(2)=indpro2_save

      NPLSTI = NPLS
      IF (ABS(INDPRO(2))>9) NPLSTI=1

      IF ((NPLS > 1) .AND. (NPLSTI == 1) .AND. LMULPL) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'TIIN IS PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(2) > 10'
        WRITE (IUNOUT,*) 'DIMENSION OF TIIN ARRAY OVERWRITTEN'
        WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE PRESENT'
        NPLSTI = NPLS
        WRITE (IUNOUT,*) ' NPLSTI = ',NPLSTI
      END IF

cdr these next 2 lines for V_IN(ipls)
      NPLSV = NPLS
      IF (MOD(ABS(INDPRO(4)),100) > 9) NPLSV = 1

      IF ((NPLS > 1) .AND. (NPLSV == 1) .AND. LMULPL) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'V_IN PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(4) > 10'
        WRITE (IUNOUT,*) 'DIMENSION OF V_IN ARRAYS OVERWRITTEN'
        WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE PRESENT'
        NPLSV = NPLS
        WRITE (IUNOUT,*) ' NPLSV = ',NPLSV
      END IF

      IF (LINCL45) THEN
        CLOSE (IUNIN)
        IUNIN = IUNIN_SAVE
        LINCL45 =.FALSE.
      END IF

C  FIND START OF NEXT INPUT BLOCK: 6

      DO
        IF (ZEILE(1:5) == '*** 6') EXIT
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR REFLECTION MODEL 600--699
C
      WRITE (iunout,*) '*** 6. GENERAL DATA FOR REFLECTION MODEL'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6665) NLTRIM
      NFR=0
      IF (NLTRIM) THEN
        READ (IUNIN,'(A72)') ZEILE
        IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').NE.0) THEN
C  PATH SPECIFICATION FOR DATABASE FOUND
          READ (IUNIN,'(A72)') ZEILE
          DO WHILE ((INDEX(ZEILE,'ON')+INDEX(ZEILE,'on')) > 0)
            NFR=NFR+1
            READ (IUNIN,'(A72)') ZEILE
          END DO
        ENDIF
      ENDIF

C  NHD6 FOR STORAGE ON ALLOCATABLE ARRAYS: number of TRIM target-projectile combinations
      IF (NFR > 0) THEN
        NHD6 = NFR
      ELSE
        NHD6 = 12
      END IF

C  FIND START OF NEXT INPUT BLOCK: 7

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR PRIMARY SOURCE  700--799
C
      WRITE (iunout,*) '*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NSTRAI
      NSTRAI_IN = NSTRAI

CDR  TRY TO SET NSTEP, THE NUMBER OF STEP FUNCTIONS FOR SOURCE SAMPLING
cdr  set nstep = highest stratum number, which receives primary source data from external code.
cdr  this must be highly case specfic. To be reconsidered !!
cpb  Step functions are only used in conjunction with source sampling.
cpb  If we find the highest stratum J that uses a step function we can be sure that there are
cpb  less than J step functions involved.
cpb  For this reason start the loop at NSTRAI and count downwards.
cpb  We imply implicitly that the numbering of step functions is in ascending order starting with 1
      NSTEP = 1
      ALLOCATE (INDSRC(NSTRAI))
      READ (IUNIN,6666) (INDSRC(J),J=1,NSTRAI)
      IF (ANY(INDSRC == 6)) THEN
        DO J=NSTRAI,1,-1
          IF (INDSRC(J) == 6) THEN
            NSTEP = J
            EXIT
          END IF
        END DO
      END IF

      READ (IUNIN,6664) RDUM1, AMPTS
      IF(AMPTS.EQ.0.0_DP) AMPTS=1.0_DP
      DO ISTRA=1,NSTRAI
        IF (INDSRC(ISTRA) == 6) CYCLE
C * ZEILE...: STRATUM NAME
        READ (IUNIN,'(A72)') ZEILE
        WRITE (IUNOUT,'(A1,A)') ' ',TRIM(ZEILE)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
cdr  number of substrata
        READ (IUNIN,6666) NSRFSI
        NSRFS = MAX(NSRFS,NSRFSI)

        DO I=1,NSRFSI
          READ (IUNIN,6666) ID, INDIM
          READ (IUNIN,6664) DUMM1, SORLIM, SORIND
          IF (INDIM == 4) NSTEP = MAX(NSTEP,NINT(SORIND))
          ISOR = NINT(SORLIM)
          DO WHILE (ISOR > 0)
cdr  here NSTEP is set to the largest step function number specified on SORIND
            ID = MOD(ISOR,10)
            IF ((ID == 4).OR.(ID==5)) NSTEP = MAX(NSTEP,NINT(SORIND))
            ISOR = ISOR / 10
          END DO
          READ (IUNIN,*)
          READ (IUNIN,*)
        END DO
        READ (IUNIN,*)
        READ (IUNIN,*)
C
      END DO
      IF (NTIME.GE.1) NSTRAI = NSTRAI + 1
      NSTRA = MAX(NSTRA,NSTRAI)

      READ (IUNIN,'(A72)') ZEILE
C
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C
      WRITE (iunout,*) '*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C
C  READ DATA FOR STATISTICS AND NON-ANALOG MODEL, 900--999
C
      WRITE (iunout,*) '*** 9. DATA FOR STATISTIC AND NON-ANALOG MODEL'

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .NE. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR STANDARD DEVIATION
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'
      READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI
      NSD = MAX(NSD,NSIGVI)
      NSDW = MAX(NSDW,NSIGSI)
      NCV = MAX(NCV,NSIGCI)

C  FIND START OF NEXT INPUT BLOCK: 10

      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C   READ DATA FOR ADDITIONAL AND SURFACE-AVERAGED TALLIES
      WRITE (iunout,*)
     . '*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION'
      WRITE (iunout,*)
     . '        ESTIMATORS AND ALGEBRAIC EXPRESSIONS'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC

      NADV = MAX(NADV,NADVI)
      NCLV = MAX(NCLV,NCLVI)
      NALV = MAX(NALV,NALVI)
      NADS = MAX(NADS,NADSI)
      NALS = MAX(NALS,NALSI)
C
      WRITE (iunout,*) '*** 10A. DATA FOR ADDITIONAL TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10B. DATA FOR COLLISION ESTIMATORS'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:2) .NE. '**')
        READ (IUNIN,'(A72)') ZEILE
      END DO

      WRITE (iunout,*) '*** 10F. DATA FOR SPECTRA'
      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3) .NE. '***') THEN
        READ (IUNIN,'(A72)') ZEILE
        DO WHILE (ZEILE(1:3) .NE. '***')
          READ (IUNIN,'(A72)') ZEILE
        END DO
      END IF
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      WRITE (iunout,*)
     .  '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT'
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      call fix_logical_input(zeile,35)
      READ (ZEILE,6665) ldum(1:35) ! in case we need
                                   ! the output switches early
C
C   READ TRCSRC (60 LOGICALS PER LINE)
      do j=0, NSTRAI, 60
        READ (IUNIN,*)
      end do

      READ (IUNIN,6666) IDUMMY(1)
      NVLPR = IDUMMY(1)
C  ERGODIC OPTION NEEDS PRINTOUT OF VOLUME, AND ONE, TWO OR THREE FURTHER TALLIES AT LEAST
      IF (NLERG) NVLPR=MAX(4,NVLPR)
      DO J=1,IDUMMY(1)
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,6666) IDUMMY(2)
      NSRPR = IDUMMY(2)
C  ERGODIC OPTION NEEDS PRINTOUT AT LEAST FROM TIME HORIZON
      IF (NLERG) NSRPR=MAX(1,NSRPR)
      DO J=1,IDUMMY(2)
        READ (IUNIN,*)
      END DO

      CALL EIRENE_ALLOC_CTRCEI(1)
      NVOLPR = IDUMMY(1)
      NSURPR = IDUMMY(2)
      TRCAMD = ldum(10)
      TRCINT = ldum(11)

C  SKIP READING ALSO POSSIBLE LINES FOR DELIBERATE DE-ACTIVATION OR RE-ACTIVATION OF TALLIES
c     to be written: allow for comment lines here

C  SEARCH START OF PLOTTING INPUT: NEXT LINE WITH F OR T
      READ (IUNIN,'(A72)') ZEILE
      CALL EIRENE_UPPERCASE(ZEILE)
      DO WHILE ((SCAN(ZEILE,'FT') == 0) .OR. (ZEILE(1:1) == '*'))
        READ (IUNIN,'(A72)') ZEILE
        CALL EIRENE_UPPERCASE(ZEILE)
      END DO
C
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C  FIRST CARD FOR (LOGICAL) PLOTTING FLAGS NOW ON 'zeile'
      READ (ZEILE(16:16),'(L1)') LRPSCUT  ! needed below for further
                                          ! storage considerations

c  there are 13 geometry plotting input cards following, before plotting of volume tallies
      DO J=1,13
        READ (IUNIN,*)
      END DO
C
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  DATA FOR PLOTS OF VOLUME-AVERAGED TALLIES
      READ (ZEILE,6666) NVOLPL
C
C
      NPLT = 1
      IF (NVOLPL > 0) THEN
C   READ PLTSRC (60 LOGICALS PER LINE)
!   backward compatible with Eirene_96
        READ (IUNIN,'(A72)') ZEILE
        IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .       INDEX(ZEILE,'t')) .ne. 0) THEN
          DO J=61, NSTRAI, 60
            READ (IUNIN,*)
          END DO
          IF (LRPSCUT) READ (IUNIN,*)
          IREAD = 0
        ELSE
          IREAD = 1
        ENDIF

        DO J=1,NVOLPL
          IF (IREAD.EQ.0 .OR. J.GT.1) READ (IUNIN,'(A72)') ZEILE
          DO WHILE (ZEILE(1:1) .EQ. '*')
            READ (IUNIN,'(A72)') ZEILE
          END DO
! For compatibility with Eirene_96
          READ (ZEILE,6666) IDUMMY(1:2)
          IF (IDUMMY(2).NE.0) THEN
            NSP = IDUMMY(2)
          ELSE
            NSP = IDUMMY(1)
          ENDIF
          NPLT = MAX(NPLT, NSP)
          READ (IUNIN,'(A72)') ZEILE
          call fix_logical_input(zeile,2)
          READ (ZEILE,6665) PLTL2D,PLTL3D
          READ (IUNIN,*)
          IF (PLTL2D) THEN
            READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
          ENDIF
          IF (PLTL3D) THEN
            READ (IUNIN,*)
            IF (IDUMMY(2).EQ.0) READ (IUNIN,*)
            DO I=1,NSP
              READ (IUNIN,*)
            END DO
            READ (IUNIN,*)
          ENDIF
        END DO
      END IF
C
C  SKIP INPUT LINES, UNTIL INPUT BLOCK 12 STARTS
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR DIAGNOSTIC MODULE  1200--1299
C
      WRITE (iunout,*) '*** 12. DATA FOR DIAGNOSTIC MODULE'
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:1) .EQ. '*')
        READ (IUNIN,'(A72)') ZEILE
      END DO
      IREAD=1

c  optional input cards: 'DEFINE_LINES'

c  read up to NUM_LINES transitions (volumetric line emissions).
c  For each transition specify NUM_COMPO different parent (donor) state components.
c  Each COMPONENT may consist of NUM_CONTRIB isotopic contributions

c  Identify the block of LINES and COMPONENTS available in this run
C  by an extra input card containing 'DEFINE_LINES'
      ULINE=ZEILE
      CALL EIRENE_UPPERCASE(ULINE)
      NADV_ADD = 0
      NLEMIS = .FALSE.
      NUM_LINES = 0

      IF (INDEX(ULINE,'DEFINE_LINES') > 0) THEN
        IREAD=0
CDR AT LEAST ONE (OR MORE) VOLUMETRIC LINE EMISSIVITY TALLY DEFINED IN INPUT BLOCK 12
cdr as additional output tally ADDV(...).
        IREAC_ADD = 0
        NLEMIS = .TRUE.
c  read number of lines, and the flag MOD_ADDV for storage mode on ADDV tallies
        READ (IUNIN,6666) NUM_LINES, MOD_ADDV
        DO ILINE=1, NUM_LINES
          READ (IUNIN,'(A80)') ZEILE
          DO WHILE (ZEILE(1:1) == '*')
            READ (IUNIN,'(A80)') ZEILE
          END DO
c
          READ (IUNIN,6666) NUM_COMPO  ! components of line ILINE
          READ (IUNIN,*)
          IF (MOD_ADDV == 0) THEN
cdr  minimal storage, but each time when a new line comes,
cdr  the emissivity profiles on ADDV must be re-calculated
            NADV_ADD = MAX(NADV_ADD, (NUM_COMPO + 1))
          ELSE
cdr  all possible emissivity profiles are kept on ADDV tallies.
            NADV_ADD =     NADV_ADD +(NUM_COMPO + 1)
          END IF
c
          DO JCOMP=1, NUM_COMPO
            READ (IUNIN,*)
            READ (IUNIN,*) NUM_CONTRIB     ! contributions to component
                                           ! JCOMP for line ILINE
c           write (iunout,*) 'num_contrib', num_contrib
            IREAC_ADD = IREAC_ADD + NUM_CONTRIB
cdr  Specify all required contributions explicitly.
cdr  In the old default this was automatically detected
cdr     from mass and charge states/numbers of hydrogenic particles.
cdr     And only one set of emission data for all contributions was used,
cdr     plus one or two population ratios.
cdr     Now we provide storage for one additional AM data set for each contribution,
cdr     plus one or two population ratios.
            DO KCONTR = 1, NUM_CONTRIB
              READ (IUNIN,6666) ISP, ITP, IRATIO
cdr do we require a QSS population ratio for this contribution?
              IF (IRATIO > 0) THEN
                READ (IUNIN,*)
cdr do we require a second QSS population ratio for this contribution?
                IF (IRATIO == 2) READ (IUNIN,*)
              END IF  !  IRATIO
            END DO    !  NUM_CONTRIB   (POSSIBLE D, H, T CONTRIBUTE
                      !                 TO GROUND STATE EMISSIVITY)
          END DO      !  NUM_COMPO     (E.G. GROUND STATE)
        END DO        !  NUM_LINES     (E.G. BA-ALPHA)

c  STORAGE FOR ADDITIONAL TALLIES NADV_ADD, AND REACTIONS IREAC_ADD (LINE EMISSIVITIES)
        NADV = NADV + NADV_ADD
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
          READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1
      ELSEIF (INDEX(ULINE,'DEFAULT_LINES') > 0) THEN
        IREAD=0
c  set old default line emissivity tallies ADDV.
        NLEMIS=.TRUE.
        READ (IUNIN,'(A80)') ZEILE
        DO WHILE (ZEILE(1:1) == '*')
            READ (IUNIN,'(A80)') ZEILE
        END DO
        IREAD=1
      END IF

      WRITE (iunout,*) '*** 12A. DATA FOR LINES OF SIGHT'
      IF (IREAD.EQ.0) READ (IUNIN,'(A80)') ZEILE
      READ (ZEILE,6666) NCHORI,NCHENI
      IREAD=0

      NCHOR = MAX(NCHOR,NCHORI)
      NCHEN = MAX(NCHEN,NCHENI)

cdr this next condition for old default: better also check for nchtal=2 ??
      NLEMIS = NLEMIS .OR. (NCHOR > 0)
      IF (NLEMIS.AND.(NUM_LINES == 0)) THEN
! USE OLD HYDROGENIC DEFAULT LINES FOR EMISSIVITY
        NUM_LINES = 6
        NUM_COMPO = 6
        MOD_ADDV = 0
        NADV=NADV +7
! USE MAXIMUM POSSIBLE NUMBER OF CONTRIBUTIONS.
cdr MAY BE FAR TOO LARGE, AS NCHAR AND NCHRG ARE NOT YET AVAILABLE
cdr tbd: check consequences for storage footprint
cdr:  corrected in "emiss branch"
        NUM_CONTRIB = NATMI + NMOLI + 2*NMOLI + 2*NMOLI + 2*NMOLI +
     .                NPLSI

c  each component of a line needs one reaction card.
c  for the hydrogen default model, there are 6 lines, with 6 components
c  each, and in total 3 density ratios (the same for all lines and components)
c  additionally 3 reaction cards are needed for ratios
        NREAC = NREAC + NUM_LINES*NUM_COMPO + 3

      END IF

C  SKIP READING REST OF THIS BLOCK
      READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO
C
C  READ DATA FOR TIME-DEPENDENT AND NONLINEAR MODE  1300--1399
C
      IREAD=0
      WRITE (iunout,*)
     .  '*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION'

C
      READ (IUNIN,'(A72)') ZEILE
      call fix_integer_input(ZEILE,3)
      READ (ZEILE,'(12I6)') NPRNLI, NINITL_READ, NPRMUL
      IF (NPRMUL > 1) THEN
        NPRNLI = NPRNLI * NPRMUL
cdr apply amplification factor for "large runs" also to census storage
      ELSEIF (AMPTS.NE.1.0_DP) THEN
        NPRNLI = INT(REAL(NPRNLI) * AMPTS)
        WRITE (iunout,*)
     .   ' CENSUS STORAGE NPRNLI ENHANCED BY FACTOR AMPTS ',AMPTS
      ENDIF
      NPRNL = MAX(NPRNL,NPRNLI)

      IF (NPRNL.LE.0.OR.NTIME.EQ.0) THEN
C  TURN OFF TIME DEP MODE IF EITHER NTIME=0 OR NPRNL=0
        IF (NPRNL.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NTIME=0'
          NPRNLI=0
          NPRNL=0
        ENDIF
        IF (NTIME.GT.0) THEN
          WRITE (IUNOUT,*) 'TIME DEP. MODE TURNED OFF, BECAUSE NPRNL=0'
          NTIME=0
        ENDIF
      ENDIF

!PB if NTIME >= 1 NSTRAI has been increased already
!PB therefore if NPRNLI <=0 reduce NSTRAI and NSTRA
      if (NTIME.GE.1.AND.NPRNLI <= 0) THEN
        NSTRAI=NSTRAI-1
        NSTRA=NSTRA-1
      ENDIF
      if ((NTIME.GE.1.AND.NPRNLI > 0).OR.NLERG) THEN
        NSTSI=NSTSI+1
        if (NLERG .AND.(NTIME .LT. 1)) NSTRAI=NSTRAI+1
      ENDIF
      NSTS = MAX(NSTS,NSTSI)
      NSTRA = MAX(NSTRA,NSTRAI)
      NLIMPS = NLIM + NSTS

      READ (IUNIN,'(A72)') ZEILE
      IF (ZEILE(1:3).EQ.'***') THEN
         IREAD = 1
      ELSE
         READ (ZEILE,'(12I6)') IDUM1, IDUM2
         READ (IUNIN,'(6E12.4)') RDUM1, RDUM2
         IREAD = 0
      END IF

C   SNAPSHOT TALLIES AND CENSUS ARRAY
      NSNVI=0
      IF (NPRNLI > 0) THEN
        IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).NE.'*') THEN
          READ (IUNIN,*)
C   READ DATA FOR SNAPSHOT TALLIES
          READ (IUNIN,*)
          WRITE (iunout,*) '*** 13A. DATA FOR SNAPSHOT TALLIES'
          READ (IUNIN,6666) NSNVI
          IREAD=0  !  fix by XB
        END IF
      END IF
      NSNV = MAX(NSNV,NSNVI)

      IF (NLERG.AND.NPRNLI.LE.0) THEN
C  NO TIME HORIZON DEFINED, DESPITE NLERG=.TRUE.
C  THEREFORE: SET A DEFAULT TIME HORIZON HERE
        IF (NTIME.EQ.0) NTIME=1
        NPRNLI=100
      ENDIF
      NSTS = MAX(NSTS,NSTSI)
      NSTRA = MAX(NSTRA,NSTRAI)
      NPRNL = MAX(NPRNL,NPRNLI)
      NLIMPS = NLIM+NSTS
      CALL EIRENE_ALLOC_CTRCEI(2)

cdr       NPRNL is only valid for writing census arrays onto fort.15
cdr  tbd: when reading fort.15 (census), the size is determined by the
cdr       size of that file, (IPRNL) not by NPRNL

C  SKIP READING REST OF THIS BLOCK
      IF (IREAD.EQ.0) READ (IUNIN,'(A72)') ZEILE
      DO WHILE (ZEILE(1:3) .NE. '***')
        READ (IUNIN,'(A72)') ZEILE
      END DO

C
C  INPUT BLOCK 14 BEGIN
C
C  READ DATA IN INTERFACING SUBROUTINE INFCOP  1400 -- 1499
C
      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'
      IF (NMODE.EQ.0) THEN
!pb     NCOPII is completely redundant
!pb     READ (IUNIN,6666) NAINI,NCOPII,NCOPIE
        READ (IUNIN,6666) NAINI,NDUMM1,NCOPIE
        NCPVI=NCOPIE
      ELSE
        NAINI=0
        NCPVI=0
        CALL EIRENE_IF0PRM(IUNIN)
      ENDIF
      DEALLOCATE (INDSRC)

cdr  some parameters may have gotten changed in IF0PRM, case-specific
      NAIN = MAX(NAIN,NAINI)
      NCPV = MAX(NCPV,NCPVI)

cdr  due to these changes there, also some derived storage parameters may have changed....
      NRAD=MAX(N1ST*N2ND*N3RD,NTRI*N3RD,NTETRA)+NADD+1 ! as in parmmod

      IF (IUNIN.NE.5) REWIND IUNIN
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'AUTOMATED STORAGE SETTING (FIND_PARAM.F)'
      CALL EIRENE_LEER(1)
C
cdr  grid size
      WRITE (iunout,'(a14,i8)') 'N1ST        = ',N1ST
      WRITE (iunout,'(a14,i8)') 'N2ND        = ',N2ND
      WRITE (iunout,'(a14,i8)') 'N3RD        = ',N3RD
      WRITE (iunout,'(a14,i8)') 'NADD        = ',NADD
      WRITE (iunout,'(a14,i8)') 'NTOR        = ',NTOR
      WRITE (iunout,'(a14,i8)') 'NRTAL       = ',NRTAL
      WRITE (iunout,'(a14,i8)') 'NLIM        = ',NLIM
      WRITE (iunout,'(a14,i8)') 'NSTS        = ',NSTS
      WRITE (iunout,'(a14,i8)') 'NPLG        = ',NPLG
      WRITE (iunout,'(a14,i8)') 'NPPART      = ',NPPART
      WRITE (iunout,'(a14,i8)') 'NKNOT       = ',NKNOT
      WRITE (iunout,'(a14,i8)') 'NTRI        = ',NTRI
      WRITE (iunout,'(a14,i8)') 'NTETRA      = ',NTETRA
      WRITE (iunout,'(a14,i8)') 'NCOORD      = ',NCOORD

      WRITE (iunout,*) ' '
      WRITE (iunout,'(a14,i8)') 'NRAD        = ',NRAD
cdr  primary source
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NSTRA       = ',NSTRA
      WRITE (iunout,'(a14,i8)') 'NSRFS       = ',NSRFS
      WRITE (iunout,'(a14,i8)') 'NSTEP       = ',NSTEP
      WRITE (iunout,'(a14,i8)') 'NPTRGT      = ',NPTRGT
cdr species
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(a14,i8)') 'NATM        = ',NATM
      WRITE (iunout,'(a14,i8)') 'NMOL        = ',NMOL
      WRITE (iunout,'(a14,i8)') 'NION        = ',NION
      WRITE (iunout,'(a14,i8)') 'NPHOT       = ',NPHOT
      WRITE (iunout,'(a14,i8)') 'NPLS        = ',NPLS

      CALL EIRENE_LEER(1)
c  additional volumetric tallies
      WRITE (iunout,'(a14,i8)') 'NADV        = ',NADV
      WRITE (iunout,'(a14,i8)') 'NCLV        = ',NCLV
      WRITE (iunout,'(a14,i8)') 'NSNV        = ',NSNV
      WRITE (iunout,'(a14,i8)') 'NALV        = ',NALV
C     WRITE (iunout,'(a14,i8)') 'NBGV        = ',NBGV ! determined later
C     WRITE (iunout,'(a14,i8)') 'NCPV        = ',NCPV ! determined later
c  additional surface-averaged tallies
      WRITE (iunout,'(a14,i8)') 'NADS        = ',NADS
      WRITE (iunout,'(a14,i8)') 'NALS        = ',NALS
c  additional input tallies
      WRITE (iunout,'(a14,i8)') 'NAIN        = ',NAIN

      CALL EIRENE_LEER(1)
c  statistical variances, covariances
      WRITE (iunout,'(a14,i8)') 'NSD         = ',NSD
      WRITE (iunout,'(a14,i8)') 'NSDW        = ',NSDW
      WRITE (iunout,'(a14,i8)') 'NCV         = ',NCV

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'MAX. NO. OF ATOMIC/MOLECULAR "REACTIONS"'
      WRITE (IUNOUT,'(a14,i8)') 'NREAC       = ',NREAC
      WRITE (IUNOUT,*) 'NREC,NREI,NRCX,NREL,NRPI: DETERMINED LATER'
      WRITE (IUNOUT,*) 'NRPH,NRBG,NROT          : DETERMINED LATER'
C     WRITE (iunout,'(a14,i8)') 'NREC        = ',NREC
C     WRITE (iunout,'(a14,i8)') 'NREI        = ',NREI
C     WRITE (iunout,'(a14,i8)') 'NRCX        = ',NRCX
C     WRITE (iunout,'(a14,i8)') 'NREL        = ',NREL
C     WRITE (iunout,'(a14,i8)') 'NRPI        = ',NRPI
C     WRITE (iunout,'(a14,i8)') 'NRBG        = ',NRBG
C     WRITE (iunout,'(a14,i8)') 'NRPH        = ',NRPH

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF STORAGE OPTIMIZATION OPTIONS'
      IF (.NOT.LDEFSTOR) THEN
C  ADJUST SOME DEFAULT STORAGE OPTIMIZATION SETTING
        NOPTIM=N1ST*N2ND*N3RD+NADD     ! = NRAD ??
C       NOPTM1=   ??
        NOPTIM_IN = NOPTIM
      ENDIF
C  TRY TO BE INTELLIGENT: AUTOMATIC STORAGE REDUCTION....
C  SWITCH OFF SUM OVER STRATA IF THERE IS ONLY ONE STRATUM TO BE CALCULATED
c  TO BE TESTED, E.G. CHECK VARIANCES, AND VARIANCES SUM OVER STRATA
      IF (NSTRAI == 1.AND.NSMSTRA.NE.0) THEN
         WRITE (iunout,*) 'NSMSTRA SET = 0, BECAUSE NSTRAI=1 '
         NSMSTRA = 0
      ENDIF
      CALL EIRENE_LEER(1)

C  OPTIONAL STORAGE/PERFORMANCE HANDLING FLAGS
      WRITE (iunout,'(a14,i8)') 'NOPTIM      = ',NOPTIM
      WRITE (iunout,'(a14,i8)') 'NOPTM1      = ',NOPTM1
      WRITE (iunout,'(a14,i8)') 'NGEOM_USR   = ',NGEOM_USR
      WRITE (iunout,'(a14,i8)') 'NCOUP_INPUT = ',NCOUP_INPUT
      WRITE (iunout,'(a14,i8)') 'NSMSTRA     = ',NSMSTRA
      WRITE (iunout,'(a14,i8)') 'NSTORAM     = ',NSTORAM
      WRITE (iunout,'(a14,i8)') 'NGSTAL      = ',NGSTAL
C
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) 'SETTING OF CENSUS STORAGE FOR T-DEP. MODE'
cdr  time-dependent options: census array size
      WRITE (iunout,'(a14,i8)') 'NPRNL       = ',NPRNL
C
      CALL EIRENE_LEER(2)
cpg
      call eirene_couple_param_consistency(nlimi,nstsi,textal,
     .                                     size(textal))
cpg
C
      RETURN
C
 6664 FORMAT (6E12.4)
 6665 FORMAT (12(5L1,1X))
 6666 FORMAT (12I6)
C
 6999 WRITE (IUNOUT,*) 'Empty input file found!'
      WRITE (IUNOUT,*)
     . 'Either remove it or replace it with a correct file.'
      CALL EIRENE_EXIT_OWN(1)
 7999 WRITE (IUNOUT,*) 'Could not open input file!'
      CALL EIRENE_EXIT_OWN(1)
      RETURN

      CONTAINS

      SUBROUTINE eirene_Ti_input(indpro2_save,l45)
      implicit none
      integer, intent(out) :: indpro2_save
      logical, intent(in) :: l45

      integer icount, nv, indpro2
cdr Nov. 2019
cdr Achieve synchronisation of input options for multi-species Ti and Vi.
cdr Return a "best guess" of what INDPRO(2) should be.

cdr Ti_IN(ipls) options have historically been just opposite to V_IN(ipls) options.
cdr This has gotten amplified to a long-lasting code mess:
cdr Increasingly inconsistent input options for TIIN(ipls)
cdr and vector V..IN(ipls) profiles
cdr and false code operation (in multi-fluid cases) in many instances.
cdr Now: We try to automatically set the TIIN flag INDPRO(2) by counting
cdr input cards in this input block 5.

      if (NPLS.eq.1) then
c  no ambiguity possible here...
        indpro2_save=iabs(indpro(2))
      endif
      write (iunout,*)
     .    'FIND_PARAM: try to infer INDPRO(2) for Ti profiles'
      write (iunout,*) 'indpro2, npls ',indpro(2),npls

c remove sign and second or third digits
      indpro2=mod(iabs(indpro(2)),10)  ! now within 1 and 9

      if (indpro2.eq.6) then
c  Ti(ipls) is from external file/code.
c  Try to find out what was meant:
        if (indpro(2).lt.0) then       ! indpro(2)=-6, -16,...
cdr  just trying, via CI
          indpro2_save=6  ! maybe this was meant in solps-iter cases?
        elseif (indpro(2).lt.10) then  ! indpro(2)=6
          indpro2_save=16 ! maybe this was meant in solps-iter cases?
        else                           ! indpro(2)=16,26,36,...106,...
cdr probably never used?
          write (iunout,*) 'unknown option for indpro(2)'
          write (iunout,*) 'use indpro(2)=6  for multispec. Ti'
          write (iunout,*) 'use indpro(2)=16 for single Ti for all'
          call eirene_exit_own(1)
        endif
c  no card counting done.
        indpro(2)=indpro2_save
        write (iunout,*) 'indpro(2) reset to: ',indpro(2)
        goto 200
      endif

C  FIND
C        START OF NEXT INPUT BLOCK: 6,
C  OR    OF OPTIONAL INPUT LINES,
C  or,   in case of external block 45, L45: of end of file
cdr  and count the remaining input lines for plasma profiles in block 5.


      icount=0
      DO
        READ (IUNIN,'(A72)',end=100) ZEILE
        IF ((ZEILE(1:3) == '***').OR.
     .      (ZEILE(1:3) == 'OPT')
     .      ) EXIT
        icount=icount+1
      END DO
      goto 101
  100 continue
      write (iunout,*) 'end of block45 found'
  101 write (iunout,*) 'cards found in block 5 ',icount

c  te:
      if (indpro(1).ne.6) icount=icount-1
c  ni
      if (indpro(3).ne.6) icount=icount-npls
c  vx,vy,vz
cdr We trust the indpro(4) input, unchanged since 1985.
cdr Only a single common flow velocity is specified.
      nv=npls
      if (indpro(4).gt.10) nv=1
      if (iabs(indpro(4)).ne.6) icount=icount-3*nv
c  B
      if (iabs(indpro(5)).ne.6) icount=icount-1
c  zi
      if (indpro(11).ne.0.and.indpro(11).ne.6) icount=icount-npls
c  Vol?
cdr optional, we do not know if such a card has been included


  200 continue
      if (npls.eq.1.or.indpro2.eq.6) return

      write (iunout,*) 'cards for Ti (+Vol?) ',icount


      IF (icount+1.lt.npls.and.icount.ge.1) then
c  we necessarily have data for only one single Ti in the input file.
        if (indpro(2).lt.10) then
c         if (indpro(2).lt.0) ...  input is in K rather than eV. Not used.
          indpro(2)=iabs(indpro(2))
          indpro(2)=10+indpro(2)
          write (iunout,*) 'apparently only one Ti profile is given'
          write (iunout,*) 'indpro(2) reset to: ',indpro(2)
        endif
      elseif ((icount.eq.npls .or. icount.eq.npls+1) .and.
     .        icount.ne.2) then
c  we necessarily have npls lines for the npls TI profiles
        if (indpro(2).lt.0.or.indpro(2).gt.10) then
c         if (indpro(2).lt.0) ...  input is in K rather than eV. Not used.
          indpro(2)=iabs(indpro(2))
          indpro(2)=mod(indpro(2),10)
          write (iunout,*) 'apparently npls Ti profiles are given'
          write (iunout,*) 'indpro(2) reset to: ',indpro(2)
        endif
      else
c  No unique decision possible. Can only happen in case NPLS=2.
c  Due to the volume tally input card (indpro(12)) being optional.
        if (npls.ne.2) then
          write (iunout,*) 'code confused wrt. multispec. Ti, Vi input'
          call eirene_exit_own(1)
        endif
        write (iunout,*) 'Unclear setting of INDPRO(2). NPLS: ',NPLS
        write (iunout,*) 'This should only occur in case NPLS=2'
        write (iunout,*) 'No. of input lines for Ti and vol? ',icount
        write (iunout,*) 'Choose INDPRO(2) > 10 (NPLSTI=1)'
        write (iunout,*) 'and assume that a VOL card is present.'
        indpro(2) = mod(indpro(2),10) + 10
        write (iunout,*) 'IF NPLSTI=2 was intended, add the VOL-card'
        write (iunout,*) 'i.e. add a card, e.g. containing just 0.0'
      endif
      indpro2_save=indpro(2)
      return
      end subroutine eirene_Ti_input

      END SUBROUTINE EIRENE_FIND_PARAM_FIXFORM
