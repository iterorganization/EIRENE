
      SUBROUTINE EIRENE_READ_JSON(NL, SAREA_SAVE, IERROR)
C
C   READ INPUT DATA AND SET DEFAULT VALUES
c   IN CASE IITER.GT.1 OR ITIMV.GT.1 : SKIP READING NEW INPUT FROM IUNIN.
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
      USE EIRMOD_MPI
      USE EIRMOD_PROFILES
      USE EIRMOD_JSON
      USE EIRMOD_PRESSURELOOP

#ifndef NO_JSON
      USE EIRMOD_CPES, ONLY: NPRS, NLIDENT,
     >    STRATEGY_UNDEFINED, STRATEGY_EMBARRASS,
     >    STRATEGY_ORIGINAL, STRATEGY_APCAS,
     >    STRATEGY_BALANCED, STRATEGY_IMPROVED,
     >    INPUT_DISTRIBUTION_STRATEGY
      USE EIRMOD_TIMEA, ONLY: EIRENE_TIMEA0
      USE EIRMOD_SECOND_OWN, ONLY: EIRENE_SECOND_OWN
      USE EIRMOD_IOUSR, ONLY: EIRENE_READ_BLOCK_11_USR
      USE EIRMOD_INFCOP, ONLY: EIRENE_IF0COP
      use json_module           !IGNORE
! this does not work  found no matching specific binding for json%get
!     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck
#endif

      IMPLICIT NONE

C
#ifndef NO_JSON
      TYPE(TEMPERATURE),POINTER :: TEMPCUR
      TYPE(DENSITY),POINTER :: DENCUR
      TYPE(VELOCITY),POINTER :: VELCUR
      TYPE(VOLUMEP),POINTER :: VOLCUR
C
      TYPE(TSURFACE), POINTER :: SURFCUR
      TYPE(REFMODEL), POINTER :: REFCUR, REFLAST
      TYPE(SPEC_REF), POINTER :: SPR, SPR_LAST

      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC, SSPEC
      TYPE(TCONTRIB) :: CNT
cdr  optional arguments for FILNAM=CONST options
      INTERFACE
        SUBROUTINE EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
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
#endif

      INTEGER, INTENT(IN) :: NL
      REAL(DP), INTENT(INOUT) :: SAREA_SAVE(NL)

C
#ifndef NO_JSON
      REAL(DP) :: AFF(3,3), AFFI(3,3)
      REAL(DP) :: ALR, ROTNRM, RPSDL, XSH, YSH, ZSH, REFNRM,
     .          DPP, SPCMN, SPCMX,SPC_SHIFT,
     .          SPCPLT_X,SPCPLT_Y,SPCPLT_SAME, SPCVX, SPCVY, SPCVZ,
     .          VNORM
#endif

      INTEGER :: IERROR, IREAD
c
C  MULTIPLIER FOR BOTH CPU TIME NTCPU AND MAX NUMBER OF MC HISTORIES NPTS, ....
#ifndef NO_JSON
      REAL(DP) :: AMPTS

      INTEGER :: II, NTLVF, NSRF, NSP, NTL, ICHORI,
     .           JSTRA, I2, NTIME0, I, J, IMSG,
     .           K, NRGEN, NTLV, MXL,
     .           IR, IPLN,
     .           ISPSRF, ISPTYP, NSPS, NSPSA, IPTYP, IPSPZ, ISRFCLL,
     .           IDIREC, IL, ISPOPT, ISPLDEG,
     .           ISTR, IADV, IDMDL,
     .           JS, KK, ELONLY

      INTEGER, SAVE :: NITER0
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LRDMLTI, LRDMLV
#else
      EXTERNAL :: EIRENE_MASAGE
#endif
      LOGICAL :: LHYDDEF, LMULPL
      EXTERNAL :: EIRENE_EXIT_OWN

C
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST ITERATION
C  STEP IN THIS RUN. IITER IS THE ACTUAL ITERATION NUMBER
C  DO NOT READ ANY INPUT, IF THIS IS NOT THE VERY FIRST TIMESTEP
C  IN THIS RUN. ITIMV IS THE ACTUAL TIMESTEP NUMBER
C
C  INITIALISE SOME DATA AND SET DEFAULTS
C
      IREAD=0
      IERROR=0
C
!pb IUNIN set in COMPRT
C
C  UNIT NUMBER FOR OUTPUT FILE: MUST BE DIFFERENT FROM: 5,8,10,11,12
C  13,14, AND 15 AND IUNIN
cdr IUNOUT has already been setin COMPRT and in subroutine EIRENE
!pb   IUNOUT=6+ifoff
C
      IUSROUT = 0

      NULLIFY(SURFLIST)
      NULLIFY(REFLIST)
C
C  SET DEFAULT DATA FOR BLOCK 13
C
      DTIMV=1.D30
      TIME0=0.
      NSNVI=0
      NTMSTP=1

      LMULPL = .FALSE.
#ifndef NO_JSON
      js = itree_num(0)
      CALL EIRENE_READ_BLOCK_0(jtrees(js),blks(0)%p)

      js = itree_num(1)
      CALL EIRENE_READ_BLOCK_1(jtrees(js),blks(1)%p)

      js = itree_num(2)
      CALL EIRENE_READ_BLOCK_2(jtrees(js),blks(2)%p)

      js = itree_num(3)
      CALL EIRENE_READ_BLOCK_3a(jtrees(js),blks(3)%p)

      js = itree_num(4)
      CALL EIRENE_READ_BLOCK_3b(jtrees(js),blks(4)%p)

      js = itree_num(5)
      CALL EIRENE_READ_BLOCK_4(jtrees(js),blks(5)%p)
      CALL EIRENE_READ_BLOCK_5(jtrees(js),blks(5)%p)

      js = itree_num(6)
      CALL EIRENE_READ_BLOCK_6(jtrees(js),blks(6)%p)

      js = itree_num(7)
      CALL EIRENE_READ_BLOCK_7(jtrees(js),blks(7)%p)

      js = itree_num(8)
      CALL EIRENE_READ_BLOCK_8(jtrees(js),blks(8)%p)

      js = itree_num(9)
      CALL EIRENE_READ_BLOCK_9(jtrees(js),blks(9)%p)

      js = itree_num(10)
      CALL EIRENE_READ_BLOCK_10(jtrees(js),blks(10)%p)

      js = itree_num(11)
      CALL EIRENE_READ_BLOCK_11(jtrees(js),blks(11)%p)

      js = itree_num(12)
      CALL EIRENE_READ_BLOCK_12(jtrees(js),blks(12)%p)

      js = itree_num(13)
      CALL EIRENE_READ_BLOCK_13(jtrees(js),blks(13)%p)
#else
      CALL EIRENE_MASAGE
     . ('JSON FORMAT NOT SUPPORTED WHEN COMPILED WITH NO_JSON OPTION')
      CALL EIRENE_EXIT_OWN(1)
      LHYDDEF = .FALSE.
#endif

C  TRY TO GENERATE INPUT FOR EIRENE CORRESPONDING TO A HYDKIN RUN.
C  NOT READY
      IF (LHYDDEF) THEN
        WRITE (IUNOUT,*)
     .    'LHYDDEF IS TRUE: OPTION NOT READY, EXIT CALLED'
        CALL EIRENE_EXIT_OWN(1)
cdr  subr. EIRENE_SETUP_HYDKIN_REACTIONS moved to folder: unfinished_business
cdr     CALL EIRENE_SETUP_HYDKIN_REACTIONS(HYDKIN_DEFAULT,CADAPT)
      ENDIF

      RETURN

#ifndef NO_JSON
      CONTAINS

!*****************************************************************

      subroutine eirene_read_block_0(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      integer :: ncm, k
      logical :: found, found_cm
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txtr
      character(kind=json_CK,len=420), dimension(:), allocatable ::
     .           cmlines
!cym/cpg end

      call json%get(p,'TXTRUN',txtr,found)
      TXTRUN = trim(txtr)
      WRITE (iunout,'(1X,A)') trim(TXTRUN)

      call json%get(p,'COMMENTS',cmlines,found_cm)
      if (found_cm) then
        ncm = size(cmlines)
        do k = 1, ncm
          WRITE (IUNOUT,'(1X,A)') TRIM(CMLINES(K))
!  STORE CM LINES FOR OUTPUTING TO JSON FILE
          call eirene_push_string_stack(cm_stack,cmlines(k))
        end do
        deallocate (cmlines)
      end if

      return
      end subroutine eirene_read_block_0

!*****************************************************************

      subroutine eirene_read_block_1(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pfile, pdb
      logical :: found, NLSPCSCL_DUM, NLSPCSCL_ON_DUM, NLSPCSPT_DUM
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: cdbh, cdbf
      character(len=:), allocatable :: cdbhc
!cym/cpg end
      integer :: ifile, nch
      external :: eirene_filepath_usr, eirene_leer, eirene_masage

      call json%get(p,'NPRLL',nprll,found)
      call json%get(p,'NMODE',nmode,found)
      call json%get(p,'NTCPU',ntcpu,found)
      call json%get(p,'NFILE',nfile,found)
      call json%get(p,'NITER0',niter0,found)
      call json%get(p,'NITER',niter,found)
      call json%get(p,'NTIME0',ntime0,found)
      call json%get(p,'NTIME',ntime,found)
      NTIME_IN = NTIME
      call json%get(p,'NEVNTMAX',nevntmax,found)
      NEVNTMAX=1000*NEVNTMAX
      call json%get(p,'NVNMXTR',nvnmxtr,found)
      call json%get(p,'NSTORAM',nstoram,found)
      call json%get(p,'NGSTAL',ngstal,found)

      CALL EIRENE_LEER(1)
      IITER=MAX0(1,NITER0)
      ITIMV=MAX0(1,NTIME0)

!  NOPTIM, NOPTM1 ... already read in find_param_json

      call json%get(p,'NLSCL',nlscl,found)
      call json%get(p,'NLTEST',nltest,found)
      call json%get(p,'NLANA',nlana,found)
      call json%get(p,'NLDRFT',nldrft,found)
      call json%get(p,'NLCRR',nlcrr,found)

      call json%get(p,'NLERG',nlerg,found)
      call json%get(p,'NLIDENT',nlident,found)
      call json%get(p,'NLONE',nlone,found)
      call json%get(p,'NLMOVIE',nlmovie,found)
      call json%get(p,'NLDFST',nldfst,found)

      call json%get(p,'NLRANMAR',nlranmar,found)
      call json%get(p,'NLCASCAD',nlcascad,found)
      call json%get(p,'NLOCTREE',nloctree,found)
      call json%get(p,'NLWRMSH',nlwrmsh,found)
      call json%get(p,'NEXVS',nexvs,found)
      call json%get(p,'NLTRIMESH',nltrimesh,found)
!pb do not overwrite
      call json%get(p,'NLSPCSCL',nlspcscl_dum,found)
      call json%get(p,'NLSPCSCL_ON',nlspcscl_on_dum,found)
      call json%get(p,'NLSOLEDGE',nlsoledge,found)
      call json%get(p,'NLSPCSPT',nlspcspt_dum,found)

      call json%get_child(p,'CFILE',pfile,found)
      if (found) then
        call json%info(pfile,n_children=nch)
!       do ifile = 1, ndbnames

        do i = 1, nch
          call json%get_child(pfile,i,pdb,found)
          call json%get(pdb,'FILE',cdbh,found)
          call json%get(pdb,'PATH',cdbf,found)
!cym/cpg
          cdbhc=cdbh
!cym/cpg end
          DO IFILE = 1,NDBNAMES
!cym/cpg
!            IF (INDEX(DBHANDLE(IFILE),cdbh) /= 0) EXIT
             IF (INDEX(DBHANDLE(IFILE),cdbhc) /= 0) EXIT
!cym/cpg end
          END DO
          IF (IFILE <= NDBNAMES) THEN
            dbhandle(ifile) = trim(cdbh)
            dbfname(ifile) = trim(cdbf)
            dbfname_in(ifile) = trim(cdbf)
            call eirene_filepath_usr(cdbf,dbfname(ifile),
     .                               1,len_trim(cdbf))
            ldbread(ifile) = .true.
            WRITE (IUNOUT,'(2A)') 'PATH SET FOR FILE ',dbhandle(ifile)
            WRITE (IUNOUT,'(2A)') 'PATH = ',trim(dbfname(ifile))
          ELSE
            WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED'
            WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',
     .                         trim(cdbh),' IGNORED'
          end if
          deallocate(cdbh)
          deallocate(cdbf)
          nullify(pdb)
        end do
      end if

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
     .  ('       EIRENE RUN AS STAND ALONE CODE')
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
        NFILEL = NFILEL - 5
        NLSHRT13 = .TRUE.
      END IF
      SELECT CASE (NFILEL)
        CASE (1)
          IF (NLSHRT13) THEN
            WRITE (iunout,*) '       EIRENE SAVES ONLY PLASMA DATA,'
            WRITE (iunout,*) '       STARTING FROM SPECIES NUMBER'
            WRITE (iunout,*) '       NPLS_FIX+1 (COUPLING'
            WRITE (iunout,*) '       VARIABLE, BLOCK 14),'
          ELSE
            WRITE (iunout,*) '       EIRENE SAVES PLASMA DATA, A&M DATA'
            WRITE (iunout,*) '       AND SOURCE DISTRIBUTION DATA'
          END IF
          WRITE (iunout,*) '       ON FILE FT13 AT END OF RUN, I.E.,'
          WRITE (iunout,*) '       AFTER LAST TIMESTEP OR ITERATION'
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
        WRITE (iunout,*) '       END OF LAST TIMESTEP ON FILE FT15'
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
c  printout:  which default storage settings have we overruled (consequences thereof?)
C
      return
      end subroutine eirene_read_block_1

!*****************************************************************

      subroutine eirene_read_block_2(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: prad, ppol, ptor, pnts, plgs, child,
     .                             pmlt, padd
      type(json_core) :: core
      logical :: found
      integer :: ilogs, i, k, nch
      integer, allocatable :: ihelp(:)
      real(dp), allocatable :: rhelp(:)
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: clab
!cym/cpg end
      external :: eirene_read_tetra, eirene_read_triang,
     .            eirene_leer, eirene_masage, eirene_exit_own

      CALL EIRENE_MASAGE
     .  ('*** 2. DATA FOR STANDARD MESH')
      CALL EIRENE_LEER(1)

      call json%get(p,'INDGRD',ihelp,found)
      indgrd(1:3) = ihelp(1:3)
      deallocate(ihelp)

      call json%get(p,'NLRAD',nlrad,found)
      call json%get(p,'NLPOL',nlpol,found)
      call json%get(p,'NLTOR',nltor,found)
      call json%get(p,'NLMLT',nlmlt,found)
      call json%get(p,'NLADD',nladd,found)
      nlpol_in = nlpol

      IF (NLRAD) THEN
C
        call json%get(p,'RADIAL_GRID',prad,found)

        call json%get(prad,'NLSLB',nlslb,found)
        call json%get(prad,'NLCRC',nlcrc,found)
        call json%get(prad,'NLELL',nlell,found)
        call json%get(prad,'NLTRI',nltri,found)
        call json%get(prad,'NLPLG',nlplg,found)
        call json%get(prad,'NLFEM',nlfem,found)
        call json%get(prad,'NLTET',nltet,found)
        call json%get(prad,'NLGEN',nlgen,found)
        nlplg_in = nlplg
        nlfem_in = nlfem

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
     .                     ' FLAGS MAY BE .TRUE. '
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

        call json%get(prad,'NR1ST',nr1st,found)
        call json%get(prad,'NRSEP',nrsep,found)
        call json%get(prad,'NRPLG',nrplg,found)
        call json%get(prad,'NPPLG',npplg,found)
        call json%get(prad,'NRKNOT',nrknot,found)
        call json%get(prad,'NCOOR',ncoor,found)
        nr1st_in = nr1st

        IF (NR1ST < 0) THEN
          NR1ST = N1ST
        ENDIF

        IF (INDGRD(1).LE.5) THEN
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            call json%get(prad,'RIA',ria,found)
            call json%get(prad,'RGA',rga,found)
            call json%get(prad,'RAA',raa,found)
            call json%get(prad,'RRA',rra,found)
            IF (NLELL.OR.NLTRI) THEN
              call json%get(prad,'EP1IN',ep1in,found)
              call json%get(prad,'EP1OT',ep1ot,found)
              call json%get(prad,'EP1CH',ep1ch,found)
              call json%get(prad,'EXEP1',exep1,found)

              call json%get(prad,'ELLIN',ellin,found)
              call json%get(prad,'ELLOT',ellot,found)
              call json%get(prad,'ELLCH',ellch,found)
              call json%get(prad,'EXELL',exell,found)
              IF (NLTRI) THEN
                call json%get(prad,'TRIIN',triin,found)
                call json%get(prad,'TRIOT',triot,found)
                call json%get(prad,'TRICH',trich,found)
                call json%get(prad,'EXTRI',extri,found)
              ENDIF
            ENDIF
          ENDIF
          IF (NLPLG) THEN
            call json%get(prad,'XPCOR',xpcor,found)
            call json%get(prad,'YPCOR',ypcor,found)
            call json%get(prad,'ZPCOR',zpcor,found)
            call json%get(prad,'PLREFL',plrefl,found)

            call json%get(prad,'POLYGON_PARTS',pnts)
            call json%info(pnts,n_children=nch)
            if (nch < npplg) then
              write (iunout,*) ' NUMBER OF POLYGON PARTS DOES',
     .             ' NOT MATCH NPPLG '
              write (iunout,*) 'NPPLG = ',npplg
              write (iunout,*) 'NCH =   ',nch
              call eirene_exit_own(1)
            end if

            do k = 1, npplg
              call core%get_child(pnts,k,child)
              call core%get(child,'NPOINT',ihelp,found)
              npoint(1:2,k) = ihelp(1:2)
              deallocate(ihelp)
              nullify(child)
            end do
            nullify(pnts)

            call json%get(prad,'POLYGONS',plgs)
            call json%info(plgs,n_children=nch)
            if (nch < nr1st) then
              write (iunout,*) ' NUMBER OF POLYGON PARTS DOES',
     .             ' NOT MATCH NR1ST '
              write (iunout,*) 'NR1ST = ',nr1st
              write (iunout,*) 'NCH =   ',nch
              call eirene_exit_own(1)
            end if
            do i = 1, nr1st
              call core%get_child(plgs,i,child)
              call core%get(child,'XPOL_YPOL',rhelp,found)
              if (.not.found) then
                write (iunout,*) 'ERROR IN READING INPUT BLOCK 2'
                write (iunout,*)
     .        'POLYGON GRID REQUESTED BUT NO DATA FOUND FOR XPOL, YPOL'
                call eirene_exit_own(1)
              end if
              if (size(rhelp) /= 2*nrplg) then
                write (iunout,*) ' NUMBER OF POLYGON NODES DOES',
     .             ' NOT MATCH NRPLG '
                write (iunout,*) 'NRPLG = ',nrplg
                write (iunout,*) 'NCH =   ',nch
                call eirene_exit_own(1)
              end if
              xpol(i,1:nrplg) = rhelp(1:2*nrplg:2)
              ypol(i,1:nrplg) = rhelp(2:2*nrplg:2)
              deallocate(rhelp)
              nullify(child)
            end do
            nullify(plgs)
            IF (PLREFL.GT.0.D0) NR1ST=NR1ST+1
          ENDIF
          IF (NLFEM .OR. NLTET) THEN
            call json%get(prad,'CASE',clab,found)
            if (.not.found) then
              WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
              WRITE (IUNOUT,*)
     .          ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
              WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
              CALL EIRENE_EXIT_OWN(1)
            else
              casename = adjustl(trim(clab))
            end if
            i2 = len_trim(casename)

            IF (NLFEM) CALL EIRENE_READ_TRIANG (CASENAME(1:I2))
            IF (NLTET) CALL EIRENE_READ_TETRA (CASENAME(1:I2))

            call json%get(prad,'XPCOR',xpcor,found)
            call json%get(prad,'YPCOR',ypcor,found)
            call json%get(prad,'ZPCOR',zpcor,found)

          END IF
          IF (NLGEN) THEN
            NRGEN=NR1ST
          ENDIF
        ELSEIF (INDGRD(1).EQ.6) THEN
!  NLPOL HAS ALREADY BEEN READ
          IF (NLSLB.OR.NLCRC.OR.NLELL.OR.NLTRI) THEN
            call json%get(prad,'RIA',ria,found)
            call json%get(prad,'RGA',rga,found)
            call json%get(prad,'RAA',raa,found)
            IREAD = 0
          ELSEIF (NLPLG.OR.NLFEM.OR.NLTET) THEN
            call json%get(prad,'XPCOR',xpcor,found)
            call json%get(prad,'YPCOR',ypcor,found)
            call json%get(prad,'ZPCOR',zpcor,found)
            IREAD = 0
          ENDIF
        ENDIF
      ENDIF

C
C  POLOIDAL MESH
C
C INPUT SUB-BLOCK 2B
C
      call json%get(p,'POLOIDAL_GRID',ppol,found)

      call json%get(ppol,'NLPLY',nlply,found)
      call json%get(ppol,'NLPLA',nlpla,found)
      call json%get(ppol,'NLPLP',nlplp,found)

      call json%get(ppol,'NP2ND',np2nd,found)
      call json%get(ppol,'NPSEP',npsep,found)
      call json%get(ppol,'NPPLA',nppla,found)
      call json%get(ppol,'NPPER',npper,found)
      np2nd_in = np2nd

      IF (INDGRD(2).LE.5) THEN
        call json%get(ppol,'YIA',yia,found)
        call json%get(ppol,'YGA',yga,found)
        call json%get(ppol,'YAA',yaa,found)
        call json%get(ppol,'YYA',yya,found)
      ELSEIF (INDGRD(2).EQ.6) THEN
      ENDIF
C
C  TOROIDAL MESH
C
C INPUT SUB-BLOCK 2C
C
!  NLTOR HAS ALREADY BEEN READ
C
      call json%get(p,'TOROIDAL_GRID',ptor,found)

      call json%get(ptor,'NLTRZ',nltrz,found)
      call json%get(ptor,'NLTRA',nltra,found)
      call json%get(ptor,'NLTRT',nltrt,found)

      call json%get(ptor,'NT3RD',nt3rd,found)
      call json%get(ptor,'NTSEP',ntsep,found)
      call json%get(ptor,'NTTRA',nttra,found)
      call json%get(ptor,'NTPER',ntper,found)
      nt3rd_in = nt3rd

      IF (INDGRD(3).LE.5) THEN
        call json%get(ptor,'ZIA',zia,found)
        call json%get(ptor,'ZGA',zga,found)
        call json%get(ptor,'ZAA',zaa,found)
        call json%get(ptor,'ZZA',zza,found)
        call json%get(ptor,'ROA',roa,found)
      ELSEIF (INDGRD(3).EQ.6) THEN
      ENDIF
C
C  MESH MULTIPLICATION
C
C INPUT SUB-BLOCK 2D
C
!  NLMLT HAS ALREADY BEEN READ
C
      IF (NLMLT) THEN
        call json%get(p,'MESH_MULTIPLICATION',pmlt,found)

        call json%get(pmlt,'NBMLT',nbmlt,found)
        call json%get(pmlt,'VOLCOR',rhelp,found)
        volcor(1:nbmlt) = rhelp(1:nbmlt)
        deallocate(rhelp)
      ELSE
        NBMLT=1
        VOLCOR(1)=1.D0
      ENDIF
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
!  NLADD HAS ALREADY BEEN READ
C
      IF (NLADD) THEN
        call json%get(p,'ADD_CELLS',padd,found)

        call json%get(padd,'NRADD',nradd,found)
        call json%get(padd,'VOLADD',rhelp,found)
        voladd(1:nradd) = rhelp(1:nradd)
        deallocate(rhelp)
      ENDIF

      nullify(prad)
      nullify(ppol)
      nullify(ptor)
      nullify(pmlt)
      nullify(padd)
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
        IF (NLPOL.AND.NRPLG.NE.NP2ND) THEN
          WRITE (iunout,*) 'ERROR IN INPUT: ' //
     .                     'NRPLG.NE.NP2ND, BUT NLPOL=TRUE'
          WRITE (iunout,*) 'NRPLG,NP2ND ',NRPLG,NP2ND
          CALL EIRENE_EXIT_OWN(1)
        END IF
      ENDIF

      return
      end subroutine eirene_read_block_2

!*****************************************************************

      subroutine eirene_read_block_3a(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pnstd, psrf
      type(json_core) :: core
      integer :: ists, nlj, idimp, IRPTA1, IRPTA2, IRPTA3,
     .           IRPTE1, IRPTE2, IRPTE3
      logical :: found
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txts
!cym/cpg end
      external :: eirene_leer, eirene_masage, eirene_exit_own
C
      IF (IREAD.EQ.0) READ (IUNIN,*)
      CALL EIRENE_MASAGE
     .  ('*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES')

      call json%get(p,'NSTSI',nstsi,found)

      WRITE (iunout,*) '        NSTSI= ',NSTSI
      CALL EIRENE_LEER(1)
      IMSG=0

      call json%get(p,'SURFACES',pnstd)

      DO ISTS=1,NSTSI
        NLJ=NLIM+ISTS
        call core%get_child(pnstd,ists,psrf)

        call json%get(psrf,'TXTSFL',txts,found)
        TXTSFL(NLJ) = trim(txts)
        deallocate(txts)

        call json%get(psrf,'IDIMP',idimp,found)
        call json%get(psrf,'INUMP',inump(ists,idimp),found)
        call json%get(psrf,'IRPTA1',irpta1,found)
        call json%get(psrf,'IRPTE1',irpte1,found)
        call json%get(psrf,'IRPTA2',irpta2,found)
        call json%get(psrf,'IRPTE2',irpte2,found)
        call json%get(psrf,'IRPTA3',irpta3,found)
        call json%get(psrf,'IRPTE3',irpte3,found)

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
          ENDIF
        ELSE ! IDIMP .EQ. 3
          IRPTA(ISTS,3) = IRPTA3
          IRPTE(ISTS,3) = IRPTE3
        ENDIF

        call eirene_read_surf_switches(json,psrf,nlj)

        call eirene_read_ref_model(json,psrf,nlj)

        IF (IMSG.EQ.1) CALL EIRENE_LEER(1)

      end do

      return
      end subroutine eirene_read_block_3a

!*****************************************************************

      subroutine eirene_read_surf_switches (json,srf,nlj)

      class(json_core),intent(inout) :: json
      type(json_value),pointer :: srf
      integer, intent(in) :: nlj
      logical :: found

      call json%get(srf,'ILIIN',iliin(nlj),found)
      call json%get(srf,'ILSIDE',ilside(nlj),found)
      call json%get(srf,'ILSWCH',ilswch(nlj),found)
      call json%get(srf,'ILEQUI',ilequi(nlj),found)
      call json%get(srf,'ILTOR',iltor(nlj),found)
      call json%get(srf,'ILCOL',ilcol(nlj),found)
      call json%get(srf,'ILFIT',ilfit(nlj),found)
      call json%get(srf,'ILCELL',ilcell(nlj),found)
      call json%get(srf,'ILBOX',ilbox(nlj),found)
      call json%get(srf,'ILPLG',ilplg(nlj),found)
      ilcol_in(nlj) = ilcol(nlj)

      IF (ABS(ILCOL(NLJ)).EQ.7) THEN
        WRITE (iunout,*) 'COLOUR FLAG ILCOL CHANGED FOR SURFACE NO. ',
     .                      NLJ
        WRITE (iunout,*) 'COLOUR NO. 7 IS RESERVED FOR "NON-ANALOGUE'
        WRITE (iunout,*)
     .      'SURFACES" (SPLITTING, R.R., WEIGHT WINDOWS,..)'
        ILCOL(NLJ)=ILCOL(NLJ)-2
        IMSG=1
      ENDIF

      return
      end subroutine eirene_read_surf_switches

!*****************************************************************

      subroutine eirene_read_ref_model (json,srf,nlj)

      class(json_core),intent(inout) :: json
      type(json_value),pointer :: srf
      integer, intent(in) :: nlj
      integer :: ilr
      logical :: found, foundm, foundr, founds(5)
      real(dp), allocatable :: rhelp(:)
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: smod
!cym/cpg end

      call json%get(srf,'SURFMOD', smod, foundm)

      if (foundm) then
!  surface model specified
        ALLOCATE(SURFCUR)
        SURFCUR%MODNAME = TRIM(ADJUSTL(SMOD(9:)))
        SURFCUR%NOSURF = NLJ
        SURFCUR%NEXT => SURFLIST
        SURFLIST => SURFCUR
        SMOD_NAME(NLJ) = SURFCUR%MODNAME
      else
!  check for local reflection model directly specified with the surface
        call json%get(srf,'ILREF',ilr,foundr)
        if (foundr) then
          ilref(1,nlj) = ilr
          call json%get(srf,'ILSPT',ilspt(nlj),found)
          call json%get(srf,'ISRS',isrs(1,nlj),found)
          call json%get(srf,'ISRC',isrc(1,nlj),found)
          call json%get(srf,'LCHSPNWL',lchspnwl(1,nlj),found)
          call json%get(srf,'ZNML',znml(nlj),found)
          call json%get(srf,'EWALL',ewall(nlj),found)
          call json%get(srf,'EWBIN',ewbin(nlj),found)
          call json%get(srf,'TRANSP',rhelp,found)
          transp(1,1:2,nlj) = rhelp(1:2)
          deallocate (rhelp)
          call json%get(srf,'FSHEAT',fsheat(nlj),found)
          call json%get(srf,'RECYCF',recycf(1,nlj),found)
          call json%get(srf,'RECYCT',recyct(1,nlj),found)
          call json%get(srf,'RECPRM',recprm(1,nlj),found)
          call json%get(srf,'EXPPL',exppl(1,nlj),found)
          call json%get(srf,'EXPEL',expel(1,nlj),found)
          call json%get(srf,'EXPIL',expil(1,nlj),found)
          call json%get(srf,'RECYCS',recycs(1,nlj),founds(1))
          call json%get(srf,'RECYCC',recycc(1,nlj),founds(2))
          call json%get(srf,'SPTPRM',sptprm(1,nlj),founds(3))
          call json%get(srf,'ESPUTS',esputs(1,nlj),founds(4))
          call json%get(srf,'ESPUTC',esputc(1,nlj),founds(5))

          IF (ALL(FOUNDS)) LCHSPNWL(1,NLJ)=1
          IF ((ILSPT(NLJ).NE.0) .AND. .NOT.ALL(FOUNDS)) THEN
            IF (NLJ > NLIM) THEN
              WRITE (iunout,*)
     .         'WARNING: SPUTTERING AT NON-DEF. SURFACE ',
     .         NLJ-NLIM+1
            ELSE
              WRITE (iunout,*)
     .         'WARNING: SPUTTERING AT ADDITIONAL SURFACE ', NLJ
            END IF
            WRITE (iunout,*)
     .        'BUT NO PARAMETERS RECYCS, RECYCC ARE READ '
            WRITE (iunout,*)
     .        'DEFAULT MODEL: "NO SPUTTERING" IS USED.'
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            ILSPT(NLJ)=0
            IMSG=1
          ENDIF

        end if
      end if

      return
      end subroutine eirene_read_ref_model


!*****************************************************************

      subroutine eirene_read_block_3b(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: padd, psrf,
     .                             pln, pscnd, pchld, ptrns
      type(transform), pointer :: trf
      integer :: nch0, nch1, nch2, ntrf, i, j, k, ih, ihelp(nlimps), nch
      real(dp), allocatable :: rh(:)
      logical :: found, found_ch0, found_ch1, found_ch2, foundt,
     .           foundl, founds
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txts
      character(kind=json_CK,len=72), dimension(:), allocatable ::
     .           ch0lines, ch1lines, ch2lines
      external :: eirene_dekey, eirene_dekeyb, eirene_rotadd,
     .            eirene_setref, eirene_setrot,
     .            eirene_xshadd, eirene_yshadd, eirene_zshadd,
     .            eirene_leer, eirene_masage, eirene_exit_own
!cym/cpg end
C
C  * 3B: READ DATA FOR ADDITIONAL SURFACES 350--399
C
      CALL EIRENE_MASAGE
     .  ('*** 3B. DATA FOR ADDITIONAL SURFACES')

      call json%get(p,'NLIMI',nlimi,found)
      WRITE (iunout,*) '        NLIMI= ',NLIMI
      CALL EIRENE_LEER(1)

      IF (NLIMI == 0) RETURN
C
!  READ CH0 LINES
      IF (NLIMI.GT.0) THEN
        IHELP(1:NLIMI) = IGJUM0(1:NLIMI)

        call json%get(p,'CH0-LINES',ch0lines,found_ch0)
        if (found_ch0) then
          nch0 = size(ch0lines)
          do k = 1, nch0
!  STORE CH0 LINES FOR OUTPUTING TO JSON FILE
            call eirene_push_string_stack(ch0_stack,ch0lines(k))
            CALL EIRENE_DEKEY (ch0lines(k)(4:72),IHELP,1,1,1,NLIMPS)
          end do
          deallocate (ch0lines)
        ENDIF

        IGJUM0(1:NLIMI) = IHELP(1:NLIMI)
      ENDIF
C
      call json%get(p,'SURFACES',padd)
      call json%info(padd,n_children=nch)
      if (nch < nlimi) then
        write (iunout,*) ' NUMBER OF ADD. SURFACES SPECIFIED DOES',
     .          ' NOT MATCH NUMBER OF SURFACES FOUND IN FILE '
        write (iunout,*) 'NLIMI = ',nlimi
        write (iunout,*) 'NCH =   ',nch
        call eirene_exit_own(1)
      end if

      DO I=1,NLIMI
        call json%get_child(padd,i,psrf)

        call json%get(psrf,'TXTSFL',txts,found)
        TXTSFL(I) = trim(txts)
        WRITE (iunout,'(1X,A)') trim(TXTSFL(I))
        IREAD=0
C
C   GENERAL SURFACE DATA
C
!   READ CH1 LINES

        call json%get(psrf,'CH1-LINES',ch1lines,found_ch1)
        if (found_ch1) then
          nch1 = size(ch1lines)
          do k = 1, nch1
!  STORE CH1 LINES FOR OUTPUTING TO JSON FILE
            call eirene_push_string_stack(ch1_stack(i),ch1lines(k))
            IF (NLIMPB >= NLIMPS) THEN
              CALL EIRENE_DEKEY
     .             (ch1lines(k)(4:72),IGJUM1,0,NLIMPS,I,NLIMPS)
            ELSE
              CALL EIRENE_DEKEYB
     .             (ch1lines(k)(4:72),IGJUM1,0,NLIMPS,I,NLIMPB,NBITS)
            END IF
          end do
          deallocate(ch1lines)
        ENDIF

!   READ CH2 LINES

        call json%get(psrf,'CH2-LINES',ch2lines,found_ch2)
        if (found_ch2) then
          nch2 = size(ch2lines)
          do k = 1, nch2
!  STORE CH2 LINES FOR OUTPUTING TO JSON FILE
            call eirene_push_string_stack(ch2_stack(i),ch2lines(k))
            IF (NLIMPB >= NLIMPS) THEN
              CALL EIRENE_DEKEY
     .             (ch2lines(k)(4:72),IGJUM2,0,NLIMPS,I,NLIMPS)
            ELSE
              CALL EIRENE_DEKEYB
     .             (ch2lines(k)(4:72),IGJUM2,0,NLIMPS,I,NLIMPB,NBITS)
            END IF
          end do
          deallocate(ch2lines)
        ENDIF

        call json%get(psrf,'RLB',rlb(i),found)
        call json%get(psrf,'SAREA',sarea_save(i),found)
        call json%get(psrf,'RLWMN',rlwmn(i),found)
        call json%get(psrf,'RLWMX',rlwmx(i),found)
        IF (SAREA_SAVE(I).LE.0.D0) SAREA_SAVE(I)=666.

        call eirene_read_surf_switches(json,psrf,i)

C  READ SURFACE COEFFICIENTS
        IF (RLB(I) < 2.) THEN
          call json%get(psrf,'A0LM',a0lm(i),found)
          call json%get(psrf,'A1LM',a1lm(i),found)
          call json%get(psrf,'A2LM',a2lm(i),found)
          call json%get(psrf,'A3LM',a3lm(i),found)
          call json%get(psrf,'A4LM',a4lm(i),found)
          call json%get(psrf,'A5LM',a5lm(i),found)
          call json%get(psrf,'A6LM',a6lm(i),found)
          call json%get(psrf,'A7LM',a7lm(i),found)
          call json%get(psrf,'A8LM',a8lm(i),found)
          call json%get(psrf,'A9LM',a9lm(i),found)
        ELSEIF (RLB(I) < 7.) THEN
! 2 point surface: straight line
          if (rlb(i) >= 2.) then
            call json%get(psrf,'P1',rh,found)
            p1(1:3,i) = rh
            deallocate(rh)
            call json%get(psrf,'P2',rh,found)
            p2(1:3,i) = rh
            deallocate(rh)
          end if
! 3 point surface: triangle
          if (rlb(i) >= 3.) then
            call json%get(psrf,'P3',rh,found)
            p3(1:3,i) = rh
            deallocate(rh)
          end if
! 4 point surface: flat quadrangle
          if (rlb(i) >= 4.) then
            call json%get(psrf,'P4',rh,found)
            p4(1:3,i) = rh
            deallocate(rh)
          end if
! 5 point surface
          if (rlb(i) >= 5.) then
            call json%get(psrf,'P5',rh,found)
            p5(1:3,i) = rh
            deallocate(rh)
          end if
! 6 point surface
          if (rlb(i) >= 6.) then
            call json%get(psrf,'P6',rh,found)
            p6(1:3,i) = rh
            deallocate(rh)
          end if
        ENDIF
C  READ BOUNDARY DATA
        IF ((RLB(I) > 0.).AND.(RLB(I) < 2.)) THEN
          call json%get(psrf,'XLIMS1',xlims1(1,i),found)
          call json%get(psrf,'YLIMS1',ylims1(1,i),found)
          call json%get(psrf,'ZLIMS1',zlims1(1,i),found)
          call json%get(psrf,'XLIMS2',xlims2(1,i),found)
          call json%get(psrf,'YLIMS2',ylims2(1,i),found)
          call json%get(psrf,'ZLIMS2',zlims2(1,i),found)
        ELSEIF (RLB(I).LE.0.D0) THEN
          IH=NINT(-RLB(I))
          ILIN(I)=EIRENE_IDEZ(IH,1,2)
          ISCN(I)=EIRENE_IDEZ(IH,2,2)
          call json%get(psrf,'LIN_BOUNDS',pln,foundl)
          if (foundl) then
            DO J=1,ILIN(I)
              call json%get_child(pln,j,pchld,found)
              call json%get(pchld,'ALIMS',alims(j,i),found)
              call json%get(pchld,'XLIMS',xlims(j,i),found)
              call json%get(pchld,'YLIMS',ylims(j,i),found)
              call json%get(pchld,'ZLIMS',zlims(j,i),found)
              nullify(pchld)
            END DO
          end if
          nullify(pln)
          call json%get(psrf,'SCN_BOUNDS',pscnd,founds)
          if (founds) then
            DO J=1,ISCN(I)
              call json%get_child(pscnd,j,pchld,found)
              call json%get(pchld,'ALIMS0',alims0(j,i),found)
              call json%get(pchld,'XLIMS1',xlims1(j,i),found)
              call json%get(pchld,'YLIMS1',ylims1(j,i),found)
              call json%get(pchld,'ZLIMS1',zlims1(j,i),found)
              call json%get(pchld,'XLIMS2',xlims2(j,i),found)
              call json%get(pchld,'YLIMS2',ylims2(j,i),found)
              call json%get(pchld,'ZLIMS2',zlims2(j,i),found)
              call json%get(pchld,'XLIMS3',xlims3(j,i),found)
              call json%get(pchld,'YLIMS3',ylims3(j,i),found)
              call json%get(pchld,'ZLIMS3',zlims3(j,i),found)
              nullify(pchld)
            ENDDO
          end if
          nullify(pscnd)
        ENDIF

!  STORE VERTICES OF ADDITIONAL SURFACES FOR OUTPUTING THEM ONTO JSON FILE
!  ORIGINAL VALUES MIGHT BE CHANGED LATER ON

        call eirene_copy_addsrf(i)
C
C READ LOCAL SURFACE INTERACTION MODEL FOR SURFACE NO. I
C
        call eirene_read_ref_model(json,psrf,i)
C
C
C  TRANSFORM FROM CONVENIENT COORDINATE SYSTEM TO EIRENE COORDINATE
C  SYSTEM. THIS IS POSSIBLE FOR ALL SURFACES READ BY NOW, I.E. FROM
C  ITINI=1 TO ITEND=I
       call json%get(psrf,'TRANSFORMATION',ptrns,foundt)
       if (foundt) then
          call json%info(ptrns,n_children=ntrf)
          do j = 1, ntrf
            call json%get_child(ptrns,j,pchld,found)
            call json%get(pchld,'ITINI',itini(i),found)
            call json%get(pchld,'ITEND',itend(i),found)
            IF (ITINI(I).LT.1) ITINI(I)=I
            IF (ITEND(I).GT.I) ITEND(I)=I
            call json%get(pchld,'XLCOR',xlcor(i),found)
            call json%get(pchld,'YLCOR',ylcor(i),found)
            call json%get(pchld,'ZLCOR',zlcor(i),found)
            call json%get(pchld,'XLREF',xlref(i),found)
            call json%get(pchld,'YLREF',ylref(i),found)
            call json%get(pchld,'ZLREF',zlref(i),found)
            call json%get(pchld,'XLROT',xlrot(i),found)
            call json%get(pchld,'YLROT',ylrot(i),found)
            call json%get(pchld,'ZLROT',zlrot(i),found)
            call json%get(pchld,'ALROT',alrot(i),found)
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
cC  SHIFT
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
C     ROTATION
            ROTNRM=XLROT(I)*XLROT(I)+YLROT(I)*YLROT(I)+ZLROT(I)*ZLROT(I)
            ALR=ALROT(I)
            IF (ALR.NE.0..AND.ROTNRM.GT.EPS10) THEN
              CALL EIRENE_SETROT(AFF,AFFI,1,XLROT(I),YLROT(I),ZLROT(I),
     .                           ALROT(I))
              CALL EIRENE_ROTADD(AFF,AFFI,ITINI(I),ITEND(I))
            ENDIF
            nullify (pchld)
          END DO
        ENDIF
        nullify (ptrns)
C
      ENDDO

      CALL EIRENE_LEER(1)

      nullify(padd)
      nullify(psrf)
      return
      end subroutine eirene_read_block_3b


!*****************************************************************

      subroutine eirene_read_block_4(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: preacs, prea, pspc, preac,
     .                             patm, pmol, pion, pphot
      logical :: found, foundr, found_crs
      integer :: nrea, ir, mp, mt, iz,
     .           jfex1mn, jfex1mx, jfex2mn, jfex2mx, iftfl, ncoef,
     .           ncrs, k, ll, m_qext, m_popesc, m_upper, m_lower, 
     .           iprftype
      real(dp) :: r1mn, r1mx, r2mn, r2mx, fp1(6), fp2(6), coef(9),
     .            RC1MIN, RC1MAX, RC2MIN, RC2MAX
      real(dp), allocatable :: ffp1(:), ffp2(:), ccoef(:),
     .                         pop_esc(:), pesc(:)
      integer, allocatable :: irow_esc(:), icol_esc(:),
     .                        irow_ext(:), ital_ext_1(:),
     .                        ital_ext_2(:), ital_ext_3(:)
      integer, allocatable :: ital_ext(:,:)
!cym/cpg : changed type
      character(kind=json_CK,len=:), allocatable :: chr
      character(kind=json_CK,len=8) :: filnam
      character(kind=json_CK,len=4) :: h123
      character(kind=json_CK,len=50) :: reac
      character(kind=json_CK,len=2) :: elname
      character(kind=json_CK,len=3) :: crc
      character(kind=json_CK,len=60) :: bundling

      character(len=8) :: filnamc
      character(len=4) :: h123c
      character(len=50) :: reacc
      character(len=3) :: crcc
      character(len=2) :: elnamec
      character(len=60) :: bundlingc
!cym/cpg end
      integer :: ndum1(1), ndum2(1), ndum4(1), mxdm
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=420), dimension(:), allocatable ::
     .           crs_lines
!cym/cpg end
      integer, allocatable :: lknd(:), jdnsl(:)
      external :: eirene_leer, eirene_masage, eirene_exit_own

      LHYDDEF = .FALSE.
      IMSG=0

      CALL EIRENE_MASAGE
     .  ('*** 4. DATA FOR SPECIES SPECIFICATION AND')
      CALL EIRENE_MASAGE
     .  ('       ATOMIC PHYSICS MODULE')
      CALL EIRENE_LEER(1)

      call json%get_child(p,'SPECIES_SPEC',pspc,found)
      call json%get_child(pspc,'REACTIONS',preac,found)
      call json%get(preac,'NREACI',nreaci,found)

      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'
      WRITE (iunout,*) '       NREACI= ',NREACI
      CALL EIRENE_LEER(1)

      call json%get(preac,'ClassicReactionStrings',crs_lines,found_crs)
      if (found_crs) then
        ncrs = size(crs_lines)
        do k = 1, ncrs
!  STORE CRS LINES FOR OUTPUTING TO JSON FILE
          call eirene_push_string_stack(crs_stack,crs_lines(k))
        end do
        deallocate (crs_lines)
      end if

      IF (NPHOTI > 0) CALL EIRENE_PH_INIT(0)

      call json%get(preac,'REAC_SPECS',preacs,foundr)
      if (foundr) then
        call json%info(preacs,n_children=nrea)
        if (nrea < nreaci) then
          write (iunout,*) ' NUMBER OF REACTIONS SPECIFIED DOES',
     .          ' NOT MATCH NUMBER OF REACTIONS FOUND IN FILE '
          write (iunout,*) 'NREACI = ',nreaci
          write (iunout,*) 'NREA =   ',nrea
          call eirene_exit_own(1)
        end if
        do il = 1, nrea ! use nrea as one reaction might have
                        ! several reaction lines
          filnam = repeat(' ',8)
          h123 = repeat(' ',4)
          reac = repeat(' ',50)
          bundling = repeat(' ',60)
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
          iftfl = 0
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

          call json%get_child(preacs,il,prea,found)
          call json%get(prea,'IR',ir,found)
          call json%get(prea,'FILNAM',chr,found)
          filnam = chr(1:8)
          deallocate(chr)
          call json%get(prea,'H123',chr,found)
          h123 = chr(1:4)
          deallocate(chr)
          call json%get(prea,'REAC',chr,found)
          reac = chr(1:50)
          deallocate(chr)
          call json%get(prea,'CRC',chr,found)
          crc = chr(1:3)
          deallocate(chr)
          call json%get(prea,'MASSP',mp,found)
          call json%get(prea,'MASST',mt,found)
          call json%get(prea,'DP',dpp,found)
          call json%get(prea,'R1MN',r1mn,found)
          call json%get(prea,'R1MX',r1mx,found)
          call json%get(prea,'R2MN',r2mn,found)
          call json%get(prea,'R2MX',r2mx,found)

          call json%get(prea,'JFEX1MN',jfex1mn,found)
          call json%get(prea,'JFEX1MX',jfex1mx,found)
          call json%get(prea,'JFEX2MN',jfex2mn,found)
          call json%get(prea,'JFEX2MN',jfex2mn,found)

          call json%get(prea,'FP1',ffp1,found)
          call json%get(prea,'FP2',ffp2,found)

          IF (INDEX(H123,json_ck_'P.').eq.0) then
c  single parametric data, or for first parameter in 2-parameteric data
            RC1MIN = -20._DP  ! lower ln(E), ln(T) default limit;
                              ! E,T in eV  (2E-9 EV)
            RC1MAX =  20._DP  ! upper ln(E), ln(T) default limit;
                              ! E,T in eV  (5E8  EV)
c  2nd parameter in 2 parametric data
            RC2MIN = -20._DP  ! lower ln(E0), ln(N) default limit;
                              ! E0 in eV, N in cm**-3
            RC2MAX =  100._DP ! upper ln(E0), ln(N) default limit;
                              ! E0 in eV, N in cm**-3
          END IF

          IF (R1MN > 0) THEN
            fp1(1:3) = ffp1(1:3)
            RC1MIN=LOG(R1MN)
            WRITE (IUNOUT,*) 'NON-DEF. R1MN FOR REACTION IR=',IR,R1MN
          ENDIF
          IF (R1MX > 0) THEN
            fp1(4:6) = ffp1(4:6)
            RC1MAX=LOG(R1MX)
            WRITE (IUNOUT,*) 'NON-DEF. R1MX FOR REACTION IR=',IR,R1MX
          ENDIF
          IF (R2MN > 0) THEN
            fp2(1:3) = ffp2(1:3)
            RC2MIN=LOG(R2MN)
            WRITE (IUNOUT,*) 'NON-DEF. R2MN FOR REACTION IR=',IR,R2MN
          ENDIF
          IF (R2MX > 0) THEN
            fp2(4:6) = ffp2(4:6)
            RC2MAX=LOG(R2MX)
            WRITE (IUNOUT,*) 'NON-DEF. R2MX FOR REACTION IR=',IR,R2MX
          ENDIF

          if (allocated(ffp1)) deallocate(ffp1)
          if (allocated(ffp2)) deallocate(ffp2)

          if (filnam == json_ck_'CONST') then
            call json%get(prea,'IFTFLG',iftfl,found)
            call json%get(prea,'NCOEF',ncoef,found)
            call json%get(prea,'COEF',ccoef,found)
            coef(1:ncoef) = ccoef(1:ncoef)
            deallocate(ccoef)
          else
            iftfl = 0
            ncoef = 0
          end if

          if (index(filnam,json_ck_'ADAS') > 0 .or.
     .        index(filnam,json_ck_'TAB2D') > 0) then
            call json%get(prea,'ELNAME',chr,found)
            elname = chr(1:2)
            deallocate(chr)
            call json%get(prea,'IZ',iz,found)
            call json%get(prea,'BUNDLING',chr,found)
            if (found) then
              ll = len_trim(trim(chr))
              bundling = chr(1:min(60,ll))
              if (ll == 0) bundling = ' '
              deallocate(chr)
            else
              bundling = repeat(' ',60)
            end if
          else
            elname = '  '
            iz = 0
            bundling = repeat(' ',60)
          end if

          if (index(filnam,json_ck_'CR') /= 0) then
            call json%get(prea,'M_POPESC',m_popesc,found)
            if (m_popesc > 0) then
              if (ALLOCATED(irow_esc)) then
                  DEALLOCATE (irow_esc)
                  DEALLOCATE (icol_esc)
                  DEALLOCATE (pop_esc)
              end if
              ALLOCATE(irow_esc(m_popesc), SOURCE = 0)
              ALLOCATE(icol_esc(m_popesc), SOURCE = 0)
              ALLOCATE(pop_esc(m_popesc), SOURCE = 1._DP)
              call json%get(prea,'IROW_ESC',irow_esc,found)
              call json%get(prea,'ICOL_ESC',icol_esc,found)
              call json%get(prea,'POP_ESC',pop_esc,found)
            end if
            call json%get(prea,'M_QEXT',m_qext,found)
            if (m_qext > 0) then
              if (ALLOCATED(irow_ext)) then
                  DEALLOCATE (irow_ext)
                  DEALLOCATE (ital_ext)
                  DEALLOCATE (ital_ext_1)
                  DEALLOCATE (ital_ext_2)
                  DEALLOCATE (ital_ext_3)
              end if
              ALLOCATE(irow_ext(m_qext), SOURCE = 0)
              ALLOCATE(ital_ext(3,m_qext), SOURCE = 0)
crc separate ital_ext due to json%get vector argument limitation
              ALLOCATE(ital_ext_1(m_qext), SOURCE = 0)
              ALLOCATE(ital_ext_2(m_qext), SOURCE = 0)
              ALLOCATE(ital_ext_3(m_qext), SOURCE = 0)
              call json%get(prea,'IROW_EXT',irow_ext,found)
              call json%get(prea,'ITAL_EXT_1',ital_ext_1,found)
              call json%get(prea,'ITAL_EXT_2',ital_ext_2,found)
              call json%get(prea,'ITAL_EXT_3',ital_ext_3,found)
crc recombine ital_ext for slreac call later
              ital_ext(1,1:m_qext) = ital_ext_1(1:m_qext)
              ital_ext(2,1:m_qext) = ital_ext_2(1:m_qext)
              ital_ext(3,1:m_qext) = ital_ext_3(1:m_qext)              
            end if
crc not used
c            call json%get(prea,'M_UPPER',m_upper,found)
c            call json%get(prea,'M_LOWER',m_lower,found)
          end if
C
C  SAVE INPUT LINES ON REACLINES,  FOR WRITING ONTO JSON-FILE
C       ONLY NEEDED FOR AUTOMATED INTERFACE TO HYDKIN DATABASE.

!pb       IL = IRLINES+1
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
           REACLINES(IL)%ITAL_EXT(1:3,1:M_QEXT) = 
     .      ITAL_EXT(1:3,1:M_QEXT)
          END IF
          REACLINES(IL)%IFTFLG = IFTFL
          REACLINES(IL)%NCOEF = NCOEF
          IF (NCOEF > 0) REACLINES(IL)%COEF(1:NCOEF) = COEF(1:NCOEF)
          IRLINES = IL

C
C  SAVE SOME OF THE INPUT FLAGS FOR LATER
C  PROCESSING (MASS SCALING, POTENTIAL ENERGY INCREMENT) IN XSTCX,XSTEI,...

          MASSP(IR)=MP
          MASST(IR)=MT
          DELPOT(IR)=DPP

!cym/cpg avoid type mismatch in call to SLREAC
!         CALL EIRENE_SLREAC (IR,FILNAM,H123,REAC,CRC,
!    .                  RC1MIN,RC1MAX,FP1,JFEX1MN,JFEX1MX, ! additional input card: asymptotics P1
!    .                  RC2MIN,RC2MAX,FP2,JFEX2MN,JFEX2MX, ! additional input card: asymptotics P2
!    .                  ELNAME,IZ,BUNDLING,                ! additional input card read for TAB2D/ADAS format
!    .                  IROW_ESC,ICOL_ESC,POP_ESC, ! (optional) additional input card read CR  internal models
!    .                  IFTFL, NCOEF, COEF)        ! (optional) for "CONST models"
          FILNAMc=FILNAM
          H123c=H123
          REACc=REAC
          CRCc=CRC
          ELNAMEc=ELNAME
          BUNDLINGC=BUNDLING
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
!cym/cpg end

          nullify(prea)
        end do  ! nreaci

        nullify (preacs)
        nullify (preac)
      end if
      IF (TRCAMD) CALL EIRENE_LEER(1)

!   set up pointers to particle types (ATOMS, MOLEUCLES, TEST_IONS, PHOTONS)
!   read the respective species numbers
      call json%get_child(pspc,'ATOMS',patm,found)
      call json%get(patm,'NATMI',natmi,found)

      call json%get_child(pspc,'MOLECULES',pmol,found)
      call json%get(pmol,'NMOLI',nmoli,found)

      call json%get_child(pspc,'TEST_IONS',pion,found)
      call json%get(pion,'NIONI',nioni,found)

      call json%get_child(pspc,'PHOTONS',pphot,found)
      call json%get(pphot,'NPHOTI',nphoti,found)

      NSPH=NPHOTI
      NSPA=NSPH+NATMI
      NSPAM=NSPH+NATMI+NMOLI
      NSPAMI=NSPH+NATMI+NMOLI+NIONI

      mxdm=max(nphot,natm,nmol,nion,1)
      allocate (lknd(mxdm))
      allocate (jdnsl(mxdm))
C
C  READ NEUTRAL ATOMS SPECIES CARDS
C

      WRITE (iunout,*)
     . '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'
      WRITE (iunout,*) '       NATMI= ',NATMI
      CALL EIRENE_LEER(1)

      lmulpl = .false.

      ityp = 1
      call eirene_read_block_4abcd
     .    (json, patm, ityp, nsph, natmi, natm, nspz,
     .     texts, 'IATM', 'A', nmassa, nchara, nprt, ndum2, isrf, isrt,
     .     nrca, nfola, ngena, nhsts, ireaca, ibulka,
     .     iscd1a, iscd2a, iscd3a, iscd4a, iscdea, iestma, ibgka,
     .     eeleca, ebulka, escd1a, freaca, edpota, lmulpl,
     .     lknd, jdnsl)
      nullify(patm)

C
C  READ NEUTRAL MOLECULES SPECIES CARDS
C
      WRITE (iunout,*)
     . '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'
      WRITE (iunout,*) '       NMOLI= ',NMOLI
      CALL EIRENE_LEER(1)

      ityp = 2
      call eirene_read_block_4abcd
     .    (json, pmol, ityp, nspa, nmoli, nmol, nspz,
     .     texts, 'IMOL', 'M', nmassm, ncharm, nprt, ndum2, isrf, isrt,
     .     nrcm, nfolm, ngenm, nhsts, ireacm, ibulkm,
     .     iscd1m, iscd2m, iscd3m, iscd4m, iscdem, iestmm, ibgkm,
     .     eelecm, ebulkm, escd1m, freacm, edpotm, lmulpl,
     .     lknd, jdnsl)
      lkindm(1:nmoli) = lknd(1:nmoli)
      nullify(pmol)

C
C  READ TEST PARTICLE IONS SPECIES CARDS
C
      WRITE (iunout,*) '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'
      WRITE (iunout,*) '       NIONI= ',NIONI
      CALL EIRENE_LEER(1)

      ityp = 3
      call eirene_read_block_4abcd
     .    (json, pion, ityp, nspam, nioni, nion, nspz,
     .     texts, 'IION', 'I', nmassi, nchari, nprt, nchrgi, isrf, isrt,
     .     nrci, nfoli, ngeni, nhsts, ireaci, ibulki,
     .     iscd1i, iscd2i, iscd3i, iscd4i, iscdei, iestmi, ibgki,
     .     eeleci, ebulki, escd1i, freaci, edpoti, lmulpl,
     .     lknd, jdnsl)
      lkindi(1:nioni) = lknd(1:nioni)
      nullify(pion)

C
C  READ NEUTRAL PHOTONS SPECIES CARDS
C
      IF (NPHOTI.GT.0) THEN
        WRITE (iunout,*)
     .   '*** 4D. NEUTRAL PHOTONS SPECIES CARDS, NPHOTI SPECIES'
        WRITE (iunout,*) '       NPHOTI= ',NPHOTI
        CALL EIRENE_LEER(1)
      END IF

      ityp = 0
      call eirene_read_block_4abcd
     .    (json, pphot, ityp, 0, nphoti, nphot, nspz,
     .     texts, 'IPHOT', 'PH', ndum1, ndum2, nprt, ndum4, isrf, isrt,
     .     nrcph, nfolph, ngenph, nhsts, ireacph, ibulkph,
     .     iscd1ph, iscd2ph, iscd3ph, iscd4ph, iscdeph, iestmph, ibgkph,
     .     eelecph, ebulkph, escd1ph, freacph, edpotph, lmulpl,
     .     lknd, jdnsl)
      nullify(pphot)

      nullify(pspc)

      deallocate (lknd)
      deallocate (jdnsl)

      return
      end subroutine eirene_read_block_4

!*****************************************************************

      subroutine eirene_read_block_4abcd
     .    (json, me, ityp, nbas, nloop, ndim, nspz,
     .     texts, cndx, cext, nmass, nchar, nprt, nchrg, isrf, isrt,
     .     nrc, nfol, ngen, nhsts, ireac, ibulk,
     .     iscd1, iscd2, iscd3, iscd4, iscde, iestm, ibgk,
     .     eelec, ebulk, escd1, freac, edpot, lmulpl,
     .     lknd, jdnsl, cdenmodel)

      class(json_core),intent(inout) :: json
      type(json_value),pointer :: me
      type(json_value),pointer :: blk, prt, rea, reas, cden
      integer, intent(in) :: ityp, nbas, nloop, ndim, nspz
      integer, intent(inout) :: nmass(ndim), nchar(ndim),
     .         isrf(nspz,0:nlimps), isrt(nspz,0:nlimps),
     .         nprt(nspz), nchrg(ndim),
     .         nrc(ndim), nfol(ndim), ngen(ndim), nhsts(nspz),
     .         ireac(ndim,*), ibulk(ndim,*), iscd1(ndim,*),
     .         iscd2(ndim,*), iscd3(ndim,*), iscd4(ndim,*),
     .         iscde(ndim,*), iestm(ndim,*), ibgk(ndim,*)
      integer, intent(out) :: lknd(ndim), jdnsl(ndim)
      real(dp), intent(inout) :: eelec(ndim,*), ebulk(ndim,*),
     .         escd1(ndim,*), freac(ndim,*), edpot(ndim,*)
      logical, intent(inout) :: lmulpl
      character(*) :: texts(*), cndx, cext
      character(*), optional :: cdenmodel(ndim)
      integer :: i, numsec, k, id, nr, np, nsc, ltxt, nre
      logical :: lden, foundb, foundr, found, foundc
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      external :: eirene_leer, eirene_exit_own

      lden = present(cdenmodel)

      lknd = 0
      jdnsl = 0

      call json%get_child(me,'SPECIES',blk,foundb)
      if (foundb) then
        call json%info(blk,n_children=np)
        if (np /= nloop) then
          write (iunout,*) ' NUMBER OF SPECIES DOES',
     .          ' NOT MATCH NUMBER OF SPECIES FOUND IN FILE '
          write (iunout,*) 'NSPECIES = ',nloop
          write (iunout,*) 'NP_FILE =  ',np
          call eirene_exit_own(1)
        end if

        do i = 1, nloop
          call json%get_child(blk,i,prt,found)
          call json%get(prt,cndx,id,found)
          ispz = nbas + id
          call json%get(prt,'SPECIES',txt,found)
          texts(ispz) = txt
          deallocate(txt)
          call json%get(prt,'NMASS'//cext,nmass(i),found)
          call json%get(prt,'NCHAR'//cext,nchar(i),found)

! DEFAULT FOR PHOTONIC AND ATOMIC SPECIES
          if (ityp <= 1) nprt(ispz) = 1

          if (ityp >= 2)
     .      call json%get(prt,'NPRT',nprt(ispz),found)

          if (ityp >= 3)
     .      call json%get(prt,'NCHRG'//cext,nchrg(i),found)

          call json%get(prt,'ISRF',isrf(ispz,1),found)
          call json%get(prt,'ISRT',isrt(ispz,1),found)
          call json%get(prt,'NUMSEC',numsec,found)
          call json%get(prt,'NRC'//cext,nrc(i),found)

          if (ityp < 4) then
            call json%get(prt,'NFOL'//cext,nfol(i),found)
            call json%get(prt,'NGEN'//cext,ngen(i),found)
          end if

          call json%get(prt,'NHSTS',nhsts(ispz),found)

          if ((ityp == 2) .or. (ityp == 3)) then
            call json%get(prt,'LKIND'//cext,lknd(i),found)
          end if

          if (ityp == 4) then
            call json%get(prt,'DENSLIM',jdnsl(i),found)
          end if

          if (ityp == 2) then
            IF (ISRT(ISPZ,1).LT.0) THEN
              WRITE (iunout,*) 'INPUT ERROR IN BLOCK 4B '
              WRITE (iunout,*) 'MOLECULAR SPECIES ',I,':'
              WRITE (iunout,*)
     .          'ISRT LT 0 OPTION IS NOT AVAILABLE ANYMORE'
              WRITE (iunout,*) 'PROBABLY YOU MEAN: ISRT= ',NLOOP+1
              call eirene_exit_own(1)
            ENDIF
          end if

          if (nrc(i) > 0) then
            call json%get_child(prt,'REACTIONS',reas,foundr)
            if (foundr) then
              call json%info(reas,n_children=nr)
              if (nr /= nrc(i)) then
                write (iunout,*) ' NUMBER OF REACTIONS FOR',
     .            cndx, i, ' DOES',
     .            ' NOT MATCH NUMBER OF REACTIONS FOUND IN FILE '
                write (iunout,*) 'NRC     = ',nrc(i)
                write (iunout,*) 'NR_FILE =  ',nr
                call eirene_exit_own(1)
              end if

              do k = 1, nrc(i)
                call json%get_child(reas,k,rea,found)
                call json%get(rea,'IREAC'//cext,ireac(i,k),found)
                call json%get(rea,'IBULK'//cext,ibulk(i,k),found)
                call json%get(rea,'ISCD1'//cext,iscd1(i,k),found)
                call json%get(rea,'ISCD2'//cext,iscd2(i,k),found)
                if (numsec > 2)
     .            call json%get(rea,'ISCD3'//cext,iscd3(i,k),found)
                if (numsec > 3)
     .            call json%get(rea,'ISCD4'//cext,iscd4(i,k),found)
                call json%get(rea,'ISCDE'//cext,iscde(i,k),found)
                if (ityp < 4) then
                  call json%get(rea,'IESTM'//cext,iestm(i,k),found)
                  call json%get(rea,'IBGK'//cext,ibgk(i,k),found)
                end if

                if ((numsec == 3) .and. (iscd3(i,k) > 0)) then
                  WRITE (iunout,*) ' WARNING: REACTION KK= ', ireac(i,k)
                  WRITE (iunout,*)
     .              ' THREE SECONDARY GROUPS USED FOR SPECIES ',
     .              TEXTS(ISPZ)
                  WRITE (iunout,*) ' ISCD1'//cext,' = ',ISCD1(I,K)
                  WRITE (iunout,*) ' ISCD2'//cext,' = ',ISCD2(I,K)
                  WRITE (iunout,*) ' ISCD3'//cext,' = ',ISCD3(I,K)
                  WRITE (iunout,*) ' ISCDE'//cext,' = ',ISCDE(I,K)
                  WRITE (iunout,*) ' IESTM'//cext,' = ',IESTM(I,K)
                  WRITE (iunout,*) ' IBGK'//cext,'  = ',IBGK(I,K)
                elseif ((numsec == 4) .and. (iscd4(i,k) > 0)) then
                  WRITE (iunout,*) ' WARNING: REACTION KK= ', ireac(i,k)
                  WRITE (iunout,*)
     .              ' FOUR SECONDARY GROUPS USED FOR SPECIES ',
     .              TEXTS(ISPZ)
                  WRITE (iunout,*) ' ISCD1'//cext,' = ',ISCD1(I,K)
                  WRITE (iunout,*) ' ISCD2'//cext,' = ',ISCD2(I,K)
                  WRITE (iunout,*) ' ISCD3'//cext,' = ',ISCD3(I,K)
                  WRITE (iunout,*) ' ISCD4'//cext,' = ',ISCD4(I,K)
                  WRITE (iunout,*) ' ISCDE'//cext,' = ',ISCDE(I,K)
                  WRITE (iunout,*) ' IESTM'//cext,' = ',IESTM(I,K)
                  WRITE (iunout,*) ' IBGK'//cext,'  = ',IBGK(I,K)
                end if

                if (ityp < 4) then
CDR  FOR THE TIME BEING:  NONLINEAR BGK REACTIONS MUST BE ELASTIC
                  kk=IREAC(I,K)
                  elonly=iswr(kk)
                  IF (IBGK(I,K).NE.0 .AND. elonly.ne.5) then
                    WRITE (IUNOUT,*)
     .                 'REACT. KK: BGK ITERAT. FOR IMOL IGNORED. NOT EL'
                    WRITE (IUNOUT,*) 'I, KK', I,KK, elonly
                    IBGK(I,K)=0
                    IMSG=1
                  endif
cdr  enforce multiple background field temperature arrays, in case of any BGK reactions.
                  LMULPL = LMULPL .OR. (IBGK(I,K) /= 0)
                end if

                call json%get(rea,'EELEC'//cext,eelec(i,k),found)
                call json%get(rea,'EBULK'//cext,ebulk(i,k),found)
                call json%get(rea,'ESCD1'//cext,escd1(i,k),found)
                call json%get(rea,'FREAC'//cext,freac(i,k),found)
                call json%get(rea,'EDPOT'//cext,edpot(i,k),found)

                NSC = 0
                IF (ISCD3(I,K) > 0) NSC = 3
                IF (ISCD4(I,K) > 0) NSC = 4
                if (nsc > 0) then
                  IF ((REACDAT(IREAC(I,K))%NOSEC > 0) .AND.
     .                (REACDAT(IREAC(I,K))%NOSEC /= NSC)) THEN
                    WRITE (IUNOUT,*) ' INCONSISTENCY FOUND CONCERNING',
     .                ' REACTION ',IREAC(I,K)
                    WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .                ' PREVIOUSLY WAS ',REACDAT(IREAC(I,K))%NOSEC
                    WRITE (IUNOUT,*) ' NUMBER OF SECONDARIES FOUND',
     .                ' NOW IS         ',NSC
                    WRITE (IUNOUT,*) ' USE',
     .                MAX(REACDAT(IREAC(I,K))%NOSEC,NSC),
     .                ' SECONDARIES '
                    REACDAT(IREAC(I,K))%NOSEC =
     .                MAX(REACDAT(IREAC(I,K))%NOSEC,NSC)
                  ELSE
                    REACDAT(IREAC(I,K))%NOSEC = NSC
                  END IF
                end if
                nullify(rea)
              end do  ! nrc
            end if  ! foundr
          end if  ! nrc > 0
          IF (IMSG.EQ.1) CALL EIRENE_LEER(1)
          IMSG=0

!  check for density models
          if ((ityp == 4) .and. lden) then
            call json%get_child(prt,'DENS_MODEL',cden,foundc)
            if (foundc) then
              call json%get(cden,'MODEL',txt,found)
              ltxt = min(len(txt),10)
              cdenmodel(i)(1:ltxt) = txt(1:ltxt)
              deallocate(txt)

              ALLOCATE (TDMPAR(I)%TDM)
              call json%get(cden,'NRE',nre,found)
              TDMPAR(I)%TDM%NRE=MAX(NRE,1)
              ALLOCATE (TDMPAR(I)%TDM%ISP(TDMPAR(I)%TDM%NRE))
              ALLOCATE (TDMPAR(I)%TDM%ITP(TDMPAR(I)%TDM%NRE))
              ALLOCATE (TDMPAR(I)%TDM%ISTR(TDMPAR(I)%TDM%NRE))
              ALLOCATE (TDMPAR(I)%TDM%IRC(TDMPAR(I)%TDM%NRE))
              SELECT CASE (CDENMODEL(I))
              CASE ('FORT.13','FTN13')
                IDMDL = IDMDL + 1
                call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(1),found)
c  default: only for bulk ions
                TDMPAR(I)%TDM%ITP(1)=4
              CASE ('FORT.10','FTN10')
                IDMDL = IDMDL + 1
                call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(1),found)
                call json%get(cden,'ITP',TDMPAR(I)%TDM%ITP(1),found)
                call json%get(cden,'ISTR',TDMPAR(I)%TDM%ISTR(1),found)
              CASE ('CONSTANT  ')
                IDMDL = IDMDL + 1
                call json%get(cden,'TVAL',TDMPAR(I)%TDM%TVAL,found)
                call json%get(cden,'DVAL',TDMPAR(I)%TDM%DVAL,found)
                call json%get(cden,'VXVAL',TDMPAR(I)%TDM%VXVAL,found)
                call json%get(cden,'VYVAL',TDMPAR(I)%TDM%VYVAL,found)
                call json%get(cden,'VZVAL',TDMPAR(I)%TDM%VZVAL,found)
              CASE ('MULTIPLY  ')
                IDMDL = IDMDL + 1
                call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(1),found)
                call json%get(cden,'ITP',TDMPAR(I)%TDM%ITP(1),found)
                call json%get(cden,'ISTR',TDMPAR(I)%TDM%ISTR(1),found)
                call json%get(cden,'DFACTOR',
     .                        TDMPAR(I)%TDM%DFACTOR,found)
                call json%get(cden,'TFACTOR',
     .                        TDMPAR(I)%TDM%TFACTOR,found)
                call json%get(cden,'VFACTOR',
     .                        TDMPAR(I)%TDM%VFACTOR,found)
                TDMPAR(I)%TDM%ITP(1)=4
              CASE ('SAHA      ')
                IDMDL = IDMDL + 1
!PB   TO BE WRITTEN
              CASE ('PLANCK')
                IDMDL = IDMDL + 1
                call json%get(cden,'G_PLANCK',
     .                        TDMPAR(I)%TDM%G_PLANCK,found)
                TDMPAR(I)%TDM%ITP(1)=4
!dr   unfinished, TO BE WRITTEN
              CASE ('BOLTZMANN ')
                IDMDL = IDMDL + 1
                call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(1),found)
                call json%get(cden,'ITP',TDMPAR(I)%TDM%ITP(1),found)
                call json%get(cden,'ISTR',TDMPAR(I)%TDM%ISTR(1),found)
                call json%get(cden,'G_BOLTZ',
     .                        TDMPAR(I)%TDM%G_BOLTZ,found)
                call json%get(cden,'DELTAE',TDMPAR(I)%TDM%DELTAE,found)
              CASE ('CORONA    ')
                IDMDL = IDMDL + 1
                call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(1),found)
                call json%get(cden,'ITP',TDMPAR(I)%TDM%ITP(1),found)
                call json%get(cden,'ISTR',TDMPAR(I)%TDM%ISTR(1),found)
                call json%get(cden,'IRC',TDMPAR(I)%TDM%IRC(1),found)
                call json%get(cden,'A_CORONA',
     .                        TDMPAR(I)%TDM%A_CORONA,found)
                IF (.NOT.REACDAT(IR)%LRTC) THEN
                  WRITE (iunout,*)
     .              ' WRONG REACTION SPECIFIED FOR CORONA MODEL '
                  WRITE (iunout,*) ' ONLY H.2 REACTIONS ARE PERMITTED '
                  WRITE (iunout,*) ' IPLS = ',I
                  CALL EIRENE_EXIT_OWN(1)
                END IF
              CASE ('COLRAD    ')
                IDMDL = IDMDL + 1
                DO J=1, TDMPAR(I)%TDM%NRE
                  call json%get(cden,'ISP',TDMPAR(I)%TDM%ISP(J),found)
                  call json%get(cden,'ITP',TDMPAR(I)%TDM%ITP(J),found)
                  call json%get(cden,'ISTR',TDMPAR(I)%TDM%ISTR(J),found)
                  call json%get(cden,'IRC',TDMPAR(I)%TDM%IRC(J),found)
                  IR = TDMPAR(I)%TDM%IRC(J)
                  IF (.NOT.REACDAT(IR)%LOTH) THEN
                    WRITE (iunout,*)
     .                ' WRONG REACTION SPECIFIED FOR COLRAD MODEL'
                    WRITE (iunout,*)
     .                ' ONLY H.11 OR H.12 (OT) REACTION DATA'
                    WRITE (IUNOUT,*) ' ARE PERMITTED'
                    WRITE (iunout,*) ' IPLS = ',I
                    CALL EIRENE_EXIT_OWN(1)
                  END IF
                END DO
              CASE DEFAULT
!PB  NOTHING TO BE DONE
              END SELECT
              nullify(cden)
            end if  ! foundc
          end if  ! ityp == 4 .and. lden

        end do  ! nloop

      end if  ! foundb

      nullify(blk)
      nullify(prt)
      nullify(reas)

      end subroutine eirene_read_block_4abcd

!*****************************************************************

      subroutine eirene_read_block_5(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pbck, pbulk, plsm, pch,
     .                             pte, pti, pdi, pvel, pbfld, pvol,
     .                             pintl, pAMD, pzi
      integer, allocatable :: ihelp(:)
      integer :: i, ndum1(1), ndum2(1), ndum3(1), ndum4(1), nch,
     .           ital, iopt, jpls
      logical :: found, foundp, foundi
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: MeshTypeCK
      character(kind=json_CK,len=:), allocatable :: AMDoutputCK
!cym/cpg end
      integer, allocatable :: lknd(:), jdnsl(:)
      external :: eirene_uppercase, eirene_leer, eirene_masage,
     .            eirene_exit_own

      CALL EIRENE_MASAGE('*** 5. DATA FOR PLASMA BACKGROUND')
      CALL EIRENE_LEER(1)

      call json%get_child(p,'BACKGROUND',pbck,found)

C
C  READ BULK IONS SPECIES CARDS
C
      call json%get_child(pbck,'BULK_IONS',pbulk,found)
      call json%get(pbulk,'NPLSI',nplsi,found)

      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'
      WRITE (iunout,*) '       NPLSI= ',NPLSI
      CALL EIRENE_LEER(1)
      NSPAMI=NSPAM+NIONI
      NSPTOT=NSPAMI+NPLSI
! counter for additional reaction cards possibly needed for "density models".
      IDMDL = 0

      allocate (lknd(npls))
      allocate (jdnsl(npls))

      ityp = 4
      call eirene_read_block_4abcd
     .    (json, pbulk, ityp, nspami, nplsi, npls, nspz,
     .     texts, 'IPLS', 'P', nmassp, ncharp, nprt, nchrgp, isrf, isrt,
     .     nrcp, ndum1, ndum2, nhsts, ireacp, ibulkp,
     .     iscd1p, iscd2p, iscd3p, iscd4p, iscdep, ndum3, ndum4,
     .     eelecp, ebulkp, escd1p, freacp, edpotp, lmulpl,
     .     lknd, jdnsl, cdenmodel)

      nullify(pbulk)
C
      DO I=1,NSPZ
        CALL EIRENE_UPPERCASE (TEXTS(I))
      END DO

!PB DENSITY LIMIT FOR PLASMA SPECIES
      DO JPLS=1,NPLSI
        ISPZ=NSPAMI+JPLS
        IF(JDNSL(JPLS).EQ.0) THEN
c         DENSLIM(JPLS)=10**(MAXEXPONENT(DENSLIM(JPLS))-1)
        DENSLIM(JPLS)=10e30
        ELSE
         DENSLIM(JPLS)=2._DP**JDNSL(JPLS)
        END IF
        WRITE(iunout,*) TEXTS(ISPZ), DENSLIM(JPLS)
      END DO
      CALL EIRENE_LEER(1)

      deallocate (lknd)
      deallocate (jdnsl)
C
      call json%get_child(pbck,'PLASMA',plsm,foundp)

      CALL EIRENE_MASAGE('*** 5B. PLASMA BACKGROUND DATA')
      CALL EIRENE_LEER(1)

cdr  indpro(2) has already been determined in find_param.
      call json%get(plsm,'INDPRO',ihelp,found)
      indpro(1:12) = ihelp(1:12)
      indpro_in(1:12) = indpro(1:12)
      deallocate(ihelp)
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

      call json%get_child(plsm,'TE',pte,found)
      call json%get_child(plsm,'TI',pti,found)
      call json%get_child(plsm,'DI',pdi,found)
      call json%get_child(plsm,'VEL',pvel,found)
      call json%get_child(plsm,'B_FIELD',pbfld,found)
      call json%get_child(plsm,'VOL',pvol,found)
      call json%get_child(plsm,'INPUT_TALLY_OPTIONS',pintl,foundi)

C  Te profile
      if ((indpro(1) <= 5) .and. (nplsi > 0)) then
        call json%get(pte,'TE0',te0,found)
        call json%get(pte,'TE1',te1,found)
        call json%get(pte,'TE2',te2,found)
        call json%get(pte,'TE3',te3,found)
        call json%get(pte,'TE4',te4,found)
        call json%get(pte,'TE5',te5,found)
      end if
      nullify (pte)

C  Ti profile(s)
cdr  default is: eV units for Ti
      NLKELV=INDPRO(2) < 0  ! Kelvin units instead of eV
      INDPRO(2)=IABS(INDPRO(2))

cdr  default is: |indpro| < 10: nplsi Ti fields, one per background species
cdr              |indpro| > 10: only one common Ti field is specified

      NPLSTI = NPLSI
      LRDMLTI=.TRUE.

cdr  only one common profile for all NPLS species?
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
        LRDMLTI=.FALSE.         !  read only one common Ti card,
                                !  despite storage for NPLS Ti profiles.
        NPLSTI = NPLSI
        WRITE (IUNOUT,*) '      NPLSTI = ',NPLSTI
        CALL EIRENE_LEER(1)
c     ELSE
c  read one Ti card for each of the NPLS Ti profiles
c       LRDMLTI=.TRUE.
      END IF

      IF (INDPRO(2) > 9) INDPRO(2) = MOD(INDPRO(2),10)

      NLMLTI = (NPLSTI > 1)

cdr  indirect addressing: IVL= MPLSTI(IPLS)   ! find Ti field for IPLS on TIIN(IVL)
      IF (NPLSTI == 1) THEN
        MPLSTI = 1                    ! find ALL Ti fields on
                                      ! IPLS=1 storage
      ELSEIF (NPLSTI == NPLS) THEN
        MPLSTI = (/ (I,I=1,NPLS) /)   ! find each individual Ti field
                                      ! on its IPLS storage
      ELSE
cdr     more general options: to be written.
cdr     E.g.: only virtual background species with different flow fields
      ENDIF

      MPLSTI=1
      IF (NLMLTI)     MPLSTI = (/ (I,I=1,NPLS) /)


      IF (INDPRO(2).LE.5.AND.NPLSI.GT.0) THEN
        call json%info(pti,n_children=nch)
        if (lrdmlti.and.(nch /= nplsi)) then
          write (iunout,*) ' ERROR READING ION TEMPERATURE PROFILES '
          write (iunout,*) ' TOO FEW TEMPERATURE PROFILES PROVIDED '
          write (iunout,*) ' NO. PROFILES FOUND ',nch
          write (iunout,*) ' NPLSI = ',NPLSI
          call eirene_exit_own(1)
        end if
        if (.not.lrdmlti .and. (nch > 1)) nch = 1
        do i = 1,nch
          call json%get_child(pti,i,pch,found)
          call json%get(pch,'TI0',ti0(i),found)
          call json%get(pch,'TI1',ti1(i),found)
          call json%get(pch,'TI2',ti2(i),found)
          call json%get(pch,'TI3',ti3(i),found)
          call json%get(pch,'TI4',ti4(i),found)
          call json%get(pch,'TI5',ti5(i),found)
          if (.not.lrdmlti) then
            DO JPLS=2,NPLSTI
              TI0(JPLS)=TI0(1)
              TI1(JPLS)=TI1(1)
              TI2(JPLS)=TI2(1)
              TI3(JPLS)=TI3(1)
              TI4(JPLS)=TI4(1)
              TI5(JPLS)=TI5(1)
            END DO
          end if
          nullify (pch)
        end do
      END IF
      nullify (pti)

c  di profiles
      IF (INDPRO(3).LE.5) THEN
        call json%info(pdi,n_children=nch)
        if (nch /= nplsi) then
          write (iunout,*) ' ERROR READING BULK ION DENSITY PROFILES '
          write (iunout,*) ' NO. OF PROFILES DOES NOT MATCH',
     .                     ' NO. OF BULK SPECIES'
          write (iunout,*) ' NO. PROFILES FOUND ',nch
          write (iunout,*) ' NPLSI = ',NPLSI
          call eirene_exit_own(1)
        end if
        do i = 1,nch
          call json%get_child(pdi,i,pch,found)
          call json%get(pch,'DI0',di0(i),found)
          call json%get(pch,'DI1',di1(i),found)
          call json%get(pch,'DI2',di2(i),found)
          call json%get(pch,'DI3',di3(i),found)
          call json%get(pch,'DI4',di4(i),found)
          call json%get(pch,'DI5',di5(i),found)
          nullify (pch)
        end do
      END IF
      nullify (pdi)

c  V_IN profile(s)
cdr  default is: cm/s units for flow field(s)
      NLMACH=INDPRO(4).LT.0  ! Mach number units instead,
                             ! rather than cm/s
      INDPRO(4)=IABS(INDPRO(4))

cdr  default is: |indpro| < 10:  nplsi flow fields, one per background species
cdr              |indpro| > 10:  only one common flow field

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
                        !  despite storage for NPLS V.IN profiles.
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
        MPLSV = 1                    ! find ALL flow fields
                                     ! on IPLS=1 storage
      ELSEIF (NPLSV == NPLSI) THEN
        MPLSV = (/ (I,I=1,NPLSI) /)  ! find each individual flow field
                                     ! on its IPLS storage
      ELSE
cdr     more general options: to be written.
cdr     E.g.: only virtual background species with different flow fields
      ENDIF
      IF ((INDPRO(4).LE.5).AND.(NPLSI.GT.0)) THEN
        call json%info(pvel,n_children=nch)
        if (LRDMLV.and.(nch /= nplsi)) then
          write (iunout,*) ' ERROR READING FLOW VELOCITY PROFILES '
          write (iunout,*) ' TOO FEW VELOCITY PROFILES PROVIDED '
          write (iunout,*) ' NO. PROFILES FOUND ',nch
          write (iunout,*) ' NPLSI = ',NPLSI
          call eirene_exit_own(1)
        end if
        if (.not.LRDMLV .and. (nch > 1)) nch = 1
        do i = 1,nch
          call json%get_child(pvel,i,pch,found)

          call json%get(pch,'VX0',vx0(i),found)
          call json%get(pch,'VX1',vx1(i),found)
          call json%get(pch,'VX2',vx2(i),found)
          call json%get(pch,'VX3',vx3(i),found)
          call json%get(pch,'VX4',vx4(i),found)
          call json%get(pch,'VX5',vx5(i),found)

          call json%get(pch,'VY0',vy0(i),found)
          call json%get(pch,'VY1',vy1(i),found)
          call json%get(pch,'VY2',vy2(i),found)
          call json%get(pch,'VY3',vy3(i),found)
          call json%get(pch,'VY4',vy4(i),found)
          call json%get(pch,'VY5',vy5(i),found)

          call json%get(pch,'VZ0',vz0(i),found)
          call json%get(pch,'VZ1',vz1(i),found)
          call json%get(pch,'VZ2',vz2(i),found)
          call json%get(pch,'VZ3',vz3(i),found)
          call json%get(pch,'VZ4',vz4(i),found)
          call json%get(pch,'VZ5',vz5(i),found)
          nullify(pch)
        end do
        if (.not.lrdmlv) then
          vx0(2:nplsti) = vx0(1)
          vx1(2:nplsti) = vx1(1)
          vx2(2:nplsti) = vx2(1)
          vx3(2:nplsti) = vx3(1)
          vx4(2:nplsti) = vx4(1)
          vx5(2:nplsti) = vx5(1)

          vy0(2:nplsti) = vy0(1)
          vy1(2:nplsti) = vy1(1)
          vy2(2:nplsti) = vy2(1)
          vy3(2:nplsti) = vy3(1)
          vy4(2:nplsti) = vy4(1)
          vy5(2:nplsti) = vy5(1)

          vz0(2:nplsti) = vz0(1)
          vz1(2:nplsti) = vz1(1)
          vz2(2:nplsti) = vz2(1)
          vz3(2:nplsti) = vz3(1)
          vz4(2:nplsti) = vz4(1)
          vz5(2:nplsti) = vz5(1)
        end if
      END IF
      nullify (pvel)

c  pitch - or B field profile
c                              !  default:
                               !  B FIELD WITH BX=0
      NLPITCH=INDPRO(5).LT.0   !  for 1D parallel B runs:
                               !  B FIELD WITH BY=0

      INDPRO(5)=IABS(INDPRO(5))
      IF (INDPRO(5).LE.5) THEN
        call json%get(pbfld,'B0',b0,found)
        call json%get(pbfld,'B1',b1,found)
        call json%get(pbfld,'B2',b2,found)
        call json%get(pbfld,'B3',b3,found)
        call json%get(pbfld,'B4',b4,found)
        call json%get(pbfld,'B5',b5,found)
      END IF
      nullify (pbfld)

c  zi profiles
      IF ((INDPRO(11) > 0) .AND. (INDPRO(3).LE.5)) THEN
        call json%info(pzi,n_children=nch)
        if (nch /= nplsi) then
          write (iunout,*) ' ERROR READING ZI PROFILES '
          write (iunout,*) ' NO. OF PROFILES DOES NOT MATCH',
     .                     ' NO. OF BULK SPECIES'
          write (iunout,*) ' NO. PROFILES FOUND ',nch
          write (iunout,*) ' NPLSI = ',NPLSI
          call eirene_exit_own(1)
        end if
        do i = 1,nch
          call json%get_child(pzi,i,pch,found)
          call json%get(pch,'ZI0',zi0(i),found)
          call json%get(pch,'ZI1',zi1(i),found)
          call json%get(pch,'ZI2',zi2(i),found)
          call json%get(pch,'ZI3',zi3(i),found)
          call json%get(pch,'ZI4',zi4(i),found)
          call json%get(pch,'ZI5',zi5(i),found)
          nullify (pch)
        end do
      END IF
      nullify (pzi)

c  cell volume profile card, OPTIONAL

      IF ((INDPRO(12) > 0) .AND. (INDPRO(12).LE.5)) THEN
        if (associated(pvol)) then
          call json%get(pvol,'VL0',vl0,found)
          call json%get(pvol,'VL1',vl1,found)
          call json%get(pvol,'VL2',vl2,found)
          call json%get(pvol,'VL3',vl3,found)
          call json%get(pvol,'VL4',vl4,found)
          call json%get(pvol,'VL5',vl5,found)
          nullify(pvol)
        else
          WRITE (iunout,*) 'ONE INPUT LINE (VOL) MISSING IN BLOCK 5'
          WRITE (iunout,*) 'AUTOMATIC CORRECTION PERFORMED'
          CALL EIRENE_LEER(1)
          vl0 = 0
        end if
      END IF

c  check for "OPTIONAL" input cards.
c  Here: explicitly switch on/off input tallies,
c        or activate gradient tallies: INTLOPTS(ITAL)
c  ITAL:  -ntali<ital<0, then: ital=iabs(ital).
c  IOPT:   turn off or add gradient tally

      if (foundi) then
cdr  read an arbitrary number of tally-option cards ITAL, IOPT
cdr  and store results on INTLOPTS(ITAL).
cdr  Above some switches LTESMO, LDESMO, etc.. have already been set from INDPRO > 100.
cdr  Possible conflicts, inconsistencies.
        call json%info(pintl,n_children=nch)
        if ((nch == 0) .or.(nch > ntali)) then
          write (iunout,*) ' ERROR READING OPTIONS FOR INPUT TALLIES'
          write (iunout,*) ' NO. OF AVAILABLE OPTIONS IS INVALID'
          write (iunout,*) ' NO. OPTIONS FOUND ',nch
          write (iunout,*) ' NTALI = ',NTALI
          call eirene_exit_own(1)
        end if
        do i = 1,nch
          call json%get_child(pintl,i,pch,found)
          call json%get(pch,'ITAL',ital,found)
          call json%get(pch,'IOPT',iopt,found)
          nullify (pch)
          IF ((ITAL >= 0) .OR. (ITAL < -NTALI)) CYCLE
          INTLOPTS(IABS(ITAL)) = IOPT
       end do
        nullify (pintl)
      end if

!db read AMDoutpars

      call json%get_child(p,'AMDoutpars',pAMD,found)
      if(found) then
!cym/cpg type mismatch
!        call json%get(pAMD,'AMDoutput',AMDOutPars%AMDoutput,found)
         call json%get(pAMD,'AMDoutput',AMDoutputCK,found)
         AMDOutPars%AMDoutput=AMDoutputCK
!cym/cpg end
         call json%get(pAMD,'NEmin',AMDOutPars%NEmin,found)
         call json%get(pAMD,'NEmax',AMDOutPars%NEmax,found)
         call json%get(pAMD,'TEmin',AMDOutPars%TEmin,found)
         call json%get(pAMD,'TEmax',AMDOutPars%TEmax,found)
         call json%get(pAMD,'NNE',AMDOutPars%NNE,found)
         call json%get(pAMD,'NTE',AMDOutPars%NTE,found)
!cym/cpg
!         call json%get(pAMD,'MeshType',AMDOutPars%MeshType,found)
         call json%get(pAMD,'MeshType',MeshTypeCK,found)
         AMDOutPars%MeshType=MeshTypeCK
         nullify (pAMD)
      endif

      nullify(plsm)
      nullify(pbck)

      return
      end subroutine eirene_read_block_5

!*****************************************************************

      subroutine eirene_read_block_6(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: ptom, psmods, pchild, pchan, pch
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable ::
     . txt, varname, spcname
      character(kind=json_CK,len=72), dimension(:), allocatable ::
     . rfilnm
!cym/cpg end
      integer :: iflr, nmods, imod, nch, i2, nfr, ideflt_sput, i, l,
     .           ispz, ico, is, ideflt_spez, jflr, j
      real(dp), allocatable :: rhelp(:)
      logical :: found, foundp, founds, foundc, lf(5)
      CHARACTER(500) :: FILE
      external :: eirene_leer, eirene_masage, eirene_exit_own

      CALL EIRENE_MASAGE
     .  ('*** 6. GENERAL DATA FOR REFLECTION MODEL')
      CALL EIRENE_LEER(1)

      call json%get(p,'NLTRIM',nltrim,found)

      nfr = 0
      if (nltrim) then

c  read TRIM reflection datasets A_on_B
        call json%get(p,'CPATH',txt,foundp)

        if (.not.foundp) then
C  NO PATH SPECIFIED FOR REFLECTION DATA BASE
          WRITE (iunout,*)
     .      ' NO PATH SPECIFIED FOR TRIM REFLECTION DATA BASE'
          WRITE (iunout,*) ' OLD TRIM DATABASE VERSION USED'
          CALL EIRENE_LEER(1)
          LTRIM_OLD=.TRUE.
C  TAKE OLD "TRIM.DAT" FILE WITH NHD6=12 TARGET-PROJECTILE COMBINATIONS
C           "TRIM.DAT" IS EXPECTED ON INPUT STREAM IUN=(21+IFOFF) (SUBR. REFDAT.F)
C  SKIP LINES CONTAINING SPECIFICATIONS FOR DATA BASES AND
C  CONTINUE READING with DATD

c  old default: all TRIM data in one single file.
c  read this (single) file "TRIM.DAT" in subr. REFDAT
          NHD6 = 12

        else
C  PATH SPECIFICATION FOR DATA BASE FOUND
          LTRIM_OLD=.FALSE.
          if (len_trim(adjustl(txt)) > 400) then
            WRITE (iunout,*) ' PATH FOR TRIM REFLECTION DATABASE IS',
     .                  ' TOO LONG !'
            CALL EIRENE_EXIT_OWN(1)
          END IF
!VK TO REMOVE \r
          I=INDEX(TXT,ACHAR(13)) !VK
          IF(I.GT.0.AND.I.LE.LEN_TRIM(TXT)) TXT(I:I)=' '
          I2=LEN(TXT)
          FILE=REPEAT(' ',420)
          FILE(1:I2)=txt
          WRITE (iunout,*) ' PATH = ',FILE(1:I2)
C  PATH FOUND. NEXT: READ ONE OR MORE CARDS FILNAM A_ON_B
          call json%get_child(p,'PROJ_ON_MAT',ptom,found)
          call json%info(ptom,n_children=nfr)
          if (nfr > 0) then
c  TRIM files for NFR target projectile combinations are requested.
c  read them one by one in subr. RDTRIM
c  RFILNM is array of strings and holds NFR filenames
            call json%get(ptom,rfilnm)
            do j=1, nfr
              I=INDEX(RFILNM(J),ACHAR(13)) !VK
              IF(I.GT.0.AND.I.LE.LEN_TRIM(RFILNM(J))) RFILNM(J)(I:I)=' ' !VK TO REMOVE \r
            end do
            nhd6 = nfr
          else
c  old default: all TRIM data in one single file.
c  read this (single) file "TRIM.DAT" in subr. REFDAT
            nhd6 = 12
          end if
        end if
        nullify(ptom)
        deallocate(txt)
      end if  ! nltrim

!  READING OF REFLECTION DATASETS 'A_ON_B' COMPLETED

      NH0=NHD1*NHD2*NHD6
      NH1=NH0*NHD3
      NH2=NH1*NHD4
      NH3=NH2*NHD5
      CALL EIRENE_ALLOC_CREF

      if (allocated(rfilnm)) then
        do iflr = 1, nfr
          file(i2+1:) = rfilnm(iflr)
          WRITE (iunout,*) ' NFR =',IFLR,' FILNAM = ',FILE(I2+1:I2+20)
          JFLR = NFR - IFLR + 1
          REFFIL(JFLR) = FILE
        end do
        deallocate (rfilnm)
        IF (NFR.GE.1) CALL EIRENE_LEER(1)
      end if

c  next: read species index sampling distributions datm, dmol, dion, dpls, and in case nphot > 0, also dphot

      call json%get(p,'DATD',rhelp,found)
      datd(1:natmi) = rhelp(1:natmi)
      deallocate(rhelp)

      call json%get(p,'DMLD',rhelp,found)
      dmld(1:nmoli) = rhelp(1:nmoli)
      deallocate(rhelp)

      call json%get(p,'DIOD',rhelp,found)
      diod(1:nioni) = rhelp(1:nioni)
      deallocate(rhelp)

      call json%get(p,'DPLD',rhelp,found)
      dpld(1:nplsi) = rhelp(1:nplsi)
      deallocate(rhelp)

      if (nphoti > 0) then   !dr try to make this more logical:
                             !dr always read dphd.
                             !dr backward compatible ?
        call json%get(p,'DPHD',rhelp,found)
        dphd(1:nphoti) = rhelp(1:nphoti)
        deallocate(rhelp)
      end if

      datd_in = datd
      dmld_in = dmld
      diod_in = diod
      dpld_in = dpld
      dphd_in = dphd

c  next: read universal surface reflection model flags

      call json%get(p,'ERMIN',ermin,found)
      call json%get(p,'ERCUT',ercut,found)
      call json%get(p,'RPROB0',rprob0,found)
      call json%get(p,'RINTEG',rinteg(1),found)
      call json%get(p,'EINTEG',einteg(1),found)
      call json%get(p,'AINTEG',ainteg(1),found)

c  read surface models identified by character string 'SURFMOD_...'

      call json%get_child(p,'SURFMODS',psmods,founds)
      if (founds) then

        call json%info(psmods,n_children=nmods)

        do imod = 1, nmods
          ALLOCATE (REFCUR)
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
          REFCUR%RCYCSR  = 0._DP
          REFCUR%RCYCCR  = 0._DP
          REFCUR%STPRMR  = 0._DP
          REFCUR%ESPTSR  = 0._DP
          REFCUR%ESPTCR  = 0._DP

          call json%get_child(psmods,imod,pchild,found)

          call json%get(pchild,'NAME',txt,found)
          REFCUR%REFNAME = TRIM(ADJUSTL(TXT(9:)))

          call json%get(pchild,'ILREF',REFCUR%JLREF(1), found)
          call json%get(pchild,'ILSPT',REFCUR%JLSPT, found)
          call json%get(pchild,'ISRS',REFCUR%JSRS(1), found)
          call json%get(pchild,'ISRC',REFCUR%JSRC(1), found)
          call json%get(pchild,'LCHSPNWL',REFCUR%JLCHSPNWL(1), found)
          call json%get(pchild,'REFCELL',REFCUR%REFCELL, found)

          call json%get(pchild,'ZNML',REFCUR%ZNMLR, found)
          call json%get(pchild,'EWALL',REFCUR%EWALLR, found)
          call json%get(pchild,'EWBIN',REFCUR%EWBINR, found)
          call json%get(pchild,'TRANSP',rhelp,found)
          REFCUR%TRANSPR(1,1:2) = rhelp(1:2)
          deallocate (rhelp)
          call json%get(pchild,'FSHEAT',REFCUR%FSHEATR, found)

          call json%get(pchild,'RECYCF',REFCUR%RCYCFR(1), found)
          call json%get(pchild,'RECYCT',REFCUR%RCYCTR(1), found)
          call json%get(pchild,'RECPRM',REFCUR%RCPRMR(1), found)
          call json%get(pchild,'EXPPL',REFCUR%EXPPLR(1), found)
          call json%get(pchild,'EXPEL',REFCUR%EXPELR(1), found)
          call json%get(pchild,'EXPIL',REFCUR%EXPILR(1), found)

          REFCUR%JLREF(2:NSPZ) = REFCUR%JLREF(1)
          REFCUR%JSRS(2:NSPZ) = REFCUR%JSRS(1)
          REFCUR%JSRC(2:NSPZ) = REFCUR%JSRC(1)
          REFCUR%JLCHSPNWL(2:NSPZ) = REFCUR%JLCHSPNWL(1)    !VK
          REFCUR%TRANSPR(2:NSPZ,1)=REFCUR%TRANSPR(1,1)
          REFCUR%TRANSPR(2:NSPZ,2)=REFCUR%TRANSPR(1,2)
          REFCUR%RCYCFR(2:NSPZ) = REFCUR%RCYCFR(1)
          REFCUR%RCYCTR(2:NSPZ) = REFCUR%RCYCTR(1)
          REFCUR%RCPRMR(2:NSPZ) = REFCUR%RCPRMR(1)
          REFCUR%EXPPLR(2:NSPZ) = REFCUR%EXPPLR(1)
          REFCUR%EXPELR(2:NSPZ) = REFCUR%EXPELR(1)
          REFCUR%EXPILR(2:NSPZ) = REFCUR%EXPILR(1)
C
C  DEFAULT SPUTER MODEL
          REFCUR%RCYCSR=RECYCS(1,0)
          REFCUR%RCYCCR=RECYCC(1,0)
          REFCUR%STPRMR=SPTPRM(1,0)
          REFCUR%ESPTSR=ESPUTS(1,0)
          REFCUR%ESPTCR=ESPUTC(1,0)

!  check for sputer model in input
          ideflt_sput = -1
          call json%get(pchild,'RECYCS',REFCUR%RCYCSR(1), lf(1))
          if (lf(1)) REFCUR%RCYCSR(2:NSPZ) = REFCUR%RCYCSR(1)

          call json%get(pchild,'RECYCC',REFCUR%RCYCCR(1), lf(2))
          if (lf(2)) REFCUR%RCYCCR(2:NSPZ) = REFCUR%RCYCCR(1)

          call json%get(pchild,'SPTPRM',REFCUR%STPRMR(1), lf(3))
          if (lf(3)) REFCUR%STPRMR(2:NSPZ) = REFCUR%STPRMR(1)

          call json%get(pchild,'ESPUTS',REFCUR%ESPTSR(1), lf(4))
          if (lf(4)) REFCUR%ESPTSR(2:NSPZ) = REFCUR%ESPTSR(1)

          call json%get(pchild,'ESPUTC',REFCUR%ESPTCR(1), lf(5))
          if (lf(5)) REFCUR%ESPTCR(2:NSPZ) = REFCUR%ESPTCR(1)

          call json%get(pchild,'REFPRESS',REFCUR%REFPRESS, found)

          if (any(lf)) then
            DO I=2,NSPZ
              REFCUR%RCYCSR(I) = REFCUR%RCYCSR(1)
              REFCUR%RCYCCR(I) = REFCUR%RCYCCR(1)
              REFCUR%STPRMR(I) = REFCUR%STPRMR(1)
              REFCUR%ESPTSR(I) = REFCUR%ESPTSR(1)
              REFCUR%ESPTCR(I) = REFCUR%ESPTCR(1)
            ENDDO
            ideflt_sput = 1
!Pressure feedback loop model
            IF (REFCUR%JLREF(1) == 4)
     .         WRITE(IUNOUT,*) "PFL with parameters: ",
     .           REFCUR%REFCELL, REFCUR%REFPRESS
          end if

!  check of non-default sputter model
          call json%get_child(pchild,'CHANGES',pchan,foundc)
          if (foundc) then
            call json%info(pchan,n_children=nch)
            do i = 1, nch
              call json%get_child(pchan,i,pch)
              call json%get(pch,'VARIABLE',varname,found)
              call json%get(pch,'SPECIES',spcname,found)

!  find species index
              ispz=-1
              ico=0
              do is=1,nspz
                if (texts(is) == spcname(1:8)) then
                  if (ico == 0) then
                    ispz = is
                    ico = 1
                  else
                    call EIRENE_leer(2)
                    write (iunout,*) ' WARNING !! '
                    write (iunout,*)
     .                ' ambigous species names found in surface ',
     .                'model "',REFCUR%REFNAME,'"'
                    write (iunout,*) varname,
     .                'is modified for species ',ispz
                  end if
                end if
              end do             ! end of "is"-loop
              if (ispz < 0) then
                write (iunout,*) ' wrong species encountered in',
     .                           ' species dep. reflection model'
                write (iunout,'(a,1x,a)') varname, spcname
                cycle
              end if

!  here we know that ispz > 0
              ideflt_spez=1     !tentatively assume: species card
              allocate(spr)
              spr%varname = varname
              spr%spcname = spcname
              if ((varname == 'ISRS') .or. (varname == 'ISRC') .or.
     .            (varname == 'LCHSPNWL')) then
                call json%get(pch,'VALUE',spr%ival,found)
              else
                call json%get(pch,'VALUE',spr%rval,found)
              end if
!             spr%next => refcur%spec_lines
!             refcur%spec_lines => spr
              nullify(spr%next)
              if (.not.associated(refcur%spec_lines)) then
                refcur%spec_lines => spr
                spr_last => spr
              else
                spr_last%next => spr
                spr_last => spr_last%next
              end if

!cym/cpg avoiding type mismatch
!  'xxxx' -> json_ck_'xxxx' consistent with varname type
              select case (varname)
              case (json_ck_'ISRS')
                REFCUR%JSRS(ispz) = spr%ival
              case (json_ck_'ISRC')
                REFCUR%JSRC(ispz) = spr%ival
              case (json_ck_'LCHSPNWL')
                REFCUR%JLCHSPNWL(ispz) = spr%ival !VK
              case (json_ck_'TRANSP1')
                REFCUR%TRANSPR(ispz,1) = spr%rval
              case (json_ck_'TRANSP2')
                REFCUR%TRANSPR(ispz,2) = spr%rval
              case (json_ck_'RECYCF')
                REFCUR%RCYCFR(ispz) = spr%rval
              case (json_ck_'RECYCT')
                REFCUR%RCYCTR(ispz) = spr%rval
              case (json_ck_'RECPRM')
                REFCUR%RCPRMR(ispz) = spr%rval
              case (json_ck_'EXPPL')
                REFCUR%EXPPLR(ispz) = spr%rval
              case (json_ck_'EXPEL')
                REFCUR%EXPELR(ispz) = spr%rval
              case (json_ck_'EXPIL')
                REFCUR%EXPILR(ispz) = spr%rval
              case (json_ck_'RECYCS')
                REFCUR%RCYCSR(ispz) = spr%rval
              case (json_ck_'RECYCC')
                REFCUR%RCYCCR(ispz) = spr%rval
              case (json_ck_'SPTPRM')
                REFCUR%STPRMR(ispz) = spr%rval
              case (json_ck_'ESPUTS')
                REFCUR%ESPTSR(ispz) = spr%rval
              case (json_ck_'ESPUTC')
                REFCUR%ESPTCR(ispz) = spr%rval
!cym/cpg end

              case default
                write (iunout,*) ' unknown variable name encountered',
     .                           ' in species dep. reflection model'
                write (iunout,*) varname
              end select

              deallocate(varname)
              deallocate(spcname)

            end do  ! nch
          end if

C all lines of this particular SURFMOD_... are read now
C
          IF (ideflt_sput.le.0.and.REFCUR%JLSPT.NE.0) THEN
            WRITE (iunout,*) 'WARNING: SPUTTERING FOR MODEL ',
     .            TRIM(REFCUR%REFNAME)
            WRITE (iunout,*)
     .            'BUT NO PARAMETERS RECYCS, RECYCC ARE READ '
            WRITE (iunout,*) 'DEFAULT MODEL: "NO SPUTTERING" IS USED. '
            WRITE (iunout,*) 'DO YOU REALLY WANT THIS?'
            REFCUR%JLSPT=0
          ENDIF
          L=LEN_TRIM(REFCUR%REFNAME)
          WRITE (iunout,*) 'SURFACE MODEL "',REFCUR%REFNAME(1:L),
     .                     '" DEFINED'
!         REFCUR%NEXT => REFLIST
!         REFLIST => REFCUR
!  add new reflection model at the end of the list
          nullify(refcur%next)
          if (.not.associated(reflist)) then
            reflist => refcur
            reflast => refcur
          else
            reflast%next => refcur
            reflast => reflast%next
          end if

        end do  ! nmods
        CALL EIRENE_LEER(1)
      end if

      return
      end subroutine eirene_read_block_6

!*****************************************************************

      subroutine eirene_read_block_7(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pstrata, pstr, psubs, psub
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      logical :: found
      integer, allocatable :: ihelp(:)
      integer :: nstr, j
      external :: eirene_leer, eirene_masage, eirene_masprm
C
C  READ DATA FOR PRIMARY SOURCE
C
      CALL EIRENE_MASAGE
     .  ('*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA')

      call json%get(p,'NSTRAI',nstrai,found)
      WRITE (iunout,*) '       NSTRAI= ',NSTRAI
      NSTRAI_IN = NSTRAI
      CALL EIRENE_LEER(1)

      call json%get(p,'INDSRC',ihelp,found)
      indsrc(1:nstrai) = ihelp(1:nstrai)

      call json%get(p,'ALLOC',alloc,found)
      call json%get(p,'AMPTS',ampts,found)

C  AMPTS: (option added 2014) common multiplier for max. allowed cpu time NTCPU,
C          and for number of histories NPTS (see below)
      IF (AMPTS.EQ.0.0_DP) AMPTS=1.0_DP
      IF (AMPTS.NE.1.0_DP) THEN
        MPTS_COMSOU=AMPTS
        NTCPU=INT(REAL(NTCPU)*AMPTS)
        WRITE (iunout,*) '       NTCPU, NPTS ENHANCED BY FACTOR AMPTS ',
     .                     AMPTS
      ELSE
        MPTS_COMSOU=1.0_DP
      END IF

      call json%get_child(p,'STRATA',pstrata,found)
      call json%info(pstrata,n_children=nstr)

      do istr = 1, nstr

        call json%get_child(pstrata,istr,pstr,found)

        call json%get(pstr,'ISTRA',jstra,found)

        call json%get(pstr,'TXTSOU',txt,found)
        TXTSOU(JSTRA) = TXT
        deallocate(txt)

        call json%get(pstr,'NLAVRP',nlavrp(jstra),found)
        call json%get(pstr,'NLAVRT',nlavrt(jstra),found)
        call json%get(pstr,'NLSYMP',nlsymp(jstra),found)
        call json%get(pstr,'NLSYMT',nlsymt(jstra),found)
!       call json%get(pstr,'NLRAY',nlray(jstra),found)
        NLRAY(JSTRA)=.FALSE.        ! CDR: UNFINISHED PROPRIETARY OPTION

        call json%get(pstr,'NPTS',npts(jstra),found)
        call json%get(pstr,'NINITL',ninitl(jstra),found)
        call json%get(pstr,'NEMODS',nemods(jstra),found)
        call json%get(pstr,'NAMODS',namods(jstra),found)

c.............................................
c  NMINPTS:   ENFORCE MINIMUM NUMBER OF HISTORIES FOR STRATUM.
cdr           currently set in various places, but not used anywhere.
        call json%get(pstr,'NMINPTS',nminpts(jstra),found)

c  NPTSDEL: =0   (VALUES NPTSDEL .GT.0 ONLY FOR INTERNAL TESTTING PROCEDURES
C                 OF MPI parallelization. Requires 2 runs with special input settings
cdr See comments in MCARLO.F
        call json%get(pstr,'NPTSDEL',nptsdel(jstra),found)
!       call json%get(pstr,'NRAYEN',nrayen(jstra),found)   !CDR  NOT IN USE
c...................................................................

C APPLY MULTIPLIER AMPTS TO SELECTED MC PARTICLE NUMBER RANGE
        IF(AMPTS.GT.0) THEN
          NPTS(JSTRA)=INT(NPTS(JSTRA)*AMPTS)
          NMINPTS(JSTRA)=INT(NMINPTS(JSTRA)*AMPTS)
        END IF
C

csw if npts < 0 --> npts = infty
        if(npts(jstra).lt.0) npts(jstra) = huge(1)-1

        call json%get(pstr,'FLUX',flux(jstra),found)
        call json%get(pstr,'SCALV',scalv(jstra),found)
        call json%get(pstr,'IVLSF',ivlsf(jstra),found)
        call json%get(pstr,'ISCLS',iscls(jstra),found)
        call json%get(pstr,'ISCLT',isclt(jstra),found)
        call json%get(pstr,'ISCL1',iscl1(jstra),found)
        call json%get(pstr,'ISCL2',iscl2(jstra),found)
        call json%get(pstr,'ISCL3',iscl3(jstra),found)
        call json%get(pstr,'ISCLB',isclb(jstra),found)
        call json%get(pstr,'ISCLA',iscla(jstra),found)

        call json%get(pstr,'NLATM',nlatm(jstra),found)
        call json%get(pstr,'NLMOL',nlmol(jstra),found)
        call json%get(pstr,'NLION',nlion(jstra),found)
        call json%get(pstr,'NLPLS',nlpls(jstra),found)
        call json%get(pstr,'NLPHOT',nlphot(jstra),found)

        call json%get(pstr,'NSPEZ',nspez(jstra),found)

        call json%get(pstr,'NLPNT',nlpnt(jstra),found)
        call json%get(pstr,'NLLNE',nllne(jstra),found)
        call json%get(pstr,'NLSRF',nlsrf(jstra),found)
        call json%get(pstr,'NLVOL',nlvol(jstra),found)
        call json%get(pstr,'NLCNS',nlcns(jstra),found)

C  SAME FOR POINT, LINE, SURFACE AND VOLUME SOURCES

        call json%get(pstr,'NSRFSI',nsrfsi(jstra),found)
        IF (NSRFSI(JSTRA).GT.NSRFS)
     .    CALL EIRENE_MASPRM(
     .         'NSRFS',5,NSRFS,'NSRFSI(I)',9,NSRFSI(I),IERROR)

        call json%get_child(pstr,'SUBSTRATA',psubs, found)

        do j = 1, nsrfsi(jstra)
          call json%get_child(psubs,j,psub,found)

          call json%get(psub,'INDIM',indim(j,jstra),found)
          call json%get(psub,'INSOR',insor(j,jstra),found)
          call json%get(psub,'INGRDA',ihelp,found)
          if (found) then
            ingrda(j,jstra,1:3) = ihelp(1:3)
            deallocate(ihelp)
          end if
          call json%get(psub,'INGRDE',ihelp,found)
          if (found) then
            ingrde(j,jstra,1:3) = ihelp(1:3)
            deallocate(ihelp)
          end if

          call json%get(psub,'SORWGT',sorwgt(j,jstra),found)
          call json%get(psub,'SORLIM',sorlim(j,jstra),found)
          call json%get(psub,'SORIND',sorind(j,jstra),found)
          call json%get(psub,'SOREXP',sorexp(j,jstra),found)
          call json%get(psub,'SORIFL',sorifl(j,jstra),found)

          call json%get(psub,'NRSOR',nrsor(j,jstra),found)
          call json%get(psub,'NPSOR',npsor(j,jstra),found)
          call json%get(psub,'NTSOR',ntsor(j,jstra),found)
          call json%get(psub,'NBSOR',nbsor(j,jstra),found)
          call json%get(psub,'NASOR',nasor(j,jstra),found)
          call json%get(psub,'NISOR',nisor(j,jstra),found)
          call json%get(psub,'ISTOR',istor(j,jstra),found)

          call json%get(psub,'SORAD1',sorad1(j,jstra),found)
          call json%get(psub,'SORAD2',sorad2(j,jstra),found)
          call json%get(psub,'SORAD3',sorad3(j,jstra),found)
          call json%get(psub,'SORAD4',sorad4(j,jstra),found)
          call json%get(psub,'SORAD5',sorad5(j,jstra),found)
          call json%get(psub,'SORAD6',sorad6(j,jstra),found)

          nullify(psub)
        end do
        nullify(psubs)

C  VELOCITY SPACE DISTRIBUTION

        call json%get(pstr,'SORENI',soreni(jstra),found)
        call json%get(pstr,'SORENE',sorene(jstra),found)
        call json%get(pstr,'SORVDX',sorvdx(jstra),found)
        call json%get(pstr,'SORVDY',sorvdy(jstra),found)
        call json%get(pstr,'SORVDZ',sorvdz(jstra),found)

        call json%get(pstr,'SORCOS',sorcos_in(jstra),found)
        call json%get(pstr,'SORMAX',sormax_in(jstra),found)
        call json%get(pstr,'SORCTX',sorctx(jstra),found)
        call json%get(pstr,'SORCTY',sorcty(jstra),found)
        call json%get(pstr,'SORCTZ',sorctz(jstra),found)
        call json%get(pstr,'RAYFRAC',rayfrac(jstra),found)

        sorcos(jstra) = sorcos_in(jstra)
        sormax(jstra) = sormax_in(jstra)

        nullify(pstr)
      end do

      nullify(pstrata)

      return
      end subroutine eirene_read_block_7

!*****************************************************************

      subroutine eirene_read_block_8(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pzones, pz, pt, pd, pv, pvol
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=72), dimension(:), allocatable ::
     .           ch3lines
!cym/cpg end
      logical :: found, foundc, foundv, foundm
      integer :: inilgj, inelgj, idion, ini, ine, nch3
      real(dp) :: tte, tti, di, vx, vy, vz, vl
      external :: eirene_dekey, eirene_dekeyb, eirene_exit_own,
     .            eirene_leer, eirene_masage
C
C     READ ADDITIONAL DATA FOR SOME SPECIFIC ZONES
C

      NZADD = 0
      NULLIFY(TEMPLIST)
      NULLIFY(DENLIST)
      NULLIFY(VELLIST)
      NULLIFY(VOLLIST)

      CALL EIRENE_MASAGE ('*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES')

      if (.not.associated(p)) return

      call json%get(p,'NZADD',nzadd,found)
      WRITE (iunout,*) '       NZADD= ',NZADD
      CALL EIRENE_LEER(1)

      if (nzadd == 0) return

      ALLOCATE (INI_ZONE(NZADD))
      ALLOCATE (INE_ZONE(NZADD))
      ALLOCATE (CH3_STACK(NZADD))
      INI_ZONE = 0
      INE_ZONE = 0

      call json%get_child(p,'ZONE',pzones,found)

      do i = 1, nzadd
        call json%get_child(pzones,i,pz,found)

        call json%get(pz,'INI',ini,found)
        call json%get(pz,'INE',ine,found)
        INI_ZONE(I) = INI
        INE_ZONE(I) = INE
        IF (INI.GT.NRAD.OR.INI.LE.0) GOTO 998
        IF (INE.GT.NRAD) GOTO 998
        IF (INE.LE.0) INE=INI

C  IGJUM3 FLAG
        NULLIFY(CH3_STACK(I)%HEAD)
        NULLIFY(CH3_STACK(I)%LAST)
        call json%get(pz,'CH3-LINES',ch3lines,foundc)
        if (foundc) then
          nch3 = size(ch3lines)
          INILGJ=MAX(1,MIN(NOPTIM,INI))
          if (inilgj.ne.ini) write (iunout,*)
     .      'CH3 option failed. ini, noptim =',ini,NOPTIM
          INELGJ=MIN(NOPTIM,INE)
          if (inelgj.ne.ine) write (iunout,*)
     .      'CH3 option failed. ine, noptim =',ine,NOPTIM
          do k = 1, nch3
            DO II=INILGJ,INELGJ
              IF (NLIMPB >= NLIMPS) THEN
                CALL EIRENE_DEKEY (ch3lines(k)(4:72),IGJUM3,
     .                             0,NOPTIM,II,NLIMPS)
              ELSE
                CALL EIRENE_DEKEYB(ch3lines(k)(4:72),IGJUM3,
     .                             0,NOPTIM,II,NLIMPB,NBITS)
              END IF
            END DO
            call eirene_push_string_stack(ch3_stack(i),ch3lines(k))
          end do
          deallocate(ch3lines)
        end if

C  TEMPERATURE
        call json%get_child(pz,'T',pt,foundc)
        if (foundc) then
          call json%get(pt,'IDION',idion,found)
          call json%get(pt,'TE',tte,found)
          call json%get(pt,'TI',tti,found)
          do ii = ini,ine
            ALLOCATE(TEMPCUR)
            TEMPCUR%II = II
            TEMPCUR%IDION = IDION
            TEMPCUR%TE = TTE
            TEMPCUR%TI = TTI
            TEMPCUR%NEXT => TEMPLIST
            TEMPLIST => TEMPCUR
          end do
        end if
        nullify(pt)

C  DENSITY
        call json%get_child(pz,'D',pd,foundc)
        if (foundc) then
          call json%get(pd,'IDION',idion,found)
          call json%get(pd,'DI',di,found)
          do ii = ini,ine
            ALLOCATE(DENCUR)
            DENCUR%II = II
            DENCUR%IDION = IDION
            DENCUR%DI = DI
            DENCUR%NEXT => DENLIST
            DENLIST => DENCUR
          end do
        end if
        nullify(pd)

C  VELOCITY (CM/SEC OR MACH)
        call json%get_child(pz,'V',pv,foundv)
        if (.not.foundv) call json%get_child(pz,'M',pv,foundm)
        if (foundv .or. foundm) then
          call json%get(pv,'IDION',idion,found)
          call json%get(pv,'VX',vx,found)
          call json%get(pv,'VY',vy,found)
          call json%get(pv,'VZ',vz,found)
          do ii = ini,ine
            ALLOCATE(VELCUR)
            VELCUR%II = II
            VELCUR%IDION = IDION
            VELCUR%VX = VX
            VELCUR%VY = VY
            VELCUR%VZ = VZ
            if (foundv) VELCUR%IZ = -1
            if (foundm) VELCUR%IZ = 1
            VELCUR%NEXT => VELLIST
            VELLIST => VELCUR
          end do
        end if
        nullify(pv)

C  VOLUME
        call json%get_child(pz,'VL',pvol,foundc)
        if (foundc) then
          call json%get(pvol,'VOL',vl,found)
          do ii = ini,ine
            ALLOCATE(VOLCUR)
            VOLCUR%II = II
            VOLCUR%VOL = VL
            VOLCUR%NEXT => VOLLIST
            VOLLIST => VOLCUR
          end do
        end if
        nullify(pvol)

        nullify(pz)
      end do
      nullify(pzones)

      return

998   CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT BLOCK FOR ADDITIONAL DATA FOR  '
      WRITE (iunout,*) 'SPECIFIC ZONES FOUND AT ZONE NO. ',I
      CALL EIRENE_EXIT_OWN(1)
      end subroutine eirene_read_block_8

!*****************************************************************

      subroutine eirene_read_block_9(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: psvi, psv, pssi, pss, psci, psc
      logical :: found
      logical, allocatable :: lhelp(:)
      integer :: num, i, ia, nvi, nsi, nci
      integer, allocatable :: iprsf(:), ihelp(:)
      real(dp), allocatable :: rhelp(:)
      external :: eirene_leer, eirene_masage, eirene_masprm,
     .            eirene_exit_own
C
C     READ DATA FOR STATISTICS AND NON-ANALOG MODEL, 900--999
C

      CALL EIRENE_MASAGE
     .  ('*** 9. DATA FOR STATISTIC AND NON-ANALOG MODEL')
      CALL EIRENE_LEER(1)

      maxlev = 0
      maxrad = 0
      maxpol = 0
      maxtor = 0
      maxadd = 0

      if (.not.associated(p)) return

      call json%get(p,'NLPRCA',lhelp,found)
      if (found) then
        nlprca(1:natmi) = lhelp(1:natmi)
        deallocate(lhelp)
      end if

      call json%get(p,'NLPRCM',lhelp,found)
      if (found) then
        nlprcm(1:nmoli) = lhelp(1:nmoli)
        deallocate(lhelp)
      end if

      call json%get(p,'NLPRCI',lhelp,found)
      if (found) then
        nlprci(1:nioni) = lhelp(1:nioni)
        deallocate(lhelp)
      end if

      call json%get(p,'NLPRCPH',lhelp,found)
      if (found) then
        nlprcph(1:nphoti) = lhelp(1:nphoti)
        deallocate(lhelp)
      end if

      call json%get(p,'NPRCSF',num,found)
      nlprcs = .false.
      if (num > 0) then
        call json%get(p,'IPRSF',iprsf,found)
        do i = 1, num
          if ((iprsf(i) >= 0) .and. (iprsf(i) <= nlimps))
     .       nlprcs(iprsf(i)) = .true.
        end do
        deallocate(iprsf)
      end if

C     DATA FOR SPLITTING AND RUSSIAN ROULETTE

      call json%get(p,'MAXLEV',maxlev,found)
      call json%get(p,'MAXRAD',maxrad,found)
      call json%get(p,'MAXPOL',maxpol,found)
      call json%get(p,'MAXTOR',maxtor,found)
      call json%get(p,'MAXADD',maxadd,found)

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

      nsspl = 0
      prmspl = 0._dp

      if (maxrad > 0) then
        call json%get(p,'NSSPL_RAD',ihelp,found)
        call json%get(p,'PRMSPL_RAD',rhelp,found)
        ia = 0
        nsspl(ia+1:ia+maxrad) = ihelp(1:maxrad)
        prmspl(ia+1:ia+maxrad) = rhelp(1:maxrad)
        deallocate(ihelp)
        deallocate(rhelp)
      end if

      if (maxpol > 0) then
        call json%get(p,'NSSPL_POL',ihelp,found)
        call json%get(p,'PRMSPL_POL',rhelp,found)
        ia = n1st
        nsspl(ia+1:ia+maxpol) = ihelp(1:maxpol)
        prmspl(ia+1:ia+maxpol) = rhelp(1:maxpol)
        deallocate(ihelp)
        deallocate(rhelp)
      end if

      if (maxtor > 0) then
        call json%get(p,'NSSPL_TOR',ihelp,found)
        call json%get(p,'PRMSPL_TOR',rhelp,found)
        ia = n1st + n2nd
        nsspl(ia+1:ia+maxtor) = ihelp(1:maxtor)
        prmspl(ia+1:ia+maxtor) = rhelp(1:maxtor)
        deallocate(ihelp)
        deallocate(rhelp)
      end if

      if (maxadd > 0) then
        call json%get(p,'NSSPL_ADD',ihelp,found)
        call json%get(p,'PRMSPL_ADD',rhelp,found)
        ia = n1st + n2nd + n3rd
        nsspl(ia+1:ia+maxadd) = ihelp(1:maxadd)
        prmspl(ia+1:ia+maxadd) = rhelp(1:maxadd)
        deallocate(ihelp)
        deallocate(rhelp)
      end if

C     DATA FOR BIAS SAMPLING
      call json%get(p,'WMINV',wminv,found)
      call json%get(p,'WMINS',wmins,found)
      call json%get(p,'WMINC',wminc,found)
      call json%get(p,'WMINL',wminl,found)
      WMINV=MAX(WMINV,EPS60)
      WMINS=MAX(WMINS,EPS60)
      WMINC=MAX(WMINC,EPS60)

      call json%get(p,'SPLPAR',splpar,found)

C  DATA FOR STANDARD DEVIATION
      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'

      call json%get(p,'NSIGVI',nsigvi,found)
      call json%get(p,'NSIGSI',nsigsi,found)
      call json%get(p,'NSIGCI',nsigci,found)
      call json%get(p,'NSIGI_SPC',nsigi_spc,found)
      WRITE (iunout,*) '       NSIGVI,NSIGSI,NSIGCI = ',
     .                         NSIGVI,NSIGSI,NSIGCI
      WRITE (iunout,*) '       NSIGI_SPC            = ',NSIGI_SPC
      CALL EIRENE_LEER(1)

      if (nsigvi > 0) then
        call json%get_child(p,'DEV_VOL_TAL',psvi,found)
        call json%info(psvi,n_children=nvi)
        if (nvi /= nsigvi) then
          write (iunout,*)
     .      ' NUMBER OF CARDS FOR STANDARD DEVIATION DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NSIGVI = ',nsigvi
          write (iunout,*) 'NVI =   ',nvi
          call eirene_exit_own(1)
        end if

        do i = 1, nsigvi
          call json%get_child(psvi,i,psv,found)
          call json%get(psv,'IGH',igh(i),found)
          call json%get(psv,'IIH',iih(i),found)
          nullify(psv)
        end do
        nullify(psvi)
      end if

      if (nsigsi > 0) then
        call json%get_child(p,'DEV_SRF_TAL',pssi,found)
        call json%info(pssi,n_children=nsi)
        if (nsi /= nsigsi) then
          write (iunout,*)
     .      ' NUMBER OF CARDS FOR STANDARD DEVIATION DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NSIGSI = ',nsigsi
          write (iunout,*) 'NSI =   ',nsi
          call eirene_exit_own(1)
        end if

        do i = 1, nsigsi
          call json%get_child(pssi,i,pss,found)
          call json%get(pss,'IGHW',ighw(i),found)
          call json%get(pss,'IIHW',iihw(i),found)
          nullify(pss)
        end do
        nullify(pssi)
      end if


      if (nsigci > 0) then
        call json%get_child(p,'CORR_COEFFS',psci,found)
        call json%info(psci,n_children=nci)
        if (nci /= nsigci) then
          write (iunout,*)
     .      ' NUMBER OF CARDS FOR CORRELATION COEFFS DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NSIGCI = ',nsigci
          write (iunout,*) 'NCI =   ',nci
          call eirene_exit_own(1)
        end if

        do i = 1, nsigci
          call json%get_child(psci,i,psc,found)
          call json%get(psc,'IGHC',ihelp,found)
          ighc(1:2,i) = ihelp(1:2)
          deallocate(ihelp)
          call json%get(psc,'IIHC',ihelp,found)
          iihc(1:2,i) = ihelp(1:2)
          deallocate(ihelp)
          nullify(pss)
        end do
        nullify(pssi)
      end if

      end subroutine eirene_read_block_9

!*****************************************************************

      subroutine eirene_read_block_10(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: padvi, pad, pclvi, pcl, palsi, pal,
     .                             pspeci, psp, palvi, padsi
      logical :: found
      integer :: i, nvi, nsi, j
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      external :: eirene_leer, eirene_masage, eirene_exit_own
C
C     READ DATA FOR ADDITIONAL AND SURFACE-AVERAGED TALLIES
C
      CALL EIRENE_MASAGE
     .  ('*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION')
      CALL EIRENE_MASAGE
     .  ('        ESTIMATORS AND ALGEBRAIC EXPRESSIONS')

      nadvi = 0
      nclvi = 0
      nalvi = 0
      nadsi = 0
      nalsi = 0
      nadspc = 0

      if (.not.associated(p)) then
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) '       NADVI,NCLVI,NALVI    = ',
     .                           NADVI,NCLVI,NALVI
        WRITE (iunout,*) '       NADSI,NALSI,NADSPC   = ',
     .                           NADSI,NALSI,NADSPC
        CALL EIRENE_LEER(1)
        ALLOCATE(ESTIML(1))
        return
      end if

      call json%get(p,'NADVI',nadvi,found)
      call json%get(p,'NCLVI',nclvi,found)
      call json%get(p,'NALVI',nalvi,found)
      call json%get(p,'NADSI',nadsi,found)
      call json%get(p,'NALSI',nalsi,found)
      call json%get(p,'NADSPC',nadspc,found)
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) '       NADVI,NCLVI,NALVI    = ',
     .                         NADVI,NCLVI,NALVI
      WRITE (iunout,*) '       NADSI,NALSI,NADSPC   = ',
     .                         NADSI,NALSI,NADSPC
      CALL EIRENE_LEER(1)

      CALL EIRENE_MASAGE('*** 10A. DATA FOR ADDITIONAL TALLIES')

      if (nadvi > 0) then
        call json%get_child(p,'ADD_VOL_TAL',padvi,found)
        call json%info(padvi,n_children=nvi)
        if (nvi /= nadvi) then
          write (iunout,*)
     .      ' NUMBER OF ADDITIONAL VOLUME AVERAGED TALLIES DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NADVI = ',nadvi
          write (iunout,*) 'NVI =   ',nvi
          call eirene_exit_own(1)
        end if

        do i = 1, nadvi
          call json%get_child(padvi,i,pad,found)
          call json%get(pad,'IADVE',iadve(i),found)
          call json%get(pad,'IADVS',iadvs(i),found)
          call json%get(pad,'IADVT',iadvt(i),found)
          call json%get(pad,'IADRC',iadrc(i),found)
          call json%get(pad,'TXTTLA',txt,found)
          TXTTLA(i) = TXT
          deallocate(txt)
          call json%get(pad,'TXTSCA',txt,found)
          TXTSCA(i) = TXT
          deallocate(txt)
          call json%get(pad,'TXTUTA',txt,found)
          TXTUTA(i) = TXT
          deallocate(txt)
          nullify(pad)
        end do
        nullify(padvi)
      end if

      CALL EIRENE_MASAGE
     .     ('*** 10B. DATA FOR COLLISION ESTIMATORS')

      if (nclvi > 0) then
        call json%get_child(p,'COLL_EST',pclvi,found)
        call json%info(pclvi,n_children=nvi)
        if (nvi /= nclvi) then
          write (iunout,*)
     .      ' NUMBER OF COLLISION ESTIMATORS DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NCLVI = ',nclvi
          write (iunout,*) 'NVI =   ',nvi
          call eirene_exit_own(1)
        end if

        do i = 1, nclvi
          call json%get_child(pclvi,i,pcl,found)
          call json%get(pcl,'ICLVE',iclve(i),found)
          call json%get(pcl,'ICLVS',iclvs(i),found)
          call json%get(pcl,'ICLVT',iclvt(i),found)
          call json%get(pcl,'ICLRC',iclrc(i),found)
          call json%get(pcl,'TXTTLC',txt,found)
          TXTTLC(i) = TXT
          deallocate(txt)
          call json%get(pcl,'TXTSCC',txt,found)
          TXTSCC(i) = TXT
          deallocate(txt)
          call json%get(pcl,'TXTUTC',txt,found)
          TXTUTC(i) = TXT
          deallocate(txt)
          nullify(pcl)
        end do
        nullify(pclvi)
      end if

      CALL EIRENE_MASAGE('*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS')

      if (nalvi > 0) then
        call json%get_child(p,'ALG_VOL_TAL',palvi,found)
        call json%info(palvi,n_children=nvi)
        if (nvi /= nalvi) then
          write (iunout,*)
     .      ' NUMBER OF ALGEBRAIC VOLUME AVERAGED TALLIES DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NALVI = ',nalvi
          write (iunout,*) 'NVI =   ',nvi
          call eirene_exit_own(1)
        end if

        do i = 1, nalvi
          call json%get_child(palvi,i,pal,found)
          call json%get(pal,'CHRTAL',txt,found)
          chrtal(i) = txt
          deallocate(txt)
          call json%get(pal,'TXTTLR',txt,found)
          TXTTLR(i) = TXT
          deallocate(txt)
          call json%get(pal,'TXTSCR',txt,found)
          TXTSCR(i) = TXT
          deallocate(txt)
          call json%get(pal,'TXTUTR',txt,found)
          TXTUTR(i) = TXT
          deallocate(txt)
          nullify(pal)
        end do
        nullify(palvi)
      end if

      CALL EIRENE_MASAGE
     .  ('*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES')

      if (nadsi > 0) then
        call json%get_child(p,'ADD_VOL_TAL',padsi,found)
        call json%info(padsi,n_children=nsi)
        if (nsi /= nadsi) then
          write (iunout,*)
     .      ' NUMBER OF ADDITIONAL SURFACE AVERAGED TALLIES DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NADSI = ',nadsi
          write (iunout,*) 'NSI =   ',nsi
          call eirene_exit_own(1)
        end if

        do i = 1, nadsi
          call json%get_child(padsi,i,pad,found)
          call json%get(pad,'IADSE',iadse(i),found)
          call json%get(pad,'IADSS',iadss(i),found)
          call json%get(pad,'IADST',iadst(i),found)
          call json%get(pad,'IADSC',iadsc(i),found)
          call json%get(pad,'TXTTLW',txt,found)
          TXTTLW(i,ntlsa) = TXT
          deallocate(txt)
          call json%get(pad,'TXTSPW',txt,found)
          TXTSPW(i,ntlsa) = TXT
          deallocate(txt)
          call json%get(pad,'TXTUNW',txt,found)
          TXTUNW(i,ntlsa) = TXT
          deallocate(txt)
          nullify(pad)
        end do
        nullify(padsi)
      end if

      CALL EIRENE_MASAGE('*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES')

      if (nalsi > 0) then
        call json%get_child(p,'ALG_SRF_TAL',palsi,found)
        call json%info(palsi,n_children=nsi)
        if (nsi /= nalsi) then
          write (iunout,*)
     .      ' NUMBER OF ALGEBRAIC SURFACE AVERAGED TALLIES DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NALSI = ',nalsi
          write (iunout,*) 'NSI =   ',nsi
          call eirene_exit_own(1)
        end if

        do i = 1, nalsi
          call json%get_child(palsi,i,pal,found)
          call json%get(pal,'CHRTLS',txt,found)
          chrtls(i) = txt
          deallocate(txt)
          call json%get(pal,'TXTTLW',txt,found)
          TXTTLW(i,ntlsr) = TXT
          deallocate(txt)
          call json%get(pal,'TXTSPW',txt,found)
          TXTSPW(i,ntlsr) = TXT
          deallocate(txt)
          call json%get(pal,'TXTUNW',txt,found)
          TXTUNW(i,ntlsr) = TXT
          deallocate(txt)
          nullify(pal)
        end do
        nullify(palsi)
      end if

      CALL EIRENE_MASAGE('*** 10F. DATA FOR SPECTRA')
      IF (NADSPC > 0) THEN
        ALLOCATE(ESTIML(NADSPC))
        IF (NSMSTRA > 0) ALLOCATE(SMESTL(NADSPC))
      ELSE
        ALLOCATE(ESTIML(1))
      END IF

      if (nadspc > 0) then
        call json%get_child(p,'SPECTRA',pspeci,found)
        call json%info(pspeci,n_children=nsp)
        if (nsp /= nadspc) then
          write (iunout,*)
     .      ' NUMBER OF SPECTRA DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NADSPC = ',nadspc
          write (iunout,*) 'NSP =   ',nsp
          call eirene_exit_own(1)
        end if

        do j = 1, nadspc
          call json%get_child(pspeci,j,psp,found)
          call json%get(psp,'ISRFCLL',isrfcll,found)
          call json%get(psp,'ISPSRF',ispsrf,found)
          call json%get(psp,'IPTYP',iptyp,found)
          call json%get(psp,'IPSPZ',ipspz,found)
          call json%get(psp,'ISPTYP',isptyp,found)
          call json%get(psp,'NSPS',nsps,found)
          call json%get(psp,'IDIREC',idirec,found)
          call json%get(psp,'ISPCOPT',ispopt,found)

          call json%get(psp,'SPCMN',spcmn,found)
          call json%get(psp,'SPCMX',spcmx,found)
          call json%get(psp,'SPC_SHIFT',spc_shift,found)
          call json%get(psp,'SPCPLT_X',spcplt_x,found)
          call json%get(psp,'SPCPLT_Y',spcplt_y,found)
          call json%get(psp,'SPCPLT_SAME',spcplt_same,found)
          SPCVX = 0._DP
          SPCVY = 0._DP
          SPCVZ = 0._DP
          IF (IDIREC /= 0) THEN
            call json%get(psp,'SPCVX',spcvx,found)
            call json%get(psp,'SPCVY',spcvy,found)
            call json%get(psp,'SPCVZ',spcvz,found)
            VNORM = SQRT(SPCVX**2+SPCVY**2+SPCVZ**2)+EPS60
            SPCVX = SPCVX / VNORM
            SPCVY = SPCVY / VNORM
            SPCVZ = SPCVZ / VNORM
            IF (ISRFCLL == 0) THEN
              WRITE (IUNOUT,*) ' SPECTRUM NUMBER ',J
              WRITE (IUNOUT,*) ' DEFINITION OF LINE OF SIGHT FOR',
     .                         ' SURFACE SPECTRUM IS NOT FORESEEN'
              WRITE (IUNOUT,*) ' SPECIFICATION OF DIRECTION IS',
     .                         ' IGNORED'
              IDIREC = 0
            END IF
          END IF

c   isrcfll=0 : surface-averaged tally
c   isrfcll=1 : volume-averaged tally, integrated over all directions
c   isrfcll=2 : volume-averaged tally, along a specific direction
c  it seems: surface-averaged directional tallies: not yet forseen
c  one may try to score coefficients of orthogonal angular expansion, per energy bin.
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
          IF (ISPOPT==2) THEN
            ESPEC%ISPLDEG = ISPLDEG
          END IF
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
          nullify(psp)
        END DO
        nullify(pspeci)
      ELSE
Cdr  No block 10F found, i.e. no spectrally resolved tallies at all.
cdr  Why do we allocate estiml in input.f and not in eirmod_cestim ?
        IF (.NOT.ALLOCATED(ESTIML)) THEN
          ALLOCATE(ESTIML(1))
        END IF
      END IF

      end subroutine eirene_read_block_10

!*****************************************************************

      subroutine eirene_read_block_11(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: tals, tal, srfs, srf,
     .                             plads, plad, plstds, plstd,
     .                             pltls, pltl, spcs, spc
      logical :: found, foundt, lrpscn, lrps3d, founds
      logical, allocatable :: lhelp(:)
      integer :: i, j, nsf
      integer, allocatable :: ihelp(:)
      real(dp), allocatable :: rhelp(:)
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      external :: eirene_leer, eirene_masage, eirene_masprm,
     .            eirene_exit_own
C
C   READ DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 1100--1199
C
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT')

! block 11a, print output
      call json%get(p,'TRCPLT',trcplt,found)
      call json%get(p,'TRCHST',trchst,found)
      call json%get(p,'TRCNAL',trcnal,found)
      call json%get(p,'TRCMOD',trcmod,found)
      call json%get(p,'TRCSIG',trcsig,found)

      call json%get(p,'TRCGRD',trcgrd,found)
      call json%get(p,'TRCSUR',trcsur,found)
      call json%get(p,'TRCREF',trcref,found)
      call json%get(p,'TRCFLE',trcfle,found)
      call json%get(p,'TRCAMD',trcamd,found)

      call json%get(p,'TRCINT',trcint,found)
      call json%get(p,'TRCLST',trclst,found)
      call json%get(p,'TRCSOU',trcsou,found)
      call json%get(p,'TRCREC',trcrec,found)
      call json%get(p,'TRCTIM',trctim,found)

      call json%get(p,'TRCBLA',trcbla,found)
      call json%get(p,'TRCBLM',trcblm,found)
      call json%get(p,'TRCBLI',trcbli,found)
      call json%get(p,'TRCBLP',trcblp,found)
      call json%get(p,'TRCBLE',trcble,found)

      call json%get(p,'TRCBLPH',trcblph,found)
      call json%get(p,'TRCTAL',trctal,found)
      call json%get(p,'TRCOCT',trcoct,found)
      call json%get(p,'TRCCEN',trccen,found)
      call json%get(p,'TRCRNF',trcrnf,found)
CVK TRACING FOR DEBUGGING, V.Kotov:  not in use in present EIRENE version
cdr  .                  TRCDBG2,TRCDBGE,TRCDBGM,TRCDBGF,TRCDBGL,
cdr  .                  TRCDBGS,TRCDBGG,TRCDBGMPI,TRCDBGC,
CPB  ACTIVATE SPECIES-RESOLVED CPU CONSUMPTION OPTION
      call json%get(p,'TRCHKTIM',trchktim,found)
      call json%get(p,'TRCSCL',trcscl,found)

      call json%get(p,'TRCSRC',lhelp,found)
      trcsrc(0:nstra) = lhelp(1:nstra+1)
      deallocate(lhelp)

! volume averaged tallies
      call json%get(p,'NVOLPR',nvolpr,found)
      call json%get(p,'NSPCPR',nspcpr,found)
      CALL EIRENE_LEER(1)
      IF(NVOLPR > NVLPR) THEN
         WRITE (iunout,*) 'NVOLPR > NVLPR, EXIT.'
         WRITE (iunout,*) ' NVOLPR, NVLPR ', NVOLPR, NVLPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
      WRITE (iunout,*) '       NVOLPR= ',NVOLPR

      call json%get_child(p,'OUT_VOL_TAL',tals,foundt)
      if (foundt) then
        call json%info(tals,n_children=ntl)
        if (ntl /= nvolpr) then
          write (iunout,*)
     .      ' NUMBER OF VOLUME AVERAGED TALLIES TO BE PRINTED DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NVOLPR = ',nvolpr
          write (iunout,*) 'NTL =   ',ntl
          call eirene_exit_own(1)
        end if

        do j = 1, nvolpr
          call json%get_child(tals,j,tal,found)
          call json%get(tal,'NTLV',ntlv,found)
          IF (NTLV.LT.-NTALI.OR.NTLV.GT.NTALV) GOTO 990
          call json%get(tal,'NFLGV',nflagv(j),found)
          call json%get(tal,'NSPEZV',ihelp,found)
          nspezv(j,1:2) = ihelp(1:2)
          deallocate (ihelp)
          call json%get(tal,'NTLVF',ntlvfl(j),found)
cdr  output stream for particular tallies: allow only ntlv>=70
          IF ((NTLVFL(J).NE.0).AND.(NTLVFL(J).LT.70)) GOTO 995
          NPRTLV(J) = NTLV
          nullify(tal)
        end do
      end if
      nullify(tals)

! output of surfaces
      call json%get(p,'NSURPR',nsurpr,found)
      IF(NSURPR > NSRPR) THEN
         WRITE (iunout,*) 'NSURPR > NSRPR, EXIT.'
         WRITE (iunout,*) ' NSURPR, NSRPR ', NSURPR, NSRPR
         CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
      WRITE (iunout,*) '       NSURPR= ',NSURPR

      call json%get_child(p,'SRF_OUTPUT',srfs,founds)
      if (founds) then
        call json%info(srfs,n_children=nsf)
        if (nsf /= nsurpr) then
          write (iunout,*)
     .      ' NUMBER OF SURFACE TALLIES TO BE PRINTED DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NSURPR = ',nsurpr
          write (iunout,*) 'NSF =   ',nsf
          call eirene_exit_own(1)
        end if

        do j = 1, nsurpr
          call json%get_child(srfs,j,srf,found)
          call json%get(srf,'NSRF',nsrf,found)
c  negative surface numbers may be used for non-default standard surfaces
          IF (NSRF.LT.0) NSRF=NLIM+IABS(NSRF)
c  zero as surface number for time horizon
          IF (NSRF.EQ.0.AND.NLIM+NSTSI.LT.NLIMPS) NSRF=NLIM+NSTSI+1

          IF (NSRF.LE.0.OR.NSRF.GT.NLIMPS) GOTO 991
          NPRSRF(J)  =NSRF
          call json%get(srf,'NTLS',nprtls(j),found)
          call json%get(srf,'NFLGS',nflags(j),found)
          call json%get(srf,'NSPEZS',ihelp,found)
          nspezs(j,1:2) = ihelp(1:2)
          deallocate (ihelp)
          call json%get(srf,'NTLSF',ntlsfl(j),found)
          nullify(srf)
        end do
      end if
      nullify(srfs)

c  overrule default switching on/off of volume-averaged tallies

      call json%get(p,'NTLVOUT',ntlvout,foundt)
      if (foundt) then
        call json%get(p,'NUMTAL_V',ihelp,found)
        if (.not.allocated(itlvout)) allocate(itlvout(ntlvout))
        itlvout(1:ntlvout) = ihelp(1:ntlvout)
        deallocate(ihelp)

        do j = 1, ntlvout
          IF ((ITLVOUT(J) /= 0) .AND. (ABS(ITLVOUT(J)) <= NTALV)) THEN
            IF (ITLVOUT(J) < 0) THEN
              LMISTALV(ABS(ITLVOUT(J))) = .TRUE. ! SWITCHED OFF
            ELSE
              LMISTALV(ITLVOUT(J)) = .FALSE. ! SWITCHED ON
            END IF
          END IF
        end do
      end if

c  overrule default switching on/off of surface-averaged tallies

      call json%get(p,'NTLSOUT',ntlsout,foundt)
      if (foundt) then
        call json%get(p,'NUMTAL_S',ihelp,found)
        if (.not.allocated(itlsout)) allocate(itlsout(ntlsout))
        itlsout(1:ntlsout) = ihelp(1:ntlsout)
        deallocate(ihelp)

        do j = 1, ntlsout
          IF ((ITLSOUT(J) /= 0) .AND. (ABS(ITLSOUT(J)) <= NTALV)) THEN
            IF (ITLSOUT(J) < 0) THEN
              LMISTALV(ABS(ITLSOUT)) = .TRUE. ! SWITCHED OFF
            ELSE
              LMISTALV(ITLSOUT(J)) = .FALSE. ! SWITCHED ON
            END IF
          END IF
        end do
      end if

C     read input block 11b (data for plotting)

C  2D GEOMETRY PLOT
      call json%get(p,'PL1ST',pl1st,found)
      call json%get(p,'PL2ND',pl2nd,found)
      call json%get(p,'PL3RD',pl3rd,found)
      call json%get(p,'PLADD',pladd,found)
      call json%get(p,'PLHST',plhst,found)

      call json%get(p,'PLCUT',lhelp,found)
      if (found) then
        plcut(1:3) = lhelp(1:3)
        deallocate(lhelp)
      end if
      call json%get(p,'PLBOX',plbox,found)
      call json%get(p,'PLSTOR',plstor,found)

      call json%get(p,'PLNUMV',plnumv,found)
      call json%get(p,'PLNUMS',plnums,found)
      call json%get(p,'PLARR',plarr,found)
      call json%get(p,'LRPSCUT',lrpscut,found)
      call json%get(p,'PLIDL',plidl,found)

      call json%get(p,'PLVTK',plvtk,found)

      call json%get(p,'NPLINR',nplinr,found)
      call json%get(p,'NPLOTR',nplotr,found)
      call json%get(p,'NPLDLR',npldlr,found)
      call json%get(p,'NPLINP',nplinp,found)
      call json%get(p,'NPLOTP',nplotp,found)
      call json%get(p,'NPLDLP',npldlp,found)
      call json%get(p,'NPLINT',nplint,found)
      call json%get(p,'NPLOTT',nplott,found)
      call json%get(p,'NPLDLT',npldlt,found)

      call eirene_read_block_11_usr(json, p)

C  3D GEOMETRY PLOT

      call json%get_child(p,'3D_ADD_SRF',plads,founds)
      if (founds) then
        do j = 1, 5
          call json%get_child(plads,j,plad,found)
          call json%get(plad,'PL3A',pl3a(j),found)
          call json%get(plad,'TEXTLA',txt,found)
          TEXTLA(J) = TXT
          deallocate(txt)
          call json%get(plad,'IPLTA',iplta(j),found)
          call json%get(plad,'IPLAA',ihelp,found)
          iplaa(j,1:iplta(j)) = ihelp(1:iplta(j))
          deallocate(ihelp)
          call json%get(plad,'IPLEA',ihelp,found)
          iplea(j,1:iplta(j)) = ihelp(1:iplta(j))
          deallocate(ihelp)
          nullify(plad)
        end do
      end if
      nullify(plads)

      call json%get_child(p,'3D_STD_SRF',plstds,founds)
      if (founds) then
        do j = 1, 3
          call json%get_child(plstds,j,plstd,found)
          call json%get(plstd,'PL3S',pl3s(j),found)
          call json%get(plstd,'TEXTLS',txt,found)
          TEXTLS(J) = TXT
          deallocate(txt)
          call json%get(plstd,'IPLTS',iplts(j),found)
          call json%get(plstd,'IPLAS',ihelp,found)
          iplas(j,1:iplts(j)) = ihelp(1:iplts(j))
          deallocate(ihelp)
          call json%get(plstd,'IPLES',ihelp,found)
          iples(j,1:iplts(j)) = ihelp(1:iplts(j))
          deallocate(ihelp)
          nullify(plstd)
        end do
      end if
      nullify(plstds)

      call json%get(p,'CH2MX',ch2mx,found)
      call json%get(p,'CH2MY',ch2my,found)
      call json%get(p,'CH2X0',ch2x0,found)
      call json%get(p,'CH2Y0',ch2y0,found)
      call json%get(p,'CH2Z0',ch2z0,found)

      call json%get(p,'CH3MX',ch3mx,found)
      call json%get(p,'CH3MY',ch3my,found)
      call json%get(p,'CH3MZ',ch3mz,found)
      call json%get(p,'CH3X0',ch3x0,found)
      call json%get(p,'CH3Y0',ch3y0,found)
      call json%get(p,'CH3Z0',ch3z0,found)

      call json%get(p,'ANGLE1',angle1,found)
      call json%get(p,'ANGLE2',angle2,found)
      call json%get(p,'ANGLE3',angle3,found)
C
C  PARTICLE HISTORY PLOTS IN 2D OR 3D GEOMETRY PLOTS
      call json%get(p,'I1TRC',i1trc,found)
      call json%get(p,'I2TRC',i2trc,found)
      call json%get(p,'ISYPLT',ihelp,found)
      if (found) then
        isyplt(1:8) = ihelp(1:8)
        deallocate(ihelp)
      end if
      call json%get(p,'ILINIE',ilinie,found)

! block 11 c, plots of tallies

      call json%get(p,'NVOLPL',nvolpl,found)
      WRITE (iunout,*) '       NVOLPL= ',NVOLPL
      CALL EIRENE_LEER(1)
C
      IF (NVOLPL.GT.NPTAL)
     .  CALL EIRENE_MASPRM('NPTAL',5,NPTAL,'NVOLPT',6,NVOLPL,IERROR)
C
      if (nvolpl > 0) then
        call json%get(p,'PLTSRC',lhelp,found)
        pltsrc=.false.
        if (found) then
          pltsrc(0:size(lhelp)-1) = lhelp
          deallocate(lhelp)
        end if

        if (lrpscut) then
          call json%get(p,'CUTPLANE',rhelp,found)
          cutplane(1:4) = rhelp(1:4)
          deallocate(rhelp)
!  TEST IF CUTPLANE IS AXIS-PARALLEL
          IF ((SUM(ABS(CUTPLANE(2:4)))-1._DP > EPS10) .OR.
     .        (CUTPLANE(2)*CUTPLANE(3) > 0._DP) .OR.
     .        (CUTPLANE(2)*CUTPLANE(4) > 0._DP) .OR.
     .        (CUTPLANE(3)*CUTPLANE(4) > 0._DP)) THEN
            WRITE (iunout,*) ' PLANE IS NOT PARALLEL TO ONE AXIS'
            WRITE (iunout,*) ' CUTTING OF TETRAHEDRA ABANDONED'
            WRITE (iunout,*) ' LRPSCUT SET TO .FALSE.'
            LRPSCUT = .FALSE.
          END IF
        end if

        IPLANE=1
        LRAPS3D=.FALSE.
        LR3DCON=.FALSE.
        RAPSDEL=-HUGE(1._DP)
C

        call json%get_child(p,'PLOT_TALLIES',pltls,found)
        call json%info(pltls,n_children=ntl)
        if (ntl /= nvolpl) then
          write (iunout,*)
     .      ' NUMBER OF TALLIES TO BE PLOTTED DOES',
     .      ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NVOLPR = ',nvolpr
          write (iunout,*) 'NTL =   ',ntl
          call eirene_exit_own(1)
        end if

        do j = 1, nvolpr
          call json%get_child(pltls,j,pltl,found)

          call json%get(pltl,'NSP',nsptal(j),found)

          call json%get(pltl,'PLTL2D',pltl2d(j),found)
          call json%get(pltl,'PLTL3D',pltl3d(j),found)
          call json%get(pltl,'PLTLLG',pltllg(j),found)
          call json%get(pltl,'PLTLER',pltler(j),found)

          call json%get(pltl,'TALZMI',talzmi(j),found)
          call json%get(pltl,'TALZMA',talzma(j),found)
          call json%get(pltl,'TALXMI',talxmi(j),found)
          call json%get(pltl,'TALXMA',talxma(j),found)
          call json%get(pltl,'TALYMI',talymi(j),found)
          call json%get(pltl,'TALYMA',talyma(j),found)

          if (pltl2d(j)) then
            call json%get(pltl,'LHIST2',lhist2(j),found)
            call json%get(pltl,'LSMOT2',lsmot2(j),found)

            call json%get_child(pltl,'SPECIES',spcs,found)
            call json%info(spcs,n_children=nsp)
            if (nsp /= nsptal(j)) then
              write (iunout,*)
     .             ' NUMBER OF SPECIES IN 2D PLOT DOES',
     .             ' NOT MATCH THE NUMBER FOUND IN FILE '
              write (iunout,*) 'PLOT NO = ',j
              write (iunout,*) 'NSPTAL = ',nsptal(j)
              write (iunout,*) 'NSP =   ',nsp
              call eirene_exit_own(1)
            end if
            do i = 1, nsptal(j)
              call json%get_child(spcs,i,spc,found)
              call json%get(spc,'ISPTAL',isptal(j,i),found)
              call json%get(spc,'NTL',ntl,found)
              call json%get(spc,'NPLIN2',nplin2(j,i),found)
              call json%get(spc,'NPLOT2',nplot2(j,i),found)
              call json%get(spc,'NPLDL2',npldl2(j,i),found)
              nullify(spc)
              IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
              NPTALI(J,I)=NTL
              NPLIN2(J,I) = MAX(NPLIN2(J,I),1)
              NPLDL2(J,I) = MAX(NPLDL2(J,I),1)
            end do
            nullify(spcs)
          end if

          if (pltl3d(j)) then
            call json%get(pltl,'LHIST3',lhist3(j),found)
            call json%get(pltl,'LCNTR3',lcntr3(j),found)
            call json%get(pltl,'LSMOT3',lsmot3(j),found)
            call json%get(pltl,'LRAPS3',lraps3(j),found)
            call json%get(pltl,'LVECT3',lvect3(j),found)
            call json%get(pltl,'LRPVC3',lrpvc3(j),found)
            call json%get(pltl,'LRPS3D',lrps3d)
            call json%get(pltl,'LRPSCN',lrpscn)

            call json%get(pltl,'LPRAD3',lprad3(j),found)
            call json%get(pltl,'LPPOL3',lppol3(j),found)
            call json%get(pltl,'LPTOR3',lptor3(j),found)

            call json%get_child(pltl,'SPECIES',spcs,found)
            call json%info(spcs,n_children=nsp)
            if (nsp /= nsptal(j)) then
              write (iunout,*)
     .             ' NUMBER OF SPECIES IN 3D PLOT DOES',
     .             ' NOT MATCH THE NUMBER FOUND IN FILE '
              write (iunout,*) 'PLOT NO = ',j
              write (iunout,*) 'NSPTAL = ',nsptal(j)
              write (iunout,*) 'NSP =   ',nsp
              call eirene_exit_own(1)
            end if
            do i = 1, nsptal(j)
              call json%get_child(spcs,i,spc,found)
              call json%get(spc,'ISPTAL',isptal(j,i),found)
              call json%get(spc,'NTL',ntl,found)
              call json%get(spc,'IPROJ3',iproj3(j,i),found)
              call json%get(spc,'NPLI13',npli13(j,i),found)
              call json%get(spc,'NPLO13',nplo13(j,i),found)
              call json%get(spc,'NPLI23',npli23(j,i),found)
              call json%get(spc,'NPLO23',nplo23(j,i),found)
              call json%add(spc,'IPLN',ipln)
              nullify(spc)
              IF (NTL.LT.-NTALI.OR.NTL.GT.NTALV.OR.NTL.EQ.0) GOTO 990
              NPTALI(J,I)=NTL
            end do
            nullify(spcs)

            call json%get(pltl,'TALW1',talw1(j),found)
            call json%get(pltl,'TALW2',talw2(j),found)
            call json%get(pltl,'FCABS1',fcabs1(j),found)
            call json%get(pltl,'FCABS2',fcabs2(j),found)
            call json%get(pltl,'RPSDL',rpsdl,found)
            IF (FCABS1(J).LE.0.0D0) FCABS1(J)=1.0D0
            IF (FCABS2(J).LE.0.0D0) FCABS2(J)=1.0D0
            LRAPS3D=LRAPS3D.OR.LRPS3D
            LR3DCON=LR3DCON.OR.LRPSCN
            IPLANE=MAX(IPLANE,IPLN)
            RAPSDEL=MAX(RAPSDEL,RPSDL)

            nullify(pltl)
          end if
        end do

        IF (NLTRA) RAPSDEL=RAPSDEL*DEGRAD

      end if

      return
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
  995 CONTINUE
      WRITE (iunout,*) 'ERROR IN INPUT: ',
     .            'WRONG UNIT NUMBER FOR OUTPUT OF TALLY SPECIFIED'
      WRITE (iunout,*) 'NTLV,NTLVF',NTLV,NTLVF
      CALL EIRENE_EXIT_OWN(1)

      end subroutine eirene_read_block_11

!*****************************************************************

      subroutine eirene_read_block_12(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: plines, pline, pchords, pchord,
     .    pcomps, pcomp, pconts, pcont, prats, prat
      logical :: found
      integer :: nlines, nchrds, iline, num_compo, ncomps, jcomp,
     .           num_contrib, irc, nconts, kcontr, nrats, irt, istchr
      integer, allocatable :: ihelp(:)
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt, key
!cym/cpg end
      external :: eirene_setup_default_emissivity,
     .            eirene_leer, eirene_masage, eirene_masprm,
     .            eirene_exit_own
C
C  READ DATA FOR DIAGNOSTIC MODULE
C
      CALL EIRENE_MASAGE('*** 12. DATA FOR DIAGNOSTIC MODULE')

c June 18: new: more general option for definition of emission profiles
cdr        in ADDV tally. Can e.g. be used in line of sight NCHTAL=2 option,
cdr        or for general additional output.
cdr  Based on population coefficients, QSS ratios, etc already read in block 4.

cdr  read further atomic/molecular data: population coefficients, QSS ratios, etc
cdr       needed for setting up volumetric line emissivity profiles
cdr       as further add. tally ADDV (additional to those already defined in block 10a)
cdr       The ADDV tallies may then be used for line of sight integration, along
cdr       the CHORDS defined further below.
cdr  Distinct from the other addv tallies from block 10a these
cdr  further ADDV tallies (beyond NADVI) have their text
cdr  (tally label, units, species) assigned in subr. emissivity.f.

cdr  check if such emissivity tallies are defined:
cdr  Search for 'DEFINE LINES' in next card
      CALL EIRENE_MASAGE
     .  ('*** 12A. DATA FOR EMISSIVITY PROFILES (ADDV)  ')
      NLEMIS = .FALSE.
      LDEF_LINES = .FALSE.

      NCHORI = 0
      NCHENI = 0

      if (.not.associated(p)) return

!  read definition of emissivity lines
      call json%get(p,'KEYWORD',key,found)
      if (key == 'DEFINE_LINES') then
        IADV = NADVI
        NLEMIS = .TRUE.

        call json%get(p,'NUM_LINES',num_lines,found)
        call json%get(p,'MOD_ADDV',mod_addv,found)
      end if

!   Option 1: Non-zero number of lines are defined
      if ((key == 'DEFINE_LINES').and.(num_lines > 0)) then
cdr read volumetric emission profile data
        LDEF_LINES = .TRUE.

        ALLOCATE (EMIS_LINES(NUM_LINES))
        EMIS_LINES%LINE_NAME = REPEAT(' ',80)
        EMIS_LINES%NUM_COMPO = 0

        call json%get_child(p,'LINES',plines,found)
        call json%info(plines,n_children=nlines)
        if (nlines /= num_lines) then
          write (iunout,*)
     .            ' NUMBER OF EMISSIVITY LINES DOES',
     .            ' NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NUM_LINES = ',num_lines
          write (iunout,*) 'NLINE =   ',nlines
          call eirene_exit_own(1)
        end if

        WRITE (IUNOUT,*) 'FURTHER EMISSIVITY PROFILES SET:'
        do iline = 1, num_lines
          call json%get_child(plines,iline,pline,found)
cdr  further below, a particular emission profile is identified
cdr  (e.g. for a chord ichori) either by its name  (and ch_line%...)
cdr  or, if that fails, by its energy (and EMIN1 flag)
          call json%get(pline,'LINE_NAME',txt,found)
          EMIS_LINES(ILINE)%LINE_NAME = txt
          deallocate(txt)
          write (iunout,*) trim(EMIS_LINES(ILINE)%LINE_NAME)

          call json%get(pline,'NUM_COMPO',num_compo,found)
          call json%get(pline,'EINSTEIN',
     .                  EMIS_LINES(ILINE)%EINSTEIN,found)
          call json%get(pline,'TRANS_EN',
     .                  EMIS_LINES(ILINE)%TRANS_EN,found)
C
          EMIS_LINES(ILINE)%NUM_COMPO = NUM_COMPO

c  Deal with ADDV storage for sum over components.
c  The storage on ADDV for individual components is done below (JCOMP loop).
cdr       IADV = IADV + 1
c  if mod_addv=0: reset storage needs for each line
c                 back to nadvi+1..nadvi+num_compo+1,
c                 i.e. addv tallies are only saved for one line at a time.
cdr         IF (MOD_ADDV == 0) IADV = NADVI + 1
          IF (MOD_ADDV == 0) IADV = NADVI

          EMIS_LINES(ILINE)%IADV_TOTAL = IADV + NUM_COMPO+1

          IF (NUM_COMPO > 0) THEN

            ALLOCATE (EMIS_LINES(ILINE)%COMPO(NUM_COMPO))
            call json%get_child(pline,'COMPONENTS',pcomps,found)
            call json%info(pcomps,n_children=ncomps)
            if (ncomps /= num_compo) then
              write (iunout,*)
     .            ' NUMBER OF COMPONENTS OF EMISSIVITY LINE DOES',
     .            ' NOT MATCH THE NUMBER FOUND IN FILE '
              write (iunout,*) 'NUM_COMPO = ',num_compo
              write (iunout,*) 'NCOMPS =   ',ncomps
              call eirene_exit_own(1)
            end if

            do jcomp = 1, num_compo
              call json%get_child(pcomps,jcomp,pcomp,found)
              call json%get(pcomp,'COMPO_NAME',txt,found)
              EMIS_LINES(ILINE)%COMPO(JCOMP)%COMPO_NAME = txt
              deallocate(txt)
              call json%get(pcomp,'NUM_CONTRIB',num_contrib,found)
              call json%get(pcomp,'COMPO_IRC',irc,found)
              EMIS_LINES(ILINE)%COMPO(JCOMP)%NUM_CONTRIB = NUM_CONTRIB
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IRC = IRC
              ALLOCATE
     .          (EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(NUM_CONTRIB))
cdr  for each new component: store emission profile on addv(iadv)
              IADV = IADV + 1
              EMIS_LINES(ILINE)%COMPO(JCOMP)%IADV = IADV

cdr  now we dwell on the contributions:
cdr  e.g. different isotopes,... but same rates, same population factors in each component
              call json%get_child(pcomp,'CONTRIBUTIONS',pconts,found)
              call json%info(pconts,n_children=nconts)
              if (nconts /= num_contrib) then
                write (iunout,*)
     .            ' NUMBER OF CONTRIBUTIONS TO COMPONENT OF EMIS ',
     .            ' LINE DOES NOT MATCH THE NUMBER FOUND IN FILE '
                write (iunout,*) 'NUM_CONTRIB = ',num_contrib
                write (iunout,*) 'NCONTS =   ',nconts
                call eirene_exit_own(1)
              end if

              do kcontr = 1, num_contrib
                CNT%ISP = -1
                CNT%ITP = -1
                CNT%IRATIO = 0
                CNT%IRC_RAT = 0
                CNT%ISP_RAT = -1
                CNT%ITP_RAT = -1
                call json%get_child(pconts,kcontr,pcont,found)
                call json%get(pcont,'ISP',cnt%isp,found)
                call json%get(pcont,'ITP',cnt%itp,found)
                call json%get(pcont,'IRATIO',cnt%iratio,found)

cdr  read QSS ratio between two densities, e.g.:  H2+/H2, if nfoli(H2+)=-1.
cdr  If density n_B of parent state (e.g. H2+) for upper level is not amongst the densities
cdr  known in this run,
cdr  but is in an (QSS) equilibrium (nfol<0) with such a density n_A (e.g. H2) instead.
                IF (CNT%IRATIO > 0) THEN
cdr  read QSS density ratio n_B/n_A(Te,ne). Then the
cdr  upper state population is n_B *pop_B(upper)= n_A * ratio * pop_B(upper)
                  call json%get_child(pcont,'RATIOS',prats,found)
                  call json%info(prats,n_children=nrats)
                  if (nrats /= cnt%iratio) then
                    write (iunout,*)
     .                    ' NUMBER OF RATIOS IN CONTRIBUTION ',
     .                    ' DOES NOT MATCH THE NUMBER FOUND IN FILE '
                    write (iunout,*) 'IRATIO = ',cnt%iratio
                    write (iunout,*) 'NRATS =   ',nrats
                    call eirene_exit_own(1)
                  end if
                  do irt = 1, cnt%iratio
                    call json%get_child(prats,irt,prat,found)
                    call json%get(prat,'IRC_RAT',
     .                            cnt%irc_rat(irt),found)

cdr   Read a second density ratio.
cdr   It may turn out that, after reading the first density ratio,
cdr   that now the new parent density n_A is still not amongst the density
cdr   known to eirene in this run. Then read a second density ratio.
cdr   Example:  n_B=n_H3+, ratio1=nH3+/nH2+. I.e. n_A =nH2+.
cdr   See routine emissivity.f for further explanations.
                    if (irt == 2) then
                      call json%get(prat,'ISP_RAT',ihelp,found)
                      cnt%isp_rat(1:2) = ihelp(1:2)
                      deallocate(ihelp)
                      call json%get(prat,'ITP_RAT',ihelp,found)
                      cnt%itp_rat(1:2) = ihelp(1:2)
                      deallocate(ihelp)
                    END IF
                    nullify(prat)
                  end do      ! irt
                  nullify(prats)
                END IF

                EMIS_LINES(ILINE)%COMPO(JCOMP)%CONTRIB(KCONTR) = CNT
                nullify(pcont)
              end do          ! kcontr

              nullify(pconts)
              nullify(pcomp)
            end do            ! jcomp

c  increase counter iadv, for next line, because sum over comp.
c                         is stored on num_compo+1
            IADV=IADV+1
            nullify(pcomps)
          end if

          nullify(pline)
        end do                ! iline
        nullify(plines)

!   Option 2: Requesting default lines, there are two ways:
!       - DEFAULT_LINES
!       - DEFINE_LINES && NUM_LINES=0  (backward compatibility)
      elseif ((key == 'DEFAULT_LINES')
     .  .or.(key == 'DEFINE_LINES').and.(num_lines <= 0)) then

        nlemis = .true.

      end if
      deallocate(key)

! no definition of emissivity lines was read in
! define OLD default emissivity model for chords, for backward compatibility.
cdr  allocate storage and fill structure EMIS_LINES
cdr  such that old default options are recovered.
cdr This is exclusive: as soon as at least one line emission profile is
cdr read from block "12.0", no default emissivities are set at all.
      IF (NLEMIS.AND..NOT.ALLOCATED(EMIS_LINES))
     .   CALL EIRENE_SETUP_DEFAULT_EMISSIVITY

c  June 18: new:  end of new code, further modifications below, for NCHTAL=2 option
      CALL EIRENE_MASAGE
     .  ('*** 12B. DATA FOR LINES OF SIGHT              ')

      call json%get(p,'NCHORI',nchori,found)
      call json%get(p,'NCHENI',ncheni,found)
      WRITE (iunout,*) '       NCHORI= ',NCHORI
      NCHOR = NCHORI
      NCHEN = NCHENI
      CALL EIRENE_LEER(1)
      WRITE (iunout,'(1x,a,2i5)') 'NCHORI,NCHENI= ',NCHORI,NCHENI
      CALL EIRENE_LEER(1)
      IF (IABS(NCHENI).GT.NCHEN)
     .    CALL EIRENE_MASPRM
     .         ('NCHEN',5,NCHEN,'IABS(NCHENI)',12,IABS(NCHENI),
     .                 IERROR)

      if (nchori > 0) then
        CALL EIRENE_ALLOC_COMSIG
        CALL EIRENE_INIT_COMSIG
        call json%get_child(p,'CHORDS',pchords,found)
        call json%info(pchords,n_children=nchrds)
        if (nchrds /= nchori) then
          write (iunout,*)
     .      ' NUMBER OF CHORDS ',
     .      ' DOES NOT MATCH THE NUMBER FOUND IN FILE '
          write (iunout,*) 'NCHORI = ',nchori
          write (iunout,*) 'NCHRDS = ',nchrds
          call eirene_exit_own(1)
        end if
        do ichori = 1, nchori
          call json%get_child(pchords,ichori,pchord,found)
          call json%get(pchord,'TXTSIG',txt,found)
          txtsig(ichori) = txt
          deallocate(txt)

          call json%get(pchord,'NCHTAL',nchtal(ichori),found)
          call json%get(pchord,'NSPSCL',nspscl(ichori),found)
          call json%get(pchord,'NSPNEW',nspnew(ichori),found)
          call json%get(pchord,'ISTCHR',istchr,found)
          nlstchr(ichori) = istchr > 0 ! automatically add directional
                                       ! cell-based spectra,
                                       ! along line of sight

          IF (NCHTAL(ICHORI) == 2) THEN
cdr  old default: neither "define_lines" nor "default_lines" found.
cdr  still:  chords with NCHTAL=2 present. Set default now, latest.
            IF (.NOT.ALLOCATED(EMIS_LINES))
     .        CALL EIRENE_SETUP_DEFAULT_EMISSIVITY
cdr  for nchtal=2: one optional extra input card may be read:
cdr                search for 'USE_LINE'
cdr                and fill CH_LINE_NAME(ICHORI) with the name of that line.
c     dr  Alternatively the energy parameters EMIN1 may be used.
            call json%get(pchord,'USE_LINE',txt,found)
            ch_line_name(ichori) = ''
            if (found) ch_line_name(ichori) = txt
            deallocate(txt)
          END IF

          call json%get(pchord,'NSPSTR',nspstr(ichori),found)
          call json%get(pchord,'NSPSPZ',nspspz(ichori),found)
! here should come: NSPTP(..), TYPE
          call json%get(pchord,'NSPINI',nspini(ichori),found)
          call json%get(pchord,'NSPEND',nspend(ichori),found)
          call json%get(pchord,'NSPBLC',nspblc(ichori),found)
          call json%get(pchord,'NSPADD',nspadd(ichori),found)

          call json%get(pchord,'EMIN1',emin1(ichori),found)
          call json%get(pchord,'EMAX1',emax1(ichori),found)
          call json%get(pchord,'ESHIFT',eshift(ichori),found)

          call json%get(pchord,'IPIVOT',ipivot(ichori),found)
          call json%get(pchord,'XPIVOT',xpivot(ichori),found)
          call json%get(pchord,'YPIVOT',ypivot(ichori),found)
          call json%get(pchord,'ZPIVOT',zpivot(ichori),found)

          call json%get(pchord,'ICHORD',ichord(ichori),found)
          call json%get(pchord,'XCHORD',xchord(ichori),found)
          call json%get(pchord,'YCHORD',ychord(ichori),found)
          call json%get(pchord,'ZCHORD',zchord(ichori),found)

          nullify(pchord)
        end do                  ! ichori
        nullify(pchords)

        call json%get(p,'PLCHOR',plchor,found)
        call json%get(p,'PLSPEC',plspec,found)
        call json%get(p,'PRSPEC',prspec,found)
        call json%get(p,'PLARGL',plargl,found)
        call json%get(p,'PRARGL',prargl,found)

      end if                    ! nchori > 0

      end subroutine eirene_read_block_12

!*****************************************************************

      subroutine eirene_read_block_13(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: psnaps, psnap
      logical :: found
      integer :: j, nsnaps
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      external :: eirene_leer, eirene_masage, eirene_exit_own
C
C  READ DATA FOR TIME-DEPENDENT AND NONLINEAR MODE  1300--1399
C
      CALL EIRENE_MASAGE
     .  ('*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION')
C
      NINITL_READ = 0
      nprnli_in = 0
      nprmul = 0

      if (associated(p)) then
        call json%get(p,'NPRNLI',nprnli_in,found)
        call json%get(p,'DEFAULT_TIME_HORIZON',ldef_time_horizon,found)
        if (nprnli_in > 0) then
          call json%get(p,'NINITL_READ',ninitl_read,found)
          call json%get(p,'NPRMUL',nprmul,found)
        else
          ninitl_read = 0
          nprmul = 0
        end if
      end if

      nprnli = nprnli_in
      IF (NPRMUL > 1) NPRNLI=NPRNLI * NPRMUL

      CALL EIRENE_LEER(1)

C     DATA FOR DEFAULT TIME HORIZON
      IF (LDEF_TIME_HORIZON) THEN
!PB     NTIME0=0
        NTIME0=MAX(0,NTIME0)
!PB     NTIME=1
        NTIME=MAX(1,NTIME)
        DTIMV=1.D0
        TIME0=0.D0
        NPTST=0
        NTMSTP=1
        NSNVI=0
        IREAD=0
      END IF

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

      ELSE                      ! read data for time horizon

!PB     if (.not.associated(p)) return
        if (ldef_time_horizon) return

        call json%get(p,'NPTST',nptst,found)
        call json%get(p,'NTMSTP',ntmstp,found)

        call json%get(p,'DTIMV',dtimv,found)
        call json%get(p,'TIME0',time0,found)

        CALL EIRENE_MASAGE('*** 13A. DATA FOR SNAPSHOT TALLIES')
        call json%get(p,'NSNVI',nsnvi,found)
        WRITE (iunout,*) '        NSNVI= ',NSNVI
        CALL EIRENE_LEER(1)

        IF (NSNVI.GT.0) THEN

          call json%get_child(p,'SNAPSHOTS',psnaps,found)
          call json%info(psnaps,n_children=nsnaps)
          if (nsnaps /= nsnvi) then
            write (iunout,*)
     .        ' NUMBER OF SNAPSHOT TALLIES ',
     .        ' DOES NOT MATCH THE NUMBER FOUND IN FILE '
            write (iunout,*) 'NSNVI = ',nsnvi
            write (iunout,*) 'NSNAPS = ',nsnaps
            call eirene_exit_own(1)
          end if
          do j = 1, nsnvi
            call json%get_child(psnaps,j,psnap,found)

            call json%get(psnap,'ISNVE',isnve(j),found)
            call json%get(psnap,'ISNVS',isnvs(j),found)
            call json%get(psnap,'ISNVT',isnvt(j),found)
            call json%get(psnap,'ISNRC',isnrc(j),found)

            call json%get(psnap,'TXTTLT',txt,found)
            txttlt(j) = txt
            deallocate(txt)
            call json%get(psnap,'TXTSCT',txt,found)
            txtsct(j) = txt
            deallocate(txt)
            call json%get(psnap,'TXTUTT',txt,found)
            txtutt(j) = txt
            deallocate(txt)
            nullify(psnap)
          end do
          nullify(psnaps)
        ENDIF

        IF (NTIME.LE.0) THEN
          WRITE (iunout,*)
     .      'ERROR IN INPUT: TIME DEP. MODE BUT NTIME.LE.0'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

      ENDIF

      end subroutine eirene_read_block_13

#endif
      END SUBROUTINE EIRENE_READ_JSON


      SUBROUTINE EIRENE_READ_BLK14_JSON(IERROR)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CINIT
      USE EIRMOD_COMPRT
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_JSON
#ifndef NO_JSON
      USE EIRMOD_INFCOP, ONLY: EIRENE_IF0COP
      use json_module           !IGNORE
!cym/cpg
!     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck
!cym/cpg end
#endif
      IMPLICIT NONE

      integer, intent(inout) :: ierror
#ifndef NO_JSON
      type(json_core) :: json
      type(json_value), pointer :: p
      type(json_value), pointer :: padds, padd
      logical :: found, founda, lshort
      integer :: j, nadds, ncopie
!cym/cpg CK-> JSON_CK, in line with not renaming types after use json_modules
      character(kind=json_CK,len=:), allocatable :: txt
!cym/cpg end
      integer :: js
      external :: eirene_leer, eirene_masage, eirene_masprm,
     .            eirene_exit_own

      js = itree_num(14)
      json = jtrees(js)
      p => blks(14)%p

      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"')
      IF (NMODE.EQ.0) THEN
C  STAND ALONE RUN, READ BLOCK *** 14 HERE
        WRITE (iunout,*) '        SUBR. INFCOP NOT CALLED.'

        call json%get(p,'NAINI',naini,found)
!PB     call json%get(p,'NCOPII',ncopii,found)
        call json%get(p,'NCOPIE',ncopie,found)

        NCPVI=NCOPIE
        WRITE (iunout,*) '        NAINI, NCPVI = ',NAINI,NCPVI
        IF (NAINI.GT.NAIN) THEN
          CALL EIRENE_MASPRM('NAIN',4,NAIN,'NAINI',5,NAINI,IERROR)
          WRITE (iunout,*) IERROR,
     .      ' INPUT OR PARAMETER ERRORS DETECTED'
          WRITE (iunout,*)
     .      ' SEE THE ERROR MESSAGES LISTED ABOVE AND CORRECT'
          WRITE (iunout,*) ' THE ERRORS BEFORE RE-EXECUTION'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        IF (NCPVI.GT.NCPV) THEN
          CALL EIRENE_MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF

        call json%get_child(p,'ADD_IN_TAL',padds,founda)
        if (founda) then
          call json%info(padds,n_children=nadds)
          if (nadds /= naini) then
            write (iunout,*)
     .        ' NUMBER OF ADDITIONAL INPUT TALLIES ',
     .        ' DOES NOT MATCH THE NUMBER FOUND IN FILE '
            write (iunout,*) 'NAINI = ',naini
            write (iunout,*) 'NADDS = ',nadds
            call eirene_exit_own(1)
          end if
          do j = 1, naini
            call json%get_child(padds,j,padd,found)

            call json%get(padd,'NAINS',nains(j),found)
            call json%get(padd,'NAINT',naint(j),found)

            call json%get(padd,'TXTPLS',txt,found)
            txtpls(j,ntaln) = txt
            deallocate(txt)
            call json%get(padd,'TXTPSP',txt,found)
            txtpsp(j,ntaln) = txt
            deallocate(txt)
            call json%get(padd,'TXTPUN',txt,found)
            txtpun(j,ntaln) = txt
            deallocate(txt)
            nullify(padd)
          end do
        end if
        nullify(padds)

      ELSEIF (NMODE.NE.0) THEN
C  COUPLED RUN, READ BLOCK *** 14 IN INTERFACING ROUTINE INFCOP (ENTRY IF0COP)
        NCPVI=0
!pb     NCOPII=0
        NCOPIE=0
        NAINI=0
C  READ BLOCK 14 AND GEOMETRY FROM EXTERNAL DATABASE (FT30),
C  also set NAINI, NCOPII, NCOPIE there
        LSHORT=.FALSE.
        CALL EIRENE_IF0COP(.FALSE.,LSHORT)
        WRITE (iunout,*) 'SUBR. INFCOP IS CALLED.'
!pb     write (iunout,*) 'NAINI, NCOPII, NCOPIE, NCPVI ',
!pb  .                    NAINI, NCOPII, NCOPIE, NCPVI
        write (iunout,*) 'NAINI, NCOPIE, NCPVI ',
     .                    NAINI, NCOPIE, NCPVI
        CALL EIRENE_LEER(1)
      END IF
#endif
      RETURN
      END SUBROUTINE EIRENE_READ_BLK14_JSON

!*****************************************************************

      subroutine eirene_read_mpi_strategy_json
      USE EIRMOD_JSON
#ifndef NO_JSON
      USE EIRMOD_COMUSR, ONLY: NPRLL
      USE EIRMOD_CPES, ONLY: NPRS,
     >    STRATEGY_UNDEFINED, STRATEGY_EMBARRASS,
     >    STRATEGY_ORIGINAL, STRATEGY_APCAS,
     >    STRATEGY_BALANCED, STRATEGY_IMPROVED,
     >    INPUT_DISTRIBUTION_STRATEGY
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      type(json_core) :: json
      type(json_value), pointer :: p
      character(80) :: strategy
      character(kind=CK,len=:), allocatable :: txt
      logical :: found
      integer :: js
      external :: eirene_leer, eirene_masage

      js = itree_num(15)
      json = jtrees(js)
      p => blks(15)%p

      strategy = repeat(' ',80)

      if (associated(p)) then
        call json%get(p,'STRATEGY',txt,found)
        if (found) then
          CALL EIRENE_MASAGE('*** INFORMATION_FOR_MPI')
          strategy(1:len_trim(txt)) = trim(txt)
        endif
        deallocate(txt)
      end if

      select case (strategy)
      case ('ORIGINAL')
        NPRLL = 1
        input_distribution_strategy = STRATEGY_ORIGINAL
        CALL EIRENE_MASAGE
     .           ('APPLYING "ORIGINAL" PARALLELIZATION STRATEGY')
      case ('APCAS')
        NPRLL = 2
        input_distribution_strategy = STRATEGY_APCAS
        CALL EIRENE_MASAGE
     .           ('APPLYING "APCAS" PARALLELIZATION STRATEGY')
      case ('BALANCED')
        NPRLL = 3
        input_distribution_strategy = STRATEGY_BALANCED
        CALL EIRENE_MASAGE
     .           ('APPLYING "BALANCED" PARALLELIZATION STRATEGY')
      case ('IMPROVED')
        NPRLL = 4
        input_distribution_strategy = STRATEGY_IMPROVED
        CALL EIRENE_MASAGE
     .           ('APPLYING "IMPROVED" PARALLELIZATION STRATEGY')
      case ('EMBARRASS')
        NPRLL = 0
        input_distribution_strategy = STRATEGY_EMBARRASS
        CALL EIRENE_MASAGE
     .           ('APPLYING "EMBARRASSINGLY PARALLEL" STRATEGY')
      case ('AUTOMATIC')
        CALL EIRENE_MASAGE
     .           ('APPLYING AUTOMATIC PARALLELIZATION STRATEGY')
#if MPI_VERSION < 3
        NPRLL = 1
        input_distribution_strategy = STRATEGY_ORIGINAL
        CALL EIRENE_MASAGE
     .           ('CURRENT AUTOMATIC ASSIGNMENT IS "ORIGINAL"')
#else
        NPRLL = 3
        input_distribution_strategy = STRATEGY_BALANCED
        CALL EIRENE_MASAGE
     .           ('CURRENT AUTOMATIC ASSIGNMENT IS "BALANCED"')
#endif
      case default
        IF (NPRS > 1 .and.
     &     input_distribution_strategy.eq.STRATEGY_UNDEFINED) THEN
          CALL EIRENE_MASAGE
     .       ('NO MPI PARALLELIZATION STRATEGY PROVIDED')
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
          CALL EIRENE_LEER(1)
        END IF
      end select
#endif
      return
      end subroutine eirene_read_mpi_strategy_json
