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
cdr             same thing: NCPVI, NCPV  (and eliminate old parameters NCOP, NCOPI)
cdr  July 17 :  lmulpl:  automatic options for multiple ion temperatures,
cdr                      multiple ion velocities in case of BGK nonlinear collisions
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
      SUBROUTINE EIRENE_FIND_PARAM_JSON
#ifndef NO_JSON
C
C   SET DIMENSIONS (PARAMETERS) FOR ALLOCATABLE ARRAYS (TALLIES, GRIDS, ETC..) AND SET SOME FURTHER DEFAULT VALUES
C
      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: N1ST, N2ND, N3RD, NADD, NADS, NADV,
     .                          NAIN, NALS, NALV, NATM, NCLV,
     .                          NCOORD, NCOUP_INPUT, NCPV, NCV,
     .                          NGEOM_USR, NGSTAL, NION, NKNOT, NLIM,
     .                          NLIMPS, NMOL, NOPTIM, NOPTM1, NPHOT,
     .                          NPLG, NPLS, NPPART, NPRNL, NPTRGT, NRAD,
     .                          NREAC, NRTAL, NSD, NSDW, NSMSTRA, NSNV,
     .                          NSRFS, NSTEP, NSTORAM, NSTRA, NSTS,
     .                          NTETRA, NTOR, NTRI, NCHEN, NCHOR, NPLT,
     .                          NSRPR, NVLPR, NHD6, NPLSTI, NPLSV,
     .                          NREAC_LINES, IFOFF, NBMAX, NGITT,
     .                          NKNOTS, NTRIS
      USE EIRMOD_CTRCEI, ONLY: TRCAMD, TRCINT,
     .                         EIRENE_ALLOC_CTRCEI
      USE EIRMOD_CINIT, ONLY: DBFNAME, DBHANDLE, NDBNAMES,
     .                        INDPRO2_SAVE, LDBREAD,
     .                        EIRENE_INIT_CINIT
      USE EIRMOD_COMUSR, ONLY: NMODE, NTIME,
     .                         NPHOTI, NATMI, NMOLI, NIONI, NPLSI,
     .                         NMASSA, NCHARA, NMASSM, NMASSI, NCHARI,
     .                         NCHRGI, LKINDI, COMUSR_FIRST_PASS
      USE EIRMOD_COMSOU, ONLY: NSTRAI
      USE EIRMOD_COMNNL, ONLY: NPRNLI
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CLOGAU, ONLY: NLWRMSH, NLSPCSCL, EIRENE_ALLOC_CLOGAU,
     .                         NLSPCSCL_ON, NLSOLEDGE
      USE EIRMOD_JSON, ONLY: jtrees, blks, itree_num,
     .                       ldef_time_horizon,nlfem_in, nlplg_in,
     .                       nlpol_in, np2nd_in, nr1st_in, nt3rd_in,
     .                       noptim_in, nrtal_in, nsmstra_in, nstrai_in,
     .                       dbfname_in,
     .                       eirene_init_input_blocks
      USE EIRMOD_COMXS, ONLY : MAXCRM, NCR_REAC, CR_KEY
!pgf      use json_module, ck => json_ck
!cym/cpg, keep using json_ck
!      use json_module           !IGNORE
      use json_module           !IGNORE
!cym/cpg

      IMPLICIT NONE

      INTERFACE
        SUBROUTINE EIRENE_IF0PRM_JSON(json,p)
        use json_module           !IGNORE
        class(json_core),intent(inout) :: json
        type(json_value), pointer, intent(in) :: p
        END SUBROUTINE EIRENE_IF0PRM_JSON
      END INTERFACE
c
      integer :: nstsi, j, nlimi
      INTEGER :: INDPRO(12)
      LOGICAL :: NLSCL, NLTEST, NLANA, NLDRFT, NLCRR, NLERG, NLIDENT,
     .           NLONE, NLMOVIE, NLCASCAD, NLDFST,
     .           NLRANMAR, NLOCTREE, NEXVS, NLTRIMESH
      logical :: found0, ldefstor
!      logical :: found0, ldefstor, lext, nlerg, lhyddef, ladapt
      logical :: lmulpl ! multiple Ti and V..IN (per species)
                        ! multiple ion velocities (per species)
cym
      character*8, allocatable :: textal(:)
cym
!pgf      character(kind=CK,len=:), allocatable :: header
!cym/cpg, keep using json_ck
      character(kind=json_CK,len=:), allocatable :: header
      EXTERNAL :: EIRENE_COUPLE_PARAM_CONSISTENCY,
     .            EIRENE_FIND_TRIANG_DIM, EIRENE_FIND_TET_DIM,
     .            EIRENE_INIT_PARAMS, EIRENE_LEER, EIRENE_EXIT_OWN

C
C  SET DEFAULT VALUES FOR STORAGE PARAMETERS
C
      CALL EIRENE_INIT_PARAMS
      LDEFSTOR=.false.

      lmulpl = .false.

      CALL EIRENE_LEER(3)

      WRITE (IUNOUT,*) 'PRINTOUT FROM EIRENE PRE-PROCESSING:'
      WRITE (IUNOUT,*) 'BROWSE INPUT FOR STORAGE NEEDS (FIND_PARAM.F)'
      WRITE (IUNOUT,*) 'FILE = ', 'eirene.input.json'
      CALL EIRENE_LEER(1)
C
c  start reading json file

      call eirene_init_input_blocks ('eirene.input.json', iunout)

C  start browsing the header

      j = itree_num(0)        ! index of the jtree that holds
                              ! the data of the input block
      call jtrees(j)%get(blks(0)%p,'TXTRUN',header,found0)
      write (iunout,'(A)') header
      CALL EIRENE_LEER(1)

!  browse block 1 "GENERAL_DATA" for parameters
      j = itree_num(1)
      call eirene_browse_block_1(jtrees(j),blks(1)%p)

!  browse block 2 "STANDARD_MESH" for parameters
      j = itree_num(2)
      call eirene_browse_block_2(jtrees(j),blks(2)%p)

      WRITE (iunout,*) '*** 3. DATA FOR BREP SURFACES'

!  browse block 3a "NONDEF_STD_SURFACES" for parameters
      j = itree_num(3)
      call eirene_browse_block_3a(jtrees(j),blks(3)%p)

!  browse block 3b "ADDITIONAL_SURFACES" for parameters
      j = itree_num(4)
      call eirene_browse_block_3b(jtrees(j),blks(4)%p)

!  browse block 4 and 5 "PHYSICS_MODEL" for parameters
      j = itree_num(5)
      call eirene_browse_block_4(jtrees(j),blks(5)%p)
      call eirene_browse_block_5(jtrees(j),blks(5)%p)

!  browse block 6 "REFLECTION_MODELS" for parameters
      j = itree_num(6)
      call eirene_browse_block_6(jtrees(j),blks(6)%p)

!  browse block 7 "SOURCES" for parameters
      j = itree_num(7)
      call eirene_browse_block_7(jtrees(j),blks(7)%p)

      WRITE (iunout,*) '*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES'

!  browse block 8 "SPECIFIC_ZONES" for parameters
!  nothing to be done
      j = itree_num(8)

!  browse block 9 "STATISTICS" for parameters
      j = itree_num(9)
      call eirene_browse_block_9(jtrees(j),blks(9)%p)

!  browse block 10 "ADDITIONAL_TALLIES" for parameters
      j = itree_num(10)
      call eirene_browse_block_10(jtrees(j),blks(10)%p)

!  browse block 11 "OUTPUT" for parameters
      j = itree_num(11)
      call eirene_browse_block_11(jtrees(j),blks(11)%p)

!  browse block 12 "DIAGNOSTICS" for parameters
      j = itree_num(12)
      call eirene_browse_block_12(jtrees(j),blks(12)%p)

!  browse block 13 "TIMEDEPENDENT_MODE" for parameters
      j = itree_num(13)
      call eirene_browse_block_13(jtrees(j),blks(13)%p)

!  browse block 14 "INTERFACING" for parameters
      if (nmode == 0) then
        j = itree_num(14)
        call eirene_browse_block_14(jtrees(j),blks(14)%p)
      else
! to be prepared
!       CALL EIRENE_IF0PRM(IUNIN)
        j = itree_num(14)
        call eirene_if0prm_json(jtrees(j),blks(14)%p)
      end if

      CALL EIRENE_ALLOC_CTRCEI(2)

C  SWITCH OFF SUM OVER STRATA IF THERE IS ONLY ONE STRATUM TO BE CALCULATED
!pb      IF (NSTRAI == 1) NSMSTRA = 0

cdr  due to these changes there, also some derived storage parmeters may have changed....
      NRAD=MAX(N1ST*N2ND*N3RD,NTRI*N3RD,NTETRA)+NADD+1 ! as in parmmod


!  output parameters for check

      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,'(A)') ' AUTOMATED STORAGE SETTING (FIND_PARAM.F)'
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
      WRITE (iunout,'(a14,i8)') 'NREAC       = ',NREAC
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
      ENDIF
C  TRY TO BE INTELLIGENT:  AUTOMATIC STORAGE REDUCTION....
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
cdr  time dependent options: census array size
      WRITE (iunout,'(a14,i8)') 'NPRNL       = ',NPRNL
C
      CALL EIRENE_LEER(2)
cpg
      call eirene_couple_param_consistency(nlimi,nstsi,textal,
     .                                     natm)
cpg
      return

      contains

!******************************************************************************


      subroutine eirene_browse_block_1 (json,p)

      implicit none
      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(inout) :: p
      type(json_value), pointer :: pfile, pdb
      integer :: i, ifile, nch
      logical :: found
!pgf      character(kind=CK,len=:), allocatable :: cdbh, cdbf
      character(kind=json_CK,len=:), allocatable :: cdbh, cdbf
      external :: eirene_filepath_usr

      WRITE (iunout,*) '*** 1. DATA FOR OPERATING MODE'

      call json%get(p,'NMODE',nmode,found)
      call json%get(p,'NTIME',ntime,found)
      call json%get(p,'NOPTIM',noptim,found)
      call json%get(p,'NOPTM1',noptm1,found)
      call json%get(p,'NGEOM_USR',ngeom_usr,found)
      call json%get(p,'NCOUP_INPUT',ncoup_input,found)
      call json%get(p,'NSMSTRA',nsmstra,found)
      call json%get(p,'NSTORAM',nstoram,found)
      call json%get(p,'NGSTAL',ngstal,found)
      call json%get(p,'NRTAL',nrtal,found)
      ldefstor = .true.
      nrtal_in = nrtal
      noptim_in = noptim
      nsmstra_in = nsmstra

      call json%get(p,'NLERG',nlerg,found)

C  NSTORAM IS REDEFINED, FINALLY EITHER =0  (A&M STORAGE SAVE MODE)
C                                    OR =9  (FULL A&M STORAGE MODE, =DEFAULT)
      NSTORAM = MIN(NSTORAM,9)
      IF (NSTORAM < 9) NSTORAM = 0
      NOPTM1 = MAX(NOPTM1,1)
C

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
      call json%get(p,'NLSPCSCL',nlspcscl,found)
      call json%get(p,'NLSPCSCL_ON',nlspcscl_on,found)
      call json%get(p,'NLSOLEDGE',nlsoledge,found)

      CALL EIRENE_INIT_CINIT

cdr scan for optional CFILE lines: path to external database files
      call json%get_child(p,'CFILE',pfile,found)
      if (found) then
        call json%info(pfile,n_children=nch)
!       do ifile = 1, ndbnames
        do i = 1, nch
          call json%get_child(pfile,i,pdb,found)
          call json%get(pdb,'FILE',cdbh,found)
          call json%get(pdb,'PATH',cdbf,found)
          DO IFILE = 1,NDBNAMES
            IF (INDEX(DBHANDLE(IFILE),cdbh) /= 0) EXIT
          END DO
          IF (IFILE <= NDBNAMES) THEN
            dbhandle(ifile) = trim(cdbh)
            dbfname_in(ifile) = trim(cdbf)
            dbfname(ifile) = trim(cdbf)
            call eirene_filepath_usr(cdbf,dbfname(ifile),
     .                               1,len_trim(cdbf))
            ldbread(ifile) = .true.
            WRITE (IUNOUT,'(2A)')
     .           'PATH SET FOR FILE ',trim(dbhandle(ifile))
            WRITE (IUNOUT,'(2A)') 'PATH = ',trim(dbfname(ifile))
          ELSE
            WRITE (IUNOUT,*) ' WRONG NAME FOR DATABASE ENTERED'
            WRITE (IUNOUT,*) ' DATABASE DEFINITION FOR ',
     .                         trim(dbhandle(ifile)),' IGNORED'
          end if
          deallocate(cdbh)
          deallocate(cdbf)
          nullify(pdb)
        end do
      end if

      return

      end subroutine eirene_browse_block_1


!******************************************************************************


      subroutine eirene_browse_block_2(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      INTEGER, allocatable :: indgrd(:)
      integer :: nr1st, nrplg, npplg, nrknot, ncoor, ntrii, ntet,
     .           np2nd, nt3rd, nttra, nbmlt, nradd, ll
      logical :: nlplg, nlfem, nltet, nlrad, nlpol, nltor, nlmlt, nladd,
     .           nltra
      logical :: found
      character(420) :: casename
!pgf      character(kind=CK,len=:), allocatable :: case
!cym/cpg, keep using json_ck
      character(kind=json_CK,len=:), allocatable :: case

      WRITE (iunout,*) '*** 2. DATA FOR VOXEL GRID GENERATION'

      call json%get(p,'INDGRD',indgrd,found)
      call json%get(p,'NLRAD',nlrad,found)
      call json%get(p,'NLPOL',nlpol,found)
      call json%get(p,'NLTOR',nltor,found)
      call json%get(p,'NLMLT',nlmlt,found)
      call json%get(p,'NLADD',nladd,found)
      nlpol_in = nlpol

C
C  RADIAL MESH
C
      call json%get(p,'RADIAL_GRID.NLPLG',nlplg,found)
      call json%get(p,'RADIAL_GRID.NLFEM',nlfem,found)
      call json%get(p,'RADIAL_GRID.NLTET',nltet,found)
      nlplg_in = nlplg
      nlfem_in = nlfem

      call json%get(p,'RADIAL_GRID.NR1ST',nr1st,found)
      nr1st_in = nr1st
      n1st = max(n1st,nr1st)

      if (nlplg) then
        call json%get(p,'RADIAL_GRID.NRPLG',nrplg,found)
        call json%get(p,'RADIAL_GRID.NPPLG',npplg,found)
        nplg = max(nplg,nrplg)
        nppart = max(nppart,npplg)
      end if

      call json%get(p,'RADIAL_GRID.NRKNOT',nrknot,found)
      call json%get(p,'RADIAL_GRID.NCOOR',ncoor,found)
      nknot = max(nknot,nrknot)
      ncoord = max(ncoord,ncoor)

      ntrii = 0
      ntet = 0
      if (indgrd(1) <= 5) then
        if (nlfem .or. nltet) then
          call json%get(p,'RADIAL_GRID.CASE',case,found)

          IF (.NOT.FOUND) THEN
            WRITE (IUNOUT,*) ' ERROR IN GEOMETRY SPECIFICATION '
            WRITE (IUNOUT,*)
     .           ' TRIANGLE OR TETRAHEDRON GRID SWITCHED ON '
            WRITE (IUNOUT,*) ' BUT NO CASENAME SPECIFIED '
            CALL EIRENE_EXIT_OWN(1)
          END IF

          LL = LEN_TRIM(CASE)
          CASENAME=REPEAT(' ',420)
          IF (LL > LEN(CASENAME)) THEN
            WRITE (IUNOUT,*) 'CASE NAME TOO LONG'
            CALL EIRENE_EXIT_OWN(1)
          END IF
          CASENAME(1:LL) = CASE(1:LL)

          IF (NLFEM) THEN
            CALL EIRENE_FIND_TRIANG_DIM(NR1ST, NTRI, NTRII, NKNOT,
     .                                  NGITT, CASENAME)
          ENDIF

          IF (NLTET) THEN
            CALL EIRENE_FIND_TET_DIM(NR1ST, NTETRA, NTET, NCOORD,
     .                               NGITT, CASENAME, IFOFF)
          ENDIF
        end if
      end if
      NTRIS=NTRII
      NKNOTS=NKNOT
      N1ST = MAX(N1ST,NR1ST)

C
C  POLOIDAL MESH
C
      call json%get(p,'POLOIDAL_GRID.NP2ND',np2nd, found)
      np2nd_in = np2nd
cdr  storage for 2nd coordinate grid only, if NLPOL=T
      if (nlpol) n2nd = max(n2nd,np2nd)

C
C  TOROIDAL MESH
C
      call json%get(p,'TOROIDAL_GRID.NLTRA',nltra,found)
      call json%get(p,'TOROIDAL_GRID.NT3RD',nt3rd,found)
      call json%get(p,'TOROIDAL_GRID.NTTRA',nttra,found)
      nt3rd_in = nt3rd
      if (nltor.and.nltra) nttra = nt3rd
cdr  storage for 3nd coordinate grid only, if NLTOR=T
      if (nltor) n3rd = max(n3rd,nt3rd)
CDR STORAGE FOR TOROIDAL EFFECTS, EVEN IF
CDR NO TOROIDAL RESOLUTION IS USED
      NTOR = MAX(NTOR,NTTRA)
C
C  MESH MULTIPLICATION
C
      call json%get(p,'MESH_MULTIPLICATION.NBMLT',nbmlt,found)
      NBMLT=MAX0(1,NBMLT)
      NBMAX=NBMLT
C
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
C
      if (nladd) then
        call json%get(p,'ADD_CELLS.NRADD',nradd,found)
        NADD = MAX(NADD,NRADD)
      end if

      return
      end subroutine eirene_browse_block_2

!******************************************************************************


      subroutine eirene_browse_block_3a(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      logical :: found

      WRITE (iunout,*) '*** 3A. DATA FOR NON-DEFAULT STANDARD SURFACES'

      call json%get(p,'NSTSI',nstsi,found)
!PB   IF NTIME >=1 NSTSI IS INCREASED in BLOCK 12
!PB   IF (NTIME.GE.1) NSTSI = NSTSI + 1
      NSTS = MAX(NSTS,NSTSI)

      return
      end subroutine eirene_browse_block_3a


!******************************************************************************


      subroutine eirene_browse_block_3b(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      logical :: found

      WRITE (iunout,*) '*** 3B. DATA FOR ADDITIONAL SURFACES'

      call json%get(p,'NLIMI',nlimi,found)
      NLIM = MAX(NLIM,NLIMI)

      return
      end subroutine eirene_browse_block_3b


!******************************************************************************


      subroutine eirene_browse_block_4(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: species, preac, pspecs, pspc,
     .                             preacs, prea
      integer :: nreaci,i,k
      logical :: found
      character(50) :: s1, s2
      character(kind=json_CK,len=:), allocatable :: zeile

      WRITE (iunout,*) '*** 4. DATA FOR SPECIES SPECIFICATION AND'
      WRITE (iunout,*) '       ATOMIC PHYSICS MODULE'

      s1 = 'SPECIES_SPEC'
      s2 = '.REACTIONS'

      WRITE (iunout,*)
     .  '       ATOMIC REACTION CARDS, NREACI DATA FIELDS'

      call json%get(p,s1//s2//'.NREACI',nreaci,found)
!PB   increase number of reactions by 1 as there are still
!PB   calls to SLREAC which use reaction number NREACI+1 (SGNAL and HE_EMISS)
      NREAC = MAX(NREAC,NREACI+1)
      
      call json%get_child(p,'SPECIES_SPEC',species,found)
      call json%get_child(species,'REACTIONS',preac,found)
      call json%get_child(preac,'REAC_SPECS',pspecs,found)
      NREAC_LINES=json%count(pspecs)

c  Count number of reaction data sets using internal coll. rad. models, currently H, He
cdr   NRC_REAC(1): no. of H_colrad reaction rates or pop-coef data,
cdr   NRC_REAC(4): no. of He_colrad reaction rates or pop-coef data
      NCR_REAC = 0
      DO k = 1, NREAC_LINES
        call json%get_child(pspecs,k,prea,found)
        call json%get(prea,'FILNAM',zeile,found)
        IF (INDEX(zeile,json_CK_'CR') > 0) THEN
!  coll. rad. model found
          IF (INDEX(zeile,json_CK_'CR_') > 0) THEN
!  coll. rad. model specified in more detail
            DO i = MAXCRM, 1, -1
              IF (LEN_TRIM(CR_KEY(I)) > 0) THEN
                IF (INDEX(zeile,TRIM(CR_KEY(I))) > 0)
     .            NCR_REAC(I) = NCR_REAC(I) + 1
              END IF
            END DO
          ELSE
!  use default coll. rad. model for hydrogen H, formulation II
          NCR_REAC(1) = NCR_REAC(1) + 1
          END IF
        ENDIF
        nullify(prea)
        deallocate(zeile)
      END DO

!  ATOMS

      WRITE (iunout,*)
     .  '*** 4A. NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES'

      call json%get_child(species,'ATOMS',pspc,found)
      call json%get(pspc,'NATMI',natmi,found)
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

      call eirene_read_spc_block (json,pspc,'A',natmi,lmulpl,
     .                            nmassa, nchara)

!  MOLECULES

      WRITE (iunout,*)
     .  '*** 4B. NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES'

      call json%get_child(species,'MOLECULES',pspc,found)
      call json%get(pspc,'NMOLI',nmoli,found)
      NMOL = MAX(NMOL,NMOLI)
      if(.not.allocated(NMASSM)) then
        allocate(NMASSM(NMOL))
        NMASSM = 0
        COMUSR_FIRST_PASS(3) = .FALSE.
      end if

      call eirene_read_spc_block (json,pspc,'M',nmoli,lmulpl,nmassm)

!  TEST IONS

      WRITE (iunout,*)
     . '*** 4C. TEST IONS SPECIES CARDS, NIONI SPECIES'

      call json%get_child(species,'TEST_IONS',pspc,found)
      call json%get(pspc,'NIONI',nioni,found)
      NION = MAX(NION,NIONI)
      if(.not.allocated(NMASSI)) then
        allocate(NMASSI(NION))
        NMASSI = 0
        allocate(NCHARI(NION))
        NCHARI = 0
        allocate(NCHRGI(NION))
        NCHRGI = 0
        allocate(LKINDI(NION))
        LKINDI = 0
        COMUSR_FIRST_PASS(4) = .FALSE.
      end if

      call eirene_read_spc_block (json,pspc,'I',nioni,lmulpl,
     .                            nmassi, nchari, nchrgi, lkindi)

!  PHOTONS

      WRITE (iunout,*)
     .  '*** 4D. PHOTONS SPECIES CARDS, NPHOTI SPECIES '

      call json%get_child(species,'PHOTONS',pspc,found)
      call json%get(pspc,'NPHOTI',nphoti,found)
      NPHOT = MAX(NPHOT,NPHOTI)

      call eirene_read_spc_block (json,pspc,'PH',nphoti,lmulpl)

      nullify(species)
      nullify(preac)
      nullify(pspecs)
      nullify(pspc)

      return
      end subroutine eirene_browse_block_4

!******************************************************************************

      subroutine eirene_read_spc_block (json,p,ch,nspc,logpl,
     .                                  nmass, nchar, nchrg, lkind)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      character(*), intent(in) :: ch
      integer, intent(in) :: nspc
      logical, intent(inout) :: logpl
      integer, intent(inout), optional :: nmass(nspc), nchar(nspc),
     .                                    nchrg(nspc), lkind(nspc)

      type(json_value), pointer :: pspecies, spchild, preac, pir
      integer :: nrc, ibgk, i, j
      logical found
!pgf      character(kind=CK,len=:), allocatable :: spname
!cym/cpg, keep using json_ck
      character(kind=json_CK,len=:), allocatable :: spname

      call json%get_child(p,'SPECIES',pspecies,found)
      call json%get_child(pspecies,spchild)
      if (.not.associated(spchild)) return

      do i = 1, nspc
        if (ch == 'A') then
          call json%get(spchild,'SPECIES',spname,found)
          if (found) then
            textal(i)(1:8) = spname(1:8)
            deallocate(spname)
          else
            textal(i)=repeat(' ',8)
          end if
        end if

        if (present(nmass)) then
          call json%get(spchild,'NMASS'//ch,nmass(i),found)
        end if
        if (present(nchar)) then
          call json%get(spchild,'NCHAR'//ch,nchar(i),found)
        end if
        if (present(nchrg)) then
          call json%get(spchild,'NCHRG'//ch,nchrg(i),found)
        end if
        if (present(lkind)) then
          call json%get(spchild,'LKIND'//ch,lkind(i),found)
        end if

        call json%get(spchild,'NRC'//ch,nrc,found)
        if (.not.found) cycle

        call json%get_child(spchild,'REACTIONS',preac,found)
        if (.not.found) cycle

        call json%get_child(preac,pir)
        do j = 1, nrc
          if (.not.associated(pir)) exit
cdr:  try to identify if there are so called BKG collisions, input flag IBGK::
cdr:  to be generalized: there may be other reactions, which require multiple Ti profiles
          call json%get(pir,'IBGK'//ch,ibgk,found)
          logpl = logpl .or. (ibgk /= 0)
          call json%get_next(pir,pir)
        end do
        call json%get_next(spchild,spchild)
        if (.not.associated(spchild)) exit
      end do

      nullify (pspecies)
      nullify (spchild)
      nullify (preac)
      nullify (pir)

      return
      end subroutine eirene_read_spc_block


!******************************************************************************


      subroutine eirene_browse_block_5(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: pback, pbulk, pspecies,
     .                             pspc, pdm, pplas
      integer :: i, ico
      INTEGER, allocatable :: indp(:)
      logical :: found

      WRITE (iunout,*) '*** 5. DATA FOR PLASMA BACKGROUND'

!  BULK IONS

      WRITE (iunout,*) '*** 5A. BULK ION SPECIES CARDS, NPLSI SPECIES'

      call json%get_child(p,'BACKGROUND',pback,found)
      call json%get_child(pback,'BULK_IONS',pbulk,found)
      call json%get(pbulk,'NPLSI',nplsi,found)
      NPLS = MAX(NPLS,NPLSI)

      call json%get_child(pbulk,'SPECIES',pspecies,found)
      call json%get_child(pspecies,pspc)

      ico = 0
      do i = 1, nplsi
        if (associated(pspc)) then
          call json%get_child(pspc,'DENS_MODEL',pdm,found)
          if (found) ico = ico + 1
        end if
        call json%get_next(pspc,pspc)
      end do

!  enhance number of reactions by 1 for use by density model
      if (ico > 0) nreac = nreac + 1

!  PLASMA

      WRITE (IUNOUT,*) '*** 5B. PLASMA BACKGROUND DATA'

      call json%get_child(pback,'PLASMA',pplas,found)
      call json%get(pplas,'INDPRO',indp,found)
      if (found) indpro(1:12) = indp(1:12)

      call eirene_Ti_input(json,pplas,indpro2_save)
      indpro(2)=indpro2_save

      NPLSTI = NPLS  ! now same as for VI.  good
      IF (INDPRO(2)<0) NPLSTI=1

      IF ((NPLS > 1) .AND. (NPLSTI == 1) .AND. LMULPL) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'TIIN PROVIDED FOR ONE SPECIES ONLY'
        WRITE (IUNOUT,*) 'DUE TO INDPRO(2) > 10'
        WRITE (IUNOUT,*) 'DIMENSION OF TIIN OVERWRITTEN'
        WRITE (IUNOUT,*) 'BECAUSE BGK REACTIONS ARE PRESENT'
        NPLSTI = NPLS
        WRITE (IUNOUT,*) ' NPLSTI = ',NPLSTI
      END IF

cdr these next 2 lines for V_IN(ipls)
      NPLSV = NPLS
      IF (MOD(ABS(INDPRO(4)),100) > 9) NPLSV = 1

      IF ((NPLS > 1) .AND. (NPLSV == 1) .AND. LMULPL) THEN
        WRITE (IUNOUT,*) 'WARNING FROM FIND_PARAM'
        WRITE (IUNOUT,*) 'V_IN PROVIDED FOR ONE SPECIES ONLY',
     .                   'DUE TO INDPRO(4) > 10'
        WRITE (IUNOUT,*) 'DIMENSION OF V_IN ARRAYS OVERWRITTEN',
     .                   'BECAUSE BGK REACTIONS ARE  PRESENT'
        NPLSV = NPLS
        WRITE (IUNOUT,*) ' NPLSV = ',NPLSV
      END IF

      deallocate(indp)
      nullify(pback)
      nullify(pbulk)
      nullify(pspecies)
      nullify(pspc)
      nullify(pdm)
      nullify(pplas)

      return
      end subroutine eirene_browse_block_5


      SUBROUTINE eirene_Ti_input(json,plsm,indpro2_save)
      implicit none
      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: plsm
      integer, intent(out) :: indpro2_save
      type(json_value), pointer :: pti

      integer :: icti, indpro2
      logical :: lfti

cdr Nov. 2019
cdr Achieve synchronisation of input options for multi-species Ti and Vi.
cdr Return a "best guess" of what INDPRO(2) should be.

cdr Ti_IN(ipls) options have historically been just opposite to V_IN(ipls) options.
cdr This has gotten amplified to a long-lasting code mess:
cdr Increasingly inconsistent input options for Ti(ipls) and V..IN(ipls) profiles
cdr and false code operation (in multi-fluid cases) in many instances.
cdr Now: We try to automatically set TI flag INDPRO(2) by counting input cards
cdr in this block.

      if (NPLS.eq.1) then
c  no ambiguity possible here...
        indpro2_save=iabs(indpro(2))
      endif

      write (iunout,*) 'indpro2, npls ',indpro(2),npls

c remove sign and second or third digits
      indpro2=mod(iabs(indpro(2)),10)  ! now within 1 and 9

      if (indpro2.eq.6) then
c  Ti(ipls) is from external file.
c  Try to find out what was meant:
        if (indpro(2).lt.0) then       ! indpro2=-6, -16,...
cdr  just trying, via CI
          indpro2_save=6  ! maybe this was meant in solps-iter cases?
        elseif (indpro(2).lt.10) then  ! indpro2=6
          indpro2_save=16 ! maybe this was meant in solps-iter cases?
        else                           ! indpro2=16,26,36,...
        endif
c  no card counting done.
c       indpro2_save= ??
      endif

C  FIND
C        START OF NEXT INPUT BLOCK: 6,
C  OR    OF OPTIONAL INPUT LINES,
C  or,   in case of external block 45, L45: of end of file
cdr  and count the remaining input lines for plasma profiles in block 5.

      call json%get_child(plsm,'TI',pti,lfti)

      icti = 0
      if (lfti) call json%info(pti,n_children=icti)

      write (iunout,*) 'cards for Ti ',icti
      if (npls.eq.1.or.indpro2.eq.6) return

      IF (icti.lt.npls.and.icti.ge.1) then
c  we necessarily have data for only one single Ti in the input file.
        if (indpro(2).lt.10) then
c         if (indpro(2).lt.0) ...  input is in K rather than eV. Not used.
          indpro(2)=iabs(indpro(2))
          indpro(2)=10+indpro(2)
          write (iunout,*) 'apparently only one Ti profile is given'
          write (iunout,*) 'indpro(2) reset to: ',indpro(2)
        endif
      elseif (icti.eq.npls) then
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
        write (iunout,*) 'No. of input lines for Ti ',icti
        write (iunout,*) 'Choose INDPRO(2) > 10 (NPLSTI=1)'
        write (iunout,*) 'and assume that a VOL card is present.'
        indpro(2) = mod(indpro(2),10) + 10
        write (iunout,*) 'IF NPLSTI=2 was intended, add the VOL-card'
        write (iunout,*) 'i.e. add a card, e.g. containing just 0.0'
      endif
      indpro2_save=indpro(2)

      nullify (pti)
      return
      end subroutine eirene_Ti_input

!******************************************************************************


      subroutine eirene_browse_block_6(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: ptom
      logical :: found, nltrim

      WRITE (iunout,*) '*** 6. GENERAL DATA FOR REFLECTION MODEL'

      call json%get(p,'NLTRIM',nltrim,found)

      if (nltrim) then
         call json%get_child(p,'PROJ_ON_MAT',ptom,found)
      else
         found = .false.
      end if

      if (found) then
        NHD6 = json%count(ptom)
      else
        NHD6 = 12
      end if

      nullify(ptom)

      return
      end subroutine eirene_browse_block_6


!******************************************************************************


      subroutine eirene_browse_block_7(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: pstra, pst, psub, psubs
      integer :: j, noc, nsrfsi, isor, ist, i,indim
      real(dp) :: sorlim, sorind
      integer, allocatable :: indsrc(:)
      logical :: found, fsub
!pgf      character(kind=CK,len=:),allocatable :: txt
!cym/cpg, keep using json_ck
      character(kind=json_CK,len=:),allocatable :: txt

      WRITE (iunout,*) '*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA'

      call json%get(p,'NSTRAI',nstrai,found)
      NSTRAI_IN = NSTRAI
      IF (NTIME.GE.1) NSTRAI = NSTRAI + 1
      NSTRA = MAX(NSTRA,NSTRAI)

CDR  TRY TO SET NSTEP, THE NUMBER OF STEP FUNCTIONS FOR SOURCE SAMPLING
cdr  set nstep = highest stratum number, which receives primary source data from external code.
cdr  this must be highly case specfic. To be reconsidered !!
cpb  Step functions are only used in conjunction withsource sampling.
cpb  If we find the highest stratum J that uses a step function we can be sure
cpb  that there are less than J step functions involved.
cpb  For this reason start the loop at NSTRAI and count downwards.
cpb  We imply implicitely that the numbering of step functions is in
cpb  ascending order starting with 1
      NSTEP = 1
      call json%get(p,'INDSRC',indsrc,found)
      IF (ANY(INDSRC == 6)) THEN
        DO J=NSTRAI_IN,1,-1
          IF (INDSRC(J) == 6) THEN
            NSTEP = J
            EXIT
          END IF
        END DO
      END IF

      call json%get_child(p,'STRATA',pstra,found)

      if (found) then
!  noc stratum blocks in input file
!  noc < nstrai if any(indsrc==6)
         noc = json%count(pstra)
         call json%get_child(pstra,pst)
         do j = 1, noc
           if (associated(pst)) then
             call json%get(pst,'TXTSOU',txt,found)
             WRITE (IUNOUT,'(A1,A)') ' ',trim(txt)
             deallocate(txt)
!  number of substrata
             call json%get(pst,'NSRFSI',nsrfsi,found)
             NSRFS = MAX(NSRFS,NSRFSI)

             call json%get_child(pst,'SUBSTRATA',psubs,fsub)

             call json%get_child(psubs,psub)
             do i = 1, nsrfsi
                if (associated(psub)) then
!  find maximum step function index
                  call json%get(psub,'INDIM',indim,found)
                  call json%get(psub,'SORLIM',sorlim,fsub)
                  call json%get(psub,'SORIND',sorind,fsub)
                  IF (INDIM == 4) NSTEP = MAX(NSTEP,NINT(SORIND))
                  isor = nint(sorlim)
                  do while (isor > 0)
                    ist =  mod(isor,10)
cdr  here NSTEP is set to the largest step function number specified on SORIND
                    if ((ist == 4).OR.(ist==5))
     .                nstep = max(nstep,nint(sorind))
                    isor = isor / 10
                  end do
                  call json%get_next(psub,psub)
                end if
             end do

             call json%get_next(pst,pst)
           end if
         end do
      end if

      deallocate (indsrc)
      nullify(pstra)
      nullify(pst)
      nullify(psub)
      nullify(psubs)

      return
      end subroutine eirene_browse_block_7


!******************************************************************************


      subroutine eirene_browse_block_9(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nsigvi, nsigsi, nsigci
      logical :: found

      WRITE (iunout,*) '*** 9. DATA FOR STATISTIC AND NON-ANALOG MODEL'

      WRITE (iunout,*) '       CARDS FOR STANDARD DEVIATION'

      if (associated(p)) then
        call json%get(p,'NSIGVI',nsigvi,found)
        call json%get(p,'NSIGSI',nsigsi,found)
        call json%get(p,'NSIGCI',nsigci,found)
      else
        nsigvi = 0
        nsigsi = 0
        nsigci = 0
      end if

      NSD = MAX(NSD,NSIGVI)
      NSDW = MAX(NSDW,NSIGSI)
      NCV = MAX(NCV,NSIGCI)

      return
      end subroutine eirene_browse_block_9


!******************************************************************************


      subroutine eirene_browse_block_10(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nadvi, nclvi, nalvi, nadsi, nalsi, nadspc
      logical :: found

      WRITE (iunout,*)
     . '*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION'
      WRITE (iunout,*)
     . '        ESTIMATORS AND ALGEBRAIC EXPRESSIONS'

      if (associated(p)) then
        call json%get(p,'NADVI',nadvi,found)
        call json%get(p,'NCLVI',nclvi,found)
        call json%get(p,'NALVI',nalvi,found)
        call json%get(p,'NADSI',nadsi,found)
        call json%get(p,'NALSI',nalsi,found)
        call json%get(p,'NADSPC',nadspc,found)
      else
        nadvi = 0
        nclvi = 0
        nalvi = 0
        nadsi = 0
        nalsi = 0
        nadspc = 0
      end if

      NADV = MAX(NADV,NADVI)
      NCLV = MAX(NCLV,NCLVI)
      NALV = MAX(NALV,NALVI)
      NADS = MAX(NADS,NADSI)
      NALS = MAX(NALS,NALSI)

      WRITE (iunout,*) '*** 10A. DATA FOR ADDITIONAL TALLIES'
      WRITE (iunout,*) '*** 10B. DATA FOR COLLISION ESTIMATORS'
      WRITE (iunout,*) '*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS'
      WRITE (iunout,*) '*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES'
      WRITE (iunout,*) '*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES'
      WRITE (iunout,*) '*** 10F. DATA FOR SPECTRA'

      return
      end subroutine eirene_browse_block_10


!******************************************************************************


      subroutine eirene_browse_block_11(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: ptals, ptal
      integer :: nvolpr, nsurpr, nvolpl, i, nsp
      logical :: found, lrpscut

      WRITE (iunout,*)
     .  '*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT'

      call json%get(p,'NVOLPR',nvolpr,found)
      call json%get(p,'NSURPR',nsurpr,found)

      NVLPR=NVOLPR
C  ERGODIC OPTION NEEDS PRINTOUT OF VOLUME, AND ONE, TWO OR THREE FURTHER TALLIES AT LEAST
      IF (NLERG) NVLPR=MAX(4,NVLPR)
      NSRPR=NSURPR
C  ERGODIC OPTION NEEDS PRINTOUT AT LEAST FROM TIME HORIZON
      IF (NLERG) NSRPR=MAX(1,NSRPR)

      CALL EIRENE_ALLOC_CTRCEI(1)

      call json%get(p,'TRCAMD',trcamd,found)
      call json%get(p,'TRCINT',trcint,found)

C  ERGODIC OPTION NEEDS PRINTOUT OF VOLUME, AND ONE, TWO OR THREE FURTHER TALLIES AT LEAST
      IF (NLERG) NVLPR=MAX(4,NVLPR)

      NSRPR=NSURPR
C  ERGODIC OPTION NEEDS PRINTOUT AT LEAST FROM TIME-HORIZON
      IF (NLERG) NSRPR=MAX(1,NSRPR)

      call json%get(p,'LRPSCUT',lrpscut,found)
      call json%get(p,'NVOLPL',nvolpl,found)

      NPLT = 1
      if (nvolpl > 0) then
        call json%get_child(p,'PLOT_TALLIES',ptals,found)
        call json%get_child(ptals,ptal)

        do i = 1, nvolpl
          if (associated(ptal)) then
            call json%get(ptal,'NSP',nsp,found)
            NPLT = MAX(NPLT, NSP)
            call json%get_next(ptal, ptal)
          end if
        end do

        nullify(ptals)
        nullify(ptal)
      end if

      return
      end subroutine eirene_browse_block_11


!******************************************************************************


      subroutine eirene_browse_block_12(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      type(json_value), pointer :: plines, pline, pcomps, pcomp,
     .                             pcnts, pcnt
!pgf      character(kind=CK,len=:), allocatable :: ckey
!cym/cpg, keep using json_ck
      character(kind=json_CK,len=:), allocatable :: ckey
      integer :: num_lines, mod_addv, i, j, k, num_compo, num_contrib,
     .           iratio, nadv_add, nchori, ncheni
      logical :: lkey, found, nlemis

      WRITE (iunout,*) '*** 12. DATA FOR DIAGNOSTIC MODULE'

      NADV_ADD = 0
      NLEMIS = .FALSE.
      NUM_LINES = 0
      NCHORI = 0
      NCHENI = 0

      if (associated(p)) then
        call json%get(p,'KEYWORD',ckey,lkey)

      !   if (lkey) then
!cym type mistmatch here
!          if (index(ckey,'DEFINE_LINES') == 0) lkey = .false.
      !     if (index(ckey,json_ck_'DEFINE_LINES') == 0) lkey = .false.
!cym end
      !   end if

        if (lkey.and.(index(ckey,json_ck_'DEFINE_LINES') > 0)) then
          NLEMIS = .TRUE.
!  read emission lines
          call json%get(p,'.NUM_LINES',num_lines,found)
          call json%get(p,'.MOD_ADDV',mod_addv, found)

          call json%get_child(p,'LINES',plines,found)
          if (associated(plines)) then
            call json%get_child(plines,pline)
          end if

!  lines
          do i = 1, num_lines
            if (associated(pline)) then
              call json%get(pline,'NUM_COMPO',num_compo,found)
              IF (MOD_ADDV == 0) THEN
                NADV_ADD = MAX(NADV_ADD,NUM_COMPO+1)
              ELSE
                NADV_ADD = NADV_ADD + NUM_COMPO + 1
              END IF

              if (.false.) then  ! read components and contributions,
                                 ! not needed
! components
              call json%get_child(pline,'COMPONENTS',pcomps,found)
              call json%get_child(pcomps,pcomp)
              do j = 1, num_compo
                if (associated(pcomp)) then
                  call json%get(pcomp,'NUM_CONTRIB',num_contrib,found)

! contributions
                  call json%get_child(pcomp,'CONTRIBUTIONS',pcnts,found)
                  call json%get_child(pcnts,pcnt)
                  do k = 1, num_contrib
                    if (associated(pcnt)) then
                      call json%get(pcnt,'IRATIO',iratio,found)
                      call json%get_next(pcnt,pcnt)
                    end if
                  end do

                  call json%get_next(pcomp,pcomp)
                end if
              end do
              end if  ! if .false.

              call json%get_next(pline,pline)
            end if
          end do

! ADD 1 FOR TOTAL  --  already done
!         IF (MOD_ADDV == 0) NADV_ADD = NADV_ADD + 1
          NADV = NADV + NADV_ADD
!         NREAC = NREAC + LINES

        elseif (lkey.and.(index(ckey,json_ck_'DEFAULT_LINES') > 0)) then
          NLEMIS = .TRUE.
        end if

        WRITE (iunout,*) '*** 12A. DATA FOR LINES OF SIGHT'

        call json%get(p,'NCHORI',nchori,found)
        call json%get(p,'NCHENI',ncheni,found)

        nullify(plines)
        nullify(pline)
        nullify(pcomps)
        nullify(pcomp)
        nullify(pcnts)
        nullify(pcnt)

      end if

      NCHOR = MAX(NCHOR,NCHORI)
      NCHEN = MAX(NCHEN,NCHENI)
      NLEMIS = NLEMIS .OR. (NCHOR > 0)
      IF (NLEMIS.AND.(NUM_LINES == 0)) THEN
! USE DEFAULT LINES FOR EMISSIVITY
        MOD_ADDV = 0
        NADV=NADV+7
        NUM_LINES = 6
        NUM_COMPO = 6
! USE MAXIMUM AS NCHAR AND NCHRG ARE NOT YET AVAILABLE
        NUM_CONTRIB = NATMI + NMOLI + 2*NMOLI + 2*NMOLI + 2*NMOLI +
     .                NPLSI
        NREAC = NREAC + NUM_LINES*NUM_COMPO + 3
      END IF

      return
      end subroutine eirene_browse_block_12


!******************************************************************************


      subroutine eirene_browse_block_13(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: nprmul, nsnvi
      logical :: found

      WRITE (iunout,*)
     .  '*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION'

      nprnli = 0
      nprmul = 1
      NSNVI = 0

      if (associated(p)) then
        call json%get(p,'NPRNLI',nprnli,found)
        call json%get(p,'DEFAULT_TIME_HORIZON',ldef_time_horizon,found)
        call json%get(p,'NPRMUL',nprmul,found)
        if (nprnli > 0) call json%get(p,'NSNVI',nsnvi,found)
      end if

      IF (NPRMUL > 1) NPRNLI = NPRNLI * NPRMUL
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

      return
      end subroutine eirene_browse_block_13


!******************************************************************************


      subroutine eirene_browse_block_14(json,p)

      class(json_core),intent(inout) :: json
      type(json_value), pointer, intent(in) :: p

      integer :: naini, ncopie, ncpvi
      logical :: found

      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'

      IF (NMODE.EQ.0) THEN
        call json%get(p,'NAINI',naini,found)
!pb     call json%get(p,'NCOPII',ncopii,found)
        call json%get(p,'NCOPIE',ncopie,found)
        NCPVI=NCOPIE
      ELSE
        NAINI=0
        NCPVI=0
!pb     CALL EIRENE_IF0PRM(IUNIN)
        j = itree_num(14)
        call eirene_if0prm_json(json,p)
      ENDIF

cdr  some parameters may have gotten changed in IF0PRM, case-specific
      NAIN = MAX(NAIN,NAINI)
      NCPV = MAX(NCPV,NCPVI)

      return
      end subroutine eirene_browse_block_14


!******************************************************************************
#else
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      EXTERNAL :: EIRENE_EXIT_OWN

      WRITE(IUNOUT,*) 'ATTEMPTING TO READ JSON FORMAT INPUT'
      WRITE(IUNOUT,*) 'BUT CODE IS NOT LINKED TO THE JSON LIBRARY!'
      CALL EIRENE_EXIT_OWN(1)
      RETURN
#endif

      END SUBROUTINE EIRENE_FIND_PARAM_JSON
