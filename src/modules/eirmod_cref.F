      MODULE EIRMOD_CREF
cdr Data for surface reflection models

cdr  june 17: remove parameter NFLR. Redundant. Was same as NHD6
cdr  may 18 : perhaps unnecessary hard-coding of dimensions in some arrays,
cdr           rather than using parameters NHD1,...NHD5 ?

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_CREF, EIRENE_DEALLOC_CREF,
     .          EIRENE_INIT_CREF, EIRENE_BROADCAST_CREF,
     .          EIRENE_DEALLOC_REFLIST,
     .          REFMODEL, SPEC_REF

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCREF(:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICREF(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R RPROB0,    ERMIN,     ERCUT,
     R ENAR(:),   DENAR(:),  WIAR(:),   DWIAR(:),  RAAR(:),
     R DRAAR(:),
     R TM(:),     TC(:),     WM(:),     WC(:),    ERDC(:),  HFTR3F(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R RINTEG(:), EINTEG(:), AINTEG(:)

      REAL(DP), ALLOCATABLE, PUBLIC :: HFTR0(:,:,:),
     .                                 HFTR1(:,:,:,:),
     .                                 HFTR2(:,:,:,:,:),
     .                                 HFTR3(:,:,:,:,:,:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I INE, INEM, INW, INWM, INR, INRM

      INTEGER, PUBLIC, SAVE :: NCREF, MCREF

      LOGICAL, PUBLIC, SAVE :: LTRIM_OLD    !  flag:
c                              true:
c                                 use old trim reflection data file,
c                                 nhd6=12 projectile-target combinations in
c                                 one single data file TRIM.DAT
c                                 NHD6=12: hard-wired in: FIND_PARAM.f
c                              false:
c                                 read NFR different trim data files A_on_B
c                                 path and names "A_on_B" as specified in input.f
c                                 set NHD6=NFR
c                                 NHD6: dimensioning in allocatable arrays
c

      CHARACTER(500), PUBLIC, ALLOCATABLE, SAVE :: REFFIL(:)

! TYPE REFMODEL MOVED HERE BECAUSE OF WRITING OF JSON-FILE
      TYPE SPEC_REF
        CHARACTER(8) :: VARNAME
        CHARACTER(8) :: SPCNAME
        INTEGER :: IVAL
        REAL(DP) :: RVAL
        TYPE(SPEC_REF), POINTER :: NEXT
      END TYPE SPEC_REF

      TYPE REFMODEL
CVK VARIABLES FOR CHEMICAL SPUTTERING: ESPTCR,JLCHSPNWL,STPRMR(:)
        CHARACTER(70) :: REFNAME
        INTEGER       :: JLSPT
        INTEGER       :: REFCELL
        INTEGER, DIMENSION(:), POINTER :: JLREF,JSRS,JSRC,JLCHSPNWL
        REAL(DP)      :: ZNMLR,EWALLR,EWBINR,FSHEATR
        REAL(DP), DIMENSION(:,:), POINTER :: TRANSPR
        REAL(DP), DIMENSION(:), POINTER ::
     .                                   RCYCFR,RCYCTR,RCPRMR,
     .                                   EXPPLR,EXPELR,EXPILR,
     .                                   RCYCSR,RCYCCR,STPRMR,
     .                                   ESPTSR,ESPTCR
        REAL(DP):: REFPRESS
        TYPE(SPEC_REF), POINTER :: SPEC_LINES
        TYPE(REFMODEL),POINTER :: NEXT
      END TYPE REFMODEL

      TYPE(REFMODEL), POINTER, PUBLIC, SAVE :: REFLIST

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_CREF

      IF (ALLOCATED(RCREF)) RETURN

C  storage for TRIM database reflection model
      NCREF = 3+12+11+7+6+5+4+5*NHD6+NHD5
      MCREF  = 6

      ALLOCATE (RCREF(NCREF))
      ALLOCATE (ICREF(MCREF))
      ALLOCATE (REFFIL(NHD6))

      ALLOCATE (HFTR0(NHD1,NHD2,NHD6))
      ALLOCATE (HFTR1(NHD1,NHD2,NHD3,NHD6))
      ALLOCATE (HFTR2(NHD1,NHD2,NHD3,NHD4,NHD6))
      ALLOCATE (HFTR3(NHD1,NHD2,NHD3,NHD4,NHD5,NHD6))

C  storage for general reflection model (both: "Behrisch-Matrix" and database models)
C  generalized to become surface-specific (allows for different values per surface)
      ALLOCATE (RINTEG(0:NLIMPS))
      ALLOCATE (EINTEG(0:NLIMPS))
      ALLOCATE (AINTEG(0:NLIMPS))

      WRITE (IUNMEM,'(A,T25,I15)')
     .       ' CREF ',(NCREF+3*(NLIMPS+1))*8 + MCREF*4 +
     .                 NHD6*500 +  ! reffile array: character-string
                                   ! for path to datafiles
     .                (NHD1*NHD2*NHD6*(1+NHD3*(1+NHD4*(1+NHD5))))*8

      RPROB0    => RCREF(1)
      ERMIN     => RCREF(2)
      ERCUT     => RCREF(3)
cdr these next arrays contain hard-coded assumptions on parameters NHDx.
cdr Why can we not use NHDx directly?
      ENAR      => RCREF( 4:15)  ! here for NHD1=12.
      DENAR     => RCREF(16:26)  ! here for NHD1=12, 12-1=11
      WIAR      => RCREF(27:33)  ! here for NHD2=7.
      DWIAR     => RCREF(34:39)  ! here for NID2=7, 7-1=6
      RAAR      => RCREF(40:44)  !              =5
      DRAAR     => RCREF(45:48)  !                  5-1=4

      TM        => RCREF(49+0*NHD6 : 48+1*NHD6)
      TC        => RCREF(49+1*NHD6 : 48+2*NHD6)
      WM        => RCREF(49+2*NHD6 : 48+3*NHD6)
      WC        => RCREF(49+3*NHD6 : 48+4*NHD6)
      ERDC      => RCREF(49+4*NHD6 : 48+5*NHD6)

      HFTR3F    => RCREF(49+5*NHD6 : 48+5*NHD6+NHD5)

      INE  => ICREF(1)
      INEM => ICREF(2)
      INW  => ICREF(3)
      INWM => ICREF(4)
      INR  => ICREF(5)
      INRM => ICREF(6)


      CALL EIRENE_INIT_CREF

      RETURN
      END SUBROUTINE EIRENE_ALLOC_CREF


      SUBROUTINE EIRENE_DEALLOC_CREF

      IF (.NOT.ALLOCATED(RCREF)) RETURN

      DEALLOCATE (RCREF)
      DEALLOCATE (ICREF)
      DEALLOCATE (REFFIL)

      DEALLOCATE (RINTEG)
      DEALLOCATE (EINTEG)
      DEALLOCATE (AINTEG)

      DEALLOCATE (HFTR0)
      DEALLOCATE (HFTR1)
      DEALLOCATE (HFTR2)
      DEALLOCATE (HFTR3)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CREF


      SUBROUTINE EIRENE_DEALLOC_REFLIST
      TYPE(REFMODEL), POINTER :: REFCUR

      REFCUR => REFLIST
      DO WHILE (ASSOCIATED(REFCUR))
        REFLIST => REFLIST%NEXT
        DEALLOCATE (REFCUR%JLREF)
        DEALLOCATE (REFCUR%JSRS)
        DEALLOCATE (REFCUR%JSRC)
        DEALLOCATE (REFCUR%JLCHSPNWL)
        DEALLOCATE (REFCUR%TRANSPR)
        DEALLOCATE (REFCUR%RCYCFR)
        DEALLOCATE (REFCUR%RCYCTR)
        DEALLOCATE (REFCUR%RCPRMR)
        DEALLOCATE (REFCUR%EXPPLR)
        DEALLOCATE (REFCUR%EXPELR)
        DEALLOCATE (REFCUR%EXPILR)
        DEALLOCATE (REFCUR%RCYCSR)
        DEALLOCATE (REFCUR%RCYCCR)
        DEALLOCATE (REFCUR%STPRMR)
        DEALLOCATE (REFCUR%ESPTSR)
        DEALLOCATE (REFCUR%ESPTCR)
        DEALLOCATE (REFCUR)
        REFCUR => REFLIST
      ENDDO
      END SUBROUTINE EIRENE_DEALLOC_REFLIST


      SUBROUTINE EIRENE_INIT_CREF

      RCREF  = 0._DP
      ICREF  = 0
      REFFIL = ' '

      RINTEG = 0._DP
      EINTEG = 0._DP
      AINTEG = 0._DP

      HFTR0 = 0._DP
      HFTR1 = 0._DP
      HFTR2 = 0._DP
      HFTR3 = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_CREF

      SUBROUTINE EIRENE_BROADCAST_CREF(ME)
      USE EIRMOD_MPI
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
      EXTERNAL :: EIRENE_CHECK_EXIT
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      IF (ME /= 0) CALL EIRENE_ALLOC_CREF

      CALL MPI_BCAST (RINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (EINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (AINTEG,NLIMPS+1,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (RCREF,NCREF,MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (ICREF,MCREF,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (LTRIM_OLD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (HFTR0,NHD1*NHD2*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR1,NHD1*NHD2*NHD3*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR2,NHD1*NHD2*NHD3*NHD4*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (HFTR3,NHD1*NHD2*NHD3*NHD4*NHD5*NHD6,
     .                MPI_REAL8,0,MPI_COMM_WORLD,ier)

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_CREF

      END MODULE EIRMOD_CREF
