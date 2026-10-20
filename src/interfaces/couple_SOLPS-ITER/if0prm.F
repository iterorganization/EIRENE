cdr called from FIND_PARAM.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)

C
cdr First pass reading from block 14
cdr Later this block 14 will be read again from interfacing routines (e.g. infcop.f),
cdr not from eirene_input.f, which only reads blocks 1 to 13.
c   this version: couple_SOLPS-ITER.
c
c   for particle balance rescaling:
c   define array: LKINDL(NPLS): assign a chemical element to the various plasma ions.
c                 e.g.:  He+ and He++ should be the same "kind",
c                        or C+,C++,...C6+ be "another kind".
c
c   and set storage for allocatable arrays:
c   NSTEP :
c   NPTRGT:
c   NAIN  :
c   NAOT  :
c   NCPV  :
c   NKNOT :
C   NTRII :
C   NCPVI : no. of special couple tallies

c  also set: NDX,NDY,NFL,NDXP,NDYP

      SUBROUTINE EIRENE_IF0PRM(IUNIN)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CINIT
      USE EIRMOD_BRAEIR
      USE EIRMOD_CLOGAU
     , , ONLY: NLTRIMESH, NLSPCSCL, NLSPCSCL_ON
      USE EIRMOD_CCOUPL
     , , ONLY: NTGPRT, NSPZE, NSPZI, EIRENE_ALLOC_CCOUPL
      use eirmod_extrab25

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IUNIN
      INTEGER :: NFLA, NCUTB, NCUTL, NDXA, NDYA, IPL, NTARGI, IT, IPRT,
     .           NAINB, IAIN, NAOTB, IAOT, NRKNOT, ITRG,
     .           NTRII
      INTEGER, SAVE :: IO
 ! Extend to read as written by uinp (jdl)
      INTEGER :: IDUMMY(0:99)=0
      REAL(DP) :: RDUMMY(0:9), RLAST
      LOGICAL :: LDUMMY(0:5), LSPRCL_LOC
      CHARACTER(72) :: ZEILE
      EXTERNAL :: FIX_LOGICAL_INPUT

C  READ INPUT BLOCK 14
      READ (IUNIN,'(A72)') ZEILE
      call fix_logical_input(zeile,6)
      READ (ZEILE,'(5L1,1X,5L1)') LDUMMY(0:5)
      IF (.NOT.NLTRIMESH) NLTRIMESH = LDUMMY(3)
      LSPRCL_LOC = LDUMMY(4)
      READ (IUNIN,'(9I6)') IDUMMY(0:8)
      IF (IDUMMY(0).EQ.0 .OR. IDUMMY(1).LT.0) THEN
! This was a junk Eirene_96 line
        READ (IUNIN,'(3I6)') NFLA,NCUTB,NCUTL
      ELSE
        NFLA  = IDUMMY(0)
        NCUTB = IDUMMY(1)
        NCUTL = IDUMMY(2)
      ENDIF

      IF(.NOT.ALLOCATED(LKINDP)) ALLOCATE(LKINDP(NPLS))
      LKINDP=0
      COMUSR_FIRST_PASS(5) = .FALSE.
      DO IPL=1,NPLS
        READ (IUNIN,'(2I6,2E12.4,I6)') IDUMMY(3:4),RDUMMY(0:1),IDUMMY(5)
        IF (IDUMMY(4).GT.0) THEN
          IF (IDUMMY(5).GT.0) THEN
            LKINDP(IPL)=IDUMMY(5)
          ELSE
            IF (IPL.EQ.1) THEN
              LKINDP(IPL)=1
              RLAST=RDUMMY(1)
            ELSE
              IF (RDUMMY(1).EQ.RLAST) THEN
                LKINDP(IPL)=LKINDP(IPL-1)
              ELSE
                LKINDP(IPL)=LKINDP(IPL-1)+1
              END IF
              RLAST=RDUMMY(1)
            END IF
          END IF
        END IF
      END DO

C  GRID SIZE IN 2D PLASMA FLUID CODE
      READ (IUNIN,'(2I6)') NDXA,NDYA
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      READ (IUNIN,'(I6)') NTARGI
C  NUMBER OF PARTS PER TARGET RECYCLING SOURCE
      IF (NTARGI.GT.0) THEN
        NSTEP=MAX(NSTEP,NTARGI)
        IF (.NOT.(IDUMMY(0).EQ.0 .OR. IDUMMY(1).LT.0)) THEN
          READ (IUNIN,'(12I6)') (IDUMMY(IT),IT=1,NTARGI)
          NPTRGT=SUM(IDUMMY(1:NTARGI))
        ELSE
          NPTRGT=NTARGI
        ENDIF
        CALL EIRENE_ALLOC_CCOUPL(1)
        IF (IDUMMY(0).EQ.0 .OR. IDUMMY(1).LT.0) THEN
          NTGPRT(1:NTARGI) = 1 ! Not provided in Eirene_96 input
        ELSE
          DO IT = 1, NTARGI
            NTGPRT(IT) = IDUMMY(IT)
          END DO
        ENDIF
        ITRG=0
        DO IT=1,NTARGI
          DO IPRT=1,NTGPRT(IT)
            ITRG=ITRG+1
  331       READ (IUNIN,'(A72)') ZEILE
            IF (ZEILE(1:1).EQ.'*') THEN
              GOTO 331
            ENDIF
            READ (ZEILE,'(12I6)') IDUMMY(0:11)
            IDUMMY(9) = MAX0(1,IDUMMY(9)) ! NSPZI
            IDUMMY(10) = MIN0(NFLA,IDUMMY(10)) ! NSPZE
            IF (IDUMMY(10).LT.IDUMMY(9)) THEN
              IDUMMY(9)=1
              IDUMMY(10)=NFLA
            ENDIF
            NSPZI(IT,IPRT)=IDUMMY(9)
            NSPZE(IT,IPRT)=IDUMMY(10)
          END DO
        END DO
#ifdef B25_EIRENE
        IF (.NOT.ALLOCATED(RCPOS_EIR)) THEN
          ALLOCATE(RCPOS_EIR(ITRG),RCBEG_EIR(ITRG),RCEND_EIR(ITRG))
          ALLOCATE(RCPRT_EIR(ITRG),RCSPI_EIR(ITRG),RCSPE_EIR(ITRG))
          ALLOCATE(RCCHR_EIR(ITRG))
        END IF
#endif
      END IF
      READ (IUNIN,*)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      READ (IUNIN,'(I6)') NAINB
      NAIN = MAX(NAIN,NAINB)
      DO IAIN=1,NAINB
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
      END DO
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      READ (IUNIN,'(I6)') NAOTB
      DO IAOT=1,NAOTB
        READ (IUNIN,*)
      END DO
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                          AND IN SUBR. COLLIDE)
      NCPVI= NPLS
      NCPV = MAX(NCPV,NCPVI)
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      NDX = NDXA
      NDY = NDYA
      NFL = NFLA
      IF (ALLOCATED(DNIB)) THEN
        NDX=UBOUND(DNIB,1)-1
        NDY=UBOUND(DNIB,2)-1
        NFL=UBOUND(DNIB,3)
      END IF
      NDXP = NDX+1
      NDYP = NDY+1
C
C
C  TRANSFER GEOMETRY
C
C  READ DATA FOR TRIANGULAR MESH
C
      NRKNOT=0
      NTRII=0
      if(.not.nltrimesh) then !{
C
        OPEN (UNIT=33,ACCESS='SEQUENTIAL',FORM='FORMATTED')
        READ(33,*,IOSTAT=IO) NRKNOT
        IF(IO.NE.0) NRKNOT=0
        WRITE(IUNOUT,'(a14,i8)') 'NRKNOT      = ',NRKNOT
        CLOSE (UNIT=33)
C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
        OPEN (UNIT=34,ACCESS='SEQUENTIAL',FORM='FORMATTED')
        READ(34,*,IOSTAT=IO) NTRII
        IF(IO.NE.0) NTRII=0
        WRITE(IUNOUT,'(a14,i8)') 'NTRII       = ',NTRII
        CLOSE (UNIT=34)
      end if !}

      NKNOT=NRKNOT
      NTRI=NTRII+1
csw 09jan2012
      IF(NOPTIM.LT.0) NOPTIM=(NTRII+1)*N3RD !VK no storage optimization
csw
cxpb Define some numbers needed by eirmod_extrab25
!cank 960623
      nnplsi=nfla
      nns=nnplsi

C  TRANSFER FLAG LSPRCL WHICH IS KNOWN ONLY IN THIS INTERFACE TO
C  GLOBALLY KNOWN VARIABLE NLSPCSCL
      NLSPCSCL = NLSPCSCL .OR. LSPRCL_LOC
      NLSPCSCL_ON = NLSPCSCL_ON .OR. LSPRCL_LOC

      RETURN
      END SUBROUTINE EIRENE_IF0PRM
