cdr called from find_param.f in initialization phase,
cdr when eirene is in "coupled mode": i.e. IF(NMODE.NE.0)
C
cdr Read block 14 from interfacing routines (not from eirene_input.f)
c   this version: couple_dummy, i.e. only dummy interfacing routines.
c
c   and set storage for allocatable arrays:
c   NSTEP :
c   NPTRGT:
c   NAIN  :
c   NAOT  :
c   NCPV  :
c   NKNOT :
C   NTRII :
C   NCPVI :  no. of special couple tallies

c  also set: NDX,NDY,NFL, NDXP, NDYP

#ifndef NO_JSON
      SUBROUTINE EIRENE_IF0PRM_JSON(json,p)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CINIT
      USE EIRMOD_COMPRT, only: iunout
      USE EIRMOD_BRAEIR
      USE EIRMOD_CLOGAU
     , , ONLY: NLTRIMESH, NLSPCSCL, NLSPCSCL_ON
      USE EIRMOD_CCOUPL
     , , ONLY: NTGPRT, NSPZE, NSPZI, EIRENE_ALLOC_CCOUPL
      use eirmod_extrab25
      use json_module           !IGNORE
     .    , lk => json_lk, rk => json_rk, ik => json_ik, ck => json_ck

      IMPLICIT NONE

      class(json_core), intent(inout) :: json
      type(json_value), pointer, intent(in) :: p
      type(json_value), pointer :: pflds, pfld, ptrgs, ptrg, prts, prt

      INTEGER :: NFLA, NCUTB, NCUTL, NDXA, NDYA, IPL, NTARGI, IT, IPRT,
     .           NAINB, NAOTB, ITRG,
     .           NTRII, NR1ST, IFL, LKP, NSPI, NSPE
      REAL(DP) :: BM, RLAST
      INTEGER, ALLOCATABLE :: NTG(:)
      LOGICAL :: FOUND, LDUMMY, FOUNDI, FOUNDE, LSPRCL_LOC
      EXTERNAL :: EIRENE_FIND_TRIANG_DIM

      WRITE (iunout,*) '*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"'

C  READ INPUT BLOCK 14
      call json%get(p,'NLTRIMESH',ldummy,found)
      IF (.NOT.NLTRIMESH.AND.FOUND) NLTRIMESH = LDUMMY

      call json%get(p,'LSPRCL',ldummy,found)
      IF (FOUND) LSPRCL_LOC = LDUMMY

      call json%get(p,'NFLA',nfla,found)
      call json%get(p,'NCUTB',ncutb,found)
      call json%get(p,'NCUTL',ncutl,found)

      IF(.NOT.ALLOCATED(LKINDP)) ALLOCATE(LKINDP(NPLS))
      LKINDP=0
      COMUSR_FIRST_PASS(5) = .FALSE.
      call json%get_child(p,'B2FLUIDS',pflds)
      DO IPL=1,NPLS
        call json%get_child(pflds,ipl,pfld)
        call json%get(pfld,'IFLB',ifl,found)
        call json%get(pfld,'BMASS',bm,found)
        call json%get(pfld,'LKINDP',lkp,found)
        if (ifl .gt. 0) then
          if (lkp .gt. 0) then
            lkindp(ipl) = lkp
          else
            if (ipl .eq. 1) then
              lkindp(ipl) = 1
              rlast = bm
            else
              if (bm .eq. rlast) then
                lkindp(ipl) = lkindp(ipl-1)
              else
                lkindp(ipl) = lkindp(ipl-1)+1
              end if
              rlast = bm
            end if
          end if
        end if
        nullify(pfld)
      END DO
      nullify(pflds)

C  GRID SIZE IN 2D PLASMA FLUID CODE
      call json%get(p,'NDXA',ndxa,found)
      call json%get(p,'NDYA',ndya,found)

C     NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      call json%get(p,'NTARGI',ntargi,found)

C     NUMBER OF PARTS PER TARGET RECYCLING SOURCE
      IF (NTARGI.GT.0) THEN
        NSTEP=MAX(NSTEP,NTARGI)
        call json%get(p,'NTGPRT',ntg,found)
        IF (FOUND) THEN
          NPTRGT=SUM(NTG)
        ELSE
          IF (.NOT.ALLOCATED(NTG)) ALLOCATE(NTG(NTARGI))
          NTG(1:NTARGI)=1
          NPTRGT=NTARGI
        END IF
        CALL EIRENE_ALLOC_CCOUPL(1)
        NTGPRT(1:NTARGI)=NTG(1:NTARGI)
        DEALLOCATE (NTG)

        call json%get_child(p,'TARGETS',ptrgs,found)
        ITRG=0
        DO IT=1,NTARGI
          call json%get_child(ptrgs,it,ptrg,found)
          call json%get_child(ptrg,'PARTS',prts,found)
          DO IPRT=1,NTGPRT(IT)
            ITRG=ITRG+1
            call json%get_child(prts,iprt,prt)
            call json%get(prt,'NSPZI',nspi,foundi)
            call json%get(prt,'NSPZE',nspe,founde)
            if (.not.foundi) nspi=1
            if (.not.founde) nspe=nfla
            NSPZI(IT,IPRT)=nspi
            NSPZE(IT,IPRT)=nspe
            nullify(prt)
          END DO
          nullify(prts)
          nullify(ptrg)
        END DO
        nullify(ptrgs)
#ifdef B25_EIRENE
        IF (.NOT.ALLOCATED(RCPOS_EIR)) THEN
          ALLOCATE(RCPOS_EIR(ITRG),RCBEG_EIR(ITRG),RCEND_EIR(ITRG))
          ALLOCATE(RCPRT_EIR(ITRG),RCSPI_EIR(ITRG),RCSPE_EIR(ITRG))
          ALLOCATE(RCCHR_EIR(ITRG))
        END IF
#endif
      END IF

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      call json%get(p,'NAINB',nainb,found)
      NAIN = MAX(NAIN,NAINB)

C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      call json%get(p,'NAOTB',naotb,found)
C
C READING BLOCK 14 FROM JSON INPUT FILE FINISHED
C
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                          AND IN SUBR. COLLIDE)
      NCPVI=NPLS
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
C  READ DATA FOR TRIANGULAR MESH
C
C
      nr1st = n1st
      call eirene_find_triang_dim(nr1st,ntri,ntrii,nknot,ngitt,' ')

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
      END SUBROUTINE EIRENE_IF0PRM_JSON
#else
      SUBROUTINE EIRENE_IF0PRM_JSON
      RETURN
      END SUBROUTINE EIRENE_IF0PRM_JSON
#endif
