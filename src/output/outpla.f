CDR SEPT.2014:  MODIFIED TO PRINT INPUT TALLIES ON UNDERLYING STRUCTURED GRID,
C               SAME AS THAT USED FOR OUTPUT TALLIES.
C               ASSUME: FINE GRID   DEFINED BY STRUCTURE NR1ST, NP2ND, NT3RD, NRADD, NSBOX
C                       COARSE GRID DEFINED BY STRUCTURE NR1TAL,NP2TAL,NT3TAL,NRADD_TAL, NSBOX_TAL
C                       MAPPING PROVIDED BY I_COARSE=NCLTAL(I_FINE)
C                       ADDITIONAL CELL REGION: NSURF+1,...    ,NSBOX (NRADD CELLS) ON FINE GRID
C                                               NSURF_TAL+1,...,NSBOX_TAL  ON UNDERLYING COARSE GRID
C                       THE GRIDS OVERLAP IN CELLS 1,....,NSURF_TAL OF THE COARSE GRID.
C               TO BE DONE: SWITCH BETWEEN OLD AND NEW OPTION.
cdr  Aug. 16: 1D grid set for printout on separate tally output streams in 1D cases
cdr  Jan. 17:  switching between grids (coarse, structured and finer, unstrucutred)
cdr            new NFLAGV option: >=0: old (default), use fine grid
cdr                               < 0: new: coarse grain by averaging onto
cdr                                         coarser structured grid,
cdr                                         and use abs(nflagv) as before nflagv.
cdr to be done: loops 121 and 122 are identical, once i_fine is set. eliminate one of them?
cdr             inttal and intvol are largely identical, remove one ?
cdr             prttal and prtvol are largely identical, remove one ?
cdr oct 18    : all input tallies selectable, also derived tallies.
cdr             also: gradient tallies of input tallies: currently no. 31--120
cdr may 19    : remove NF=NFRSTP(ITAL) (unused, meaning ?), comments...
cdr jan 22    : bug fix: missing weighting function for PMOM_VEC tally added
cdr             as long as fine-grid data are identical to coarse grid
cdr             data, as e.g. in B2.5 interfaces, this bug made no difference.
cdr mar 23    : account for bgk virt. species as if they were "densmodel" species
C
      SUBROUTINE EIRENE_OUTPLA(ICAL)
C  This routine prints background tallies as requested in input block 11.
C  ical=0: called after initialization phase, prior to Monte Carlo loop.
C  ical=1: called in postprocessing phase after internal iteration loop
C          Some background parameters
C          may have been modified due to iteration loop.
C          Print only the modified tallies.
C
C  PRINT INPUT TALLIES ONTO OUTPUT FILE IUNOUT
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS, ONLY: npbgkp
      USE EIRMOD_CSPEI
      USE EIRMOD_CINIT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICAL
      REAL(DP), ALLOCATABLE :: HELPP(:),HELPW(:),HELPS(:),X1D(:)
      INTEGER,  ALLOCATABLE :: NCLTPR(:)
      REAL(DP) :: TALTYP(NTALI)
      REAL(DP) :: TALAV, HELPI, TALTOT, TOTALW
      INTEGER :: IR, IP, IT, I, I_FINE, NBLCKA, IB, IPRV, ITAL,
     .           NXM, NYM, NZM, NR1PR, NP2PR, NT3PR, NSBPR, NFLGPR,
     .           ITALI, K, NFTI, NFTE, KG,
     .           KS, KK, ICO   ! indexing in vectorial tallies
      LOGICAL :: LPPOST  ! POST PROCESSED INPUT TALLY ?
      EXTERNAL :: EIRENE_INTTAL, EIRENE_INTVOL, EIRENE_PRTTAL,
     .            EIRENE_PRTVOL, EIRENE_HEADNG, EIRENE_LEER,
     .            EIRENE_MASAGE, EIRENE_MASR1, EIRENE_MASR2

C  INDICATOR FOR THE TALLIES THAT MAY HAVE BEEN MODIFIED IN POSTPROCESSING
C  CURRENTLY: BULK ION TEMP (-2), BULK ION DENSITY (-4), AND BULK ION DRIFT VELOCITY (-5,-6,-7)
      INTEGER :: JPRTAL(5) = (/-2,-4,-5,-6,-7/)
C
cdr: extensive or intensive quantities? Needed for averaging....
C  TYPE OF TALLY: TALTYP=0: #              (#-UNITS)
C                 TALTYP=1: # DENSITY      (#-UNITS/CM**3)
C                 TALTYP=2: VOLUME         (CM**3)
C                 TALTYP=3: DIMENSIONLESS  (1)
C                 TALTYP=4: UNKNOWN        (?)
      TALTYP(1)=0
      TALTYP(2)=0
      TALTYP(3)=1
      TALTYP(4)=1
      TALTYP(5)=0
      TALTYP(6)=0
      TALTYP(7)=0
      TALTYP(8)=3
      TALTYP(9)=3
      TALTYP(10)=3
      TALTYP(11)=3
      TALTYP(12)=4
      TALTYP(13)=0
      TALTYP(14)=2   ! CELL VOLUME
      TALTYP(15)=3
      TALTYP(16)=3
      TALTYP(17)=3
      TALTYP(18)=0
      TALTYP(19)=0
      TALTYP(20)=0
      TALTYP(21)=0
      TALTYP(22)=0
      TALTYP(23)=0  ! bv_vec, units: cm/s
      TALTYP(24)=0  ! pmom_vec, units: g*cm/s
      TALTYP(25)=0  ! psi, units: Tesla*m
      TALTYP(26)=3  ! zi

      TALTYP(27)=0  ! free27 units ??
      TALTYP(28)=0  ! free28 units ??
      TALTYP(29)=0  ! free29 units ??
      TALTYP(30)=0  ! free30 units ??

cdr to be done: weighting function for gradient tallies. Tentatively set =0
      TALTYP(31:NTALI)=0

      IF (ICAL == 1) THEN

cdr "density model" is jargon for setting field particle parameters di, ti and vi_vec
cdr from other field particle parameters plus some further data, e.g. certain rate coefficients.

!  IS ANY "DENSITY MODEL" DEFINED IN THIS RUN AT ALL?
        IF (ALL(CDENMODEL == REPEAT(' ',LEN(CDENMODEL)))
cdr  tbd: .and. all: npbgkp(ipls,1) .eq.0     also BGK virt. species should be printed.
cdr   since BGK field species should be just special cases of "cdenmodel" cases
     .     ) RETURN

!  IS OUTPUT OF SUCH AN "DENSITY-MODEL" INPUT TALLY ASKED FOR?
        DO IPRV=1,NVOLPR
          ITAL=NPRTLV(IPRV)
cdr further below we check, if the printout species range NSPEZI, NSPEZE
cdr includes any IPLS with a "density model"
          IF (ANY(JPRTAL == ITAL)) THEN
            CALL EIRENE_LEER (2)
            CALL EIRENE_HEADNG
     .        ('BACKGROUND TALLIES CHANGED IN POSTPROCESSING',44)
            EXIT
          END IF
        END DO
!  NO OUTPUT OF INPUT TALLIES REQUIRED
        IF (IPRV > NVOLPR) RETURN
      END IF

!  IF NVOLPR <= 0 NOTHING TO BE DONE
      IF (NVOLPR <= 0) RETURN

      ALLOCATE (HELPP(NRAD))
      ALLOCATE (HELPW(NRAD))
      ALLOCATE (HELPS(NRAD))
      ALLOCATE (NCLTPR(NRAD))
      ALLOCATE (X1D(NRAD))

!  1D STRUCTURED GRID FOR PRINTOUT
      X1D=0.D0
      IF (LEVGEO.LE.3.AND.NP2ND.EQ.1.AND.NT3RD.EQ.1) THEN
        DO I=1,NR1ST
          X1D(I)=RHOZNE(I)
        ENDDO
      ENDIF

C

C
C  PRINT THOSE INPUT VOLUME-AVERAGED TALLIES, WHICH HAVE BEEN SELECTED

C  PUT TALLY ITAL ONTO ARRAY HELPP, THIS IS RESOLVED WRT. FINEST GRID
C
      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)
      DO 100 IPRV=1,NVOLPR
        ITAL=NPRTLV(IPRV)
        IF ((ICAL == 1) .AND. (ALL(JPRTAL .NE. ITAL))) CYCLE
c  negative tally numbers ital: background (input) tallies, printed here
c  positive tally numbers ital: output tallies, printed from OUTEIR.
        IF (ITAL.LT.0) THEN
          ITALI=-ITAL
!  TALLY SWITCHED OFF ?
          IF (.NOT.LIVTALI(ITALI)) THEN
            WRITE (iunout,*) ' TALLY NOT AVAILABLE (OUTPLA)',
     .                       ' ITAL = ', ITAL
            CALL EIRENE_LEER(1)
            CYCLE
          END IF

          NFTI=1
cdr  be careful with indirect species index addressing here.
cdr  NFTE=NPLSI for Ti and Vin_VEC profiles.
cdr            not: NFTE = NPLSTI or NPLSV, resp.
cdr  due to different density weighting.
cdr  Full physical species range is needed even if
cdr  "intensive" tallies Ti or Vin_VEC are made identical
cdr  for different IPLS, e.g. even IF NFRSTP(ITALI)=1

          NFTE=NFSTPI(ITALI)
          IF (NSPEZV(IPRV,1).GT.0) THEN
c  print tally only for selected range of species indices
            NFTI=NSPEZV(IPRV,1)
            NFTE=MAX(NFTI,NSPEZV(IPRV,2))
          ENDIF
c
          DO 119 K=NFTI,NFTE
c  check for valid range of tally ITALI
            IF (K.GT.NFSTPI(ITALI)) THEN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*)
     .          'SPECIES INDEX OUT OF RANGE IN SUBR. OUTPLA'
              WRITE (iunout,*) 'ITALI, K, ',ITALI,K
              CALL EIRENE_LEER(1)
              GOTO 119
            ENDIF

c  HELPP, HELPW: for weighting when averaging tallies across several cells
cdr  if ical=1: print only post processed field particles
cdr  print only species ipls=k with either "cdenmodel", or npbgkp(k,1) .ne.0
            lppost=ical .eq.1 .and. (VERIFY(CDENMODEL(K),' ') == 0) 
     .             .and. (npbgkp(k,1) == 0)

            SELECT CASE (ITALI)
            CASE (1)
              HELPP(1:NSBOX) = TEIN(1:NSBOX)
            CASE (2)
cdr was missing here: verify cdenmodel(k) ?
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = TIIN(MPLSTI(K),1:NSBOX)
            CASE (3)
              HELPP(1:NSBOX) = DEIN(1:NSBOX)
            CASE (4)
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = DIIN(K,1:NSBOX)
            CASE (5)
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = VXIN(MPLSV(K),1:NSBOX)
            CASE (6)
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = VYIN(MPLSV(K),1:NSBOX)
            CASE (7)
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = VZIN(MPLSV(K),1:NSBOX)
            CASE (8)
              HELPP(1:NSBOX) = BXIN(1:NSBOX)
            CASE (9)
              HELPP(1:NSBOX) = BYIN(1:NSBOX)
            CASE (10)
              HELPP(1:NSBOX) = BZIN(1:NSBOX)
            CASE (11)
              HELPP(1:NSBOX) = BFIN(1:NSBOX)
            CASE (12)
              HELPP(1:NSBOX) = ADIN(K,1:NSBOX)
            CASE (13)
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = EDRIFT(K,1:NSBOX)
            CASE (14)
              HELPP(1:NSBOX) = VOL(1:NSBOX)
            CASE (15)
              HELPP(1:NSBOX) = WGHT(K,1:NSBOX)
            CASE (16)
              HELPP(1:NSBOX) = BXPERP(1:NSBOX)
            CASE (17)
              HELPP(1:NSBOX) = BYPERP(1:NSBOX)
            CASE (18)
              HELPP(1:NSBOX) = EXIN(1:NSBOX)
            CASE (19)
              HELPP(1:NSBOX) = EYIN(1:NSBOX)
            CASE (20)
              HELPP(1:NSBOX) = EZIN(1:NSBOX)
            CASE (21)
              HELPP(1:NSBOX) = EFIN(1:NSBOX)
            CASE (22)
              HELPP(1:NSBOX) = POT(1:NSBOX)
cdr  next two tallies are vectorial tallies
            CASE (23)
cdr K ranges from 1 to 3*NPLS, KS: 1:NPLS
cdr BV_VEC(1:3*NPLSV)
              KS=MOD(K,NPLS)
              IF (KS == 0 ) KS=NPLS
              ICO=(K-1)/NPLS  ! = 0,1,2
              KK=ICO*NPLSV
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = BV_VEC(KK+MPLSV(KS),1:NSBOX)
            CASE (24)
cdr  K ranges from 1 to 3*NPLS, KS: 1:NPLS
cdr  PMOM_VEC(1:3*NPLS)
              KS=MOD(K,NPLS)
              IF (KS == 0 ) KS=NPLS
              IF (lppost) CYCLE
              HELPP(1:NSBOX) = PMOM_VEC(K,1:NSBOX)
            CASE (25)
              HELPP(1:NSBOX) = PSI(1:NSBOX)
            CASE (26)
              HELPP(1:NSBOX) = ZIIN(K,1:NSBOX)
            CASE (27)
              HELPP(1:NSBOX) = FREE27(1:NSBOX)
            CASE (28)
              HELPP(1:NSBOX) = FREE28(1:NSBOX)
            CASE (29)
              HELPP(1:NSBOX) = FREE29(1:NSBOX)
            CASE (30)
              HELPP(1:NSBOX) = FREE30(1:NSBOX)
c
            CASE (31:120)  ! ntali=120, constant required here
!  GRADIENTS
              KG = NADDP(ITALI)+K
              HELPP(1:NSBOX) = PLSTLS(KG,1:NSBOX)
            CASE DEFAULT
              WRITE (iunout,*)
     .          ' WRONG TALLY NUMBER (OUTPLA), ITAL = ',ITAL
              WRITE (iunout,*) ' NO OUTPUT PERFORMED '
              CALL EIRENE_LEER(1)
              GOTO 100
            END SELECT

C  IN CASE OF TWO GRIDS, (COARSE-GRAINING)
C  SWITCH BETWEEN PRINTOUT OF TALLY ITALI EITHER ON UNDERLYING FINE GRID (OLD DEFAULT)
C  OR OF THE AVERAGED (COARSE GRAINED) TALLY ON COARSE GRID
C  IF NCLTAL(I_FINE)==I_COARSE, EVERYWHERE, THEN THERE IS ONLY ONE GRID

cdr  Apr.22  coarse graining may be incorrect for Ti, VXIN,... and BVIN,
cdr          i.e. for tallies with indirect species index addressing
cdr  Apr.24  If no coarser grid was defined, e.g. in interfacing routines,
cdr          then the next two options should be identical.
            if (nflagv(iprv).lt.0) then
c  before printing: reset (by weighted averaging) input tallies
c                   onto coarser (structured) grid
              nr1pr  =nr1tal
              np2pr  =np2tal
              nt3pr  =nt3tal
              nsbpr  =nsbox_tal
              nflgpr =-nflagv(iprv)
              do i_fine=1,nsbox
                ncltpr(i_fine) = ncltal(i_fine)
              enddo
            else
c  no coarse graining of input tallies (default):
c  print on underlying fine grid. Also if there is only one single grid
              nr1pr  =nr1st
              np2pr  =np2nd
              nt3pr  =nt3rd
              nsbpr  =nsbox
              nflgpr =nflagv(iprv)
              do i_fine=1,nsbox
                ncltpr(i_fine) = i_fine
              enddo
            endif

C
C  SET WEIGHTING FUNCTION HELPW FOR INPUT TALLY ITALI
C  CURRENTLY INPUT TALLIES ARE GIVEN ON FINE "GEOMETRY MESH", STRUCTURE NR1ST,NP2ND,....
C  USE NCLTPR(I) = NCELL ARRAY TO COARSE GRAIN ONTO STRUCTURED GRID.
C  The same cell indexing is also used to score output tallies directly
C  on coarse grained grid.
C
            TOTALW=0.D0
            HELPW = 0.D0
            HELPS = HELPP ! ORIGINAL INPUT TALLY ON FINE GRID:
                          ! MOVED TO HELPS
            HELPP = 0.D0  ! BUILD NEW TALLY NOW, ON UNDERLYING
                          ! COARSER GRID, IF WANTED
c                         ! otherwise: helpp=helps
            DO 121 IB=1,NBMLT
             NBLCKA=NSTRD*(IB-1)
             DO IR=1,NXM
             DO IP=1,NYM
             DO IT=1,NZM
              I_FINE=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
C  COARSE-GRAINING OF INPUT TALLY ITAL ONTO GRID DEFINED BY NCLTPR(I-FINE),
C  STRUCTURE NR1TAL,NP2TAL,...
C  IF THERE IS ONLY ONE SINGLE GRID, THEN NCLTPR(I)==I,
C  AND NO COARSE-GRAINING IS DONE
              I=NCLTPR(I_FINE)
              SELECT CASE (ITALI)

cdr  we should use TALTYP(ITALI) for the weights HELP,HELPW,
cdr  rather than explicit coding here for each tally.

              CASE (1)
C  1: ELECTR. TEMPERATURE: NE*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+HELPS(I_FINE)*DEIN(I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DEIN(I_FINE)*VOL(I_FINE)
              CASE (2)
C  2: ION TEMPERATURE: NI(K)*VOLUME-WEIGHTED AVERAGES
cdr  here K also must range from 1 to NPLSI,  not only to NPLSTI
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (3:4)
C  3,4: PARTICLE DENSITY PROFILES: VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+HELPS(I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+VOL(I_FINE)
              CASE (5:7)
C  5,6,7: ION DRIFT VELOCITY: NI(K)*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (8:11)
C  8,9,10,11: B FIELD UNIT VECTOR, B FIELD STRENGTH "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (16:17)
C  16,17: B_PERP FIELD: "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (14)
C  14: CELL VOLUME  = UNWEIGHTED SUM
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=1.D0
              CASE (12,15)
C  12: ADDITIONAL TALLY (NO.12)
C  15: WEIGHT WINDOW FUNCTION  (NO.15)
C  "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (13)
C  13: ION DRIFT ENERGY: NI(K)*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (18:21)
C  18,19,20,21: E FIELD UNIT VECTOR, E FIELD STRENGTH
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (22)
C  22: (ELECTRIC) POTENTIAL
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (23)
C  23: FLOW VELOCITY WRT. B
                KS=MOD(K,NPLSI)
                IF (KS == 0 ) KS=NPLSI
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(KS,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(KS,I_FINE)*VOL(I_FINE)
              CASE (24)
C  24: MOMENTUM FLOW VECTOR, : NI(K)*VOLUME-WEIGHTED AVERAGES
cdr for tally 24 K ranges from 1 to 3*NPLS
cdr Jan.22:  weighting was missing. Now corrected for this tally
                KS=MOD(K,NPLSI)
                IF (KS == 0 ) KS=NPLSI
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(KS,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(KS,I_FINE)*VOL(I_FINE)
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (25)
C  25) PSI
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (26)
C  26) ZI
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (27)
C  27) FREE27
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (28)
C  28) FREE28
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (29)
C  29) FREE29
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (30)
C  30) FREE30
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0

              CASE (31:120)  ! ntali=120, constant required here
C  (31 .. NTALI) GRADIENTS
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              END SELECT
             END DO
             END DO
             END DO
  121       CONTINUE
C
C  SAME LOOP AGAIN, (IDENTICAL CODE INSIDE LOOP) OVER ADDITIONAL CELL REGION
            DO 122 I_FINE=NSURF+1,NSURF+NRADD
C             I_FINE=I_FINE
C  COARSE-GRAINING OF INPUT TALLY ITAL ONTO GRID DEFINED BY NCLTPR(I-FINE),
C  STRUCTURE NR1TAL,NP2TAL,....
C  WHEN THERE IS ONLY ONE SINGLE GRID, THEN NCLTPR(I)==I, AND NO COARSE-GRAINING IS DONE
              I=NCLTPR(I_FINE)
              SELECT CASE (ITALI)
              CASE (1)
C  ELECTR. TEMPERATURE: NE*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+HELPS(I_FINE)*DEIN(I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DEIN(I_FINE)*VOL(I_FINE)
              CASE (2)
C  ION TEMPERATURE: NI(K)*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (3:4)
C  PARTICLE DENSITY PROFILES: VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+HELPS(I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+VOL(I_FINE)
              CASE (5:7)
C  ION DRIFT VELOCITY: NI(K)*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (8:11)
C  B FIELD UNIT VECTOR, B FIELD STRENGTH "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (16:17)
C  B_PERP FIELD: "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (14)
C  CELL VOLUME  = UNWEIGHTED SUM
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=1.D0
              CASE (12,15)
C  ADDITIONAL TALLY (NO.12)
C  WEIGHT WINDOW FUNCTION  (NO.15)
C  "1 - WEIGHTED" AVERAGES, = ARITHM. MEAN
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (13)
C  ION DRIFT ENERGY: NI(K)*VOLUME-WEIGHTED AVERAGES
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(K,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(K,I_FINE)*VOL(I_FINE)
              CASE (18:21)
C  E FIELD UNIT VECTOR, E FIELD STRENGTH
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (22)
C  (ELECTRIC) POTENTIAL
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (23)
C  FLOW VELOCITY VECTOR wrt B FIELD
cdr K ranges from 1 to 3*NPLS
                KS=MOD(K,NPLS)
                IF (KS == 0 ) KS=NPLS
cdr weighting is the same for all 3 vector components of BV_VEC
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(KS,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(KS,I_FINE)*VOL(I_FINE)
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (24)
C  MOMENTUM FLOW VECTOR wrt. B field
cdr for tally 24 K ranges from 1 to 3*NPLS
cdr Jan.22:  weighting was missing, Now corrected for this tally
                KS=MOD(K,NPLS)
                IF (KS == 0 ) KS=NPLS
                HELPP(I)=HELPP(I)+
     .                   HELPS(I_FINE)*DIIN(KS,I_FINE)*VOL(I_FINE)
                HELPW(I)=HELPW(I)+DIIN(KS,I_FINE)*VOL(I_FINE)
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (25)
C  PSI
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (26)
C  ZI
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (27)
C  27) FREE27
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (28)
C  28) FREE28
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (29)
C  29) FREE29
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (30)
C  30) FREE30
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              CASE (31:120)  ! ntali=120, constant required here
C  GRADIENTS
                HELPP(I)=HELPP(I)+HELPS(I_FINE)
                HELPW(I)=HELPW(I)+1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              END SELECT
  122       CONTINUE

C  COARSE-GRAINING: NORMALIZE WEIGHTED SUMS BY THEIR WEIGHT
C  TOTAL OF WEIGHT, FOR TOTAL TALLY AVERAGE TALAV
            DO I=1,NSBPR
              TOTALW=TOTALW+HELPW(I)
              HELPP(I)=HELPP(I)/(HELPW(I)+EPS60)
            ENDDO

C  INTEGRATE AVERAGES ALONG COORDINATES, USE WEIGHTING HELPW
            IF (ITALI.NE.NTALO) THEN
              CALL EIRENE_INTTAL (HELPP,HELPW,1,1,NSBPR,HELPI,
     .                            NR1PR,NP2PR,NT3PR,NBMLT)
C  SAME FOR TALLY NTALO: CELL VOLUME
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL EIRENE_INTVOL (HELPP,      1,1,NSBPR,HELPI,
     .                            NR1PR,NP2PR,NT3PR,NBMLT)
            ENDIF
            TALTOT=HELPI
            TALAV=HELPI/(TOTALW+EPS60)
C
            IF (TALTOT.EQ.0.D0) GOTO 118

C   PRINT TALLY.
C   IN CASE NFLAGV<0: USING UNDERLYING COARSER GRID
            IF (ITALI.NE.NTALO) THEN
              CALL EIRENE_PRTTAL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,X1D,
     .                    NR1PR,NP2PR,NT3PR,NBMLT,NSBPR,
     .                    NFLGPR,NTLVFL(IPRV))
C   SPECIAL TREATMENT FOR PRINTING CELL VOLUMES (INPUT TALLY ITALI=NTALO=14)
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL EIRENE_PRTVOL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,X1D,
     .                    NR1PR,NP2PR,NT3PR,NBMLT,NSBPR,
     .                    NFLGPR,NTLVFL(IPRV))
            ENDIF
            CALL EIRENE_LEER(2)
            IF (TALTYP(ITALI).EQ.0) THEN
              CALL EIRENE_MASAGE
     .        ('WEIGHTED MEAN VALUE ("UNITS")               ')
              CALL EIRENE_MASR1 ('MEAN:   ',TALAV)
              CALL EIRENE_LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.1) THEN
              CALL EIRENE_MASAGE
     .        ('WEIGHTED TOTAL ("UNITS*CM**3"), MEAN ("UNITS")')
              CALL EIRENE_MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL EIRENE_LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.2) THEN
              CALL EIRENE_MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS")                         ')
              CALL EIRENE_MASR1 ('TOTAL:  ',TALTOT)
              CALL EIRENE_LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.3) THEN
              CALL EIRENE_MASAGE
     .        ('ARITHMETIC MEAN ("UNITS")                          ')
              CALL EIRENE_MASR1 ('MEAN:   ',TALAV)
              CALL EIRENE_LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.4) THEN
              CALL EIRENE_MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS"), MEAN ("UNITS")')
              CALL EIRENE_MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL EIRENE_LEER(3)
            ENDIF
            GOTO 119
C
  118       CONTINUE
C   PRINT ONLY THE HEADER FOR TALLY, BECAUSE TALLY IDENTICALLY ZERO
            CALL EIRENE_PRTTAL
     .                     (TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                      TXTPUN(K,ITALI),
     .                      HELPP,X1D,
     .                      NR1PR,NP2PR,NT3PR,NBMLT,NSBPR,-1,0)
            CALL EIRENE_MASAGE
     .              ('IDENTICALLY ZERO, NOT PRINTED')
            CALL EIRENE_LEER(2)
  119     CONTINUE
        ENDIF
  100 CONTINUE
      CALL EIRENE_LEER(2)

      DEALLOCATE (HELPP)
      DEALLOCATE (HELPW)
      DEALLOCATE (HELPS)
      DEALLOCATE (NCLTPR)
      DEALLOCATE (X1D)
C
      RETURN
      END SUBROUTINE EIRENE_OUTPLA
