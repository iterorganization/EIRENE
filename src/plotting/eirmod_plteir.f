      MODULE EIRMOD_PLTEIR

cdr  30.4.04:  call plttly for spectra corrected.
cdr            first bin (no. 0) and last bin (no. nsts+1) contain
cdr            the fluxes outside the range of spectra.
cpb  30.7.04:  deal with switched off tallies
cdr  10.6.05:  further modifications of plot for spectra (text,
c              total, plot vs. wavelength, plot 2 spectra into same frame)
!pb  18.12.06: general checking of XMCP removed to allow plots of
!              input tallies even if no Monte Carlo particle has been followed
!    10.01.07: SUBROUTINE PLTEIR_REINIT added for reinitialization of EIRENE
!    30.01.09: TEXT CORRECTED FOR PLOTS OF SPECTRA (INPUT BLOCK 10F)
!    30.01.09: additional wavelength unit plots only for photon spectra.
C              Turned off for all other particle types
cdr  Oct.14  : bug fix re. 'l_same', make sure that first spectra plot is on own frame,
cdr            even if other (volumetric) output tallies have already been plotted
cdr            from same stratum in same call to plteir.
cdr  Aug.15  : scaling of spectrum tallies: hard-wired options. To be done !
cdr  Mai 19  : Plotting of spectra incomplete. Tbd: use velocity scale too,
cdr            Allow also plotting for volumetric spectra.
cdr            Careful: log energy scale option, and combination with negative energies?
cdr            Directional spectra
cdr            Sum over species of intensive input tallies is certainly nonsense !
cdr            proper weighting is missing here.
cpb  Oct 22  : I0 and INDX removed
C

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRPTL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLMSK
      USE EIRMOD_CPLOT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMSOU
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_PLTEIR, EIRENE_PLTEIR_REINIT

      INTEGER, SAVE :: IFIRST=0

      CONTAINS

C
      SUBROUTINE EIRENE_PLTEIR (ISTRA)
C
C  ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA
C  PLOT PLASMA TALLIES ONLY ONCE, BUT OUTPUT TALLIES FOR ALL STRATA AS REQUESTED.
C
C
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA

      REAL(DP), ALLOCATABLE :: VECTOR(:,:), VECSAV(:,:), VSDVI(:,:)
      REAL(DP), ALLOCATABLE :: XSPEC(:), YSPEC(:,:), VSPEC(:,:),
     .          WLSPEC(:), YSPECWL(:,:), VSPECWL(:,:)
      REAL(DP), ALLOCATABLE :: DUMMY(:)
      REAL(DP) :: XXP3D_DUM(1), YYP3D_DUM(1)
      REAL(DP), ALLOCATABLE :: YMN2(:), YMX2(:), YMNLG2(:), YMXLG2(:)
      REAL(DP) :: XMI, XMA, TMIN, TMAX, XI, XE, DEL, OUTAUI,
     .            SPCAN, SPC00, WL00, DE, DW
      INTEGER :: IR1(NPLT), IR2(NPLT), IRS(NPLT)
      INTEGER :: IXXE, IXXI, IYYE, IYYI, K, ISPC, NSPS, INULL,
     .           NF, NFT, I, IA, N, IXSET2, ISPZ, IALG, N1SDVI, ISAVE,
     .           IALV, ITL, JTAL, IBLD, ICURV, IE, IXSET3, IS,
     .           IERR, ICINC, IYSET3, IX, I2M, J, IRAD, I1, I2, IT,
     .           ITT, ITP,
     .           KS, KK, ICO, JJTAL
      LOGICAL :: LPLOT2(NPLT), LSDVI(NPLT), LINLOG, L_SAME
      CHARACTER(24) :: TXUNIT(NPLT), TXSPEC(NPLT)
      CHARACTER(24) :: TXUNT1, TXSPC1
      CHARACTER(72) :: TXTALL(NPLT)
      CHARACTER(72) :: TXTLL1
      CHARACTER(72) :: HEAD,  HEAD0, HEAD1, HEAD2, HEAD3, HEAD4,
     .                 HEAD5, HEAD6, HEAD7, HEAD8, HEAD9, HEAD10, TXHEAD
      EXTERNAL :: EIRENE_ALGTAL, EIRENE_INTTAL, EIRENE_ISOLNE,
     .            EIRENE_PL3DPG, EIRENE_PLOT3D, EIRENE_PLTTLY,
     .            EIRENE_RPSCOL, EIRENE_RPSVEC, EIRENE_RSTRT,
     .            EIRENE_SYMET, EIRENE_VECLNE,
     .            EIRENE_LEER, EIRENE_EXIT_OWN
C
C
      SAVE
      IF (IFIRST.EQ.0) THEN
        ISAVE=ISTRA
        IFIRST=1
      ENDIF
C
      IF (TRCPLT) THEN
        WRITE (iunout,*) 'PLTEIR CALLED, ISTRA, XMCP: ',
     .                                   ISTRA,XMCP(ISTRA)
        IF (XMCP(ISTRA).EQ.0.0) THEN
          WRITE (iunout,*) 'PLOTTING ABANDONED FOR ALL OUTPUT TALLIES'
        ENDIF
      ENDIF
C
C  prepare plot frame.

C  NULLPUNKT AUF DEM PAPIER

      X0PL=10.
      Y0PL=3.
C  ACHSENLAENGEN
      LENX=25.
      LENY=20.
C  ACHSENUNTERTEILUNG VORGEGEBEN?
C  NEIN!
      STPSZX=0.
      STPSZY=0.
      INTNRX=0
      INTNRY=0
C  ACHSE LOGARITHMISCH?
      LOGX=.FALSE.
C     LOGY VIA INPUT
C  LOG. ACHSE MIN
      MINLY=0
C  LOG. ACHSE MAX
C     MAXLY WERDEN BERECHNET IN ANPSGL
C  ZEICHNE NETZLINIEN EIN
      GRIDX=.TRUE.
      GRIDY=.TRUE.
C  MACHE GRADE GRENZEN, X-ACHSE (Y ACHSE, NUR WENN TALZMI=TALZMA=666.)
      FITX=.TRUE.
C  NEW FRAME FOR EACH PICTURE IN PLTTLY
      L_SAME=.FALSE.
C
C  prepare plotting output tallies for present stratum ISTRA

      IF (XMCP(ISTRA).EQ.0.0) GOTO 10

C  PROVIDE EIRENE OUTPUT TALLIES FOR SELECTED STRATUM ISTRA
      IF (IESTR.EQ.ISTRA) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTRA
        IF (TRCFLE) WRITE (IUNOUT,*) 'FROM PLTEIR:'
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,
     .             NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .             TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRAD,NR1ST,NP2ND,NT3RD,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IESTR=ISTRA
        IF (TRCFLE) WRITE (IUNOUT,*) 'FROM PLTEIR:'
        CALL EIRENE_RSTRT(ISTRA,NSTRAI,
     .             NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .             TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRAD,NR1ST,NP2ND,NT3RD,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN PLTEIR: DATA FOR STRATUM ISTRA= ',
     .                    ISTRA
        WRITE (iunout,*) 'ARE NOT AVAILABLE. PLOTS ABANDONED'
        RETURN
      ENDIF

   10 CONTINUE
C
      IF (ISTRA.EQ.0) HEAD=
     . 'SUM OVER STRATA                                             '//
     . '          '
      IF (ISTRA.NE.0) THEN
        HEAD=
     .   'STRATUM NO.                                               '//
     .   '            '
        WRITE (HEAD(13:15),'(I3)') ISTRA
      ENDIF
C
      HEAD0=
     . 'VOLUME-AVERAGED BACKGROUND TALLY, INPUT                    '//
     . '           '
      HEAD1=
     . 'DEFAULT VOLUME-AVERAGED TALLY, TRACKLENGTH-ESTIMATED       '//
     . '           '
      HEAD2=
     . 'ADDITIONAL VOLUME-AVERAGED TALLY, TRACKLENGTH-ESTIMATED    '//
     . '           '
      HEAD3=
     . 'ADDITIONAL VOLUME-AVERAGED TALLY, COLLISION-ESTIMATED      '//
     . '           '
      HEAD4=
     . 'VOLUME-AVERAGED TALLY, SNAPSHOT-ESTIMATED                  '//
     . '           '
      HEAD5=
     . 'VOLUME-AVERAGED TALLY, FOR COUPLING TO PLASMA CODE         '//
     . '           '
      HEAD6=
     . 'BGK TALLY                                                  '//
     . '           '
      HEAD7=
     . 'ALGEBRAIC FUNCTION OF VOLUME-AVERAGED TALLIES              '//
     . '           '
      HEAD8=
     . 'RELATIVE STANDARD DEVIATION                                '//
     . '           '
      HEAD9=
     . 'SPECTRUM (VS. ENERGY, EV)                                  '//
     . '           '
      HEAD10=
     . 'SPECTRUM (VS. WAVELENGTH, NM)                              '//
     . '           '
C
      IALG=0
C
C  .......................................
C
C   LOOP OVER NVOLPL
C  .......................................
C
      IF (NVOLPL > 0) THEN
        ALLOCATE (VECTOR(NRAD,NPLT))
        IF (ANY(PLTL2D(1:NVOLPL).AND.PLTL3D(1:NVOLPL)))
     .     ALLOCATE (VECSAV(NRAD,NPLT))
        IF (ANY(PLTLER(1:NVOLPL))) THEN
          ALLOCATE (VSDVI(NRAD,NPLT))
          N1SDVI = NRAD
        ELSE
          ALLOCATE (VSDVI(1,1))
          N1SDVI = 1
        END IF
      END IF

      ALLOCATE(DUMMY(NRTAL))
      ALLOCATE(YMN2(NPLT), YMX2(NPLT), YMNLG2(NPLT), YMXLG2(NPLT))
      DO 10000 IBLD=1,NVOLPL
C
        IF (PLTL2D(IBLD).OR.PLTL3D(IBLD)) THEN
C
          DO 110 ICURV=1,NSPTAL(IBLD)
            JTAL=NPTALI(IBLD,ICURV)
C  REDO ALGEBRAIC TALLY IN CASE NFILEN=2 OR NFILEN=7
            IF (JTAL.EQ.NTALR.AND.IALG.EQ.0.AND.
     .          (NFILEN.EQ.2.OR.NFILEN.EQ.7)) THEN
              CALL EIRENE_ALGTAL
              IALG=1
              DO 105 IALV=1,NALVI
                DUMMY(1:NSBOX_TAL) = ALGV(IALV,1:NSBOX_TAL)
                CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,NSBOX_TAL,
     .                       ALGVI(IALV,ISTRA),
     .                       NR1TAL,NP2TAL,NT3TAL,NBMLT)
                ALGV(IALV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
  105         CONTINUE
            ENDIF

            ITL=IABS(JTAL)
C  PLOT OUTPUT TALLIES ONLY FOR STRATA WITH TWO OR MORE HISTORIES
            IF (JTAL.GT.0.AND.XMCP(ISTRA).LE.1) GOTO 10000
C  PLOT INPUT TALLIES ONLY ONCE PER ITERATION
            IF (JTAL.LT.0.AND.ISAVE.NE.ISTRA) GOTO 10000
c
            TXHEAD=HEAD0
            IF (JTAL.GT.0)     TXHEAD=HEAD1
            IF (JTAL.EQ.NTALA) TXHEAD=HEAD2
            IF (JTAL.EQ.NTALC) TXHEAD=HEAD3
            IF (JTAL.EQ.NTALT) TXHEAD=HEAD4
            IF (JTAL.EQ.NTALM) TXHEAD=HEAD5
            IF (JTAL.EQ.NTALB) TXHEAD=HEAD6
            IF (JTAL.EQ.NTALR) TXHEAD=HEAD7
            IF (TRCPLT) THEN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'PLOT REQUESTED FOR TALLY NO. ',JTAL
            ENDIF
C
C
C .............................................
C
C  PUT TALLY ONTO ARRAY: VECTOR
C .............................................
C
C
            LSDVI(ICURV)=.FALSE.
            LPLOT2(ICURV)=.FALSE.
            ISPZ=ISPTAL(IBLD,ICURV)

            IF (JTAL.LT.0.) THEN
cdr  here we deal with input tallies (and gradients thereof)
cdr  ITL = IABS(JTAL)
cdr  physical species index range: 1:NF, independent of possible indirect addressing
              NF=NFSTPI(ITL)
              VECTOR(:,ICURV)=0.

!  INPUT TALLY SWITCHED OFF ?
              IF (.NOT.LIVTALI(ITL)) THEN
                WRITE (iunout,*) TXTPLS(1,ITL)
                WRITE (iunout,*) ' TALLY NOT AVAILABLE (PLTEIR)',
     .                           ' JTAL = ', JTAL
                CALL EIRENE_LEER(1)
                CYCLE
              END IF

              IF (ISPZ.EQ.0) THEN
cdr  sum over species: this is nonsense in case of intensive quantities,
cdr                    such as Ti,V_in,
cdr                    and also in case of derivatives.
cdr  tbd: summing with proper weighting, as in outtal.f
                if (NF .gt. 1) then
                  write (iunout,*) 'DR: wrong code in plteir '
                  write (iunout,*) 'ITAL, ISPZ ',itl,ispz
                  write (iunout,*) 'plot abandonned for safety '
                  goto 110
                endif
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
cdr  this makes no sense. Ti cannot be summed.
                  VECTOR(1:NSBOX,ICURV) = SUM(TIIN(1:NF,1:NSBOX),1)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = SUM(DIIN(1:NF,1:NSBOX),1)
                CASE (5)
cdr  summing should use density weighting.
cdr  clearly wrong.
                  VECTOR(1:NSBOX,ICURV) = SUM(VXIN(1:NF,1:NSBOX),1)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = SUM(VYIN(1:NF,1:NSBOX),1)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = SUM(VZIN(1:NF,1:NSBOX),1)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = SUM(ADIN(1:NF,1:NSBOX),1)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = SUM(EDRIFT(1:NF,1:NSBOX),1)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = SUM(WGHT(1:NF,1:NSBOX),1)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE (18)
                  VECTOR(1:NSBOX,ICURV) = EXIN(1:NSBOX)
                CASE (19)
                  VECTOR(1:NSBOX,ICURV) = EYIN(1:NSBOX)
                CASE (20)
                  VECTOR(1:NSBOX,ICURV) = EZIN(1:NSBOX)
                CASE (21)
                  VECTOR(1:NSBOX,ICURV) = EFIN(1:NSBOX)
                CASE (22)
                  VECTOR(1:NSBOX,ICURV) = POT(1:NSBOX)
                CASE (23)
cdr vectorial tally: NF=3*NPLSI

cdr  so far: only parallel component can be plotted here
cdr  weighting ? unfinished...
                  VECTOR(1:NSBOX,ICURV) =
     .              SUM(BV_VEC(1:NPLSV,1:NSBOX),1)
                CASE (24)
cdr vectorial tally: NF=3*NPLSI

cdr  so far: only parallel component can be plotted here
cdr  weighting ? unfinished
                  VECTOR(1:NSBOX,ICURV) =
     .              SUM(PMOM_VEC(1:NPLSI,1:NSBOX),1)
                CASE (25)
                  VECTOR(1:NSBOX,ICURV) = PSI(1:NSBOX)
                CASE (26)
                  VECTOR(1:NSBOX,ICURV) = SUM(ZIIN(1:NF,1:NSBOX),1)
! INPUT TALLIES 27 -- 30:  CURRENTLY UNUSED (FREE)
! GRADIENTS OF INPUT TALLIES
                CASE (31:120)     ! ntali=120, constant required here
                  KK = NADDP(ITL)
                  VECTOR(1:NSBOX,ICURV) =
     .                   SUM(PLSTLS(KK+1:KK+NF,1:NSBOX),1)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED'
                  CALL EIRENE_LEER(1)
                  GOTO 10000
                END SELECT
cdr  done with sum over species index

              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NF) THEN
cdr  individual species indices
                SELECT CASE (ITL)
                CASE (1)
                  VECTOR(1:NSBOX,ICURV) = TEIN(1:NSBOX)
                CASE (2)
                  VECTOR(1:NSBOX,ICURV) = TIIN(MPLSTI(ISPZ),1:NSBOX)
                CASE (3)
                  VECTOR(1:NSBOX,ICURV) = DEIN(1:NSBOX)
                CASE (4)
                  VECTOR(1:NSBOX,ICURV) = DIIN(ISPZ,1:NSBOX)
                CASE (5)
                  VECTOR(1:NSBOX,ICURV) = VXIN(MPLSV(ISPZ),1:NSBOX)
                CASE (6)
                  VECTOR(1:NSBOX,ICURV) = VYIN(MPLSV(ISPZ),1:NSBOX)
                CASE (7)
                  VECTOR(1:NSBOX,ICURV) = VZIN(MPLSV(ISPZ),1:NSBOX)
                CASE (8)
                  VECTOR(1:NSBOX,ICURV) = BXIN(1:NSBOX)
                CASE (9)
                  VECTOR(1:NSBOX,ICURV) = BYIN(1:NSBOX)
                CASE (10)
                  VECTOR(1:NSBOX,ICURV) = BZIN(1:NSBOX)
                CASE (11)
                  VECTOR(1:NSBOX,ICURV) = BFIN(1:NSBOX)
                CASE (12)
                  VECTOR(1:NSBOX,ICURV) = ADIN(ISPZ,1:NSBOX)
                CASE (13)
                  VECTOR(1:NSBOX,ICURV) = EDRIFT(ISPZ,1:NSBOX)
                CASE (14)
                  VECTOR(1:NSBOX,ICURV) = VOL(1:NSBOX)
                CASE (15)
                  VECTOR(1:NSBOX,ICURV) = WGHT(ISPZ,1:NSBOX)
                CASE (16)
                  VECTOR(1:NSBOX,ICURV) = BXPERP(1:NSBOX)
                CASE (17)
                  VECTOR(1:NSBOX,ICURV) = BYPERP(1:NSBOX)
                CASE (18)
                  VECTOR(1:NSBOX,ICURV) = EXIN(1:NSBOX)
                CASE (19)
                  VECTOR(1:NSBOX,ICURV) = EYIN(1:NSBOX)
                CASE (20)
                  VECTOR(1:NSBOX,ICURV) = EZIN(1:NSBOX)
                CASE (21)
                  VECTOR(1:NSBOX,ICURV) = EFIN(1:NSBOX)
                CASE (22)
                  VECTOR(1:NSBOX,ICURV) = POT(1:NSBOX)
                CASE (23)
cdr physical species index range: 3*nplsi, even if indirect addressing
                  KS = MOD(ISPZ,NPLSI)
                  IF (KS == 0) KS=NPLSI   ! =1,,..,NPLSI
                  ICO=(ISPZ-1)/NPLSI      ! =0,1,2
                  KK=ICO*NPLSV
                  VECTOR(1:NSBOX,ICURV) = BV_VEC(KK+MPLSV(KS),1:NSBOX)
                CASE (24)
cdr  ispz= 1...3*nplsi
                  VECTOR(1:NSBOX,ICURV) = PMOM_VEC(ISPZ,1:NSBOX)
                CASE (25)
                  VECTOR(1:NSBOX,ICURV) = PSI(1:NSBOX)
                CASE (26)
                  VECTOR(1:NSBOX,ICURV) = ZIIN(ISPZ,1:NSBOX)
! INPUT TALLIES 27 -- 30:  CURRENTLY UNUSED (FREE)
! GRADIENTS
                CASE (31:120)     ! ntali=120, constant required here
                  KK = NADDP(ITL)+ISPZ
                  VECTOR(1:NSBOX,ICURV) = PLSTLS(KK,1:NSBOX)
                CASE DEFAULT
                  WRITE (iunout,*) ' WRONG TALLY NUMBER IN PLTEIR',
     .                        ' JTAL = ',JTAL
                  WRITE (iunout,*) ' NO PLOT PERFORMED'
                  CALL EIRENE_LEER(1)
                  CYCLE
                END SELECT

              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE'
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*)
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF'
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF

            ELSEIF (JTAL.GE.0) THEN
cdr  plot output tallies
              NFT=NFSTVI(ITL)
              NF=NFIRST(ITL)
              VECTOR(:,ICURV)=0.

              IF (.NOT.LIVTALV(JTAL)) THEN
                WRITE (iunout,*) TXTTAL(1,JTAL)
                WRITE (iunout,*) 'TALLY SWITCHED OFF'
                WRITE (iunout,*) 'ALL PLOTS FOR THIS TALLY TURNED OFF'
                CYCLE
              END IF

              IF (ISPZ.EQ.0) THEN
c  sum over species
c  output tallies are also intensive quantities, must be volume weighted.
c  But this cancels here. No density weighting as e.g. for Ti Vi input tallies.
                DO 122 K=1,NFT
                  DO I=1,NRAD
                    VECTOR(I,ICURV)=VECTOR(I,ICURV)+
     .                              ESTIMV(NADDV(ITL)+K,NCLTAL(I))
                  END DO
  122           CONTINUE
              ELSEIF (ISPZ.GT.0.AND.ISPZ.LE.NFT) THEN
                DO 125 I=1,NRAD
                  VECTOR(I,ICURV)=ESTIMV(NADDV(ITL)+ISPZ,NCLTAL(I))
  125           CONTINUE
              ELSE
                IF (TRCPLT) THEN
                  WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE'
                  WRITE (iunout,*) 'ICURV,ISPTAL(IBLD,ICURV) ',
     .                              ICURV,ISPZ
                  WRITE (iunout,*)
     .              'ALL PLOTS FOR THIS TALLY TURNED OFF'
                ENDIF
                PLTL2D(IBLD)=.FALSE.
                PLTL3D(IBLD)=.FALSE.
                GOTO 110
              ENDIF
C
              IF (PLTLER(IBLD)) THEN
C  CHECK IF STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY
                DO 126 N=1,NSIGVI
                  IF (IIH(N).NE.JTAL) GOTO 126
                  IF (IGH(N).NE.ISPZ.AND.IGH(N).NE.0) GOTO 126
                  LSDVI(ICURV)=.TRUE.
                  DO 127 I=1,NRAD
                    VSDVI(I,ICURV)=SIGMA(N,NCLTAL(I))
  127             CONTINUE
  126           CONTINUE
              ENDIF
C
            ENDIF
C
            IF (PLTL2D(IBLD) .AND. PLTL3D(IBLD)) THEN
              DO 129 I=1,NRAD
                VECSAV(I,ICURV)=VECTOR(I,ICURV)
  129         CONTINUE
            END IF
C
  110     CONTINUE
C
C ...................................
C                                   .
C    VECTOR(IC,ICURV) IS SET NOW    .
C ...................................
C
C
          LOGY=PLTLLG(IBLD)
C
          IF (PLTL2D(IBLD)) THEN
C
C  ............................
C
C   SET ABSCISSA FOR 2D PLOT
C  ............................
C
C  SET ABSCISSA FROM GRID DATA, OR USE ONE OF THE INPUT OPTIONS:
            IXSET2=0
            IF (TALXMI(IBLD).NE.0..OR.TALXMA(IBLD).NE.0.) THEN
              XMI=TALXMI(IBLD)
              XMA=TALXMA(IBLD)
              IA= 100000000
              IE=-100000000
              DO ICURV=1,NSPTAL(IBLD)
                IA=MIN(IA,NPLIN2(IBLD,ICURV))
                IE=MAX(IE,NPLOT2(IBLD,ICURV))
              ENDDO
              DEL=IE-IA
              IF (XMI.LT.XMA) THEN
C  EQUIDISTANT IN LIN SCALE
                DO I=IA,IE
                  XXP2D(I)=XMI+(I-IA)/DEL*(XMA-XMI)
                ENDDO
              ELSEIF (XMI.GT.XMA.AND.XMI.GT.0.AND.XMA.GT.0) THEN
C  EQUIDISTANT IN LOG SCALE
                XI=LOG(XMA)
                XE=LOG(XMI)
                DO I=IA,IE
                  XXP2D(I)=EXP(XI+(I-IA)/DEL*(XE-XI))
                ENDDO
C  USER-DEFINED ABSCISSA, XXP2D_USR
              ELSEIF (XMI.GT.XMA.AND.(XMI.LE.0.OR.XMA.LE.0)) THEN
                DO I=IA,IE
                  XXP2D(I)=XXP2D_USR(I,IBLD)
                ENDDO
              ENDIF
              IXSET2=1
              XMI=XXP2D(IA)*(1.+1.E-6)
              XMA=XXP2D(IE)/(1.+1.E-6)
              GOTO 139
            ENDIF
C
C  TRY DEFAULT OPTION TO SET PLOT GRID FROM 1ST (RADIAL) GRID
C
            IXSET2=0
            IF (LEVGEO.EQ.1.OR.LEVGEO.EQ.2) THEN
C   USE RADIAL SURFACE-CENTERED GRID "RHOSRF"
C   ...SAME FOR EACH Y- OR POLOIDAL, IF APPLICABLE
              DO 130 I=1,NR1ST
                XXP2D(I)=RHOSRF(I)
  130         CONTINUE
              DO 131 J=2,NP2ND*NT3RD*NBMLT
                DO I=1,NR1ST
                  XXP2D(I+(J-1)*NR1ST)=XXP2D(I)
                END DO
  131         CONTINUE
              XXP2D(NSURF+1:NRAD)=0.
              IXSET2=1
            ELSEIF (LEVGEO.EQ.3) THEN
C   USE PERPEND. ARCLENGTH "BGLP" IN CASE OF POLYGON GRID,
C   ...FOR EACH POLOIDAL AND TOROIDAL POSITION, IF APPLICABLE
              DO 133 I=1,NR1ST
                DO J=1,NP2ND
                  DO K=1,NT3RD
                    IRAD=I+((J-1)+(K-1)*NP2T3)*NR1P2
                    XXP2D(IRAD)=BGLP(I,J)
                  END DO
                END DO
  133         CONTINUE
              XXP2D(NSURF+1:NRAD)=0.
              IXSET2=1
            ELSE
C   NO 2D PLOT OPTIONS AVAILABLE
            ENDIF
            XMI=XXP2D(NPLIN2(IBLD,1))*(1.+1.E-6)
            XMA=XXP2D(NPLOT2(IBLD,1))/(1.+1.E-6)
C
  139       CONTINUE
C
            IF (IXSET2.NE.1) THEN
              WRITE (iunout,*) ' NO GRID SET FOR 2D PLOTTING '
              WRITE (iunout,*) ' IXSET2 = ',IXSET2
              WRITE (iunout,*) ' NO 2D PLOTTING IS DONE '
              GOTO 1000
            ENDIF
C
C  IN CASE OF LSMOT2, SET ZONE-CENTERED ABSCISSA
C  GRID FROM SURFACE-CENTERED GRID "X"
C
            IF (LSMOT2(IBLD)) THEN
              DO 137 J=1,NRAD-1
                XXP2D(J)=(XXP2D(J)+XXP2D(J+1))*0.5
  137         CONTINUE
            ENDIF
C
            DO 140 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IR1(ICURV)=NPLIN2(IBLD,ICURV)
              IR2(ICURV)=NPLOT2(IBLD,ICURV)
              IRS(ICURV)=NPLDL2(IBLD,ICURV)
C
              YMNLG2(ICURV)=1.D60
              YMXLG2(ICURV)=-1.D60
              IF (JTAL.GT.0.) THEN
                JJTAL = JTAL
                IF (NEXTVI(JTAL) > 0) THEN
                  KK = 0
                  DO K=1,JTAL
                    KK = KK + 1
                    IF ((K > NEXTVI(JTAL) .AND.
     .                  (MOD(K,NEXTVI(JTAL)) == 1))) KK = KK + 1
                  END DO
                  JJTAL = KK
                END IF
                CALL EIRENE_FETCH_OUTAUI (OUTAUI,JJTAL,ISPZ,ISTRA,
     .           IUNOUT)
                IF (OUTAUI.EQ.0.) THEN
                  IF (TRCPLT) THEN
                    WRITE (iunout,*) 'TALLY NO. ',JTAL,
     .                               ' CURVE NO. ',ICURV
                    WRITE (iunout,*) 'NOT PLOTTED BECAUSE'
                    WRITE (iunout,*) 'ZERO INTEGRAL (OUTAUI=0.)'
                  ENDIF
                  YMN2(ICURV)=0.
                  YMX2(ICURV)=0.
                  YMNLG2(ICURV)=0.
                  YMXLG2(ICURV)=0.
                  GOTO 140
                ENDIF
              ENDIF
              LPLOT2(ICURV)=.TRUE.
              I1=IR1(ICURV)
              I2=IR2(ICURV)
              I2M=I2-1
              IS=IRS(ICURV)
C
C YMNLG2, YMXLG2: REAL MAX/MIN, FOR LEGENDE ON 2D PLOT ONLY
              DO 141 I=I1,I2M,IS
                YMNLG2(ICURV)=MIN(YMNLG2(ICURV),VECTOR(I,ICURV))
  141         CONTINUE
              DO 142 I=I1,I2M,IS
                YMXLG2(ICURV)=MAX(YMXLG2(ICURV),VECTOR(I,ICURV))
  142         CONTINUE
C
C YMN2, YMX2: FOR AXIS
              FITY=.TRUE.
              IF (TALZMI(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMN2(ICURV)=TALZMI(IBLD)
                DO 143 I=1,NRAD
                  VECTOR(I,ICURV)=MAX(YMN2(ICURV),VECTOR(I,ICURV))
  143           CONTINUE
                IF (LOGY) YMN2(ICURV)=YMN2(ICURV)*(1.+1.E-6)
              ELSE
                YMN2(ICURV)=YMNLG2(ICURV)
              ENDIF
C
              IF (TALZMA(IBLD).NE.666.) THEN
                IF (.NOT.LOGY) FITY=.FALSE.
                YMX2(ICURV)=TALZMA(IBLD)
                DO 144 I=1,NRAD
                  VECTOR(I,ICURV)=MIN(YMX2(ICURV),VECTOR(I,ICURV))
  144           CONTINUE
                IF (LOGY) YMX2(ICURV)=YMX2(ICURV)/(1.+1.E-6)
              ELSE
                YMX2(ICURV)=YMXLG2(ICURV)
              ENDIF
  140       CONTINUE
C
C  PLOT ALL CURVES REQUESTED FROM THIS TALLY INTO ONE PICTURE
            DO 150 ICURV=1,NSPTAL(IBLD)
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPEC(ICURV)='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNIT(ICURV)=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTALL(ICURV)=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNIT(ICURV)=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTALL(ICURV)=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTALL(ICURV)=TXTPLS(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTPSP(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTALL(ICURV)=TXTTAL(ISPZ,ITL)
                  TXSPEC(ICURV)=TXTSPC(ISPZ,ITL)
                  TXUNIT(ICURV)=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
  150       CONTINUE
            IERR=0
            L_SAME=.FALSE.
            CALL EIRENE_PLTTLY (XXP2D,VECTOR,VSDVI,YMN2,YMX2,
     .             IR1,IR2,IRS,
     .             NSPTAL(IBLD),TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .             LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,LHIST2(IBLD),IERR,
     .             N1SDVI,NRAD,L_SAME)
            IF (TRCPLT) THEN
              IF (IERR.GT.0) THEN
                WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,
     .                           ' ABANDONED'
                WRITE (iunout,*) 'ERROR CODE FROM SUBR. PLTTLY: ',IERR
                WRITE (iunout,*) 'XMI,XMA ',XMI,XMA
                GOTO 1000
              ENDIF
              WRITE (iunout,*) '2D PLOT FOR TALLY NO. ',JTAL,' DONE'
              WRITE (iunout,*) 'XMIN= ',XMI,' XMAX= ',XMA
              DO 160 ICURV=1,NSPTAL(IBLD)
                IF (LPLOT2(ICURV))
     .            WRITE (iunout,*) 'ICURV= ',ICURV,
     .                        ' YMIN= ',YMNLG2(ICURV),
     .                        ' YMAX= ',YMXLG2(ICURV),
     .                        ' LSDVI= ',LSDVI(ICURV)
  160         CONTINUE
            ENDIF
C
          ENDIF
C
 1000     CONTINUE
C
C   3D PLOT GRID
C
          IF (PLTL3D(IBLD)) THEN
C
            DO 1040 ICURV=1,NSPTAL(IBLD)
              IF (PLTL2D(IBLD)) THEN
                DO 1035 I=1,NRAD
                   VECTOR(I,ICURV)=VECSAV(I,ICURV)
 1035           CONTINUE
              END IF
C  SYMMETRY CONDITION AT POLAR ANGLE THETA=YIA AND THETA=2*PI+YIA
C  NOT READY: IXTL3 NOT DEFINED HERE. ENFORCE SYMMETRY AUTOMATICALLY EARLIER
C             IF (LEVGEO.EQ.2.AND.IYTL3.EQ.NP2ND) THEN
C               DO 1036 I=1,IXTL3
C                 VECTOR(I+NP2NDM*NR1ST,ICURV)=VECTOR(I,ICURV)
C 1036          CONTINUE
C             ENDIF
 1040       CONTINUE
C
C SET QUASIRECTANGULAR PLOT GRIDS XXP3D (IX), IX=1,IXTL3
C                             AND YYP3D (IY), IY=1,IYTL3
C FOR EACH 3D PICTURE
C
            IXSET3=0
            IYSET3=0
            IF (LEVGEO.EQ.1) THEN

              IF (NLPOL.AND.NLTOR.AND.NLTRZ.AND.LPRAD3(IBLD)) THEN
c  set a y-z grid, by abuse of notation on xxp3d,yyp3d
                IXTL3=NP2ND
                DO I=1,IXTL3
                  XXP3D(I)=PSURF(I)
                END DO
                IXSET3=1
                IYTL3=NT3RD
                DO I=1,IYTL3
                  YYP3D(I)=ZSURF(I)
                END DO
                DO I=1,NP2ND
                  DO J=1,NT3RD
                    XPOL(I,J) = PSURF(I)
                    YPOL(I,J) = ZSURF(J)
                  END DO
                END DO
                IYSET3=1
              END IF

              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
c  at this point: either lppol3 or lptor3 must be true
c  set a x-y or a x-z grid, by abuse of notation on xxp3d,yyp3d
                IXTL3=NR1ST
                XXP3D(1:IXTL3)=RHOSRF(1:IXTL3)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
c  at this point: lppol3 must be true, i.e. we need x-z grid
                IYTL3=NT3RD
                YYP3D(1:IYTL3)=ZSURF(1:IYTL3)
                DO I=1,NR1ST
                  DO J=1,NT3RD
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=ZSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
c  at this point: lptor3 must be true, i.e. we need x-y grid
                IYTL3=NP2ND
                YYP3D(1:IYTL3)=PSURF(1:IYTL3)
                DO I=1,NR1ST
                  DO J=1,NP2ND
                    XPOL(I,J)=RHOSRF(I)
                    YPOL(I,J)=PSURF(J)
                  ENDDO
                ENDDO
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.2) THEN
C
              IF (NLRAD.AND..NOT.LPRAD3(IBLD)) THEN
                IXTL3=NR1ST
                XXP3D(1:IXTL3)=RHOSRF(1:IXTL3)
                IXSET3=1
              ENDIF
              IF (NLTOR.AND.NLTRZ.AND..NOT.LPTOR3(IBLD)) THEN
                IYTL3=NT3RD
                YYP3D(1:IYTL3)=ZSURF(1:IYTL3)
                IYSET3=1
              ENDIF
              IF (NLPOL.AND..NOT.LPPOL3(IBLD)) THEN
                IYTL3=NP2ND
                YYP3D(1:IYTL3-1)=0.5*(PSURF(2:IYTL3)+PSURF(1:IYTL3-1))
                YYP3D(NP2ND)=PSURF(1)+PI2A
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.3) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXTL3=NR1ST
                DO 228 IX=1,IXTL3
                  XXP3D(IX)=IX
  228           CONTINUE
                IXSET3=1
                IYTL3=NP2ND
                DO 230 IX=1,IYTL3
                  YYP3D(IX)=IX
  230           CONTINUE
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.4) THEN
C
              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
              ENDIF
C
            ELSEIF (LEVGEO.EQ.5) THEN
C
!              IF (LPTOR3(IBLD)) THEN
                IXSET3=1
                IYSET3=1
!              ENDIF
            ELSE
              WRITE (iunout,*) '3D PLOT OPTION TO BE WRITTEN, LEVGEO '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C
C  LOOP ICURV=1,....
C
            ICINC=1
            IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) ICINC=2
            DO 1160 ICURV=1,NSPTAL(IBLD),ICINC
              JTAL=NPTALI(IBLD,ICURV)
              ITL=IABS(JTAL)
              ISPZ=ISPTAL(IBLD,ICURV)
              IF (ISPZ.EQ.0) THEN
                TXSPC1='SUM OVER SPECIES        '
                IF (JTAL.LT.0) TXUNT1=TXTPUN(1,ITL)
                IF (JTAL.LT.0) TXTLL1=TXTPLS(1,ITL)
                IF (JTAL.GE.0) TXUNT1=TXTUNT(1,ITL)
                IF (JTAL.GE.0) TXTLL1=TXTTAL(1,ITL)
              ELSE
                IF (JTAL.LT.0) THEN
                  TXTLL1=TXTPLS(ISPZ,ITL)
                  TXSPC1=TXTPSP(ISPZ,ITL)
                  TXUNT1=TXTPUN(ISPZ,ITL)
                ELSE
                  TXTLL1=TXTTAL(ISPZ,ITL)
                  TXSPC1=TXTSPC(ISPZ,ITL)
                  TXUNT1=TXTUNT(ISPZ,ITL)
                ENDIF
              ENDIF
C
              IF (IXSET3+IYSET3.LT.2) THEN
                WRITE (iunout,*) ' NO GRIDS SET FOR 3D PLOTTING '
                WRITE (iunout,*) ' IXSET3,IYSET3 = ',IXSET3,IYSET3
                WRITE (iunout,*) ' NO 3D PLOTTING IS DONE '
                GOTO 10000
              ENDIF
C
              LINLOG=PLTLLG(IBLD)
              TMIN=TALZMI(IBLD)
              TMAX=TALZMA(IBLD)
C
 1200         CONTINUE
C
C
C  CONTOUR PLOTS
              IF (LCNTR3(IBLD)) THEN
                CALL EIRENE_ISOLNE (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS PLOTS
              ELSEIF (LRAPS3(IBLD)) THEN
                CALL EIRENE_RPSCOL (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D HISTOGRAM
              ELSEIF (LHIST3(IBLD)) THEN
                CALL EIRENE_PL3DPG (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  3D SURFACE PLOTS, IN CUBE
              ELSEIF (LSMOT3(IBLD)) THEN
                CALL EIRENE_PLOT3D (VECTOR(1,ICURV),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D,YYP3D,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,TALW1(IBLD),TALW2(IBLD),
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  VECTOR FIELD PLOT
              ELSEIF (LVECT3(IBLD)) THEN
                IXXI=NPLI13(IBLD,ICURV)
                IXXE=NPLO13(IBLD,ICURV)
                IYYI=NPLI23(IBLD,ICURV)
                IYYE=NPLO23(IBLD,ICURV)
                CALL EIRENE_VECLNE (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXXI,IXXE,IYYI,IYYE,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
C  WRITE FILES FOR RAPS VECTOR PLOTS
              ELSEIF (LRPVC3(IBLD)) THEN
                CALL EIRENE_RPSVEC (VECTOR(1,ICURV),
     .                       VECTOR(1,ICURV+1),IBLD,ICURV,
     .                       IXTL3,IYTL3,XXP3D_DUM,YYP3D_DUM,
     .                       TXTLL1,TXSPC1,TXUNT1,
     .                       LINLOG,TMAX,TMIN,
     .                       HEAD,TXTRUN,TXHEAD,TRCPLT)
              ELSE
                WRITE (iunout,*) 'NO 3D PLOT OPTION FOR IBLD= ',IBLD
                GOTO 10000
              ENDIF
C
C   PLOT STANDARD DEVIATION PROFILE FOR THIS TALLY, IF REQUESTED
C
              IF (LVECT3(IBLD).OR.LRPVC3(IBLD)) GOTO 1160
              IF (PLTLER(IBLD).AND.LSDVI(ICURV)) THEN
                TXUNT1='%                       '
                TXHEAD=HEAD8
                INULL=0
                DO 1222 I=1,NRAD
                  VECTOR(I,ICURV)=VSDVI(I,ICURV)
                  IF (LRAPS3(IBLD).AND.
     .               (ABS(VECTOR(I,ICURV)) < EPS30)) THEN
                    VECTOR(I,ICURV)=101._DP
                    INULL = INULL + 1
                  END IF
 1222           CONTINUE
                LINLOG=.FALSE.
                TMIN=0.
                TMAX=100.
                IF (LRAPS3(IBLD).AND.INULL.GT.0) THEN
                  TMAX=101.
                  WRITE (IUNOUT,*) 'RAPS GRAPHICS FOR STD. DEVIATION:'
                  WRITE (IUNOUT,*) 'IBLD, ICURV ',IBLD, ICURV
                  WRITE (IUNOUT,*)  INULL, ' CELLS WITH 0 HISTORIES'
                  WRITE (IUNOUT,*) 'STD. DEV. SET = 101% IN THESE CELLS'
                  WRITE (IUNOUT,*) 'TO PERMIT SPECIAL CHOICE OF COLOUR'
                END IF
                LSDVI(ICURV)=.FALSE.
                GOTO 1200
              ENDIF
C
 1160       CONTINUE
C  LOOP ICURV FINISHED
          ENDIF
C
        ELSEIF (NSPTAL(IBLD).GT.0) THEN
          WRITE (iunout,*) 'PLOT REQUEST FOR TALLY NO. ',JTAL,
     .                     ' BUT NEITHER'
          WRITE (iunout,*)
     .       'PLTL2D NOR PLTL3D TRUE. NO PLOT FOR THIS TALLY'
        ENDIF
C
10000 CONTINUE
      DEALLOCATE(DUMMY)

C  LOOP IBLD FINISHED, NO PICTURE PRODUCED IN CASE XMCP=0 AND OUTPUT TALLY REQUESTED
C
C  NEXT: PLOT ENERGY (WAVELENGTH) SPECTRA, IF ANY HAVE BEEN SCORED
C        PLOTTING IS NOT YET CONDITIONED BY FLAGS
C        ALL PLOTS FOR ALL SPECTRA ARE ALWAYS DONE
C
      IF (XMCP(ISTRA).LE.1.0) GOTO 20000
C
      DO ISPC=1,NADSPC
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
C  THESE ARE EQUALLY SPACED LINEARLY OR LOGARITHMICALLY
C       LOGX=ESTIML(ISPC)%LOG
        NSPS=ESTIML(ISPC)%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        SPCAN=ESTIML(ISPC)%SPCMIN
        SPC00=ESTIML(ISPC)%ESP_00
        if (logx) then
c  to be written: log x scale
        endif
C  x axis: ENERGY-BIN FACES
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%SPCDEL-SPC00
        END DO
C  y axis: ENERGY BIN AVERAGES (approx: value at energy-bin centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%SGM(I)
        END DO

        YMN2(1)=MINVAL(YSPEC(1:NSPS,1))
        YMX2(1)=MAXVAL(YSPEC(1:NSPS,1))
        IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
        YMNLG2(1)=YMN2(1)
        YMXLG2(1)=YMX2(1)

CDR:  NOT READY: ABUSE SPCPLT FOR MIN MAX ON PLOT, ALWAYS: LIN-LOG SCALE
        IF (ESTIML(ISPC)%SPC_XPLT.NE.666.) THEN
          YMN2(1)=ESTIML(ISPC)%SPC_XPLT
        ENDIF
        IF (ESTIML(ISPC)%SPC_YPLT.NE.666.) THEN
          YMX2(1)=ESTIML(ISPC)%SPC_YPLT
        ENDIF

        LSDVI(1)=NSIGI_SPC > 0
        LPLOT2(1)=.TRUE.
        IR1(1)=1
        IR2(1)=NSPS+1
        IRS(1)=1
        XMI=XSPEC(1)
        XMA=XSPEC(NSPS+1)

C  LINEAR OR LOGARITHMIC Y SCALE ?
        LOGY=.FALSE.
        FITY=.FALSE.

CDR  PLOT SPECTRUM LOGARITHMICALLY ???
CDR     IF (ESTIML(ISPC)%SPC_YPLT.GT.0.0) THEN
        IF (.TRUE.) THEN
          LOGY=.TRUE.
          FITY=.TRUE.
        ENDIF
CDR

        IF (ESTIML(ISPC)%ISRFCLL == 0) THEN
         TXTALL(1)='SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .             'SPECIES                '
        ELSE
         TXTALL(1)='SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .             'SPECIES                '
        ENDIF
        WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%ISPCSRF
        WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%IPRTYP
        WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%IPRSP
        IT = ESTIML(ISPC)%ISPCTYP
        ITT= ESTIML(ISPC)%ISRFCLL
        TXSPEC=REPEAT(' ',24)
        TXUNIT=REPEAT(' ',24)
        IF (ITT.EQ.0.AND.IT == 1) TXUNIT='AMP/BIN(EV)             '
        IF (ITT.EQ.0.AND.IT == 2) TXUNIT='WATT/BIN(EV)            '
        IF (ITT.EQ.1.AND.IT == 1) TXUNIT='#/CM**3/BIN(EV)         '
        IF (ITT.EQ.1.AND.IT == 2) TXUNIT='EV/CM**3/BIN(EV)        '
cdr  itt=2 was still missing....  units probably: (TO BE CHECKED)
cdr  June20: probably IDIREC controls directional vs. nondirectional spectrum
        IF (ITT.EQ.2.AND.IT == 1) TXUNIT='#/CM**3/BIN(EV)/STERAD  '
        IF (ITT.EQ.2.AND.IT == 2) TXUNIT='EV/CM**3/BIN(EV)/STERAD '
        TXHEAD=REPEAT(' ',72)
        TXHEAD(1:30)=HEAD9(1:30)
        TXHEAD(32:42)='INTEGRAL: '
        WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%SPCS
        IERR=0
C  MANY SPECTRA INTO ONE PICTURE
        L_SAME=ESTIML(ISPC)%SPC_SAME .NE. 1.D0
C  ENFORCE NEW FRAME FOR 1ST SPECTRUM
        IF (ISPC.EQ.1) L_SAME=.FALSE.
        CALL EIRENE_PLTTLY (XSPEC,YSPEC,VSPEC,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

      END DO

C  NOW REPEAT SAME PLOTS, BUT VS. WAVELENGTH
      IF (NPHOTI > 0) THEN
CDR TO BE DONE: DISTINGUISH BETWEEN PHOTON AND PARTICLE SPECTRA

      DO ISPC=1,NADSPC
        ITP=ESTIML(ISPC)%IPRTYP
        IF (ITP.NE.0) CYCLE
C  THERE ARE NSPS BINS, AND NSPS+1 ENERGY BIN BOUNDARIES
        NSPS=ESTIML(ISPC)%NSPC
        ALLOCATE (XSPEC(NSPS+1))
        ALLOCATE (YSPEC(NSPS+1,1))
        ALLOCATE (VSPEC(NSPS+1,1))
        ALLOCATE (WLSPEC(NSPS+1))
        ALLOCATE (YSPECWL(NSPS+1,1))
        ALLOCATE (VSPECWL(NSPS+1,1))
        SPCAN=ESTIML(ISPC)%SPCMIN
        SPC00=ESTIML(ISPC)%ESP_00
C  x axis: ENERGY-BIN FACES
        DO I=1,NSPS+1
          XSPEC(I)=SPCAN+(I-1)*ESTIML(ISPC)%SPCDEL
        END DO
C  y axis: ENERGY BIN AVERAGES (approx: value at energy-bin centres)
        DO I=1,NSPS
          YSPEC(I,1)=ESTIML(ISPC)%SPC(I)
          IF (NSIGI_SPC > 0) VSPEC(I,1)=ESTIML(ISPC)%SGM(I)
        END DO

C  PLOT ALSO VS. WAVELENGTH (NM)
        WL00            =HPCL/MAX(1.E-6_DP,SPC00)*1.E7_DP
C  x axis: cell faces
        DO I=1,NSPS+1
          WLSPEC(NSPS-I+1+1)=HPCL/MAX(1.E-6_DP,XSPEC(I))*1.E7_DP
          WLSPEC(NSPS-I+1+1)=WLSPEC(NSPS-I+1+1)-WL00
        END DO
C  y axis: cell averages (approx: cell centres)
        DO I=1,NSPS
          YSPECWL(NSPS-I+1,1)=YSPEC(I,1)
          IF (NSIGI_SPC > 0) VSPECWL(NSPS-I+1,1)=VSPEC(I,1)
        END DO
C  rescaling: flux/ev to flux/nm
        DO I=1,NSPS
          DE=XSPEC(NSPS-I+1+1)-XSPEC(NSPS-I+1)
          DW=WLSPEC(I+1)-WLSPEC(I)
          YSPECWL(I,1) = YSPECWL(I,1)*DE/DW
        END DO

        DEALLOCATE (XSPEC)
        DEALLOCATE (YSPEC)
        DEALLOCATE (VSPEC)

        YMN2(1)=MINVAL(YSPECWL(1:NSPS,1))
        YMX2(1)=MAXVAL(YSPECWL(1:NSPS,1))
        IF (ABS(YMX2(1)-YMN2(1)) < EPS30) YMX2(1) = YMN2(1) + 1._dp
        YMNLG2(1)=YMN2(1)
        YMXLG2(1)=YMX2(1)
        LSDVI(1)=NSIGI_SPC > 0
        LPLOT2(1)=.TRUE.
        IR1(1)=1
        IR2(1)=NSPS+1
        IRS(1)=1
        XMI=WLSPEC(1)
        XMA=WLSPEC(NSPS+1)
        LOGY=.TRUE.
        FITY=.TRUE.
        IF (ESTIML(ISPC)%ISRFCLL == 0) THEN
          TXTALL(1)=
     .    'SPECTRUM FOR SURFACE        PARTICLE TYPE        '//
     .    'SPECIES                '
        ELSE
          TXTALL(1)=
     .    'SPECTRUM FOR CELL           PARTICLE TYPE        '//
     .    'SPECIES                '
        END IF
        WRITE (TXTALL(1)(22:27),'(I6)') ESTIML(ISPC)%ISPCSRF
        WRITE (TXTALL(1)(43:48),'(I6)') ESTIML(ISPC)%IPRTYP
        WRITE (TXTALL(1)(58:63),'(I6)') ESTIML(ISPC)%IPRSP
        IT = ESTIML(ISPC)%ISPCTYP
        ITT= ESTIML(ISPC)%ISRFCLL
        TXSPEC=REPEAT(' ',24)
        TXUNIT=REPEAT(' ',24)
        IF (ITT.EQ.0.AND.IT == 1) TXUNIT='AMP/BIN(NM)             '
        IF (ITT.EQ.0.AND.IT == 2) TXUNIT='WATT/BIN(NM)            '
        IF (ITT.EQ.1.AND.IT == 1) TXUNIT='#/CM**3/BIN(NM)         '
        IF (ITT.EQ.1.AND.IT == 2) TXUNIT='EV/CM**3/BIN(NM)        '
cdr  itt=2 was still missing....  units probably: (TO BE CHECKED)
        IF (ITT.EQ.2.AND.IT == 1) TXUNIT='#/CM**3/BIN(EV)/STERAD  '
        IF (ITT.EQ.2.AND.IT == 2) TXUNIT='EV/CM**3/BIN(EV)/STERAD '
        TXHEAD=REPEAT(' ',72)
        TXHEAD(1:30)=HEAD10(1:30)
        TXHEAD(32:42)='INTEGRAL: '
        WRITE (TXHEAD(43:55),'(ES12.4)') ESTIML(ISPC)%SPCS
        IERR=0
        L_SAME=.TRUE.
        IF (ISPC.EQ.1) L_SAME=.FALSE.
        L_SAME=ESTIML(ISPC)%SPC_SAME .NE. 1.D0
        CALL EIRENE_PLTTLY (WLSPEC,YSPECWL,VSPECWL,YMN2,YMX2,
     .       IR1,IR2,IRS,
     .       1,TXTALL,TXSPEC,TXUNIT,TXTRUN,TXHEAD,
     .       LSDVI,XMI,XMA,YMNLG2,YMXLG2,LPLOT2,.TRUE.,IERR,
     .       NSPS+1,NSPS+1,L_SAME)
        DEALLOCATE (WLSPEC)
        DEALLOCATE (YSPECWL)
        DEALLOCATE (VSPECWL)
      END DO  !  LOOP OVER PHOTON SPECTRA ENDS HERE

      END IF

20000 CONTINUE

      DEALLOCATE(YMN2,YMX2,YMNLG2,YMXLG2)
      IF (ALLOCATED(VECTOR)) DEALLOCATE(VECTOR)
      IF (ALLOCATED(VECSAV)) DEALLOCATE(VECSAV)
      IF (ALLOCATED(VSDVI))  DEALLOCATE(VSDVI)
      RETURN
      END SUBROUTINE EIRENE_PLTEIR

C     the following SUBROUTINE is for reinitialization of EIRENE (DMH)

      SUBROUTINE EIRENE_PLTEIR_REINIT
      IMPLICIT NONE
      IFIRST = 0
      return
      END SUBROUTINE EIRENE_PLTEIR_REINIT

      END MODULE EIRMOD_PLTEIR
