C  aug. 05: corrected electron energy loss rate for default rec. rate
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
! 25.03.07: check of mass conservation only for up to two secondaries
! 2013    : DSUB (RESCALING OF DENSITY IN H.4 FITS) REMOVED, NOW DONE IN RATE_COEFF.F
! 2013    : DENSITY LIMIT 1E8 SET FOR POLYNOMIAL FITS (ARRAY PLS).
cdr  oct.14:  pls made allocatable, plus minor synchronisation with other xsect... routines
cdr  Nov.14:  reaction scaling factor removed from Bremsstrahlung.
CDR           bremsstrahlung: new function eirene_brems, replaces Gaunt factor function
cdr  June 15:  added: default He+ --> He(1S) + rad  model. same analytic form of rate as for H+ default model.
cdr  April 16:  typo re TABRC1 for default He recombination corrected. Correction by SOLPS-ITER group
cdr             should not have had any effect, on any run, so far,
cdr             since this reaction did not exist in EIRENE at all until June 15
cdr  Jan 18  :  call energy_rate_coeff with lexp=true, because internal colrad (ifit=5)
cdr             option is now available.
cdr  May 18 :  still missing. low Te cut-off (should be done as in xstei, there:
cdr            0.1 eV, until asymptotics from database are fully implemented.
cdr            DEIMIN density cut-off now redundant, due to defaults read from AMJUEL

C
      SUBROUTINE EIRENE_XSECTP
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_PHOTON

      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: PLS(:), TEPLS(:)
      REAL(DP) :: DELE, FCTKKL, ZX, DEIMIN, RMASS2, FACTKK,
     .            RMASS2_2, CORSUM, COU, EIRENE_RATE_COEFF,
     .            EIRENE_ENERGY_RATE_COEFF,
     .            BREMS, Z, eirene_brems,
     .            DEIMAX, TEIMIN
      INTEGER :: IIRC, IION3, IPLS3, IATM3, IMOL3, KK, NRC, IATM,
     .           IRRC, J, IDSC, IPLS, NSERC5, KREAD, MODC, IATM1,
     .           ITYP, ISPZ, ITYP2, ISPZ2, IPHOT3
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LEXP, LADAS
      EXTERNAL :: EIRENE_ENERGY_RATE_COEFF, EIRENE_RATE_COEFF,
     .            EIRENE_BREMS, EIRENE_XSTRC,
     .            EIRENE_LEER, EIRENE_MASBOX, EIRENE_EXIT_OWN
      SAVE

      ALLOCATE (PLS(NSTORDR))
      ALLOCATE (TEPLS(NSTORDR))

cdr: set hard-wired lower density for H.4, H.10 type fits from AMJUEL: 1e8 cm**-3
cdr: at this lower limit density the fits are produced such
cdr: that they collapse to the corona limit values.
      DEIMIN=LOG(1.D8)
      DEIMAX=LOG(1.D16) !VK
      TEIMIN=log(0.1d0) !csw
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))
          PLS(J)=MIN(DEIMAX,PLS(J))        !VK
          TEPLS(J)=MAX(TEIMIN,TEINL(J))
   10   CONTINUE
      END IF

C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
C
        IDSC=0
        LGPRC(IPLS,0)=0
C
        DO NRC=1,NRCP(IPLS)
          KK=IREACP(IPLS,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.7) GOTO 994
        ENDDO
C
        IF (NRCP(IPLS).EQ.0) THEN   ! default model for bulk ion IPLS
C
          IF (NCHARP(IPLS).EQ.1 .AND. NCHRGP(IPLS).EQ.1) THEN
! this is now H+, or D+, or T+
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL, for capture on all levels of H
C  HYDR. RECOMBINATION RATE COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            DO 52 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM)) THEN
C
                IDSC=IDSC+1
                NRRCI=NRRCI+1
                IF (NRRCI.GT.NREC) GOTO 992
                IRRC=NRRCI
                LGPRC(IPLS,IDSC)=IRRC

                IF (NSTORDR >= NRAD) THEN
                  DO 51 J=1,NSBOX
                    ZX=EIONH/MAX(1.E-5_DP,TEIN(J))
C  rate = rate coeff: <sig v> times electr. density, 1/s per ion
                    TABRC1(IRRC,J)=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(J)
C  maxw. electron energy loss rate due to recombination
c                   corsum=0._dp  ! old default: 1.5*Te
C  correction due to energy dependence in rec. cross-section
C  corsum=d(ln<sig v>)/d(ln Te)
c  corsum approx -0.5 for Te --> 0
c  corsum approx  0.0 for Te approx 11.43
c  corsum approx +1.0 for Te --> infinity
                    corsum=(-0.5_dp*zx+0.59)/(zx+0.59)
                    EELRC1(IRRC,J)=-(1.5+CORSUM)*TEIN(J)*TABRC1(IRRC,J)
   51             CONTINUE
cdr  this setting kk=-1 is confusing. It may work, but
cdr  kk=-1 is already reserved for "minimal" default H+p charge exchange process
                  NREARC(IRRC) = -1
                  NELRRC(IRRC) = -1
                ELSE          ! storage saving mode:
                              ! tabrc1, eelrc1 to be found "on the fly"
                  NREARC(IRRC) = -1
                  NELRRC(IRRC) = -1
                END IF
                IATM1=IATM
                NATPRC(IRRC)=IATM1
                NIOPRC(IRRC)=0
                NPLPRC(IRRC)=0
                NMLPRC(IRRC)=0
C
                MODCOL(6,2,IRRC)=1
                MODCOL(6,4,IRRC)=1
              ENDIF
   52       CONTINUE
C
            NPRCI(IPLS)=IDSC

cdr Sept. 19 added: nprt=1, to avoid confusing D2+ and He+ virtual background here
          ELSEIF (NCHARP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1
     .      .AND. NPRT(ISPZ).EQ.1) THEN  ! this is now He+
C
C  DEFAULT HELIUM + RADIATIVE RECOMBINATION MODEL
C  HELIUM ION (HE+) RECOMBINATION RATE COEFFICIENT (1/S/CCM) E + HE+ --> HE(1S) + RAD.
C  JANEV ET. AL. FORMULA H.2. 2.3.13, BASED ON SOBELMAN 1979
C  (BORN-COULOMB APPROXIMATION), SIMILAR EXPRESSION AS FOR HYDROGEN DEFAULT RECOMBINATION MODEL
C
            DO 54 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM)) THEN
C
                IDSC=IDSC+1
                NRRCI=NRRCI+1
                IF (NRRCI.GT.NREC) GOTO 992
                IRRC=NRRCI
                LGPRC(IPLS,IDSC)=IRRC

                IF (NSTORDR >= NRAD) THEN
                  DO 53 J=1,NSBOX
                    ZX=EIONHE/MAX(1.E-5_DP,TEIN(J))
C  rate = [rate coeff <sig v>] times [electr. density], 1/s per ion
c    1.96e-14*sqrt(EionHe/Ry) = 3.5487E-14
                    TABRC1(IRRC,J)=3.5487E-14*ZX**1.5/(ZX+0.35)*DEIN(J)
C  maxw. electron energy loss rate due to recombination
c                   corsum=0._dp  ! old default: 1.5*Te
C  correction due to energy dependence in rec. cross-section
C  corsum=d(ln<sig v>)/d(ln Te)
c  corsum approx -0.5 for Te --> 0
c  corsum approx  0.0 for Te approx 11.5
c  corsum approx +1.0 for Te --> infinity
                    corsum=(-0.5_dp*zx+0.35)/(zx+0.35)
                    EELRC1(IRRC,J)=-(1.5+CORSUM)*TEIN(J)*TABRC1(IRRC,J)
   53             CONTINUE
cdr  this setting kk=-2 is confusing. It may work, but
cdr  kk=-2 is already reserved for other "minimal" default processes
                  NREARC(IRRC) = -2
                  NELRRC(IRRC) = -2
                ELSE          ! storage saving mode:
                              ! tabrc1, eelrc1 to be found "on the fly"
                  NREARC(IRRC) = -2
                  NELRRC(IRRC) = -2
                END IF
                IATM1=IATM
                NATPRC(IRRC)=IATM1
                NIOPRC(IRRC)=0
                NPLPRC(IRRC)=0
                NMLPRC(IRRC)=0
C
                MODCOL(6,2,IRRC)=1
                MODCOL(6,4,IRRC)=1
              ENDIF
   54       CONTINUE
C
            NPRCI(IPLS)=IDSC
          ENDIF
C
C  NON-DEFAULT MODEL: 240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
csw check photonic process
            if(iswr(kk)==7) then    ! PH Processes
              idsc=idsc+1
              nrrci=nrrci+1
              IF (NRRCI.GT.NREC) GOTO 992
cdr  here should come call to xstph or xstot
              call EIRENE_XSTRC(ipls,nrc,idsc,nrrci)
              cycle
csw end branch
            ELSEIF (ISWR(KK).EQ.6) THEN  ! RC Processes
C
              FACTKK=FREACP(IPLS,NRC)
              IF (FACTKK.EQ.0.D0) FACTKK=1.
C  RECOMBINATION MODEL FOR BULK IONS
              IDSC=IDSC+1
              NRRCI=NRRCI+1
              IF (NRRCI.GT.NREC) GOTO 992

              IRRC=NRRCI
              LGPRC(IPLS,IDSC)=IRRC
cdr  for notational consistency: here should come a call to routine xstrc,
cdr  for RC type processes
cdr  as already in case of xsecta, xsectm, xsecti, etc...
cdr  There for the corresponding ei,el,cx and pi processes
cdr  this next stuff should go into xstrc.f
              ITYP=EIRENE_IDEZ(ISCD1P(IPLS,NRC),1,3)
              ISPZ=EIRENE_IDEZ(ISCD1P(IPLS,NRC),3,3)

              IF ((ISPZ < 1) .OR. (ISPZ > MAXSPC(ITYP))) GOTO 995
              IF (ITYP.EQ.3) THEN
                NIOPRC(IRRC)=ISPZ
                RMASS2=RMASSI(ISPZ)
              ELSEIF (ITYP.EQ.4) THEN
                NPLPRC(IRRC)=ISPZ
                RMASS2=RMASSP(ISPZ)
              ELSEIF (ITYP.EQ.1) THEN
                NATPRC(IRRC)=ISPZ
                RMASS2=RMASSA(ISPZ)
              ELSEIF (ITYP.EQ.2) THEN
                NMLPRC(IRRC)=ISPZ
                RMASS2=RMASSM(ISPZ)
              ELSEIF (ITYP.EQ.0) THEN
                NPHPRC(IRRC)=ISPZ
                RMASS2=0.
              ENDIF

              ITYP2=EIRENE_IDEZ(ISCD2P(IPLS,NRC),1,3)
              ISPZ2=EIRENE_IDEZ(ISCD2P(IPLS,NRC),3,3)
              IF ((ISCD2P(IPLS,NRC) /= 0) .AND.
     .           ((ISPZ2 < 1) .OR. (ISPZ2 > MAXSPC(ITYP2)))) GOTO 995
              IF (ITYP2.EQ.4) THEN
                NPLPRC_2(IRRC)=ISPZ2
                RMASS2_2      =RMASSP(ISPZ2)
              ELSE
                RMASS2_2=0._DP
              ENDIF
C  CHECK MASS CONSERVATION
              IF (REACDAT(KK)%NOSEC < 3) THEN
                IF (RMASSP(IPLS).NE.(RMASS2+RMASS2_2)) GOTO 993
              END IF
C
C  1.) CROSS-SECTION(TE)
C           NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3) --> RATE (1/S) per Ion
C
C  2.A) RATE COEFFICIENT = CONST.
C           TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
              IF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.1) THEN
                IF (NSTORDR >= NRAD) THEN

cdr  lexp should not be set from mod(iftflg), that has completely different meaning !!!!!
                  LEXP = .NOT. (MOD(IFTFLG(KK,2),100) == 10)
                  DO J=1,NSBOX
!pb                 IF (LGVAC(J,IPLS)) CYCLE
cdr  a density-independent rate can exist also in a vacuum cell.
cdr  e.g. spontaneous emission of a line, also treated as "recombination" event
cdr       by analogy.
                    IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                    COU = EIRENE_RATE_COEFF(KK,J,TEINL(J),0._DP,LEXP,0)
                    TABRC1(IRRC,J)=COU*FACTKK
                    IF (IFTFLG(KK,2) < 100)
     .                TABRC1(IRRC,J)=TABRC1(IRRC,J)*DEIN(J)
                  END DO
                  NREARC(IRRC) = KK
                ELSE
C  DO NOT STORE DATA, BUT COMPUTE THEM WHEN NEEDED
                  NREARC(IRRC) = KK
                END IF
                MODCOL(6,2,IRRC)=1
C             ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
              ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
                IF (NSTORDR >= NRAD) THEN
                  DO J=1,NSBOX
!pb                 IF (LGVAC(J,IPLS)) CYCLE
                    IF (LGVAC(J,NPLS+1).AND.IFTFLG(KK,2) < 100) CYCLE
                    COU = EIRENE_RATE_COEFF(KK,J,TEINL(J),PLS(J),
     .                   .TRUE.,1)
                    TABRC1(IRRC,J)=COU*FACTKK
                    IF (IFTFLG(KK,2) < 100)
     .                TABRC1(IRRC,J)=TABRC1(IRRC,J)*DEIN(J)
                  END DO

                  NREARC(IRRC) = KK
                ELSE
C  DO NOT STORE DATA, BUT COMPUTE THEM WHEN NEEDED
                  NREARC(IRRC) = KK
                END IF
                MODCOL(6,2,IRRC)=1
              ENDIF  ! (MODCLF(KK),3,5) options

              FACRRC(IRRC,1) = FACTKK
              FACRRC(IRRC,2) = LOG(FACTKK)
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ELECTRON ENERGY LOSS RATE  eV/s per ion
C  flags: NELRRC, and for storage saving mode: additionally JELRRC
C
              NSERC5=EIRENE_IDEZ(ISCDEP(IPLS,NRC),5,5)

              IF (NSERC5.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATE COEFF.
                IF (NSTORDR >= NRAD) THEN
                  DO 101 J=1,NSBOX
                    EELRC1(IRRC,J)=-EELECP(IPLS,NRC)*TABRC1(IRRC,J)
  101             CONTINUE
                  NELRRC(IRRC) = 0
                ELSE
                  NELRRC(IRRC) = 0
                  JELRRC(IRRC) = -1
                  EELRC1(IRRC,1)=EELECP(IPLS,NRC)
                END IF
                MODCOL(6,4,IRRC)=1

              ELSEIF (NSERC5.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ELECTRON = -1.5*TE*RATE COEFF.
                IF (NSTORDR >= NRAD) THEN
                  DO 102 J=1,NSBOX
                    EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
  102             CONTINUE
                  NELRRC(IRRC) = 0
                ELSE
                  NELRRC(IRRC) = 0
                  JELRRC(IRRC) = -2
                END IF
                MODCOL(6,4,IRRC)=1
C
              ELSEIF (NSERC5.EQ.3) THEN

                KREAD=NINT(EELECP(IPLS,NRC))
                IF ((KREAD < 1) .OR. (KREAD > NREACI)) GOTO 996
                MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
c  special treatment in case bremsstrahlung is contained in energy loss rate
c  as e.g. the case in ADAS ADF11- PRB files
                LADAS = EIRENE_IS_RTCEW_TAB2D(KREAD)
C  4.C)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >= NRAD) THEN
                    DO J = 1, NSBOX
                      IF (LGVAC(J,NPLS+1)) CYCLE
C   CAREFUL: EELRC1 IS TO BE TAKEN NEGATIVE, IF IT IS A LOSS!
                      EELRC1(IRRC,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,J,
     .                               TEINL(J),0._DP,.TRUE.,0)
                      EELRC1(IRRC,J)=-EELRC1(IRRC,J)*DEIN(J)*FACTKK
C  SUBTRACT BREMSSTRAHLUNG, if it was included in recombination energy loss rate
c  (since eelrc1 is taken negative, add the bremsstrahlung)
cnh 28.10.2019
                      IF (LADAS) THEN
                        IF (LGVAC(J,IPLS)) CYCLE
                        IF (NCHRGP(IPLS)==0) THEN
                          BREMS = 0._DP
                        ELSE
c                         Charge
                          IF(ZIIN(IPLS,J).NE.ZVAC) THEN
                            Z = ZIIN(IPLS,J)
                          ELSE
                            Z = DBLE(NCHRGP(IPLS))
                          ENDIF
                          BREMS =EIRENE_BREMS(TEIN(J),DEIN(J),Z)/ELCHA  !eV/s/ion
                        END IF
                        EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                      ENDIF
c  bremsstrahlung correction done.

                    END DO
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                  ELSE  ! storage saving mode: only eelrc1(irrc,1)
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                  END IF
                  MODCOL(6,4,IRRC)=1
C  4.D)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE,EBEAM)
C               ELSEIF (MODC.EQ.2) THEN
C        IRRELEVANT
C                 MODCOL(6,4,IRRC)=2
C  4.E)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE,NE), eV/s/ion
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    DO J = 1, NSBOX
                      IF (LGVAC(J,NPLS+1)) CYCLE
C  change logical from false to true, to avoid log(erate), with erate negative
C  as it may result from internal CR code H_COL,...., when used with delpot=0.0
                      EELRC1(IRRC,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,J,
cdr  .                               TEINL(J),PLS(J),.FALSE.,1)
     .                               TEINL(J),PLS(J),.TRUE.,1)
                      EELRC1(IRRC,J)=-EELRC1(IRRC,J)*DEIN(J)*FACTKK
cdr  old code, for log(e_rate) return. Not possible with h_colrad, due to sign change
cdr                   EEMX=MAX(-100._DP,EELRC1(IRRC,J)+DEINL(J))+FCTKKL
cdr                   EELRC1(IRRC,J)=-EXP(EEMX)

C  SUBTRACT BREMSSTRAHLUNG, if it was included in recombination energy loss rate
c  (since eelrc1 is taken negative, add the bremsstrahlung)
                      IF (LADAS) THEN
                        IF (LGVAC(J,IPLS)) CYCLE
                        IF (NCHRGP(IPLS)==0) THEN
                          BREMS = 0._DP
                        ELSE
c                         Charge
                          IF(ZIIN(IPLS,J).NE.ZVAC) THEN
                            Z = ZIIN(IPLS,J)
                          ELSE
                            Z = DBLE(NCHRGP(IPLS))
                          ENDIF
                          BREMS =EIRENE_BREMS(TEIN(J),DEIN(J),Z)/ELCHA  !eV/s/ion
                        END IF
                        EELRC1(IRRC,J) = EELRC1(IRRC,J) + BREMS
                      ENDIF
c  bremsstrahlung correction done.

                    END DO  ! nsbox
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  ELSE  ! STORAGE SAVING MODE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                  END IF

                  MODCOL(6,4,IRRC)=1
                ENDIF   ! MODC =1, 2 or 3

                FACRRC(IRRC,1) = FACTKK
                FACRRC(IRRC,2) = LOG(FACTKK)
C  EELRC1: NEGATIVE SIGN: LOSS FOR ELECTRONS
C
C  SHIFT ELECTRON COOLING RATE BY DELE * TABRC
c  DELE= +IONISATION POTENTIAL TURNS A RADIATION LOSS COMPONENT
C         INTO ELECTRON ENERGY LOSS/GAIN (SIGN CHANGE POSSIBLE)
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    DO 110 J=1,NSBOX
                      EELRC1(IRRC,J)=EELRC1(IRRC,J)+
     .                               DELE*TABRC1(IRRC,J)
  110               CONTINUE
c  STORAGE SAVING MODE AND DELPOT NE 0.0
c                 ELSE  ! ??
                  END IF
C
                ENDIF   ! DELPOT
              ENDIF  ! NSERC5
            ELSE
              GOTO 997
            ENDIF  ! iswr(kk)
C
   82     CONTINUE
          NPRCI(IPLS)=IDSC
C
C  NO MODEL DEFINED
        ELSE
          NPRCI(IPLS)=0
        ENDIF
C
        NPRCIM(IPLS)=NPRCI(IPLS)-1
        LGPRC(IPLS,0)=NPRCI(IPLS)
C
        IF (TRCAMD) THEN

          CALL EIRENE_MASBOX
     .    ('BULK ION SPECIES IPLS = '//TEXTS(NSPAMI+IPLS))
          CALL EIRENE_LEER(1)
C
          IF (LGPRC(IPLS,0).EQ.0) THEN

            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO RECOMBINATION '
            CALL EIRENE_LEER(1)
          ELSE
            DO 220 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              WRITE (iunout,*) 'RECOMBINATION NO. IRRC= ',IRRC
              WRITE (iunout,*) 'RECOMBINATION INTO SPECIES:'
              IION3=NIOPRC(IRRC)
              IF (IION3.NE.0) WRITE (iunout,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC(IRRC)
              IF (IPLS3.NE.0) WRITE (iunout,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC(IRRC)
              IF (IATM3.NE.0) WRITE (iunout,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC(IRRC)
              IF (IMOL3.NE.0) WRITE (iunout,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
              IPHOT3=NPHPRC(IRRC)
              IF (IPHOT3.NE.0) WRITE (iunout,*) 'PHOTON  IPHOT= ',
     .                                     TEXTS(IPHOT3)
C  and, possibly, a second secondary
              IION3=NIOPRC_2(IRRC)
              IF (IION3.NE.0) WRITE (iunout,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC_2(IRRC)
              IF (IPLS3.NE.0) WRITE (iunout,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC_2(IRRC)
              IF (IATM3.NE.0) WRITE (iunout,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC_2(IRRC)
              IF (IMOL3.NE.0) WRITE (iunout,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
              IPHOT3=NPHPRC_2(IRRC)
              IF (IPHOT3.NE.0) WRITE (iunout,*) 'PHOTON  IPHOT= ',
     .                                     TEXTS(IPHOT3)
C
C             WRITE (iunout,*) 'ELECTRONS: PELPRC,EELRC1'
C             IF (NSTORDR >= NRAD) THEN
C               WRITE (iunout,*) 'EL      ',1.,EELRC1(IRRC,1)
C             ELSE
C               WRITE (iunout,*) 'EL      ',1.,FEELRC1(IRRC,1)
C             END IF

              CALL EIRENE_LEER(1)
              WRITE (IUNOUT,*) 'COLLISION MODEL: '
              WRITE (iunout,*) 'PROCESS NO. KK ',NREARC(IRRC)
              WRITE (IUNOUT,*) 'MODCOL         ',
     .                  MODCOL(6,1,IRRC),MODCOL(6,2,IRRC),
     .                  MODCOL(6,3,IRRC),MODCOL(6,4,IRRC)
              WRITE (IUNOUT,'(1X,A15,1(1PE12.4))') 'SCALING FACTOR ',
     .                     FACRRC(IRRC,1)
              CALL EIRENE_LEER(1)
  220       CONTINUE   !irrc for ipls
          ENDIF
          CALL EIRENE_LEER(1)

        ENDIF  !trcamd
C
 1000 CONTINUE

      DEALLOCATE (PLS)
      DEALLOCATE (TEPLS)
C
      RETURN
C
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED'
      WRITE (iunout,*) 'NREC TOO SMALL, CHECK PARAMETER STATEMENTS'
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED'
      WRITE (iunout,*) 'MASS CONSERVATION VIOLATED, IPLS,IRRC ',
     .                  IPLS,IRRC
      CALL EIRENE_EXIT_OWN(1)
  994 CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTP.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE'
      WRITE (iunout,*) 'IPLS = ',IPLS
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED'
      CALL EIRENE_EXIT_OWN(1)
  995 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED'
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: EXIT CALLED'
      WRITE (iunout,*)
     .  'WRONG REACTION INDEX SPECIFIED FOR KREAD IN REACTION KK'
      WRITE (iunout,*) 'KK ',KK
      WRITE (IUNOUT,*) 'KREAD ',KREAD
      CALL EIRENE_EXIT_OWN(1)
  997 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTP: ISCDE FLAG'
      WRITE (iunout,*) 'IRRC, EFLAG, KK, ISWR ',IRRC,NSERC5,KK,ISWR(KK)
      CALL EIRENE_EXIT_OWN(1)
C
      END SUBROUTINE EIRENE_XSECTP
