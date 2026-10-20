!pb  28.06.06: bug fix for NSTORAM=0 and MODC=1
!pb            FACREA=FACTKK instead of FACREA=log(FACTKK)
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!dr  30.01.07: if lgvac(..,npls+1)  cycle (do not evaluate rates in vacuum)
!pb  20.04.07: allow for third and fourth secondary

!pb  June  07: LHCOL: indicate: direct coupling to col-rad code, rather
!pb            than reading fits or data tables.
cdr  Jan. 2014:
!dr  to be done: also for recombination, and generalize to other species (He,...)
!dr             currently: label H.4 2.1.5 or H.10 2.1.5 are not used in case LHCOL?
!   23.02.14:   nomenclature changed IPL --> IPP to provide consistency with XSTPI.f
!   02.02.15:   ONLY COMMENTS ADDED
!dr Jan   16:   EPLDS: SPECIES INDEX ADDED. Old EPLDS is now EPLEI(..,0,..)

!pb APR   16:   pplds -> pplei
!pb APR   16:   patds -> patei, eatds -> eatei
!pb APR   16:   pmlds -> pmlei, emlds -> emlei
!pb APR   16:   piods -> pioei, eiods -> eioei
!pb APR   16:   pelds -> pelei, eelds -> eelei
!pb MAY   16:   tabds1 -> tabei1
!pb JUL   16:   ehvds1 -> ehvei1


cdr Aug.16  :   minor synchronisation with xstpi.f.
cdr Sept.16 :   Started to implement H.3 rate coeff.
cdr             for high E0, low Te cases. Needs to be added: TABEI3
cdr May 17  :   TABEI1: safety cut-off for TEE at 0.1 eV,
cdr 	        (rather than the eirene default TVAC (=0.02 eV), added in more cases
cdr         :   tbd: to be replaced by a proper Arrhenius form extrapolation
!pb Juli 17 :   LHCOL removed
cdr Sept 18 :   JELREI: flags for electron energy loss, in storage save mode only.
cdr             i.e. only needed in FEELEI1

C
      SUBROUTINE EIRENE_XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,ITHRD,IFRTH,
     .                 EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
!
!  Set rate coefficients and reaction kinetics for e + ISP collisions.
!
!
c   rmass: mass of incident test particle
c   irei:  counting index for this particular electron impact collision
c   isp:   incident test particle species identifier
c   PLS:   precomputed log of electron density


C
C  SET NON-DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMXS
      use EIRMOD_ctrcei, only: trcamd

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, EELEC, FACTKK,CHRDF0
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ITHRD, IFRTH,
     .                       ISCDE, IESTM, KK
      REAL(DP) :: CF(9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL,
     .          EE, TB, TEE,
     .          ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV, COU,
     .          EIRENE_RATE_COEFF,
     .          EIRENE_ENERGY_RATE_COEFF,
     .          DELE,
     .          FP1(6),FP2(6)
      INTEGER :: MODC, KREAD, J, IPP, IP,
     .           III, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      REAL(DP), PARAMETER :: TMINL=-2.3_DP
      type(poly_data), pointer :: rp
      type(fit_forms), pointer :: rt
      EXTERNAL :: EIRENE_DBL_POLY, EIRENE_LEER, EIRENE_EXIT_OWN,
     .            EIRENE_RATE_COEFF, EIRENE_ENERGY_RATE_COEFF

      ITYP=EIRENE_IDEZ(IFRST,1,3)
      INUM=EIRENE_IDEZ(IFRST,2,3)
      ISPE=EIRENE_IDEZ(IFRST,3,3)
C
C ACCUMULATED MASS OF SECONDARIES: ACCMAS (AMU)
      ACCMAS=0.D0
      ACCMSA=0.D0
      ACCMSM=0.D0
      ACCMSI=0.D0
      ACCMSP=0.D0
      ACCINV=0.D0
      ACCINA=0.D0
      ACCINM=0.D0
      ACCINI=0.D0
      ACCINP=0.D0
C
      ICOUNT=1
   85 CONTINUE

      IF ((ISPE < 1) .OR. (ISPE > MAXSPC(ITYP))) GOTO 994

!  ACCMAS: accumulated mass of all secondaries (all types)
!  ACCINV: accumulated inverse mass of all secondaries (all types)

!  ACCMSA: accumulated mass of ATOMIC secondaries (type ITYP=1)
!  ACCINA: accumulated inverse mass of ATOMIC secondaries (type ITYP=1)
!  analogously for molecule, test ion and bulk secondaries
      IF (ITYP.EQ.1) THEN
        IAT=ISPE
        IAA=NSPH+IAT
        PATEI(IREI,IAT)=PATEI(IREI,IAT)+INUM
        P2ND(IREI,IAA)=P2ND(IREI,IAA)+INUM
        ACCMAS=ACCMAS+INUM*RMASSA(IAT)
        ACCMSA=ACCMSA+INUM*RMASSA(IAT)
        ACCINV=ACCINV+INUM/RMASSA(IAT)
        ACCINA=ACCINA+INUM/RMASSA(IAT)
        EATEI(IREI,IAT,1)=RMASSA(IAT)
        EATEI(IREI,IAT,2)=1./RMASSA(IAT)
      ELSEIF (ITYP.EQ.2) THEN
        IML=ISPE
        IMM=NSPA+IML
        PMLEI(IREI,IML)=PMLEI(IREI,IML)+INUM
        P2ND(IREI,IMM)=P2ND(IREI,IMM)+INUM
        ACCMAS=ACCMAS+INUM*RMASSM(IML)
        ACCMSM=ACCMSM+INUM*RMASSM(IML)
        ACCINV=ACCINV+INUM/RMASSM(IML)
        ACCINM=ACCINM+INUM/RMASSM(IML)
        EMLEI(IREI,IML,1)=RMASSM(IML)
        EMLEI(IREI,IML,2)=1./RMASSM(IML)
      ELSEIF (ITYP.EQ.3) THEN
        IIO=ISPE
        III=NSPAM+IIO
        PIOEI(IREI,IIO)=PIOEI(IREI,IIO)+INUM
        P2ND(IREI,III)=P2ND(IREI,III)+INUM
        ACCMAS=ACCMAS+INUM*RMASSI(IIO)
        ACCMSI=ACCMSI+INUM*RMASSI(IIO)
        ACCINV=ACCINV+INUM/RMASSI(IIO)
        ACCINI=ACCINI+INUM/RMASSI(IIO)
        EIOEI(IREI,IIO,1)=RMASSI(IIO)
        EIOEI(IREI,IIO,2)=1./RMASSI(IIO)
      ELSEIF (ITYP.EQ.4) THEN
        IPP=ISPE
        PPLEI(IREI,IPP)=PPLEI(IREI,IPP)+INUM
        ACCMAS=ACCMAS+INUM*RMASSP(IPP)
        ACCMSP=ACCMSP+INUM*RMASSP(IPP)
        ACCINV=ACCINV+INUM/RMASSP(IPP)
        ACCINP=ACCINP+INUM/RMASSP(IPP)
        EPLEI(IREI,IPP,1)=RMASSP(IPP)
        EPLEI(IREI,IPP,2)=1./RMASSP(IPP)
      ENDIF
C
      IF (ISCND.NE.0.AND.ICOUNT.EQ.1) THEN
        ITYP=EIRENE_IDEZ(ISCND,1,3)
        INUM=EIRENE_IDEZ(ISCND,2,3)
        ISPE=EIRENE_IDEZ(ISCND,3,3)
        ICOUNT=2
        GOTO 85
      ELSEIF (ITHRD.NE.0.AND.ICOUNT.EQ.2) THEN
        ITYP=EIRENE_IDEZ(ITHRD,1,3)
        INUM=EIRENE_IDEZ(ITHRD,2,3)
        ISPE=EIRENE_IDEZ(ITHRD,3,3)
        ICOUNT=3
        GOTO 85
      ELSEIF (IFRTH.NE.0.AND.ICOUNT.EQ.3) THEN
        ITYP=EIRENE_IDEZ(IFRTH,1,3)
        INUM=EIRENE_IDEZ(IFRTH,2,3)
        ISPE=EIRENE_IDEZ(IFRTH,3,3)
        ICOUNT=4
        GOTO 85
      ENDIF
C
      IF (ABS(ACCMAS-RMASS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTEI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT TEST SPECIES ',TEXTS(ISP)
        WRITE (iunout,*)
     .    'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
        CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
      DO IAT=1,NATMI
        EATEI(IREI,IAT,1)=EATEI(IREI,IAT,1)/ACCMAS
        EATEI(IREI,IAT,2)=EATEI(IREI,IAT,2)/ACCINV
      ENDDO
      EATEI(IREI,0,1)=ACCMSA/ACCMAS
      EATEI(IREI,0,2)=ACCINA/ACCINV
      DO IML=1,NMOLI
        EMLEI(IREI,IML,1)=EMLEI(IREI,IML,1)/ACCMAS
        EMLEI(IREI,IML,2)=EMLEI(IREI,IML,2)/ACCINV
      ENDDO
      EMLEI(IREI,0,1)=ACCMSM/ACCMAS
      EMLEI(IREI,0,2)=ACCINM/ACCINV
      DO IIO=1,NIONI
        EIOEI(IREI,IIO,1)=EIOEI(IREI,IIO,1)/ACCMAS
        EIOEI(IREI,IIO,2)=EIOEI(IREI,IIO,2)/ACCINV
      ENDDO
      EIOEI(IREI,0,1)=ACCMSI/ACCMAS
      EIOEI(IREI,0,2)=ACCINI/ACCINV
      DO IPP=1,NPLSI
        EPLEI(IREI,IPP,1)=EPLEI(IREI,IPP,1)/ACCMAS
        EPLEI(IREI,IPP,2)=EPLEI(IREI,IPP,2)/ACCINV
      ENDDO
      EPLEI(IREI,0,1)=ACCMSP/ACCMAS
      EPLEI(IREI,0,2)=ACCINP/ACCINV
C
      CHRDIF=CHRDF0
      DO 83 IIO=1,NIONI
        CHRDIF=CHRDIF+PIOEI(IREI,IIO)*NCHRGI(IIO)
   83 CONTINUE

      DO 84 IP=1,NPLSI
        CHRDIF=CHRDIF+PPLEI(IREI,IP)*NCHRGP(IP)
   84 CONTINUE
      PELEI(IREI)=PELEI(IREI)+CHRDIF
C
C
C  1.) CROSS-SECTION(TE) : NOT NEEDED
C
C
C..................................................................
C  2.) RATE COEFFICIENT  (CM**3/S) * ELECTRON DENSITY (CM**-3)
C..................................................................
C

C
C  2.A) RATE COEFFICIENT = CONST.
C     TO BE WRITTEN

      IF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.1) THEN
C  2.B) RATE COEFFICIENT(TE)  (FIXED, OR: LOW DENSITY (CORONA) RATE COEFF.
C                              E0 FIXED (E.G. E0=0.0)
        IF (NSTORDR >= NRAD) THEN
C  RATE: (1/S) =
C  RATE COEFFICIENT: (CM^3/S) * DENSITY (CM^3)
          DO J=1,NSBOX
            IF (LGVAC(J,NPLS+1)) CYCLE
            TEE=TEINL(J)
cdr  safety cut-off at TE= 0.1 eV. (TVAC=0.02)
            TEE = max(tminl,TEE)
            COU = EIRENE_RATE_COEFF(KK,J,TEE,0._DP,.TRUE.,0)
            TABEI1(IREI,J)=COU*FACTKK
C  IS TABEI1 A RATE COEFFICIENT OR ALREADY A RATE ?
            IF (IFTFLG(KK,2) < 100)
     .        TABEI1(IREI,J)=TABEI1(IREI,J)*DEIN(J)
          END DO
          NREAEI(IREI) = KK
        ELSE ! NOT SUFFICIENT STORAGE ON TABEI1
          NREAEI(IREI) = KK
        ENDIF
        MODCOL(1,2,IREI)=1

      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C       NEND=9
C  TO BE WRITTEN
        if (.true.) goto 996

        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)
          rt => reacdat(kk)%rtc
          fp1(1:3) = rt%fp1l
          fp1(4:6) = rt%fp1r
          fp2(1:3) = rt%fp2b
          fp2(4:6) = rt%fp2t
          DO J=1,NSBOX
            IF (LGVAC(J,NPLS+1)) CYCLE
              TEE=TEINL(J)
cdr  safety cut-off at TE= 0.1 eV. (note: TVAC=0.02)
              TEE = max(tminl,TEE)
c  evaluate 2 parametric fit,
c  collapse this to a one parameter fit CF for EB dependence, evaluated at TEE.
              rp => reacdat(KK)%rtc%poly
              call EIRENE_dbl_poly (rp%dblpol,tee,0._dp,cou,cf,
     .               rt%rc1min, rt%rc1max, fp1, rt%jfex1mn, rt%jfex1mx,
     .               rt%rc2min, rt%rc2max, fp2, rt%jfex2mn, rt%jfex2mx,
     .               trcamd)
cdr  not ready, tabei1 --> tabei3 to be done.
C             TABEI3(IREI,J,1:9) = CF(1:9)
C             TABEI3(IREI,J,1)=TABEI3(IREI,J,1)+DEINL(J)+FCTKKL
          END DO
        ELSE ! NOT SUFFICIENT STORAGE ON TABEI3
C  STORAGE SAVE MODE NOT READY FOR THIS OPTION ??

        ENDIF
        MODCOL(1,2,IREI)=2

      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)

          DO J=1,NSBOX
            IF (LGVAC(J,NPLS+1)) CYCLE
            TEE=TEINL(J)
cdr  safety cut-off at Te= 0.1 eV. (note: TVAC=0.02)
            TEE = max(tminl,TEE)
cdr  safety cut-off at ne= 1e8 cm**-3 already in PLS(..) from calling program. DVAC=1.0e2)
            COU = EIRENE_RATE_COEFF(KK,J,TEE,PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DEINL(J)
            TB=MAX(-100._DP,TB)
            TABEI1(IREI,J)=EXP(TB)
          END DO
          NREAEI(IREI) = KK
        ELSE ! NOT SUFFICIENT STORAGE ON TABEI1
C  WHAT DO WE DO IN CASE NSTORDR < NRAD  ?
          NREAEI(IREI) = KK
        ENDIF
        MODCOL(1,2,IREI)=1 !  indicate: rate coefficient as fct.
                           !  of local plasma conditions only
      ELSE
        GOTO 996
      ENDIF

      FACREI(IREI,1) = FACTKK
      FACREI(IREI,2) = LOG(FACTKK)
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ENERGY RATES
C
C  4.A: ELECTRON ENERGY LOSS RATES
C
      EFLAG=EIRENE_IDEZ(ISCDE,5,5)
      IF (EFLAG.EQ.0) THEN
C  4.A1) ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATE COEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 101 J=1,NSBOX
                  EELEI1(IREI,J)=EELEC
  101           CONTINUE
                NELREI(IREI)=0
              ELSE ! NOT SUFFICIENT STORAGE ON EELEI1
                NELREI(IREI)=0
                JELREI(IREI)=-1
                EELEI1(IREI,1)=EELEC
              END IF
              MODCOL(1,4,IREI)=1

      ELSEIF (EFLAG.EQ.1) THEN
C  4.A2) ENERGY LOSS RATE OF IMP. ELECTRON = 1.5*TE*RATE COEFF
              IF (NSTORDR >= NRAD) THEN
                DO 103 J=1,NSBOX
                  IF (LGVAC(J,NPLS+1)) CYCLE
                  EELEI1(IREI,J)=-1.5*TEIN(J)
  103           CONTINUE
                NELREI(IREI)=0
              ELSE ! NOT SUFFICIENT STORAGE ON EELEI1
                JELREI(IREI)=-2
                NELREI(IREI)=0
              END IF
              MODCOL(1,4,IREI)=1

      ELSEIF (EFLAG.EQ.3) THEN
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE), NO. KREAD
                KREAD=NINT(EELEC)
                IF ((KREAD < 1) .OR. (KREAD > NREACI)) GOTO 998
                MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    DO 102 J=1,NSBOX
                      IF (LGVAC(J,NPLS+1)) CYCLE
                      EELEI1(IREI,J)=-EIRENE_ENERGY_RATE_COEFF(KREAD,J,
     .                                TEINL(J),
     .                                0._DP,.TRUE.,0)*DEIN(J)*FACTKK/
     .                                (TABEI1(IREI,J)+EPS60)
  102               CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ENDIF
                  MODCOL(1,4,IREI)=1
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A5) ENERGY LOSS RATE OF IMP. ELECTRON = EN.-WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    DO J = 1, NSBOX
                      IF (LGVAC(J,NPLS+1)) CYCLE
cdr return erate, rather than ln(erate), because a recombination
cdr reaction (with delpot .ne. 0) may be used as EI reaction too,
cdr e.g. when trace ions are followed and recombine.
cdr But delpot .ne.0 with ln(erate) causes trouble with internal CR models.
cdr Had already been taken care of similarly in xstrc.f
                      EE = EIRENE_ENERGY_RATE_COEFF(KREAD,J,TEINL(J),
cdr  .                                              PLS(J),.FALSE.,1) ! to be removed
     .                                              PLS(J),.TRUE.,1)
crc check if EE is a rate or a rate coefficient
crc if it is a rate coefficient, it should be multiplied by DEIN(J)
                      IF (IFTFLG(KREAD,2) < 100) EE = EE*DEIN(J)
                      EELEI1(IREI,J)=-EE*FACTKK/
     .                               (TABEI1(IREI,J)+EPS60)
cdr                   EE = MAX(-100._DP,EE+FCTKKL+DEINL(J))     ! to be removed
cdr                   EELEI1(IREI,J)=-EXP(EE)/(TABEI1(IREI,J)+EPS60)   ! to be removed
                    END DO
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  END IF
                  MODCOL(1,4,IREI)=1
                ENDIF
                FACREI(IREI,1)=FACTKK
                FACREI(IREI,2)=LOG(FACTKK)
C  NEGATIVE SIGN: LOSS FOR ELECTRONS
C  POSITIVE SIGN: GAIN FOR ELECTRONS
C  SHIFT ELECTRON COOLING RATE BY DELE * TABEI
c  DELE= -IONISATION POTENTIAL TURNS EELEI INTO A RADIATION LOSS COMPONENT ONLY
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    DO J=1,NSBOX
                      EELEI1(IREI,J)=EELEI1(IREI,J)+DELE
                    END DO
                  END IF
                ENDIF
      ELSE
        GOTO 997
      ENDIF
C
C  4.B: HEAVY PARTICLE ENERGY GAIN RATE, "KINETIC ENERGY RELEASE" KER
C
      EFLAG=EIRENE_IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.B1)  RATE = CONST.*RATE COEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVEI1(IREI,J)=EHEAVY
  201     CONTINUE
          NHVREI(IREI)=0
        ELSE
          NHVREI(IREI)=0
          EHVEI1(IREI,1)=EHEAVY
        END IF

C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION

      ELSEIF (EFLAG.EQ.3) THEN
C  4.B3)  ENERGY RATE = EN.-WEIGHTED RATE(TE)
        KREAD=NINT(EHEAVY)
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        FACREI(IREI,1)=FACTKK
        FACREI(IREI,2)=LOG(FACTKK)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
              IF (LGVAC(J,NPLS+1)) CYCLE
              EHVEI1(IREI,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,J,TEINL(J),
     .             0._DP,.TRUE.,0)*DEIN(J)*FACTKK/(TABEI1(IREI,J)+EPS60)
  202       CONTINUE
            NHVREI(IREI)=KREAD
          ELSE
            NHVREI(IREI)=KREAD
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTEI: MODC=EFLAG '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ELSE
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEI(IREI,1)=EIRENE_IDEZ(IESTM,1,3)
      IESTEI(IREI,2)=EIRENE_IDEZ(IESTM,2,3)
      IESTEI(IREI,3)=EIRENE_IDEZ(IESTM,3,3)
C
      IF (IESTEI(IREI,1).NE.0) THEN
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART. BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        CALL EIRENE_LEER(1)
        IESTEI(IREI,1)=0
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM. BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        CALL EIRENE_LEER(1)
        IESTEI(IREI,2)=0
      ENDIF
      RETURN
C
C-----------------------------------------------------------------------
C
  994 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: EXIT CALLED'
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (iunout,*) IREI,KK
      CALL EIRENE_EXIT_OWN(1)
  997 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (iunout,*) 'IREI, EFLAG ',IREI,EFLAG
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: INVALID KREAD'
      WRITE (iunout,*) IREI,KREAD
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_XSTEI
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIRENE_XSTEI_1(IREI)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IREI
      INTEGER :: IA, IPP, IO, ISPZ, IAT, IM, IIO, IML
      REAL(DP) :: P2N
C
C  SET TOTAL NUMBER OF SECONDARIES BY TYPE OF SECONDARY: P..DS(IREI,0)
C  AND
C  CONVERT SECONDARY SPECIES DISTRIBUTION P2ND(IREI)  INTO
C  CUMULATIVE DISTRIBUTION (NOT YET NORMALIZED, THIS IS DONE BELOW)

C  ATOM SECONDARIES
      DO 510 IAT=1,NATMI
        IA=NSPH+IAT
        PATEI(IREI,0)=PATEI(IREI,0)+
     +                      PATEI(IREI,IAT)
        P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                      P2ND(IREI,IA)
  510 CONTINUE
C  MOLECULE SECONDARIES
      DO 520 IML=1,NMOLI
        IM=NSPA+IML
        PMLEI(IREI,0)=PMLEI(IREI,0)+
     +                      PMLEI(IREI,IML)
        P2ND(IREI,IM)=P2ND(IREI,IM-1)+
     +                      P2ND(IREI,IM)
  520 CONTINUE
C  TEST ION SECONDARIES
      DO 530 IIO=1,NIONI
        IO=NSPAM+IIO
        PIOEI(IREI,0)=PIOEI(IREI,0)+
     +                      PIOEI(IREI,IIO)
        P2ND(IREI,IO)=P2ND(IREI,IO-1)+
     +                      P2ND(IREI,IO)
  530 CONTINUE
C  BULK SECONDARIES (NOT ON P2ND)
      DO 540 IPP=1,NPLSI
        PPLEI(IREI,0)=PPLEI(IREI,0)+
     +                      PPLEI(IREI,IPP)
  540 CONTINUE
C
C  TOTAL NUMBER OF SECONDARIES
      P2NEI(IREI)=PATEI(IREI,0)+PMLEI(IREI,0)+
     .            PIOEI(IREI,0)

C  FINALY: NORMALIZE SECONDARY TEST PARTICLE SPECIES DISTRIBUTION P2ND
C          SUCH THAT IT BECOMES A CUMULATIVE SAMPLING DISTRIBUTION
C          FOR TEST PARTICLE SECONDARIES
C          NORMALIZATION DOES NOT EXTEND OVER SECONDARY BULK PARTICLES
      P2N=P2ND(IREI,NSPAMI)
      DO 550 ISPZ=NSPH+1,NSPAMI
        IF (P2N.GT.0.D0)
     .  P2ND(IREI,ISPZ)=P2ND(IREI,ISPZ)/P2N
  550 CONTINUE
C
      RETURN
      END SUBROUTINE EIRENE_XSTEI_1
C
C-----------------------------------------------------------------------
C
      SUBROUTINE EIRENE_XSTEI_2(IREI)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMXS
CTK      use EIRMOD_ctrcei, only: trcamd

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IREI
      REAL(DP) :: EN, EA, EI, EIRENE_FEHVEI1, EIRENE_FEELEI1
      INTEGER :: IPP, IP, IRAD, IO, IIO, IMIN, IMAX, IAT, IA, IML, IM
      EXTERNAL :: EIRENE_LEER, EIRENE_FEHVEI1, EIRENE_FEELEI1
C
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL EIRENE_LEER(1)

      EI=1.D30
      EA=-1.D30
      imin=0
      imax=0
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELEI1(IREI,IRAD)
        ELSE
          EN=EIRENE_FEELEI1(IREI,IRAD)
        END IF
        if (en < ei) imin=irad
        if (en > ea) imax=irad
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
  875 CONTINUE

      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'

      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (iunout,*) 'ELECTRONS: PELEI, CONSTANT ENERGY: EEL'
        WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELEI(IREI),EI
      ELSE
        WRITE (iunout,*)
     .    'ELECTRONS: PELEI, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (iunout,'(1X,A8,3(1PE12.4))') 'EL      ',PELEI(IREI),EI,EA
      ENDIF
c     write (iunout,*) ' imin = ', imin, ' imax = ',imax
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVEI1(IREI,IRAD)
        ELSE
          EN=EIRENE_FEHVEI1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
  876 CONTINUE
      IF (PPLEI(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLEI'
        DO 874 IPP=1,NPLSI
          IP=NSPAMI+IPP
          IF (PPLEI(IREI,IPP).NE.0.D0)
     .      WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLEI(IREI,IPP)
  874   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EPLEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLEI(IREI,0,1),
     .                                 ' * E0 + ',EPLEI(IREI,0,2)*EI
        ELSEIF (EI.NE.1.D30) THEN
          WRITE (iunout,*) 'ENERGY: EPLEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLEI(IREI,0,1),
     .                                 ' * E0 + ',EPLEI(IREI,0,2),
     .                                 ' * EHEAVY'
C  IN CASE OF EI PROCESSES: COM IS SET EQ. E0
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (iunout,*) 'BULK IONS: NONE'
      ENDIF
      CALL EIRENE_LEER(1)
C
      WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NEI(IREI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL EIRENE_LEER(1)
        GOTO 880
      ENDIF
C
      IF (PATEI(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATEI'
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATEI(IREI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATEI(IREI,IAT)
  871   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATEI(IREI,0,1),
     .                                 ' * E0 + ',EATEI(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATEI(IREI,0,1),
     .                                 ' * E0 + ',EATEI(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLEI(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLEI'
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLEI(IREI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLEI(IREI,IML)
  872   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLEI(IREI,0,1),
     .                                 ' * E0 + ',EMLEI(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLEI(IREI,0,1),
     .                                 ' * E0 + ',EMLEI(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIOEI(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIOEI'
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIOEI(IREI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIOEI(IREI,IIO)
  873   CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIOEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIOEI(IREI,0,1),
     .                                 ' * E0 + ',EIOEI(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIOEI'
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIOEI(IREI,0,1),
     .                                 ' * E0 + ',EIOEI(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF

  880 CONTINUE

      CALL EIRENE_LEER(1)

      IF (IESTEI(IREI,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART. BALANCE '
      IF (IESTEI(IREI,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM. BALANCE '
      IF (IESTEI(IREI,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN. BALANCE '
      CALL EIRENE_LEER(1)

      WRITE (IUNOUT,*) 'COLLISION MODEL: '
      WRITE (IUNOUT,*) 'PROCESS NO. KK ',NREAEI(IREI)
      WRITE (IUNOUT,*) 'MODCOL         ',
     .                  MODCOL(1,1,IREI),MODCOL(1,2,IREI),
     .                  MODCOL(1,3,IREI),MODCOL(1,4,IREI)
      WRITE (IUNOUT,'(1X,A15,1(1PE12.4))') 'SCALING FACTOR ',
     .                  FACREI(IREI,1)
      CALL EIRENE_LEER(1)

      RETURN
      END SUBROUTINE EIRENE_XSTEI_2
