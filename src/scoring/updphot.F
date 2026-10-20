C  27.6.05 updphot: iadd removed
C  21.01.06: photon background for test atoms: removed
C  18.04.06: test ions and atoms: synchronized
C            bug fix: V0_para  --> parmom_0 for elastic momentum source
C                                  contribution from atoms.
C  10.01.07: parallel momentum exchange tallies MAPL_VEC, MMPL_VEC, MIPL_VEC
C            included as default EIRENE tallies.
C            Before these tallies have been updated in problem
C            specific section UPTCOP, as COPV tallies.
C  12.02.07: Add user-supplied B field and plasma flow option indpro=8
C            to evaluation of parallel momentum sources,
C            Do not use BV_VEC, PMOM_VEC arrays in this case, because they
C            may not have been initialized in subr. PLASMA_DERIV
C            for these options.
C            Use vsig_parp und val_parp instead.
C  25.04.07 update of tallies because of PI reactions revised
C  07.08.07 collision estimators vollstaendig fuer atom, mol und iion.
C           entries: atm, mol, ion voll synchronisiert.
C  28.8.07: esigpi(...,4) --> PL, esigpi(...,5)--> EL
c  oct.14:  some intermediate scoring of additional tally ADDV removed
c  06.08.15 arguments added to vecusr
c  24.08.15 comments and documention wrt. BGK collision treatment
cdr dec.15: tracklength estimators for heavy test particle post-collision energies
cdr         in PI processes added. For A, M, I incident test particles.
cdr dec.15: further corrections, lea --> leio, and other logical flags for turning on-off estimators

cdr nov.15: tracklength estimators for eapl,empl,eipl: species ipl resolved.
cdr apr.16: bug fix J.Lore re index in lgiel. This part of code is still unused,
cdr         so no effect on any result.  Few further comments corrected
cdr fix ph7 PHV_LGPHOT --> lgphph
cdr feb.23: prepare for merging into subr. update (unified treatment)


C
      SUBROUTINE EIRENE_UPDPHOT (XSTOR2,XSTORV2,IFLAG)
C
C ESTIMATORS ARE UPDATED FOR EACH TRACK TAKING T/VEL SEC.
C T (CM) IS STORED ON CLPD ARRAY FOR ONE OR MORE CELLS, THAT HAVE
C BEEN CROSSED WITHOUT COLLISION.

C
C  NCOU: NUMBER OF PIECES OF TRACK IN DIFFERENT CELLS SCORED IN THIS PRESENT CALL (BUT FIXED NRCELL)
C     I: INDIVIDUAL TRACK, I=1,NCOU
C  IRDO: TRACK IS IN (FINE) GEOMETRY CELL IRDO (=NRCELL+NUPC(I)*NR1P2+NBLCKA)
C  IRD:  ESTIMATORS ARE UPDATED IN (COARSE) SCORING CELL IRD  (=NCLTAL(IRDO))
C
C  IFLAG: CURRENTLY ONLY USED FOR PHOTON TALLIES, TO AVOID CANCELLATION OF TERMS

C  IFLAG=1:
C  IFLAG=2:
C  IFLAG=3:
C  IFLAG=4: CALLED FROM WITHIN STATIC LOOP  (PATH LENGTH SET TO MFP), OR CALLED AT POINT OF COLLISION
C  IFLAG=5:

C  SPECIAL TREATMENT OF "BGK" COLLISIONS (= ELASTIC COLLISIONS WITH VIRTUAL BACKGROUND SPECIES)
C
C  A) NPBGK..(ITEST) : IF GT 0, THE CORRESPONDING PARTICLE (IATM, IMOL OR IION) IS A SO-CALLED "BGK" SPECIES
C                      IF, ADDITIONALLY, LBGKV = T, THEN ADDITIONAL BGK TALLIES ARE SCORED VIA A CALL TO UPTBGK
C  B) SIGBGK         : TOTAL RATE OF BGK TYPE COLLISIONS. INCIDENT ATOM AND ITS ENERGY IS NOT LOST
C  C) NPBGKP (IPLS,1): IREL ELASTIC COLLISION CONTRIBUTIONS WITH BULK COLLISION PARTNERS WITH NPBGKP(IPLS,1)>0
C          ARE NOT INCLUDED IN SOURCE/SINK TALLIES.

C          IN CASE OF EAPL THIS IS IMPORTANT, IN ORDER NOT TO MIX ENERGY SOURCES FOR REAL BACKGROUND
C          IONS WITH ENERGY SOURCES FOR VIRTUAL BACKGROUND "IONS"  (MISSING SPECIES INDEX)
C          BUT: CURRENTLY MISSING IN EAAT: CONTRIBUTIONS OF ENERGY EXCHANGE DUE TO BGK COLLISIONS
C          (BOTH SOURCE (DUE TO C) AND SINK (DUE TO B)


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_CSDVI
      USE EIRMOD_COMXS
      USE EIRMOD_CCONA
      USE EIRMOD_PHOTON
      USE EIRMOD_CINIT
      USE EIRMOD_CLOGAU
      USE EIRMOD_UPTUSR

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN OUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                            XSTORV2(NSTORV,N2ND+N3RD)
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: WTRSIG, DIST, WTR, WTRE0, WV, CNDYNPH, WTRV
      INTEGER :: IRD, I, IRDO, INUM
C     INTEGER :: NPBGK
C SECONDARY SPECIES IDENTIFIERS
      INTEGER :: IAT1,IAT2,IML1,IML2,IIO1,IIO2,IPH1,IPH2,IPL1,IPL2
C PH PROCESSES
      INTEGER :: IXPH,IRPH,IPHPH
C    .          ,UPDF        ! out, something for stim. emiss ?
      EXTERNAL :: EIRENE_UPTCOP

C
C  ESTIMATORS FOR PHOTONS
C

cdr: summer 2016
cdr  this entire section has not been looked at (nor used) since long.
cdr  before running photons again: check code carefully.
C
C
      WV=WEIGHT/VEL
c     NPBGK=NPBGKPH(IPHOT)
C
      IF (NADVI.GT.0) CALL EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
      IF (NCPVI.GT.0) CALL EIRENE_UPTCOP(XSTOR2,XSTORV2,WV)
cdr generalise flag NPBGK to mean: model collision term for bi-linear collision
cdr only in this case: set background radiation intensity profiles.
cdr   IF (NPBGK.GT.0) CALL ....
C
      IF (IUPDTE == 2) RETURN

      CNDYNPH = EV_TO_ERG/CLIGHT
C
      DO 131 I=1,NCOU
        DIST=CLPD(I)
        WTR=WV*DIST
        WTRE0=WTR*E0
        WTRV=WTR*VEL*CNDYNPH
        IRDO=NRCELL+NUPC(I)*NR1P2+NBLCKA
        IRD=NCLTAL(IRDO)

C  FOR STANDARD DEVIATION: INDICATE CELLS THAT HAVE BEEN MET BY THE PRESENT MC HISTORY
        IF (IMETCL(IRD) == 0) THEN
          NCLMT = NCLMT+1
          ICLMT(NCLMT) = IRD
          IMETCL(IRD) = NCLMT
        END IF
C
C  PARTICLE, MOMENTUM AND ENERGY DENSITY ESTIMATORS
C
        IF (LEDENPH) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          EDENPH(IPHOT,IRD)=EDENPH(IPHOT,IRD)+WTRE0
        ENDIF
        IF (LPDENPH) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          PDENPH(IPHOT,IRD)=PDENPH(IPHOT,IRD)+WTR
        ENDIF
        IF (LEDENPH.OR.LPDENPH) LMETSP(IPHOT)=.TRUE.

        IF (LVXDENPH) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          VXDENPH(IPHOT,IRD)=VXDENPH(IPHOT,IRD)+WTRV*VELX
        ENDIF
        IF (LVYDENPH) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          VYDENPH(IPHOT,IRD)=VYDENPH(IPHOT,IRD)+WTRV*VELY
        ENDIF
        IF (LVZDENPH) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          VZDENPH(IPHOT,IRD)=VZDENPH(IPHOT,IRD)+WTRV*VELZ
        ENDIF
        IF (LVXDENPH.OR.LVYDENPH.OR.LVZDENPH) LMETSP(IPHOT)=.TRUE.
C
C  ESTIMATORS FOR SOURCES AND SINKS
C  NEGATIVE SIGN MEANS: LOSS FOR PARTICLES
C  POSITIVE SIGN MEANS: GAIN FOR PARTICLES
C
        IF (LGVAC(IRDO,0)) GOTO 131
C
        if (ncou.gt.1) then
          XSTOR(:,:) = XSTOR2(:,:,I)
          XSTORV(:)  = XSTORV2(:,I)
        endif
C
C  PRE-COLLISION RATES, ASSUME: TEST PARTICLES (AND THEIR ENERGY) ARE LOST
C
cdr     IF ((LAST_EVENT%IFLAG == 1) .AND.
cdr past event option turned off. 
cdr biased results in 8 cell tests
cdr but perfect balances in 1 cell tests. ???
cdr i.e. something still not correct with this option.
        IF ((LAST_EVENT%IFLAG == -1) .AND.  ! turns this option off.
     .      (LAST_EVENT%NCELL == IRD)) THEN
cdr  particle was born, is now on first flight
cdr  particle track is still in same cell

! use collision estimator for first cell ("brick") along the track.
! In case UPDATE is called at a collision point: sample 1 (score the whole weight)
cdr This perfectly cancels the source rate. This is an "idle" Monte Carlo particle.
! In case UPDATE is called at any other instance (no collision in cell IRD): score with a zero weight (do not score)
          IF ((IFLAG == 4).OR.(IFLAG == 5)) THEN
            IF (LPPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)-WEIGHT
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPHT2(IPHOT,IPHOT)=PPHPHT2(IPHOT,IPHOT)-WEIGHT
                LMETSP2(1:NPHOT,0:NPHOT) => LMETSP(NTS_IPH+1:NTS_PHPH)
                LMETSP2(IPHOT,0)=.TRUE.
                LMETSP2(IPHOT,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPHT(IPHOT+IPHOT*NPHOT,IRD) =
     .           PPHPHT(IPHOT+IPHOT*NPHOT,IRD)-WEIGHT
                LMETSP(NTS_IPH+IPHOT)=.TRUE.
                LMETSP(NTS_IPH+IPHOT+IPHOT*NPHOT)=.TRUE.
#endif
              END IF
            ENDIF
            IF (LEPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              EPHPHT(IRD)=EPHPHT(IRD)-WEIGHT*E0
            ENDIF
          ELSE
!  nothing to be done, sample a 0
          END IF

        ELSE
! use tracklength estimator
          WTRSIG=WTR*(SIGTOT-SIGBGK)
          IF (LPPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)-WTRSIG
            IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PPHPHT2(IPHOT,IPHOT)=PPHPHT2(IPHOT,IPHOT)-WTRSIG
              LMETSP2(1:NPHOT,0:NPHOT) => LMETSP(NTS_IPH+1:NTS_PHPH)
              LMETSP2(IPHOT,0)=.TRUE.
              LMETSP2(IPHOT,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PPHPHT(IPHOT+IPHOT*NPHOT,IRD) =
     .         PPHPHT(IPHOT+IPHOT*NPHOT,IRD)-WTRSIG
              LMETSP(NTS_IPH+IPHOT)=.TRUE.
              LMETSP(NTS_IPH+IPHOT+IPHOT*NPHOT)=.TRUE.
#endif
            END IF
          ENDIF
          IF (LEPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EPHPHT(IRD)=EPHPHT(IRD)-WTRSIG*E0
          ENDIF
        END IF
C
C..........................................................................
C
C  PHOTONIC (PH) CONTRIBUTION for photons
cdr  in analogy to CX processes. better: analogy to PI ? wg. stim emission ?
C
        IF (LGPHPH(IPHOT,0,0).EQ.0) GOTO 133
C  DEFAULT TRACKLENGTH ESTIMATOR
        DO IPHPH=1,NPHPHI(IPHOT)
          IRPH=LGPHPH(IPHOT,IPHPH,0)
          IPLS=LGPHPH(IPHOT,IPHPH,1)
cdr  check IRPH=0 in initialisation phase only once
cdr       if(IRPH == 0) cycle

cdr       UPDF=LGPHPH(IPHOT,IPHPH,4)
          LOGPLS(IPLS,ISTRA)=.TRUE.
C
          WTRSIG=WTR*SIGVPH(IRPH)
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE-COLLISION RATES HERE
C
          IF (IESTPH(IRPH,1).NE.0) THEN
            IF (LPPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)+WTRSIG
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPHT2(IPHOT,IPHOT)=PPHPHT2(IPHOT,IPHOT)+WTRSIG
                LMETSP2(1:NPHOT,0:NPHOT) => LMETSP(NTS_IPH+1:NTS_PHPH)
                LMETSP2(IPHOT,0)=.TRUE.
                LMETSP2(IPHOT,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPHT(IPHOT+IPHOT*NPHOT,IRD) =
     .           PPHPHT(IPHOT+IPHOT*NPHOT,IRD)+WTRSIG
                LMETSP(NTS_IPH+IPHOT)=.TRUE.
                LMETSP(NTS_IPH+IPHOT+IPHOT*NPHOT)=.TRUE.
#endif
              END IF
            ENDIF
cdr         if(updf==1) PPHPHT(IPHOT,IRD)=PPHPHT(IPHOT,IRD)+WTRSIG !prob. wrong
          ELSE
C
C  PRE-COLLISION RATES, BULK IONS
C
cdr  if(ipls > 0) then
cdr  do this check in initialisation, only once
            IF (LPPHPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PPHPL(IPLS,IRD)=PPHPL(IPLS,IRD)-WTRSIG
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPL2(1:NPLS,0:NPHOT) => PPHPL(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPL2(IPLS,IPHOT)=PPHPL2(IPLS,IPHOT)-WTRSIG
                LMETSP2(1:NPLS,0:NPHOT) => LMETSP(NTS_PHPH+1:NTS_PPH)
                LMETSP2(IPLS,0)=.TRUE.
                LMETSP2(IPLS,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PPHPL(IPLS+IPHOT*NPLS,IRD) =
     .           PPHPL(IPLS+IPHOT*NPLS,IRD)-WTRSIG
                LMETSP(NTS_PHPH+IPLS)=.TRUE.
                LMETSP(NTS_PHPH+IPLS+IPHOT*NPLS)=.TRUE.
#endif
              END IF
              LMETSP(NSPAMI+IPLS)=.TRUE.
            END IF
C
C  POST-COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)
C
cdr-test    goto 133
!pb         IF (N1STPH(IRPH,3).NE.0) THEN
!pb N1STPH(IRPH,3) does not include bulk
C  FIRST SECONDARY:
              IF (N1STPH(IRPH,1).EQ.1) THEN
                IAT1=N1STPH(IRPH,2)
                INUM=N1STPH(IRPH,3)
                LOGATM(IAT1,ISTRA)=.TRUE.
                IF (LPPHAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHAT(IAT1,IRD)= PPHAT(IAT1,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHAT2(1:NATM,0:NPHOT) => PPHAT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHAT2(IAT1,IPHOT)=PPHAT2(IAT1,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NATM,0:NPHOT) => LMETSP(NTS_PI+1:NTS_APH)
                    LMETSP2(IAT1,0)=.TRUE.
                    LMETSP2(IAT1,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHAT(IAT1+IPHOT*NATM,IRD) =
     .               PPHAT(IAT1+IPHOT*NATM,IRD)+WTRSIG*INUM
                    LMETSP(NTS_PI+IAT1)=.TRUE.
                    LMETSP(NTS_PI+IAT1+IPHOT*NATM)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPH+IAT1)=.TRUE.
                END IF
              ELSEIF (N1STPH(IRPH,1).EQ.2) THEN
                IML1=N1STPH(IRPH,2)
                INUM=N1STPH(IRPH,3)
                LOGMOL(IML1,ISTRA)=.TRUE.
                IF (LPPHML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHML(IML1,IRD)= PPHML(IML1,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHML2(1:NMOL,0:NPHOT) => PPHML(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHML2(IML1,IPHOT)=PPHML2(IML1,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NMOL,0:NPHOT) => LMETSP(NTS_APH+1:NTS_MPH)
                    LMETSP2(IML1,0)=.TRUE.
                    LMETSP2(IML1,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHML(IML1+IPHOT*NMOL,IRD) =
     .               PPHML(IML1+IPHOT*NMOL,IRD)+WTRSIG*INUM
                    LMETSP(NTS_APH+IML1)=.TRUE.
                    LMETSP(NTS_APH+IML1+IPHOT*NMOL)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPA+IML1)=.TRUE.
                END IF
              ELSEIF (N1STPH(IRPH,1).EQ.3) THEN
                IIO1=N1STPH(IRPH,2)
                INUM=N1STPH(IRPH,3)
                LOGION(IIO1,ISTRA)=.TRUE.
                IF (LPPHIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHIO(IIO1,IRD)= PPHIO(IIO1,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHIO2(1:NION,0:NPHOT) => PPHIO(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHIO2(IIO1,IPHOT)=PPHIO2(IIO1,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NION,0:NPHOT) => LMETSP(NTS_MPH+1:NTS_IPH)
                    LMETSP2(IIO1,0)=.TRUE.
                    LMETSP2(IIO1,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHIO(IIO1+IPHOT*NION,IRD) =
     .               PPHIO(IIO1+IPHOT*NION,IRD)+WTRSIG*INUM
                    LMETSP(NTS_MPH+IIO1)=.TRUE.
                    LMETSP(NTS_MPH+IIO1+IPHOT*NION)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPAM+IIO1)=.TRUE.
                END IF
              ELSEIF (N1STPH(IRPH,1).EQ.4) THEN
                IPL1=N1STPH(IRPH,2)
C               INUM=N1STPH(IRPH,3)
                INUM=1
                LOGPLS(IPL1,ISTRA)=.TRUE.
                IF (LPPHPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHPL(IPL1,IRD)= PPHPL(IPL1,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHPL2(1:NPLS,0:NPHOT) => PPHPL(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPL2(IPL1,IPHOT)=PPHPL2(IPL1,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NPLS,0:NPHOT)=>LMETSP(NTS_PHPH+1:NTS_PPH)
                    LMETSP2(IPL1,0)=.TRUE.
                    LMETSP2(IPL1,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPL(IPL1+IPHOT*NPLS,IRD) =
     .               PPHPL(IPL1+IPHOT*NPLS,IRD)+WTRSIG*INUM
                    LMETSP(NTS_PHPH+IPL1)=.TRUE.
                    LMETSP(NTS_PHPH+IPL1+IPHOT*NPLS)=.TRUE.
#endif
                  END IF
csw added updf check (stim.em)
cdr: not ready
cdr  stim emission: am besten: 2 secondaries in group 1. hier jedoch:
cdr  dazu PI-process vervollstandigen.
cdr  test photon + ipls --> ipl2 (bulk particle)
cdr  if this bulk is a photon (same as iphot), dann noch eins rauf auf tally
                  LMETSP(NSPAMI+IPL1)=.TRUE.
                END IF
csw added branch
              ELSEIF (N1STPH(IRPH,1).EQ.0) then
                iph1=N1STPH(IRPH,2)
c               inum=N1STPH(IRPH,3)
                inum=1                   
cdr  test iph1 > 0 only once, in initialisation. here: removed
                if ((inum > 0) .and. (iph1 > 0)) then
                  logphot(iph1,istra)=.true.
                  IF (LPPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPHT(iph1,ird)=PPHPHT(iph1,ird)+wtrsig*inum
                    IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                      PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                      PPHPHT2(iph1,iphot)=PPHPHT2(iph1,iphot)+
     .                                    wtrsig*inum
                      LMETSP2(1:NPHOT,0:NPHOT)=>
     .                       LMETSP(NTS_IPH+1:NTS_PHPH)
                      LMETSP2(IPH1,0)=.TRUE.
                      LMETSP2(IPH1,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPHT(IPH1+IPHOT*NPHOT,IRD) =
     .               PPHPHT(IPH1+IPHOT*NPHOT,IRD)+WTRSIG*INUM
                    LMETSP(NTS_IPH+IPH1)=.TRUE.
                    LMETSP(NTS_IPH+IPH1+IPHOT*NPHOT)=.TRUE.
#endif
                    END IF
                    LMETSP(IPH1)=.TRUE.
                  END IF
                end if
              ENDIF
!pb         END IF

!pb         IF (N2NDPH(IRPH,3).NE.0) THEN
!pb N2NDPH(IRPH,3) does not include bulk
C  SECOND SECONDARY:
              IF (N2NDPH(IRPH,1).EQ.1) THEN
                IAT2=N2NDPH(IRPH,2)
c               INUM=N2NDPH(IRPH,3)
                inum=1   
                LOGATM(IAT2,ISTRA)=.TRUE.
                IF (LPPHAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHAT(IAT2,IRD)= PPHAT(IAT2,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHAT2(1:NATM,0:NPHOT) => PPHAT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHAT2(IAT2,IPHOT)=PPHAT2(IAT2,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NATM,0:NPHOT)=>LMETSP(NTS_PI+1:NTS_APH)
                    LMETSP2(IAT2,0)=.TRUE.
                    LMETSP2(IAT2,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHAT(IAT2+IPHOT*NATM,IRD) =
     .               PPHAT(IAT2+IPHOT*NATM,IRD)+WTRSIG*INUM
                    LMETSP(NTS_PI+IAT2)=.TRUE.
                    LMETSP(NTS_PI+IAT2+IPHOT*NATM)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPH+IAT2)=.TRUE.
                END IF
              ELSEIF (N2NDPH(IRPH,1).EQ.2) THEN
                IML2=N2NDPH(IRPH,2)
c               INUM=N2NDPH(IRPH,3)
                inum=1 
                LOGMOL(IML2,ISTRA)=.TRUE.
                IF (LPPHML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHML(IML2,IRD)= PPHML(IML2,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHML2(1:NMOL,0:NPHOT) => PPHML(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHML2(IML2,IPHOT)=PPHML2(IML2,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NMOL,0:NPHOT)=>LMETSP(NTS_APH+1:NTS_MPH)
                    LMETSP2(IML2,0)=.TRUE.
                    LMETSP2(IML2,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHML(IML2+IPHOT*NMOL,IRD) =
     .               PPHML(IML2+IPHOT*NMOL,IRD)+WTRSIG*INUM
                    LMETSP(NTS_APH+IML2)=.TRUE.
                    LMETSP(NTS_APH+IML2+IPHOT*NMOL)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPA+IML2)=.TRUE.
                END IF
              ELSEIF (N2NDPH(IRPH,1).EQ.3) THEN
                IIO2=N2NDPH(IRPH,2)
c               INUM=N2NDPH(IRPH,3)
                inum=1 
                LOGION(IIO2,ISTRA)=.TRUE.
                IF (LPPHIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHIO(IIO2,IRD)= PPHIO(IIO2,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHIO2(1:NION,0:NPHOT) => PPHIO(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHIO2(IIO2,IPHOT)=PPHIO2(IIO2,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NION,0:NPHOT)=>LMETSP(NTS_MPH+1:NTS_IPH)
                    LMETSP2(IIO2,0)=.TRUE.
                    LMETSP2(IIO2,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHIO(IIO2+IPHOT*NION,IRD) =
     .               PPHIO(IIO2+IPHOT*NION,IRD)+WTRSIG*INUM
                    LMETSP(NTS_MPH+IIO2)=.TRUE.
                    LMETSP(NTS_MPH+IIO2+IPHOT*NION)=.TRUE.
#endif
                  END IF
                  LMETSP(NSPAM+IIO2)=.TRUE.
                END IF
              ELSEIF (N2NDPH(IRPH,1).EQ.4) THEN
                IPL2=N2NDPH(IRPH,2)
C               INUM=N2NDPH(IRPH,3)
                INUM=1
                LOGPLS(IPL2,ISTRA)=.TRUE.
                IF (LPPHPL) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  PPHPL(IPL2,IRD)= PPHPL(IPL2,IRD)+WTRSIG*INUM
                  IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                    PPHPL2(1:NPLS,0:NPHOT) => PPHPL(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPL2(IPL2,IPHOT)=PPHPL2(IPL2,IPHOT)+WTRSIG*INUM
                    LMETSP2(1:NPLS,0:NPHOT)=>LMETSP(NTS_PHPH+1:NTS_PPH)
                    LMETSP2(IPL2,0)=.TRUE.
                    LMETSP2(IPL2,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPL(IPL2+IPHOT*NPLS,IRD) =
     .               PPHPL(IPL2+IPHOT*NPLS,IRD)+WTRSIG*INUM
                    LMETSP(NTS_PHPH+IPL2)=.TRUE.
                    LMETSP(NTS_PHPH+IPL2+IPHOT*NPLS)=.TRUE.
#endif
                  END IF
csw added updf check (stim.em)
cdr  stim emission: am besten: 2 secondaries in group 2.
cdr  dazu PI-process vervollstandigen.
cdr  test photon + ipls --> ipl2 (bulk particle)
cdr  if this bulk is a photon (same as iphot), dann noch eins rauf auf tally
                  LMETSP(NSPAMI+IPL2)=.TRUE.
                END IF
csw added branch
              ELSEIF (N2NDPH(IRPH,1).EQ.0) then
                iph2=N2NDPH(IRPH,2)
c               inum=N2NDPH(IRPH,3)
                inum=1  
cdr test iph2 > 0 removed, to be done only once in initialisation
                if ((inum > 0) .and. (iph2 > 0)) then
                  logphot(iph2,istra)=.true.
                  IF (LPPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    PPHPHT(iph2,ird)=PPHPHT(iph2,ird)+wtrsig*inum
                    IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                      PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,IRD)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                      PPHPHT2(iph2,IPHOT)=PPHPHT2(iph2,IPHOT)+
     .                                    wtrsig*inum
                      LMETSP2(1:NPHOT,0:NPHOT)=>
     .                       LMETSP(NTS_IPH+1:NTS_PHPH)
                      LMETSP2(IPH2,0)=.TRUE.
                      LMETSP2(IPH2,IPHOT)=.TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                      PPHPHT(IPH2+IPHOT*NPHOT,IRD) =
     .                 PPHPHT(IPH2+IPHOT*NPHOT,IRD)+WTRSIG*INUM
                      LMETSP(NTS_IPH+IPH2)=.TRUE.
                      LMETSP(NTS_IPH+IPH2+IPHOT*NPHOT)=.TRUE.
#endif
                    END IF
                    LMETSP(iph2)=.true.
                  END IF
                END IF
              ENDIF
!pb         ENDIF
          ENDIF
C
C  PARTICLE ESTIMATORS DONE. NEXT: ENERGY ESTIMATORS
C
C  COLLISION ESTIMATOR IN SUBR. COLLIDE ?
C  COMPENSATE PRE-COLLISION RATES HERE
C
          IF (LEPHPHT.AND.(IESTPH(IRPH,3).NE.0)) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            EPHPHT(IRD)=EPHPHT(IRD)          +WTRSIG*E0
cdr         if(updf==1) EPHPHT(IRD)=EPHPHT(IRD)+WTRSIG*E0 ! verm. falsch
          ELSE
C
C  PRE-COLLISION RATES, BULK IONS
C
cdr if(ipls > 0) then : this check only once in initialisation. removed
cdr         IF (LEPHEL) EPHPL(IRD)    =EPHPL(IRD)  -WTRSIG*E0  vermutl. falsch
cdr  gibt es schon ESIGPH ? ist dann IRPH das richtige argument ?
cdr         IF (LEPHPL) EPHPL(IRD)    =EPHPL(IRD)  -WTRSIG*ESIGPH(IRPH)
C
C  POST-COLLISION RATES, ALL SECONDARIES (TEST AND BULK PARTICLES)

!dr       IF (N1STPH(IRPH,3).NE.0) THEN
!dr N1STPH(IRPH,3) does not include bulk
C  FIRST SECONDARY:
cdr: erstmal raus. vermutlich muss hier die energie des
cdr  neuen schwerteilchens stehen, nicht die E0 des IPHOT.
cdr         IF (N1STPH(IRPH,1).EQ.1) THEN
cdr           IAT1=N1STPH(IRPH,2)
cdr           LOGATM(IAT1,ISTRA)=.TRUE.
cdr           IF (LEPHAT) EPHAT(IRD)     = EPHAT(IRD)  +WTRSIG*E0
C           ELSEIF (N1STPH(IRPH,1).EQ.2) THEN
C             IML1=N1STPH(IRPH,2)
C             LOGMOL(IML1,ISTRA)=.TRUE.
C             IF (LEPHML) EPHML(IRD)     = EPHML(IRD)     +WTRSIG*E0
cdr         ELSEIF (N1STPH(IRPH,1).EQ.3) THEN
cdr           IIO1=N1STPH(IRPH,2)
cdr           LOGION(IIO1,ISTRA)=.TRUE.
cdr           IF (LEPHIO) EPHIO(IRD)     = EPHIO(IRD)     +WTRSIG*E0
cdr         ELSEIF (N1STPH(IRPH,1).EQ.4) THEN
cdr           IPL1=N1STPH(IRPH,2)
cdr           LOGPLS(IPL1,ISTRA)=.TRUE.
cdr           IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)   +WTRSIG*E0
csw added updf check (stim.em)
cdr           if(LEPHPL.AND.(updf==1.and.phv_is_plsphot(ipl1)>0)) then
cdr              EPHPL(IRD)     = EPHPL(IRD)+WTRSIG*E0
cdr           endif
csw added branch
cdr: photon scattering, e0_in = e0_out, also: delta scattering,
cdr  mit eventuell teilchenmultiplikation
cdr         ELSEIF (N1STPH(IRPH,1).EQ.0) THEN
              IF (N1STPH(IRPH,1).EQ.0) THEN
                IPH1=N1STPH(IRPH,2)
                INUM=N1STPH(IRPH,2)
cdr  if(iph1 > 0) then abfrage hier raus, nur in initialisation
                LOGPHOT(IPH1,ISTRA)=.TRUE.
                IF (LEPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                  EPHPHT(IRD)=EPHPHT(IRD) +WTRSIG*E0*INUM
                ENDIF
              ENDIF
!dr         ENDIF

!dr     IF (N2NDPH(IRPH,3).NE.0) THEN
!dr N2NDPH(IRPH,3) does not include bulk
C  SECOND SECONDARY:
cdr : gleiches problem wie oben E0_out ? bei test photon impact mit E0 ?
cdr         IF (N2NDPH(IRPH,1).EQ.1) THEN
cdr           IAT2=N2NDPH(IRPH,2)
cdr           LOGATM(IAT2,ISTRA)=.TRUE.
cdr           IF (LEPHAT) EPHAT(IRD)     = EPHAT(IRD)     +WTRSIG*E0
C           ELSEIF (N2NDPH(IRPH,1).EQ.2) THEN
C             IML2=N2NDPH(IRPH,2)
C             LOGMOL(IML2,ISTRA)=.TRUE.
C             IF (LEPHML) EPHML(IRD)     = EPHML(IRD)     +WTRSIG*E0
cdr         ELSEIF (N2NDPH(IRPH,1).EQ.3) THEN
cdr           IIO2=N2NDPH(IRPH,2)
cdr           LOGION(IIO2,ISTRA)=.TRUE.
cdr           IF (LEPHIO) EPHIO(IRD)     = EPHIO(IRD)     +WTRSIG*E0
cdr         ELSEIF (N2NDPH(IRPH,1).EQ.4) THEN
cdr           IPL2=N2NDPH(IRPH,2)
cdr           LOGPLS(IPL2,ISTRA)=.TRUE.
cdr           IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)     +WTRSIG*E0
csw added updf check (stim.em)
cdr           if(updf==1.and.phv_is_plsphot(ipl2)>0) then
cdr              IF (LEPHPL) EPHPL(IRD)     = EPHPL(IRD)     +WTRSIG*E0
cdr           endif
csw added branch
cdr         ELSEIF (N2NDPH(IRPH,1).EQ.0) THEN
              IF (N2NDPH(IRPH,1).EQ.0) THEN
                IPH2=N2NDPH(IRPH,2)
c               INUM=N2NDPH(IRPH,3)
                inum=1 
cdr if(iph2 > 0) then  ! dieser test nur in initialisation phase
                if ((inum > 0) .and. (iph2 > 0)) then
                  LOGPHOT(IPH2,ISTRA)=.TRUE.
                  IF (LEPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                    EPHPHT(IRD) = EPHPHT(IRD)+WTRSIG*E0*INUM
                  ENDIF
                end if
              ENDIF
!dr         ENDIF
          ENDIF
C
        ENDDO
  133   CONTINUE
  131 CONTINUE
      RETURN

      END SUBROUTINE EIRENE_UPDPHOT
