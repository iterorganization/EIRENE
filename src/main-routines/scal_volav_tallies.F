!pb  31.07.09: scaling for snapshot tallies corrected
!              only for stationary calculations snapshot tallies
!              need to be scaled by timestep DTIMV
cdr  dec. 15 : species index added in volumetric energy tallies for bulk ions
cdr            eapl,empl,eipl,ephpl,eppl (38, 44, 50, 56 and 84)
cdr  March 18: X.B.: bug fix re scaling of addv, copv, snapv, copv,
cdr                  in case scltal =4 this was erroneously overwritten with 0.
cdr May 19:    scoring of volumetric spectra: ISRFCLL unclear. Use VOL or VOLTAL ?
cdr            scaling with FATM, FMOL,... done elsewhere?

      SUBROUTINE EIRENE_SCAL_VOLAV_TALLIES (ISTR, ZWW, ZW,
     .                               ZVOLIN, ZVOLIW, SCLTAL, N1DIM)
cdr
c  scaling of volume-averaged tallies:
c  zvolin : source strength(amp)/vol(cell volume cm^-3)/elementary-charge(amp/(1/s)),
c           --> e.g. converts cumulated flight times [s] to densities [cm^-3]
c  zvoliw : source strength(amp)/vol(cell volume cm^-3),
c           --> e.g. converts counted events per cell [1] to source rates (amp/cm^-3)
c                or  also: cumulated (flight times/reaction time) [1] to source rates (amp/cm^-3)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CGRID
      USE EIRMOD_COUTAU
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEZ
      USE EIRMOD_COMNNL
      USE EIRMOD_CLOGAU

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTR, N1DIM
      REAL(DP), INTENT(IN) :: ZWW, ZW
      REAL(DP), INTENT(INOUT) :: SCLTAL(N1DIM,*)
      REAL(DP), INTENT(IN) :: ZVOLIN(*), ZVOLIW(*)

      REAL(DP) :: ZFAC, FACDT
      REAL(DP) :: DEL, DELI, ELEFT, ERIGHT
      INTEGER :: IATM, I, J, IMOL, IION, IPHOT, IPLS, IADV, ICLV,
     .           ISNV, ICPV, IBGV, ISPC, ICL,
     .           KK, ICO
#if ( defined(G95) || defined(LEGACYCOMP) )
      INTEGER :: JATM, JMOL, JION, JPHOT
#endif
C
C
C*****VOLUME-AVERAGED TALLIES  220 - 239
C
C  ATOMIC PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 220 IATM=1,NATMI
        IF (LOGATM(IATM,ISTR)) THEN
          DO 221 J=1,NSBOX_TAL
            IF (LPDENA) PDENA(IATM,J)=PDENA(IATM,J)*ZVOLIN(J)
            IF (LEDENA) EDENA(IATM,J)=EDENA(IATM,J)*ZVOLIN(J)
            IF (LPAAT) THEN
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAAT2(1:NATM,0:NATM) => PAAT(:,J)
                PAAT2(IATM,:) = PAAT2(IATM,:) *ZVOLIW(J)
#else
                DO JATM = 0, NATM
                  PAAT(IATM+JATM*NATM,J) =
     .             PAAT(IATM+JATM*NATM,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PAAT(IATM,J) = PAAT(IATM,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPMAT) THEN
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMAT2(1:NATM,0:NMOL) => PMAT(:,J)
                PMAT2(IATM,:) = PMAT2(IATM,:) *ZVOLIW(J)
#else
                DO JMOL = 0, NMOL
                  PMAT(IATM+JMOL*NATM,J) =
     .             PMAT(IATM+JMOL*NATM,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PMAT(IATM,J) = PMAT(IATM,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPIAT) THEN
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIAT2(1:NATM,0:NION) => PIAT(:,J)
                PIAT2(IATM,:) = PIAT2(IATM,:) *ZVOLIW(J)
#else
                DO JION = 0, NION
                  PIAT(IATM+JION*NATM,J) =
     .             PIAT(IATM+JION*NATM,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PIAT(IATM,J) = PIAT(IATM,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPPHAT) THEN
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHAT2(1:NATM,0:NPHOT) => PPHAT(:,J)
                PPHAT2(IATM,:) = PPHAT2(IATM,:) *ZVOLIW(J)
#else
                DO JPHOT = 0, NPHOT
                  PPHAT(IATM+JPHOT*NATM,J) =
     .             PPHAT(IATM+JPHOT*NATM,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
               PPHAT(IATM,J) = PPHAT(IATM,J)*ZVOLIW(J)
              END IF
            END IF
            IF (LPGENA) PGENA(IATM,J)=PGENA(IATM,J)*ZVOLIW(J)
            IF (LEGENA) EGENA(IATM,J)=EGENA(IATM,J)*ZVOLIW(J)
            IF (LVGENA) VGENA(IATM,J)=VGENA(IATM,J)*ZVOLIW(J)
            IF (LPPAT)  PPAT(IATM,J) =PPAT(IATM,J) *ZVOLIW(J)
            IF (LVXDENA) VXDENA(IATM,J)=VXDENA(IATM,J)*ZVOLIN(J)
            IF (LVYDENA) VYDENA(IATM,J)=VYDENA(IATM,J)*ZVOLIN(J)
            IF (LVZDENA) VZDENA(IATM,J)=VZDENA(IATM,J)*ZVOLIN(J)
            IF (LRAEL)  RAEL(IATM,J) =RAEL(IATM,J) *ZVOLIW(J)
  221     CONTINUE
        ENDIF
        SCLTAL(IATM,1)=1
        SCLTAL(IATM,5)=1
        SCLTAL(IATM,10)=3
        SCLTAL(IATM,16)=3
        SCLTAL(IATM,22)=3
        SCLTAL(IATM,28)=3
        SCLTAL(IATM,63)=3
        SCLTAL(IATM,67)=3
        SCLTAL(IATM,71)=3
        SCLTAL(IATM,75)=3
        SCLTAL(IATM,85)=1
        SCLTAL(IATM,89)=1
        SCLTAL(IATM,93)=1
        SCLTAL(IATM,101)=3
  220 CONTINUE
C
C  MOLECULAR PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 222 IMOL=1,NMOLI
        IF (LOGMOL(IMOL,ISTR)) THEN
          DO 223 J=1,NSBOX_TAL
            IF (LPDENM) PDENM(IMOL,J)=PDENM(IMOL,J)*ZVOLIN(J)
            IF (LEDENM) EDENM(IMOL,J)=EDENM(IMOL,J)*ZVOLIN(J)
            IF (LPAML) THEN
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAML2(1:NMOL,0:NATM) => PAML(:,J)
                PAML2(IMOL,:) = PAML2(IMOL,:) *ZVOLIW(J)
#else
                DO JATM = 0, NATM
                  PAML(IMOL+JATM*NMOL,J) =
     .             PAML(IMOL+JATM*NMOL,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PAML(IMOL,J) = PAML(IMOL,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPMML) THEN
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMML2(1:NMOL,0:NMOL) => PMML(:,J)
                PMML2(IMOL,:) = PMML2(IMOL,:) *ZVOLIW(J)
#else
                DO JMOL = 0, NMOL
                  PMML(IMOL+JMOL*NMOL,J) =
     .             PMML(IMOL+JMOL*NMOL,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PMML(IMOL,J) = PMML(IMOL,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPIML) THEN
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIML2(1:NMOL,0:NION) => PIML(:,J)
                PIML2(IMOL,:) = PIML2(IMOL,:) *ZVOLIW(J)
#else
                DO JION = 0, NION
                  PIML(IMOL+JION*NMOL,J) =
     .             PIML(IMOL+JION*NMOL,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PIML(IMOL,J) = PIML(IMOL,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPPHML) THEN
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHML2(1:NMOL,0:NPHOT) => PPHML(:,J)
                PPHML2(IMOL,:) = PPHML2(IMOL,:) *ZVOLIW(J)
#else
                DO JPHOT = 0, NPHOT
                  PPHML(IMOL+JPHOT*NMOL,J) =
     .             PPHML(IMOL+JPHOT*NMOL,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PPHML(IMOL,J)= PPHML(IMOL,J)*ZVOLIW(J)
              END IF
            END IF
            IF (LPGENM) PGENM(IMOL,J)=PGENM(IMOL,J)*ZVOLIW(J)
            IF (LEGENM) EGENM(IMOL,J)=EGENM(IMOL,J)*ZVOLIW(J)
            IF (LVGENM) VGENM(IMOL,J)=VGENM(IMOL,J)*ZVOLIW(J)
            IF (LPPML)  PPML(IMOL,J) =PPML(IMOL,J) *ZVOLIW(J)
            IF (LVXDENM) VXDENM(IMOL,J)=VXDENM(IMOL,J)*ZVOLIN(J)
            IF (LVYDENM) VYDENM(IMOL,J)=VYDENM(IMOL,J)*ZVOLIN(J)
            IF (LVZDENM) VZDENM(IMOL,J)=VZDENM(IMOL,J)*ZVOLIN(J)
            IF (LRMEL)  RMEL(IMOL,J) =RMEL(IMOL,J) *ZVOLIW(J)
  223     CONTINUE
        ENDIF
        SCLTAL(IMOL,2)=1
        SCLTAL(IMOL,6)=1
        SCLTAL(IMOL,11)=3
        SCLTAL(IMOL,17)=3
        SCLTAL(IMOL,23)=3
        SCLTAL(IMOL,29)=3
        SCLTAL(IMOL,64)=3
        SCLTAL(IMOL,68)=3
        SCLTAL(IMOL,72)=3
        SCLTAL(IMOL,76)=3
        SCLTAL(IMOL,86)=1
        SCLTAL(IMOL,90)=1
        SCLTAL(IMOL,94)=1
        SCLTAL(IMOL,102)=3
  222 CONTINUE
C
C  TEST ION PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO 225 IION=1,NIONI
        IF (LOGION(IION,ISTR)) THEN
          DO 226 J=1,NSBOX_TAL
            IF (LPDENI) PDENI(IION,J)=PDENI(IION,J)*ZVOLIN(J)
            IF (LEDENI) EDENI(IION,J)=EDENI(IION,J)*ZVOLIN(J)
            IF (LPAIO) THEN
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAIO2(1:NION,0:NATM) => PAIO(:,J)
                PAIO2(IION,:) = PAIO2(IION,:) *ZVOLIW(J)
#else
                DO JATM = 0, NATM
                  PAIO(IION+JATM*NION,J) =
     .             PAIO(IION+JATM*NION,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PAIO(IION,J) = PAIO(IION,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPMIO) THEN
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMIO2(1:NION,0:NMOL) => PMIO(:,J)
                PMIO2(IION,:) = PMIO2(IION,:) *ZVOLIW(J)
#else
                DO JMOL = 0, NMOL
                  PMIO(IION+JMOL*NION,J) =
     .             PMIO(IION+JMOL*NION,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PMIO(IION,J) = PMIO(IION,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPIIO) THEN
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIIO2(1:NION,0:NION) => PIIO(:,J)
                PIIO2(IION,:) = PIIO2(IION,:) *ZVOLIW(J)
#else
                DO JION = 0, NION
                  PIIO(IION+JION*NION,J) =
     .             PIIO(IION+JION*NION,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PIIO(IION,J) =PIIO(IION,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPPHIO) THEN
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHIO2(1:NION,0:NPHOT) => PPHIO(:,J)
                PPHIO2(IION,:) = PPHIO2(IION,:) *ZVOLIW(J)
#else
                DO JPHOT = 0, NPHOT
                  PPHIO(IION+JPHOT*NION,J) =
     .             PPHIO(IION+JPHOT*NION,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PPHIO(IION,J)=PPHIO(IION,J)*ZVOLIW(J)
              END IF
            END IF
            IF (LPGENI) PGENI(IION,J)=PGENI(IION,J)*ZVOLIW(J)
            IF (LEGENI) EGENI(IION,J)=EGENI(IION,J)*ZVOLIW(J)
            IF (LVGENI) VGENI(IION,J)=VGENI(IION,J)*ZVOLIW(J)
            IF (LPPIO)  PPIO(IION,J) =PPIO(IION,J) *ZVOLIW(J)
            IF (LVXDENI) VXDENI(IION,J)=VXDENI(IION,J)*ZVOLIN(J)
            IF (LVYDENI) VYDENI(IION,J)=VYDENI(IION,J)*ZVOLIN(J)
            IF (LVZDENI) VZDENI(IION,J)=VZDENI(IION,J)*ZVOLIN(J)
            IF (LRIEL)  RIEL(IION,J) =RIEL(IION,J)*ZVOLIW(J)
  226     CONTINUE
        ENDIF
        SCLTAL(IION,3)=1
        SCLTAL(IION,7)=1
        SCLTAL(IION,12)=3
        SCLTAL(IION,18)=3
        SCLTAL(IION,24)=3
        SCLTAL(IION,30)=3
        SCLTAL(IION,65)=3
        SCLTAL(IION,69)=3
        SCLTAL(IION,73)=3
        SCLTAL(IION,77)=3
        SCLTAL(IION,87)=1
        SCLTAL(IION,91)=1
        SCLTAL(IION,95)=1
        SCLTAL(IION,103)=3
  225 CONTINUE
C
C  PHOTON PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C
      DO IPHOT=1,NPHOTI
        IF (LOGPHOT(IPHOT,ISTR)) THEN
          DO J=1,NSBOX_TAL
            IF (LPDENPH) PDENPH(IPHOT,J)=PDENPH(IPHOT,J)*ZVOLIN(J)
            IF (LEDENPH) EDENPH(IPHOT,J)=EDENPH(IPHOT,J)*ZVOLIN(J)
            IF (LPAPHT) THEN
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAPHT2(1:NPHOT,0:NATM) => PAPHT(:,J)
                PAPHT2(IPHOT,:) = PAPHT2(IPHOT,:) *ZVOLIW(J)
#else
                DO JATM = 0, NATM
                  PAPHT(IPHOT+JATM*NPHOT,J) =
     .             PAPHT(IPHOT+JATM*NPHOT,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PAPHT(IPHOT,J) = PAPHT(IPHOT,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPMPHT) THEN
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMPHT2(1:NPHOT,0:NMOL) => PMPHT(:,J)
                PMPHT2(IPHOT,:) = PMPHT2(IPHOT,:) *ZVOLIW(J)
#else
                DO JMOL = 0, NMOL
                  PMPHT(IPHOT+JMOL*NPHOT,J) =
     .             PMPHT(IPHOT+JMOL*NPHOT,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PMPHT(IPHOT,J) = PMPHT(IPHOT,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPIPHT) THEN
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIPHT2(1:NPHOT,0:NION) => PIPHT(:,J)
                PIPHT2(IPHOT,:) = PIPHT2(IPHOT,:) *ZVOLIW(J)
#else
                DO JION = 0, NION
                  PIPHT(IPHOT+JION*NPHOT,J) =
     .             PIPHT(IPHOT+JION*NPHOT,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PIPHT(IPHOT,J) = PIPHT(IPHOT,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPPHPHT) THEN
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPHT2(1:NPHOT,0:NPHOT) => PPHPHT(:,J)
                PPHPHT2(IPHOT,:) = PPHPHT2(IPHOT,:) *ZVOLIW(J)
#else
                DO JPHOT = 0, NPHOT
                  PPHPHT(IPHOT+JPHOT*NPHOT,J) =
     .             PPHPHT(IPHOT+JPHOT*NPHOT,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PPHPHT(IPHOT,J) = PPHPHT(IPHOT,J)*ZVOLIW(J)
              END IF
            END IF
            IF (LPGENPH) PGENPH(IPHOT,J)=PGENPH(IPHOT,J)*ZVOLIW(J)
            IF (LEGENPH) EGENPH(IPHOT,J)=EGENPH(IPHOT,J)*ZVOLIW(J)
            IF (LVGENPH) VGENPH(IPHOT,J)=VGENPH(IPHOT,J)*ZVOLIW(J)
            IF (LPPPHT)  PPPHT(IPHOT,J) =PPPHT(IPHOT,J) *ZVOLIW(J)
            IF (LVXDENPH) VXDENPH(IPHOT,J)=VXDENPH(IPHOT,J)*ZVOLIN(J)
            IF (LVYDENPH) VYDENPH(IPHOT,J)=VYDENPH(IPHOT,J)*ZVOLIN(J)
            IF (LVZDENPH) VZDENPH(IPHOT,J)=VZDENPH(IPHOT,J)*ZVOLIN(J)
          END DO
        ENDIF
        SCLTAL(IPHOT,4)=1
        SCLTAL(IPHOT,8)=1
        SCLTAL(IPHOT,13)=3
        SCLTAL(IPHOT,19)=3
        SCLTAL(IPHOT,25)=3
        SCLTAL(IPHOT,31)=3
        SCLTAL(IPHOT,66)=3
        SCLTAL(IPHOT,70)=3
        SCLTAL(IPHOT,74)=3
        SCLTAL(IPHOT,78)=3
        SCLTAL(IPHOT,88)=1
        SCLTAL(IPHOT,92)=1
        SCLTAL(IPHOT,96)=1
      END DO
C
C  BULK ION PARTICLE SPECIES LOOP FOR THE STRATUM ISTRA
C  ONLY SOURCE RATES, NO DENSITIES
C
      DO 227 IPLS=1,NPLSI
        IF (LOGPLS(IPLS,ISTR)) THEN
          DO 228 J=1,NSBOX_TAL
            IF (LPAPL) THEN
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAPL2(1:NPLS,0:NATM) => PAPL(:,J)
                PAPL2(IPLS,:) = PAPL2(IPLS,:) *ZVOLIW(J)
#else
                DO JATM = 0, NATM
                  PAPL(IPLS+JATM*NPLS,J) =
     .             PAPL(IPLS+JATM*NPLS,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PAPL(IPLS,J) = PAPL(IPLS,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPMPL) THEN
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMPL2(1:NPLS,0:NMOL) => PMPL(:,J)
                PMPL2(IPLS,:) = PMPL2(IPLS,:) *ZVOLIW(J)
#else
                DO JMOL = 0, NMOL
                  PMPL(IPLS+JMOL*NPLS,J) =
     .             PMPL(IPLS+JMOL*NPLS,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PMPL(IPLS,J) = PMPL(IPLS,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPIPL) THEN
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIPL2(1:NPLS,0:NION) => PIPL(:,J)
                PIPL2(IPLS,:) = PIPL2(IPLS,:) *ZVOLIW(J)
#else
                DO JION = 0, NION
                  PIPL(IPLS+JION*NPLS,J) =
     .             PIPL(IPLS+JION*NPLS,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PIPL(IPLS,J) = PIPL(IPLS,J) *ZVOLIW(J)
              END IF
            END IF
            IF (LPPHPL) THEN
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPL2(1:NPLS,0:NPHOT) => PPHPL(:,J)
                PPHPL2(IPLS,:) = PPHPL2(IPLS,:) *ZVOLIW(J)
#else
                DO JPHOT = 0, NPHOT
                  PPHPL(IPLS+JPHOT*NPLS,J) =
     .             PPHPL(IPLS+JPHOT*NPLS,J)*ZVOLIW(J)
                END DO
#endif
              ELSE
                PPHPL(IPLS,J)= PPHPL(IPLS,J)*ZVOLIW(J)
              END IF
            END IF
            IF (LPPPL)  PPPL(IPLS,J) =PPPL(IPLS,J) *ZVOLIW(J)

            IF (LEAPL)  EAPL(IPLS,J) =EAPL(IPLS,J) *ZVOLIW(J)
            IF (LEMPL)  EMPL(IPLS,J) =EMPL(IPLS,J) *ZVOLIW(J)
            IF (LEIPL)  EIPL(IPLS,J) =EIPL(IPLS,J) *ZVOLIW(J)
            IF (LEPHPL) EPHPL(IPLS,J)=EPHPL(IPLS,J)*ZVOLIW(J)
            IF (LEPPL)  EPPL(IPLS,J) =EPPL(IPLS,J) *ZVOLIW(J)
            DO ICO=0,2
              KK=ICO*NPLS
              IF (LMAPL_VEC)  MAPL_VEC(KK+IPLS,J) =
     .                        MAPL_VEC(KK+IPLS,J) *ZVOLIW(J)
              IF (LMMPL_VEC)  MMPL_VEC(KK+IPLS,J) =
     .                        MMPL_VEC(KK+IPLS,J) *ZVOLIW(J)
              IF (LMIPL_VEC)  MIPL_VEC(KK+IPLS,J) =
     .                        MIPL_VEC(KK+IPLS,J) *ZVOLIW(J)
              IF (LMPHPL_VEC) MPHPL_VEC(KK+IPLS,J)=
     .                        MPHPL_VEC(KK+IPLS,J)*ZVOLIW(J)
            ENDDO
  228     CONTINUE
          SCLTAL(IPLS,14)=3
          SCLTAL(IPLS,20)=3
          SCLTAL(IPLS,26)=3
          SCLTAL(IPLS,32)=3
          SCLTAL(IPLS,38)=3
          SCLTAL(IPLS,44)=3
          SCLTAL(IPLS,50)=3
          SCLTAL(IPLS,56)=3
          SCLTAL(IPLS,79)=3
          SCLTAL(IPLS,84)=3
          SCLTAL(IPLS,97)=3
          SCLTAL(IPLS,98)=3
          SCLTAL(IPLS,99)=3
          SCLTAL(IPLS,100)=3
        ENDIF
  227 CONTINUE
C
C  ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES FOR THE STRATUM ISTRA
C  TALLY ADDV; NO. NTALA
C  # is a placeholder for the units of the detector function used
C  for scoring these additional tallies (see manual, block 10)
C
      IF (LADDV) THEN
        DO 230 IADV=1,NADVI
          IF (IADVE(IADV).EQ.1) THEN
C  SCALE # PER VOLUME
            DO 231 J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*ZVOLIN(J)
  231       CONTINUE
            SCLTAL(IADV,NTALA)=1
          ELSEIF (IADVE(IADV).EQ.2) THEN
C  SCALE # PER CELL
            DO 232 J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*ZW
  232       CONTINUE
            SCLTAL(IADV,NTALA)=2
C  SCALE # AMP/S PER VOLUME
          ELSEIF (IADVE(IADV).EQ.3) THEN
            DO 233 J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*ZVOLIW(J)
  233       CONTINUE
            SCLTAL(IADV,NTALA)=3
C  SCALE # AMP/S PER CELL
          ELSEIF (IADVE(IADV).EQ.4) THEN
            DO 234 J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*ZWW
  234       CONTINUE
            SCLTAL(IADV,NTALA)=4
          ELSE
C  DO NOT SCALE AT ALL
            SCLTAL(IADV,NTALA)=0
          ENDIF
  230   CONTINUE
      END IF
C
C  ADDITIONAL COLLISION-ESTIMATED TALLIES FOR THE STRATUM ISTRA
C  TALLY COLV; NO. NTALC
C
      IF (LCOLV) THEN
        DO 235 ICLV=1,NCLVI
          IF (ICLVE(ICLV).EQ.1) THEN
C  SCALE # PER VOLUME
            DO 236 J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*ZVOLIN(J)
  236       CONTINUE
            SCLTAL(ICLV,NTALC)=1
          ELSEIF (ICLVE(ICLV).EQ.2) THEN
C  SCALE # PER CELL
            DO 237 J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*ZW
  237       CONTINUE
            SCLTAL(ICLV,NTALC)=2
          ELSEIF (ICLVE(ICLV).EQ.3) THEN
C  SCALE # AMP/S PER VOLUME
            DO 238 J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*ZVOLIW(J)
  238       CONTINUE
            SCLTAL(ICLV,NTALC)=3
          ELSEIF (ICLVE(ICLV).EQ.4) THEN
C  SCALE # AMP/S PER CELL
            DO 239 J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*ZWW
  239       CONTINUE
            SCLTAL(ICLV,NTALC)=4
          ELSE
C  DO NOT SCALE AT ALL
            SCLTAL(ICLV,NTALC)=0
          ENDIF
  235   CONTINUE
      END IF
C
C  ADDITIONAL SNAPSHOT-ESTIMATED TALLIES, FOR THE STRATUM ISTRA
C  TALLY SNAPV; NO. NTALT,
C  FIRST: SNAPV=SNAPV*DTIMV, THEN: SCALING
C
      IF (LSNAPV) THEN
        FACDT = 1._DP
        IF (NTMSTP < 0) FACDT = DTIMV
        DO 245 ISNV=1,NSNVI
          IF (ISNVE(ISNV).EQ.1) THEN
C  SCALE # PER VOLUME
            DO 246 J=1,NSBOX_TAL
              SNAPV(ISNV,J)=SNAPV(ISNV,J)*FACDT*ZVOLIN(J)
  246       CONTINUE
            SCLTAL(ISNV,NTALT)=1
          ELSEIF (ISNVE(ISNV).EQ.2) THEN
C  SCALE # PER CELL
            DO 247 J=1,NSBOX_TAL
              SNAPV(ISNV,J)=SNAPV(ISNV,J)*FACDT*ZW
  247       CONTINUE
            SCLTAL(ISNV,NTALT)=2
          ELSEIF (ISNVE(ISNV).EQ.3) THEN
C  SCALE # AMP/S PER VOLUME
            DO 248 J=1,NSBOX_TAL
              SNAPV(ISNV,J)=SNAPV(ISNV,J)*FACDT*ZVOLIW(J)
  248       CONTINUE
            SCLTAL(ISNV,NTALT)=3
          ELSEIF (ISNVE(ISNV).EQ.4) THEN
C  SCALE # AMP/S PER CELL
            DO 249 J=1,NSBOX_TAL
              SNAPV(ISNV,J)=SNAPV(ISNV,J)*FACDT*ZWW
  249       CONTINUE
            SCLTAL(ISNV,NTALT)=4
          ELSE
C  DO NOT SCALE AT ALL
            DO J=1,NSBOX_TAL
              SNAPV(ISNV,J)=SNAPV(ISNV,J)*FACDT
            ENDDO
            SCLTAL(ISNV,NTALT)=0
          ENDIF
  245   CONTINUE
      END IF
C
C  TALLIES FOR FLUID CODE COUPLING FOR THE STRATUM ISTRA
C  TALLY NO. NTALM
C
      IF (LCOPV) THEN
        DO 255 ICPV=1,NCPVI
          IF (ICPVE(ICPV).EQ.1) THEN
C  SCALE # PER VOLUME
            DO 256 J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*ZVOLIN(J)
  256       CONTINUE
            SCLTAL(ICPV,NTALM)=1
          ELSEIF (ICPVE(ICPV).EQ.2) THEN
C  SCALE # PER CELL
            DO 257 J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*ZW
  257       CONTINUE
            SCLTAL(ICPV,NTALM)=2
          ELSEIF (ICPVE(ICPV).EQ.3) THEN
C  SCALE # AMP/S PER VOLUME
            DO 258 J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*ZVOLIW(J)
  258       CONTINUE
            SCLTAL(ICPV,NTALM)=3
          ELSEIF (ICPVE(ICPV).EQ.4) THEN
C  SCALE # AMP/S PER CELL
            DO 259 J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*ZWW
  259       CONTINUE
            SCLTAL(ICPV,NTALM)=4
          ELSE
C  DO NOT SCALE AT ALL
            SCLTAL(ICPV,NTALM)=0
          ENDIF
  255   CONTINUE
      END IF
C
C  TALLIES FOR ITERATIVE MODE, E.G.: BGK OPTION. STRATUM ISTRA
C  TALLY NO. NTALB
C
      IF (LBGKV) THEN
        DO 265 IBGV=1,NBGVI
          IF (IBGVE(IBGV).EQ.1) THEN
C  SCALE # PER VOLUME
            DO 266 J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*ZVOLIN(J)
  266       CONTINUE
            SCLTAL(IBGV,NTALB)=1
          ELSEIF (IBGVE(IBGV).EQ.2) THEN
C  SCALE # PER CELL
            DO 267 J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*ZW
  267       CONTINUE
            SCLTAL(IBGV,NTALB)=2
          ELSEIF (IBGVE(IBGV).EQ.3) THEN
C  SCALE # AMP/S PER VOLUME
            DO 268 J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*ZVOLIW(J)
  268       CONTINUE
            SCLTAL(IBGV,NTALB)=3
          ELSEIF (IBGVE(IBGV).EQ.4) THEN
C  SCALE # AMP/S PER CELL
            DO 269 J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*ZWW
  269       CONTINUE
            SCLTAL(IBGV,NTALB)=4
          ELSE
C  DO NOT SCALE AT ALL
            SCLTAL(IBGV,NTALB)=0
          ENDIF
  265   CONTINUE
      END IF
C
C  OTHER TALLIES ESTIMATED FROM HISTORIES, NO FIRST (SPECIES) INDEX
C
      DO 270 J=1,NSBOX_TAL
        IF (LPAEL)  PAEL(J) =PAEL(J) *ZVOLIW(J)
        IF (LEAEL)  EAEL(J) =EAEL(J) *ZVOLIW(J)
        IF (LEAAT)  EAAT(J) =EAAT(J) *ZVOLIW(J)
        IF (LEAML)  EAML(J) =EAML(J) *ZVOLIW(J)
        IF (LEAIO)  EAIO(J) =EAIO(J) *ZVOLIW(J)
        IF (LEAPHT) EAPHT(J)=EAPHT(J)*ZVOLIW(J)
C
        IF (LPMEL)  PMEL(J) =PMEL(J) *ZVOLIW(J)
        IF (LEMEL)  EMEL(J) =EMEL(J) *ZVOLIW(J)
        IF (LEMAT)  EMAT(J) =EMAT(J) *ZVOLIW(J)
        IF (LEMML)  EMML(J) =EMML(J) *ZVOLIW(J)
        IF (LEMIO)  EMIO(J) =EMIO(J) *ZVOLIW(J)
        IF (LEMPHT) EMPHT(J)=EMPHT(J)*ZVOLIW(J)
C
        IF (LPIEL)  PIEL(J) =PIEL(J) *ZVOLIW(J)
        IF (LEIEL)  EIEL(J) =EIEL(J) *ZVOLIW(J)
        IF (LEIAT)  EIAT(J) =EIAT(J) *ZVOLIW(J)
        IF (LEIML)  EIML(J) =EIML(J) *ZVOLIW(J)
        IF (LEIIO)  EIIO(J) =EIIO(J) *ZVOLIW(J)
        IF (LEIPHT) EIPHT(J)=EIPHT(J)*ZVOLIW(J)
C
        IF (LPPHEL)  PPHEL(J) =PPHEL(J) *ZVOLIW(J)
        IF (LEPHEL)  EPHEL(J) =EPHEL(J) *ZVOLIW(J)
        IF (LEPHAT)  EPHAT(J) =EPHAT(J) *ZVOLIW(J)
        IF (LEPHML)  EPHML(J) =EPHML(J) *ZVOLIW(J)
        IF (LEPHIO)  EPHIO(J) =EPHIO(J) *ZVOLIW(J)
        IF (LEPHPHT) EPHPHT(J)=EPHPHT(J)*ZVOLIW(J)

        IF (LEPAT)  EPAT(J) =EPAT(J) *ZVOLIW(J)
        IF (LEPML)  EPML(J) =EPML(J) *ZVOLIW(J)
        IF (LEPIO)  EPIO(J) =EPIO(J) *ZVOLIW(J)
        IF (LEPPHT) EPPHT(J)=EPPHT(J)*ZVOLIW(J)
  270 CONTINUE
      SCLTAL(1,9)=3
      SCLTAL(1,15)=3
      SCLTAL(1,21)=3
      SCLTAL(1,27)=3
      SCLTAL(1,33)=3
      SCLTAL(1,34)=3
      SCLTAL(1,35)=3
      SCLTAL(1,36)=3
      SCLTAL(1,37)=3
      SCLTAL(1,39)=3
      SCLTAL(1,40)=3
      SCLTAL(1,41)=3
      SCLTAL(1,42)=3
      SCLTAL(1,43)=3
      SCLTAL(1,45)=3
      SCLTAL(1,46)=3
      SCLTAL(1,47)=3
      SCLTAL(1,48)=3
      SCLTAL(1,49)=3
      SCLTAL(1,51)=3
      SCLTAL(1,52)=3
      SCLTAL(1,53)=3
      SCLTAL(1,54)=3
      SCLTAL(1,55)=3
      SCLTAL(1,80)=3
      SCLTAL(1,81)=3
      SCLTAL(1,82)=3
      SCLTAL(1,83)=3

C  CHECK: ALL TALLIES 1 -- 100 SCALED ?   TBD.


C   SCALE AND INTEGRATE VOLUMETRIC SPECTRA

      DO ISPC=1,NADSPC
        IF (ESTIML(ISPC)%ISRFCLL /= 0) THEN
          ICL = ESTIML(ISPC)%ISPCSRF
          IF (ESTIML(ISPC)%ISRFCLL == 1) THEN
!  scoring cell (coarse grid)
            ZFAC = ZVOLIN(ICL)     ! ZW/VOL or ZW/VOLTAL ?
          ELSE IF (ESTIML(ISPC)%ISRFCLL == 2) THEN
!  geometry cell (fine grid)
            ZFAC = ZW / VOL(ICL)   ! ZW/VOL or ZW/VOLTAL ?
          ELSE
            ZFAC = 1._DP
            WRITE (IUNOUT,*) ' NO SCALING PERFORMED FOR SPECTRUM NO. ',
     .                         ISPC
          END IF
          ESTIML(ISPC)%SPCS = 0._DP
          DO I = 0, ESTIML(ISPC)%NSPC+1
            IF (I.EQ.0) THEN
              ELEFT = MIN(ESTIML(ISPC)%SPCMIN,
     .                   ESTIML(ISPC)%ESP_MIN)
              ERIGHT= ESTIML(ISPC)%SPCMIN
            ELSE IF (I.EQ.ESTIML(ISPC)%NSPC+1) THEN
              ELEFT = ESTIML(ISPC)%SPCMAX
              ERIGHT= MAX(ESTIML(ISPC)%SPCMAX,
     .                   ESTIML(ISPC)%ESP_MAX)
            ELSE
              ELEFT = ESTIML(ISPC)%SPCMIN+
     .               ESTIML(ISPC)%SPCDEL*(I-1)
              ERIGHT= ESTIML(ISPC)%SPCMIN+
     .               ESTIML(ISPC)%SPCDEL*I
            END IF
            IF (ESTIML(ISPC)%LOG) THEN
              DEL = 10._DP**ERIGHT-10._DP**ELEFT+EPS60
            ELSE
              DEL = ERIGHT-ELEFT+EPS60
            END IF
cdr 1/DE, DE= energy increment for bin no. I.
            DELI = 1._DP/(DEL+EPS60)

C  SCALE: FROM SCORING TALLY UNITS # PER ENERGY BIN --> TALLY UNITS # PER EV
            ESTIML(ISPC)%SPC(I) =
     .      ESTIML(ISPC)%SPC(I)*ZFAC*DELI
C  INTEGRATE--> TALLY UNITS
cdr  Test tbd: in case of total (not directional) spectrum, i.e. for IDIREC=0, this
cdr            integral must coincide with the particle density PDEN.. or the energy density EDEN..,
cdr            in the selected cell, depending on ISPTYP=1,
cdr            or ISPTYP=2, respectively.
cdr  See text in scale_surf_tallies: this test should
cdr  also work for directional resolved spectra.
            ESTIML(ISPC)%SPCS = ESTIML(ISPC)%SPCS +
     .                          ESTIML(ISPC)%SPC(I)*DEL
          END DO
        END IF
      END DO

      RETURN

      END SUBROUTINE EIRENE_SCAL_VOLAV_TALLIES
