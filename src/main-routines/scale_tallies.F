cdr   jan.2023   spumpi total (spumpi(0,istra) corrected.
cdr              spumpi is currently unused.
cdr    dec. 15:  added species index ipls, for volumetric energy source tallies for bulk ions
cdr              eapl,empl,eipl,ephpl
cdr   24.09.14:  scaling of new sputter tallies with fatm, fmol,fion,nphot: corrected
c  spring 2014:  new sputter tallies introduced: emitted species-resolved


      SUBROUTINE EIRENE_SCALE_TALLIES (ISTRA)
C
C  RESCALE TRACKLENGTH-ESTIMATED VOLUME-AVERAGED TALLIES
C  AS WELL AS SURFACE-FLUX TALLIES
C  WITH PARTICLE BALANCE CORRECTION FACTORS FATM,FMOL,FION,FPHOT
C  TO ENFORCE PERFECT GLOBAL PARTICLE BALANCE
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      USE EIRMOD_CLOGAU
      USE EIRMOD_COUTAU
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      USE EIRMOD_CESTIM
      USE EIRMOD_CSPEZ
      USE EIRMOD_CCOUPL

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), ALLOCATABLE :: FATM(:), FMOL(:), FION(:), FPHOT(:)
      REAL(DP) :: FADD
      INTEGER :: IATM, IMOL, IION, IPLS, IPHOT, IADV, ICLV, IBGV, ICPV,
     .           JATM, J, ISPC,
#if ( defined(G95) || defined(LEGACYCOMP) )
     .           JMOL, JION, JPHOT,
#endif
     .           KK, ICO
      CHARACTER*24 HLP_FRM
      EXTERNAL :: EIRENE_GETSCL4, EIRENE_GETSCLS,
     .            EIRENE_TALLY_SANITY_CHECK,
     .            EIRENE_LEER, EIRENE_MASR4, EIRENE_MASRR1
C
      ALLOCATE(FATM(0:NATMI))
      ALLOCATE(FMOL(0:NMOLI))
      ALLOCATE(FION(0:NIONI))
      ALLOCATE(FPHOT(0:NPHOTI))
cdr  always find particle balance scaling factors, per type
      CALL EIRENE_GETSCL4 (ISTRA,FATM(0),FMOL(0),FION(0),FPHOT(0))
cdr  optional: also finer rescaling per species
      IF (NLSPCSCL) THEN
        CALL EIRENE_TALLY_SANITY_CHECK(ISTRA)
        CALL EIRENE_GETSCLS (ISTRA,NATMI,NMOLI,NIONI,NPHOTI,
     .                       FATM,FMOL,FION,FPHOT)
      ELSE
        FATM(1:NATMI) = FATM(0)
        FMOL(1:NMOLI) = FMOL(0)
        FION(1:NIONI) = FION(0)
        FPHOT(1:NPHOTI) = FPHOT(0)
      END IF
C
      IF (.NOT.NLSCL) THEN

        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO RESCALING DONE (NLSCL=FALSE)'
        CALL EIRENE_LEER(2)

      ELSEIF (NLSCL) THEN

        IF (NLSPCSCL) THEN
          FASCL(:,ISTRA) = FATM
          FMSCL(:,ISTRA) = FMOL
          FISCL(:,ISTRA) = FION
          FPHSCL(:,ISTRA)= FPHOT
        ELSE
          FASCL(:,ISTRA) = FATM(0)
          FMSCL(:,ISTRA) = FMOL(0)
          FISCL(:,ISTRA) = FION(0)
          FPHSCL(:,ISTRA)= FPHOT(0)
        END IF
        call eirene_leer(1)
        write (iunout,*) 'ISTRA = ',ISTRA
        write (hlp_frm,'(a,i2,a)') '(a,',NATM+1,'f15.8)'
        write (iunout,hlp_frm) 'FASCL = ',FASCL(:,ISTRA)
        write (hlp_frm,'(a,i2,a)') '(a,',NMOL+1,'f15.8)'
        write (iunout,hlp_frm) 'FMSCL = ',FMSCL(:,ISTRA)
        write (hlp_frm,'(a,i2,a)') '(a,',NION+1,'f15.8)'
        write (iunout,hlp_frm) 'FISCL = ',FISCL(:,ISTRA)
        write (hlp_frm,'(a,i2,a)') '(a,',NPHOT+1,'f15.8)'
        write (iunout,hlp_frm) 'FPHSCL= ',FPHSCL(:,ISTRA)
        CALL EIRENE_LEER(2)
C
C  CARRY OUT SCALING OF VOLUME AND SURFACE TALLIES, RESP.
C
C  ATOM TALLIES
C
c  volumetric atomic tallies
        DO 2101 IATM=1,NATMI
          if (.not.logatm(iatm,istra)) cycle
          DO 111 J=1,NSBOX_TAL
            IF (LPDENA) PDENA(IATM,J)=PDENA(IATM,J)*FATM(IATM)
            IF (LEDENA) EDENA(IATM,J)=EDENA(IATM,J)*FATM(IATM)
            IF (LPAAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAAT2(1:NATM,LB_ATM:NATM) => PAAT(:,J)
                PAAT2(IATM,1:NATMI) = PAAT2(IATM,1:NATMI)*FATM(1:NATMI)
                PAAT(IATM,J) = SUM(PAAT2(IATM,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PAAT(IATM+(JATM-LB_ATM)*NATM,J) =
     .             PAAT(IATM+(JATM-LB_ATM)*NATM,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PAAT(IATM,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PAAT(IATM,J) = PAAT(IATM,J) +
     .             PAAT(IATM+(JATM-LB_ATM)*NATM,J)
                END DO
#endif
              ELSE
                PAAT(IATM,J)=PAAT(IATM,J)*FATM(0)
              END IF
            END IF
            IF (LPMAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMAT2(1:NATM,LB_MOL:NMOL) => PMAT(:,J)
                PMAT2(IATM,1:NMOLI) = PMAT2(IATM,1:NMOLI)*FMOL(1:NMOLI)
                PMAT(IATM,J) = SUM(PMAT2(IATM,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PMAT(IATM+(JMOL-LB_MOL)*NATM,J) =
     .             PMAT(IATM+(JMOL-LB_MOL)*NATM,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PMAT(IATM,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PMAT(IATM,J) = PMAT(IATM,J) +
     .             PMAT(IATM+(JMOL-LB_MOL)*NATM,J)
                END DO
#endif
              ELSE
                PMAT(IATM,J)=PMAT(IATM,J)*FMOL(0)
              END IF
            END IF
            IF (LPIAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIAT2(1:NATM,LB_ION:NION) => PIAT(:,J)
                PIAT2(IATM,1:NIONI) = PIAT2(IATM,1:NIONI)*FION(1:NIONI)
                PIAT(IATM,J) = SUM(PIAT2(IATM,1:NIONI))
#else
                DO JION = 1, NIONI
                  PIAT(IATM+(JION-LB_ION)*NATM,J) =
     .             PIAT(IATM+(JION-LB_ION)*NATM,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PIAT(IATM,J) = 0.0_DP
                DO JION = 1, NIONI
                  PIAT(IATM,J) = PIAT(IATM,J) +
     .             PIAT(IATM+(JION-LB_ION)*NATM,J)
                END DO
#endif
              ELSE
                PIAT(IATM,J)=PIAT(IATM,J)*FION(0)
              END IF
            END IF
            IF (LPPHAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHAT2(1:NATM,LB_PHOT:NPHOT) => PPHAT(:,J)
                PPHAT2(IATM,1:NPHOTI) = PPHAT2(IATM,1:NPHOTI)*
     .                                  FPHOT(1:NPHOTI)
                PPHAT(IATM,J) = SUM(PPHAT2(IATM,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J) =
     .             PPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PPHAT(IATM,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PPHAT(IATM,J) = PPHAT(IATM,J) +
     .             PPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J)
                END DO
#endif
              ELSE
                PPHAT(IATM,J)=PPHAT(IATM,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENA) PGENA(IATM,J)=PGENA(IATM,J)*FATM(IATM)
            IF (LEGENA) EGENA(IATM,J)=EGENA(IATM,J)*FATM(IATM)
            IF (LVGENA) VGENA(IATM,J)=VGENA(IATM,J)*FATM(IATM)
            IF (LVXDENA) VXDENA(IATM,J)=VXDENA(IATM,J)*FATM(IATM)
            IF (LVYDENA) VYDENA(IATM,J)=VYDENA(IATM,J)*FATM(IATM)
            IF (LVZDENA) VZDENA(IATM,J)=VZDENA(IATM,J)*FATM(IATM)
            IF (LRAEL) RAEL(IATM,J)=RAEL(IATM,J)*FATM(IATM)
  111     CONTINUE
c  atomic surface tallies, outgoing species: iatm
          DO 310 J=1,NLMPGS
            IF (LPOTAT) POTAT(IATM,J)=POTAT(IATM,J)*FATM(IATM)
            IF (LPRFAAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAAT2(1:NATM,LB_ATM:NATM) => PRFAAT(:,J)
                PRFAAT2(IATM,1:NATMI) = PRFAAT2(IATM,1:NATMI)*
     .                                  FATM(1:NATMI)
                PRFAAT(IATM,J) = SUM(PRFAAT2(IATM,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PRFAAT(IATM+(JATM-LB_ATM)*NATM,J) =
     .             PRFAAT(IATM+(JATM-LB_ATM)*NATM,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PRFAAT(IATM,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PRFAAT(IATM,J) = PRFAAT(IATM,J) +
     .             PRFAAT(IATM+(JATM-LB_ATM)*NATM,J)
                END DO
#endif
              ELSE
                PRFAAT(IATM,J)=PRFAAT(IATM,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMAT2(1:NATM,LB_MOL:NMOL) => PRFMAT(:,J)
                PRFMAT2(IATM,1:NMOLI) = PRFMAT2(IATM,1:NMOLI)*
     .                                  FMOL(1:NMOLI)
                PRFMAT(IATM,J) = SUM(PRFMAT2(IATM,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PRFMAT(IATM+(JMOL-LB_MOL)*NATM,J) =
     .             PRFMAT(IATM+(JMOL-LB_MOL)*NATM,J)*FMOL(JMOL)
                END DO
               IF (LB_MOL.EQ.0) PRFMAT(IATM,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PRFMAT(IATM,J) = PRFMAT(IATM,J) +
     .             PRFMAT(IATM+(JMOL-LB_MOL)*NATM,J)
                END DO
#endif
              ELSE
                PRFMAT(IATM,J)=PRFMAT(IATM,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIAT2(1:NATM,LB_ION:NION) => PRFIAT(:,J)
                PRFIAT2(IATM,1:NIONI) = PRFIAT2(IATM,1:NIONI)*
     .                                  FION(1:NIONI)
                PRFIAT(IATM,J) = SUM(PRFIAT2(IATM,1:NIONI))
#else
                DO JION = 1, NIONI
                  PRFIAT(IATM+(JION-LB_ION)*NATM,J) =
     .             PRFIAT(IATM+(JION-LB_ION)*NATM,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PRFIAT(IATM,J) = 0.0_DP
                DO JION = 1, NIONI
                  PRFIAT(IATM,J) = PRFIAT(IATM,J) +
     .             PRFIAT(IATM+(JION-LB_ION)*NATM,J)
                END DO
#endif
              ELSE
                PRFIAT(IATM,J)=PRFIAT(IATM,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHAT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHAT2(1:NATM,LB_PHOT:NPHOT) => PRFPHAT(:,J)
                PRFPHAT2(IATM,1:NPHOTI) = PRFPHAT2(IATM,1:NPHOTI)*
     .                                    FPHOT(1:NPHOTI)
                PRFPHAT(IATM,J) = SUM(PRFPHAT2(IATM,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PRFPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J) =
     .             PRFPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PRFPHAT(IATM,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PRFPHAT(IATM,J) = PRFPHAT(IATM,J) +
     .             PRFPHAT(IATM+(JPHOT-LB_PHOT)*NATM,J)
                END DO
#endif
              ELSE
                PRFPHAT(IATM,J)=PRFPHAT(IATM,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTAT) EOTAT(IATM,J)=EOTAT(IATM,J)*FATM(IATM)
            IF (LERFAAT) ERFAAT(IATM,J)=ERFAAT(IATM,J)*FATM(IATM)
            IF (LERFMAT) ERFMAT(IATM,J)=ERFMAT(IATM,J)*FMOL(0)
            IF (LERFIAT) ERFIAT(IATM,J)=ERFIAT(IATM,J)*FION(0)
            IF (LERFPHAT) ERFPHAT(IATM,J)=ERFPHAT(IATM,J)*FPHOT(0)
cdr sputter tallies
            IF (LSPTAAT) SPTAAT(IATM,J)=SPTAAT(IATM,J)*FATM(IATM)
            IF (LSPTMAT) SPTMAT(IATM,J)=SPTMAT(IATM,J)*FMOL(0)
            IF (LSPTIAT) SPTIAT(IATM,J)=SPTIAT(IATM,J)*FION(0)
            IF (LSPTPHAT) SPTPHAT(IATM,J)=SPTPHAT(IATM,J)*FPHOT(0)
cdr  ?? scaling with bulk flux?
            IF (LSPTPAT) SPTPAT(IATM,J)=SPTPAT(IATM,J)*FATM(IATM)
            IF (LSPUMP) SPUMP(NSPH+IATM,J)=SPUMP(NSPH+IATM,J)*FATM(IATM)
  310     CONTINUE

 2101   CONTINUE

        if (lsptatot) sptatot = sptatot * fatm(0)

c  integrated atomic tallies, both volumetric and surface-averaged
        DO 2111 IATM=1,NATMI
          if (.not.logatm(iatm,istra)) cycle
          PDENAI(IATM,ISTRA)=PDENAI(IATM,ISTRA)*FATM(IATM)
          EDENAI(IATM,ISTRA)=EDENAI(IATM,ISTRA)*FATM(IATM)
          PGENAI(IATM,ISTRA)=PGENAI(IATM,ISTRA)*FATM(IATM)
          EGENAI(IATM,ISTRA)=EGENAI(IATM,ISTRA)*FATM(IATM)
          VGENAI(IATM,ISTRA)=VGENAI(IATM,ISTRA)*FATM(IATM)
          VXDENAI(IATM,ISTRA)=VXDENAI(IATM,ISTRA)*FATM(IATM)
          VYDENAI(IATM,ISTRA)=VYDENAI(IATM,ISTRA)*FATM(IATM)
          VZDENAI(IATM,ISTRA)=VZDENAI(IATM,ISTRA)*FATM(IATM)

          POTATI(IATM,ISTRA)=POTATI(IATM,ISTRA)*FATM(IATM)
          EOTATI(IATM,ISTRA)=EOTATI(IATM,ISTRA)*FATM(IATM)
          ERFAAI(IATM,ISTRA)=ERFAAI(IATM,ISTRA)*FATM(IATM)
          ERFMAI(IATM,ISTRA)=ERFMAI(IATM,ISTRA)*FMOL(0)
          ERFIAI(IATM,ISTRA)=ERFIAI(IATM,ISTRA)*FION(0)
          ERFPHAI(IATM,ISTRA)=ERFPHAI(IATM,ISTRA)*FPHOT(0)

          SPTAATI(IATM,ISTRA)=SPTAATI(IATM,ISTRA)*FATM(IATM)
          SPTMATI(IATM,ISTRA)=SPTMATI(IATM,ISTRA)*FMOL(0)
          SPTIATI(IATM,ISTRA)=SPTIATI(IATM,ISTRA)*FION(0)
          SPTPHATI(IATM,ISTRA)=SPTPHATI(IATM,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPATI(IATM,ISTRA)=SPTPATI(IATM,ISTRA)*FATM(IATM)
          SPUMPI(NSPH+IATM,ISTRA)=SPUMPI(NSPH+IATM,ISTRA)*FATM(IATM)
          RAELI(IATM,ISTRA)=RAELI(IATM,ISTRA)*FATM(IATM)

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAATI2(0:NATM,LB_ATM:NATM) => PAATI(:,ISTRA)
            PAATI2(IATM,1:NATMI) = PAATI2(IATM,1:NATMI)*FATM(1:NATMI)
            PAATI(IATM,ISTRA) = SUM(PAATI2(IATM,1:NATMI))
            PRFAAI2(0:NATM,LB_ATM:NATM) => PRFAAI(:,ISTRA)
            PRFAAI2(IATM,1:NATMI) = PRFAAI2(IATM,1:NATMI)*FATM(1:NATMI)
            PRFAAI(IATM,ISTRA) = SUM(PRFAAI2(IATM,1:NATMI))
#else
            DO JATM = 1, NATMI
              PAATI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA) =
     .         PAATI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA)*FATM(JATM)
              PRFAAI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA) =
     .         PRFAAI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA)*FATM(JATM)
            END DO
            IF (LB_ATM.EQ.0) PAATI(IATM,ISTRA) = 0.0_DP
            IF (LB_ATM.EQ.0) PRFAAI(IATM,ISTRA) = 0.0_DP
            DO JATM = 1, NATMI
              PAATI(IATM,ISTRA) = PAATI(IATM,ISTRA) +
     .         PAATI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA)
              PRFAAI(IATM,ISTRA) = PRFAAI(IATM,ISTRA) +
     .         PRFAAI(IATM+(JATM-LB_ATM)*(NATM+1),ISTRA)
            END DO
#endif
          ELSE
            PAATI(IATM,ISTRA)=PAATI(IATM,ISTRA)*FATM(0)
            PRFAAI(IATM,ISTRA)=PRFAAI(IATM,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMATI2(0:NATM,LB_MOL:NMOL) => PMATI(:,ISTRA)
            PMATI2(IATM,1:NMOLI) = PMATI2(IATM,1:NMOLI)*FMOL(1:NMOLI)
            PMATI(IATM,ISTRA) = SUM(PMATI2(IATM,1:NMOLI))
            PRFMAI2(0:NATM,LB_MOL:NMOL) => PRFMAI(:,ISTRA)
            PRFMAI2(IATM,1:NMOLI) = PRFMAI2(IATM,1:NMOLI)*FMOL(1:NMOLI)
            PRFMAI(IATM,ISTRA) = SUM(PRFMAI2(IATM,1:NMOLI))
#else
            DO JMOL = 1, NMOLI
              PMATI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA) =
     .         PMATI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA)*FMOL(JMOL)
              PRFMAI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA) =
     .         PRFMAI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA)*FMOL(JMOL)
            END DO
            IF (LB_MOL.EQ.0) PMATI(IATM,ISTRA) = 0.0_DP
            IF (LB_MOL.EQ.0) PRFMAI(IATM,ISTRA) = 0.0_DP
            DO JMOL = 1, NMOLI
              PMATI(IATM,ISTRA) = PMATI(IATM,ISTRA) +
     .         PMATI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA)
              PRFMAI(IATM,ISTRA) = PRFMAI(IATM,ISTRA) +
     .         PRFMAI(IATM+(JMOL-LB_MOL)*(NATM+1),ISTRA)
            END DO
#endif
          ELSE
            PMATI(IATM,ISTRA)=PMATI(IATM,ISTRA)*FMOL(0)
            PRFMAI(IATM,ISTRA)=PRFMAI(IATM,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIATI2(0:NATM,LB_ION:NION) => PIATI(:,ISTRA)
            PIATI2(IATM,1:NIONI) = PIATI2(IATM,1:NIONI)*FION(1:NIONI)
            PIATI(IATM,ISTRA) = SUM(PIATI2(IATM,1:NIONI))
            PRFIAI2(0:NATM,LB_ION:NION) => PRFIAI(:,ISTRA)
            PRFIAI2(IATM,1:NIONI) = PRFIAI2(IATM,1:NIONI)*FION(1:NIONI)
            PRFIAI(IATM,ISTRA) = SUM(PRFIAI2(IATM,1:NIONI))
#else
            DO JION = 1, NIONI
              PIATI(IATM+(JION-LB_ION)*(NATM+1),ISTRA) =
     .         PIATI(IATM+(JION-LB_ION)*(NATM+1),ISTRA)*FION(JION)
              PRFIAI(IATM+(JION-LB_ION)*(NATM+1),ISTRA) =
     .         PRFIAI(IATM+(JION-LB_ION)*(NATM+1),ISTRA)*FION(JION)
            END DO
            IF (LB_ION.EQ.0) PIATI(IATM,ISTRA) = 0.0_DP
            IF (LB_ION.EQ.0) PRFIAI(IATM,ISTRA) = 0.0_DP
            DO JION = 1, NIONI
              PIATI(IATM,ISTRA) = PIATI(IATM,ISTRA) +
     .         PIATI(IATM+(JION-LB_ION)*(NATM+1),ISTRA)
              PRFIAI(IATM,ISTRA) = PRFIAI(IATM,ISTRA) +
     .         PRFIAI(IATM+(JION-LB_ION)*(NATM+1),ISTRA)
            END DO
#endif
          ELSE
            PIATI(IATM,ISTRA)=PIATI(IATM,ISTRA)*FION(0)
            PRFIAI(IATM,ISTRA)=PRFIAI(IATM,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHATI2(0:NATM,LB_PHOT:NPHOT) => PPHATI(:,ISTRA)
            PPHATI2(IATM,1:NPHOT) = PPHATI2(IATM,1:NPHOT)*
     .                              FPHOT(1:NPHOTI)
            PPHATI(IATM,ISTRA) = SUM(PPHATI2(IATM,1:NPHOT))
            PRFPHAI2(0:NATM,0:NPHOT) => PRFPHAI(:,ISTRA)
            PRFPHAI2(IATM,1:NPHOTI) = PRFPHAI2(IATM,1:NPHOTI)*
     .                                FPHOT(1:NPHOTI)
            PRFPHAI(IATM,ISTRA) = SUM(PRFPHAI2(IATM,1:NPHOTI))
#else
            DO JPHOT = 1, NPHOTI
              PPHATI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA) =
     .         PPHATI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA)*FPHOT(JPHOT)
              PRFPHAI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA) =
     .         PRFPHAI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA)*FPHOT(JPHOT)
            END DO
            IF (LB_PHOT.EQ.0) PPHATI(IATM,ISTRA) = 0.0_DP
            IF (LB_PHOT.EQ.0) PRFPHAI(IATM,ISTRA) = 0.0_DP
            DO JPHOT = 1, NPHOTI
              PPHATI(IATM,ISTRA) = PPHATI(IATM,ISTRA) +
     .         PPHATI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA)
              PRFPHAI(IATM,ISTRA) = PRFPHAI(IATM,ISTRA) +
     .         PRFPHAI(IATM+(JPHOT-LB_PHOT)*(NATM+1),ISTRA)
            END DO
#endif
          ELSE
            PPHATI(IATM,ISTRA)=PPHATI(IATM,ISTRA)*FPHOT(0)
            PRFPHAI(IATM,ISTRA)=PRFPHAI(IATM,ISTRA)*FPHOT(0)
          END IF

 2111   CONTINUE

        IF (LOGATM(0,ISTRA)) THEN
          PDENAI(0,ISTRA)=SUM(PDENAI(1:NATMI,ISTRA))
          EDENAI(0,ISTRA)=SUM(EDENAI(1:NATMI,ISTRA))
          PGENAI(0,ISTRA)=SUM(PGENAI(1:NATMI,ISTRA))
          EGENAI(0,ISTRA)=SUM(EGENAI(1:NATMI,ISTRA))
          VGENAI(0,ISTRA)=SUM(VGENAI(1:NATMI,ISTRA))
          VXDENAI(0,ISTRA)=SUM(VXDENAI(1:NATMI,ISTRA))
          VYDENAI(0,ISTRA)=SUM(VYDENAI(1:NATMI,ISTRA))
          VZDENAI(0,ISTRA)=SUM(VZDENAI(1:NATMI,ISTRA))

          POTATI(0,ISTRA)=SUM(POTATI(1:NATMI,ISTRA))
          EOTATI(0,ISTRA)=SUM(EOTATI(1:NATMI,ISTRA))
          ERFAAI(0,ISTRA)=SUM(ERFAAI(1:NATMI,ISTRA))
          ERFMAI(0,ISTRA)=SUM(ERFMAI(1:NATMI,ISTRA))
          ERFIAI(0,ISTRA)=SUM(ERFIAI(1:NATMI,ISTRA))
          ERFPHAI(0,ISTRA)=SUM(ERFPHAI(1:NATMI,ISTRA))

          SPTAATI(0,ISTRA)=SUM(SPTAATI(1:NATMI,ISTRA))
          SPTMATI(0,ISTRA)=SUM(SPTMATI(1:NATMI,ISTRA))
          SPTIATI(0,ISTRA)=SUM(SPTIATI(1:NATMI,ISTRA))
          SPTPHATI(0,ISTRA)=SUM(SPTPHATI(1:NATMI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPATI(0,ISTRA)=SUM(SPTPATI(1:NATMI,ISTRA))
          RAELI(0,ISTRA)=SUM(RAELI(1:NATMI,ISTRA))

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAATI2(0:NATM,LB_ATM:NATM) => PAATI(:,ISTRA)
            PRFAAI2(0:NATM,LB_ATM:NATM) => PRFAAI(:,ISTRA)
            DO JATM = LB_ATM,NATMI
              PAATI2(0,JATM) = SUM(PAATI2(1:NATMI,JATM))
              PRFAAI2(0,JATM) = SUM(PRFAAI2(1:NATMI,JATM))
            END DO
#else
            DO JATM = LB_ATM,NATMI
              PAATI(JATM*(NATM+1),ISTRA) =
     .         SUM(PAATI(1+JATM*(NATM+1):NATMI+JATM*(NATM+1),ISTRA))
              PRFAAI(JATM*(NATM+1),ISTRA) =
     .         SUM(PRFAAI(1+JATM*(NATM+1):NATMI+JATM*(NATM+1),ISTRA))
            END DO
#endif
          ELSE
            PAATI(0,ISTRA)=SUM(PAATI(1:NATMI,ISTRA))
            PRFAAI(0,ISTRA)=SUM(PRFAAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMATI2(0:NATM,LB_MOL:NMOL) => PMATI(:,ISTRA)
            PRFMAI2(0:NATM,LB_MOL:NMOL) => PRFMAI(:,ISTRA)
            DO IMOL=LB_MOL,NMOLI
              PMATI2(0,IMOL) = SUM(PMATI2(1:NATMI,IMOL))
              PRFMAI2(0,IMOL) = SUM(PRFMAI2(1:NATMI,IMOL))
            END DO
#else
            DO IMOL = LB_MOL, NMOLI
              PMATI(IMOL*(NATM+1),ISTRA) =
     .         SUM(PMATI(1+IMOL*(NATM+1):NATMI+IMOL*(NATM+1),ISTRA))
              PRFMAI(IMOL*(NATM+1),ISTRA) =
     .         SUM(PRFMAI(1+IMOL*(NATM+1):NATMI+IMOL*(NATM+1),ISTRA))
            END DO
#endif
          ELSE
            PMATI(0,ISTRA)=SUM(PMATI(1:NATMI,ISTRA))
            PRFMAI(0,ISTRA)=SUM(PRFMAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIATI2(0:NATM,LB_ION:NION) => PIATI(:,ISTRA)
            PRFIAI2(0:NATM,LB_ION:NION) => PRFIAI(:,ISTRA)
            DO IION=LB_ION,NIONI
              PIATI2(0,IION) = SUM(PIATI2(1:NATMI,IION))
              PRFIAI2(0,IION) = SUM(PRFIAI2(1:NATMI,IION))
            END DO
#else
            DO IION = LB_ION, NIONI
              PIATI(IION*(NATM+1),ISTRA) =
     .         SUM(PIATI(1+IION*(NATM+1):NATMI+IION*(NATM+1),ISTRA))
              PRFIAI(IION*(NATM+1),ISTRA) =
     .         SUM(PRFIAI(1+IION*(NATM+1):NATMI+IION*(NATM+1),ISTRA))
            END DO
#endif
          ELSE
            PIATI(0,ISTRA)=SUM(PIATI(1:NATMI,ISTRA))
            PRFIAI(0,ISTRA)=SUM(PRFIAI(1:NATMI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHATI2(0:NATM,LB_PHOT:NPHOT) => PPHATI(:,ISTRA)
            PRFPHAI2(0:NATM,LB_PHOT:NPHOT) => PRFPHAI(:,ISTRA)
            DO IPHOT=LB_PHOT,NPHOTI
              PPHATI2(0,IPHOT) = SUM(PPHATI2(1:NATMI,IPHOT))
              PRFPHAI2(0,IPHOT) = SUM(PRFPHAI2(1:NATMI,IPHOT))
            END DO
#else
            DO IPHOT = LB_PHOT, NPHOTI
              PPHATI(IPHOT*(NATM+1),ISTRA) =
     .         SUM(PPHATI(1+IPHOT*(NATM+1):NATMI+IPHOT*(NATM+1),ISTRA))
              PRFPHAI(IPHOT*(NATM+1),ISTRA) =
     .         SUM(PRFPHAI(1+IPHOT*(NATM+1):NATMI+IPHOT*(NATM+1),ISTRA))
            END DO
#endif
          ELSE
            PPHATI(0,ISTRA)=SUM(PPHATI(1:NATMI,ISTRA))
            PRFPHAI(0,ISTRA)=SUM(PRFPHAI(1:NATMI,ISTRA))
          END IF

        END IF

        sptatti(istra) = sptatti(istra)*fatm(0)

        DO 2112 J=1,NSBOX_TAL
          IF (LEAAT) EAAT(J)=EAAT(J)*FATM(0)
          IF (LEMAT) EMAT(J)=EMAT(J)*FMOL(0)
          IF (LEIAT) EIAT(J)=EIAT(J)*FION(0)
          IF (LEPHAT) EPHAT(J)=EPHAT(J)*FPHOT(0)
 2112   CONTINUE
        EAATI(ISTRA)=EAATI(ISTRA)*FATM(0)
        EMATI(ISTRA)=EMATI(ISTRA)*FMOL(0)
        EIATI(ISTRA)=EIATI(ISTRA)*FION(0)
        EPHATI(ISTRA)=EPHATI(ISTRA)*FPHOT(0)
C
C  MOLECULE TALLIES
C
        DO 2115 IMOL=1,NMOLI
          if (.not.logmol(imol,istra)) cycle
          DO 115 J=1,NSBOX_TAL
            IF (LPDENM) PDENM(IMOL,J)=PDENM(IMOL,J)*FMOL(IMOL)
            IF (LEDENM) EDENM(IMOL,J)=EDENM(IMOL,J)*FMOL(IMOL)
            IF (LPAML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAML2(1:NMOL,LB_ATM:NATM) => PAML(:,J)
                PAML2(IMOL,1:NATMI) = PAML2(IMOL,1:NATMI)*FATM(1:NATMI)
                PAML(IMOL,J) = SUM(PAML2(IMOL,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PAML(IMOL+(JATM-LB_ATM)*NMOL,J) =
     .             PAML(IMOL+(JATM-LB_ATM)*NMOL,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PAML(IMOL,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PAML(IMOL,J) = PAML(IMOL,J) +
     .             PAML(IMOL+(JATM-LB_ATM)*NMOL,J)
                END DO
#endif
              ELSE
                PAML(IMOL,J)=PAML(IMOL,J)*FATM(0)
              END IF
            END IF
            IF (LPMML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMML2(1:NMOL,LB_MOL:NMOL) => PMML(:,J)
                PMML2(IMOL,1:NMOLI) = PMML2(IMOL,1:NMOLI)*FMOL(1:NMOLI)
                PMML(IMOL,J) = SUM(PMML2(IMOL,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PMML(IMOL+(JMOL-LB_MOL)*NMOL,J) =
     .             PMML(IMOL+(JMOL-LB_MOL)*NMOL,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PMML(IMOL,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PMML(IMOL,J) = PMML(IMOL,J) +
     .             PMML(IMOL+(JMOL-LB_MOL)*NMOL,J)
                END DO
#endif
              ELSE
                PMML(IMOL,J)=PMML(IMOL,J)*FMOL(0)
              END IF
            END IF
            IF (LPIML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIML2(1:NMOL,LB_ION:NION) => PIML(:,J)
                PIML2(IMOL,1:NIONI) = PIML2(IMOL,1:NIONI)*FION(1:NIONI)
                PIML(IMOL,J) = SUM(PIML2(IMOL,1:NIONI))
#else
                DO JION = 1, NIONI
                  PIML(IMOL+(JION-LB_ION)*NMOL,J) =
     .             PIML(IMOL+(JION-LB_ION)*NMOL,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PIML(IMOL,J) = 0.0_DP
                DO JION = 1, NIONI
                  PIML(IMOL,J) = PIML(IMOL,J) +
     .             PIML(IMOL+(JION-LB_ION)*NMOL,J)
                END DO
#endif
              ELSE
                PIML(IMOL,J)=PIML(IMOL,J)*FION(0)
              END IF
            END IF
            IF (LPPHML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHML2(1:NMOL,LB_PHOT:NPHOT) => PPHML(:,J)
                PPHML2(IMOL,1:NPHOTI) = PPHML2(IMOL,1:NPHOTI)*
     .                                  FPHOT(1:NPHOTI)
                PPHML(IMOL,J) = SUM(PPHML2(IMOL,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J) =
     .             PPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PPHML(IMOL,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PPHML(IMOL,J) = PPHML(IMOL,J) +
     .             PPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J)
                END DO
#endif
              ELSE
                PPHML(IMOL,J)=PPHML(IMOL,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENM) PGENM(IMOL,J)=PGENM(IMOL,J)*FMOL(IMOL)
            IF (LEGENM) EGENM(IMOL,J)=EGENM(IMOL,J)*FMOL(IMOL)
            IF (LVGENM) VGENM(IMOL,J)=VGENM(IMOL,J)*FMOL(IMOL)
            IF (LVXDENM) VXDENM(IMOL,J)=VXDENM(IMOL,J)*FMOL(IMOL)
            IF (LVYDENM) VYDENM(IMOL,J)=VYDENM(IMOL,J)*FMOL(IMOL)
            IF (LVZDENM) VZDENM(IMOL,J)=VZDENM(IMOL,J)*FMOL(IMOL)
            IF (LRMEL) RMEL(IMOL,J)=RMEL(IMOL,J)*FMOL(IMOL)
  115     CONTINUE
          DO 315 J=1,NLMPGS
            IF (LPOTML) POTML(IMOL,J)=POTML(IMOL,J)*FMOL(IMOL)
            IF (LPRFAML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAML2(1:NMOL,LB_ATM:NATM) => PRFAML(:,J)
                PRFAML2(IMOL,1:NATMI) = PRFAML2(IMOL,1:NATMI)*
     .                                  FATM(1:NATMI)
                PRFAML(IMOL,J) = SUM(PRFAML2(IMOL,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PRFAML(IMOL+(JATM-LB_ATM)*NMOL,J) =
     .             PRFAML(IMOL+(JATM-LB_ATM)*NMOL,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PRFAML(IMOL,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PRFAML(IMOL,J) = PRFAML(IMOL,J) +
     .             PRFAML(IMOL+(JATM-LB_ATM)*NMOL,J)
                END DO
#endif
              ELSE
                PRFAML(IMOL,J)=PRFAML(IMOL,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMML2(1:NMOL,LB_MOL:NMOL) => PRFMML(:,J)
                PRFMML2(IMOL,1:NMOLI) = PRFMML2(IMOL,1:NMOLI)*
     .                                  FMOL(1:NMOLI)
                PRFMML(IMOL,J) = SUM(PRFMML2(IMOL,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PRFMML(IMOL+(JMOL-LB_MOL)*NMOL,J) =
     .             PRFMML(IMOL+(JMOL-LB_MOL)*NMOL,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PRFMML(IMOL,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PRFMML(IMOL,J) = PRFMML(IMOL,J) +
     .             PRFMML(IMOL+(JMOL-LB_MOL)*NMOL,J)
                END DO
#endif
              ELSE
                PRFMML(IMOL,J)=PRFMML(IMOL,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIML2(1:NMOL,LB_ION:NION) => PRFIML(:,J)
                PRFIML2(IMOL,1:NIONI) = PRFIML2(IMOL,1:NIONI)*
     .                                  FION(1:NIONI)
                PRFIML(IMOL,J) = SUM(PRFIML2(IMOL,1:NIONI))
#else
                DO JION = 1, NIONI
                  PRFIML(IMOL+(JION-LB_ION)*NMOL,J) =
     .             PRFIML(IMOL+(JION-LB_ION)*NMOL,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PRFIML(IMOL,J) = 0.0_DP
                DO JION = 1, NIONI
                  PRFIML(IMOL,J) = PRFIML(IMOL,J) +
     .             PRFIML(IMOL+(JION-LB_ION)*NMOL,J)
                END DO
#endif
              ELSE
                PRFIML(IMOL,J)=PRFIML(IMOL,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHML) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHML2(1:NMOL,LB_PHOT:NPHOT) => PRFPHML(:,J)
                PRFPHML2(IMOL,1:NPHOTI) = PRFPHML2(IMOL,1:NPHOTI)*
     .                                    FPHOT(1:NPHOTI)
                PRFPHML(IMOL,J) = SUM(PRFPHML2(IMOL,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PRFPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J) =
     .             PRFPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PRFPHML(IMOL,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PRFPHML(IMOL,J) = PRFPHML(IMOL,J) +
     .             PRFPHML(IMOL+(JPHOT-LB_PHOT)*NMOL,J)
                END DO
#endif
              ELSE
                PRFPHML(IMOL,J)=PRFPHML(IMOL,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTML) EOTML(IMOL,J)=EOTML(IMOL,J)*FMOL(IMOL)
            IF (LERFAML) ERFAML(IMOL,J)=ERFAML(IMOL,J)*FATM(0)
            IF (LERFMML) ERFMML(IMOL,J)=ERFMML(IMOL,J)*FMOL(IMOL)
            IF (LERFIML) ERFIML(IMOL,J)=ERFIML(IMOL,J)*FION(0)
            IF (LERFPHML) ERFPHML(IMOL,J)=ERFPHML(IMOL,J)*FPHOT(0)

            IF (LSPTAML) SPTAML(IMOL,J)=SPTAML(IMOL,J)*FATM(0)
            IF (LSPTMML) SPTMML(IMOL,J)=SPTMML(IMOL,J)*FMOL(IMOL)
            IF (LSPTIML) SPTIML(IMOL,J)=SPTIML(IMOL,J)*FION(0)
            IF (LSPTPHML) SPTPHML(IMOL,J)=SPTPHML(IMOL,J)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPML) SPTPML(IMOL,J)=SPTPML(IMOL,J)*FMOL(IMOL)
            IF (LSPUMP) SPUMP(NSPA+IMOL,J)=SPUMP(NSPA+IMOL,J)*FMOL(IMOL)
  315     CONTINUE
 2115   CONTINUE

        if (lsptmtot) sptmtot = sptmtot*fmol(0)

        DO 2116 IMOL=1,NMOLI
          if (.not.logmol(imol,istra)) cycle
          PDENMI(IMOL,ISTRA)=PDENMI(IMOL,ISTRA)*FMOL(IMOL)
          EDENMI(IMOL,ISTRA)=EDENMI(IMOL,ISTRA)*FMOL(IMOL)
          PGENMI(IMOL,ISTRA)=PGENMI(IMOL,ISTRA)*FMOL(IMOL)
          EGENMI(IMOL,ISTRA)=EGENMI(IMOL,ISTRA)*FMOL(IMOL)
          VGENMI(IMOL,ISTRA)=VGENMI(IMOL,ISTRA)*FMOL(IMOL)
          VXDENMI(IMOL,ISTRA)=VXDENMI(IMOL,ISTRA)*FMOL(IMOL)
          VYDENMI(IMOL,ISTRA)=VYDENMI(IMOL,ISTRA)*FMOL(IMOL)
          VZDENMI(IMOL,ISTRA)=VZDENMI(IMOL,ISTRA)*FMOL(IMOL)

          POTMLI(IMOL,ISTRA)=POTMLI(IMOL,ISTRA)*FMOL(IMOL)
          EOTMLI(IMOL,ISTRA)=EOTMLI(IMOL,ISTRA)*FMOL(IMOL)
          ERFAMI(IMOL,ISTRA)=ERFAMI(IMOL,ISTRA)*FATM(0)
          ERFMMI(IMOL,ISTRA)=ERFMMI(IMOL,ISTRA)*FMOL(IMOL)
          ERFIMI(IMOL,ISTRA)=ERFIMI(IMOL,ISTRA)*FION(0)
          ERFPHMI(IMOL,ISTRA)=ERFPHMI(IMOL,ISTRA)*FPHOT(0)

          SPTAMLI(IMOL,ISTRA)=SPTAMLI(IMOL,ISTRA)*FATM(0)
          SPTMMLI(IMOL,ISTRA)=SPTMMLI(IMOL,ISTRA)*FMOL(IMOL)
          SPTIMLI(IMOL,ISTRA)=SPTIMLI(IMOL,ISTRA)*FION(0)
          SPTPHMLI(IMOL,ISTRA)=SPTPHMLI(IMOL,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPMLI(IMOL,ISTRA)=SPTPMLI(IMOL,ISTRA)*FMOL(IMOL)
          SPUMPI(NSPA+IMOL,ISTRA)=SPUMPI(NSPA+IMOL,ISTRA)*FMOL(IMOL)
          RMELI(IMOL,ISTRA)=RMELI(IMOL,ISTRA)*FMOL(IMOL)

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAMLI2(0:NMOL,LB_ATM:NATM) => PAMLI(:,ISTRA)
            PAMLI2(IMOL,1:NATMI) = PAMLI2(IMOL,1:NATMI)*FATM(1:NATMI)
            PAMLI(IMOL,ISTRA) = SUM(PAMLI2(IMOL,1:NATMI))
            PRFAMI2(0:NMOL,LB_ATM:NATM) => PRFAMI(:,ISTRA)
            PRFAMI2(IMOL,1:NATMI) = PRFAMI2(IMOL,1:NATMI)*FATM(1:NATMI)
            PRFAMI(IMOL,ISTRA) = SUM(PRFAMI2(IMOL,1:NATMI))
#else
            DO JATM = 1, NATMI
              PAMLI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA) =
     .         PAMLI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA)*FATM(JATM)
              PRFAMI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA) =
     .         PRFAMI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA)*FATM(JATM)
            END DO
            IF (LB_ATM.EQ.0) PAMLI(IMOL,ISTRA) = 0.0_DP
            IF (LB_ATM.EQ.0) PRFAMI(IMOL,ISTRA) = 0.0_DP
            DO JATM = 1, NATMI
              PAMLI(IMOL,ISTRA) = PAMLI(IMOL,ISTRA) +
     .         PAMLI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA)
              PRFAMI(IMOL,ISTRA) = PRFAMI(IMOL,ISTRA) +
     .         PRFAMI(IMOL+(JATM-LB_ATM)*(NMOL+1),ISTRA)
            END DO
#endif
          ELSE
            PAMLI(IMOL,ISTRA)=PAMLI(IMOL,ISTRA)*FATM(0)
            PRFAMI(IMOL,ISTRA)=PRFAMI(IMOL,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMMLI2(0:NMOL,LB_MOL:NMOL) => PMMLI(:,ISTRA)
            PMMLI2(IMOL,1:NMOLI) = PMMLI2(IMOL,1:NMOLI)*FMOL(1:NMOLI)
            PMMLI(IMOL,ISTRA) = SUM(PMMLI2(IMOL,1:NMOLI))
            PRFMMI2(0:NMOL,LB_MOL:NMOL) => PRFMMI(:,ISTRA)
            PRFMMI2(IMOL,1:NMOLI) = PRFMMI2(IMOL,1:NMOLI)*FMOL(1:NMOLI)
            PRFMMI(IMOL,ISTRA) = SUM(PRFMMI2(IMOL,1:NMOLI))
#else
            DO JMOL = 1, NMOLI
              PMMLI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA) =
     .         PMMLI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA)*FMOL(JMOL)
              PRFMMI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA) =
     .         PRFMMI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA)*FMOL(JMOL)
            END DO
            IF (LB_MOL.EQ.0) PMMLI(IMOL,ISTRA) = 0.0_DP
            IF (LB_MOL.EQ.0) PRFMMI(IMOL,ISTRA) = 0.0_DP
            DO JMOL = 1, NMOLI
              PMMLI(IMOL,ISTRA) = PMMLI(IMOL,ISTRA) +
     .         PMMLI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA)
              PRFMMI(IMOL,ISTRA) = PRFMMI(IMOL,ISTRA) +
     .         PRFMMI(IMOL+(JMOL-LB_MOL)*(NMOL+1),ISTRA)
            END DO
#endif
          ELSE
            PMMLI(IMOL,ISTRA)=PMMLI(IMOL,ISTRA)*FMOL(0)
            PRFMMI(IMOL,ISTRA)=PRFMMI(IMOL,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIMLI2(0:NMOL,LB_ION:NION) => PIMLI(:,ISTRA)
            PIMLI2(IMOL,1:NIONI) = PIMLI2(IMOL,1:NIONI)*FION(1:NIONI)
            PIMLI(IMOL,ISTRA) = SUM(PIMLI2(IMOL,1:NIONI))
            PRFIMI2(0:NMOL,LB_ION:NION) => PRFIMI(:,ISTRA)
            PRFIMI2(IMOL,1:NIONI) = PRFIMI2(IMOL,1:NIONI)*FION(1:NIONI)
            PRFIMI(IMOL,ISTRA) = SUM(PRFIMI2(IMOL,1:NIONI))
#else
            DO JION = 1, NIONI
              PIMLI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA) =
     .         PIMLI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA)*FION(JION)
              PRFIMI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA) =
     .         PRFIMI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA)*FION(JION)
            END DO
            IF (LB_ION.EQ.0) PIMLI(IMOL,ISTRA) = 0.0_DP
            IF (LB_ION.EQ.0) PRFIMI(IMOL,ISTRA) = 0.0_DP
            DO JION = 1, NIONI
              PIMLI(IMOL,ISTRA) = PIMLI(IMOL,ISTRA) +
     .         PIMLI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA)
              PRFIMI(IMOL,ISTRA) = PRFIMI(IMOL,ISTRA) +
     .         PRFIMI(IMOL+(JION-LB_ION)*(NMOL+1),ISTRA)
            END DO
#endif
          ELSE
            PIMLI(IMOL,ISTRA)=PIMLI(IMOL,ISTRA)*FION(0)
            PRFIMI(IMOL,ISTRA)=PRFIMI(IMOL,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHMLI2(0:NMOL,LB_PHOT:NPHOT) => PPHMLI(:,ISTRA)
            PPHMLI2(IMOL,1:NPHOTI) = PPHMLI2(IMOL,1:NPHOTI)*
     .                               FPHOT(1:NPHOTI)
            PPHMLI(IMOL,ISTRA) = SUM(PPHMLI2(IMOL,1:NPHOTI))
            PRFPHMI2(0:NMOL,LB_PHOT:NPHOT) => PRFPHMI(:,ISTRA)
            PRFPHMI2(IMOL,1:NPHOTI) = PRFPHMI2(IMOL,1:NPHOTI)*
     .                                FPHOT(1:NPHOTI)
            PRFPHMI(IMOL,ISTRA) = SUM(PRFPHMI2(IMOL,1:NPHOTI))
#else
            DO JPHOT = 1, NPHOTI
              PPHMLI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA) =
     .         PPHMLI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA)*FPHOT(JPHOT)
              PRFPHMI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA) =
     .         PRFPHMI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA)*FPHOT(JPHOT)
            END DO
            IF (LB_PHOT.EQ.0) PPHMLI(IMOL,ISTRA) = 0.0_DP
            IF (LB_PHOT.EQ.0) PRFPHMI(IMOL,ISTRA) = 0.0_DP
            DO JPHOT = 1, NPHOTI
              PPHMLI(IMOL,ISTRA) = PPHMLI(IMOL,ISTRA) +
     .         PPHMLI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA)
              PRFPHMI(IMOL,ISTRA) = PRFPHMI(IMOL,ISTRA) +
     .         PRFPHMI(IMOL+(JPHOT-LB_PHOT)*(NMOL+1),ISTRA)
            END DO
#endif
          ELSE
            PPHMLI(IMOL,ISTRA)=PPHMLI(IMOL,ISTRA)*FPHOT(0)
            PRFPHMI(IMOL,ISTRA)=PRFPHMI(IMOL,ISTRA)*FPHOT(0)
          END IF

 2116   CONTINUE

        if (logmol(0,istra)) then
          PDENMI(0,ISTRA)=SUM(PDENMI(1:NMOLI,ISTRA))
          EDENMI(0,ISTRA)=SUM(EDENMI(1:NMOLI,ISTRA))
          PGENMI(0,ISTRA)=SUM(PGENMI(1:NMOLI,ISTRA))
          EGENMI(0,ISTRA)=SUM(EGENMI(1:NMOLI,ISTRA))
          VGENMI(0,ISTRA)=SUM(VGENMI(1:NMOLI,ISTRA))
          VXDENMI(0,ISTRA)=SUM(VXDENMI(1:NMOLI,ISTRA))
          VYDENMI(0,ISTRA)=SUM(VYDENMI(1:NMOLI,ISTRA))
          VZDENMI(0,ISTRA)=SUM(VZDENMI(1:NMOLI,ISTRA))

          POTMLI(0,ISTRA)=SUM(POTMLI(1:NMOLI,ISTRA))
          EOTMLI(0,ISTRA)=SUM(EOTMLI(1:NMOLI,ISTRA))
          ERFAMI(0,ISTRA)=SUM(ERFAMI(1:NMOLI,ISTRA))
          ERFMMI(0,ISTRA)=SUM(ERFMMI(1:NMOLI,ISTRA))
          ERFIMI(0,ISTRA)=SUM(ERFIMI(1:NMOLI,ISTRA))
          ERFPHMI(0,ISTRA)=SUM(ERFPHMI(1:NMOLI,ISTRA))

          SPTAMLI(0,ISTRA)=SUM(SPTAMLI(1:NMOLI,ISTRA))
          SPTMMLI(0,ISTRA)=SUM(SPTMMLI(1:NMOLI,ISTRA))
          SPTIMLI(0,ISTRA)=SUM(SPTIMLI(1:NMOLI,ISTRA))
          SPTPHMLI(0,ISTRA)=SUM(SPTPHMLI(1:NMOLI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPMLI(0,ISTRA)=SUM(SPTPMLI(1:NMOLI,ISTRA))
          RMELI(0,ISTRA)=SUM(RMELI(1:NMOLI,ISTRA))

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAMLI2(0:NMOL,LB_ATM:NATM) => PAMLI(:,ISTRA)
            PRFAMI2(0:NMOL,LB_ATM:NATM) => PRFAMI(:,ISTRA)
            DO IATM = LB_ATM,NATMI
              PAMLI2(0,IATM) = SUM(PAMLI2(1:NMOLI,IATM))
              PRFAMI2(0,IATM) = SUM(PRFAMI2(1:NMOLI,IATM))
            END DO
#else
            DO IATM = LB_ATM,NATMI
              PAMLI(IATM*(NMOL+1),ISTRA) =
     .         SUM(PAMLI(1+IATM*(NMOL+1):NMOLI+IATM*(NMOL+1),ISTRA))
              PRFAMI(IATM*(NMOL+1),ISTRA) =
     .         SUM(PRFAMI(1+IATM*(NMOL+1):NMOLI+IATM*(NMOL+1),ISTRA))
            END DO
#endif
          ELSE
            PAMLI(0,ISTRA)=SUM(PAMLI(1:NMOLI,ISTRA))
            PRFAMI(0,ISTRA)=SUM(PRFAMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMMLI2(0:NMOL,LB_MOL:NMOL) => PMMLI(:,ISTRA)
            PRFMMI2(0:NMOL,LB_MOL:NMOL) => PRFMMI(:,ISTRA)
            DO IMOL = LB_MOL,NMOLI
              PMMLI2(0,IMOL) = SUM(PMMLI2(1:NMOLI,IMOL))
              PRFMMI2(0,IMOL) = SUM(PRFMMI2(1:NMOLI,IMOL))
            END DO
#else
            DO IMOL = LB_MOL, NMOLI
              PMMLI(IMOL*(NMOL+1),ISTRA) =
     .         SUM(PMMLI(1+IMOL*(NMOL+1):NMOLI+IMOL*(NMOL+1),ISTRA))
              PRFMMI(IMOL*(NMOL+1),ISTRA) =
     .         SUM(PRFMMI(1+IMOL*(NMOL+1):NMOLI+IMOL*(NMOL+1),ISTRA))
            END DO
#endif
          ELSE
            PMMLI(0,ISTRA)=SUM(PMMLI(1:NMOLI,ISTRA))
            PRFMMI(0,ISTRA)=SUM(PRFMMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIMLI2(0:NMOL,LB_ION:NION) => PIMLI(:,ISTRA)
            PRFIMI2(0:NMOL,LB_ION:NION) => PRFIMI(:,ISTRA)
            DO IION = LB_ION,NIONI
              PIMLI2(0,IION) = SUM(PIMLI2(1:NMOLI,IION))
              PRFIMI2(0,IION) = SUM(PRFIMI2(1:NMOLI,IION))
            END DO
#else
            DO IION = LB_ION, NIONI
              PIMLI(IION*(NMOL+1),ISTRA) =
     .         SUM(PIMLI(1+IION*(NMOL+1):NMOLI+IION*(NMOL+1),ISTRA))
              PRFIMI(IION*(NMOL+1),ISTRA) =
     .         SUM(PRFIMI(1+IION*(NMOL+1):NMOLI+IION*(NMOL+1),ISTRA))
            END DO
#endif
          ELSE
            PIMLI(0,ISTRA)=SUM(PIMLI(1:NMOLI,ISTRA))
            PRFIMI(0,ISTRA)=SUM(PRFIMI(1:NMOLI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHMLI2(0:NMOL,LB_PHOT:NPHOT) => PPHMLI(:,ISTRA)
            PRFPHMI2(0:NMOL,LB_PHOT:NPHOT) => PRFPHMI(:,ISTRA)
            DO IPHOT = LB_PHOT,NPHOTI
              PPHMLI2(0,IPHOT) = SUM(PPHMLI2(1:NMOLI,IPHOT))
              PRFPHMI2(0,IPHOT) = SUM(PRFPHMI2(1:NMOLI,IPHOT))
            END DO
#else
            DO IPHOT = LB_PHOT, NPHOTI
              PPHMLI(IPHOT*(NMOL+1),ISTRA) =
     .         SUM(PPHMLI(1+IPHOT*(NMOL+1):NMOLI+IPHOT*(NMOL+1),ISTRA))
              PRFPHMI(IPHOT*(NMOL+1),ISTRA) =
     .         SUM(PRFPHMI(1+IPHOT*(NMOL+1):NMOLI+IPHOT*(NMOL+1),ISTRA))
            END DO
#endif
          ELSE
            PPHMLI(0,ISTRA)=SUM(PPHMLI(1:NMOLI,ISTRA))
            PRFPHMI(0,ISTRA)=SUM(PRFPHMI(1:NMOLI,ISTRA))
          END IF

        end if

        sptmtti(istra) = sptmtti(istra)*fmol(0)

        DO 2117 J=1,NSBOX_TAL
          IF (LEAML) EAML(J)=EAML(J)*FATM(0)
          IF (LEMML) EMML(J)=EMML(J)*FMOL(0)
          IF (LEIML) EIML(J)=EIML(J)*FION(0)
          IF (LEPHML) EPHML(J)=EPHML(J)*FPHOT(0)
 2117   CONTINUE
        EAMLI(ISTRA)=EAMLI(ISTRA)*FATM(0)
        EMMLI(ISTRA)=EMMLI(ISTRA)*FMOL(0)
        EIMLI(ISTRA)=EIMLI(ISTRA)*FION(0)
        EPHMLI(ISTRA)=EPHMLI(ISTRA)*FPHOT(0)
C
C  TEST ION TALLIES
C
        DO 420 IION=1,NIONI
          if (.not.logion(iion,istra)) cycle
          DO 421 J=1,NSBOX_TAL
            IF (LPDENI) PDENI(IION,J)=PDENI(IION,J)*FION(IION)
            IF (LEDENI) EDENI(IION,J)=EDENI(IION,J)*FION(IION)
            IF (LPAIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAIO2(1:NION,LB_ATM:NATM) => PAIO(:,J)
                PAIO2(IION,1:NATMI) = PAIO2(IION,1:NATMI)*FATM(1:NATMI)
                PAIO(IION,J) = SUM(PAIO2(IION,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PAIO(IION+(JATM-LB_ATM)*NION,J) =
     .             PAIO(IION+(JATM-LB_ATM)*NION,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PAIO(IION,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PAIO(IION,J) = PAIO(IION,J) +
     .             PAIO(IION+(JATM-LB_ATM)*NION,J)
                END DO
#endif
              ELSE
                PAIO(IION,J)=PAIO(IION,J)*FATM(0)
              END IF
            END IF
            IF (LPMIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMIO2(1:NION,LB_MOL:NMOL) => PMIO(:,J)
                PMIO2(IION,1:NMOLI) = PMIO2(IION,1:NMOLI)*FMOL(1:NMOLI)
                PMIO(IION,J) = SUM(PMIO2(IION,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PMIO(IION+(JMOL-LB_MOL)*NION,J) =
     .             PMIO(IION+(JMOL-LB_MOL)*NION,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PMIO(IION,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PMIO(IION,J) = PMIO(IION,J) +
     .             PMIO(IION+(JMOL-LB_MOL)*NION,J)
                END DO
#endif
              ELSE
                PMIO(IION,J)=PMIO(IION,J)*FMOL(0)
              END IF
            END IF
            IF (LPIIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIIO2(1:NION,LB_ION:NION) => PIIO(:,J)
                PIIO2(IION,1:NIONI) = PIIO2(IION,1:NIONI)*FION(1:NIONI)
                PIIO(IION,J) = SUM(PIIO2(IION,1:NIONI))
#else
                DO JION = 1, NIONI
                  PIIO(IION+(JION-LB_ION)*NION,J) =
     .             PIIO(IION+(JION-LB_ION)*NION,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PIIO(IION,J) = 0.0_DP
                DO JION = 1, NIONI
                  PIIO(IION,J) = PIIO(IION,J) +
     .             PIIO(IION+(JION-LB_ION)*NION,J)
                END DO
#endif
              ELSE
                PIIO(IION,J)=PIIO(IION,J)*FION(0)
              END IF
            END IF
            IF (LPPHIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHIO2(1:NION,LB_PHOT:NPHOT) => PPHIO(:,J)
                PPHIO2(IION,1:NPHOTI) = PPHIO2(IION,1:NPHOTI)*
     .                                  FPHOT(1:NPHOTI)
                PPHIO(IION,J) = SUM(PPHIO2(IION,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PPHIO(IION+(JPHOT-LB_PHOT)*NION,J) =
     .             PPHIO(IION+(JPHOT-LB_PHOT)*NION,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PPHIO(IION,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PPHIO(IION,J) = PPHIO(IION,J) +
     .             PPHIO(IION+(JPHOT-LB_PHOT)*NION,J)
                END DO
#endif
              ELSE
                PPHIO(IION,J)=PPHIO(IION,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENI) PGENI(IION,J)=PGENI(IION,J)*FION(IION)
            IF (LEGENI) EGENI(IION,J)=EGENI(IION,J)*FION(IION)
            IF (LVGENI) VGENI(IION,J)=VGENI(IION,J)*FION(IION)
            IF (LVXDENI) VXDENI(IION,J)=VXDENI(IION,J)*FION(IION)
            IF (LVYDENI) VYDENI(IION,J)=VYDENI(IION,J)*FION(IION)
            IF (LVZDENI) VZDENI(IION,J)=VZDENI(IION,J)*FION(IION)
            IF (LRIEL) RIEL(IION,J)=RIEL(IION,J)*FION(IION)
  421     CONTINUE
          DO 422 J=1,NLMPGS
            IF (LPOTIO) POTIO(IION,J)=POTIO(IION,J)*FION(IION)
            IF (LPRFAIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAIO2(1:NION,LB_ATM:NATM) => PRFAIO(:,J)
                PRFAIO2(IION,1:NATMI) = PRFAIO2(IION,1:NATMI)*
     .                                  FATM(1:NATMI)
                PRFAIO(IION,J) = SUM(PRFAIO2(IION,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PRFAIO(IION+(JATM-LB_ATM)*NION,J) =
     .             PRFAIO(IION+(JATM-LB_ATM)*NION,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PRFAIO(IION,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PRFAIO(IION,J) = PRFAIO(IION,J) +
     .             PRFAIO(IION+(JATM-LB_ATM)*NION,J)
                END DO
#endif
              ELSE
                PRFAIO(IION,J)=PRFAIO(IION,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMIO2(1:NION,LB_MOL:NMOL) => PRFMIO(:,J)
                PRFMIO2(IION,1:NMOLI) = PRFMIO2(IION,1:NMOLI)*
     .                                  FMOL(1:NMOLI)
                PRFMIO(IION,J) = SUM(PRFMIO2(IION,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PRFMIO(IION+(JMOL-LB_MOL)*NION,J) =
     .             PRFMIO(IION+(JMOL-LB_MOL)*NION,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PRFMIO(IION,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PRFMIO(IION,J) = PRFMIO(IION,J) +
     .             PRFMIO(IION+(JMOL-LB_MOL)*NION,J)
                END DO
#endif
              ELSE
                PRFMIO(IION,J)=PRFMIO(IION,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIIO2(1:NION,LB_ION:NION) => PRFIIO(:,J)
                PRFIIO2(IION,1:NIONI) = PRFIIO2(IION,1:NIONI)*
     .                                  FION(1:NIONI)
                PRFIIO(IION,J) = SUM(PRFIIO2(IION,1:NIONI))
#else
                DO JION = 1, NIONI
                  PRFIIO(IION+(JION-LB_ION)*NION,J) =
     .             PRFIIO(IION+(JION-LB_ION)*NION,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PRFIIO(IION,J) = 0.0_DP
                DO JION = 1, NIONI
                  PRFIIO(IION,J) = PRFIIO(IION,J) +
     .             PRFIIO(IION+(JION-LB_ION)*NION,J)
                END DO
#endif
              ELSE
                PRFIIO(IION,J)=PRFIIO(IION,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHIO) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHIO2(1:NION,LB_PHOT:NPHOT) => PRFPHIO(:,J)
                PRFPHIO2(IION,1:NPHOTI) = PRFPHIO2(IION,1:NPHOTI)*
     .                                    FPHOT(1:NPHOTI)
                PRFPHIO(IION,J) = SUM(PRFPHIO2(IION,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PRFPHIO(IION+(JPHOT-LB_PHOT)*NION,J) =
     .             PRFPHIO(IION+(JPHOT-LB_PHOT)*NION,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PRFPHIO(IION,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PRFPHIO(IION,J) = PRFPHIO(IION,J) +
     .             PRFPHIO(IION+(JPHOT-LB_PHOT)*NION,J)
                END DO
#endif
              ELSE
                PRFPHIO(IION,J)=PRFPHIO(IION,J)*FPHOT(0)
              END IF
            END IF
            IF (LEOTIO) EOTIO(IION,J)=EOTIO(IION,J)*FION(IION)
            IF (LERFAIO) ERFAIO(IION,J)=ERFAIO(IION,J)*FATM(0)
            IF (LERFMIO) ERFMIO(IION,J)=ERFMIO(IION,J)*FMOL(0)
            IF (LERFIIO) ERFIIO(IION,J)=ERFIIO(IION,J)*FION(IION)
            IF (LERFPHIO) ERFPHIO(IION,J)=ERFPHIO(IION,J)*FPHOT(0)
            IF (LSPTAIO) SPTAIO(IION,J)=SPTAIO(IION,J)*FATM(0)
            IF (LSPTMIO) SPTMIO(IION,J)=SPTMIO(IION,J)*FMOL(0)
            IF (LSPTIIO) SPTIIO(IION,J)=SPTIIO(IION,J)*FION(IION)
            IF (LSPTPHIO) SPTPHIO(IION,J)=SPTPHIO(IION,J)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPIO) SPTPIO(IION,J)=SPTPIO(IION,J)*FION(IION)
            IF (LSPUMP) SPUMP(NSPAM+IION,J)=
     .                             SPUMP(NSPAM+IION,J)*FION(IION)
  422     CONTINUE
  420   CONTINUE

        if (lsptitot) sptitot = sptitot*fion(0)

        DO 431 IION=1,NIONI
          if (.not.logion(iion,istra)) cycle
          PDENII(IION,ISTRA)=PDENII(IION,ISTRA)*FION(IION)
          EDENII(IION,ISTRA)=EDENII(IION,ISTRA)*FION(IION)
          PGENII(IION,ISTRA)=PGENII(IION,ISTRA)*FION(IION)
          EGENII(IION,ISTRA)=EGENII(IION,ISTRA)*FION(IION)
          VGENII(IION,ISTRA)=VGENII(IION,ISTRA)*FION(IION)
          VXDENII(IION,ISTRA)=VXDENII(IION,ISTRA)*FION(IION)
          VYDENII(IION,ISTRA)=VYDENII(IION,ISTRA)*FION(IION)
          VZDENII(IION,ISTRA)=VZDENII(IION,ISTRA)*FION(IION)

          POTIOI(IION,ISTRA)=POTIOI(IION,ISTRA)*FION(IION)
          EOTIOI(IION,ISTRA)=EOTIOI(IION,ISTRA)*FION(IION)
          ERFAII(IION,ISTRA)=ERFAII(IION,ISTRA)*FATM(0)
          ERFMII(IION,ISTRA)=ERFMII(IION,ISTRA)*FMOL(0)
          ERFIII(IION,ISTRA)=ERFIII(IION,ISTRA)*FION(IION)
          ERFPHII(IION,ISTRA)=ERFPHII(IION,ISTRA)*FPHOT(0)
          SPTAIOI(IION,ISTRA)=SPTAIOI(IION,ISTRA)*FATM(0)
          SPTMIOI(IION,ISTRA)=SPTMIOI(IION,ISTRA)*FMOL(0)
          SPTIIOI(IION,ISTRA)=SPTIIOI(IION,ISTRA)*FION(IION)
          SPTPHIOI(IION,ISTRA)=SPTPHIOI(IION,ISTRA)*FPHOT(0)
cdr  ?? scaling with bulk flux ??
          SPTPIOI(IION,ISTRA)=SPTPIOI(IION,ISTRA)*FION(IION)
          SPUMPI(NSPAM+IION,ISTRA)=SPUMPI(NSPAM+IION,ISTRA)*FION(IION)
          RIELI(IION,ISTRA)=RIELI(IION,ISTRA)*FION(IION)

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAIOI2(0:NION,LB_ATM:NATM) => PAIOI(:,ISTRA)
            PAIOI2(IION,1:NATMI) = PAIOI2(IION,1:NATMI)*FATM(1:NATMI)
            PAIOI(IION,ISTRA) = SUM(PAIOI2(IION,1:NATMI))
            PRFAII2(0:NION,LB_ATM:NATM) => PRFAII(:,ISTRA)
            PRFAII2(IION,1:NATMI) = PRFAII2(IION,1:NATMI)*FATM(1:NATMI)
            PRFAII(IION,ISTRA) = SUM(PRFAII2(IION,1:NATMI))
#else
            DO JATM = 1, NATMI
              PAIOI(IION+(JATM-LB_ATM)*(NION+1),ISTRA) =
     .         PAIOI(IION+(JATM-LB_ATM)*(NION+1),ISTRA)*FATM(JATM)
              PRFAII(IION+(JATM-LB_ATM)*(NION+1),ISTRA) =
     .         PRFAII(IION+(JATM-LB_ATM)*(NION+1),ISTRA)*FATM(JATM)
            END DO
            IF (LB_ATM.EQ.0) PAIOI(IION,ISTRA) = 0.0_DP
            IF (LB_ATM.EQ.0) PRFAII(IION,ISTRA) = 0.0_DP
            DO JATM = 1, NATMI
              PAIOI(IION,ISTRA) = PAIOI(IION,ISTRA) +
     .         PAIOI(IION+(JATM-LB_ATM)*(NION+1),ISTRA)
              PRFAII(IION,ISTRA) = PRFAII(IION,ISTRA) +
     .         PRFAII(IION+(JATM-LB_ATM)*(NION+1),ISTRA)
            END DO
#endif
          ELSE
            PAIOI(IION,ISTRA)=PAIOI(IION,ISTRA)*FATM(0)
            PRFAII(IION,ISTRA)=PRFAII(IION,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMIOI2(0:NION,LB_MOL:NMOL) => PMIOI(:,ISTRA)
            PMIOI2(IION,1:NMOLI) = PMIOI2(IION,1:NMOLI)*FMOL(1:NMOLI)
            PMIOI(IION,ISTRA) = SUM(PMIOI2(IION,1:NMOLI))
            PRFMII2(0:NION,LB_MOL:NMOL) => PRFMII(:,ISTRA)
            PRFMII2(IION,1:NMOLI) = PRFMII2(IION,1:NMOLI)*FMOL(1:NMOLI)
            PRFMII(IION,ISTRA) = SUM(PRFMII2(IION,1:NMOLI))
#else
            DO JMOL = 1, NMOLI
              PMIOI(IION+(JMOL-LB_MOL)*(NION+1),ISTRA) =
     .         PMIOI(IION+(JMOL-LB_MOL)*(NION+1),ISTRA)*FMOL(JMOL)
              PRFMII(IION+(JMOL-LB_MOL)*(NION+1),ISTRA) =
     .         PRFMII(IION+(JMOL-LB_MOL)*(NION+1),ISTRA)*FMOL(JMOL)
            END DO
            IF (LB_MOL.EQ.0) PMIOI(IION,ISTRA) = 0.0_DP
            IF (LB_MOL.EQ.0) PRFMII(IION,ISTRA) = 0.0_DP
            DO JMOL = 1, NMOLI
              PMIOI(IION,ISTRA) = PMIOI(IION,ISTRA) +
     .         PMIOI(IION+(JMOL-LB_MOL)*(NION+1),ISTRA)
              PRFMII(IION,ISTRA) = PRFMII(IION,ISTRA) +
     .         PRFMII(IION+(JMOL-LB_MOL)*(NION+1),ISTRA)
            END DO
#endif
          ELSE
            PMIOI(IION,ISTRA)=PMIOI(IION,ISTRA)*FMOL(0)
            PRFMII(IION,ISTRA)=PRFMII(IION,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIIOI2(0:NION,LB_ION:NION) => PIIOI(:,ISTRA)
            PIIOI2(IION,1:NIONI) = PIIOI2(IION,1:NIONI)*FION(1:NIONI)
            PIIOI(IION,ISTRA) = SUM(PIIOI2(IION,1:NIONI))
            PRFIII2(0:NION,LB_ION:NION) => PRFIII(:,ISTRA)
            PRFIII2(IION,1:NIONI) = PRFIII2(IION,1:NIONI)*FION(1:NIONI)
            PRFIII(IION,ISTRA) = SUM(PRFIII2(IION,1:NIONI))
#else
            DO JION = 1, NIONI
              PIIOI(IION+(JION-LB_ION)*(NION+1),ISTRA) =
     .         PIIOI(IION+(JION-LB_ION)*(NION+1),ISTRA)*FION(JION)
              PRFIII(IION+(JION-LB_ION)*(NION+1),ISTRA) =
     .         PRFIII(IION+(JION-LB_ION)*(NION+1),ISTRA)*FION(JION)
            END DO
            IF (LB_ION.EQ.0) PIIOI(IION,ISTRA) = 0.0_DP
            IF (LB_ION.EQ.0) PRFIII(IION,ISTRA) = 0.0_DP
            DO JION = 1, NIONI
              PIIOI(IION,ISTRA) = PIIOI(IION,ISTRA) +
     .         PIIOI(IION+(JION-LB_ION)*(NION+1),ISTRA)
              PRFIII(IION,ISTRA) = PRFIII(IION,ISTRA) +
     .         PRFIII(IION+(JION-LB_ION)*(NION+1),ISTRA)
            END DO
#endif
          ELSE
            PIIOI(IION,ISTRA)=PIIOI(IION,ISTRA)*FION(0)
            PRFIII(IION,ISTRA)=PRFIII(IION,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHIOI2(0:NION,LB_PHOT:NPHOT) => PPHIOI(:,ISTRA)
            PPHIOI2(IION,1:NPHOTI) = PPHIOI2(IION,1:NPHOTI)*
     .                               FPHOT(1:NPHOTI)
            PPHIOI(IION,ISTRA) = SUM(PPHIOI2(IION,1:NPHOTI))
            PRFPHII2(0:NION,LB_PHOT:NPHOT) => PRFPHII(:,ISTRA)
            PRFPHII2(IION,1:NPHOTI) = PRFPHII2(IION,1:NPHOTI)*
     .                                FPHOT(1:NPHOTI)
            PRFPHII(IION,ISTRA) = SUM(PRFPHII2(IION,1:NPHOTI))
#else
            DO JPHOT = 1, NPHOTI
              PPHIOI(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA) =
     .         PPHIOI(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA)*FPHOT(JPHOT)
              PRFPHII(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA) =
     .         PRFPHII(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA)*FPHOT(JPHOT)
            END DO
            IF (LB_PHOT.EQ.0) PPHIOI(IION,ISTRA) = 0.0_DP
            IF (LB_PHOT.EQ.0) PRFPHII(IION,ISTRA) = 0.0_DP
            DO JPHOT = 1, NPHOTI
              PPHIOI(IION,ISTRA) = PPHIOI(IION,ISTRA) +
     .         PPHIOI(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA)
              PRFPHII(IION,ISTRA) = PRFPHII(IION,ISTRA) +
     .         PRFPHII(IION+(JPHOT-LB_PHOT)*(NION+1),ISTRA)
            END DO
#endif
          ELSE
            PPHIOI(IION,ISTRA)=PPHIOI(IION,ISTRA)*FPHOT(0)
            PRFPHII(IION,ISTRA)=PRFPHII(IION,ISTRA)*FPHOT(0)
          END IF

  431   CONTINUE

        if (logion(0,istra)) then
          PDENII(0,ISTRA)=SUM(PDENII(1:NIONI,ISTRA))
          EDENII(0,ISTRA)=SUM(EDENII(1:NIONI,ISTRA))
          PGENII(0,ISTRA)=SUM(PGENII(1:NIONI,ISTRA))
          EGENII(0,ISTRA)=SUM(EGENII(1:NIONI,ISTRA))
          VGENII(0,ISTRA)=SUM(VGENII(1:NIONI,ISTRA))
          VXDENII(0,ISTRA)=SUM(VXDENII(1:NIONI,ISTRA))
          VYDENII(0,ISTRA)=SUM(VYDENII(1:NIONI,ISTRA))
          VZDENII(0,ISTRA)=SUM(VZDENII(1:NIONI,ISTRA))

          POTIOI(0,ISTRA)=SUM(POTIOI(1:NIONI,ISTRA))
          EOTIOI(0,ISTRA)=SUM(EOTIOI(1:NIONI,ISTRA))
          ERFAII(0,ISTRA)=SUM(ERFAII(1:NIONI,ISTRA))
          ERFMII(0,ISTRA)=SUM(ERFMII(1:NIONI,ISTRA))
          ERFIII(0,ISTRA)=SUM(ERFIII(1:NIONI,ISTRA))
          ERFPHII(0,ISTRA)=SUM(ERFPHII(1:NIONI,ISTRA))
          SPTAIOI(0,ISTRA)=SUM(SPTAIOI(1:NIONI,ISTRA))
          SPTMIOI(0,ISTRA)=SUM(SPTMIOI(1:NIONI,ISTRA))
          SPTIIOI(0,ISTRA)=SUM(SPTIIOI(1:NIONI,ISTRA))
          SPTPHIOI(0,ISTRA)=SUM(SPTPHIOI(1:NIONI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPIOI(0,ISTRA)=SUM(SPTPIOI(1:NIONI,ISTRA))
          RIELI(0,ISTRA)=SUM(RIELI(1:NIONI,ISTRA))

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAIOI2(0:NION,LB_ATM:NATM) => PAIOI(:,ISTRA)
            PRFAII2(0:NION,LB_ATM:NATM) => PRFAII(:,ISTRA)
            DO IATM = LB_ATM,NATMI
              PAIOI2(0,IATM) = SUM(PAIOI2(1:NIONI,IATM))
              PRFAII2(0,IATM) = SUM(PRFAII2(1:NIONI,IATM))
            END DO
#else
            DO IATM = LB_ATM,NATMI
              PAIOI(IATM*(NION+1),ISTRA) =
     .         SUM(PAIOI(1+IATM*(NION+1):NIONI+IATM*(NION+1),ISTRA))
              PRFAII(IATM*(NION+1),ISTRA) =
     .         SUM(PRFAII(1+IATM*(NION+1):NIONI+IATM*(NION+1),ISTRA))
            END DO
#endif
          ELSE
            PAIOI(0,ISTRA)=SUM(PAIOI(1:NIONI,ISTRA))
            PRFAII(0,ISTRA)=SUM(PRFAII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMIOI2(0:NION,LB_MOL:NMOL) => PMIOI(:,ISTRA)
            PRFMII2(0:NION,LB_MOL:NMOL) => PRFMII(:,ISTRA)
            DO IMOL = LB_MOL,NMOLI
              PMIOI2(0,IMOL) = SUM(PMIOI2(1:NIONI,IMOL))
              PRFMII2(0,IMOL) = SUM(PRFMII2(1:NIONI,IMOL))
            END DO
#else
            DO IMOL = LB_MOL, NMOLI
              PMIOI(IMOL*(NION+1),ISTRA) =
     .         SUM(PMIOI(1+IMOL*(NION+1):NIONI+IMOL*(NION+1),ISTRA))
              PRFMII(IMOL*(NION+1),ISTRA) =
     .         SUM(PRFMII(1+IMOL*(NION+1):NIONI+IMOL*(NION+1),ISTRA))
            END DO
#endif
          ELSE
            PMIOI(0,ISTRA)=SUM(PMIOI(1:NIONI,ISTRA))
            PRFMII(0,ISTRA)=SUM(PRFMII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIIOI2(0:NION,LB_ION:NION) => PIIOI(:,ISTRA)
            PRFIII2(0:NION,LB_ION:NION) => PRFIII(:,ISTRA)
            DO IION = LB_ION,NIONI
              PIIOI2(0,IION) = SUM(PIIOI2(1:NIONI,IION))
              PRFIII2(0,IION) = SUM(PRFIII2(1:NIONI,IION))
            END DO
#else
            DO IION = LB_ION, NIONI
              PIIOI(IION*(NION+1),ISTRA) =
     .         SUM(PIIOI(1+IION*(NION+1):NIONI+IION*(NION+1),ISTRA))
              PRFIII(IION*(NION+1),ISTRA) =
     .         SUM(PRFIII(1+IION*(NION+1):NIONI+IION*(NION+1),ISTRA))
            END DO
#endif
          ELSE
            PIIOI(0,ISTRA)=SUM(PIIOI(1:NIONI,ISTRA))
            PRFIII(0,ISTRA)=SUM(PRFIII(1:NIONI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHIOI2(0:NION,LB_PHOT:NPHOT) => PPHIOI(:,ISTRA)
            PRFPHII2(0:NION,LB_PHOT:NPHOT) => PRFPHII(:,ISTRA)
            DO IPHOT = LB_PHOT,NPHOTI
              PPHIOI2(0,IPHOT) = SUM(PPHIOI2(1:NIONI,IPHOT))
              PRFPHII2(0,IPHOT) = SUM(PRFPHII2(1:NIONI,IPHOT))
            END DO
#else
            DO IPHOT = LB_PHOT, NPHOTI
              PPHIOI(IPHOT*(NION+1),ISTRA) =
     .         SUM(PPHIOI(1+IPHOT*(NION+1):NIONI+IPHOT*(NION+1),ISTRA))
              PRFPHII(IPHOT*(NION+1),ISTRA) =
     .         SUM(PRFPHII(1+IPHOT*(NION+1):NIONI+IPHOT*(NION+1),ISTRA))
            END DO
#endif
          ELSE
            PPHIOI(0,ISTRA)=SUM(PPHIOI(1:NIONI,ISTRA))
            PRFPHII(0,ISTRA)=SUM(PRFPHII(1:NIONI,ISTRA))
          END IF

        end if

        sptitti(istra) = sptitti(istra)*fion(0)

        DO 432 J=1,NSBOX_TAL
          IF (LEAIO) EAIO(J)=EAIO(J)*FATM(0)
          IF (LEMIO) EMIO(J)=EMIO(J)*FMOL(0)
          IF (LEIIO) EIIO(J)=EIIO(J)*FION(0)
          IF (LEPHIO) EPHIO(J)=EPHIO(J)*FPHOT(0)
  432   CONTINUE
        EAIOI(ISTRA)=EAIOI(ISTRA)*FATM(0)
        EMIOI(ISTRA)=EMIOI(ISTRA)*FMOL(0)
        EIIOI(ISTRA)=EIIOI(ISTRA)*FION(0)
        EPHIOI(ISTRA)=EPHIOI(ISTRA)*FPHOT(0)
C
C  PHOTON TALLIES
C
        DO IPHOT=1,NPHOTI
          if (.not.logphot(iphot,istra)) cycle
          DO J=1,NSBOX_TAL
            IF (LPDENPH) PDENPH(IPHOT,J)=PDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LEDENPH) EDENPH(IPHOT,J)=EDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LPAPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAPHT2(1:NPHOT,LB_ATM:NATM) => PAPHT(:,J)
                PAPHT2(IPHOT,1:NATMI) = PAPHT2(IPHOT,1:NATMI)*
     .                                  FATM(1:NATMI)
                PAPHT(IPHOT,J) = SUM(PAPHT2(IPHOT,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J) =
     .             PAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PAPHT(IPHOT,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PAPHT(IPHOT,J) = PAPHT(IPHOT,J) +
     .             PAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J)
                END DO
#endif
              ELSE
                PAPHT(IPHOT,J)=PAPHT(IPHOT,J)*FATM(0)
              END IF
            END IF
            IF (LPMPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMPHT2(1:NPHOT,LB_MOL:NMOL) => PMPHT(:,J)
                PMPHT2(IPHOT,1:NMOLI) = PMPHT2(IPHOT,1:NMOLI)*
     .                                  FMOL(1:NMOLI)
                PMPHT(IPHOT,J) = SUM(PMPHT2(IPHOT,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J) =
     .             PMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PMPHT(IPHOT,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PMPHT(IPHOT,J) = PMPHT(IPHOT,J) +
     .             PMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J)
                END DO
#endif
              ELSE
                PMPHT(IPHOT,J)=PMPHT(IPHOT,J)*FMOL(0)
              END IF
            END IF
            IF (LPIPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIPHT2(1:NPHOT,LB_ION:NION) => PIPHT(:,J)
                PIPHT2(IPHOT,1:NIONI) = PIPHT2(IPHOT,1:NIONI)*
     .                                  FION(1:NIONI)
                PIPHT(IPHOT,J) = SUM(PIPHT2(IPHOT,1:NIONI))
#else
                DO JION = 1, NIONI
                  PIPHT(IPHOT+(JION-LB_ION)*NPHOT,J) =
     .             PIPHT(IPHOT+(JION-LB_ION)*NPHOT,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PIPHT(IPHOT,J) = 0.0_DP
                DO JION = 1, NIONI
                  PIPHT(IPHOT,J) = PIPHT(IPHOT,J) +
     .             PIPHT(IPHOT+(JION-LB_ION)*NPHOT,J)
                END DO
#endif
              ELSE
                PIPHT(IPHOT,J)=PIPHT(IPHOT,J)*FION(0)
              END IF
            END IF
            IF (LPPHPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPHT2(1:NPHOT,LB_PHOT:NPHOT) => PPHPHT(:,J)
                PPHPHT2(IPHOT,1:NPHOTI) = PPHPHT2(IPHOT,1:NPHOTI)*
     .                                    FPHOT(1:NPHOTI)
                PPHPHT(IPHOT,J) = SUM(PPHPHT2(IPHOT,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J) =
     .             PPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PPHPHT(IPHOT,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PPHPHT(IPHOT,J) = PPHPHT(IPHOT,J) +
     .             PPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J)
                END DO
#endif
              ELSE
                PPHPHT(IPHOT,J)=PPHPHT(IPHOT,J)*FPHOT(0)
              END IF
            END IF
            IF (LPGENPH) PGENPH(IPHOT,J)=PGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LEGENPH) EGENPH(IPHOT,J)=EGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVGENPH) VGENPH(IPHOT,J)=VGENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVXDENPH) VXDENPH(IPHOT,J)=VXDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVYDENPH) VYDENPH(IPHOT,J)=VYDENPH(IPHOT,J)*FPHOT(IPHOT)
            IF (LVZDENPH) VZDENPH(IPHOT,J)=VZDENPH(IPHOT,J)*FPHOT(IPHOT)
          END DO
C  SURFACE-AVERAGED TALLIES
          DO J=1,NLMPGS
            IF (LPOTPHT)  POTPHT (IPHOT,J)=POTPHT (IPHOT,J)*FPHOT(IPHOT)
            IF (LPRFAPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAPHT2(1:NPHOT,LB_ATM:NATM) => PRFAPHT(:,J)
                PRFAPHT2(IPHOT,1:NATMI) = PRFAPHT2(IPHOT,1:NATMI)*
     .                                    FATM(1:NATMI)
                PRFAPHT(IPHOT,J) = SUM(PRFAPHT2(IPHOT,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PRFAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J) =
     .             PRFAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PRFAPHT(IPHOT,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PRFAPHT(IPHOT,J) = PRFAPHT(IPHOT,J) +
     .             PRFAPHT(IPHOT+(JATM-LB_ATM)*NPHOT,J)
                END DO
#endif
              ELSE
                PRFAPHT(IPHOT,J)=PRFAPHT(IPHOT,J)*FATM(0)
              END IF
            END IF
            IF (LPRFMPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMPHT2(1:NPHOT,LB_MOL:NMOL) => PRFMPHT(:,J)
                PRFMPHT2(IPHOT,1:NMOLI) = PRFMPHT2(IPHOT,1:NMOLI)*
     .                                    FMOL(1:NMOLI)
                PRFMPHT(IPHOT,J) = SUM(PRFMPHT2(IPHOT,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PRFMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J) =
     .             PRFMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PRFMPHT(IPHOT,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PRFMPHT(IPHOT,J) = PRFMPHT(IPHOT,J) +
     .             PRFMPHT(IPHOT+(JMOL-LB_MOL)*NPHOT,J)
                END DO
#endif
               ELSE
                PRFMPHT(IPHOT,J)=PRFMPHT(IPHOT,J)*FMOL(0)
              END IF
            END IF
            IF (LPRFIPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIPHT2(1:NPHOT,LB_ION:NION) => PRFIPHT(:,J)
                PRFIPHT2(IPHOT,1:NIONI) = PRFIPHT2(IPHOT,1:NIONI)*
     .                                    FION(1:NIONI)
                PRFIPHT(IPHOT,J) = SUM(PRFIPHT2(IPHOT,1:NIONI))
#else
                DO JION = 1, NIONI
                  PRFIPHT(IPHOT+(JION-LB_ION)*NPHOT,J) =
     .             PRFIPHT(IPHOT+(JION-LB_ION)*NPHOT,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PRFIPHT(IPHOT,J) = 0.0_DP
                DO JION = 1, NIONI
                  PRFIPHT(IPHOT,J) = PRFIPHT(IPHOT,J) +
     .             PRFIPHT(IPHOT+(JION-LB_ION)*NPHOT,J)
                END DO
#endif
               ELSE
                PRFIPHT(IPHOT,J)=PRFIPHT(IPHOT,J)*FION(0)
              END IF
            END IF
            IF (LPRFPHPHT) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHPHT2(1:NPHOT,LB_PHOT:NPHOT) => PRFPHPHT(:,J)
                PRFPHPHT2(IPHOT,1:NPHOTI) = FPHOT(1:NPHOTI)
                PRFPHPHT(IPHOT,J) = SUM(PRFPHPHT2(IPHOT,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PRFPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J) =
     .             PRFPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PRFPHPHT(IPHOT,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PRFPHPHT(IPHOT,J) = PRFPHPHT(IPHOT,J) +
     .             PRFPHPHT(IPHOT+(JPHOT-LB_PHOT)*NPHOT,J)
                END DO
#endif
               ELSE
                PRFPHPHT(IPHOT,J)=PRFPHPHT(IPHOT,J)*FPHOT(0)
              END IF
            END IF
            IF (LPRFPPHT) PRFPPHT(IPHOT,J)=PRFPPHT(IPHOT,J)*FPHOT(IPHOT)

            IF (LEOTPHT)  EOTPHT (IPHOT,J)=EOTPHT (IPHOT,J)*FPHOT(IPHOT)
            IF (LERFAPHT) ERFAPHT(IPHOT,J)=ERFAPHT(IPHOT,J)*FATM(0)
            IF (LERFMPHT) ERFMPHT(IPHOT,J)=ERFMPHT(IPHOT,J)*FMOL(0)
            IF (LERFIPHT) ERFIPHT(IPHOT,J)=ERFIPHT(IPHOT,J)*FION(0)
            IF (LERFPHPHT) ERFPHPHT(IPHOT,J)=
     .                                    ERFPHPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LERFPPHT) ERFPPHT(IPHOT,J)=ERFPPHT(IPHOT,J)*FPHOT(IPHOT)

            IF (LSPTAPHT) SPTAPHT(IPHOT,J)=SPTAPHT(IPHOT,J)*FATM(0)
            IF (LSPTMPHT) SPTMPHT(IPHOT,J)=SPTMPHT(IPHOT,J)*FMOL(0)
            IF (LSPTIPHT) SPTIPHT(IPHOT,J)=SPTIPHT(IPHOT,J)*FION(0)
            IF (LSPTPHPHT) SPTPHPHT(IPHOT,J)=
     .                                    SPTPHPHT(IPHOT,J)*FPHOT(IPHOT)
cdr  ?? scaling with bulk flux ??
            IF (LSPTPPHT) SPTPPHT(IPHOT,J)=SPTPPHT(IPHOT,J)*FPHOT(IPHOT)
            IF (LSPUMP) SPUMP(IPHOT,J)=SPUMP(IPHOT,J)*FPHOT(IPHOT)
          END DO
        END DO

        if (lsptphtot) sptphtot = sptphtot*fphot(0)

        DO IPHOT=1,NPHOTI
          if (.not.logphot(iphot,istra)) cycle
          PDENPHI(IPHOT,ISTRA)=PDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EDENPHI(IPHOT,ISTRA)=EDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          PGENPHI(IPHOT,ISTRA)=PGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EGENPHI(IPHOT,ISTRA)=EGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VGENPHI(IPHOT,ISTRA)=VGENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VXDENPHI(IPHOT,ISTRA)=VXDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VYDENPHI(IPHOT,ISTRA)=VYDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)
          VZDENPHI(IPHOT,ISTRA)=VZDENPHI(IPHOT,ISTRA)*FPHOT(IPHOT)

          POTPHTI(IPHOT,ISTRA)=POTPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          EOTPHTI(IPHOT,ISTRA)=EOTPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)

          ERFAPHTI(IPHOT,ISTRA)=ERFAPHTI(IPHOT,ISTRA)*FATM(0)
          ERFMPHTI(IPHOT,ISTRA)=ERFMPHTI(IPHOT,ISTRA)*FMOL(0)
          ERFIPHTI(IPHOT,ISTRA)=ERFIPHTI(IPHOT,ISTRA)*FION(0)
          ERFPHPHTI(IPHOT,ISTRA)=ERFPHPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)

          SPTAPHTI(IPHOT,ISTRA)=SPTAPHTI(IPHOT,ISTRA)*FATM(0)
          SPTMPHTI(IPHOT,ISTRA)=SPTMPHTI(IPHOT,ISTRA)*FMOL(0)
          SPTIPHTI(IPHOT,ISTRA)=SPTIPHTI(IPHOT,ISTRA)*FION(0)
          SPTPHPHTI(IPHOT,ISTRA)=SPTPHPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
cdr  ?? scaling with bulk flux ??
          SPTPPHTI(IPHOT,ISTRA)=SPTPPHTI(IPHOT,ISTRA)*FPHOT(IPHOT)
          SPUMPI(IPHOT,ISTRA)=SPUMPI(IPHOT,ISTRA)*FPHOT(IPHOT)

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAPHTI2(0:NPHOT,LB_ATM:NATM) => PAPHTI(:,ISTRA)
            PAPHTI2(IPHOT,1:NATMI) = PAPHTI2(IPHOT,1:NATMI)*
     .                               FATM(1:NATMI)
            PAPHTI(IPHOT,ISTRA) = SUM(PAPHTI2(IPHOT,1:NATMI))
            PRFAPHTI2(0:NPHOT,LB_ATM:NATM) => PRFAPHTI(:,ISTRA)
            PRFAPHTI2(IPHOT,1:NATMI) = PRFAPHTI2(IPHOT,1:NATMI)*
     .                                 FATM(1:NATMI)
            PRFAPHTI(IPHOT,ISTRA) = SUM(PRFAPHTI2(IPHOT,1:NATMI))
#else
            DO JATM = 1, NATMI
              PAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA) =
     .         PAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA)*FATM(JATM)
              PRFAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA) =
     .         PRFAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA)*FATM(JATM)
            END DO
            IF (LB_ATM.EQ.0) PAPHTI(IPHOT,ISTRA) = 0.0_DP
            IF (LB_ATM.EQ.0) PRFAPHTI(IPHOT,ISTRA) = 0.0_DP
            DO JATM = 1, NATMI
              PAPHTI(IPHOT,ISTRA) = PAPHTI(IPHOT,ISTRA) +
     .         PAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA)
              PRFAPHTI(IPHOT,ISTRA) = PRFAPHTI(IPHOT,ISTRA) +
     .         PRFAPHTI(IPHOT+(JATM-LB_ATM)*(NPHOT+1),ISTRA)
            END DO
#endif
          ELSE
            PAPHTI (IPHOT,ISTRA)=PAPHTI (IPHOT,ISTRA)*FATM(0)
            PRFAPHTI(IPHOT,ISTRA)=PRFAPHTI(IPHOT,ISTRA)*FATM(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMPHTI2(0:NPHOT,LB_MOL:NMOL) => PMPHTI(:,ISTRA)
            PMPHTI2(IPHOT,1:NMOLI) = PMPHTI2(IPHOT,1:NMOLI)*
     .                               FMOL(1:NMOLI)
            PMPHTI(IPHOT,ISTRA) = SUM(PMPHTI2(IPHOT,1:NMOLI))
            PRFMPHTI2(0:NPHOT,LB_MOL:NATM) => PRFMPHTI(:,ISTRA)
            PRFMPHTI2(IPHOT,1:NMOLI) = PRFMPHTI2(IPHOT,1:NMOLI)*
     .                                 FMOL(1:NMOLI)
            PRFMPHTI(IPHOT,ISTRA) = SUM(PRFMPHTI2(IPHOT,1:NMOLI))
#else
            DO JMOL = 1, NMOLI
              PMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA) =
     .         PMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA)*FMOL(JMOL)
              PRFMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA) =
     .         PRFMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA)*FMOL(JMOL)
            END DO
            IF (LB_MOL.EQ.0) PMPHTI(IPHOT,ISTRA) = 0.0_DP
            IF (LB_MOL.EQ.0) PRFMPHTI(IPHOT,ISTRA) = 0.0_DP
            DO JMOL = 1, NMOLI
              PMPHTI(IPHOT,ISTRA) = PMPHTI(IPHOT,ISTRA) +
     .         PMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA)
              PRFMPHTI(IPHOT,ISTRA) = PRFMPHTI(IPHOT,ISTRA) +
     .         PRFMPHTI(IPHOT+(JMOL-LB_MOL)*(NPHOT+1),ISTRA)
            END DO
#endif
          ELSE
            PMPHTI (IPHOT,ISTRA)=PMPHTI (IPHOT,ISTRA)*FMOL(0)
            PRFMPHTI(IPHOT,ISTRA)=PRFMPHTI(IPHOT,ISTRA)*FMOL(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIPHTI2(0:NPHOT,LB_ION:NION) => PIPHTI(:,ISTRA)
            PIPHTI2(IPHOT,1:NIONI) = PIPHTI2(IPHOT,1:NIONI)*
     .                               FION(1:NIONI)
            PIPHTI(IPHOT,ISTRA) = SUM(PIPHTI2(IPHOT,1:NIONI))
            PRFIPHTI2(0:NPHOT,LB_ION:NION) => PRFIPHTI(:,ISTRA)
            PRFIPHTI2(IPHOT,1:NIONI) = PRFIPHTI2(IPHOT,1:NIONI)*
     .                                 FION(1:NIONI)
            PRFIPHTI(IPHOT,ISTRA) = SUM(PRFIPHTI2(IPHOT,1:NIONI))
#else
            DO JION = 1, NIONI
              PIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA) =
     .         PIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA)*FION(JION)
              PRFIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA) =
     .         PRFIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA)*FION(JION)
            END DO
            IF (LB_ION.EQ.0) PIPHTI(IPHOT,ISTRA) = 0.0_DP
            IF (LB_ION.EQ.0) PRFIPHTI(IPHOT,ISTRA) = 0.0_DP
            DO JION = 1, NIONI
              PIPHTI(IPHOT,ISTRA) = PIPHTI(IPHOT,ISTRA) +
     .         PIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA)
              PRFIPHTI(IPHOT,ISTRA) = PRFIPHTI(IPHOT,ISTRA) +
     .         PRFIPHTI(IPHOT+(JION-LB_ION)*(NPHOT+1),ISTRA)
            END DO
#endif
          ELSE
            PIPHTI (IPHOT,ISTRA)=PIPHTI (IPHOT,ISTRA)*FION(0)
            PRFIPHTI(IPHOT,ISTRA)=PRFIPHTI(IPHOT,ISTRA)*FION(0)
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PPHPHTI(:,ISTRA)
            PPHPHTI2(IPHOT,1:NPHOTI) = PPHPHTI2(IPHOT,1:NPHOTI)*
     .                                 FPHOT(1:NPHOTI)
            PPHPHTI(IPHOT,ISTRA) = SUM(PPHPHTI2(IPHOT,1:NPHOTI))
            PRFPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PRFPHPHTI(:,ISTRA)
            PRFPHPHTI2(IPHOT,1:NPHOTI) = PRFPHPHTI2(IPHOT,1:NPHOTI)*
     .                                   FPHOT(1:NPHOTI)
            PRFPHPHTI(IPHOT,ISTRA) = SUM(PRFPHPHTI2(IPHOT,1:NPHOTI))
#else
            DO JPHOT = 1, NPHOTI
              PPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA) =
     .         PPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA)*
     .         FPHOT(JPHOT)
              PRFPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA) =
     .         PRFPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA)*
     .         FPHOT(JPHOT)
            END DO
            IF (LB_PHOT.EQ.0) PPHPHTI(IPHOT,ISTRA) = 0.0_DP
            IF (LB_PHOT.EQ.0) PRFPHPHTI(IPHOT,ISTRA) = 0.0_DP
            DO JPHOT = 1, NPHOTI
              PPHPHTI(IPHOT,ISTRA) = PPHPHTI(IPHOT,ISTRA) +
     .         PPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA)
              PRFPHPHTI(IPHOT,ISTRA) = PRFPHPHTI(IPHOT,ISTRA) +
     .         PRFPHPHTI(IPHOT+(JPHOT-LB_PHOT)*(NPHOT+1),ISTRA)
            END DO
#endif
          ELSE
            PPHPHTI(IPHOT,ISTRA)=PPHPHTI(IPHOT,ISTRA)*FPHOT(0)
            PRFPHPHTI(IPHOT,ISTRA)=PRFPHPHTI(IPHOT,ISTRA)*FPHOT(0)
          END IF

        END DO

        if (logphot(0,istra)) then
          PDENPHI(0,ISTRA)=SUM(PDENPHI(1:NPHOTI,ISTRA))
          EDENPHI(0,ISTRA)=SUM(EDENPHI(1:NPHOTI,ISTRA))
          PGENPHI(0,ISTRA)=SUM(PGENPHI(1:NPHOTI,ISTRA))
          EGENPHI(0,ISTRA)=SUM(EGENPHI(1:NPHOTI,ISTRA))
          VGENPHI(0,ISTRA)=SUM(VGENPHI(1:NPHOTI,ISTRA))
          VXDENPHI(0,ISTRA)=SUM(VXDENPHI(1:NPHOTI,ISTRA))
          VYDENPHI(0,ISTRA)=SUM(VYDENPHI(1:NPHOTI,ISTRA))
          VZDENPHI(0,ISTRA)=SUM(VZDENPHI(1:NPHOTI,ISTRA))

          POTPHTI(0,ISTRA)=SUM(POTPHTI(1:NPHOTI,ISTRA))
          EOTPHTI(0,ISTRA)=SUM(EOTPHTI(1:NPHOTI,ISTRA))
          ERFAPHTI(0,ISTRA)=SUM(ERFAPHTI(1:NPHOTI,ISTRA))
          ERFMPHTI(0,ISTRA)=SUM(ERFMPHTI(1:NPHOTI,ISTRA))
          ERFIPHTI(0,ISTRA)=SUM(ERFIPHTI(1:NPHOTI,ISTRA))
          ERFPHPHTI(0,ISTRA)=SUM(ERFPHPHTI(1:NPHOTI,ISTRA))
          SPTAPHTI(0,ISTRA)=SUM(SPTAPHTI(1:NPHOTI,ISTRA))
          SPTMPHTI(0,ISTRA)=SUM(SPTMPHTI(1:NPHOTI,ISTRA))
          SPTIPHTI(0,ISTRA)=SUM(SPTIPHTI(1:NPHOTI,ISTRA))
          SPTPHPHTI(0,ISTRA)=SUM(SPTPHPHTI(1:NPHOTI,ISTRA))
cdr  ?? scaling with bulk flux ??
          SPTPPHTI(0,ISTRA)=SUM(SPTPPHTI(1:NPHOTI,ISTRA))

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAPHTI2(0:NPHOT,LB_ATM:NATM) => PAPHTI(:,ISTRA)
            PRFAPHTI2(0:NPHOT,LB_ATM:NATM) => PRFAPHTI(:,ISTRA)
            DO IATM = LB_ATM,NATMI
              PAPHTI2(0,IATM) = SUM(PAPHTI2(1:NPHOTI,IATM))
              PRFAPHTI2(0,IATM) = SUM(PRFAPHTI2(1:NPHOTI,IATM))
            END DO
#else
            DO IATM = LB_ATM,NATMI
              PAPHTI(IATM*(NPHOT+1),ISTRA) =
     .         SUM(PAPHTI(1+IATM*(NPHOT+1):NPHOTI+IATM*(NPHOT+1),ISTRA))
              PRFAPHTI(IATM*(NION+1),ISTRA) =
     .         SUM(PRFAPHTI(1+IATM*(NPHOT+1):
     .                 NPHOTI+IATM*(NPHOT+1),ISTRA))
            END DO
#endif
          ELSE
            PAPHTI(0,ISTRA)=SUM(PAPHTI(1:NPHOTI,ISTRA))
            PRFAPHTI(0,ISTRA)=SUM(PRFAPHTI(1:NPHOTI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMPHTI2(0:NPHOT,LB_MOL:NMOL) => PMPHTI(:,ISTRA)
            PRFMPHTI2(0:NPHOT,LB_MOL:NMOL) => PRFMPHTI(:,ISTRA)
            DO IMOL = LB_MOL,NMOLI
              PMPHTI2(0,IMOL) = SUM(PMPHTI2(1:NPHOTI,IMOL))
              PRFMPHTI2(0,IMOL) = SUM(PRFMPHTI2(1:NPHOTI,IMOL))
            END DO
#else
            DO IMOL = LB_MOL, NMOLI
              PMPHTI(IMOL*(NPHOT+1),ISTRA) =
     .         SUM(PMPHTI(1+IMOL*(NPHOT+1):NPHOTI+IMOL*(NPHOT+1),ISTRA))
              PRFMPHTI(IMOL*(NPHOT+1),ISTRA) =
     .         SUM(PRFMPHTI(1+IMOL*(NPHOT+1):
     .                 NPHOTI+IMOL*(NPHOT+1),ISTRA))
            END DO
#endif
          ELSE
            PMPHTI(0,ISTRA)=SUM(PMPHTI(1:NPHOTI,ISTRA))
            PRFMPHTI(0,ISTRA)=SUM(PRFMPHTI(1:NPHOTI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIPHTI2(0:NPHOT,LB_ION:NION) => PIPHTI(:,ISTRA)
            PRFIPHTI2(0:NPHOT,LB_ION:NION) => PRFIPHTI(:,ISTRA)
            DO IION = LB_ION,NIONI
              PIPHTI2(0,IION) = SUM(PIPHTI2(1:NPHOTI,IION))
              PRFIPHTI2(0,IION) = SUM(PRFIPHTI2(1:NPHOTI,IION))
            END DO
#else
            DO IION = LB_ION, NIONI
              PIPHTI(IION*(NPHOT+1),ISTRA) =
     .         SUM(PIPHTI(1+IION*(NPHOT+1):NPHOTI+IION*(NPHOT+1),ISTRA))
              PRFIPHTI(IION*(NPHOT+1),ISTRA) =
     .         SUM(PRFIPHTI(1+IION*(NPHOT+1):
     .                 NPHOTI+IION*(NPHOT+1),ISTRA))
            END DO
#endif
          ELSE
            PIPHTI(0,ISTRA)=SUM(PIPHTI(1:NPHOTI,ISTRA))
            PRFIPHTI(0,ISTRA)=SUM(PRFIPHTI(1:NPHOTI,ISTRA))
          END IF

          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PPHPHTI(:,ISTRA)
            PRFPHPHTI2(0:NPHOT,LB_PHOT:NPHOT) => PRFPHPHTI(:,ISTRA)
            DO IPHOT = LB_PHOT,NPHOTI
              PPHPHTI2(0,IPHOT) = SUM(PPHPHTI2(1:NPHOTI,IPHOT))
              PRFPHPHTI2(0,IPHOT) = SUM(PRFPHPHTI2(1:NPHOTI,IPHOT))
            END DO
#else
            DO IPHOT = LB_PHOT, NPHOTI
              PPHPHTI(IPHOT*(NPHOT+1),ISTRA) =
     .         SUM(PPHPHTI(1+IPHOT*(NPHOT+1):
     .                NPHOTI+IPHOT*(NPHOT+1),ISTRA))
              PRFPHPHTI(IPHOT*(NPHOT+1),ISTRA) =
     .         SUM(PRFPHPHTI(1+IPHOT*(NPHOT+1):
     .                  NPHOTI+IPHOT*(NPHOT+1),ISTRA))
            END DO
#endif
          ELSE
            PPHPHTI(0,ISTRA)=SUM(PPHPHTI(1:NPHOTI,ISTRA))
            PRFPHPHTI(0,ISTRA)=SUM(PRFPHPHTI(1:NPHOTI,ISTRA))
          END IF

        END IF

        SPUMPI(0,ISTRA) = SUM(SPUMPI(1:NSPAM+NIONI,ISTRA))

        sptphtti(istra) = sptphtti(istra)*fphot(0)

        DO J=1,NSBOX_TAL
          IF (LEAPHT) EAPHT(J)=EAPHT(J)*FATM(0)
          IF (LEMPHT) EMPHT(J)=EMPHT(J)*FMOL(0)
          IF (LEIPHT) EIPHT(J)=EIPHT(J)*FION(0)
          IF (LEPHPHT) EPHPHT(J)=EPHPHT(J)*FPHOT(0)
        END DO
        EAPHTI(ISTRA)=EAPHTI(ISTRA)*FATM(0)
        EMPHTI(ISTRA)=EMPHTI(ISTRA)*FMOL(0)
        EIPHTI(ISTRA)=EIPHTI(ISTRA)*FION(0)
        EPHPHTI(ISTRA)=EPHPHTI(ISTRA)*FPHOT(0)

C  ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES
C
        IF (LADDV) THEN
          DO 423 IADV=1,NADVI
            IF (IADRC(IADV).EQ.0) THEN
              FADD=FPHOT(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.1) THEN
              FADD=FATM(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.2) THEN
              FADD=FMOL(IADVS(IADV))
            ELSEIF (IADRC(IADV).EQ.3) THEN
              FADD=FION(IADVS(IADV))
            ELSE
              CYCLE
            ENDIF
            DO J=1,NSBOX_TAL
              ADDV(IADV,J)=ADDV(IADV,J)*FADD
            ENDDO
            ADDVI(IADV,ISTRA)=ADDVI(IADV,ISTRA)*FADD
  423     CONTINUE
        END IF
C
C  ADDITIONAL COLLISION-ESTIMATED TALLIES
C
        IF (LCOLV) THEN
          DO 426 ICLV=1,NCLVI
            IF (ICLRC(ICLV).EQ.0) THEN
              FADD=FPHOT(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.1) THEN
              FADD=FATM(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.2) THEN
              FADD=FMOL(ICLVS(ICLV))
            ELSEIF (ICLRC(ICLV).EQ.3) THEN
              FADD=FION(ICLVS(ICLV))
            ELSE
              CYCLE
            ENDIF
            DO J=1,NSBOX_TAL
              COLV(ICLV,J)=COLV(ICLV,J)*FADD
            ENDDO
            COLVI(ICLV,ISTRA)=COLVI(ICLV,ISTRA)*FADD
  426     CONTINUE
        END IF
C
C  SNAPSHOT TALLIES
C
C  TO BE WRITTEN
C
C  TALLIES FOR COUPLING TO FLUID PLASMA CODE
C
        IF (LCOPV) THEN
          DO 435 ICPV=1,NCPVI
            IF (ICPRC(ICPV).EQ.1) THEN
              FADD=FATM(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.2) THEN
              FADD=FMOL(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.3) THEN
              FADD=FION(ICPVS(ICPV))
            ELSEIF (ICPRC(ICPV).EQ.0) THEN
              FADD=FPHOT(ICPVS(ICPV))
            ELSE
              CYCLE
            ENDIF
            DO J=1,NSBOX_TAL
              COPV(ICPV,J)=COPV(ICPV,J)*FADD
            ENDDO
            COPVI(ICPV,ISTRA)=COPVI(ICPV,ISTRA)*FADD
  435     CONTINUE
        END IF
C
C  TALLIES FOR BGK SELF-COLLISION ITERATIONS
C
        IF (LBGKV) THEN
          DO 437 IBGV=1,NBGVI
            IF (IBGRC(IBGV).EQ.1) THEN
              FADD=FATM(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.2) THEN
              FADD=FMOL(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.3) THEN
              FADD=FION(IBGVS(IBGV))
            ELSEIF (IBGRC(IBGV).EQ.0) THEN
              FADD=FPHOT(IBGVS(IBGV))
            ELSE
              CYCLE
            ENDIF
            DO J=1,NSBOX_TAL
              BGKV(IBGV,J)=BGKV(IBGV,J)*FADD
            ENDDO
            BGKVI(IBGV,ISTRA)=BGKVI(IBGV,ISTRA)*FADD
  437     CONTINUE
        END IF
C
C  BULK ION TALLIES
C
        DO IPLS=1,NPLSI
          if (.not.logpls(ipls,istra)) cycle
          DO J=1,NSBOX_TAL
            IF (LPAPL) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PAPL2(1:NPLS,LB_ATM:NATM) => PAPL(:,J)
                PAPL2(IPLS,1:NATMI) = PAPL2(IPLS,1:NATMI)*FATM(1:NATMI)
                PAPL(IPLS,J) = SUM(PAPL2(IPLS,1:NATMI))
#else
                DO JATM = 1, NATMI
                  PAPL(IPLS+(JATM-LB_ATM)*NPLS,J) =
     .             PAPL(IPLS+(JATM-LB_ATM)*NPLS,J)*FATM(JATM)
                END DO
                IF (LB_ATM.EQ.0) PAPL(IPLS,J) = 0.0_DP
                DO JATM = 1, NATMI
                  PAPL(IPLS,J) = PAPL(IPLS,J) +
     .             PAPL(IPLS+(JATM-LB_ATM)*NPLS,J)
                END DO
#endif
              ELSE
                PAPL(IPLS,J)=PAPL(IPLS,J)*FATM(0)
              END IF
            END IF
            IF (LPMPL) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PMPL2(1:NPLS,LB_MOL:NMOL) => PMPL(:,J)
                PMPL2(IPLS,1:NMOLI) = PMPL2(IPLS,1:NMOLI)*FMOL(1:NMOLI)
                PMPL(IPLS,J) = SUM(PMPL2(IPLS,1:NMOLI))
#else
                DO JMOL = 1, NMOLI
                  PMPL(IPLS+(JMOL-LB_MOL)*NPLS,J) =
     .             PMPL(IPLS+(JMOL-LB_MOL)*NPLS,J)*FMOL(JMOL)
                END DO
                IF (LB_MOL.EQ.0) PMPL(IPLS,J) = 0.0_DP
                DO JMOL = 1, NMOLI
                  PMPL(IPLS,J) = PMPL(IPLS,J) +
     .             PMPL(IPLS+(JMOL-LB_MOL)*NPLS,J)
                END DO
#endif
              ELSE
                PMPL(IPLS,J)=PMPL(IPLS,J)*FMOL(0)
              END IF
            END IF
            IF (LPIPL) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PIPL2(1:NPLS,LB_ION:NION) => PIPL(:,J)
                PIPL2(IPLS,1:NIONI) = PIPL2(IPLS,1:NIONI)*FION(1:NIONI)
                PIPL(IPLS,J) = SUM(PIPL2(IPLS,1:NIONI))
#else
                DO JION = 1, NIONI
                  PIPL(IPLS+(JION-LB_ION)*NPLS,J) =
     .             PIPL(IPLS+(JION-LB_ION)*NPLS,J)*FION(JION)
                END DO
                IF (LB_ION.EQ.0) PIPL(IPLS,J) = 0.0_DP
                DO JION = 1, NIONI
                  PIPL(IPLS,J) = PIPL(IPLS,J) +
     .             PIPL(IPLS+(JION-LB_ION)*NPLS,J)
                END DO
#endif
              ELSE
                PIPL(IPLS,J)=PIPL(IPLS,J)*FION(0)
              END IF
            END IF
            IF (LPPHPL) THEN
              IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PPHPL2(1:NPLS,LB_PHOT:NPHOT) => PPHPL(:,J)
                PPHPL2(IPLS,1:NPHOTI) = PPHPL2(IPLS,1:NPHOTI)*
     .                                  FPHOT(1:NPHOTI)
                PPHPL(IPLS,J) = SUM(PPHPL2(IPLS,1:NPHOTI))
#else
                DO JPHOT = 1, NPHOTI
                  PPHPL(IPLS+(JPHOT-LB_PHOT)*NPLS,J) =
     .             PPHPL(IPLS+(JPHOT-LB_PHOT)*NPLS,J)*FPHOT(JPHOT)
                END DO
                IF (LB_PHOT.EQ.0) PPHPL(IPLS,J) = 0.0_DP
                DO JPHOT = 1, NPHOTI
                  PPHPL(IPLS,J) = PPHPL(IPLS,J) +
     .             PPHPL(IPLS+(JPHOT-LB_PHOT)*NPLS,J)
                END DO
#endif
              ELSE
                PPHPL(IPLS,J)=PPHPL(IPLS,J)*FPHOT(0)
              END IF
            END IF
            IF (LEAPL) EAPL(IPLS,J)=EAPL(IPLS,J)*FATM(0)
            IF (LEMPL) EMPL(IPLS,J)=EMPL(IPLS,J)*FMOL(0)
            IF (LEIPL) EIPL(IPLS,J)=EIPL(IPLS,J)*FION(0)
            IF (LEPHPL) EPHPL(IPLS,J)=EPHPL(IPLS,J)*FPHOT(0)
            DO ICO=0,2
              KK=ICO*NPLS
              IF (LMAPL_VEC)
     .          MAPL_VEC(KK+IPLS,J)=MAPL_VEC(KK+IPLS,J)*FATM(0)
              IF (LMMPL_VEC)
     .          MMPL_VEC(KK+IPLS,J)=MMPL_VEC(KK+IPLS,J)*FMOL(0)
              IF (LMIPL_VEC)
     .          MIPL_VEC(KK+IPLS,J)=MIPL_VEC(KK+IPLS,J)*FION(0)
              IF (LMPHPL_VEC)
     .          MPHPL_VEC(KK+IPLS,J)=MPHPL_VEC(KK+IPLS,J)*FPHOT(0)
            ENDDO
          END DO
        END DO

cdr  sum over species IPLS
        DO IPLS=1,NPLSI
          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PAPLI2(0:NPLS,LB_ATM:NATM) => PAPLI(:,ISTRA)
            PAPLI2(IPLS,1:NATMI) = PAPLI2(IPLS,1:NATMI)*FATM(1:NATMI)
            PAPLI(IPLS,ISTRA) = SUM(PAPLI2(IPLS,1:NATMI))
#else
            DO JATM = 1, NATMI
              PAPLI(IPLS+(JATM-LB_ATM)*(NPLS+1),ISTRA) =
     .         PAPLI(IPLS+(JATM-LB_ATM)*(NPLS+1),ISTRA)*FATM(JATM)
            END DO
            IF (LB_ATM.EQ.0) PAPLI(IPLS,ISTRA) = 0.0_DP
            DO JATM = 1, NATMI
              PAPLI(IPLS,ISTRA) = PAPLI(IPLS,ISTRA) +
     .         PAPLI(IPLS+(JATM-LB_ATM)*(NPLS+1),ISTRA)
            END DO
#endif
          ELSE
            PAPLI(IPLS,ISTRA)=PAPLI(IPLS,ISTRA)*FATM(0)
          END IF
          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PMPLI2(0:NPLS,LB_MOL:NMOL) => PMPLI(:,ISTRA)
            PMPLI2(IPLS,1:NMOLI) = PMPLI2(IPLS,1:NMOLI)*FMOL(1:NMOLI)
            PMPLI(IPLS,ISTRA) = SUM(PMPLI2(IPLS,1:NMOLI))
#else
            DO JMOL = 1, NMOLI
              PMPLI(IPLS+(JMOL-LB_MOL)*(NPLS+1),ISTRA) =
     .         PMPLI(IPLS+(JMOL-LB_MOL)*(NPLS+1),ISTRA)*FMOL(JMOL)
            END DO
            IF (LB_MOL.EQ.0) PMPLI(IPLS,ISTRA) = 0.0_DP
            DO JMOL = 1, NMOLI
              PMPLI(IPLS,ISTRA) = PMPLI(IPLS,ISTRA) +
     .         PMPLI(IPLS+(JMOL-LB_MOL)*(NPLS+1),ISTRA)
            END DO
#endif
          ELSE
            PMPLI(IPLS,ISTRA)=PMPLI(IPLS,ISTRA)*FMOL(0)
          END IF
          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PIPLI2(0:NPLS,LB_ION:NION) => PIPLI(:,ISTRA)
            PIPLI2(IPLS,1:NIONI) = PIPLI2(IPLS,1:NIONI)*FION(1:NIONI)
            PIPLI(IPLS,ISTRA) = SUM(PIPLI2(IPLS,1:NIONI))
#else
            DO JION = 1, NIONI
              PIPLI(IPLS+(JION-LB_ION)*(NPLS+1),ISTRA) =
     .         PIPLI(IPLS+(JION-LB_ION)*(NPLS+1),ISTRA)*FION(JION)
            END DO
            IF (LB_ION.EQ.0) PIPLI(IPLS,ISTRA) = 0.0_DP
            DO JION = 1, NIONI
              PIPLI(IPLS,ISTRA) = PIPLI(IPLS,ISTRA) +
     .         PIPLI(IPLS+(JION-LB_ION)*(NPLS+1),ISTRA)
            END DO
#endif
          ELSE
            PIPLI(IPLS,ISTRA)=PIPLI(IPLS,ISTRA)*FION(0)
          END IF
          IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
            PPHPLI2(0:NPLS,LB_PHOT:NPHOT) => PPHPLI(:,ISTRA)
            PPHPLI2(IPLS,1:NPHOTI) = PPHPLI2(IPLS,1:NPHOTI)*
     .                               FPHOT(1:NPHOTI)
            PPHPLI(IPLS,ISTRA) = SUM(PPHPLI2(IPLS,1:NPHOTI))
#else
            DO JPHOT = 1, NPHOTI
              PPHPLI(IPLS+(JPHOT-LB_PHOT)*(NPLS+1),ISTRA) =
     .         PPHPLI(IPLS+(JPHOT-LB_PHOT)*(NPLS+1),ISTRA)*FPHOT(JPHOT)
            END DO
            IF (LB_PHOT.EQ.0) PPHPLI(IPLS,ISTRA) = 0.0_DP
            DO JPHOT = 1, NPHOTI
              PPHPLI(IPLS,ISTRA) = PPHPLI(IPLS,ISTRA) +
     .         PPHPLI(IPLS+(JPHOT-LB_PHOT)*(NPLS+1),ISTRA)
            END DO
#endif
          ELSE
            PPHPLI(IPLS,ISTRA)=PPHPLI(IPLS,ISTRA)*FPHOT(0)
          END IF
          EAPLI(IPLS,ISTRA)=EAPLI(IPLS,ISTRA)*FATM(0)
          EMPLI(IPLS,ISTRA)=EMPLI(IPLS,ISTRA)*FMOL(0)
          EIPLI(IPLS,ISTRA)=EIPLI(IPLS,ISTRA)*FION(0)
          EPHPLI(IPLS,ISTRA)=EPHPLI(IPLS,ISTRA)*FPHOT(0)
          DO ICO=0,2
            KK=ICO*NPLS
            MAPLI_VEC(KK+IPLS,ISTRA)=MAPLI_VEC(KK+IPLS,ISTRA)*FATM(0)
            MMPLI_VEC(KK+IPLS,ISTRA)=MMPLI_VEC(KK+IPLS,ISTRA)*FMOL(0)
            MIPLI_VEC(KK+IPLS,ISTRA)=MIPLI_VEC(KK+IPLS,ISTRA)*FION(0)
            MPHPLI_VEC(KK+IPLS,ISTRA)=MPHPLI_VEC(KK+IPLS,ISTRA)*FPHOT(0)
          ENDDO
        END DO

        IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
          PAPLI2(0:NPLS,LB_ATM:NATM) => PAPLI(:,ISTRA)
          DO IATM = LB_ATM,NATMI
            PAPLI2(0,IATM) = SUM(PAPLI2(1:NPLSI,IATM))
          END DO
#else
          DO IATM = LB_ATM,NATMI
            PAPLI(IATM*(NPLS+1),ISTRA) =
     .       SUM(PAPLI(1+IATM*(NPLS+1):NPLSI+IATM*(NPLS+1),ISTRA))
          END DO
#endif
        ELSE
          PAPLI(0,ISTRA)=SUM(PAPLI(1:NPLSI,ISTRA))
        END IF

        IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
          PMPLI2(0:NPLS,LB_MOL:NMOL) => PMPLI(:,ISTRA)
          DO IMOL = LB_MOL,NMOLI
            PMPLI2(0,IMOL) = SUM(PMPLI2(1:NPLSI,IMOL))
          END DO
#else
          DO IMOL = LB_MOL,NMOLI
            PMPLI(IMOL*(NPLS+1),ISTRA) =
     .       SUM(PMPLI(1+IMOL*(NPLS+1):NPLSI+IMOL*(NPLS+1),ISTRA))
          END DO
#endif
        ELSE
          PMPLI(0,ISTRA)=SUM(PMPLI(1:NPLSI,ISTRA))
        END IF

        IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
          PIPLI2(0:NPLS,LB_ION:NION) => PIPLI(:,ISTRA)
          DO IION = LB_ION,NIONI
            PIPLI2(0,IION) = SUM(PIPLI2(1:NPLSI,IION))
          END DO
#else
          DO IION = LB_ION,NIONI
            PIPLI(IION*(NPLS+1),ISTRA) =
     .       SUM(PIPLI(1+IION*(NPLS+1):NPLSI+IION*(NPLS+1),ISTRA))
          END DO
#endif
        ELSE
          PIPLI(0,ISTRA)=SUM(PIPLI(1:NPLSI,ISTRA))
        END IF

        IF (NLSPCSCL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
          PPHPLI2(0:NPLS,LB_PHOT:NPHOT) => PPHPLI(:,ISTRA)
          DO IPHOT = LB_PHOT,NPHOTI
            PPHPLI2(0,IPHOT) = SUM(PPHPLI2(1:NPLSI,IPHOT))
          END DO
#else
          DO IPHOT = LB_PHOT,NPHOTI
            PPHPLI(IPHOT*(NPLS+1),ISTRA) =
     .       SUM(PPHPLI(1+IPHOT*(NPLS+1):NPLSI+IPHOT*(NPLS+1),ISTRA))
          END DO
#endif
        ELSE
          PPHPLI(0,ISTRA)=SUM(PPHPLI(1:NPLSI,ISTRA))
        END IF

        EAPLI(0,ISTRA)=SUM(EAPLI(1:NPLSI,ISTRA))
        EMPLI(0,ISTRA)=SUM(EMPLI(1:NPLSI,ISTRA))
        EIPLI(0,ISTRA)=SUM(EIPLI(1:NPLSI,ISTRA))
        EPHPLI(0,ISTRA)=SUM(EPHPLI(1:NPLSI,ISTRA))
cdr  parallel component
        MAPLI_VEC(-1,ISTRA) =SUM(MAPLI_VEC(1:NPLSI,ISTRA))
        MMPLI_VEC(-1,ISTRA) =SUM(MMPLI_VEC(1:NPLSI,ISTRA))
        MIPLI_VEC(-1,ISTRA) =SUM(MIPLI_VEC(1:NPLSI,ISTRA))
        MPHPLI_VEC(-1,ISTRA)=SUM(MPHPLI_VEC(1:NPLSI,ISTRA))
cdr  perp. component
        MAPLI_VEC(-2,ISTRA) =SUM(MAPLI_VEC(NPLS+1:2*NPLSI,ISTRA))
        MMPLI_VEC(-2,ISTRA) =SUM(MMPLI_VEC(NPLS+1:2*NPLSI,ISTRA))
        MIPLI_VEC(-2,ISTRA) =SUM(MIPLI_VEC(NPLS+1:2*NPLSI,ISTRA))
        MPHPLI_VEC(-2,ISTRA)=SUM(MPHPLI_VEC(1:NPLS+1:2*NPLSI,ISTRA))
cdr  diamagn. component
        MAPLI_VEC(-3,ISTRA) =SUM(MAPLI_VEC(2*NPLS+1:3*NPLSI,ISTRA))
        MMPLI_VEC(-3,ISTRA) =SUM(MMPLI_VEC(2*NPLS+1:3*NPLSI,ISTRA))
        MIPLI_VEC(-3,ISTRA) =SUM(MIPLI_VEC(2*NPLS+1:3*NPLSI,ISTRA))
        MPHPLI_VEC(-3,ISTRA)=SUM(MPHPLI_VEC(2*NPLS+1:3*NPLSI,ISTRA))

cdr  so far: sum over vectorial components. Makes no sense. Use vector norm
cdr          but: unused anyway.
        MAPLI_VEC(0,ISTRA)=SUM(MAPLI_VEC(-3:-1,ISTRA))
        MMPLI_VEC(0,ISTRA)=SUM(MMPLI_VEC(-3:-1,ISTRA))
        MIPLI_VEC(0,ISTRA)=SUM(MIPLI_VEC(-3:-1,ISTRA))
        MPHPLI_VEC(0,ISTRA)=SUM(MPHPLI_VEC(-3:-1,ISTRA))
C
C  ELECTRON TALLIES
C
        DO 551 J=1,NSBOX_TAL
          IF (LPAEL) PAEL(J)=PAEL(J)*FATM(0)
          IF (LPMEL) PMEL(J)=PMEL(J)*FMOL(0)
          IF (LPIEL) PIEL(J)=PIEL(J)*FION(0)
          IF (LPPHEL) PPHEL(J)=PPHEL(J)*FPHOT(0)
  551   CONTINUE
        PAELI(ISTRA)=PAELI(ISTRA)*FATM(0)
        PMELI(ISTRA)=PMELI(ISTRA)*FMOL(0)
        PIELI(ISTRA)=PIELI(ISTRA)*FION(0)
        PPHELI(ISTRA)=PPHELI(ISTRA)*FPHOT(0)
        DO 552 J=1,NSBOX_TAL
          IF (LEAEL) EAEL(J)=EAEL(J)*FATM(0)
          IF (LEMEL) EMEL(J)=EMEL(J)*FMOL(0)
          IF (LEIEL) EIEL(J)=EIEL(J)*FION(0)
          IF (LEPHEL) EPHEL(J)=EPHEL(J)*FPHOT(0)
  552   CONTINUE
        EAELI(ISTRA)=EAELI(ISTRA)*FATM(0)
        EMELI(ISTRA)=EMELI(ISTRA)*FMOL(0)
        EIELI(ISTRA)=EIELI(ISTRA)*FION(0)
        EPHELI(ISTRA)=EPHELI(ISTRA)*FPHOT(0)
C
C  SPECTRUM TALLIES
C
        DO ISPC=1,NADSPC
          SELECT CASE (ESTIML(ISPC)%IPRTYP)
          CASE (0)
            FADD = FPHOT(ESTIML(ISPC)%IPRSP)
          CASE (1)
            FADD = FATM(ESTIML(ISPC)%IPRSP)
          CASE (2)
            FADD = FMOL(ESTIML(ISPC)%IPRSP)
          CASE (3)
            FADD = FION(ESTIML(ISPC)%IPRSP)
          CASE DEFAULT
            FADD = 1._DP
          END SELECT
          ESTIML(ISPC)%SPC = ESTIML(ISPC)%SPC * FADD
          ESTIML(ISPC)%SPCS = ESTIML(ISPC)%SPCS * FADD
        END DO
C
        IF (NLSPCSCL_ON) THEN
          CALL EIRENE_LEER(1)
          WRITE (iunout,*) 'RESCALING OF TRACKLENGTH TALLIES COMPLETED'
          WRITE (iunout,*) 'RESCALING FACTORS:'
          CALL EIRENE_MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .                        FATM(0), FMOL(0), FION(0), FPHOT(0))
          IF (NLSPCSCL) THEN
            IF(NATMI.GE.1)
     .       CALL EIRENE_MASRR1 ('FATM       ', FATM(1:NATMI), NATMI, 6)
            IF(NMOLI.GE.1)
     .       CALL EIRENE_MASRR1 ('FMOL       ', FMOL(1:NMOLI), NMOLI, 6)
            IF(NIONI.GE.1)
     .       CALL EIRENE_MASRR1 ('FION       ', FION(1:NIONI), NIONI, 6)
            IF(NPHOTI.GE.1)
     .       CALL EIRENE_MASRR1 ('FPHOT      ',FPHOT(1:NPHOTI),NPHOTI,6)
          END IF
          CALL EIRENE_LEER(2)
        ENDIF
C
      ENDIF

      DEALLOCATE(FATM,FMOL,FION,FPHOT)

      RETURN
      END SUBROUTINE EIRENE_SCALE_TALLIES
