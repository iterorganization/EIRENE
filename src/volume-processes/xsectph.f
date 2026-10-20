Cdr  analogue to  xsecta, xsectm, xsecti and xsectp
cdr  Prepare processes for photon type test particles.
cdr  Printout to be sync.
cdr  And: call xstph for non-default photonic reaction options.
cdr  Unfinished code, for photons.
c
Cdr  in xsectp, there is only a call to xstrc.f for iswr=7 (irph).
cdr             needed there: a call to xstrc for iswr=6  (irrc)
cdr  
C
      SUBROUTINE EIRENE_XSECTPH
C
C  TABLE FOR CROSS-SECTION AND REACTION RATES FOR PHOTONS
C
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CTRCEI
      USE EIRMOD_PHOTON
      IMPLICIT NONE
c
cdr   PHOTON COLLISIONS, PH - type (separate from OT processes, which are
cdr                                 of H.11, H.12 type, popul. ratios)
c
      integer :: kk,iphot,idsc,nrc,ipl0,ipl1,ipl2,ityp1,ityp2,ifnd,
     .           updf,mode,irph
      EXTERNAL :: EIRENE_LEER, EIRENE_MASBOX

      DO IPHOT=1,NPHOTI
        IDSC=0
        LGPHPH(IPHOT,0,0)=0
        LGPHPH(IPHOT,0,1)=0
C
C   AT PRESENT NO DEFAULT (MINIMAL) MODEL
C
        IF (NRCPH(IPHOT).EQ.0) THEN
          NPHPHI(IPHOT)=0
C
C  NON-DEFAULT "PH" MODEL:
C
        ELSEIF(NRCPH(IPHOT) > 0) THEN
          DO NRC=1,NRCPH(IPHOT)
            KK=IREACPH(IPHOT,NRC)
            IF (ISWR(KK).NE.7) CYCLE
cdr  This reaction no. NRC is a photonic reaction for IPHOT indeed.
            IDSC=IDSC+1  ! count per line photon IPHOT
            NRPHI=NRPHI+1 ! count all ph processes. 
            IRPH=NRPHI
            LGPHPH(IPHOT,IDSC,0)=IRPH

            NREAPH(IRPH) = KK
            CALL EIRENE_XSTPH (IPHOT,NRC,IRPH,IDSC)
          ENDDO
          NPHPHI(IPHOT)=IDSC
C  NO "PH" MODEL DEFINED
        ELSE
          NPHPHI(IPHOT)=0
        ENDIF

        NPHPHIM(IPHOT)=NPHPHI(IPHOT)-1
        LGPHPH(IPHOT,0,0)=NPHPHI(IPHOT)

      ENDDO

cdr sync terminolgy with rest of code
cdr unfinished

      DO IPHOT=1,NPHOTI
C
        IF (TRCAMD) THEN
          CALL EIRENE_MASBOX ('PHOTON SPECIES IPHOT = '//TEXTS(IPHOT))
          CALL EIRENE_LEER(1)
C
          IF(NPHPHI(iphot).eq.0) then
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'NO PHOTONIC REACTION WITH BULK PARTICLES'
            CALL EIRENE_LEER(1)
          ELSE
            DO IDSC=1,NPHPHI(IPHOT)

              irph=LGPHPH(iphot,idsc,0)
              ipl0=LGPHPH(iphot,idsc,1)
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'PHOTONIC REACTION NO. IRPH= ',IRPH
              CALL EIRENE_LEER(1)

              ityp1=N1STph(irph,1)
              ipl1= N1STph(irph,2)
              ityp2=N2NDph(irph,1)
              ipl2= N2NDph(irph,2)

              write (iunout,*) 'irph,ipl0,kk'
              write (iunout,*)  irph,ipl0,kk
              write (iunout,*) 'ityp1,ipl1,ityp2,ipl2'
              write (iunout,*)  ityp1,ipl1,ityp2,ipl2
              call EIRENE_leer(1)
            enddo
          endif
        endif
      enddo

      RETURN
      END SUBROUTINE EIRENE_XSECTPH
