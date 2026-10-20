cdr  formerly: ph_xsectph
cdr  now  sync with other routines for preparation of rates
cdr  here: ph type processes


      SUBROUTINE EIRENE_XSTPH(ipht,nrc,irph,idsc)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS
      USE EIRMOD_CCONA
      USE eirmod_cgrid, ONLY: NSBOX
      USE EIRMOD_COMPRT, ONLY: IUNOUT



      IMPLICIT NONE
      integer, intent(in) :: ipht,nrc,irph,idsc
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,
     .    nseot4,ierr,ipl0ti
      real(dp) :: factkk, ebulk

c external
      integer, external :: eirene_idez

      kk=ireacph(ipht,nrc)
      factkk=freacph(ipht,nrc)
      if(factkk == 0.) factkk=1.

      IPL0 =eirene_IDEZ(IBULKPH(ipht,nrc),3,3)
      IPL1 =eirene_IDEZ(ISCD1PH(ipht,nrc),3,3)
      IPL2 =eirene_IDEZ(ISCD2PH(ipht,nrc),3,3)
      ITYP0=eirene_IDEZ(IBULKPH(ipht,nrc),1,3)  ! always: =4, bulk
      ITYP1=eirene_IDEZ(ISCD1PH(ipht,nrc),1,3)
      ITYP2=eirene_IDEZ(ISCD2PH(ipht,nrc),1,3)

c collect atomic level information data, N0, N1, N2
crc save bulk species
      LGPHPH(ipht,idsc,1)=ipl0

CDR  1ST SECONDARY
      N1STph(irph,1) = ityp1
      N1STph(irph,2) = ipl1
      N1STph(irph,3) = 0
      IF (ityp1 < 4)
     .  N1STph(irph,3) = eirene_IDEZ(ISCD1PH(ipht,nrc),2,3)
CDR  2ND SECONDARY
      N2NDph(irph,1) = ityp2
      N2NDph(irph,2) = ipl2
      N2NDph(irph,3) = 0
      IF (ityp2 < 4)
     .  N2NDph(irph,3) = eirene_IDEZ(ISCD2PH(ipht,nrc),2,3)

cdr  CROSS-SECTION: HERE: BEAM-BEAM, NO DOPPLER FROM THERMAL MOTION
      MODCOL(7,1,IRPH)=KK
cdr  COLLISION MODEL 2: BEAM-MAXWELL, i.e. Doppler broadening included
      MODCOL(7,2,IRPH)=2
cdr  COLLISION MODEL 4: BEAM-BEAM, i.e. no Doppler broadening
      MODCOL(7,2,IRPH)=4
c
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
cdr   NSEOT4=eirene_IDEZ(ISCDE,4,5)
cdr bulk energy loss not ready
      nseot4=0
      ebulk=0.
cdr
      IF (NSEOT4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATE COEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        write (iunout,*) ' in xstph, nseot4=0 '
        IF (EBULK.LE.0.D0) THEN
          write (iunout,*) ' in xstph, nseot4=0, ebulk <= 0 '
          IF (NSTORDR >= NRAD) THEN
            write (iunout,*) ' in xstph, nstordr>nrad ',irph
            IPL0TI=MPLSTI(IPL0)
            EPLPH3(Irph,1:NSBOX,1)=1.5*TIIN(IPL0TI,1:NSBOX)
            IF (LEDRIFT) EPLPH3(Irph,1:NSBOX,1)=
     .                   EPLPH3(Irph,1:NSBOX,1)+EDRIFT(IPL0,1:NSBOX)
            NELRPH(Irph) = -3
          ELSE
            write (iunout,*) ' in xstph, nstordr<nrad '
            NELRPH(Irph) = -3
          END IF
        ELSE
          write (iunout,*) ' in xstph, nseot4=0, ebulk > 0 '
          IF (NSTORDR >= NRAD) THEN
            write (iunout,*) ' in xstph, nstordr>nrad '
            EPLPH3(Irph,1:NSBOX,1)=EBULK
            IF (LEDRIFT) EPLPH3(Irph,1:NSBOX,1)=
     .                   EPLPH3(Irph,1:NSBOX,1)+EDRIFT(IPL0,1:NSBOX)
            NELRPH(Irph) = -2
          ELSE
            NELRPH(Irph) = -2
            EPLPH3(Irph,1,1)=EBULK
            write (iunout,*) ' in xstph, nstordr<nrad '
          END IF
        ENDIF
        MODCOL(7,4,IRPH)=3
        write (iunout,*) ' in xstph, Modcol(7,4,irph) ',
     .                MODCOL(7,4,Irph)
      ELSEIF (NSEOT4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATE COEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        write (iunout,*) ' in xstph, nseot4=1 '
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            IPL0TI=MPLSTI(IPL0)
            EPLPH3(Irph,1:NSBOX,1)=1.5*TIIN(IPL0TI,1:NSBOX)
            IF (LEDRIFT) EPLPH3(Irph,1:NSBOX,1)=
     .                   EPLPH3(Irph,1:NSBOX,1)+EDRIFT(IPL0,1:NSBOX)
            NELRPH(Irph) = -3
          ELSE
            NELRPH(Irph) = -3
          END IF
        ELSE
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTPH '
          WRITE (iunout,*) 'MODIFIED TREATMENT OF photon collision '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THAN WITH T = TIIN '
          CALL EIRENE_LEER(1)
          IF (NSTORDR >= NRAD) THEN
            EPLPH3(Irph,1:NSBOX,1)=EBULK
            IF (LEDRIFT) EPLPH3(Irph,1:NSBOX,1)=
     .                   EPLPH3(Irph,1:NSBOX,1)+EDRIFT(IPL0,1:NSBOX)
            NELRPH(Irph) = -2
          ELSE
            NELRPH(Irph) = -2
            EPLPH3(Irph,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(7,4,Irph)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
c     ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.-WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.-WEIGHTED RATE
      ELSE
        IERR=5
        GOTO 996
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION

      IESTph(irph,1) = eirene_IDEZ(IESTMPH(IPHT,nrc),1,3)
      IESTph(irph,2) = eirene_IDEZ(IESTMPH(IPHT,nrc),2,3)
      IESTph(irph,3) = eirene_IDEZ(IESTMPH(IPHT,nrc),3,3)
C
cdr  not ready
c
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       WRITE (iunout,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART. BALANCE '
c       WRITE (iunout,*) 'IRCX = ',IRCX
c       WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       CALL EIRENE_LEER(1)
c       IESTCX(IRCX,1)=0
c     ENDIF
c     IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       WRITE (iunout,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM. BALANCE '
c       WRITE (iunout,*) 'IRCX = ',IRCX
c       WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       CALL EIRENE_LEER(1)
c       IESTCX(IRCX,2)=0
c     ENDIF
c     IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
c       WRITE (iunout,*) 'WARNING: COLL.EST NOT AVAILABLE FOR EN. BALANCE '
c       WRITE (iunout,*) 'IRCX = ',IRCX
c       WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
c       CALL EIRENE_LEER(1)
c       IESTCX(IRCX,3)=0
c     ENDIF
      return
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPH: EXIT CALLED'
      WRITE (iunout,*) 'NO CROSS-SECTION AVAILABLE FOR NON-DEFAULT PH'
      WRITE (iunout,*) 'KK,IPHT,IPL0 ',KK,IPHT,IPL0
      WRITE (iunout,*) 'EITHER PROVIDE CROSS-SECTION OR USE DIFFERENT'
      WRITE (iunout,*) 'POST-COLLISION SAMPLING FLAG ISCDEA'
      CALL EIRENE_EXIT_OWN(1)
      RETURN
      END SUBROUTINE EIRENE_XSTPH

C
      SUBROUTINE EIRENE_XSTPH_2(Irph,IPL)
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      integer, intent(in) :: irph,ipl
C
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'Photon REACTION NO. IRPH= ',IRPH
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'Collision WITH BULK IONS IPLS:'
      WRITE (iunout,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'
c     ITYP1=N1STX(IRCX,1)
c     ITYP2=N2NDX(IRCX,1)
c     ISPZ1=N1STX(IRCX,2)
c     ISPZ2=N2NDX(IRCX,2)
c     IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
c     IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
c     IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
c     IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
c     IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
c     IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
c     IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
c     IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
c     WRITE (iunout,*) 'IPLS= ',TEXTS(NSPAMI+IPL),'I2ND1= ',TEXTS1,
c    .                    'I2ND2= ',TEXTS2
c     CALL EIRENE_LEER(1)
      RETURN
C
      END SUBROUTINE EIRENE_XSTPH_2
