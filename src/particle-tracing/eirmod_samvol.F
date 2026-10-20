      MODULE EIRMOD_SAMVOL
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_PHOTON
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_SAMUSR, ONLY: EIRENE_SAMUSR, EIRENE_SAMUSR_INIT

      IMPLICIT NONE

      PUBLIC :: EIRENE_SAMVL0, EIRENE_SAMVL1,
     .          EIRENE_SAMVOL_REINIT

      REAL(DP), ALLOCATABLE, SAVE :: FREC(:,:,:), VSOURC(:,:), VSMXI(:)
      REAL(DP), ALLOCATABLE, SAVE :: RQ21(:), PS21(:)
      REAL(DP), ALLOCATABLE, SAVE :: ASIMP(:,:)
      INTEGER, ALLOCATABLE, SAVE  :: ISOURC(:,:), ICMX(:),
     .                               IFREC(:)
      INTEGER, SAVE :: ISTROLD=-1
      LOGICAL, ALLOCATABLE, SAVE :: LPLSSR(:), LATSSR(:)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(ISTROLD)
!$OMP THREADPRIVATE(FREC,VSOURC,VSMXI,RQ21,PS21,ASIMP)
!$OMP THREADPRIVATE(ISOURC,ICMX,IFREC)
!$OMP THREADPRIVATE(LPLSSR,LATSSR)
#endif

      CONTAINS

cdr Aug.  22 : notation, and bugfix: sumn,sumnt vs. sum (only affecting
cdr            diagnostic output)
cdr June  22 : vectorial momentum source rates.
cdr            as in collide, update, etc...
cdr            These momentum source tallies follow the
cdr            different sign convention that came with
cdr            this extension.
cdr Sep.  21 : a bit more and corrected documentation
cdr Nov.  19 : Nested do 6 loop: erroneous exit from loop
cdr            corrected (300919): Possible significant
cdr            effect on diagnostic output (EIO loss)
cdr Sept. 19 : added: nprt(ispz)=1 condition for bremsstrahlung (exclude molec. ions)
chf Nov.  18 : samvol_usr added, for levgeo=10 option
cdr Jan   18 : only notational change, to distinguish surface substrata from volume substrata
cdr  5.14.15 : vecusr called with ncell, and 0,0,0 (center of gravity)
cdr  2.11.14 : new function eirene_brems: bremsstrahlung in W per ion
cdr            replaces explicit expression.
cdr 21.10.14 : bug fix: spectral cut-off flag ICCT set to zero for default vol.rec (KK=0)
cdr           -->now runs again on eirene default vol.rec model.
cdr 30.10.14 : lplssr true even if npts=0, to allow setting up volume source tallies,
cdr            even if npts=0 for the vol-rec stratum

cdr  1111.07: "istep out of range" error message removed once again.
!pb  2203.07: LEVGEO=6 --> LEVGEO=10
!pb  2710.06: use flux set by user-defined sampling routine
!pb  1001.06: SUBROUTINE SAMVOL_REINIT added for reinitialsation of Eirene
!pb  1812.06: calculate bremsstrahlung
!pb  2408.06: set output values for DIWL and SHWL
cdr  2008.06: tiwl(*), ... instead of tiwl(npls),... to unify code.
cdr  0604.06: check "istep out of range" moved to correct place
c    0311.05: iplsti moved after check of validity of ipls, to produce legal exit
c             rather than code crash
C  JET 2005, PATCH 1: NEW ARGUMENTS EFWL AND SHWL IN PARAMETER LIST
c                     FOR SUBROUTINES SMVOL1 AND SMUSR1
C
C  SAMVL0:
cdr  IN CASE: NLVOL .AND. NLPLS:
cdr  scan over all IPLS=1,NPLSI and for those who have assigned RC processes to them:
C    DEFINE THE CUMULATIVE DISTRIBUTION FUNCTION
C    FREC(IPLS,IRRC,ICELL) FOR EACH VOLUME SOURCE DISTRIBUTION TABRC1, FOR SAMPLING
C    THE CELL INDEX ICELL OF THE VOLUME SOURCE PARTICLE.

C    Also the total (volume integrated) particle rate SREC, (from tabrc1(irrc))
C         the field particle energy rate EIO, (from 1.5 Ti(ipls,:) + EDRIFT(ipls,:)
C         and electron energy rate EEL, (from eelrc1(irrc))
C         are computed for diagnostics.
C
C
C    A FEW GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING
C    OF THE STARTING POINTS IN EACH CELL ARE PRE-COMPUTED
C
C    THE SOURCE STRENGTH FLUX(ISTRA) IS MODIFIED FOR THE
C    STRATA WITH NLVOL(ISTRA)=.TRUE.
C
C  SAMVL1:
C    THE INITIAL COORDINATES OF A TEST FLIGHT ARE SAMPLED,
C    AND THE CELL NUMBERS ARE COMPUTED
C
      SUBROUTINE EIRENE_SAMVL0
      IMPLICIT NONE
      REAL(DP) :: ADD, EIRENE_FTABRC1, CDYN, REC, XC, YC, ZC,
     .            ADD_VEC(0:2), BX(0:2), BY(0:2), BZ(0:2), BF,
     .            VX, VY, VZ, EIRENE_FEELRC1,
     .            VPARA, VPERA, VDIAA, EELRC,
     .            MOMPARA, MOMPERA, MOMDIAA,
     .            TOT_BREMS(NPLS), Z,
     .            BREMS, EIRENE_BREMS,
cdr  sum over sub-strata
     .            SUMNT, SUMEIT, SUMEI, SUMN,
     .            X1, Y1, X2, Y2, X3, Y3
      INTEGER :: ISTR, MXREC, MXPLS, IVOLSI, IVL,
     .           IFPLS, IIRC, IRRC, I, J, ICCT, IPLSTI, IPLSV,
     .           IR1, IR2, IP1, IP2, IT1, IT2, IR, IP, IT,
     .           ICC, IRC, ISTEP, IFRC, IND, JPLS, ISTRAI,
     .           KK, ICO
      EXTERNAL :: EIRENE_LEER, EIRENE_MASAJR, EIRENE_MASJ2R,
     .            EIRENE_MASJ3, EIRENE_MASJR2, EIRENE_MASAJR3,
     .            EIRENE_EXIT_OWN, EIRENE_BFIELD, EIRENE_VECUSR,
     .            EIRENE_FTABRC1, EIRENE_FEELRC1, EIRENE_BREMS

#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP PARALLEL DEFAULT(none)
!$OMP& SHARED(NATMI,NPLS,NPLSI,NSTRAI,NPRCI,NRAD,NSRFS,N1ST,N2ND,NSPEZ)
!$OMP& SHARED(LEVGEO,FLUX,NLVOL,NLPLS)
!$OMP& PRIVATE(MXREC,MXPLS)
#endif
      IF (.NOT.ALLOCATED(FREC)) THEN

cdr  some preparatory work, for plasma (field particle) sources (strata with NLPLS=T),
cdr  termed: "recombination", which is sometimes by abuse of language

C  LPLSSR(IPLS):
C  IDENTIFY THOSE FIELD PARTICLE SPECIES IPLS TO WHICH A
C                 VOLUME SOURCE DISTRIBUTION IS ASSIGNED

        ALLOCATE (LPLSSR(NPLSI))
        ALLOCATE (LATSSR(NATMI))
        LPLSSR = .FALSE.
        LATSSR = .FALSE.
        DO ISTR=1,NSTRAI
          IF (NLVOL(ISTR) .AND. NLPLS(ISTR)
     .        .AND. (FLUX(ISTR) > 0._DP)) THEN
            IPLS = NSPEZ(ISTR)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
c  nspez out of range for at least one stratum: Set volumetric sources for ALL species
              LPLSSR = .TRUE.
            ELSE
              LPLSSR(IPLS) = .TRUE.
            END IF
          END IF
        END DO

        MXREC=MAXVAL(NPRCI(1:NPLSI))
        MXPLS=COUNT(LPLSSR(1:NPLSI))
c  cumulated distribution of cell number index NCELL
        ALLOCATE (FREC(0:MXPLS,0:MXREC,0:NRAD))
c        ALLOCATE (FREC(0:NPLSI,0:MXREC,0:NRAD))
c  same, but compressed: only non-zero entries thereof, for faster NCELL sampling by inversion method
        ALLOCATE (VSOURC(NSRFS,0:NRAD))
        ALLOCATE (VSMXI(NSRFS))
        IF (LEVGEO == 2) ALLOCATE (RQ21(N1ST))
        IF (LEVGEO == 2) ALLOCATE (PS21(N2ND))
        IF (LEVGEO == 3) ALLOCATE (ASIMP(2,NRAD))
        ALLOCATE (ISOURC(NSRFS,0:NRAD))
        ALLOCATE (ICMX(NSRFS))
        ALLOCATE (IFREC(NPLS))
      END IF
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP END PARALLEL
#endif
C
      FREC=0.
      SREC=0.
      IFREC=0
C
cdr  prepare volumetric source distribution
cdr  for primary field particles ITYP_PRIM = 4 with LPLSSR(IPLS) = T

cdr  select case ityp_prim
cdr  case = 4

      IFPLS=0
      DO 2 JPLS=1,NPLSI
        IPLS=JPLS
        IF (LGPRC(IPLS,0).EQ.0) GOTO 2
        IF (.NOT.LPLSSR(IPLS)) GOTO 2
        IFPLS=IFPLS+1
        IFREC(IPLS)=IFPLS
        DO 3 IIRC=1,NPRCI(IPLS)
          IRRC=LGPRC(IPLS,IIRC)
          KK=NREARC(IRRC)
          ICCT=0
C  SPECTRAL CUT-OFF FOR SOURCE RATE: ONLY FOR PHOTONS SO FAR.
          IF (KK.GT.0) THEN
            ICCT=NREACT(KK)
          ENDIF
          DO J=1,NSBOX
            ADD=0.
C  EXCLUDE DEAD CELLS (GRID CUTS, ISOLATED CELLS FROM COUPLE_.., ETC)
C  EXCLUDE IPLS-VACUUM CELLS

cdr
C  tabrc1(irrc,:) is a volumetric rate per ion  (1/s)
c  turn this into rate in (amp per ion), factor 'elch'
c  and then into source rate: amp per cell, factor di(ipls) * vol
c  This then is the (discrete) cell number distribution FREC (cumulated)
c  for this process irrc
cdr
            IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
              IF (NSTORDR >= NRAD) THEN
                ADD=TABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              ELSE
                ADD=EIRENE_FTABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              END IF
            END IF
C  SPECTRAL CUT-OFF FOR SOURCE RATE (ONLY USED FOR PHOTONS SO FAR)
            IF (ICCT > 0)
     .        ADD = ADD*(XINTLEFT(ICCT,J) +
     .                   XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

            FREC(IFPLS,IIRC,J)  =FREC(IFPLS,IIRC,J-1)+ADD
            SREC(IPLS,IRRC)     =SREC(IPLS,IRRC)+ADD
          END DO
    3   CONTINUE
    2 CONTINUE

C  SUM OVER SPECIES AND RECOMBINATION TYPE INDICES
      DO 4 JPLS=1,NPLSI
        IF (LGPRC(JPLS,0).EQ.0) GOTO 4
        IF (.NOT.LPLSSR(JPLS)) GOTO 4
        IFPLS=IFREC(JPLS)
        DO 5 IIRC=1,NPRCI(JPLS)
          IRRC=LGPRC(JPLS,IIRC)
          SREC(JPLS,0)=SREC(JPLS,0)+SREC(JPLS,IRRC)
          SREC(0,IRRC)=SREC(0,IRRC)+SREC(JPLS,IRRC)
          SREC(0,0)   =SREC(0,0)   +SREC(JPLS,IRRC)
          DO J=1,NSBOX
            FREC(IFPLS,0,J)=FREC(IFPLS,0,J)+FREC(IFPLS,IIRC,J)
          END DO
    5   CONTINUE
    4 CONTINUE
C
C
      IF (TRCSOU.AND.IFPLS.GT.0) THEN
cdr  for diagnostics only: further global primary volume source rates
        EIO=0.
        EEL=0.
        MOM_VEC=0.
C
        DO 7 JPLS=1,NPLSI
          IPLS=JPLS
          CDYN=CNDYNP(IPLS)
          IF (LGPRC(IPLS,0).EQ.0) GOTO 7
          IF (.NOT.LPLSSR(IPLS)) GOTO 7
          IFPLS=IFREC(IPLS)
          IPLSTI = MPLSTI(IPLS)
          IPLSV = MPLSV(IPLS)
          DO 6 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)

            ICCT=0
            IF (KK.GT.0) THEN
              ICCT=NREACT(KK)
            ENDIF
            DO J=1,NSBOX
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
c  FREC is in Amp, so ADD is in: eV * Amp = Watt
                REC=FREC(IFPLS,IIRC,J)-FREC(IFPLS,IIRC,J-1)
                IF (REC.LE.0.D0) CYCLE
                ADD=1.5*TIIN(IPLSTI,J)*REC
                IF (LEDRIFT) ADD=ADD+EDRIFT(IPLS,J)*REC
C  SPECTRAL CUT-OFF, CURRENTLY ONLY FOR PHOTONS
                IF (ICCT > 0)
     .            ADD = ADD*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

                EIO(IPLS,IRRC)=EIO(IPLS,IRRC)-ADD
                EIO(IPLS,0)   =EIO(IPLS,0   )-ADD

CDR  position x0,y0,z0 is not yet known here
cdr  take center of gravity in cell, if needed (last parameter (logical) in bfield.f)
                xc=0.
                yc=0.
                zc=0.
                mompara=0.
                mompera=0.
                momdiaa=0.
                CALL EIRENE_BFIELD (J, XC,YC,ZC, BX,BY,BZ, BF,.FALSE.)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR(2,J,XC,YC,ZC,VX,VY,VZ,IPLS,.FALSE.)
                  VPARA=VX*BX(0)+VY*BY(0)+VZ*BZ(0)
                  MOMPARA=VPARA*CDYN
                  VPERA=VX*BX(1)+VY*BY(1)+VZ*BZ(1)
                  MOMPERA=VPERA*CDYN
                  VDIAA=VX*BX(2)+VY*BY(2)+VZ*BZ(2)
                  MOMDIAA=VDIAA*CDYN
                ELSE IF (INDPRO(5) == 8) THEN
                  VX = VXIN(IPLSV,J)
                  VY = VYIN(IPLSV,J)
                  VZ = VZIN(IPLSV,J)
                  VPARA=VX*BX(0)+VY*BY(0)+VZ*BZ(0)
                  MOMPARA=VPARA*CDYN
                  VPERA=VX*BX(1)+VY*BY(1)+VZ*BZ(1)
                  MOMPERA=VPERA*CDYN
                  VDIAA=VX*BX(2)+VY*BY(2)+VZ*BZ(2)
                  MOMDIAA=VDIAA*CDYN
                ELSE
cdr  check if storage for pmom_vec(:,:) is available, else: use zero
                  IF (LPMOM_VEC) THEN
                    MOMPARA=PMOM_VEC(IPLS,J)
                    MOMPERA=PMOM_VEC(IPLS+NPLS,J)
                    MOMDIAA=PMOM_VEC(IPLS+2*NPLS,J)
                  ENDIF
                ENDIF
                ADD_VEC(0)=MOMPARA*REC
                ADD_VEC(1)=MOMPERA*REC
                ADD_VEC(2)=MOMDIAA*REC

                IF (ICCT > 0)
     .            ADD_VEC(0:2) = ADD_VEC(0:2)*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))
                DO ICO=0,2
                  KK=ICO*NPLS+IPLS
                  MOM_VEC(KK,IRRC)=MOM_VEC(KK,IRRC)-ADD_VEC(ICO)
                  MOM_VEC(KK,0)   =MOM_VEC(KK,0   )-ADD_VEC(ICO)
                ENDDO
              ENDIF
            END DO
    6     CONTINUE
C
C  associated electron cooling/heating rate: eelrc: EV *CM**3/S
          DO 8 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)

            ICCT=0
            IF (KK.GT.0) THEN
              ICCT=NREACT(KK)
            ENDIF
            DO J=1,NSBOX
              ADD=0.D0
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
                IF (NSTORDR >= NRAD) THEN
                  EELRC = EELRC1(IRRC,J)
                ELSE
                  EELRC = EIRENE_FEELRC1(IRRC,J)
                END IF
c  Turn eV/s/particle into Watt/cell
                ADD=EELRC*DIIN(IPLS,J)*VOL(J)*ELCHA
C  SPECTRAL CUT-OFF (PHOTONS ONLY)
                IF (ICCT > 0)
     .            ADD = ADD*(XINTLEFT(ICCT,J) +
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))

              ENDIF

              EEL(IPLS,IRRC)=EEL(IPLS,IRRC)+ADD
              EEL(IPLS,0   )=EEL(IPLS,0   )+ADD
            ENDDO   !  nsbox loop
    8     CONTINUE  !  irrc loop

cdr  testing internal CR model, using amjuel and h_colrad rates, nrrc=2,
cdr  with scaling factor 0.5 each. ....TEST OK, FEB 18, out again.
cdr   if (ipls.eq.1) then
c          do j=1,nsbox
c            write (iunout,*) j,eelrc1(1,j),eelrc1(2,j),tein(j),
c    .         dein(j),lgvac(j,1),nstgrd(j)
c          enddo
cdr   endif

    7   CONTINUE    !  npls loop

C  BREMSSTRAHLUNG ORIGINATING FROM IONS IPLS, CHARGE Z=NCHRGP(IPLS)
cdr needed for diagnostics only, and also for electron energy loss rates,
cdr in which bremsstrahlung might be either lumped into eelrc1 (ADAS) or not (ELSE), depending on data source.

C  only: atomic ions. Exclude here, for the time being: molecular ions
        TOT_BREMS = 0._DP
        DO JPLS=1,NPLSI
          IPLS=JPLS
          ISPZ=NSPAMI+IPLS
          IF (NCHRGP(IPLS) == 0.OR.NPRT(ISPZ).GT.1) CYCLE
          DO J = 1, NSBOX
cnh         28.10.2019
            IF (ZIIN(IPLS,J).NE.ZVAC) THEN
              Z = ZIIN(IPLS,J)
            ELSE
              Z = DBLE(NCHRGP(IPLS))
            ENDIF
            IF (LGVAC(J,NPLS+1).OR.LGVAC(J,IPLS)) CYCLE
            BREMS=EIRENE_BREMS(TEIN(J),DEIN(J),Z)  ! Watt per ion
            TOT_BREMS(IPLS) = TOT_BREMS(IPLS) +
     .                        BREMS*DIIN(IPLS,J)*VOL(J) ! Watt per cell
          END DO
        END DO
C
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'DIAGNOSTICS FROM SUBR. SAMVL0: '
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'VOLUME RECOMBINATION RATES INTEGRATED OVER'
        WRITE (iunout,*) 'ENTIRE COMPUTATIONAL GRID '
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION PARTICLE LOSS (AMP): '
        ITYP=4
        DO 10 JPLS=1,NPLSI
          IPLS=JPLS
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 11 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,SREC          ',
     .                   TEXTS(ISPZ),IRRC,-SREC(IPLS,IRRC))
   11     CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,SREC(IPLS,0)  ',
     .                   TEXTS(ISPZ),0   ,-SREC(IPLS,0))
          ENDIF
   10   CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION ENERGY LOSS (WATT): '
        DO 12 JPLS=1,NPLSI
          IPLS=JPLS
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 13 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,EIO           ',
     .                   TEXTS(ISPZ),IRRC,EIO(IPLS,IRRC))
   13     CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,EIO(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EIO(IPLS,0))
          ENDIF
   12   CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ELECTRON ENERGY LOSS (WATT): '
        DO 14 JPLS=1,NPLSI
          IPLS=JPLS
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 15 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR('IPLS,IRRC,EEL           ',
     .                   TEXTS(ISPZ),IRRC,EEL(IPLS,IRRC))
   15     CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR('IPLS,TOT.,EEL(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EEL(IPLS,0))
          ENDIF
   14   CONTINUE
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'RECOMBINATION PARALLEL MOMENTUM LOSS : '
        DO 16 JPLS=1,NPLSI
          IPLS=JPLS
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 17 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL EIRENE_MASAJR3('IPLS,IRRC,MOM_VEC       ',
     .               TEXTS(ISPZ),IRRC,MOM_VEC(IPLS,IRRC),
     .               MOM_VEC(IPLS+NPLS,IRRC),MOM_VEC(IPLS+2*NPLS,IRRC))
   17     CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL EIRENE_MASAJR3('IPLS,TOT.,MOM_VEC(IPL,0)',
     .               TEXTS(ISPZ),0   ,MOM_VEC(IPLS,0),
     .               MOM_VEC(IPLS+NPLS,0),MOM_VEC(IPLS+2*NPLS,0))
          ENDIF
   16   CONTINUE
        CALL EIRENE_LEER(1)

        WRITE (iunout,*) 'BREMSSTRAHLUNG (WATT): '
        DO JPLS=1,NPLSI
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,JPLS)
          CALL EIRENE_MASAJR('IPLS,TOT.BREMSSTRAHLUNG ',
     .                 TEXTS(ISPZ),JPLS,TOT_BREMS(JPLS))
        END DO

      ENDIF    !trcsou
C
C  SET TOTAL SOURCE STRENGTH FOR STRATA WITH NLVOL(ISTRA)=.TRUE.,
C
      DO 50 ISTRAI=1,NSTRAI
        ISTRA=ISTRAI
        IF (NLVOL(ISTRA).AND.NLPLS(ISTRA).AND.NPTS(ISTRA).GT.0) THEN
          IPLS=NSPEZ(ISTRA)
          IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
            WRITE (iunout,*) 'SOURCE SPECIES INDEX NSPEZ OUT OF RANGE'
            WRITE (iunout,*) 'ISTRA, NSPEZ(ISTRA)   ',
     .                        ISTRA, NSPEZ(ISTRA)
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          IPLSTI = MPLSTI(IPLS)

          SUMNT=0.D0
          SUMEIT=0.D0
C  VOLUMETRIC SUB-STRATA
          DO 53 IVOLSI=1,NSRFSI(ISTRA)
            IVL=IVOLSI
            SUMN=0.D0
            SUMEI=0.D0
            IF (SORLIM(IVL,ISTRA).LT.0) THEN
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR USER SPECIFIED VOLUME SOURCE
              CALL EIRENE_SAMUSR_INIT(IVL,ISTRA,
     .                    SORAD1(IVL,ISTRA),SORAD2(IVL,ISTRA),
     .                    SORAD3(IVL,ISTRA),SORAD4(IVL,ISTRA),
     .                    SORAD5(IVL,ISTRA),SORAD6(IVL,ISTRA))
!pb assume flux is set in samusr
cdr April 22 bug fix. do not overwrite total flux sumnt (old: summ) here.
              SUMN=FLUX(ISTRA)
            ELSE
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR EXTERNAL VOLUMETRIC SOURCES (use the rates: TABRC1(irrc,:))
C  ACCOUNT FOR INGRDA(IVOLSI,ISTRA,...), INGRDE(IVOLSI,ISTRA,...) section of comp. grid
              I=ISTRA
              ICC=0
              IRC=-1
              IF (NR1ST.GT.1) THEN
               IF (INGRDA(IVL,I,1).LE.0 .OR.
     .             INGRDE(IVL,I,1).LE.0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,1),INGRDE(.,.,1)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
                INGRDA(IVL,I,1)=1
                INGRDE(IVL,I,1)=MAX0(1,NR1ST)
                CALL EIRENE_LEER(1)
               ENDIF
              ENDIF
              IF (NP2ND.GT.1) THEN
               IF (INGRDA(IVL,I,2).LE.0 .OR.
     .             INGRDE(IVL,I,2).LE.0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,2),INGRDE(.,.,2)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
                INGRDA(IVL,I,2)=1
                INGRDE(IVL,I,2)=MAX0(1,NP2ND)
                CALL EIRENE_LEER(1)
               ENDIF
              ENDIF
              IF (NT3RD.GT.1) THEN
               IF (INGRDA(IVL,I,3).LE.0 .OR.
     .             INGRDE(IVL,I,3).LE.0) THEN
                CALL EIRENE_LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*)
     .            'NEW INPUT FOR INGRDA(.,.,3),INGRDE(.,.,3)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
                INGRDA(IVL,I,3)=1
                INGRDE(IVL,I,3)=MAX0(1,NT3RD)
                CALL EIRENE_LEER(1)
               ENDIF
              ENDIF
              IF (NPRCI(IPLS).EQ.0) THEN
                WRITE (iunout,*) 'NO VOLUMETRIC SOURCE DISTRIBUTION'
                WRITE (iunout,*) 'DEFINED. SUBSTRATUM TURNED OFF'
                WRITE (iunout,*) 'IPLS,IVOLSI,ISTRA ',IPLS,IVOLSI,ISTRA
                SORWGT(IVL,ISTRA)=0.D0
                GOTO 53
              ENDIF
              IF (NLRAD) THEN
                IR1=MAX0(1,INGRDA(IVL,ISTRA,1))
                IR2=MIN0(NR1ST,INGRDE(IVL,ISTRA,1))
              ELSE
                IR1=1
                IR2=2
              ENDIF
              IF (NLPOL) THEN
                IP1=MAX0(1,INGRDA(IVL,ISTRA,2))
                IP2=MIN0(NP2ND,INGRDE(IVL,ISTRA,2))
              ELSE
                IP1=1
                IP2=2
              ENDIF
              IF (NLTOR) THEN
                IT1=MAX0(1,INGRDA(IVL,ISTRA,3))
                IT2=MIN0(NT3RD,INGRDE(IVL,ISTRA,3))
              ELSE
                IT1=1
                IT2=2
              ENDIF

              ISTEP=NINT(SORIND(IVL,ISTRA))
              IFPLS=IFREC(IPLS)
              DO 52 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF (ISTEP.EQ.IRRC) THEN
C  ONE SINGLE VOLUMETRIC PROCESS RATE IRRC IDENTIFIED
                  IRC=IRRC
                  IFRC=IIRC
                ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL RECOMBINATION PROCESSES FOR SPECIES IPLS
                  IRC=0
                  IFRC=0
                  ISTEP=-1
                ELSE
C  TRY OTHER RECOMBINATION PROCESS ASSIGNED TO IPLS
                  GOTO 52
                ENDIF
                DO IR=IR1,IR2-1
                  DO IP=IP1,IP2-1
                    DO IT=IT1,IT2-1
                      NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      REC=FREC(IFPLS,IFRC,NCELL)-
     .                    FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                      IF (REC.GT.0.D0) THEN
                        ICC=ICC+1
                        SUMN=SUMN+REC
                        SUMEI=SUMEI-1.5*TIIN(IPLSTI,NCELL)*REC
                        IF (LEDRIFT) SUMEI=SUMEI-EDRIFT(IPLS,NCELL)*REC
                      ENDIF
                    END DO
                  END DO
                END DO
   52         CONTINUE   ! summing over irrc
c
              IF (SUMN.EQ.0.D0) THEN
                WRITE (IUNOUT,*) 'NO VOLUMETRIC SOURCE RATE FOR:'
                WRITE (IUNOUT,*) 'ISTRA, IVOLSI, IPLS, ISTEP ',
     .                            ISTRA, IVL   , IPLS, ISTEP, 
     .                            SUM(DIIN(IPLS,1:NSBOX))
                WRITE (IUNOUT,*) 'EITHER: ISTEP OUT OF RANGE IN SAMVOL'
                WRITE (IUNOUT,*) 'OR: DENSITY OF "RECOMBINING" IPLS = 0'
                SORWGT(IVL,ISTRA)=0.D0
                GOTO 53
              ENDIF

              SORWGT(IVL,ISTRA)=SUMN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'SUB-STRATUM WEIGHT REDEFINED'
              CALL EIRENE_MASJ2R
     .          ('IVOLSI,ISTRA,SORWGT     ',IVOLSI,ISTRA,SUMN)
              IF (TRCSOU) THEN
                CALL EIRENE_MASJ3 ('IRRC,IPLS,ICMX          ',
     .                              IRC ,IPLS,ICC)
                CALL EIRENE_LEER(1)
              ENDIF

cdr  sum over sub-strata
              SUMNT=SUMNT+SUMN
              SUMEIT=SUMEIT+SUMEI
            ENDIF
   53     CONTINUE
C
          IF (SUMNT.GT.0.D0) THEN
            FLUX(ISTRA)=SUMNT
            WRITE (iunout,*) 'SOURCE STRENGTH REDEFINED'
            CALL EIRENE_MASJR2('ISTRA, FLUX, EIFLUX     ',
     .                          ISTRA, FLUX(ISTRA), SUMEIT)
            CALL EIRENE_LEER(1)
          ELSE
crc now requires 2 iterations before stratum is off, for bgk with 1 iteration
crc only relevant for coupled runs
crc better solutions probably exist, quick fix for now
            IF (FLUX(ISTRA).GE.1.0D0) THEN
                WRITE (iunout,*) 'SOURCE STRENGTH ZERO ON ',
     .            'FIRST PASS, SETTING TO 0.5'
                CALL EIRENE_LEER(1)
                FLUX(ISTRA)=0.5D0
            ELSE
                FLUX(ISTRA)=0.D0
                WRITE (iunout,*) 'SOURCE ISTRA= ',ISTRA,' TURNED OFF '
                CALL EIRENE_LEER(1)
            END IF    
          ENDIF
        ENDIF
   50 CONTINUE
C
C  PREPARE SOME GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING IN STANDARD MESH CELLS
      select case (LEVGEO)
      case (2)
        IF (NLPOL) THEN
          DO 54 IP=1,NP2NDM
            PS21(IP)=PSURF(IP+1)-PSURF(IP)
   54     CONTINUE
        ENDIF
        DO 55 IR=1,NR1STM
          RQ21(IR)=RQ(IR+1)-RQ(IR)
   55   CONTINUE

      case (3)
c  split quadrangle into two triangles,
c  then 1st: sample triangle according to its relative area,
c  then 2nd: sample uniformly within this triangle
        IT=1
        DO IR=1,NR1ST-1
         DO IP=1,NP2ND-1
          IND=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          X1=XPOL(IR,IP)
          X2=XPOL(IR,IP+1)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR,IP)
          Y2=YPOL(IR,IP+1)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(1,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
          X1=XPOL(IR+1,IP)
          X2=XPOL(IR,IP)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR+1,IP)
          Y2=YPOL(IR,IP)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(2,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
         END DO
        END DO
      end select
C
      RETURN
      END SUBROUTINE EIRENE_SAMVL0
C
C  AT THIS POINT: CALLED FROM PARTICLE LOOP TO INITIALIZE TEST FLIGHT
C
      SUBROUTINE EIRENE_SAMVL1
     . (ISTRAI,NVLM,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,ZIWL,WEISPZ)
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: TEWL, SHWL, TIWL(*), DIWL(*),
     .                         VXWL(*), VYWL(*), VZWL(*),
     .                         EFWL(*), ZIWL(*), WEISPZ(*)
      INTEGER, INTENT(IN) :: ISTRAI, NVLM
      REAL(DP) :: ADD, X1, Y1, X2, Y2, X3, Y3, ZEP1, ZRM1, WINK,
     .            RR, RRI, RRD, RRN, EPR, ELR, Z1, Z2, Z3, ZZ,
     .            X4, Y4, Z4, X01, CNORM
      INTEGER :: IFPLS, IIRC, IRRC, IVOLSI, IVL, ICC,
     .           IR1, IR2, IP1, IP2, IT1, IT2, IR, IP, IT, ISTEP,
     .           IFRC, IAUSR, IBUSR, IRUSR, ITUSR, IPUSR,
     .           IC1, IC2, IL, IU, IM, ICELL, IN, JSPZ
      EXTERNAL :: EIRENE_FPOLYT_3, EIRENE_FPOLYT_4, EIRENE_FZRTRI,
     .            EIRENE_NCELLN, EIRENE_SAMVOL_USR, EIRENE_EXIT_OWN

C  USER-SUPPLIED SOURCE
C
      IF (SORLIM(NVLM,ISTRAI).LT.0) THEN
        CALL EIRENE_SAMUSR(NVLM,X0,Y0,Z0,
     .              SORAD1(NVLM,ISTRAI),SORAD2(NVLM,ISTRAI),
     .              SORAD3(NVLM,ISTRAI),SORAD4(NVLM,ISTRAI),
     .              SORAD5(NVLM,ISTRAI),SORAD6(NVLM,ISTRAI),
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,ZIWL,
     .              WEISPZ)
        NRCELL=IRUSR
        NPCELL=IPUSR
        NTCELL=ITUSR
        NACELL=IAUSR
        NBLOCK=IBUSR
        NBLCKA=NSTRD*(NBLOCK-1)+NACELL
        NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C
        MTSURF=0
        NLSRFZ=.FALSE.
        MPSURF=0
        NLSRFY=.FALSE.
        MRSURF=0
        NLSRFX=.FALSE.
        EFWL(1:NPLSI)=0._DP
        RETURN
      ENDIF
C
C  VOLUME RECOMBINATION SOURCE
C
C  TENTATIVELY ASSUME: A BULK ION WILL BE GENERATED
      LGPART=.TRUE.
      ITYP=4
C
      IF (.NOT.NLPLS(ISTRAI)) GOTO 999

      IF (ISTROLD /= ISTRAI) THEN
        ISTROLD=ISTRAI
        IPLS=NSPEZ(ISTRAI)
        DO IVOLSI=1,NSRFSI(ISTRAI)
          IVL=IVOLSI
          ICC=0
          VSOURC(IVL,0)=0.D0
          IF (NLRAD) THEN
            IR1=MAX0(1,INGRDA(IVL,ISTRAI,1))
            IR2=MIN0(NR1ST,INGRDE(IVL,ISTRAI,1))
          ELSE
            IR1=1
            IR2=2
          ENDIF
          IF (NLPOL) THEN
            IP1=MAX0(1,INGRDA(IVL,ISTRAI,2))
            IP2=MIN0(NP2ND,INGRDE(IVL,ISTRAI,2))
          ELSE
            IP1=1
            IP2=2
          ENDIF
          IF (NLTOR) THEN
            IT1=MAX0(1,INGRDA(IVL,ISTRAI,3))
            IT2=MIN0(NT3RD,INGRDE(IVL,ISTRAI,3))
          ELSE
            IT1=1
            IT2=2
          ENDIF

          ISTEP=NINT(SORIND(IVL,ISTRAI))
          IFPLS=IFREC(IPLS)
          DO IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            IF (ISTEP.EQ.IRRC) THEN
              IFRC=IIRC
            ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL EXTERNAL VOLUMETRIC SOURCE PROCESSES FOR SPECIES IPLS
              IFRC=0
              ISTEP=-1
            ELSE
              CYCLE
            ENDIF
            DO IR=IR1,IR2-1
              DO IP=IP1,IP2-1
                DO IT=IT1,IT2-1
                  NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  ADD=FREC(IFPLS,IFRC,NCELL)-
     .                FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                  IF (ADD.GT.0.D0) THEN
                    ICC=ICC+1
                    ISOURC(IVL,ICC)=NCELL
                    VSOURC(IVL,ICC)=VSOURC(IVL,ICC-1)+ADD
                  ENDIF
                END DO
              END DO
            END DO
          END DO ! IIRC
          ICMX(IVL)=ICC
          VSMXI(IVL) = 1._DP / VSOURC(IVL,ICC)
        END DO ! IVOLSI
C
      END IF
C
C  FIND CELL NUMBER: NCELL
C
      IF (INDIM(NVLM,ISTRAI) .GE. 0) THEN
cdr analog sampling, no weighting
!PB  choose cell according to cell contribution to total source strength
        IC1=0
        IC2=ICMX(NVLM)
        ZEP1=RANF_EIRENE()*VSOURC(NVLM,IC2)

        IL=IC1
        IU=IC2

c  binary search
        DO WHILE (IU-IL.gt.1)
          IM=INT((IU+IL)*0.5)
          IF(ZEP1.GE.VSOURC(NVLM,IM)) THEN
            IL=IM
          ELSE
            IU=IM
          ENDIF
        END DO
c
        ICELL=IU

        NCELL=ISOURC(NVLM,ICELL)
      ELSE

cdr non-analog sampling.
cdr Here use uniform distribution of cell indices and weighting
cdr tbd: correlated sampling: use previous (reference) distribution and weighting
cdr      rather than uniform sampling.
        IC1=0
        IC2=ICMX(NVLM)

        ICELL = MIN(INT(1+RANF_EIRENE()*(IC2-1)),IC2)+IC1
        NCELL = ISOURC(NVLM,ICELL)

        WEIGHT=(VSOURC(NVLM,ICELL)-VSOURC(NVLM,ICELL-1))*VSMXI(NVLM)*
     .         (IC2-IC1)+IC1

      END IF
C
      IF (NCELL.GT.NSURF) GOTO 991
C
C  A CELL NUMBER NCELL HAS NOW BEEN SAMPLED
C
      CALL EIRENE_NCELLN(NCELL,NRCELL,NPCELL,NTCELL,NACELL,NBLOCK,
     .            NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
C
C  FIND TOROIDAL COORDINATE IN NTCELL
C
      IF (.NOT.NLTOR) THEN
C       NTCELL=1
        IPERID=1
        IF (NLTRZ) THEN
          Z0=0.
        ELSEIF (NLTRA) THEN
C  TACTICALLY ASSUME: PARTICLE STARTS IN LOCAL TOR. BASIS CELL NO.1
          ZRM1=ZSURF(1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
          IPERID=1
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
C         IPERID=LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMVOL      ')
        ELSEIF (NLTRT) THEN
          GOTO 999
        ENDIF
      ELSEIF (NLTOR) THEN
        IPERID=NTCELL
C  SAMPLE IN CELL NTCELL
        IF (NLTRZ) THEN
          Z0=ZSURF(NTCELL)+RANF_EIRENE()*(ZSURF(NTCELL+1)-ZSURF(NTCELL))
        ELSEIF (NLTRT) THEN
          PHI=ZSURF(NTCELL)+RANF_EIRENE()*
     .        (ZSURF(NTCELL+1)-ZSURF(NTCELL))
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ELSEIF (NLTRA) THEN
          ZRM1=ZFULL*(NTCELL-1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ENDIF
      ENDIF
C
C  FIND RADIAL AND POLOIDAL COORDINATE
C
      select case (LEVGEO)
      case (1)
        X0=RSURF(NRCELL)+RANF_EIRENE()*(RSURF(NRCELL+1)-RSURF(NRCELL))
        IF (NLPOL) THEN
          Y0=PSURF(NPCELL)+RANF_EIRENE()*(PSURF(NPCELL+1)-PSURF(NPCELL))
        ELSE
          Y0=YIA+RANF_EIRENE()*(YAA-YIA)
        END IF
C..........................................................................
      case (2)
        IF (NLCRC) THEN
C  POLOIDAL COORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL COORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          X0=RR*COS(WINK)
          Y0=RR*SIN(WINK)
        ELSEIF (NLELL) THEN
CDR NOT READY. STRICTLY, THETA AND R ARE CORRELATED. USE
CDR            MARGINAL AND CONDITIONAL DISTRIBUTION F1(R) AND
CDR            F2(PHI, GIVEN R)
C  POLOIDAL COORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL COORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          RRI=RSURF(NRCELL)
          RRD=RSURF(NRCELL+1)-RRI
          RRN=(RR-RRI)/RRD
C
          ELR=ELL(NRCELL)+RRN*(ELL(NRCELL+1)-ELL(NRCELL))
          EPR=EP1(NRCELL)+RRN*(EP1(NRCELL+1)-EP1(NRCELL))
          X0=RR*COS(WINK)+EPR
          Y0=RR*SIN(WINK)*ELR
        ELSEIF (NLTRI) THEN
          GOTO 999
        ENDIF
C...................................................................
      case (3)
        IF (.NOT.NLPOL) THEN
          GOTO 999
        ENDIF
        IN = NRCELL + (NPCELL-1)*NR1ST
        ZEP1=AREA(IN)*RANF_EIRENE()
        IF (ZEP1.LE.ASIMP(1,IN)) THEN
C   POINT TO BE SAMPLED WITHIN TRIANGLE 1
          X1=XPOL(NRCELL,NPCELL)
          X2=XPOL(NRCELL,NPCELL+1)
          X3=XPOL(NRCELL+1,NPCELL+1)
          Y1=YPOL(NRCELL,NPCELL)
          Y2=YPOL(NRCELL,NPCELL+1)
          Y3=YPOL(NRCELL+1,NPCELL+1)
        ELSE
C   POINT TO BE SAMPLED WITHIN TRIANGLE 2
          X1=XPOL(NRCELL+1,NPCELL)
          X2=XPOL(NRCELL,NPCELL)
          X3=XPOL(NRCELL+1,NPCELL+1)
          Y1=YPOL(NRCELL+1,NPCELL)
          Y2=YPOL(NRCELL,NPCELL)
          Y3=YPOL(NRCELL+1,NPCELL+1)
        ENDIF
        IPOLG=NPCELL
        Z1=0.
        Z2=0.
        Z3=0.
        CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,ZZ)

C...................................................................
      case (4)
        X1=XTRIAN(NECKE(1,NCELL))
        X2=XTRIAN(NECKE(2,NCELL))
        X3=XTRIAN(NECKE(3,NCELL))
        Y1=YTRIAN(NECKE(1,NCELL))
        Y2=YTRIAN(NECKE(2,NCELL))
        Y3=YTRIAN(NECKE(3,NCELL))
        Z1=0.
        Z2=0.
        Z3=0.
        CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,ZZ)
C.................................................................
      case (5)
        X1=XTETRA(NTECK(1,NCELL))
        Y1=YTETRA(NTECK(1,NCELL))
        Z1=ZTETRA(NTECK(1,NCELL))
        X2=XTETRA(NTECK(2,NCELL))
        Y2=YTETRA(NTECK(2,NCELL))
        Z2=ZTETRA(NTECK(2,NCELL))
        X3=XTETRA(NTECK(3,NCELL))
        Y3=YTETRA(NTECK(3,NCELL))
        Z3=ZTETRA(NTECK(3,NCELL))
        X4=XTETRA(NTECK(4,NCELL))
        Y4=YTETRA(NTECK(4,NCELL))
        Z4=ZTETRA(NTECK(4,NCELL))
        CALL EIRENE_FPOLYT_4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     .                       X0,Y0,Z0)
C....................................................................
      case (10)
chf added Nov. 2018
        CALL EIRENE_SAMVOL_USR(NCELL,X0,Y0,Z0)
        NRCELL=NCELL
      end select
C
      IF (NLTRA) THEN
C  FIND Z0 FROM X01,PHI IN LOCAL TOROIDAL CELL NTCELL
        X01=X0+RMTOR
        CALL EIRENE_FZRTRI(X0,Z0,NTCELL,X01,PHI,NTCELL)
      ENDIF
C
      MTSURF=0
      NLSRFZ=.FALSE.
      MPSURF=0
      NLSRFY=.FALSE.
      MRSURF=0
      NLSRFX=.FALSE.
C
C  NEXT: ANALOG SPECIES INDEX DISTRIBUTION: WEISPZ(IPL)
C
      DO 630 JSPZ=1,NSPZ
        WEISPZ(JSPZ)=-1.
  630 CONTINUE
C
C  NOT IN USE ANYMORE
C  CURRENTLY: ONLY SINGLE SPECIES VOLUME SOURCES POSSIBLE
C  MULTI SPECIES VOL-SOURCES HAVE TO BE TREATED BY STRATIFIED SAMPLING
C     IF (NSPEZ(ISTRAI).LE.0) THEN
C       IF (NCELL.EQ.1) THEN
C         DO 640 IPL=1,NPLSI
C           IREC=0
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,1))/
C    .                  (FREC(0,  0,1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
C 640     CONTINUE
C       ELSE
C         DO 645 IPL=1,NPLSI
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,NCELL)-FREC(IFPLS,0,NCELL-1))/
C     .                 (FREC(0,    0,NCELL)-FREC(0,    0,NCELL-1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
C 645     CONTINUE
C       ENDIF
C     ENDIF
C
      CRTX=SORAD4(NVLM,ISTRAI)
      CRTY=SORAD5(NVLM,ISTRAI)
      CRTZ=SORAD6(NVLM,ISTRAI)
      CNORM=SQRT(CRTX**2+CRTY**2+CRTZ**2)+EPS60
      CRTX=CRTX/CNORM
      CRTY=CRTY/CNORM
      CRTZ=CRTZ/CNORM
!PB
      TEWL=TEIN(NCELL)
      TIWL(1:NPLSI)=TIIN(MPLSTI(1:NPLSI),NCELL)
      DIWL(1:NPLSI)=DIIN(1:NPLSI,NCELL)
      VXWL(1:NPLSI)=VXIN(MPLSV(1:NPLSI),NCELL)
      VYWL(1:NPLSI)=VYIN(MPLSV(1:NPLSI),NCELL)
      VZWL(1:NPLSI)=VZIN(MPLSV(1:NPLSI),NCELL)
      IF (maxval(ZIIN(1:NPLSI,NCELL)).GT.ZVAC) THEN
        ZIWL(1:NPLSI)=ZIIN(1:NPLSI,NCELL)
      ELSE
        ZIWL(1:NPLSI)=DBLE(NCHRGP(1:NPLSI))
      END IF
      EFWL(1:NPLSI)=0._DP
      SHWL=0._DP

      RETURN
C
C
  991 CONTINUE
      WRITE (iunout,*) 'SAMPLING ERROR IN SAMVOL'
      WRITE (iunout,*) 'NCELL,NSURF,NSBOX ',NCELL,NSURF,NSBOX
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*) 'UNWRITTEN OPTION IN SAMVOL'
      CALL EIRENE_EXIT_OWN(1)


      END SUBROUTINE EIRENE_SAMVL1

C     the following SUBROUTINE is for reinitialization of EIRENE (DMH)

      SUBROUTINE EIRENE_SAMVOL_REINIT
      IMPLICIT NONE

#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP PARALLEL DEFAULT(none)
!$OMP& SHARED(LEVGEO)
#endif
      ISTROLD = -1

      DEALLOCATE (LPLSSR)
      DEALLOCATE (LATSSR)
      DEALLOCATE (FREC)
      DEALLOCATE (VSOURC)
      DEALLOCATE (VSMXI)
      IF (LEVGEO == 2) DEALLOCATE (RQ21)
      IF (LEVGEO == 2) DEALLOCATE (PS21)
      IF (LEVGEO == 3) DEALLOCATE (ASIMP)
      DEALLOCATE (ISOURC)
      DEALLOCATE (ICMX)
      DEALLOCATE (IFREC)
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP END PARALLEL
#endif

      return

      END SUBROUTINE EIRENE_SAMVOL_REINIT

      END MODULE EIRMOD_SAMVOL
