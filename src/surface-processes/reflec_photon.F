CDR  NOV 17 : lemtspw arguments corrected
cdr  jan 18 : start to implement bi-directional reflectance functions
cdr  dec.18 : remove unfinished "hollmann database model"
cpb  oct.19 : set ISPZO=ISPZ for default reflection model
cdr  oct 21 : ilref=1 model: partial absorbing, cosine, specular
C
C
C  REFLECT ESCAPING PHOTONS
C  INPUT:
C       ILREF = 0  PERFECT ABSORPTION
C       ILREF = 1  PARTIAL REFLECTION, SPECULAR PLUS COSINE COMPONENTS
C       ILREF = 2  BI-DIRECTIONAL REFLECTANCE FUNCTION: tbd.  currently: absorb.
C
C       ITYP  = 0  INCIDENT PHOTON
C  OUTPUT:
C     LGPART= TRUE AND:
C       ITYP = 0  PHOTON IPHOT IS RETURNED TO CALLING PROGRAM
C     LGPART= FALSE  NO PARTICLE IS RETURNED (ABSORPTION)
C       ITYP = -1
C
C
C  INITIALIZE SURFACE REFLECTION MODELS FOR PHOTONS
C
      SUBROUTINE EIRENE_REFLC0_PHOTON
      USE EIRMOD_CTRCEI
      IMPLICIT NONE
C
      INTEGER, SAVE :: IFIRST=0
      EXTERNAL :: EIRENE_LEER

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(IFIRST)
#endif

      IF (IFIRST.EQ.1) RETURN
      IFIRST=1
C
C  INITIALIZE HOLLMAN REFLECTANCE FUNCTION; OUT
C     formerly: CALL EIRENE_INIT_REFL_HLM
C
      IF (TRCREF) THEN
        CALL EIRENE_LEER(2)
      ENDIF
C
C
C  PRINTOUT REFLECTION PROPERTIES OF SURFACES
C  ALREADY DONE FROM REFLC0
C
C
      RETURN
      END SUBROUTINE EIRENE_REFLC0_PHOTON
C
      SUBROUTINE EIRENE_REFLC1_PHOTON (WMIN,XMP,XCP,NPRIN,IGASF,IGAST)
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CTRCEI
      USE EIRMOD_CADGEO
      USE EIRMOD_CLGIN
      USE EIRMOD_COMUSR
      USE EIRMOD_CLOGAU
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_CRAND
      USE EIRMOD_CREF
      USE EIRMOD_CCONA
      USE EIRMOD_PHOTON
      USE EIRMOD_CPES
      USE EIRMOD_CSDVI
      USE EIRMOD_RANF, ONLY: RANF_EIRENE, RANSET_EIRENE, RANGET_EIRENE
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: WMIN, XMP, XCP
      INTEGER, INTENT(IN) :: NPRIN
      INTEGER, INTENT(INOUT) :: IGASF, IGAST
      INTEGER, SAVE :: NPANOLD=0
      INTEGER :: MODREF, ISPZO,
     .           MSS, IDUMMY
      REAL(DP) :: DUMMY, XMW, XCW, COSI2,
     .            EXPP, EXPE, EXPI, RINTG, EINTG, AINTG, COSIN,
     .            THETA_OUT=0._DP, ALPHA_OUT=0._DP,
     .            FR1, ZCPHI, ZSPHI, ZCTHET, ZSTHET, VX, VY, VZ,
     .            RPROB, PRFCF, PRFCT, PABS, PLAMBERT, WABS,
     .            ZEP1
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      EXTERNAL :: EIRENE_ROTATE, EIRENE_ROTATF

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(NPANOLD)
#endif
C
C  RE-SYNCHRONIZE RANDOM NUMBERS (CORRELATED SAMPLING): out.
C  tbd: introduce counter of max. possible random numbers used up to this point
C  see report Kalos, 1955, on correlated sampling, synchronisation...
C
      IF (NLCRR.AND.(NPANU.NE.NPANOLD).AND..FALSE.) THEN
C  RE-INITIALIZE RANDOM NUMBERS FOR EACH PARTICLE, TO GENERATE CORRELATION

        idummy=ranset_eirene(iseedr)
        DUMMY=RANF_EIRENE()
        ISEEDR=ranget_eirene(iseedr)
        ISEEDR=INTMAX-ISEEDR
        NPANOLD=NPANU
      END IF
C..................................................................
C
C  SURFACE NUMBER  : MSURF (MSURF=0: DEFAULT MODEL)
C  SPECIES INDEX   : ISPZ
C
cdr fix ph2
      ISPZO=ISPZ
      MODREF=ILREF(ISPZ,MSURF)
      XMW=ZNML(MSURF)
      XCW=ZNCL(MSURF)
C
C
C   MODREF=0: "PERFECTLY ABSORBING SURFACE", DEFAULT
C   MODREF=1: "PARTIALLY ABSORBING, SPECULAR OR COSINE"
C

!    RPROB: probability of the photon to be reflected
      RPROB = 0._dp

      IF (MODREF.EQ.1.or.MODREF.EQ.2) THEN
CDR MODREF=1 simple specular, or cosine or absorb:  recycf--> recyct --> 1.0
CDR modreF=2 BDRF MODEL, WITH GAUSSIAN LOBE IN SPECULAR PART AND LAMBERTIAN IN THERMAL PART

cdr reflected (specular) fraction
        PRFCF=max(0._DP,(min(1._DP,RECYCF(ISPZ,MSURF))))
cdr prfct: total reemitted fraction.
cdr pabs : (1-recyct) is the absorbed fraction
        PRFCT=max(0._DP,(min(1._DP,RECYCT(ISPZ,MSURF))))
        PABS=1.-PRFCT
cdr prfct is cumulated SUM of emitted (Lambertian) and specular reflected (BDRF) part
cdr cumulative distribution of (spec-ref)-(lambert)-(absorb) fractions: prfcf,prfct,1.0
        prfcf=min(prfcf,prfct)
        plambert=prfct-prfcf

cdr parameters for bi-directional reflectance function
        EXPP=EXPPL(ISPZ,MSURF)
        EXPE=EXPEL(ISPZ,MSURF)
        EXPI=EXPIL(ISPZ,MSURF)

        RINTG=RINTEG(MSURF)
        EINTG=EINTEG(MSURF)
        AINTG=AINTEG(MSURF)
        ISPZO=ISPZ
      ENDIF
C
C
      IF (MODREF.EQ.1) THEN
cdr: tbd.: suppression of absorption. Coding: See below
C   TENTATIVELY ASSUME REFLECTION
        LGPART=.TRUE.

C   CHECK FOR ABSORPTION, LAMBERTIAN OR SPECULAR COMPONENT
        ZEP1=ranf_eirene()
        IF (ZEP1.GT.PRFCT) THEN  ! (PRFCT < ZEP1 <= 1)
C   ABSORPTION
          GOTO 700
        ELSEIF (ZEP1.GT.PRFCF) THEN
C   LAMBERTIAN (COSINE) DISTRIBUTION
          GOTO 600
        ELSE                     ! (0 < ZEP1 <= PRFCF)
C   COSINE OF ANGLE OF INCIDENCE against outer normal
          COSIN=VELX*CRTX+VELY*CRTY+VELZ*CRTZ
          IF (COSIN.LT.0.D0) GOTO 993
C   PURELY SPECULAR REFLECTION
          COSI2=-(COSIN+COSIN)
          VELX=VELX+COSI2*CRTX
          VELY=VELY+COSI2*CRTY
          VELZ=VELZ+COSI2*CRTZ
          RETURN
        ENDIF
      ENDIF

      IF (MODREF.GE.2) THEN
c  SPECULAR LOBE, BDRF MODEL. UNFINISHED
c  because ilref=2 is the default, for all test particles, here
c  we use absorption as default for photons.
        GOTO 700
      ENDIF
C
C  UNFINISHED DATABASE MODEL REMOVED HERE. OLD STATEMENT LABELS: 100 TO 130.
C
cdr  currently this code part from here to statement label 600 cannot be reached.
!    RPROB: probability of the photon to be reflected
      RPROB = 0._dp

C  DECIDE IF PARTICLE IS TO BE REFLECTED OR ABSORBED
C  (NO THERMAL RE-EMISSION MODEL FOR INCIDENT PHOTONS)
C
      IF (WEIGHT.GT.WMIN) THEN
C  WITH SUPPRESSION OF ABSORPTION
        WABS=WEIGHT*(1.D0-RPROB)
        IF (WABS.GT.0.D0) THEN
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WABS
            IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              SPUMP(ISPZO,MSURFG)=SPUMP(ISPZO,MSURFG)+WABS
            END IF
            LMETSPW(ISPZO) = .TRUE.
          ENDIF
        ENDIF
        WEIGHT=WEIGHT-WABS
        IF (WEIGHT.LE.EPS30) GOTO 700
      ELSE
C  NO SUPPRESSION OF ABSORPTION
        FR1=RANF_EIRENE( )
        IF (FR1.GT.RPROB) GOTO 700
      ENDIF
C
C  SPECIES OF REFLECTED PARTICLE
      IF (IGASF.LT.1.OR.IGASF.GT.NPHOTI) GOTO 992
      IPHOT=IGASF
      ISPZ=IPHOT
      ITYP=0
C
C  ENERGY (WAVELENGTH): NOT MODIFIED
C
  600 CONTINUE
C
C   LAMBERTIAN (COSINE DISTRIBUTION)
C
      IF (INIV4.LE.0) CALL EIRENE_FCOSIN
      VX=FC1(INIV4)
      VY=FC2(INIV4)
      VZ=FC3(INIV4)
      INIV4=INIV4-1
c  /vx,vy,vz/ is cosine distributed around surface normal vector/-1,0,0/
      CALL EIRENE_ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
      RETURN
C
Cdr  unfinished option.... rc: which one?
C  POLAR ANGLE OF REFLECTION (THETA)
      ZCTHET=COS(THETA_OUT)
C  LIMIT COSINE OF POLAR ANGLE TO 85. DEGREES
C  (I.E., 5 DEGREES AGAINST SURFACE TANGENTIAL PLANE)
      ZCTHET=MIN(0.999999_DP,MAX(0.08716_DP,ZCTHET))
      ZSTHET=SQRT(1.-ZCTHET*ZCTHET)
C
C  AZIMUTHAL ANGLE OF REFLECTION (PHI)
C
      ZCPHI=COS(ALPHA_OUT)
      ZCPHI=MAX(-.999999_DP,MIN(0.999999_DP,ZCPHI))
      ZSPHI=SQRT(1.-ZCPHI*ZCPHI)
      ZSPHI=ZSPHI*SIGN(1._DP,(RANF_EIRENE( )-0.5_DP))
C
      VX=-ZCTHET
      VY=ZSTHET*ZSPHI
      VZ=ZSTHET*ZCPHI
      IF (COSIN.GT.0.999999) THEN
        CALL EIRENE_ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
      ELSE
        CALL EIRENE_ROTATE
     .  (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ,COSIN)
      ENDIF
      RETURN
C
C  ABSORB PARTICLE AT THIS SURFACE
C
  700 CONTINUE
      IF (LSPUMP) THEN
        IF (MSURF.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZO,MSURF)=SPUMP(ISPZO,MSURF)+WEIGHT
        ENDIF
        IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZO,MSURFG)=SPUMP(ISPZO,MSURFG)+WEIGHT
        END IF
        LMETSPW(ISPZO) = .TRUE.
      END IF
      LGPART=.FALSE.
      WEIGHT=0.
      ITYP=-1
      RETURN
C
C  ERROR MESSAGES FROM SUBR. REFLEC_PHOTON
C
C
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (iunout,*) 'MSURF = ',MSS
      WRITE (iunout,*) 'IGASF,IGAST,MODREF ?? '
      WRITE (iunout,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
c
  993 CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. REFLEC_PHOTON '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (iunout,*) 'MSURF = ',MSS
      WRITE (iunout,*) 'COSIN.LT.0. ', COSIN
      WRITE (iunout,*) 'STOP HISTORY NO. NPANU= ',NPANU
      GOTO 999
C
  999 IF (NLTRC)  CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
      LGPART=.FALSE.
      WEIGHT=0.
      RETURN
C
      END SUBROUTINE EIRENE_REFLC1_PHOTON
