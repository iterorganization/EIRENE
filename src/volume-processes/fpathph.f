c  25.11.05: option modcol(3,4...)=3 added
c            (first implemented in fpatha)
c            CX rate option 4 added (adopted from fpatha)
C               added: jcou,ncou
!pb  30.08.06:  data structure for reaction data redefined
!pb  12.10.06:  modcol revised
!pb  22.11.06:  flag for shift of first parameter to rate_coeff introduced
!pb  28.11.06:  initialization of XSTOR reactivated because of trouble in
!pb             BGK iteration
!pb  22.03.07:  PI reactions revised
cdr  oct.14  :  ftabcx3 added. Full tests still to be done
cdr  oct.14  :  synchronized with fpathm, fpathi

cdr 31.10.14 :  speedup of final cut-off evaluations

cdr note:       sgnl_poly evaluations are just the 8th-order polynomial,
cdr             plus rcmin,rcmax consideration.
cdr             unless rcmin,rcmax are set (as it is the case currently here),
cdr             there is no need to call  --> move to in-line
cdr 06.08.15 :  arguments added to vecusr
cdr 13.08.15 :  cflag(4,1) changed from 2 to 1 (as it was in fpatha).  Is that correct ??

cdr dec. 15:    missing: ftabel3
cdr jan. 16:    call to ftabcx3 added and tested for modcol=1 option

cdr aug. 16:    bug fix re EXPO in PI branch
cdr sept.16:    PI process: use v0/vth >> 1. to switch to beam-rate coeff
cdr             EI process: started to check for H.3, H.1 options for EI processes
cdr                         according to v0/vth >> 1. criteria
cdr Nov. 16:    cflag(7,mstor0) rather than cflag(6,3), see comments

cdr Jan. 18:    This entire routine is largely unfinished. Photon transport
cdr             with eirene currently not possible.
cdr             Started to prepare re-activating this option: for now: comments only
cdr Nov. 18:    notational cleanup: OT processes --> PH processes, to avoid confusion
cdr             with OT ("other") processes of type H.11, H.12, population ratios.
cdr

C
      FUNCTION EIRENE_FPATHPH (K,CFLAG,JCOU,NCOU)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR PHOTON
C   "BEAM PARTICLES" , SPECIES IPHOT, OF VELOCITY (E0,VEL_X,Y,Z) IN DRIFTING MAXWELLIAN BACKGROUND MEDIUM
C   IN CELL K

C
C   INPUT:
C   IPHOT     : PHOTON LINE SPECIES INDEX (INPUT VIA COMMON)
C   K         : CURRENT GRID CELL
C   JCOU, NCOU: THERE WILL BE NCOU CALLS TO FPATH, FOR SAME TEST PARTICLE
C               COORDINATES WITH DIFFERENT CELL NUMBER K.
C               THIS CURRENT CALL IS CALL NO. JCOU.

C   OUTPUT: COMMON COMLCA
C           CFLAG: FLAG FOR SAMPLING OF POST-COLLISION STATES
C           CFLAG(1,...): EI
C           CFLAG(2,...): NOT IN USE, was DS process class in very old versions
C           CFLAG(3,...): CX
C           CFLAG(4,...): PI
C           CFLAG(5,...): EL
C           CFLAG(6,...): RC
c           CFLAG(7,...): PH (photonic processes, formerly: OT)
C
C   FLAG FOR POST-COLLISION DISTRIBUTION IN VELOCITY SPACE
C  CFLAG(...,IRCL),  IRCL: IREI,..., IRCX,IRPI,IREL,IRRC,IRPH
C      =0: VI: DELTA COLLISION IN VELOCITY SPACE (BUT DIFFERENT
C                                                 SPECIES ALLOWED)
C      =1: VI: MONOENERGETIC AND ISOTROPIC IN FRAME MOVING WITH BULK SPECIES
C      =2: VI: DRIFTING MAXWELLIAN
C      =3: VI: SIGMA-V-WEIGHTED MAXWELLIAN IN FRAME MOVING WITH BULK SPECIES
C      =X  VI: DELTA COLLISION IN VELOCITY SPACE: VI=V0 (BUT DIFFERENT SPECIES ALLOWED)
C              TO BE WRITTEN
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CZT1
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_PHOTON

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(7,MSTOR0)
      INTEGER, INTENT(IN) :: K,JCOU,NCOU

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLSV)
      REAL(DP) :: EIRENE_FPATHPH,
     .            sigmax, sigv,
     .            DENEL, VX, VY, VZ, PVELQ0, fac,
     .            XC, YC, ZC,
cdr  functions for 'on the fly' evaluation of a&m data
     .            EIRENE_FEPLPH3
      INTEGER :: J, KK, irph, ipph, IL, jpls,
     .           IPLSV
      EXTERNAL :: EIRENE_VECUSR, EIRENE_EXIT_OWN,
     .            EIRENE_FEPLPH3
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTORV=0.D0
!pb   IF (NCOU.GT.1) THEN
        XSTOR=0.D0
!pb   ENDIF
      EIRENE_FPATHPH=1.D10
      SIGMAX=0.D0
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)

      DO 2 JPLS=1,NPLSI
        ZTI(JPLS)=ZT1(JPLS,K)
        DENIO(JPLS)=DIIN(JPLS,K)
    2 CONTINUE
C
C  TRANSFORM TEST PARTICLE VELOCITY TO FRAME MOVING WITH BULK SPECIES IPLS
C            PVELQ(IPLSV) IS THE VELOCITY IN THESE REFERENCE FRAMES, SQUARED 
C
      PVELQ0=VEL*VEL
      DO 3 JPLS=1,NPLS
        IPLSV=MPLSV(JPLS)
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            XC=0.
            YC=0.
            ZC=0.
            CALL EIRENE_VECUSR (2,K,XC,YC,ZC,VX,VY,VZ,JPLS,.FALSE.)
          ELSE
            VX=VXIN(IPLSV,K)
            VY=VYIN(IPLSV,K)
            VZ=VZIN(IPLSV,K)
          END IF
          PVELQ(IPLSV)=(VELX*VEL-VX)**2+
     .                 (VELY*VEL-VY)**2+
     .                 (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLSV)=PVELQ0
        ENDIF
    3 CONTINUE
C
c  PH processes (photonic reactions)
c
csw
      if(LGPHPH(iphot,0,0) == 0) goto 70
      do 61 ipph=1,NPHPHI(iphot)
!  -->  lgxph, with x=ph, IRPH corresponds to: irei, ircx, ....
        irph =LGPHPH(iphot,ipph,0)
!  -->  ipls: bulk, with iteration, just like others
        ipls =LGPHPH(iphot,ipph,1)
!  -->  kk=NREAPH(IRPH), analogous to KK=NREAPI(IRPI) for PI processes
        kk   =NREAPH(IRPH)

        IF (LGVAC(K,IPLS)) GOTO 61
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(7,2,IRPH).EQ.1) THEN
          GOTO 997
        ELSEIF (MODCOL(7,2,IRPH).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE. FULL ACCOUNT FOR DOPPLER SHIFT

cdr       kk   = nreaph(IRPH)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr  getcoeff liefert nun maxw. average ueber Ti(ipls) (background neutrals), z.b. voigt, ....
crc check if diin > 0
          IF (diin(ipls,k) > 0.D0) THEN
              call EIRENE_PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
              sigv=sigv*diin(ipls,k)
          ELSE
              sigv=0.D0
          ENDIF
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVPH(IRPH)=SIGV
          GOTO 997
        ELSEIF (MODCOL(7,2,   IRPH).EQ.4) THEN
C  MODEL 4:
C  BEAM-BEAM RATE. IGNORE DOPPLER SHIFT DUE TO THERMAL MOTION,
C                  INCLUDE DOPPLER SHIFT DUE TO DIRECTED MOTION
cdr       kk   = nreaph(IRPH)
cdr   effective energy e0_eff due to doppler shift from directed motion
cdr       e0_eff=
cdr       ireac=modcol(7,1,IRPH)
cdr  ireac entspricht "typ" in Getcoeff - cross-section (lorentz, vdw, ...)
cdr  allerdings kann hier der "querschnitt" von hintergrundparametern abhaengen
crc check if diin > 0
          IF (diin(ipls,k) > 0.D0) THEN
              call EIRENE_PH_GETCOEFF(kk,iphot,0,k,ipls,fac,sigv)
              sigv=sigv*diin(ipls,k)
          ELSE
              sigv=0.D0
          ENDIF
cdr  besser: sig vc. E0_effective, ohne v = vrel = c, dann
cdr  dann:   sigv=sig * c * diin
cdr  denn:   sigv enthaelt hier keine faltung ueber background maxw. vel-verteil
          if(phv_muldens .EQ. 0) then
cdr  hier in fpathph kann es keine spontanen raten (1/s) geben.
cdr  spaeter: allgemein raten (1/s) auch fuer testteilchen (fpatha, fpathm, fpat
cdr           als neue option einfuehren, analog Aik in xsectp.
            GOTO 997
          endif
          SIGVPH(IRPH)=sigv
cdr       ESIGPH(IRPH,1)=e0*sigv   ziemlich sicher falsch
        ELSE
          GOTO 997
        ENDIF

        SIGMAX=MAX(SIGMAX,SIGVPH(IRPH))
        SIGPHT=SIGPHT+SIGVPH(IRPH)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(7,4,IRPH).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
C  (ONLY NEEDED FOR TRACKLENGTH ESTIMATOR)
          IF (NSTORDR >= NRAD) THEN
            ESIGPH(IRPH,1)=EPLPH3(IRPH,K,1)
          ELSE
            ESIGPH(IRPH,1)=EIRENE_FEPLPH3(IRPH,K)
          END IF
          CFLAG(7,IRPH)=2
        ELSEIF (MODCOL(7,4,IRPH).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGPH(IRPH,1)=EPLPH3(IRPH,K,1)
          ELSE
            ESIGPH(IRPH,1)=EIRENE_FEPLPH3(IRPH,K)
          END IF
          CFLAG(7,IRPH)=1
        ELSE
          GOTO 997
        ENDIF
   61 CONTINUE

   70 CONTINUE
c
C     TOTAL
C
C
C  CUT OFF RESIDUAL RATES, WHICH SHOULD STRICTLY BE ZERO
C  TO AVOID SPURIOUS ENTRIES TO COLLISION RATE TALLIES
C  CURRENTLY: CUT-OFF AT 1E-10 TIMES SIGMAX
C
      IF (SIGPHT.GT.0._DP) THEN
        DO IRPH=1,NRPH
          IF (SIGVPH(IRPH) .LE. SIGMAX*1.D-10) THEN
            SIGPHT=SIGPHT-SIGVPH(IRPH)
            SIGVPH(IRPH) = 0.D0
          END IF
        END DO
      END IF

      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT+SIGPHT
      IF (SIGTOT.GT.1.D-20) THEN
        EIRENE_FPATHPH=VEL/SIGTOT
        ZMFPI=1./EIRENE_FPATHPH
      ENDIF
C
      RETURN
  997 CONTINUE
      WRITE (iunout,*)
     .  'ERROR IN FPATHPH: INCONSISTENT PHOTON COLL. DATA'
      WRITE (iunout,*) 'ITYP,IPHOT,IRPH,MODCOL(7,J,IRPH),J=1,4 '
      WRITE (iunout,*)  ITYP,IPHOT,IRPH,(MODCOL(7,J,IRPH),J=1,4)
      CALL EIRENE_EXIT_OWN(1)
      END FUNCTION EIRENE_FPATHPH
