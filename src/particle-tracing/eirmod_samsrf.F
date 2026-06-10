      MODULE EIRMOD_SAMSRF
C
C  THIS SUBROUTINE SAMPLES THE INITIAL COORDINATES OF A PARTICLE
C  HISTORY, WHICH STARTS ON A SURFACE
C
c  eirene_samsf0:
c  eirene_samsf1:
c  eirene_samsf2:  deallocate temporary arrays

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CINIT
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_CLGIN
      USE EIRMOD_CTRIG
      USE EIRMOD_LEARC1, ONLY: EIRENE_LEARC1
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_SAMUSR, ONLY: EIRENE_SAMUSR_INIT, EIRENE_SAMUSR
      USE EIRMOD_SHEATH, ONLY: EIRENE_SHEATH
      USE EIRMOD_COMXS, ONLY: L_MACRO

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_SAMSF0, EIRENE_SAMSF1,
     .          EIRENE_SAMSF2

      REAL(DP), ALLOCATABLE, SAVE ::
     .        ALEFT(:,:,:), BRGHT(:,:,:), XI(:,:,:), XE(:,:,:)
      REAL(DP), SAVE :: FL,VX,VY,VZ,XC,YC,ZC
      INTEGER, ALLOCATABLE, SAVE :: INDTEC(:,:)
      INTEGER, SAVE :: ISTEP_SPEZ, ISTEP, IS1, IPLSTI, IPLSV

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(FL,VX,VY,VZ,XC,YC,ZC,
!$OMP& ISTEP_SPEZ,ISTEP,IS1,IPLSTI,IPLSV)
#endif

      CONTAINS

C  6.4.04: include torl in default step function
C  nov. 11.05  jet-2005 patch 1: new parameter shwl at call to smsrf1
C              added: print ekflx, eshflx: total kinetic and sheath fluxes
C              from newly added step functions  elstep, shstep
C              now return shwl, efwl(ipls) for use in locate.f
C              in locate: new options nemod1=8,9 for using efwl(ipls)
C              efwl(ipls) is the kinetic ion energy flux (W) at sheath entrance
C  nov.14.05: bug fix, v.k., istra --> istep for aleft(1,....) (twice)
C  may 06: sampling on 2nd-order additional surface modified:
C          sign of determinant is now sampled. rejection with respect to bounding box
c  aug 06: new option SORIND.gt.100 for step functions introducted. 3rd digit is now
c          pointer to species of step function
!pb 22.03.07: LEVGEO=6 --> LEVGEO=10
C  june 13: some comments, clean-up, further checks and error exits.
C           DELR REMOVED FROM LOOP 3010
C           ELSTEP and SHSTEP added in definition of default step function
C  oct. 14: some preliminary options for correlation sampling removed,
C           back to development branch
c  may  15: argument in first call to vecusr: ipls, rather than iplsv, now everywhere
c  aug. 15: arguments in vecusr added: ncell, x0, y0, z0
C  NOV. 15: INDSRF: SURFACE NUMBER FOR SHEATH MODEL, ONLY IN CASE OF STEP FUNCTION OPTION
CDR         now: default is ALWAYS set. INDSRF is e.g. argument in call to fct. SHEATH(...)
cdr nov.16: istra --> istrai, ispz -->jspz, and a bit more info on diagnostic printout
cdr dec.17: cleanup, comments
cdr may 20: Default step function: projection to flux tube in levgeo=4 (non-orthodonal
cdr         grids) was turned off at some point in time. This must be unphysical.
cdr jun.22: step function species range was the IPLS range, for the
cdr         default plasma recycling step functions set here.
cdr         This may (and did) cause
cdr         wrong physics and even code crashes when used with atm, mol or test ion
cdr         sources with step function distributions, e.g. gas puffs.
c   sep 22: Addressed the issue above for step functions related to atoms and molecules
cnh nov 23: The solution of the issue above gave problems for the hybrid fluid-kinetic model
c           where particles are transferred from the fluid atom to the kinetic atom
c           population and sent into the void regions. Hence, the fix is not done for an
c           atomic source in combination with the hybrid model (L_macro > 0). Hence, issues
c           will still be present for some multi-species gas puffs with the hybrid model.
C
      SUBROUTINE EIRENE_SAMSF0

C  INITIALIZE DATA FOR SURFACE SAMPLING FOR STRATUM NO. ISTRAI
C
      IMPLICIT NONE

cym variable defined in the module removed from the list
      REAL(DP):: TORL(NSTEP,NGITT), GAMMA, CUR,
     .           RANDIF,
     .           DELR, TESH, CTHETA, DELRR, CS, V_MEAN, ZZ,
     .           EIRENE_STEP0, EIRENE_STEP
cym end
      REAL(DP):: FLX(NPLS),EKFLX(NPLS),ESHFLX(NPLS),
     .           DISH(NPLS),VPSH(NPLS),ZISH(NPLS)
      INTEGER :: ISPZD(NPLS), ISP, ISPZ
      INTEGER :: ISTRAI, IERROR, ISRFS, ISOR, ISORFL, INDSRF, ISTR, ISR,
     .           NL3J, NL2J, NL1J, IP, ISTS, IT, KAN, KEN,
     .           JATM, JMOL, JPLS, JSPZ,
     .           K, NBIN, NSMX, NANZ, IS, ITRI, ITET,
     .           ISGRD1, IS2, ISGRD2,
     .           ISGRD3, INS,
     .           MSTEP(NSTEP),MMSTEP(NSTEP)
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      EXTERNAL :: EIRENE_MASPRM, EIRENE_MASJ2R, EIRENE_MASJR2,
     .            EIRENE_MASR2, EIRENE_LEER, EIRENE_EXIT_OWN,
     .            EIRENE_STEP, EIRENE_STEP0, EIRENE_TET_STEP,
     .            EIRENE_VECUSR

      MSTEP=0      ! REDUCE DUPLICATED DIAGNOSTIC PRINTOUT
                   ! FOR STEP FUNCTION ISTEP
      MMSTEP=ISTUF ! STEP FUNCTIONS ALREADY SET EXTERNALLY

      DO ISTRAI=1,NSTRAI

      IF (.NOT.NLSRF(ISTRAI)) CYCLE

      IERROR=0

      IF (.NOT.ALLOCATED(INDTEC)) THEN
#ifndef LEGACYCOMP
        ALLOCATE (INDTEC(3*NSRFS,NSTRA), SOURCE = 0)
        ALLOCATE (ALEFT(3,NSRFS,NSTRA), SOURCE = 0._DP)
        ALLOCATE (BRGHT(3,NSRFS,NSTRA), SOURCE = 0._DP)
        ALLOCATE (XI(3,NSRFS,NSTRA), SOURCE = 0._DP)
        ALLOCATE (XE(3,NSRFS,NSTRA), SOURCE = 0._DP)
#else
        ALLOCATE (INDTEC(3*NSRFS,NSTRA))
        ALLOCATE (ALEFT(3,NSRFS,NSTRA))
        ALLOCATE (BRGHT(3,NSRFS,NSTRA))
        ALLOCATE (XI(3,NSRFS,NSTRA))
        ALLOCATE (XE(3,NSRFS,NSTRA))
        INDTEC = 0
        ALEFT = 0._DP
        BRGHT = 0._DP
        XI = 0._DP
        XE = 0._DP
#endif
      END IF
C
C  LOOP OVER SOURCE SURFACES: ISRFS
C
      DO 1 ISRFS=1,NSRFSI(ISTRAI)
C
c  sampling distribution, for all three coordinates (and time):
c                         4 digits: TZYX
        ISOR=NINT(SORLIM(ISRFS,ISTRAI))
c  initial birth point flags (ifpath,....) for particle tracing
        ISORFL=EIRENE_IDEZ(NINT(SORIFL(ISRFS,ISTRAI)),4,4)
c  number of surface for current surface source segment
        INDSRF=INSOR(ISRFS,ISTRAI)
        IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)

C  HAS THIS SURFACE SOURCE A PRE-PROGRAMMED DISTRIBUTION
C
        IF (ISOR.LE.0) THEN
C  NO. USER-DEFINED SOURCE SAMPLING. NOW INITIALIZE USER-SUPPLIED SOURCE SAMPLING
          ISR=ISRFS
          ISTR=ISTRAI
          CALL EIRENE_SAMUSR_INIT(ISR,ISTR,
     .                SORAD1(ISR,ISTR),SORAD2(ISR,ISTR),
     .                SORAD3(ISR,ISTR),SORAD4(ISR,ISTR),
     .                SORAD5(ISR,ISTR),SORAD6(ISR,ISTR))
          GOTO 1
        ENDIF
C
        IF (INDIM(ISRFS,ISTRAI) .EQ. 1) THEN

c  source is on radial (x-) grid surface x= const. r= const, etc...
c  the poloidal range of source region should be on ingrd..(...,2),
c  not on ingrd..(...,1)
c  sample 2nd and 3rd coordinate, evaluate 1st coordinate
          IF (INGRDA(ISRFS,ISTRAI,1).NE.INGRDE(ISRFS,ISTRAI,1)) THEN
            WRITE (iunout,*) 'WARNING FROM SAMSF0, ISTRAI= ',ISTRAI
            WRITE (iunout,*) 'NEW INPUT FOR INGRDA,INGRDE....'
            WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT'
            INGRDA(ISRFS,ISTRAI,2)=INGRDA(ISRFS,ISTRAI,1)
            INGRDE(ISRFS,ISTRAI,2)=INGRDE(ISRFS,ISTRAI,1)
            INGRDA(ISRFS,ISTRAI,1)=INSOR(ISRFS,ISTRAI)
            INGRDE(ISRFS,ISTRAI,1)=INSOR(ISRFS,ISTRAI)
          ENDIF
        ENDIF

C  set sampling distributions for 2 of the three coordinates x,y,z
c  2 coordinates are sampled,
c  the third coordinate is evaluated from surface equation
        NL1J=ISRFS
        NL2J=NL1J+NSRFS
        NL3J=NL2J+NSRFS
c  FLAG FOR SPATIAL SAMPLING IS CODED ON INDTEC
c  (ORIGINALLY: THIS WAS INPUT FLAG SORLIM --> ISOR)
        INDTEC(NL1J,ISTRAI)=EIRENE_IDEZ(ISOR,1,4)
        INDTEC(NL2J,ISTRAI)=EIRENE_IDEZ(ISOR,2,4)
        INDTEC(NL3J,ISTRAI)=EIRENE_IDEZ(ISOR,3,4)
C
C  IS A STEP FUNCTION REQUESTED? (indtec=4 for one of the coordinates?)
C
        ISTEP=0
        ISTEP_SPEZ=0

        IF (INDTEC(NL1J,ISTRAI).NE.4.AND.INDTEC(NL2J,ISTRAI).NE.4.AND.
     .      INDTEC(NL3J,ISTRAI).NE.4) GOTO 7
C
C  YES. CHECK INPUT DATA AND STORAGE
C
C  ISTEP      IS 1ST AND 2ND DIGIT "BA" OF REAL FLAG SORIND (="CBA.0")
C  ISTEP_SPEZ IS 3RD DIGIT "C" OF REAL FLAG SORIND
        ISTEP=MOD(NINT(REAL(SORIND(ISRFS,ISTRAI),DP)),100)
        ISTEP_SPEZ=INT(SORIND(ISRFS,ISTRAI)/100)
C
        IF (ISTEP.EQ.0) THEN
          WRITE (iunout,*) 'ERROR IN PRIMARY SOURCE DATA '
          WRITE (iunout,*) 'STEP FUNCTION REQUESTED FOR SOURCE SURFACE '
          WRITE (iunout,*) 'NO. ', INSOR(ISRFS,ISTRAI),
     .      ' BUT SORIND.EQ.0..'
          CALL EIRENE_EXIT_OWN(1)
        ELSEIF (ISTEP.GT.NSTEP) THEN
          CALL EIRENE_MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  HAS THIS STEP FUNCTION NO. ISTEP STILL TO BE INITIALIZED?
C
        IF (ISTUF(ISTEP).EQ.0) THEN
C  YES
          IF (TRCSOU) THEN
            CALL EIRENE_LEER(1)
            WRITE (IUNOUT,*) 'FOR STRATUM ISTRA: ',ISTRAI
            WRITE (iunout,*) 'SET DEFAULT STEP FUNCTION ISTEP: ',ISTEP
          ENDIF

C
C  INITIALIZE STEP FUNCTION NO. ISTEP BY DEFAULT MODEL:
C  DEFAULT MODEL FOR ISTEP: NPLSI STEP FUNCTIONS FROM BULK ION FLUXES
C
C  STEP FUNCTIONS FOR SAMPLING 1ST COORDINATE (X) FROM "RADIAL" DISTRIBUTIONS ONLY:
C       source surfaces must hence be either x-y or x-z surface
c       IN CASE OF UNSTRUCTURED GRIDS: LEVGEO 4 AND LEVGEO 5,
C       ALL NON-DEFAULT SURFACES ARE REGARDED as x-s surface, t=const
C STEP FUNCTIONS IN ALL OTHER CASES MUST BE SET ELSEWHERE (E.G. IN USER OR INTERFACE ROUTINES)
C
          IF (INDTEC(NL1J,ISTRAI).EQ.4) THEN
C
C  USE X-OR RADIAL DISTRIBUTION OF ION FLUX 0.5*NI(R,Y0,Z0)*CS(R,Y0,Z0)*(DELTA-Z)
C  WITH: CS = COMMON ION ACOUSTIC SPEED, PARALLEL TO B FIELD.
C  TAKE RADIAL (X) PROFILE OF PLASMA DATA ON A POLOIDAL (Y = CONST) SURFACE:
C  AT SOME GIVEN POLOIDAL (Y) POSITION IP
C  AT SOME GIVEN TOROIDAL (Z) POSITION IT
C  SCALE FLUX DENSITY WITH A TOROIDAL LENGTH, I.E.
C  WITH: DELTA-Z = LENGTH IN TOROIDAL OR Z-DIRECTION: EITHER "ZDF" OR "2 PI R"

!  INDSRF IS THE SURFACE NUMBER OF THIS SURFACE SOURCE SEGMENT.
!         IT CAN BE EITHER AN ADDITIONAL SURFACE (1,...NLIM)
!         OR A NON-DEFAULT STANDARD SURFACE (THEN  NLIM+1...NLIM+NSTS)

            INDSRF = 0
C  Y0:
            IF (INDIM(ISRFS,ISTRAI).NE.2 .OR. .NOT.NLPOL) THEN
              IP=1
            ELSEIF (INDIM(ISRFS,ISTRAI).EQ.2 .AND. NLPOL) THEN
              IP=INSOR(ISRFS,ISTRAI)
              do ists=1,nstsi
                if ((inump(ists,2) == ip) .and.
     .              (irpta(ists,1) <= ingrda(isrfs,istrai,1)) .and.
     .              (irpte(ists,1) >= ingrde(isrfs,istrai,1)))
     .             indsrf=nlim+ists
              end do
            ENDIF
C  Z0:
            IF (INDIM(ISRFS,ISTRAI).NE.3 .OR. .NOT.NLTOR) THEN
              IT=1
            ELSEIF (INDIM(ISRFS,ISTRAI).EQ.3 .AND. NLTOR) THEN
              IT=INSOR(ISRFS,ISTRAI)
              WRITE (iunout,*) 'DEF. STEP FUNCTIONS ONLY AVAILABLE FOR'
              WRITE (iunout,*) 'SAMPLING OF RADIAL (OR X-) COORDINATE'
              WRITE (iunout,*) 'ON A POLOIDAL (OR Y-) SURFACE '
              WRITE (iunout,*) 'DEFINE STEP FUNCTION ELSEWHERE, ',
     .                         'EG. IN SAMUSR'
              WRITE (iunout,*) 'CALL EIRENE_EXIT '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
C
C  CURRENTLY: IN STANDARD GEOMETRY BLOCK IBLOCK=1
C

            select case (LEVGEO)
            case (1)
              KAN=1
              KEN=NR1STM
              DO K=1,NR1STM
                IRSTEP(ISTEP,K)=K
                IASTEP(ISTEP,K)=0
                IBSTEP(ISTEP,K)=1
                IPSTEP(ISTEP,K)=IP
                IF (ISORFL == 1) IPSTEP(ISTEP,K)=IP-1
                ITSTEP(ISTEP,K)=IT
                RRSTEP(ISTEP,K)=RSURF(K)
                IF (NLTRZ) THEN
                  TORL(ISTEP,K)=ZDF
                ELSEIF (NLTRA.OR.NLTRT) THEN
                  TORL(ISTEP,K)=(RSURF(K+1)+RSURF(K))/2._DP
                  TORL(ISTEP,K)=TORL(ISTEP,K)*2._DP*PIA
                ENDIF
              ENDDO
              NBIN=NR1ST
            case (2:3)
cdr May 2020 to be done:
cdr Possible conflict here:  levgeo=2, nlcrc, then: XPOL,YPOL and BGLP may not exist.
cdr The y-surface is: pol. angle=const.
cdr Furthermore: we assume here that the "radial surface" is normal to the
cdr plasma flow direction, i.e. plasma flows in y-z-plane only.
cdr Therefore here no projection of target surface to flux tube is done.
cdr For inclinded targets (e.g. levgeo=4) this can be different.
cdr Generally: what about pitch ? Also needed?
              KAN=1
              KEN=NR1STM
              DO K=1,NR1STM
                IRSTEP(ISTEP,K)=K
                IASTEP(ISTEP,K)=0
                IBSTEP(ISTEP,K)=1
                IPSTEP(ISTEP,K)=IP
                IF (ISORFL == 1) IPSTEP(ISTEP,K)=IP-1
                ITSTEP(ISTEP,K)=IT
cdr may 2020
                if (levgeo.eq.2 .and. nlcrc) then
                  RRSTEP(ISTEP,K)=RSURF(K)
                else
                  RRSTEP(ISTEP,K)=BGLP(K,IP)
                endif
                IF (NLTRZ) THEN
                  TORL(ISTEP,K)=ZDF
                ELSEIF (NLTRA.OR.NLTRT) THEN
cdr may 2020
                  IF (levgeo.eq.2 .and. NLCRC) THEN
                    TORL(ISTEP,K)=(RSURF(K+1)+RSURF(K))/2._DP
                  ELSE
                    TORL(ISTEP,K)=(XPOL(K+1,IP)+XPOL(K,IP))/2._DP
                  ENDIF
                  TORL(ISTEP,K)=TORL(ISTEP,K)*2._DP*PIA
                ENDIF
              ENDDO
              NBIN=NR1ST
            case (4)
C  TRIANGULAR GRID. RRSURF IS INTEGRATED ALONG A SET OF TRIANGLE SIDES.
              K=0
              RRSTEP(ISTEP,1) = 0._DP
              INDSRF=INSOR(ISRFS,ISTRAI)
              IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)
              DO ITRI=1,NTRII
                DO IS=1,3
                  IF (INMTI(IS,ITRI) == INDSRF) THEN
                    K=K+1
                    IS1=IS+1
                    IF (IS1.GT.3) IS1=1
                    IRSTEP(ISTEP,K)=ITRI
                    IPSTEP(ISTEP,K)=IS
                    ITSTEP(ISTEP,K)=IT
                    IASTEP(ISTEP,K)=0
                    IBSTEP(ISTEP,K)=1
!pb  projection to B field switched off!!!
!pb  to allow for step functions on surfaces perpendicular to magnetic field
cdr  ??
!pb                 BABS=SQRT(BXIN(ITRI)**2+BYIN(ITRI)**2+BZIN(ITRI)**2)
!pb                 CTHETA=ABS((PTRIX(IS,ITRI)*BXIN(ITRI) +
!pb     .                       PTRIY(IS,ITRI)*BYIN(ITRI))/BABS)
                    CTHETA = 1._DP
cdr  May 2020:  I do not think that is is correct.
cdr  theta is the angle between parallel (to B) plasma flux and the surface normal.
cdr  in levgeo=4 this angle can be nonzero, i.e. cos(theta) ne.1., even
cdr  ctheta=0 (on surfaces parallel to B) is possible.

                    DELRR=SQRT(
     .              (XTRIAN(NECKE(IS,ITRI))-XTRIAN(NECKE(IS1,ITRI)))**2+
     .              (YTRIAN(NECKE(IS,ITRI))-YTRIAN(NECKE(IS1,ITRI)))**2)

                    RRSTEP(ISTEP,K+1)=RRSTEP(ISTEP,K) + CTHETA*DELRR
                    IF (NLTRZ) THEN
                      TORL(ISTEP,K)=ZDF
                    ELSEIF (NLTRA.OR.NLTRT) THEN
                      TORL(ISTEP,K)=
     .                  (XTRIAN(NECKE(IS,ITRI))+XTRIAN(NECKE(IS1,ITRI)))
     .                  /2._DP
                      TORL(ISTEP,K)=TORL(ISTEP,K)*2._DP*PIA
                    ENDIF
                  END IF
                END DO
              END DO
              KAN=1
              KEN=K
              NBIN=K+1
            case (5)
C  GRID OF TETRAHEDRA. RRSTEP IS A CUMULATED SURFACE AREA, INTEGRATING
C                       OVER A SET OF TETRAHEDRON SIDES (=TRIANGLES)
cdr  i.e. no DELTA_Z factors here.
              K=0
              RRSTEP(ISTEP,1) = 0._DP
              INDSRF=INSOR(ISRFS,ISTRAI)
              IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)
              DO ITET=1,NTET
                DO IS=1,4
                  IF (INMTIT(IS,ITET) == INDSRF) THEN
cdr set irstep, ipstep and cumulated variable RRSTEP, purely geometrical.
cdr no projection to flux tube cross-section is done.
                    CALL EIRENE_TET_STEP (ISTEP,ITET,IS,K)
                  END IF
                END DO
              END DO
              KAN=1
              KEN=K
C  toroidal length: already included in RRSTEP, which is a surface area
              TORL(ISTEP,KAN:KEN)=1._DP
              NBIN=K+1
            case (10)
              WRITE (iunout,*) 'DEF. STEP FUNCTIONS NOT AVAILABLE FOR'
              WRITE (iunout,*) 'USER-SUPPLIED GEOMETRY OPTION '
              WRITE (iunout,*) 'DEFINE STEP FUNCTION ELSEWHERE, ',
     .                         'EG. IN SAMUSR'
              WRITE (iunout,*) 'CALL EXIT '
              CALL EIRENE_EXIT_OWN(1)
            end select
C
C  NOW SET THE FLUX DISTRIBUTION FLSTEP, AS WELL AS TE, TI, V-PLASMA, NI, at surface.
cdr also set sheath potential SHSTEP, and energy flux ELSTEP.
cdr Set these STEP functions in default species index range ispz=1,npls.
cdr This made sense only for strata with NLPLS=.TRUE
cdr July 22: extend code towards nlatm, nlmol, etc. sources
C

            DO K=KAN,KEN
              NBLCKA=NSTRD*(IBSTEP(ISTEP,K)-1)+IASTEP(ISTEP,K)
              select case (LEVGEO)
              case (4:5)
                NCELL=IRSTEP(ISTEP,K)+
     .            (ITSTEP(ISTEP,K)-1)*NP2T3*
     .             NR1P2+NBLCKA
              case default
                NCELL=IRSTEP(ISTEP,K)+
     .            ((IPSTEP(ISTEP,K)-1)+(ITSTEP(ISTEP,K)-1)*NP2T3)*
     .             NR1P2+NBLCKA
              end select

              TESTEP(ISTEP,K)=TEIN(NCELL)
              IF (NLPLS(ISTRAI)) THEN
cdr  the code below sets bulk ion (field particle) fluxes in the full species range 1:npls
               DO 2 JPLS=1,NPLSI
                IPLS=JPLS
                IPLSTI = MPLSTI(IPLS)
                IPLSV = MPLSV(IPLS)
                TISTEP(IPLSTI,ISTEP,K)=TIIN(IPLSTI,NCELL)
                IF (INDPRO(4) == 8) THEN
cdr  no spatial information x0,y0,z0 of birth point available here for vecusr.
c    set drift velocities at cell center
                  XC=0.
                  YC=0.
                  ZC=0.
                  CALL EIRENE_VECUSR(2,NCELL,XC,YC,ZC,VX,VY,VZ,IPLS,
     .                               .FALSE.)

                  VXSTEP(IPLSV,ISTEP,K)=VX
                  VYSTEP(IPLSV,ISTEP,K)=VY
                  VZSTEP(IPLSV,ISTEP,K)=VZ
                ELSE
                  VXSTEP(IPLSV,ISTEP,K)=VXIN(IPLSV,NCELL)
                  VYSTEP(IPLSV,ISTEP,K)=VYIN(IPLSV,NCELL)
                  VZSTEP(IPLSV,ISTEP,K)=VZIN(IPLSV,NCELL)
                END IF
                DISTEP(IPLS,ISTEP,K)=DIIN(IPLS,NCELL)
                IF (ZIIN(IPLS,NCELL).NE.ZVAC) THEN
                  ZISTEP(IPLS,ISTEP,K)=ZIIN(IPLS,NCELL)
                ELSE
                  ZISTEP(IPLS,ISTEP,K)=DBLE(NCHRGP(IPLS))
                END IF
c  isothermal ion acoustic speed
                CS=CVEL2A*SQRT((TIIN(IPLSTI,NCELL)+TEIN(NCELL))/
     .                          RMASSP(IPLS))
c  default: sonic flux density FL = 1/2 n_u cs = n_t cs  (AMPS)
                FLSTEP(IPLS,ISTEP,K)=CS*DIIN(IPLS,NCELL)*
     .                               TORL(ISTEP,K)*ELCHA
                IF (INDSRF == 0) THEN
                   SHSTEP(ISTEP,K) = 0._DP
                ELSE
                   SHSTEP(ISTEP,K)=FSHEAT(INDSRF)
                END IF
c  drifting maxwellian at sheath entrance: (not including sheath acceleration)
c  drift speed = isothermal ion sound speed
                ELSTEP(IPLS,ISTEP,K)=(3._DP*TISTEP(IPLSTI,ISTEP,K) +
     .                               0.5_DP*TESTEP(ISTEP,K)) *
     .                               FLSTEP(IPLS,ISTEP,K)
    2          CONTINUE

              ELSEIF (NLMOL(ISTRAI)) THEN
cdr  the code below sets molecular fluxes in the full species range 1:nmol
               DO 3 JMOL=1,NMOLI
                IMOL=JMOL

                TISTEP(IMOL,ISTEP,K)=0.026 ! room temperature kT, eV
                VXSTEP(IMOL,ISTEP,K)=0.0
                VYSTEP(IMOL,ISTEP,K)=0.0
                VZSTEP(IMOL,ISTEP,K)=0.0

                DISTEP(IMOL,ISTEP,K)=1.0  ! uniform density, a.u.
                ZISTEP(IMOL,ISTEP,K)=0.0
c  default: thermal flux density FL = 1/4 n v_mean  (AMPS)
                V_mean=CVEL2A*SQRT(8.0/PIA)*
     .                 SQRT(TISTEP(IMOL,ISTEP,K)/RMASSM(IMOL))
                FLSTEP(IMOL,ISTEP,K)=1./4.*V_mean*DISTEP(IMOL,ISTEP,K)*
     .                                            TORL(ISTEP,K)*ELCHA

                SHSTEP(ISTEP,K)=0.0
cdr half sided stationary Maxwellian energy flux= 2 T FL
                ELSTEP(IMOL,ISTEP,K)=2._DP*TISTEP(IMOL,ISTEP,K)*
     .                                     FLSTEP(IMOL,ISTEP,K)
    3          CONTINUE

              ELSEIF (NLATM(ISTRAI)) THEN
cdr  the code below sets atomic fluxes in the full species range 1:natm
               DO JATM=1,NATMI
                IATM=JATM

                TISTEP(IATM,ISTEP,K)=0.026 ! room temperature kT, eV
                VXSTEP(IATM,ISTEP,K)=0.0
                VYSTEP(IATM,ISTEP,K)=0.0
                VZSTEP(IATM,ISTEP,K)=0.0

                DISTEP(IATM,ISTEP,K)=1.0  ! uniform density, a.u.
                ZISTEP(IATM,ISTEP,K)=0.0
c  default: thermal flux density FL = 1/4 n v_mean  (AMPS)
                V_mean=CVEL2A*SQRT(8.0/PIA)*
     .                 SQRT(TISTEP(IATM,ISTEP,K)/RMASSA(IATM))
                FLSTEP(IATM,ISTEP,K)=1./4.*V_mean*DISTEP(IATM,ISTEP,K)*
     .                                            TORL(ISTEP,K)*ELCHA

                SHSTEP(ISTEP,K)=0.0
cdr half sided stationary Maxwellian energy flux= 2 T FL
                ELSTEP(IATM,ISTEP,K)=2._DP*TISTEP(IATM,ISTEP,K)*
     .                                     FLSTEP(IATM,ISTEP,K)
               ENDDO
              ELSE
                WRITE (IUNOUT,*)
     .           'DEFAULT STEP FUNCTION, BUT UNEXPECTED PARTICLE TYPE'
                CALL EIRENE_EXIT_OWN(1)
              ENDIF

            END DO
C
C  LAST INTERVAL BOUNDARY FOR SAMPLING DISTRIBUTION
C
            select case (LEVGEO)
            case (1)
              RRSTEP(ISTEP,NR1ST)=RSURF(NR1ST)
            case (2:3)
              RRSTEP(ISTEP,NR1ST)=BGLP(NR1ST,IP)
            case (4:)
C             RRSTEP(ISTEP,NBIN) ALREADY SET ABOVE
            end select

cdr finalize step function preparation:
            if (nlatm(istrai)) FL=EIRENE_STEP(1,NATMI,NBIN,ISTEP,1)
            if (nlmol(istrai)) FL=EIRENE_STEP(1,NMOLI,NBIN,ISTEP,2)
            if (nlion(istrai)) FL=EIRENE_STEP(1,NIONI,NBIN,ISTEP,3)
            if (nlpls(istrai)) FL=EIRENE_STEP(1,NPLSI,NBIN,ISTEP,4)
C
c  NO DEFAULT STEP FUNCTIONS AVAILABLE FOR Y OR Z SURFACE SOURCES
C
          ELSEIF (INDTEC(NL2J,ISTRAI).EQ.4) THEN
            WRITE (iunout,*) 'DEFAULT STEP FUNCTIONS ONLY AVAILABLE FOR'
            WRITE (iunout,*) 'SAMPLING OF RADIAL COORDINATE '
            WRITE (iunout,*)
     .        'DEFINE STEP FUNCTION ELSEWHERE, EG. IN SAMUSR'
            WRITE (iunout,*) 'CALL EXIT '
            CALL EIRENE_EXIT_OWN(1)
C
          ELSEIF (INDTEC(NL3J,ISTRAI).EQ.4) THEN
            WRITE (iunout,*) 'DEFAULT STEP FUNCTIONS ONLY AVAILABLE FOR'
            WRITE (iunout,*) 'SAMPLING OF RADIAL COORDINATE '
            WRITE (iunout,*)
     .        'DEFINE STEP FUNCTION ELSEWHERE, EG. IN SAMUSR'
            WRITE (iunout,*) 'CALL EXIT '
            CALL EIRENE_EXIT_OWN(1)
          ENDIF

        ELSEIF (ISTUF(ISTEP).NE.0) THEN
C  NO, ISTEP IS ALREADY DONE.
C  TBD: CHECK FOR SAME NLPLS,NLMOL,, ETC... AS IN EARLIER ISTEP SETTING
          IF (TRCSOU) THEN
            CALL EIRENE_LEER(1)
            WRITE (IUNOUT,*) 'FOR STRATUM ISTRA: ',ISTRAI
            WRITE (iunout,*) 'STEP FUNCTION ISTEP: ',ISTEP
            WRITE (IUNOUT,*) 'WAS ALREADY SET '
            CALL EIRENE_LEER(1)
          ENDIF

        ENDIF  !dr  Now: istuf=0 --> istuf=1

C  DEFAULT INITIALIZATION OF STEP FUNCTION ISTEP IN SAMSRF IS COMPLETE NOW.
c  ALTERNATIVELY: STEP FUNCTIONS MAY HAVE BEEN DEFINED FROM EXTERNAL MODULES.

C  AVAILABLE SPECIES RANGE FOR THIS ISTEP: NSPSTI,NSPSTE = 1, NPLSI IS SET
C  INDICATOR: ISTUF(ISTEP)=1 IS SET.

C  FUNCTION STEP(ISTEP,...) CAN BE USED FOR PRESENT STRATUM ISTRAI,
C  BUT ALSO FOR LATER STRATA.
c
c  PREPARE SOME STEP FUNCTION DIAGNOSTIC OUTPUT (SAME FOR DEFAULT OR EXTERNAL STEP FUNCTIONS)
C
        IF (TRCSOU.AND.ISTEP.GT.0 .AND.
cdr  printout only once for each step function ISTEP
cdr  Hidden link: istep and istrai are now not independent.
cdr  tbd:  make sure that for a given ISTEP
cdr        the istrai variables: NLPLS, NLMOL, NLATM,...etc. are identical
     .      (MSTEP(ISTEP).EQ.0 .OR. MMSTEP(ISTEP).EQ.1)) THEN
          MSTEP(ISTEP)=1
          MMSTEP(ISTEP)=0

          NSMX=NSMAX(ISTEP)
C  IDENTIFY SPECIES WITH NONZERO FLUX
          FLX=0.
          EKFLX=0.
          ESHFLX=0.
          DO JSPZ=NSPSTI(ISTEP),NSPSTE(ISTEP)
            ISPZ=JSPZ
            DO K=1,NSMX-1
              DELR=RRSTEP(ISTEP,K+1)-RRSTEP(ISTEP,K)
              FL=FLSTEP(ISPZ,ISTEP,K)*DELR
              FLX(ISPZ)=FLX(ISPZ)+FL
              EKFLX(ISPZ)=EKFLX(ISPZ)+ELSTEP(ISPZ,ISTEP,K)*DELR
c  sheath factor given on step function shstep along target?
              IF (SHSTEP(ISTEP,K) > 0) THEN
c  this can happen with either nlpls or nlion:
                ZZ=0.
                if (nlpls(istrai)) ZZ=ZISTEP(ISPZ,ISTEP,K)
                if (nlion(istrai)) ZZ=nchrgi(ispz)
                ESHFLX(ISPZ)=ESHFLX(ISPZ)+
     .               FL*SHSTEP(ISTEP,K)*TESTEP(ISTEP,K)*
     .               ZZ
              ELSEIF (NLPLS(ISTRAI)) THEN
c  employ default eirene sheath model.
c  ISPZ=IPLS now.
c  To be done. plasma flow velocity v..step should first be projected
c  towards surface normal.
                GAMMA=0.
                CUR=0.
                TESH=TESTEP(ISTEP,K)
CDR  THIS NEXT LOOP CAN GO OUT: IT IS NEEDED ONLY ONCE, NOT FOR EACH IPLS.
                DO IP=1,NPLSI
                  IPLSV = MPLSV(ISPZ)
                  VPSH(IP)=SQRT(VXSTEP(IPLSV,ISTEP,K)**2
     .                         +VYSTEP(IPLSV,ISTEP,K)**2
     .                         +VZSTEP(IPLSV,ISTEP,K)**2)
                  DISH(IP) = DISTEP(IP,ISTEP,K)
                  ZISH(IP) = ZISTEP(IP,ISTEP,K)
                END DO
CDR
                ESHFLX(ISPZ)=ESHFLX(ISPZ)+
     .             FL*ZISH(ISPZ)*EIRENE_SHEATH(TESH,DISH,VPSH,
     .                              ZISH,GAMMA,CUR,NPLSI,INDSRF)
              END IF  ! NO SHEATH
            ENDDO  ! ISPZ
          ENDDO

          NANZ=COUNT(FLX(:).GT.0.D0)
          ISPZD(1:NANZ)=PACK((/(ISP,ISP=NSPSTI(ISTEP),NSPSTE(ISTEP))/),
     .                  FLX(NSPSTI(ISTEP):NSPSTE(ISTEP)).GT.0.D0)
C
          WRITE (iunout,*) 'FUNCTION STEP NO. ',ISTEP,': '
          WRITE (IUNOUT,*) 'FLUXES IN AMP/CM**2 '
          IF (NANZ.LE.5) THEN
            WRITE (iunout,
     .       '(1X,A4,A12,5(2X,A7,I2,A1))')
     .       '   K','  RRSTEP    ',('FLSTEP(',ISPZD(ISP),')',
     .                               ISP=1,NANZ)
          ELSE
            WRITE (iunout,
     .       '(1X,A4,A12,5(2X,A7,I2,A1)/(17x,5(2X,A7,I2,A1)))')
     .       '   K','  RRSTEP    ',('FLSTEP(',ISPZD(ISP),')',
     .                               ISP=1,NANZ)
          END IF
          DO 4 K=1,NSMX-1
            IF (NANZ.LE.5) THEN
              WRITE (iunout,'(1X,I4,1P,6E12.4)')
     .               K,RRSTEP(ISTEP,K),
     .               (FLSTEP(ISPZD(ISP),ISTEP,K),ISP=1,NANZ)
            ELSE
              WRITE (iunout,'(1X,I4,1P,6E12.4/(17x,1P,5E12.4))')
     .               K,RRSTEP(ISTEP,K),
     .               (FLSTEP(ISPZD(ISP),ISTEP,K),ISP=1,NANZ)
            END IF
    4     CONTINUE

          WRITE (iunout,'(1X,I4,1P,2E12.4)') NSMX,RRSTEP(ISTEP,NSMX)
          CALL EIRENE_LEER(1)

          WRITE (iunout,*) 'FLUXES: PART., KINET., SHEATH; INTEGRATED:'
          DO 5 ISP=1,NANZ
            CALL EIRENE_MASJ2R('ISTEP,ISPZ,FLUX      [A]',
     .                          ISTEP,ISPZD(ISP),FLX(ISPZD(ISP)))
            CALL EIRENE_MASJ2R('ISTEP,ISPZ,EKIN-FLUX [W]',
     .                          ISTEP,ISPZD(ISP),EKFLX(ISPZD(ISP)))
            CALL EIRENE_MASJ2R('ISTEP,ISPZ,ESH-FLUX  [W]',
     .                          ISTEP,ISPZD(ISP),ESHFLX(ISPZD(ISP)))
    5     CONTINUE
          CALL EIRENE_LEER(2)
        ENDIF

C  NEXT: CONTINUE INITIALIZATION FOR STRATUM ISTRAI.
C
C
    7   CONTINUE
C
C  DEFINE LEFT AND RIGHT BOUNDARY OF SAMPLING INTERVALS FOR STRATUM ISTRAI.
C
C  FLAG ISTEP INDICATES: STEP FUNCTION ISTEP IS TO BE USED FOR THIS STRATUM
cdr  next 2 lines: same as already done above.
        ISTEP=MOD(NINT(REAL(SORIND(ISRFS,ISTRAI),DP)),100)
        ISTEP_SPEZ=INT(SORIND(ISRFS,ISTRAI)/100)
C
        IF (INDIM(ISRFS,ISTRAI).EQ.1) THEN
C  SOURCE ON RADIAL (X=CONST) SURFACE:
C                  SAMPLE SECOND AND THIRD COORDINATE,
C                  COMPUTE FIRST COORDINATE
C
          INDTEC(NL1J,ISTRAI)=0
C
          ISGRD2=INGRDA(ISRFS,ISTRAI,2)+INGRDE(ISRFS,ISTRAI,2)
          ISGRD3=INGRDA(ISRFS,ISTRAI,3)+INGRDE(ISRFS,ISTRAI,3)

          IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.ISGRD2.GT.0) THEN
            IS1=MAX0(1,INGRDA(ISRFS,ISTRAI,2))
            IS2=MIN0(NRPLG,INGRDE(ISRFS,ISTRAI,2))
            ALEFT(2,ISRFS,ISTRAI)=BGL(INSOR(ISRFS,ISTRAI),IS1)
            BRGHT(2,ISRFS,ISTRAI)=BGL(INSOR(ISRFS,ISTRAI),IS2)
          ELSEIF (LEVGEO.EQ.4.OR.LEVGEO.EQ.5) THEN
            IF (ISTEP == 0) THEN
              select case (LEVGEO)
              case (4)
                INS = INSOR(ISRFS,ISTRAI)
                IF (INS < 0) INS = ABS(INS) + NLIM
                ALEFT(2,ISRFS,ISTRAI)=SURF_TRIAN(INS)%BGLT(1)
                BRGHT(2,ISRFS,ISTRAI)=SURF_TRIAN(INS)%
     .                               BGLT(SURF_TRIAN(INS)%NUMTR+1)
              case (5)
                WRITE (iunout,*) ' ERROR IN SAMSRF '
                WRITE (iunout,*)
     .            ' SAMPLING ON NON-DEFAULT STANDARD X-SURFACE '
                WRITE (iunout,*) ' IS NOT FORESEEN FOR THIS LEVGEO '
                WRITE (iunout,*) ' EXCEPT FOR STEP FUNCTION OPTIONS '
                CALL EIRENE_EXIT_OWN(1)
              end select
c  STEP FUNCTIONS ON RADIAL SURFACE MAY HAVE BEEN DEFINED EXTERNALLY
C  TAKE FULL RANGE STEP FUNCTION. NO USE OF INGRDA, INGRDE FORESEEN
            ELSE IF (ISTEP.NE.0) THEN

              IF (NSMAX(ISTEP).NE.0) THEN
                ALEFT(2,ISRFS,ISTRAI)=RRSTEP(ISTEP,1)
                BRGHT(2,ISRFS,ISTRAI)=RRSTEP(ISTEP,NSMAX(ISTEP))
              ELSE
                WRITE (iunout,*) ' ERROR IN SAMSRF '
                WRITE (iunout,*)
     .            ' CANNOT FIND RANGE FOR STEP FUNCTION SAMPLING '
                WRITE (IUNOUT,*) ' BECAUSE NSMAX(ISTEP) = 0  '
                CALL EIRENE_EXIT_OWN(1)
              END IF
            END IF
C
          ELSEIF (LEVGEO == 1) THEN
            ALEFT(2,ISRFS,ISTRAI)=SORAD3(ISRFS,ISTRAI)
            BRGHT(2,ISRFS,ISTRAI)=SORAD4(ISRFS,ISTRAI)
          ENDIF
          ALEFT(3,ISRFS,ISTRAI)=SORAD5(ISRFS,ISTRAI)
          BRGHT(3,ISRFS,ISTRAI)=SORAD6(ISRFS,ISTRAI)
        ENDIF
C
        IF (INDIM(ISRFS,ISTRAI).EQ.2) THEN
C  SOURCE ON POLOIDAL (Y0= CONST) SURFACE:
C            SAMPLE FIRST AND THIRD COORDINATE,
C            COMPUTE SECOND COORDINATE
C
          INDTEC(NL2J,ISTRAI)=0
C
          ISGRD1=INGRDA(ISRFS,ISTRAI,1)+INGRDE(ISRFS,ISTRAI,1)
          ISGRD3=INGRDA(ISRFS,ISTRAI,3)+INGRDE(ISRFS,ISTRAI,3)
          IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.ISGRD1.GT.0) THEN
            IS1=MAX0(1,INGRDA(ISRFS,ISTRAI,1))
            IS2=MIN0(NR1ST,INGRDE(ISRFS,ISTRAI,1))
            ALEFT(1,ISRFS,ISTRAI)=BGLP(IS1,INSOR(ISRFS,ISTRAI))
            BRGHT(1,ISRFS,ISTRAI)=BGLP(IS2,INSOR(ISRFS,ISTRAI))
          ELSE
            ALEFT(1,ISRFS,ISTRAI)=SORAD1(ISRFS,ISTRAI)
            BRGHT(1,ISRFS,ISTRAI)=SORAD2(ISRFS,ISTRAI)
          ENDIF
          ALEFT(3,ISRFS,ISTRAI)=SORAD5(ISRFS,ISTRAI)
          BRGHT(3,ISRFS,ISTRAI)=SORAD6(ISRFS,ISTRAI)
        ENDIF
C
        IF (INDIM(ISRFS,ISTRAI).EQ.3) THEN
C  SOURCE ON TOROIDAL (Z=CONST) SURFACE:
C                  SAMPLE FIRST AND SECOND COORDINATE,
C                  COMPUTE THIRD COORDINATE
C
          INDTEC(NL3J,ISTRAI)=0
C
          ISGRD1=INGRDA(ISRFS,ISTRAI,1)+INGRDE(ISRFS,ISTRAI,1)
          ISGRD2=INGRDA(ISRFS,ISTRAI,2)+INGRDE(ISRFS,ISTRAI,2)
C         IF (LEVGEO.EQ.3.AND.ISGRD1.GT.0) THEN
C         ELSE
            ALEFT(1,ISRFS,ISTRAI)=SORAD1(ISRFS,ISTRAI)
            BRGHT(1,ISRFS,ISTRAI)=SORAD2(ISRFS,ISTRAI)
C         ENDIF
          ALEFT(2,ISRFS,ISTRAI)=SORAD3(ISRFS,ISTRAI)
          BRGHT(2,ISRFS,ISTRAI)=SORAD4(ISRFS,ISTRAI)
        ENDIF
C
        IF (INDIM(ISRFS,ISTRAI).EQ.4) THEN
C  SOURCE ON COMBINATION OF RADIAL AND POLOIDAL SURFACE SEGMENTS:
C  ONLY FOR EXTERNALLY DEFINED STEP FUNCTIONS
C      SAMPLE ARC-LENGTH AND THIRD COORDINATE,
C      COMPUTE FIRST AND SECOND COORDINATE FROM ARCLENGTH
C
          INDTEC(NL2J,ISTRAI)=0
C
          IF (ISTEP.NE.0) THEN
            IF (NSMAX(ISTEP).NE.0) THEN
              ALEFT(1,ISRFS,ISTRAI)=RRSTEP(ISTEP,1)
              BRGHT(1,ISRFS,ISTRAI)=RRSTEP(ISTEP,NSMAX(ISTEP))
            ELSE
              WRITE (iunout,*) ' ERROR IN SAMSRF '
              WRITE (iunout,*)
     .           ' INDIM=4 ONLY FORESEEN WITH STEP FUNCTION '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
          ELSE
            WRITE (iunout,*) ' ERROR IN SAMSRF '
            WRITE (iunout,*) ' INDIM=4 ONLY FORESEEN WITH STEP FUNCTION'
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
          ALEFT(3,ISRFS,ISTRAI)=SORAD5(ISRFS,ISTRAI)
          BRGHT(3,ISRFS,ISTRAI)=SORAD6(ISRFS,ISTRAI)
        ENDIF
C
        IF (INDIM(ISRFS,ISTRAI).EQ.0) THEN
C  SOURCE ON ADDITIONAL SURFACE
          ALEFT(1,ISRFS,ISTRAI)=SORAD1(ISRFS,ISTRAI)
          BRGHT(1,ISRFS,ISTRAI)=SORAD2(ISRFS,ISTRAI)
          ALEFT(2,ISRFS,ISTRAI)=SORAD3(ISRFS,ISTRAI)
          BRGHT(2,ISRFS,ISTRAI)=SORAD4(ISRFS,ISTRAI)
          ALEFT(3,ISRFS,ISTRAI)=SORAD5(ISRFS,ISTRAI)
          BRGHT(3,ISRFS,ISTRAI)=SORAD6(ISRFS,ISTRAI)
        ENDIF
C
C
C
C  STRATUM INITIALISATION FINISHED.

C  (IF ISTEP GT. 0: STEP(ISTEP,...) WAS REQUESTED FOR STRATUM ISTRAI,
C  AND THIS STEP FUNCTION IS ALSO FULLY INITIALIZED.
C
C  NEXT: STORE SOME DATA TO SPEED UP SAMPLING FROM STEP FUNCTION
C  IN STRATUM ISTRAI, FOR THE SELECTED (RANGE OF) SPECIES.
C
        ISPZ=0
        IF (NSPEZ(ISTRAI).GT.0) ISPZ=NSPEZ(ISTRAI)

        IF (ISTEP_SPEZ.GT.0) ISPZ=ISTEP_SPEZ
C
        IF (ISTEP.GT.0) THEN
          IF (ISPZ.GT.0.AND.ISPZ.LT.NSPSTI(ISTEP)) GOTO 991
          IF (              ISPZ.GT.NSPSTE(ISTEP)) GOTO 991
C
          IF (INDTEC(NL1J,ISTRAI).EQ.4) THEN
            XI(1,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,ALEFT(1,ISRFS,ISTRAI))
            XE(1,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,BRGHT(1,ISRFS,ISTRAI))
            IF (TRCSOU) THEN
              WRITE (IUNOUT,*) 'ISPZ, STEP FUNCTION SAMPLING INTERVAL'
              CALL EIRENE_MASJR2('ISPZ, ALEFT, BRGHT      ', ISPZ,
     .                     ALEFT(1,ISRFS,ISTRAI),BRGHT(1,ISRFS,ISTRAI))
              WRITE (IUNOUT,*) 'SUB-RANGE FOR UNIFORM RANDOM NUMBERS '
              CALL EIRENE_MASR2('XI,XE           ',
     .                       XI(1,ISRFS,ISTRAI),XE(1,ISRFS,ISTRAI))
              RANDIF=XE(1,ISRFS,ISTRAI)-XI(1,ISRFS,ISTRAI)
              IF (ABS(RANDIF-1.0).GE.1.0D-4) THEN
                DO ISP=1,NANZ
                  WRITE (IUNOUT,*)
     .              'TRUNCATED FLUX ONTO SOURCE SURFACE [A] ',
     .                   ISPZD(ISP),RANDIF*FLX(ISPZD(ISP))
                END DO
              ENDIF
            END IF
          ELSEIF (INDTEC(NL2J,ISTRAI).EQ.4) THEN
            XI(2,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,ALEFT(2,ISRFS,ISTRAI))
            XE(2,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,BRGHT(2,ISRFS,ISTRAI))
            IF (TRCSOU) THEN
              WRITE (IUNOUT,*) 'ISPZ, STEP FUNCTION SAMPLING INTERVAL'
              CALL EIRENE_MASJR2('ISPZ, ALEFT, BRGHT      ', ISPZ,
     .                     ALEFT(2,ISRFS,ISTRAI),BRGHT(2,ISRFS,ISTRAI))
              WRITE (IUNOUT,*) 'SUB-RANGE FOR UNIFORM RANDOM NUMBERS '
              CALL EIRENE_MASR2('XI,XE           ',
     .                       XI(2,ISRFS,ISTRAI),XE(2,ISRFS,ISTRAI))
              DO ISP=1,NANZ
                WRITE (IUNOUT,*) 'FLUX ONTO SOURCE SURFACE [A] ',ISP,
     .                (XE(2,ISRFS,ISTRAI)-XI(2,ISRFS,ISTRAI))*
     .                 FLX(ISPZD(ISP))
              END DO
            END IF
          ELSEIF (INDTEC(NL3J,ISTRAI).EQ.4) THEN
            XI(3,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,ALEFT(3,ISRFS,ISTRAI))
            XE(3,ISRFS,ISTRAI)=
     .           EIRENE_STEP0(ISPZ,ISTEP,BRGHT(3,ISRFS,ISTRAI))
            IF (TRCSOU) THEN
              WRITE (IUNOUT,*) 'ISPZ, STEP FUNCTION SAMPLING INTERVAL'
              CALL EIRENE_MASJR2('ISPZ, ALEFT, BRGHT      ', ISPZ,
     .                     ALEFT(3,ISRFS,ISTRAI),BRGHT(3,ISRFS,ISTRAI))
              WRITE (IUNOUT,*) 'SUB-RANGE FOR UNIFORM RANDOM NUMBERS '
              CALL EIRENE_MASR2('XI,XE           ',
     .                       XI(3,ISRFS,ISTRAI),XE(3,ISRFS,ISTRAI))
              DO ISP=1,NANZ
                WRITE (IUNOUT,*) 'FLUX ONTO SOURCE SURFACE [A] ',ISP,
     .                (XE(3,ISRFS,ISTRAI)-XI(3,ISRFS,ISTRAI))*
     .                 FLX(ISPZD(ISP))
              END DO
            END IF
          ENDIF
        ENDIF
C
    1 CONTINUE
C
      IF (TRCSOU) CALL EIRENE_LEER(2)

      END DO    ! STRATA LOOP

      RETURN

  991 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF '
      WRITE (iunout,*) 'SPECIES INDEX DISTRIBUTION NOT KNOWN'
      WRITE (iunout,*) 'ISPZ, STRATUM NO: ISTRAI= ',ISPZ,ISTRAI
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_SAMSF0
C
      SUBROUTINE EIRENE_SAMSF1
     .      (NLSF,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,ZIWL,WEISPZ)
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: TEWL, SHWL, VXWL(*), VYWL(*), VZWL(*),
     .                         TIWL(*), DIWL(*), EFWL(*), WEISPZ(*),
     .                         ZIWL(*)
      INTEGER, INTENT(IN) :: NLSF
      REAL(DP) :: ZZ(3)
      REAL(DP) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, ELLZZ1, EP1ZZ1,
     .          WINK, XR, DET,
     .          X0TEST, Y0TEST, Z0TEST,
     .          S, AN, P, Q, VVX, D, VVI, VVY, BL, PH,
     .          RNF, ZH, DELTA, ZM, XLAMDA,
     .          EIRENE_STEP1
      INTEGER :: ISID, IDUM, NDUM, EIRENE_LEARC2, NT,
     .           IEN, IAN, ITRI, ITET,
     .           EIRENE_LEARCA,
     .           ICOUNT, IPLG, I, ILTR, IAUSR,
     .           IBUSR, IRUSR, IPUSR, ITUSR, IK, J, JCALC, IINDEX,
     .           JPLS, JSPZ
      LOGICAL :: LOGTST
      INTEGER :: ITSIDE(3,4)
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_FZRTOR, EIRENE_FZRTRI,
     .            EIRENE_FPOLYT_3, EIRENE_SURTST, EIRENE_VECUSR,
     .            EIRENE_MASJ3, EIRENE_MASR3, EIRENE_MASR4,
     .            EIRENE_LEARCA, EIRENE_LEARC2, EIRENE_STEP1
C
      ICOUNT=0
  100 CONTINUE
      ICOUNT=ICOUNT+1
C
C   NLSF=SURFACE INDEX IN (NSRFS) SOURCE ARRAYS
C
      JCALC=0
      ISTEP=0
C
      WEISPZ(1:NSPZ)=-1.
C
C   USER-SUPPLIED GENERATOR FOR X0,Y0,Z0. ONLY ONE CALL
C   FOR ALL 3 COORDINATES. SUBR. SURTST IS NOT CALLED!
C
      IF (SORLIM(NLSF,ISTRA).LT.0.D0) THEN
        CALL EIRENE_SAMUSR(NLSF,X0,Y0,Z0,
     .              SORAD1(NLSF,ISTRA),SORAD2(NLSF,ISTRA),
     .              SORAD3(NLSF,ISTRA),SORAD4(NLSF,ISTRA),
     .              SORAD5(NLSF,ISTRA),SORAD6(NLSF,ISTRA),
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,ZIWL,
     .              WEISPZ)
        ISTEP=-1
        ZZ(1)=X0
        ZZ(2)=Y0
        ZZ(3)=Z0
        GOTO 1500
      ENDIF
C
      ZZ = 0._DP
      DO J=1,3
        IK=NLSF+(J-1)*NSRFS
        select case (INDTEC(IK,ISTRA))
        case (1)
C   DELTA DISTRIBUTION AT CENTER OF INTERVAL
          ZZ(J)=(ALEFT(J,NLSF,ISTRA)+BRGHT(J,NLSF,ISTRA))*0.5
        case (2)
C   UNIFORM DISTRIBUTION IN THIS COORDINATE
          ZZ(J)=RANF_EIRENE( )*(BRGHT(J,NLSF,ISTRA)-
     .          ALEFT(J,NLSF,ISTRA))+ALEFT(J,NLSF,ISTRA)
        case (3)
C   TRUNCATED EXPONENTIAL DECAY WITH LENGTH XLAMDA, FOR ONE COORDINATE ONLY
C   PARAMETER: SOREXP
C   METHOD: COVEYOU-TRICK (SPANIER-GELBARD, ADDISON WESLEY, P 35)
          DELTA=BRGHT(J,NLSF,ISTRA)-ALEFT(J,NLSF,ISTRA)
          XLAMDA=SOREXP(NLSF,ISTRA)
          ZM=DELTA/XLAMDA
          ZH=MOD(-LOG(RANF_EIRENE( )),ZM)
          ZZ(J)=XLAMDA*ZH+ALEFT(J,NLSF,ISTRA)
        case (4)
C   STEP FUNCTION NO. ISTEP, FOR ONE COORDINATE ONLY
C   PARAMETER: SORIND
          ISTEP=MOD(NINT(REAL(SORIND(NLSF,ISTRA),DP)),100)
          ISTEP_SPEZ = INT(SORIND(NLSF,ISTRA)/100)
          ISPZ=NSPEZ(ISTRA)
          IF (ISTEP_SPEZ.GT.0) ISPZ=ISTEP_SPEZ
c  RNF: uniform in spatial sampling interval
          RNF=XI(J,NLSF,ISTRA)+RANF_EIRENE( )*
     .        (XE(J,NLSF,ISTRA)-XI(J,NLSF,ISTRA))
          ZZ(J)=EIRENE_STEP1(IINDEX,ISTEP,RNF,ISPZ)
        case default
C   ZZ(JCALC) IS TO BE CALCULATED FROM SURFACE EQUATION
          IF (JCALC.NE.0) GOTO 997
          JCALC=J
        end select
      END DO
C
      IPOLG=1
C
 1500 CONTINUE
C
C  1500---2000
C  IN CASE OF DEFAULT SAMPLING:
C    TWO OUT OF THE THREE COORDINATES ZZ(1),ZZ(2),ZZ(3) ARE
C    SAMPLED.
C    FIND 3RD COORDINATE FROM SURFACE EQUATION
C    CONVERT ZZ(1), ZZ(2), ZZ(3) TO CARTESIAN COORDINATES X0,Y0,Z0
C  IN CASE OF USER-SUPPLIED SAMPLING:
C    THE THREE CARTESIAN COORDINATES X0,Y0,Z0 ARE GIVEN NOW
C    ONLY SURFACE FLAGS NLSRFX,...MRSURF,.... MUST BE SET
C
      IF (INDIM(NLSF,ISTRA).EQ.0) THEN
C  BIRTH POINT ON ADDITIONAL SURFACE, ZZ( ) ARE CARTESIAN
        NLSRFX=.FALSE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
        MTSURF=0
        MPSURF=0
        MRSURF=0
        MASURF=INSOR(NLSF,ISTRA)
        ILTR=ILTOR(MASURF)
        IF (NLTRZ.OR.     ILTR.NE.0) THEN
C Z0 IS GIVEN IN LOCAL COORDINATE SYSTEM ILTOR
          Z0=ZZ(3)
        ELSEIF (NLTRA.AND.ILTR.EQ.0) THEN
          PHI=ZZ(3)*DEGRAD
C Z0 IS IRRELEVANT HERE, AND IS DETERMINED LATER FROM PHI
        ELSEIF (NLTRT.AND.ILTR.EQ.0) THEN
          PHI=ZZ(3)*DEGRAD
C Z0 IS IRRELEVANT HERE, AND IS DETERMINED LATER FROM PHI
        ENDIF
c
        IF (SORLIM(NLSF,ISTRA).LT.0.D0) GOTO 2000
C
        IF (JCALC.EQ.1) THEN
C  FIND X COORDINATE X0 FROM Y=Y0 AND Z=Z0 ON SURFACE NO. MASURF
          Y0=ZZ(2)
          IF (NLTRA.AND.ILTR.EQ.0) THEN
            IF(.NOT.(A3LM(MASURF).EQ.0.0_DP .AND.
     &               A6LM(MASURF).EQ.0.0_DP .AND.
     &               A8LM(MASURF).EQ.0.0_DP .AND.
     &               A9LM(MASURF).EQ.0.0_DP)) THEN
              WRITE (iunout,*) 'Z0 IN SAMSRF FOR JCALC=1 ?? '
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
          ENDIF
          IF (JUMLIM(MASURF).NE.0) THEN
            IF (ABS(A1LM(MASURF)).LE.EPS12) GOTO 9931
            AN=-1./A1LM(MASURF)
            X0=(A0LM(MASURF)+Y0*A2LM(MASURF)+Z0*A3LM(MASURF))*AN
          ELSE
            S=A4LM(MASURF)
            Q=A0LM(MASURF)+(A2LM(MASURF)+A5LM(MASURF)*Y0)*Y0+
     .        (A3LM(MASURF)+A6LM(MASURF)*Z0+A9LM(MASURF)*Y0)*Z0
            P=A1LM(MASURF)+A7LM(MASURF)*Y0+A8LM(MASURF)*Z0
            IF (ABS(S).GT.EPS12) THEN
C  SECOND ORDER IN X
              AN=1./S
              Q=Q*AN
              P=P*AN
              PH=-P*0.5
              DET=PH*PH-Q
              IF (DET.LT.0.D0) THEN
                GOTO 999
              ELSEIF (DET.EQ.0.D0) THEN
                X0=PH
              ELSEIF (DET.GT.0.D0) THEN
                X0=PH+SIGN(1._DP,RANF_EIRENE()-0.5_DP)*SQRT(DET)
                CALL EIRENE_SURTST(X0,Y0,Z0,MASURF,LOGTST)
                IF (.NOT.LOGTST) THEN
                  IF (ICOUNT.LE.1000) GOTO 100
                  WRITE (IUNOUT,*) 'PROBABLY ILL-DEFINED SURFACE '
                  WRITE (IUNOUT,*) 'SAMPLING, MASURF = ',MASURF
                  LGPART=.FALSE.
                  GOTO 998
                ENDIF
              ENDIF
            ELSEIF (ABS(P).GT.EPS12) THEN
C  FIRST ORDER IN X
              X0=Q/P
            ELSE
C  INDEPENDENT OF X
              GOTO 9931
            ENDIF
          ENDIF
C  CARRY OUT RANGE TEST FOR X0?
          X0TEST=ABS(ALEFT(1,NLSF,ISTRA)-BRGHT(1,NLSF,ISTRA))
          IF (X0TEST.LT.1.D-10) GOTO 2000
          IF (ALEFT(1,NLSF,ISTRA).GT.X0 .OR.
     .        BRGHT(1,NLSF,ISTRA).LT.X0) THEN
            IF (ICOUNT.LT.1000) GOTO 100
            WRITE (iunout,*)
     .        'WARNING FROM SAMSRF FROM X0TEST, ICOUNT=1000 '
            LGPART=.FALSE.
            GOTO 998
          ENDIF

        ELSEIF (JCALC.EQ.2) THEN
C    FIND Y COORDINATE Y0 FROM X=X0 AND Z=Z0 ON SURFACE NO. MASURF
          X0=ZZ(1)
          IF (NLTRA.AND.ILTR.EQ.0) THEN
C  FIND Z0, NT,  FROM X0,PHI
            IF (PHI.LT.ZSURF(1)) PHI=PHI+PI2A
            IF (PHI.GT.ZSURF(NTTRA)) PHI=PHI-PI2A
            NT=EIRENE_LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMSRF      ')
            IF (NT.LE.0.OR.NT.GT.NTTRAM) THEN
              WRITE (iunout,*) 'NT OUT OF RANGE IN SAMSRF '
              WRITE (iunout,*) PHI,ZFULL
              CALL EIRENE_EXIT_OWN(1)
            ENDIF
            X01=X0+RMTOR
            CALL EIRENE_FZRTRI(X0,Z0,NT,X01,PHI,NT)
          ENDIF
          IF (JUMLIM(MASURF).NE.0) THEN
            IF (ABS(A2LM(MASURF)).LE.EPS12) GOTO 9931
            AN=-1./A2LM(MASURF)
            Y0=(A0LM(MASURF)+X0*A1LM(MASURF)+Z0*A3LM(MASURF))*AN
          ELSE
            S=A5LM(MASURF)
            Q=A0LM(MASURF)+(A1LM(MASURF)+A4LM(MASURF)*X0)*X0+
     .        (A3LM(MASURF)+A6LM(MASURF)*Z0+A8LM(MASURF)*X0)*Z0
            P=A2LM(MASURF)+A7LM(MASURF)*X0+A9LM(MASURF)*Z0
            IF (ABS(S).GT.EPS12) THEN
C  SECOND ORDER IN Y
              AN=1./S
              Q=Q*AN
              P=P*AN
              PH=-P*0.5
              DET=PH*PH-Q
              IF (DET.LT.0.D0) THEN
                GOTO 999
              ELSEIF (DET.EQ.0.D0) THEN
                Y0=PH
              ELSEIF (DET.GT.0.D0) THEN
                Y0=PH+SIGN(1._DP,RANF_EIRENE()-0.5_DP)*SQRT(DET)
                CALL EIRENE_SURTST(X0,Y0,Z0,MASURF,LOGTST)
                IF (.NOT.LOGTST) THEN
                  IF (ICOUNT.LE.1000) GOTO 100
                  WRITE (IUNOUT,*) 'PROBABLY ILL-DEFINED SURFACE '
                  WRITE (IUNOUT,*) 'SAMPLING, MASURF = ',MASURF
                  LGPART=.FALSE.
                  GOTO 998
                ENDIF
              ENDIF
            ELSEIF (ABS(P).GT.EPS12) THEN
C  FIRST ORDER IN Y
              Y0=Q/P
            ELSE
C  INDEPENDENT OF Y
              GOTO 9931
            ENDIF
          ENDIF
C  CARRY OUT RANGE TEST FOR Y0?
          Y0TEST=ABS(ALEFT(2,NLSF,ISTRA)-BRGHT(2,NLSF,ISTRA))
          IF (Y0TEST.LT.1.D-10) GOTO 2000
          IF (ALEFT(2,NLSF,ISTRA).GT.Y0 .OR.
     .        BRGHT(2,NLSF,ISTRA).LT.Y0) THEN
            IF (ICOUNT.LT.1000) GOTO 100
            WRITE (iunout,*)
     .        'WARNING FROM SAMSRF FROM Y0TEST, ICOUNT=1000 '
            LGPART=.FALSE.
            GOTO 998
          ENDIF

        ELSEIF (JCALC.EQ.3) THEN
C    FIND Z COORDINATE Z0 FROM X=X0 AND Y=Y0 ON SURFACE NO. MASURF
          X0=ZZ(1)
          Y0=ZZ(2)
          IF ((NLTRA.OR.NLTRT).AND.ILTR.EQ.0) GOTO 9931
          IF (JUMLIM(MASURF).NE.0) THEN
            IF (ABS(A3LM(MASURF)).LE.EPS12) GOTO 9931
            AN=-1./A3LM(MASURF)
            Z0=(A0LM(MASURF)+X0*A1LM(MASURF)+Y0*A2LM(MASURF))*AN
          ELSE
            S= A6LM(MASURF)
            Q= A0LM(MASURF)+(A1LM(MASURF)+A4LM(MASURF)*X0)*X0+
     .        (A2LM(MASURF)+A5LM(MASURF)*Y0+A7LM(MASURF)*X0)*Y0
            P= A3LM(MASURF)+A8LM(MASURF)*X0+A9LM(MASURF)*Y0
            IF (ABS(S).GT.EPS12) THEN
C  SECOND ORDER IN Z
              AN=1./S
              Q=Q*AN
              P=P*AN
              PH=-P*0.5
              DET=PH*PH-Q
              IF (DET.LT.0.D0) THEN
                GOTO 999
              ELSEIF (DET.EQ.0.D0) THEN
                Z0=PH
              ELSEIF (DET.GT.0.D0) THEN
                Z0=PH+SIGN(1._DP,RANF_EIRENE()-0.5_DP)*SQRT(DET)
                CALL EIRENE_SURTST(X0,Y0,Z0,MASURF,LOGTST)
                IF (.NOT.LOGTST) THEN
                  IF (ICOUNT.LE.1000) GOTO 100
                  WRITE (IUNOUT,*) 'PROBABLY ILL-DEFINED SURFACE '
                  WRITE (IUNOUT,*) 'SAMPLING, MASURF = ',MASURF
                  LGPART=.FALSE.
                  GOTO 998
                ENDIF
              ENDIF
            ELSEIF (ABS(P).GT.EPS12) THEN
C  FIRST ORDER IN Z
              Z0=Q/P
            ELSE
C  INDEPENDENT OF Z
              GOTO 9931
            ENDIF
          ENDIF
C  CARRY OUT RANGE TEST FOR Z0?
          Z0TEST=ABS(ALEFT(3,NLSF,ISTRA)-BRGHT(3,NLSF,ISTRA))
          IF (Z0TEST.LT.1.D-10) GOTO 2000
          IF (ALEFT(3,NLSF,ISTRA).GT.Z0 .OR.
     .        BRGHT(3,NLSF,ISTRA).LT.Z0) THEN
            IF (ICOUNT.LT.1000) GOTO 100
            WRITE (iunout,*) 'WARNING FROM SAMSRF, Z0TEST, ICOUNT=1000 '
            LGPART=.FALSE.
            GOTO 998
          ENDIF
C
        ELSE
          X0=ZZ(1)
          Y0=ZZ(2)
          Z0=ZZ(3)
        ENDIF
C
        CALL EIRENE_SURTST(X0,Y0,Z0,MASURF,LOGTST)
        IF (.NOT.LOGTST) GOTO 998
C
      ELSEIF (INDIM(NLSF,ISTRA).EQ.1) THEN
C  BIRTH POINT ON STANDARD RADIAL SURFACE  MRSURF
C  Y- OR POLOIDAL AND Z- OR TOROIDAL COORDINATES MUST BE SAMPLED
C  X-OR RADIAL COORDINATE IS COMPUTED FROM SURFACE EQUATION
C  SPECIAL CASE: LEVGEO=5: SAMPLE ON SURFACE IN 3D, DISCRETISED BY
C                SET OF TRIANGLES (SELECTED SIDES OF TETRAHEDRA)
        NLSRFX=.TRUE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
        MRSURF=INSOR(NLSF,ISTRA)
        MPSURF=0
        MTSURF=0
        MASURF=0
        ILTR=0
        IF (NLTRZ) THEN
          Z0=ZZ(3)
        ELSEIF (NLTRA) THEN
          PHI=ZZ(3)*DEGRAD
        ELSEIF (NLTRT) THEN
          PHI=ZZ(3)*DEGRAD
        ENDIF
        IF (SORLIM(NLSF,ISTRA).LT.0.D0) GOTO 2000
        IF (JCALC.EQ.2.OR.JCALC.EQ.3) GOTO 993
C
        select case (LEVGEO)
        case (1)
          X0=RSURF(MRSURF)
          Y0=ZZ(2)
        case (2)
          X0=RSURF(MRSURF)*COS(ZZ(2)*DEGRAD)+EP1(MRSURF)
          Y0=RSURF(MRSURF)*SIN(ZZ(2)*DEGRAD)*ELL(MRSURF)
        case (3)
C  SAMPLING ON A RADIAL POLYGONIAL SURFACE
          BL=ZZ(2)
          DO 1501 I=1,NPPLG
            DO J=NPOINT(1,I),NPOINT(2,I)
              IF (BL.LE.BGL(MRSURF,J)) GOTO 1502
            END DO
 1501     CONTINUE
          GOTO 996
 1502     CONTINUE
          IPLG=J-1
          D=BL-BGL(MRSURF,IPLG)
          VVX=VPLX(MRSURF,IPLG)
          VVY=VPLY(MRSURF,IPLG)
          VVI=1./SQRT(VVX*VVX+VVY*VVY)
          X0=XPOL(MRSURF,IPLG)+D*VVX*VVI
          Y0=YPOL(MRSURF,IPLG)+D*VVY*VVI
        case (4)
C  SAMPLING ON SURFACE COMPOSED OF TRIANGLE SIDES
          BL=ZZ(2)
          IF (MRSURF < 0) MRSURF=ABS(MRSURF) + NLIM
          IF (ISTEP.LE.0) THEN
C  NO STEP FUNCTION,
            DO I=1,SURF_TRIAN(MRSURF)%NUMTR
              IF (BL.LE.SURF_TRIAN(MRSURF)%BGLT(I+1)) GOTO 1503
            END DO
            GOTO 996
 1503       CONTINUE
            ITRI = SURF_TRIAN(MRSURF)%ITRIAS(I)
            IPLG = SURF_TRIAN(MRSURF)%ITRISI(I)
            D = BL - SURF_TRIAN(MRSURF)%BGLT(I)
            MRSURF=ITRI
C  USE STEP FUNCTION ISTEP
          ELSE
            ITRI=IRSTEP(ISTEP,IINDEX)
            IPLG=IPSTEP(ISTEP,IINDEX)
            MRSURF=ITRI
            D=BL-RRSTEP(ISTEP,IINDEX)
          END IF
          VVX=VTRIX(IPLG,ITRI)
          VVY=VTRIY(IPLG,ITRI)
          VVI=1./SQRT(VVX*VVX+VVY*VVY)
          IF (IREVERS(IPLG,ITRI) > 0) THEN
            X0=XTRIAN(NECKE(IPLG,ITRI))+D*VVX*VVI
            Y0=YTRIAN(NECKE(IPLG,ITRI))+D*VVY*VVI
          ELSE
            IS1=IPLG+1
            IF (IS1 > 3) IS1 = 1
            X0=XTRIAN(NECKE(IS1,ITRI))-D*VVX*VVI
            Y0=YTRIAN(NECKE(IS1,ITRI))-D*VVY*VVI
          END IF
        case (5)
C  SAMPLING ON A SURFACE COMPOSED OF TETRAHEDRON SIDES
          IF (ISTEP.LE.0) GOTO 995
          ITET=IRSTEP(ISTEP,IINDEX)
          ISID=IPSTEP(ISTEP,IINDEX)
          MRSURF=ITET
          X1=XTETRA(NTECK(ITSIDE(1,ISID),ITET))
          Y1=YTETRA(NTECK(ITSIDE(1,ISID),ITET))
          Z1=ZTETRA(NTECK(ITSIDE(1,ISID),ITET))
          X2=XTETRA(NTECK(ITSIDE(2,ISID),ITET))
          Y2=YTETRA(NTECK(ITSIDE(2,ISID),ITET))
          Z2=ZTETRA(NTECK(ITSIDE(2,ISID),ITET))
          X3=XTETRA(NTECK(ITSIDE(3,ISID),ITET))
          Y3=YTETRA(NTECK(ITSIDE(3,ISID),ITET))
          Z3=ZTETRA(NTECK(ITSIDE(3,ISID),ITET))
          CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,Z0)
        case (10)
          WRITE (iunout,*) 'ERROR EXIT FROM SAMSRF. NLPOL ',LEVGEO
          CALL EIRENE_EXIT_OWN(1)
        end select
C
      ELSEIF (INDIM(NLSF,ISTRA).EQ.2) THEN
C  BIRTH POINT ON STANDARD POLOIDAL SURFACE MPSURF
C  X-OR RADIAL COORDINATE AND Z-OR TOROIDAL COORDINATE MUST BE SAMPLED
C  Y-OR POLOIDAL COORDINATE IS COMPUTED FROM SURFACE EQUATION
        NLSRFX=.FALSE.
        NLSRFY=.TRUE.
        NLSRFZ=.FALSE.
        MRSURF=0
        MPSURF=INSOR(NLSF,ISTRA)
        MTSURF=0
        MASURF=0
        ILTR=0
        IF (NLTRZ) THEN
          Z0=ZZ(3)
        ELSEIF (NLTRA) THEN
          PHI=ZZ(3)*DEGRAD
        ELSEIF (NLTRT) THEN
          PHI=ZZ(3)*DEGRAD
        ENDIF
        IF (SORLIM(NLSF,ISTRA).LT.0.D0) GOTO 2000
        IF (JCALC.EQ.1.OR.JCALC.EQ.3) GOTO 993
C
        IF (LEVGEO.EQ.1) THEN
          X0=ZZ(1)
          Y0=PSURF(MPSURF)
        ELSEIF (LEVGEO.EQ.2.AND.NLCRC) THEN
          X0=ZZ(1)*COSPH(MPSURF)
          Y0=ZZ(1)*SINPH(MPSURF)
        ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
          BL=ZZ(1)
          DO 1601 J=2,NR1ST
            IF (BL.LE.BGLP(J,MPSURF)) GOTO 1602
 1601     CONTINUE
          GOTO 996
 1602     CONTINUE
          IPLG=J-1
          D=BL-BGLP(IPLG,MPSURF)
          VVX=VVTX(IPLG,MPSURF)
          VVY=VVTY(IPLG,MPSURF)
          VVI=1./SQRT(VVX*VVX+VVY*VVY)
          X0=XPOL(IPLG,MPSURF)+D*VVX*VVI
          Y0=YPOL(IPLG,MPSURF)+D*VVY*VVI
        ELSE
          WRITE (iunout,*) 'ERROR IN SAMSRF: LEVGEO AND INDIM? '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
      ELSEIF (INDIM(NLSF,ISTRA).EQ.3) THEN
C  BIRTH POINT ON STANDARD TOROIDAL SURFACE
C  X-OR RADIAL COORDINATE AND Y-OR TOROIDAL COORDINATE MUST BE SAMPLED
C  Z-OR TOROIDAL COORDINATE IS COMPUTED FROM SURFACE EQUATION
        NLSRFX=.FALSE.
        NLSRFY=.FALSE.
        NLSRFZ=.TRUE.
        MRSURF=0
        MPSURF=0
        MTSURF=INSOR(NLSF,ISTRA)
        MASURF=0
        ILTR=0
        IF (SORLIM(NLSF,ISTRA).LT.0.D0) GOTO 2000
        IF (JCALC.EQ.1.OR.JCALC.EQ.2) GOTO 993
C
        IF (LEVGEO.EQ.1) THEN
          X0=ZZ(1)
          Y0=ZZ(2)
          Z0=ZSURF(MTSURF)
          IF (NLTRA) PHI=ZSURF(MTSURF)
        ELSEIF (LEVGEO.EQ.2.AND.NLCRC) THEN
          EP1ZZ1=0.0
          ELLZZ1=1.0
          X0=ZZ(1)*COS(ZZ(2)*DEGRAD)+EP1ZZ1
          Y0=ZZ(1)*SIN(ZZ(2)*DEGRAD)*ELLZZ1
          Z0=ZSURF(MTSURF)
          IF (NLTRA) PHI=ZSURF(MTSURF)
          IF (.NOT.NLCRC) GOTO 992
        ELSE
C  TO BE WRITTEN
          WRITE (iunout,*) 'ERROR IN SAMPLE, SOURCE ON TOR. SURFACE'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
      ELSEIF (INDIM(NLSF,ISTRA).EQ.4) THEN
C  LEVGEO=3,4:
C     BIRTH POINT ON MIX OF STANDARD RADIAL AND POLOIDAL SURFACES
C     ARC-LENGTH COORDINATE IS SAMPLED FROM STEP FUNCTION
C  LEVGEO=5:
C     TRIANGULAR SURFACE SEGMENT IS SAMPLED FROM STEP FUNCTION
C     UNIFORM DISTRIBUTION WITHIN SURFACE SEGMENT
        IF (ISTEP.LE.0) GOTO 995
        NLSRFZ=.FALSE.
        MTSURF=0
        MASURF=0
        ILTR=0
        IF (NLTRZ) THEN
          Z0=ZZ(3)
        ELSEIF (NLTRA) THEN
          PHI=ZZ(3)*DEGRAD
        ELSEIF (NLTRT) THEN
          PHI=ZZ(3)*DEGRAD
        ENDIF
        IF (SORLIM(NLSF,ISTRA).LT.0.D0) GOTO 2000
        IF (JCALC.EQ.1.OR.JCALC.EQ.3) GOTO 993
C
        BL=ZZ(1)
        DO 1701 J=2,NSMAX(ISTEP)
          IF (BL.LE.RRSTEP(ISTEP,J)) GOTO 1702
 1701   CONTINUE
        GOTO 996
 1702   CONTINUE
        IPLG=J-1
        D=BL-RRSTEP(ISTEP,IPLG)
        IF (LEVGEO.EQ.4) THEN
          IF (ISTEP.LE.0) GOTO 995
          NLSRFX=.TRUE.
          NLSRFY=.FALSE.
          ITRI=IRSTEP(ISTEP,IINDEX)
          ISID=IPSTEP(ISTEP,IINDEX)
          MRSURF=ITRI
          VVX=VTRIX(ISID,ITRI)
          VVY=VTRIY(ISID,ITRI)
          VVI=1./SQRT(VVX*VVX+VVY*VVY)
          IF (IREVERS(ISID,ITRI) > 0) THEN
            X0=XTRIAN(NECKE(ISID,ITRI))+D*VVX*VVI
            Y0=YTRIAN(NECKE(ISID,ITRI))+D*VVY*VVI
          ELSE
            IS1=ISID+1
            IF (IS1 > 3) IS1 = 1
            X0=XTRIAN(NECKE(IS1,ITRI))-D*VVX*VVI
            Y0=YTRIAN(NECKE(IS1,ITRI))-D*VVY*VVI
          END IF
        ELSEIF (LEVGEO.EQ.5) THEN
          IF (ISTEP.LE.0) GOTO 995
          NLSRFX=.TRUE.
          NLSRFY=.FALSE.
          ITET=IRSTEP(ISTEP,IINDEX)
          ISID=IPSTEP(ISTEP,IINDEX)
          MRSURF=ITET
          X1=XTETRA(NTECK(ITSIDE(1,ISID),ITET))
          Y1=YTETRA(NTECK(ITSIDE(1,ISID),ITET))
          Z1=ZTETRA(NTECK(ITSIDE(1,ISID),ITET))
          X2=XTETRA(NTECK(ITSIDE(2,ISID),ITET))
          Y2=YTETRA(NTECK(ITSIDE(2,ISID),ITET))
          Z2=ZTETRA(NTECK(ITSIDE(2,ISID),ITET))
          X3=XTETRA(NTECK(ITSIDE(3,ISID),ITET))
          Y3=YTETRA(NTECK(ITSIDE(3,ISID),ITET))
          Z3=ZTETRA(NTECK(ITSIDE(3,ISID),ITET))
          CALL EIRENE_FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,Z0)
        ELSEIF (LEVGEO.EQ.3) THEN
          IF (IGSTEP(ISTEP,IPLG).LT.200000) THEN
C  RADIAL PART
            NLSRFX=.TRUE.
            NLSRFY=.FALSE.
            MPSURF=0
            MRSURF=IGSTEP(ISTEP,IPLG)-100000
            VVX=VPLX(MRSURF,IPSTEP(ISTEP,IPLG))
            VVY=VPLY(MRSURF,IPSTEP(ISTEP,IPLG))
            VVI=1./SQRT(VVX*VVX+VVY*VVY)
            X0=XPOL(MRSURF,IPSTEP(ISTEP,IPLG))+D*VVX*VVI
            Y0=YPOL(MRSURF,IPSTEP(ISTEP,IPLG))+D*VVY*VVI
          ELSE
C  POLOIDAL PART
            NLSRFX=.FALSE.
            NLSRFY=.TRUE.
            MRSURF=0
            MPSURF=IGSTEP(ISTEP,IPLG)-200000
            VVX=VVTX(IRSTEP(ISTEP,IPLG),MPSURF)
            VVY=VVTY(IRSTEP(ISTEP,IPLG),MPSURF)
            VVI=1./SQRT(VVX*VVX+VVY*VVY)
            X0=XPOL(IRSTEP(ISTEP,IPLG),MPSURF)+D*VVX*VVI
            Y0=YPOL(IRSTEP(ISTEP,IPLG),MPSURF)+D*VVY*VVI
          ENDIF
        ELSE
          GOTO 992
        ENDIF
      ENDIF
C
 2000 CONTINUE
C
C  X0,Y0,Z0 AND PHI ARE GIVEN NOW.
C  MAKE SURE, THAT THESE ARE IN PROPER PERIODICITY BLOCK
C
      IF (NLTRA.AND.ILTR.GT.0) THEN
C  Z0 IS A CARTESIAN COORDINATE IN LOCAL SYSTEM ILTR
C  BIRTHPOINT IS ON ADDITIONAL SURFACE (BECAUSE ILTR.NE.0)
C  FIND PHI,NT AND THEN X0,Z0 IN CELL NT
        CALL EIRENE_FZRTOR(X0,Z0,ILTR,XR,PHI,NT,.FALSE.,0)
        IF (NT.NE.ILTR)
     .    CALL EIRENE_FZRTRI(X0,Z0,NT,XR,PHI,NT)
        X01=X0+RMTOR
      ELSEIF (NLTRA.AND.ILTR.EQ.0) THEN
C  PHI IS THE TOROIDAL ANGLE (RADIANS)
C  BIRTHPOINT IS EITHER ON ADDITIONAL OR ON STANDARD GRID SURFACE
C  FIND Z0, NT,  FROM X0,PHI
C
C  MOVE PHI AWAY FROM TOROIDAL PERIODICITY SURFACE
        IF (PHI.EQ.0.0.AND.ABS(ZSURF(1)).LT.EPS10) PHI=0.01_DP
C  DONE
        IF (PHI.LT.ZSURF(1)) PHI=PHI+PI2A
        IF (PHI.GT.ZSURF(NTTRA)) PHI=PHI-PI2A
        NT=EIRENE_LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMSRF      ')
        IF (NT.LE.0.OR.NT.GT.NTTRAM) THEN
          WRITE (iunout,*) 'NT OUT OF RANGE IN SAMSRF '
          WRITE (iunout,*) PHI,ZFULL
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
        X01=X0+RMTOR
        CALL EIRENE_FZRTRI(X0,Z0,NT,X01,PHI,NT)
      ENDIF
C
C  FIND ZONE NUMBER NRCELL, NACELL ,NBLOCK AND IPOLG
C
      IF (NASOR(NLSF,ISTRA).GT.0.AND.NRSOR(NLSF,ISTRA).EQ.0) THEN
C  NRCELL=0
C  IPOLG=0
C  NBLOCK=NBMLTP
C  NACELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NASOR
        NACELL=NASOR(NLSF,ISTRA)
        NRCELL=0
        NBLOCK=NBMLTP
        IPOLG=0
      ELSEIF (NRSOR(NLSF,ISTRA).GT.0.AND.NASOR(NLSF,ISTRA).EQ.0) THEN
C  NRCELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NRSOR
C  IPOLG  IS EXPLICITLY DEFINED BY INPUT VARIABLE NISOR
C  NBLOCK IS EXPLICITLY DEFINED BY INPUT VARIABLE NBSOR
C  NACELL=0
        NRCELL=NRSOR(NLSF,ISTRA)
        NBLOCK=NBSOR(NLSF,ISTRA)
        IPOLG=NISOR(NLSF,ISTRA)
        NACELL=0
      ELSEIF (NRSOR(NLSF,ISTRA).EQ.0.AND.NASOR(NLSF,ISTRA).EQ.0) THEN
C  NRCELL IS COMPUTED IN STANDARD MESH
C  IPOLG IS COMPUTED IN STANDARD MESH
C  NBLOCK IS EXPLICITLY DEFINED BY INPUT VARIABLE NBSOR
C  NACELL=0
        IF (NLSRFX) THEN
          NRCELL=MIN0(NR1STM,MRSURF)
          IAN=MRSURF
          IEN=MRSURF
          NDUM=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,IAN,IEN,NLSRFX,NLSRFY,NPANU,
     .        'SAMSRF      ')
        ELSEIF (NLSRFY) THEN
          IPOLG=MIN0(NP2NDM,MPSURF)
          IAN=MPSURF
          IEN=MPSURF
          NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IDUM,IAN,IEN,NLSRFX,NLSRFY,
     .                  NPANU,'SAMSRF      ')
        ELSE
          NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,NLSRFX,NLSRFY,
     .                  NPANU,'SAMSRF      ')
        ENDIF
        NACELL=0
        NBLOCK=NBSOR(NLSF,ISTRA)
      ELSEIF (NRSOR(NLSF,ISTRA).LT.0.AND.ISTEP.GT.0) THEN
C  NRCELL IS FOUND FROM STEP FUNCTION DATA
C  IPOLG IS FOUND FROM STEP FUNCTION DATA
C  NBLOCK IS FOUND FROM STEP FUNCTION DATA
C  NACELL IS FOUND FROM STEP FUNCTION DATA
        NRCELL=IRSTEP(ISTEP,IINDEX)
        IPOLG =IPSTEP(ISTEP,IINDEX)
        NBLOCK=IBSTEP(ISTEP,IINDEX)
        NACELL=IASTEP(ISTEP,IINDEX)
      ELSEIF (NRSOR(NLSF,ISTRA).LT.0.AND.ISTEP.LT.0) THEN
C  NRCELL IS TRANSFERRED FROM SUBR. SAMUSR
C  IPOLG IS TRANSFERRED FROM SUBR. SAMUSR
C  NBLOCK IS TRANSFERRED FROM SUBR. SAMUSR
C  NACELL IS TRANSFERRED FROM SUBR. SAMUSR
        NRCELL=IRUSR
        NACELL=IAUSR
        NBLOCK=IBUSR
        IPOLG =IPUSR
      ELSE
        GOTO 999
      ENDIF
C
      NTCELL=1
      IPERID=1
      IF (NLTOR.AND.NACELL.EQ.0) THEN
        IF (NLTRZ) THEN
          IF (NTSOR(NLSF,ISTRA).GT.0) THEN
C  NTCELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NTSOR
            NTCELL=NTSOR(NLSF,ISTRA)
          ELSEIF (NTSOR(NLSF,ISTRA).LT.0.AND.ISTEP.GT.0) THEN
C  NTCELL IS FOUND FROM STEP FUNCTION DATA
            NTCELL=ITSTEP(ISTEP,IINDEX)
          ELSEIF (NTSOR(NLSF,ISTRA).LT.0.AND.ISTEP.LT.0) THEN
C  NTCELL IS TRANSFERRED FROM SUBR. SAMUSR
            NTCELL=ITUSR
          ELSEIF (NTSOR(NLSF,ISTRA).EQ.0) THEN
C  NTCELL IS COMPUTED IN STANDARD MESH
            NTCELL=EIRENE_LEARCA(Z0,ZSURF,1,NT3RD,1,'SAMSRF      ')
          ELSE
            GOTO 999
          ENDIF
!pb
          NNTCLL = NTCELL
        ELSEIF (NLTRA) THEN
C  NTSOR NOT AVAILABLE FOR NLTRA OPTION
          NTCELL=NT
          IPERID=NTCELL
!pb
          NNTCLL = NTCELL
        ELSEIF (NLTRT) THEN
          WRITE (iunout,*) 'NLTRT: TO BE WRITTEN IN SAMSRF '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ELSEIF (.NOT.NLTOR.OR.NACELL.GT.0) THEN
C       IF (NLTRZ) THEN
C  ALL DONE
        IF (NLTRA) THEN
C  NTSOR NOT AVAILABLE FOR NLTRA OPTION
          IPERID=NT
!pb
          NNTCLL = NTCELL
        ELSEIF (NLTRT) THEN
          WRITE (iunout,*) 'NLTRT: TO BE WRITTEN IN SAMSRF '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ENDIF
C
      NPCELL=1
      IF (NLPOL.AND.NACELL.EQ.0) THEN
        IF (NPSOR(NLSF,ISTRA).GT.0) THEN
C  NPCELL IS EXPLICITLY DEFINED BY INPUT VARIABLE NPSOR
          NPCELL=NPSOR(NLSF,ISTRA)
        ELSEIF (NPSOR(NLSF,ISTRA).LT.0.AND.ISTEP.GT.0) THEN
C  NPCELL IS FOUND FROM STEP FUNCTION DATA
          IF ((LEVGEO == 4) .OR. (LEVGEO == 5)) THEN
            NPCELL = 1
          ELSE
            NPCELL=IPSTEP(ISTEP,IINDEX)
          END IF
        ELSEIF (NPSOR(NLSF,ISTRA).LT.0.AND.ISTEP.LT.0) THEN
C  NPCELL IS TRANSFERRED FROM SUBR. SAMUSR
          NPCELL=IPUSR
        ELSEIF (NPSOR(NLSF,ISTRA).EQ.0) THEN
C  NPCELL IS COMPUTED IN STANDARD MESH
          IF (LEVGEO.EQ.1) THEN
            NPCELL=EIRENE_LEARCA(Y0,PSURF,1,NP2ND,1,'SAMSRF')
          ELSEIF (LEVGEO.EQ.2.AND.NLCRC) THEN
            WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
            NPCELL=EIRENE_LEARCA(WINK,PSURF,1,NP2ND,1,'SAMSRF')
          ELSEIF (LEVGEO.EQ.2.OR.LEVGEO.EQ.3) THEN
            IF (NLSRFY) THEN
              NPCELL=MIN0(NP2NDM,MPSURF)
            ELSE
              NPCELL=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'SAMSRF')
            ENDIF
          ELSE
            WRITE (iunout,*) 'ERROR EXIT FROM SAMSRF. NLPOL ',LEVGEO
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
        ELSE
          GOTO 999
        ENDIF
      ENDIF
C
      NBLOCK=MAX0(1,NBLOCK)
      NBLOCK=MIN0(NBLOCK,NBMLT)
      IF (NRCELL.GT.0) NACELL=0
      IF (NACELL.GT.0) NBLOCK=NBMLTP
      NBLCKA=NSTRD*(NBLOCK-1)+NACELL
      NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C
      IF (ISTEP.GT.0) THEN
C  TAKE BACKGROUND MEDIUM DATA AT PLACE OF BIRTH FROM STEP FUNCTION ISTEP
        TEWL=TESTEP(ISTEP,IINDEX)
        SHWL=SHSTEP(ISTEP,IINDEX)
        DO 3010 JPLS=1,MAX(NPLSI,NATMI,NMOLI)
          IF (NLPLS(ISTRA)) THEN
            IF (JPLS.GT.NPLSI) CYCLE
            IPLSTI=MPLSTI(JPLS)
            IPLSV=MPLSV(JPLS)
          ELSE IF (NLATM(ISTRA)) THEN
            IF (JPLS.GT.NATMI.AND.L_MACRO.LT.EPS12) CYCLE
            IPLSTI=JPLS
            IPLSV=JPLS
          ELSE IF (NLMOL(ISTRA)) THEN
            IF (JPLS.GT.NMOLI) CYCLE
            IPLSTI=JPLS
            IPLSV=JPLS
          END IF
          TIWL(JPLS)=TISTEP(IPLSTI,ISTEP,IINDEX)
          VXWL(JPLS)=VXSTEP(IPLSV,ISTEP,IINDEX)
          VYWL(JPLS)=VYSTEP(IPLSV,ISTEP,IINDEX)
          VZWL(JPLS)=VZSTEP(IPLSV,ISTEP,IINDEX)
          DIWL(JPLS)=DISTEP(JPLS,ISTEP,IINDEX)
                  FL=FLSTEP(JPLS,ISTEP,IINDEX)
          EFWL(JPLS)=ELSTEP(JPLS,ISTEP,IINDEX)/(FL+EPS30)
          ZIWL(JPLS)=ZISTEP(JPLS,ISTEP,IINDEX)
 3010   CONTINUE
      ELSEIF (ISTEP.EQ.0) THEN
C  TAKE BACKGROUND MEDIUM DATA AT PLACE OF BIRTH FROM LOCAL BULK PLASMA DATA
C                              IN SAMPLED CELL NCELL
        TEWL=TEIN(NCELL)
        SHWL=0.
        DO 3020 JPLS=1,NPLSI
          IPLS=JPLS
          IPLSTI = MPLSTI(IPLS)
          IPLSV = MPLSV(IPLS)
          TIWL(IPLS)=TIIN(IPLSTI,NCELL)
          IF (INDPRO(4) == 8) THEN
cdr  no spatial information x0,y0,z0 of birth point available here for vecusr.
c    set drift velocities at cell center
            XC=0.
            YC=0.
            ZC=0.
            CALL EIRENE_VECUSR (2,NCELL,XC,YC,ZC,VX,VY,VZ,IPLS,
     .                          .FALSE.)
            VXWL(IPLS)=VX
            VYWL(IPLS)=VY
            VZWL(IPLS)=VZ
          ELSE
            VXWL(IPLS)=VXIN(IPLSV,NCELL)
            VYWL(IPLS)=VYIN(IPLSV,NCELL)
            VZWL(IPLS)=VZIN(IPLSV,NCELL)
          END IF
          DIWL(IPLS)=DIIN(IPLS,NCELL)
          EFWL(IPLS)=0.
 3020   CONTINUE
      ELSEIF (ISTEP.LT.0) THEN
C  TEWL, TIWL, .... SHWL ALREADY DEFINED IN SAMUSR
C  NOTHING MORE TO BE DONE HERE
      ENDIF
      IF(ISTEP.LE.0 .OR. sum(ZIWL(1:NPLSI))==0.0) THEN
cnh   28.11.2019
      DO JPLS=1,NPLSI
        IF (ZIIN(JPLS,NCELL).NE.ZVAC) THEN
          ZIWL(JPLS) = ZIIN(JPLS,NCELL)
        ELSE
          ZIWL(JPLS) = DBLE(NCHRGP(JPLS))
        ENDIF
      ENDDO
      ENDIF
C
C  SET ANALOG SPECIES INDEX DISTRIBUTION WEISPZ
C
      IF (NSPEZ(ISTRA).LE.0) THEN
        IF (ISTEP.GT.0) THEN
C  WEISPZ FROM STEP FUNCTION
          DO 4100 JSPZ=NSPSTI(ISTEP),NSPSTE(ISTEP)
            WEISPZ(JSPZ)=FLSTEP(JSPZ,ISTEP,IINDEX)/
     .                   FLSTEP(0,   ISTEP,IINDEX)
 4100     CONTINUE
          DO JSPZ=1,NSPSTI(ISTEP)-1
            WEISPZ(JSPZ)=0.D0
          ENDDO
          DO JSPZ=NSPSTE(ISTEP)+1,NSPZ
            WEISPZ(JSPZ)=0.D0
          ENDDO
        ELSE
C  WEISPZ NOT DEFINED
        ENDIF
      ENDIF
C
      RETURN
  992 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF '
      WRITE (iunout,*) 'INDIM, LEVGEO: ',INDIM(NLSF,ISTRA),LEVGEO
      WRITE (iunout,*) 'THIS CASE IS TO BE WRITTEN '
      CALL EIRENE_EXIT_OWN(1)
  993 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      WRITE (iunout,*) 'SOURCE ON STANDARD SURFACE BUT INVALID "JCALC"'
      WRITE (iunout,*) 'CHANGE INPUT FLAG "INDTEC"'
      CALL EIRENE_EXIT_OWN(1)
 9931 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      WRITE (iunout,*)
     .  'SOURCE ON ADDITIONAL SURFACE BUT INVALID "JCALC"'
      WRITE (iunout,*) 'CHANGE INPUT FLAG "INDTEC"'
      CALL EIRENE_EXIT_OWN(1)
  995 CONTINUE
      WRITE (iunout,*) ' ERROR IN SAMSRF '
      WRITE (iunout,*) 'INDIM, LEVGEO: ',INDIM(NLSF,ISTRA),LEVGEO
      WRITE (iunout,*)
     .  ' AND NO STEP FUNCTION FOR COMPUTING THE BIRTH PLACE'
      WRITE (iunout,*) ' THIS CASE IS TO BE WRITTEN '
      CALL EIRENE_EXIT_OWN(1)
  996 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      WRITE (iunout,*)
     . 'ARCLENGTH BL GREATER THAN THE LENGTH OF THE POLYGON'
      WRITE (iunout,*) 'BL = ',BL
      CALL EIRENE_EXIT_OWN(1)
  997 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      WRITE (iunout,*) 'CHECK INPUT OF SORLIM, BECAUSE EIRENE TRIES TO'
      WRITE (iunout,*)
     .  'COMPUTE MORE THAN ONE VARIABLE FROM THE SURFACE'
      WRITE (iunout,*) 'EQUATION'
      CALL EIRENE_EXIT_OWN(1)
  998 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      CALL EIRENE_MASR3('X0, Y0, Z0              ',X0,Y0,Z0)
      WRITE (iunout,*) 'BUT OUTSIDE THE VALID AREA OF THIS SURFACE'
      CALL EIRENE_EXIT_OWN(1)
  999 CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMSRF. EXIT CALLED'
      WRITE (iunout,*) 'INITIAL CELL NUMBER INVALID OR DET NEGATIVE'
      CALL EIRENE_MASR4('X0, Y0, Z0, DET                 ',X0,Y0,Z0,DET)
      WRITE (iunout,*) 'ISTEP ',ISTEP
      CALL EIRENE_MASJ3('NBLOCK,NACELL,NRCELL    ',NBLOCK,NACELL,NRCELL)
      CALL EIRENE_MASJ3('NPCELL,NTCELL,IPOLG     ',NPCELL,NTCELL,IPOLG)
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_SAMSF1

      SUBROUTINE EIRENE_SAMSF2
      IMPLICIT NONE

      IF (ALLOCATED(INDTEC)) THEN
        DEALLOCATE (INDTEC)
        DEALLOCATE (ALEFT)
        DEALLOCATE (BRGHT)
        DEALLOCATE (XI)
        DEALLOCATE (XE)
      END IF

      RETURN
      END SUBROUTINE EIRENE_SAMSF2

      END MODULE EIRMOD_SAMSRF
