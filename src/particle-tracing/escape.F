cdr march 18  : bug fix re semi-transp. surfaces.
cdr             This intermediate bug was introduced in jan 18 commit
cdr jan. 18   : outgoing flux tallies scored in eirene_update_surface(ind=1).
cdr             Semi-transp fluxes: score only incident and emitted current fractions
cdr             for which surfaces are NOT transparent.
cdr             update_sptflx: different meaning of flag IND. More consistent
cdr             now with IND-flag in other surface scoring routines.

cdr nov. 17   : lmetspw arguments corrected
cdr sept.17   : no ion sheath orbit correction at mirror surfaces (=symmetry BC)
cdr aug.17    : bug fix. cond exp. estimator, on purely absorbing surface.
c               return 3, if icol=1, even for purely absorbing surfaces.
c               plus some minor clean-up, commenting.
c 06.08.15    : arguments added to vecusr
c   aug.15    : periodicity and icol=1, return iret=3 rather than return iret=2
C   OCT.14    : ARGUMENTS IN VELOCS: WEIGHT AND VWL
C   OCT.14    : SPUTTERING: SCORE FLUXES ALSO IN CASE SPUTTERED PARTICLES ARE NOT FOLLOWED
c               new meaning of isrs, isrc=0: sputter, score fluxes, but do not follow.
C               IGASP_OLD = 0 OPTION. INTRODUCED: ITYP_OLD, IGASP_OLD, IGASC_OLD

CDR 21.04.2010: ILIIN=2: DO NOT SPUTTER (CRTX,... ETC. IS NOT DEFINED).
!   20.06.07: maximum number of stored splitting particle MAXLEVEL introduced
CDR 07.12.06: comments, correction for ILIIN=-3 option:
c             old: prf... and erf... tallies have also been filled in
c                  case of iliin=-3 and "negative" incident direction, due to
C                  typing error. but they have not been printed, so no erroneous
c                  output was produced from outflx.f
c                  The correct "net" fluxes are on POT.... and EOT.... tallies.
c             now: do not fill prf... and erf... tallies at all in case ILIIN=-3
!pb 08.11.06: definition of splitting arrays changed
!             RSPLST(NLEVEL,1:NPARTC) --> RSPLST(1:NPARTC,NLEVEL)
!             ISPLST(NLEVEL,1:MPARTC) --> ISPLST(1:MPARTC,NLEVEL)
!pb 21.08.06: calls to update_surface introduced
c   20.01.06: sheath repulsion added for negative ions (or
c             positive ions and positive sheath potentials).
C
      SUBROUTINE EIRENE_ESCAPE(PR,SG,IRET)
C
C  PROCESS ESCAPING PARTICLES
C  INPUT:
C        PR: Probability to reach surface, times: sign relative to surf. normal
C            i.e.: WPR=WEIGHT*PR is positive or negative incident flux of test
C                  particle, sign depends upon surface orientation.
C        SG: sign relative to surf. normal. SG=1.0 OR SG=-1.0
C        LGPART=.TRUE.  UPDATE TALLIES FOR INCIDENT PARTICLES,
C                       THEN CALL SURFACE MODEL (SPUTER, REFLEC,...)
C                       THEN UPDATE TALLIES FOR EMITTED PARTICLES
C        LGPART=.FALSE. UPDATE TALLIES FOR INCIDENT PARTICLES, THEN STOP.
C        ICOL:  =1:  CONDITIONAL EXPECTATION ESTIMATOR, AND AN EARLIER COLLISION
C                    ALONG THIS TRACK IS STORED
cdr    IND flag used internally:  IND=1 score on OT tallies (going "onto" surface)
cdr                               IND=2 score on RF tallies (coming "from" surface)
C
C  RETURN  : STOP TRACK OF THIS PARTICLE TYPE. RETURN TO SUBR. MCARLO
C  RETURN IRET=1: START NEW TRACK OF SAME TYPE IN CALLING PROGRAM
C                 (I.E. IN FOLNEUT OR FOLION)
C  RETURN IRET=2: CONTINUE THIS TRACK IN CALLING PROGRAM (TRANSP. SURFACE)
C  RETURN IRET=3: RESTORE COLLISION DATA, CONDITIONAL EXPECTATION WAS USED
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CINIT
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSPL
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_CTRIG
      USE EIRMOD_CUPD
      USE EIRMOD_CFPLK
      USE EIRMOD_CSDVI
      USE EIRMOD_RANF, ONLY: RANF_EIRENE
      USE EIRMOD_SWITCH_PARTINFO, ONLY: EIRENE_SWITCH_PARTINFO
      USE EIRMOD_SPUTER, ONLY: EIRENE_SPUTR1
      USE EIRMOD_REFLEC, ONLY: EIRENE_REFLC1
      USE EIRMOD_REFUSR, ONLY: EIRENE_RF2USR
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC
      USE EIRMOD_SHEATH, ONLY: EIRENE_SHEATH

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: PR, SG
      INTEGER, INTENT(OUT) :: IRET
      REAL(DP) :: DIWL(NPLS), VPWL(NPLS), ZIWL(NPLS)
      REAL(DP) :: VXSPTC, VYSPTC, VZSPTC, VSPTC, ESPTC, VXSPTP,
     .          VYSPTP, VZSPTP, ESPTP, VELXS, VELYS, VELZS, VSPTP, TW,
     .          E0TERM, FR2, COSI2, ZVZ, WABS,
     .          CUR, GAMMA, TEWL, VX, VY, VZ, FCHAR, WPR, FMASS,
     .          FLX, YIELD1, YIELD2, VELS, WEIGHS, E0S, ESHET,EVCQ,
     .          VSHETQ, V, VELSH, VC, VCQ, VC2,
     .          SPLFLG,
     .          VXR, VYR, VZR, VWL, WGHTVS, RATR, YSPTWL
      INTEGER :: ISG, ISPZS, I, J, IDIM, MS, IC, IP, ISTS,
     .           ISSPTP, ISSPTC, IPV, MODREF, MOL_DEFAULT,
     .           ITYP_OLD, IGASP_OLD, IGASC_OLD, IOLD, IND
      LOGICAL :: NLSPUT, LTRANS
      EXTERNAL :: EIRENE_FCOSIN, EIRENE_REFLC1_PHOTON, EIRENE_ROTATF,
     .            EIRENE_UPDATE_SPECTRUM,
     .            EIRENE_UPDATE_SPTFLX, EIRENE_UPDATE_SURFACE,
     .            EIRENE_UPSUSR, EIRENE_VECUSR, EIRENE_VELOCS,
     .            EIRENE_LEER, EIRENE_EXIT_OWN
C
      IRET = 0

      IF (IMETWL(MSURF) == 0) THEN
        NWLMT = NWLMT+1
        IWLMT(NWLMT) = MSURF
        IMETWL(MSURF) = NWLMT
      END IF

C
C  .................
C  .               .
C  .  PERIODICITY: .
C  .................
C
C  CURRENTLY: ONLY IN CASE (LEVGEO=1, NLTRZ). HENCE: VEL_OLD=VEL_NEW
C  TO BE WRITTEN: TOROIDICITY AS SPECIAL CASE OF PERIODICITY
C  ALL FIDDLING WITH VELOCITIES AT PERIODIC SURFACES MUST HAVE BEEN DONE
C  ALREADY IN CALLING PROGRAMS, E.G. STDNOR.F
C  CURRENTLY: NO SURFACE TALLIES AT PERIODICITY SURFACES

      IF (ILIIN(MSURF).GE.4) THEN
cdr: unfinished option: store trajectories for later postprocessing
cdr: unused
c       NLTRJ = .FALSE.
c       TRAJ(ITRJ)%TRJ%NO_SURF = MSURF
C  CONDITIONAL EXPECTATION ESTIMATOR: HAS THIS PARTICLE COLLIDED IN THE VOLUME,
C  BEFORE IT HIT THE WALL?
        IF (ICOL.EQ.1) then
          colflag = .true.  ! probably unused, perhaps in TIM?
          IRET = 3
          RETURN
        ENDIF
C  NO, PARTICLE HAS ARRIVED AT PERIODICITY SURFACE MSURF
        IF (.NOT.LGPART) THEN
          WRITE (IUNOUT,*) 'ERROR AT PERIODICITY SURFACE, LGPART=FALSE '
          IRET = 0
          RETURN
        ENDIF
        IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
          CALL EIRENE_CHCTRC(X0,Y0,Z0,0,11)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        ENDIF
        IRET = 1
        RETURN
      ENDIF
C
C  .............................
C  .                           .
C  .  UPDATE INCIDENT FLUXES:  .
C  .............................
C
      ITYP_OLD=ITYP
      WPR=WEIGHT*PR

      select case (ITYP_OLD)
      case (0)
        IOLD = IPHOT
      case (1)
        IOLD = IATM
      case (2)
        IOLD = IMOL
      case (3)
        IOLD = IION
      case (4)
        IOLD = IPLS
      case default
        IOLD = 0
      end select
C
C  UPDATE PARTICLE OUTGOING FLUX ONTO SURFACE MSURF
C  UPDATE ENERGY OUTGOING FLUX ONTO SURFACE MSURF
C
C  SPATIAL RESOLUTION ON NON-DEFAULT STANDARD SURFACE?
      IF (MSURF.GT.NLIM.AND.NLMPGS.GT.NLIMPS) THEN
        select case (LEVGEO)
        case (:3)
          ISTS=MSURF-NLIM
          MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
          FLX=FLXOUT(MSURFG)
        case (4)
          MSURFG=NLIM+NSTS+INSPAT(IPOLGN,MRSURF)
          FLX=FLXOUT(MSURFG)
        case (5)
cdr  to be written
c         MSURFG=NLIM+NSTS+INSPAT(IPOLGN,MRSURF)
c         FLX=FLXOUT(MSURFG)
        case default
          MSURFG=0
          FLX=FLXOUT(MSURF)
        end select
      ELSEIF (MSURF.GT.0) THEN
        MSURFG=0
        FLX=FLXOUT(MSURF)
      ELSE
        WRITE (iunout,*) 'ERROR IN ESCAPE: MSURF=0. KILL PARTICLE '
        IRET = 0
        RETURN
      ENDIF
C
C  FLUXES ARE UPDATED HERE IN THE FOLLOWING CASES:
C
C  FOR NON-TRANSPARENT SURFACES: SCORE INCIDENT FLUX
C  FOR     TRANSPARENT SURFACES: SCORE ONE-SIDED FLUX, ONLY POSITIVE COMPONENT
C  FOR     TRANSPARENT SURFACES: SCORE NET FLUX IN CASE OF ILIIN(MSURF)=-3
C
      IF ((ILIIN(MSURF).LT.0).AND.(SG.LT.0.D0).AND.(ILIIN(MSURF).NE.-3))
     .GOTO 10
C
C   HERE: EITHER: ILIIN .GE.0,    NOT TRANSPARENT, SCORE OUTGOING FLUX (WPR >=0 always)
c             OR: ILIIN.EQ.-3,    TRANSPARENT,     SCORE NET FLUX (WPR contains sign)
C             OR: SG .GT.0        SCORE ONE-SIDED POSITIVE CURRENT ONLY, EVEN FOR TRANSPARENT SURF.
c
C                                         NEGATIVE CURRENTS ARE TAKEN CARE OF FURTHER BELOW

      IF (ITYP.EQ.0) THEN
C  INCIDENT PHOTONS
        FMASS=0._dp
        FCHAR=0._dp
      ELSEIF (ITYP.EQ.1) THEN
C  INCIDENT ATOMS
        FMASS=DBLE(NMASSA(IATM))
        FCHAR=DBLE(NCHARA(IATM))
      ELSEIF (ITYP.EQ.2) THEN
C  INCIDENT MOLECULES
        FMASS=DBLE(NMASSM(IMOL))
        FCHAR=DBLE(NCHARM(IMOL))
      ELSEIF (ITYP.EQ.3) THEN

C  SPECIAL CASE: INCIDENT TEST IONS. SHEATH ACCELERATION.
!pb ispz calculated for check of semitransparency
cdr: sheath ion orbit part: make separate routine. tbd: E.g. Lindner-formula.
cdr: iliin=3:  mirror     : no sheath acceleration
cdr: iliin>=4: periodicity: no sheath acceleration
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
        IF (
     .      (ILIIN(MSURF).GT.0)                                    .AND.
     .    (abs(transp(ispz,1,msurf))+abs(transp(ispz,2,msurf))==0) .AND.
     .      (ILIIN(MSURF).LT.3)                  ) then
C  SURFACE MSURF IS NOT MADE (PARTIALLY) TRANSPARENT FOR TEST ION SPECIES IION (ISPZ)
          ESHET=0.D0
C  ACCOUNT FOR ELECTROSTATIC SHEATH AT SURFACE FOR TEST IONS
          IF (FSHEAT(MSURF).LE.0.D0) THEN
            GAMMA=0.
            CUR=0.
            IC=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
            IF (.NOT.LGVAC(IC,NPLS+1)) THEN
              TEWL=TEIN(IC)
              DO 30 IP=1,NPLSI
                IPV=MPLSV(IP)
                IF (INDPRO(4) == 8) THEN
                  CALL EIRENE_VECUSR (2,IC,X0,Y0,Z0,VX,VY,VZ,IP,
     .                                .TRUE.)
                ELSE
                  VX=VXIN(IPV,IC)
                  VY=VYIN(IPV,IC)
                  VZ=VZIN(IPV,IC)
                ENDIF
                VPWL(IP)=SQRT(VX**2+VY**2+VZ**2)
                DIWL(IP)=DIIN(IP,IC)
cnh      28.10.2019
                IF (ZIIN(IP,IC).NE.ZVAC) THEN
                  ZIWL(IP) = ZIIN(IP,IC)
                ELSE
                  ZIWL(IP) = DBLE(NCHRGP(IP))
                ENDIF
   30         CONTINUE
cnh      02.11.2019   NCHRGP --> ZIWL
              ESHET=NCHRGI(IION)*EIRENE_SHEATH(TEWL,DIWL,VPWL,
     .                                  ZIWL,GAMMA,CUR,NPLSI,MSURF)
            ENDIF
          ELSE
            IC=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
            IF (.NOT.LGVAC(IC,NPLS+1)) THEN
              TEWL=TEIN(IC)
              ESHET=FSHEAT(MSURF)*TEWL*NCHRGI(IION)
            ENDIF
          ENDIF
C  SIMPLE SHEATH MODEL:
C    SHEATH IS NOT RESOLVED SPATIALLY.
C    CASE A: ACCELERATING SHEATH (ESHET.GT.0):
C      THE VELOCITY COMPONENT
C      NORMAL TO THE SURFACE IS ENHANCED BY A "SHEATH-VELOCITY"
C      VELSH (CM/S), SUCH THAT THE ENERGY INCREASES BY ESHET(EV)
C    CASE B: REPULSIVE SHEATH (ESHET.LT.0):
C      CASE B1: PARTICLE STILL OVERCOMES SHEATH, BUT WITH REDUCED ENERGY
C      CASE B2: PARTICLE IS REFLECTED BY SHEATH
C
C  VSHETQ SHEATH ENERGY IN UNITS OF CM^2/S^2
C     NEGATIVE FOR REPULSIVE POTENTIAL
C     POSITIVE FOR ATTRACTIVE POTENTIAL
          VSHETQ=ESHET*RSQDVI(IION)*RSQDVI(IION)
C  VC = NORMAL COMPONENT OF VELOCITY VECTOR:
          VC=VEL*(VELX*CRTX+VELY*CRTY+VELZ*CRTZ)
          VCQ=VC*VC
C
cdr       eshet=0.  !  no sheath acceleration,  for testing only
c
          IF (ESHET.GT.0.D0) THEN
C  POSITIVELY CHARGED  IONS, ATTRACTIVE SHEATH POTENTIAL
C  ADD VELOCITY DUE TO SHEATH ACCELERATION
            VELSH=-VC+SQRT(VCQ+VSHETQ)
            VX=VEL*VELX+VELSH*CRTX
            VY=VEL*VELY+VELSH*CRTY
            VZ=VEL*VELZ+VELSH*CRTZ
            V=SQRT(VX*VX+VY*VY+VZ*VZ)
            VELX=VX/V
            VELY=VY/V
            VELZ=VZ/V
            E0=E0+ESHET
            VEL=SQRT(E0)*RSQDVI(IION)
            EELFI(IION,ISTRA)=EELFI(IION,ISTRA)+ESHET*WPR
          ELSEIF (ESHET.LT.0.D0) THEN
C  NEGATIVELY CHARGED  IONS
            IF (VCQ.GT.-VSHETQ) THEN
C  CASE B1
C  REDUCE NORMAL VELOCITY DUE TO SHEATH REPULSION
              VELSH=-VC+SQRT(VCQ+VSHETQ)
              VX=VEL*VELX+VELSH*CRTX
              VY=VEL*VELY+VELSH*CRTY
              VZ=VEL*VELZ+VELSH*CRTZ
              V=SQRT(VX*VX+VY*VY+VZ*VZ)
              VELX=VX/V
              VELY=VY/V
              VELZ=VZ/V
              E0=E0+ESHET
              VEL=SQRT(E0)*RSQDVI(IION)
              EELFI(IION,ISTRA)=EELFI(IION,ISTRA)+ESHET*WPR
            ELSEIF (VCQ.LT.-VSHETQ) THEN
C  CASE B2: REFLECTION AT SHEATH, NO CHANGE IN ENERGY
C           NO SURFACE TALLY SCORING, RETURN TO CALLING PROGRAM
              VC2=VC+VC
              VX=VEL*VELX-VC2*CRTX
              VY=VEL*VELY-VC2*CRTY
              VZ=VEL*VELZ-VC2*CRTZ
              V=SQRT(VX*VX+VY*VY+VZ*VZ)
              VELX=VX/V
              VELY=VY/V
              VELZ=VZ/V
C             E0=E0
C             VEL=SQRT(E0)*RSQDVI(IION)

cdr unfinished option: store trajectories
cdr           NLTRJ = .FALSE.
cdr           TRAJ(ITRJ)%TRJ%NO_SURF = MSURF

              IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
                WRITE (IUNOUT,*) 'REFLECT FROM SHEATH AT MSURF= ',MSURF
                WRITE (IUNOUT,*) 'NEW SPEED UNIT VECTOR: VELX,VELY,VELZ'
                WRITE (IUNOUT,*)  VELX,VELY,VELZ
                EVCQ=VCQ/(RSQDVI(IION)*RSQDVI(IION))
                WRITE (IUNOUT,*)  'VEL,ESHEAT,EVCQ,E0 ',
     .                             VEL,ESHET ,EVCQ,E0
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
              ENDIF
              IRET = 1
              RETURN
            ENDIF  !  REFLECT AT SHEATH OR TRANSMISSION ?
          ENDIF !  SHEATH > < 0 ?

        ENDIF ! ILIIN > 0 ?

        FMASS=DBLE(NMASSI(IION))
        FCHAR=DBLE(NCHARI(IION))

      ENDIF
cdr sheath done

      IND=1
      CALL EIRENE_UPDATE_SURFACE(ITYP_OLD,IOLD,WPR,IND)
C
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

!  PARTICLE TYPE AND SPECIES MAY HAVE CHANGED
!  PREPARE POINTER FOR UNIFIED SUBROUTINES
      CALL EIRENE_SWITCH_PARTINFO   !dr  HERE NOT NEEDED?

C
   10 CONTINUE
C
C  ADDITIONAL OUTGOING SURFACE FLUX TALLIES
      IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WPR,1)
      IF (NADSPC_S.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WPR,1,0)
C
C  STOP TRAJECTORY, WITHOUT SPUTTERING, FOR SOME REASON IN SUBROUTINE ADDCOL OR STDCOL
C
      IF (.NOT.LGPART.AND.WPR.LE.0.0.AND.ICOL.EQ.0) THEN
        IRET = 0
        RETURN
      END IF

C  PURELY ABSORBING SURFACE. ALSO: NO SPUTTERING HERE
      IF (ILIIN(MSURF).EQ.2) GOTO 50
C
C   ........................
C   .                      .
C   .  CALL SPUTTER MODEL  .
C   ........................
C
C
C  TENTATIVELY ASSUME: NO SPUTTERED PARTICLES WILL BE FOLLOWED
      NLSPUT=.FALSE.
C
      WGHTSP=0.
      WGHTSC=0.
      YIELD1=0.
      YIELD2=0.
      ISSPTP=0
      ISSPTC=0
C
      IF (ILSPT(MSURF).NE.0) THEN
C
C  AT THIS POINT: ILIIN.EQ.1
C  SAVE PARAMETERS OF INCIDENT PARTICLE
        E0S=E0
        WEIGHS=WEIGHT
        VELS=VEL
        VELXS=VELX
        VELYS=VELY
        VELZS=VELZ
        ISPZS=ISPZ
        CALL EIRENE_SPUTR1(WMINS,FMASS,FCHAR,FLX,
     .              ISRS(ISPZ,MSURF),
     .              YIELD1,
     .              ISSPTP,ESPTP,VSPTP,VXSPTP,VYSPTP,VZSPTP,
     .              ISRC(ISPZ,MSURF),
     .              YIELD2,
     .              ISSPTC,ESPTC,VSPTC,VXSPTC,VYSPTC,VZSPTC,
     .              YSPTWL)
        NLSPUT=YIELD1.GT.0..OR.YIELD2.GT.0.
        WGHTSP=WPR*YIELD1
        WGHTSC=WPR*YIELD2
C
C  UPDATE SPUTTERED FLUX IF AVAILABLE. SORTED BY INCIDENT PARTICLE TYPE
C
C  SCORING SHIFTED TO SUBROUTINE EIRENE_UPDATE_SPTFLX, CALLED SEPARATELY
C  FOR PHYSICAL AND CHEMICAL SPUTTERING RESP.
C  at this point: ityp, iatm,.... incident particle, not the sputtered particle
C
        ITYP_OLD = ITYP
        IGASP_OLD = ISRS(ISPZ,MSURF)
        IGASC_OLD = ISRC(ISPZ,MSURF)
C
      ENDIF ! SPUTTER MODEL DONE, SO FAR.
CDR SCORING, SPLITTING ETC. WILL BE DONE FURTHER BELOW
C
   50 CONTINUE
C
C   ...................................................................
C   .                                                                 .
C   .  "SEMI-TRANSPARENT" SURFACE, MADE
c   .        TRANSPARENT FOR A FRACTION "TRANSP" OF THE FLUX          .
C   ...................................................................
C
C
      LTRANS=.FALSE.
      IF (TRANSP(ISPZ,1,MSURF).GT.0.D0.OR.
     .    TRANSP(ISPZ,2,MSURF).GT.0.D0) THEN
C
C  AT THIS POINT: ILIIN(MSURF).GT.0
C
        IF (SG.GT.0) ISG=1
        IF (SG.LT.0) ISG=2
        RATR=RANF_EIRENE( )
        LTRANS=RATR.LE.TRANSP(ISPZ,ISG,MSURF)

        IF (LTRANS) THEN
C  A NON-TRANSPARENT SURFACE IS MADE TRANSPARENT FOR THIS
C  PARTICULAR PARTICLE
C  STANDARD OR ADDITIONAL SURFACE?
          MS=MSURF
          IF (MSURF.GT.NLIM) THEN
            ISTS=MSURF-NLIM
            MS=-ISTS
            IF (INUMP(ISTS,1).NE.0) IDIM=1
            IF (INUMP(ISTS,2).NE.0) IDIM=2
            IF (INUMP(ISTS,3).NE.0) IDIM=3
C  CELL NUMBER SWITCHES LIKE A TRANSPARENT DEFAULT STANDARD SURFACE IN STDCOL
            IF (IDIM.EQ.1) NRCELL=NRCELL+NINCX
            IF (IDIM.EQ.2) THEN
              NPCELL=NPCELL+NINCY
              IPOLG=MPSURF
            END IF
            IF (IDIM.EQ.3) NTCELL=NTCELL+NINCZ
          ENDIF
          IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
            CALL EIRENE_LEER(1)
            WRITE (iunout,*) 'SURFACE MSURF= ',MS,' IS MADE TRANSPARENT'
            WRITE (iunout,*) 'ORIENTATION, SPECIES: ',ISG,ISPZ
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
          ENDIF
        ENDIF
C
      ENDIF
C
C
C   ...................................
C   .                                 .
C   .  PERFECTLY ABSORBING SURFACES   .
C   ...................................
C
C
C  NOTHING ELSE TO BE DONE, RETURN
C
      IF (ILIIN(MSURF).EQ.2.AND..NOT.LTRANS) THEN
        IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WPR
          IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WPR
          END IF
          LMETSPW(ISPZ) = .TRUE.
        ENDIF

cdr  apparently an unfinished option ?
cdr     NLTRJ = .FALSE.
cdr     TRAJ(ITRJ)%TRJ%NO_SURF = MSURF

        IF (ICOL.EQ.1) THEN
          IRET = 3
          RETURN                ! conditional expectation estimator.
                                ! Continue.
        END IF
C
        WEIGHT=0.D0
        LGPART=.FALSE.

        IRET = 0
        RETURN
      ENDIF
C
C   ................................................
C   .                                              .
C   .   MIRROR                                     .
C   .   REEMITTED FLUX=INCOMING FLUX AND RETURN    .
C   ................................................
C
      IF (ILIIN(MSURF).EQ.3.AND..NOT.LTRANS) THEN
C
C ITYP_OLD=ITNEW=ITYP
cdr this next stuff can go into update_surface(wpr,2)
C
        IF (ITYP.EQ.1) THEN
          IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          ENDIF
          IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT2(IATM,IATM)=PRFAAT2(IATM,IATM)+WPR
              LMETSPW2(1:NATM,0:NATM) => LMETSPW(NSPZTOTWS+1:NTWS_AA)
              LMETSPW2(IATM,0) = .TRUE.
              LMETSPW2(IATM,IATM) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM+IATM*NATM,MSURF)=
     .         PRFAAT(IATM+IATM*NATM,MSURF)+WPR
              LMETSPW(NSPZTOTWS+IATM) = .TRUE.
              LMETSPW(NSPZTOTWS+IATM+IATM*NATM) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFAAT .OR. LPRFAAT) LMETSPW(NSPH+IATM) =.TRUE.
        ELSEIF (ITYP.EQ.2) THEN
          IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          ENDIF
          IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML2(IMOL,IMOL)=PRFMML2(IMOL,IMOL)+WPR
              LMETSPW2(1:NMOL,0:NMOL) => LMETSPW(NTWS_MA+1:NTWS_MM)
              LMETSPW2(IMOL,0) = .TRUE.
              LMETSPW2(IMOL,IMOL) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL+IMOL*NMOL,MSURF)=
     .         PRFMML(IMOL+IMOL*NMOL,MSURF)+WPR
              LMETSPW(NTWS_MA+IMOL) = .TRUE.
              LMETSPW(NTWS_MA+IMOL+IMOL*NMOL) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFMML .OR. LPRFMML) LMETSPW(NSPA+IMOL) =.TRUE.
        ELSEIF (ITYP.EQ.3) THEN
          IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          ENDIF
          IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO2(IION,IION)=PRFIIO2(IION,IION)+WPR
              LMETSPW2(1:NION,0:NION) => LMETSPW(NTWS_IM+1:NTWS_II)
              LMETSPW2(IION,0) = .TRUE.
              LMETSPW2(IION,IION) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION+IION*NION,MSURF)=
     .         PRFIIO(IION+IION*NION,MSURF)+WPR
              LMETSPW(NTWS_IM+IION) = .TRUE.
              LMETSPW(NTWS_IM+IION+IION*NION) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFIIO .OR. LPRFIIO) LMETSPW(NSPAM+IION) =.TRUE.
        ELSEIF (ITYP.EQ.0) THEN
          IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WPR
          ENDIF
          IF (LPRFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WPR
            IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT2(IPHOT,IPHOT)=PRFPHPHT2(IPHOT,IPHOT)+WPR
              LMETSPW2(1:NPHOT,0:NPHOT) => LMETSPW(NTWS_PHM+1:NTWS_PHPH)
              LMETSPW2(IPHOT,0) = .TRUE.
              LMETSPW2(IPHOT,IPHOT) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURF)=
     .         PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURF)+WPR
              LMETSPW(NTWS_PHM+IPHOT) = .TRUE.
              LMETSPW(NTWS_PHM+IPHOT+IPHOT*NPHOT) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFPHPHT .OR. LPRFPHPHT) LMETSPW(IPHOT) =.TRUE.
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT2(IATM,IATM)=PRFAAT2(IATM,IATM)+WPR
                LMETSPW2(1:NATM,0:NATM) => LMETSPW(NSPZTOTWS+1:NTWS_AA)
                LMETSPW2(IATM,0) = .TRUE.
                LMETSPW2(IATM,IATM) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT(IATM+IATM*NATM,MSURFG)=
     .           PRFAAT(IATM+IATM*NATM,MSURFG)+WPR
                LMETSPW(NSPZTOTWS+IATM) = .TRUE.
                LMETSPW(NSPZTOTWS+IATM+IATM*NATM) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFAAT .OR. LPRFAAT) LMETSPW(NSPH+IATM) =.TRUE.
          ELSEIF (ITYP.EQ.2) THEN
            IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML2(IMOL,IMOL)=PRFMML2(IMOL,IMOL)+WPR
                LMETSPW2(1:NMOL,0:NMOL) => LMETSPW(NTWS_MA+1:NTWS_MM)
                LMETSPW2(IMOL,0) = .TRUE.
                LMETSPW2(IMOL,IMOL) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML(IMOL+IMOL*NMOL,MSURFG)=
     .           PRFMML(IMOL+IMOL*NMOL,MSURFG)+WPR
                LMETSPW(NTWS_MA+IMOL) = .TRUE.
                LMETSPW(NTWS_MA+IMOL+IMOL*NMOL) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFMML .OR. LPRFMML) LMETSPW(NSPA+IMOL) =.TRUE.
          ELSEIF (ITYP.EQ.3) THEN
            IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO2(IION,IION)=PRFIIO2(IION,IION)+WPR
                LMETSPW2(1:NION,0:NION) => LMETSPW(NTWS_IM+1:NTWS_II)
                LMETSPW2(IION,0) = .TRUE.
                LMETSPW2(IION,IION) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO(IION+IION*NION,MSURFG)=
     .           PRFIIO(IION+IION*NION,MSURFG)+WPR
                LMETSPW(NTWS_IM+IION) = .TRUE.
                LMETSPW(NTWS_IM+IION+IION*NION) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFIIO .OR. LPRFIIO) LMETSPW(NSPAM+IION) =.TRUE.
          ELSEIF (ITYP.EQ.0) THEN
            IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WPR
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT2(IPHOT,IPHOT)=PRFPHPHT2(IPHOT,IPHOT)+WPR
                LMETSPW2(1:NPHOT,0:NPHOT)=>LMETSPW(NTWS_PHM+1:NTWS_PHPH)
                LMETSPW2(IPHOT,0) = .TRUE.
                LMETSPW2(IPHOT,IPHOT) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURFG)=
     .           PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURFG)+WPR
                LMETSPW(NTWS_PHM+IPHOT) = .TRUE.
                LMETSPW(NTWS_PHM+IPHOT+IPHOT*NPHOT) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFPHPHT .OR. LPRFPHPHT) LMETSPW(IPHOT) =.TRUE.
          ENDIF
        ENDIF
      ENDIF
C
C  EITHER: SEMI-TRANSPARENT SURFACE....
C
      IF (LTRANS) THEN
C  CONTINUE WITH UNMODIFIED VELOCITY.

C  COMPENSATE INCIDENT SURFACE FLUX TALLY CONTRIBUTIONS
C  SCORED ABOVE.
        IND=1
        CALL EIRENE_UPDATE_SURFACE (ITYP_OLD,IOLD,-WPR,IND)
cdr  for additional surface tallies: this is treated as re-emitted flux.
        IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WPR,2)
        IF (NADSPC.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WPR,2,0)
        COLFLAG = .TRUE.
        IRET = 2
        RETURN
C
C  ... OR: PERFECT SPECULAR REFLECTION, MIRROR,
C
      ELSEIF (ILIIN(MSURF).EQ.3) THEN   ! (and here: .not ltrans)
cdr  why no call to update_surface for emitted fluxes?
cdr  Because already done above, under: if (iliin=3 and not ltransp)
        COSI2=-2.*(VELX*CRTX+VELY*CRTY+VELZ*CRTZ)
        VELX=VELX+COSI2*CRTX
        VELY=VELY+COSI2*CRTY
        VELZ=VELZ+COSI2*CRTZ
        IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WPR,2)
        IF (NADSPC.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WPR,2,0)

cdr  redundant code from an unfinished option?
C       NLTRJ = .FALSE.
C       TRAJ(ITRJ)%TRJ%NO_SURF = MSURF
cdr
        IF (ICOL.EQ.1) THEN
          IRET = 3
          RETURN
        END IF
        IRET = 1
        RETURN
      ENDIF
C
C   .........................
C   .                       .
C   .  TRANSPARENT SURFACE  .
C   .........................
C
      IF (ILIIN(MSURF).LT.0) THEN
C
C  ONE-SIDED FLUX: UPDATE ONLY NEGATIVE COMPONENT HERE,
C                  POSITIVE COMPONENT, SG.GT.0, WAS ALREADY ON "OT-TALLIES"
C  IN CASE ILIIN=-3: NET FLUXES HAVE ALREADY BEEN UPDATED ABOVE ON "OT-TALLIES".
C                    NEED NOT BE UPDATED AGAIN HERE.
C
        IF ((SG.GT.0.D0).OR.(ILIIN(MSURF).EQ.-3)) GOTO 90
C
C  HERE: ILIIN NE -3, AND SG LE 0, SCORE ONE-SIDED "NEGATIVE" CURRENTS (WPR <=0)
C
C ITYP_OLD=ITNEW=ITYP
C
        IF (ITYP.EQ.1) THEN
          IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFAAT(IATM,MSURF)=ERFAAT(IATM,MSURF)+E0*WPR
          ENDIF
          IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFAAT(IATM,MSURF)=PRFAAT(IATM,MSURF)+WPR
            IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT2(IATM,IATM)=PRFAAT2(IATM,IATM)+WPR
              LMETSPW2(1:NATM,0:NATM)=>LMETSPW(NSPZTOTWS+1:NTWS_AA)
              LMETSPW2(IATM,0) = .TRUE.
              LMETSPW2(IATM,IATM) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM+IATM*NATM,MSURF)=
     .         PRFAAT(IATM+IATM*NATM,MSURF)+WPR
              LMETSPW(NSPZTOTWS+IATM) = .TRUE.
              LMETSPW(NSPZTOTWS+IATM+IATM*NATM) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFAAT .OR. LPRFAAT) LMETSPW(NSPH+IATM) =.TRUE.
        ELSEIF (ITYP.EQ.2) THEN
          IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFMML(IMOL,MSURF)=ERFMML(IMOL,MSURF)+E0*WPR
          ENDIF
          IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFMML(IMOL,MSURF)=PRFMML(IMOL,MSURF)+WPR
            IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML2(IMOL,IMOL)=PRFMML2(IMOL,IMOL)+WPR
              LMETSPW2(1:NMOL,0:NMOL)=>LMETSPW(NTWS_MA+1:NTWS_MM)
              LMETSPW2(IMOL,0) = .TRUE.
              LMETSPW2(IMOL,IMOL) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL+IMOL*NMOL,MSURF)=
     .         PRFMML(IMOL+IMOL*NMOL,MSURF)+WPR
              LMETSPW(NTWS_MA+IMOL) = .TRUE.
              LMETSPW(NTWS_MA+IMOL+IMOL*NMOL) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFMML .OR. LPRFMML) LMETSPW(NSPA+IMOL) =.TRUE.
        ELSEIF (ITYP.EQ.3) THEN
          IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFIIO(IION,MSURF)=ERFIIO(IION,MSURF)+E0*WPR
          ENDIF
          IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFIIO(IION,MSURF)=PRFIIO(IION,MSURF)+WPR
            IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO2(IION,IION)=PRFIIO2(IION,IION)+WPR
              LMETSPW2(1:NION,0:NION)=>LMETSPW(NTWS_IM+1:NTWS_II)
              LMETSPW2(IION,0) = .TRUE.
              LMETSPW2(IION,IION) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION+IION*NION,MSURF)=
     .         PRFIIO(IION+IION*NION,MSURF)+WPR
              LMETSPW(NTWS_IM+IION) = .TRUE.
              LMETSPW(NTWS_IM+IION+IION*NION) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFIIO .OR. LPRFIIO) LMETSPW(NSPAM+IION) =.TRUE.
        ELSEIF (ITYP.EQ.0) THEN
          IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            ERFPHPHT(IPHOT,MSURF)=ERFPHPHT(IPHOT,MSURF)+E0*WPR
          ENDIF
          IF (LPRFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            PRFPHPHT(IPHOT,MSURF)=PRFPHPHT(IPHOT,MSURF)+WPR
            IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
              PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURF)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT2(IPHOT,IPHOT)=PRFPHPHT2(IPHOT,IPHOT)+WPR
              LMETSPW2(1:NPHOT,0:NPHOT)=>LMETSPW(NTWS_PHI+1:NTWS_PHPH)
              LMETSPW2(IPHOT,0) = .TRUE.
              LMETSPW2(IPHOT,IPHOT) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURF)=
     .         PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURF)+WPR
              LMETSPW(NTWS_PHM+IPHOT) = .TRUE.
              LMETSPW(NTWS_PHM+IPHOT+IPHOT*NPHOT) = .TRUE.
#endif
            END IF
          ENDIF
          IF (LERFPHPHT .OR. LPRFPHPHT) LMETSPW(IPHOT) =.TRUE.
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (ITYP.EQ.1) THEN
            IF (LERFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFAAT(IATM,MSURFG)=ERFAAT(IATM,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFAAT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFAAT(IATM,MSURFG)=PRFAAT(IATM,MSURFG)+WPR
              IF (NLSPCSCL_ATM) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFAAT2(1:NATM,0:NATM) => PRFAAT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT2(IATM,IATM)=PRFAAT2(IATM,IATM)+WPR
                LMETSPW2(1:NATM,0:NATM)=>LMETSPW(NSPZTOTWS+1:NTWS_AA)
                LMETSPW2(IATM,0) = .TRUE.
                LMETSPW2(IATM,IATM) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFAAT(IATM+IATM*NATM,MSURFG)=
     .           PRFAAT(IATM+IATM*NATM,MSURFG)+WPR
                LMETSPW(NSPZTOTWS+IATM) = .TRUE.
                LMETSPW(NSPZTOTWS+IATM+IATM*NATM) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFAAT .OR. LPRFAAT) LMETSPW(NSPH+IATM) =.TRUE.
          ELSEIF (ITYP.EQ.2) THEN
            IF (LERFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFMML(IMOL,MSURFG)=ERFMML(IMOL,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFMML) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFMML(IMOL,MSURFG)=PRFMML(IMOL,MSURFG)+WPR
              IF (NLSPCSCL_MOL) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFMML2(1:NMOL,0:NMOL) => PRFMML(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML2(IMOL,IMOL)=PRFMML2(IMOL,IMOL)+WPR
                LMETSPW2(1:NMOL,0:NMOL)=>LMETSPW(NTWS_MA+1:NTWS_MM)
                LMETSPW2(IMOL,0) = .TRUE.
                LMETSPW2(IMOL,IMOL) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFMML(IMOL+IMOL*NMOL,MSURFG)=
     .           PRFMML(IMOL+IMOL*NMOL,MSURFG)+WPR
                LMETSPW(NTWS_MA+IMOL) = .TRUE.
                LMETSPW(NTWS_MA+IMOL+IMOL*NMOL) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFMML .OR. LPRFMML) LMETSPW(NSPA+IMOL) =.TRUE.
          ELSEIF (ITYP.EQ.3) THEN
            IF (LERFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFIIO(IION,MSURFG)=ERFIIO(IION,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFIIO) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFIIO(IION,MSURFG)=PRFIIO(IION,MSURFG)+WPR
              IF (NLSPCSCL_ION) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFIIO2(1:NION,0:NION) => PRFIIO(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO2(IION,IION)=PRFIIO2(IION,IION)+WPR
                LMETSPW2(1:NION,0:NION)=>LMETSPW(NTWS_IM+1:NTWS_II)
                LMETSPW2(IION,0) = .TRUE.
                LMETSPW2(IION,IION) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFIIO(IION+IION*NION,MSURFG)=
     .           PRFIIO(IION+IION*NION,MSURFG)+WPR
                LMETSPW(NTWS_IM+IION) = .TRUE.
                LMETSPW(NTWS_IM+IION+IION*NION) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFIIO .OR. LPRFIIO) LMETSPW(NSPAM+IION) =.TRUE.
          ELSEIF (ITYP.EQ.0) THEN
            IF (LERFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              ERFPHPHT(IPHOT,MSURFG)=ERFPHPHT(IPHOT,MSURFG)+E0*WPR
            ENDIF
            IF (LPRFPHPHT) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              PRFPHPHT(IPHOT,MSURFG)=PRFPHPHT(IPHOT,MSURFG)+WPR
              IF (NLSPCSCL_PHOT) THEN
#if ( !defined(G95) && !defined(LEGACYCOMP) )
                PRFPHPHT2(1:NPHOT,0:NPHOT) => PRFPHPHT(:,MSURFG)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT2(IPHOT,IPHOT)=PRFPHPHT2(IPHOT,IPHOT)+WPR
                LMETSPW2(1:NPHOT,0:NPHOT)=>LMETSPW(NTWS_PHI+1:NTWS_PHPH)
                LMETSPW2(IPHOT,0) = .TRUE.
                LMETSPW2(IPHOT,IPHOT) = .TRUE.
#else
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURFG)=
     .           PRFPHPHT(IPHOT+IPHOT*NPHOT,MSURFG)+WPR
                LMETSPW(NTWS_PHM+IPHOT) = .TRUE.
                LMETSPW(NTWS_PHM+IPHOT+IPHOT*NPHOT) = .TRUE.
#endif
              END IF
            ENDIF
            IF (LERFPHPHT .OR. LPRFPHPHT) LMETSPW(IPHOT) =.TRUE.
          ENDIF
        ENDIF
C
   90   CONTINUE
        IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WPR,2)
        IF (NADSPC.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WPR,2,0)
        colflag = .true.
        IRET = 2
        RETURN
      ENDIF
C
C
C   .............................
C   .                           .
C   .  NON-TRANSPARENT SURFACE  .
C   .............................
C
      XGENER=0.
C
C  ALL INCIDENT SURFACE TALLIES, AND
C  INCLUDING CONDITIONAL EXPECTATION CORRECTIONS, ARE UPDATED.
C  TRANSPARENT, MIRROR AND ABSORBING SURFACES ARE ALSO DONE
C
C  CONDITIONAL EXPECTATION ESTIMATOR: HAS THIS PARTICLE COLLIDED IN THE VOLUME,
C  BEFORE IT HIT THE WALL?
C
      IF (ICOL.EQ.1) THEN
cdr     NLTRJ = .FALSE.
cdr     TRAJ(ITRJ)%TRJ%NO_SURF = MSURF
        IRET = 3
        RETURN
      END IF
C
C  ..........................................................................
C
C  NOW DEAL WITH REFLECTED AND/OR SPUTTERED PARTICLES.
C  ..........................................................................
C
C  REFLECTION FROM SURFACE
C
C  ......................................................
C  .                                                    .
C  .  REFLECTION MODEL 600--699 FOR INCIDENT MOLECULES  .
C  ......................................................
C
cdr  tbd: make a universal "thermal" reflection model.
cdr       called from here for molecules.
cdr       called from reflec1 for atoms etc.
cdr       called from reflec_photons for photons...
C
      IF (ITYP.EQ.2) THEN
C
C  TEST FOR REFLECTION
C
        IF (WEIGHT.GE.WMINS) THEN
C  WITH SUPPRESSION OF ABSORPTION
          WABS=WEIGHT*(1.D0-RECYCT(ISPZ,MSURF))
          IF (WABS.GT.0.D0) THEN
            IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WABS
              IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
                SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WABS
              END IF
              LMETSPW(ISPZ) = .TRUE.
            ENDIF
            WEIGHT=WEIGHT-WABS
          ENDIF
          IF (WEIGHT.GT.EPS30) GOTO 610
          LGPART=.FALSE.
          IRET = 0
          RETURN
        ELSEIF (RECYCT(ISPZ,MSURF).NE.1.D0) THEN
C  NO SUPPRESSION OF ABSORPTION
          ZVZ=RANF_EIRENE( )
          IF (ZVZ.LT.RECYCT(ISPZ,MSURF)) GOTO 610
C  ABSORB THIS PARTICLE
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
            IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WEIGHT
            END IF
            LMETSPW(ISPZ) = .TRUE.
          ENDIF
          LGPART=.FALSE.
          IRET = 0
          RETURN
        ENDIF
C
  610   CONTINUE
C
C  NEW SPECIES: AGAIN MOLECULE
C
C       ITYP=2
        IMOL = ISRT(ISPZ,MSURF)
        MOL_DEFAULT = IMOL
        MODREF=ILREF(ISPZ,MSURF)
        IF (MODREF.GE.9) CALL EIRENE_RF2USR (IMOL,MOL_DEFAULT)
        IF (IMOL.GT.NMOLI) THEN
          FR2=RANF_EIRENE( )
          DO 621 I=1,NMOLI
            IMOL=I
            IF (FR2.LE.DMOL(IMOL)) GOTO 622
  621     CONTINUE
          GOTO 995
  622     CONTINUE
        ELSEIF (IMOL.EQ.0) THEN
C  NO THERMAL EMISSION, ABSORB INSTEAD
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
            IF (MSURFG.GT.0) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
              SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WEIGHT
            END IF
            LMETSPW(ISPZ) = .TRUE.
          ENDIF
          LGPART=.FALSE.
          IRET = 0
          RETURN
        ELSEIF (IMOL.LT.0) THEN
C  NOT IN USE
          GOTO 995
        ENDIF
C
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

!  PARTICLE TYPE AND SPECIES MAY HAVE CHANGED
!  IF SO, THEN PREPARE POINTERS FOR UNIFIED SUBROUTINES
cdr only if ispz .ne. ispz_old:
        CALL EIRENE_SWITCH_PARTINFO
C
        E0TERM=EWALL(MSURF)
        IF (E0TERM.GT.0.0_DP) THEN
C  MONOENERGETIC DISTRIBUTION
          E0=E0TERM
          VEL=RSQDVM(IMOL)*SQRT(E0)
C   AZIMUTHAL ANGLE: EQUIDISTRIBUTION
C   POLAR ANGLE: COSINE
          IF (INIV4.EQ.0) CALL EIRENE_FCOSIN
          VX=FC1(INIV4)
          VY=FC2(INIV4)
          VZ=FC3(INIV4)
          INIV4=INIV4-1
          CALL EIRENE_ROTATF (VELX,VELY,VELZ,VX,VY,VZ,CRTX,CRTY,CRTZ)
        ELSEIF (E0TERM.LT.0.0_DP) THEN
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
          TW=-E0TERM
c these variables are necessary as the corresponding arguments in velocs
c are INTENT(IN) !
          VXR = 0.0_DP
          VYR = 0.0_DP
          VZR = 0.0_DP
          VWL = 0.0_DP !  INDICATE: NON-DRIFTING MAXWELLIAN FLUX,
                       !            WEIGHT IS NOT MODIFIED
          WGHTVS= WEIGHT
          CALL EIRENE_VELOCS(WGHTVS,
     .       TW,0._DP,VWL,VXR,VYR,VZR,RSQDVM(IMOL),
     .                      CVRSSM(IMOL),
     .                     -CRTX,-CRTY,-CRTZ,
     .                      E0,VELX,VELY,VELZ,VEL)
        ELSE
          WRITE (iunout,*) 'EWALL = 0: MODEL NOT AVAILABLE'
          WRITE (iunout,*) 'ERROR IN ESCAPE, EXIT CALLED'
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
C
C  ...............................................
C  .                                             .
C  .  REFLECTION MODEL FOR ATOMS OR TEST IONS    .
C  ...............................................
C
      ELSEIF (ITYP.EQ.1.OR.ITYP.EQ.3) THEN
C
C
C  THESE INCIDENT PARTICLES MAY HAVE SPUTTERED AT THIS  SURFACE
C
        SPLFLG=0.
C  UPDATE TOTAL SPUTTERED FLUXES IN CASE OF UNKNOWN SPECIES INDEX FOR SPUTTERED PARTICLE
C  IF ISSPT... GT.0 ('EMITTED SPECIES KNOWN') THEN SPUTTERED FLUXES WILL BE UPDATED BELOW,
C                                             AND SPUTTERED PARTICLES MAY BE FOLLOWED
        IF (WGHTSP.GT.0..AND.ISSPTP.EQ.0)
     .    CALL EIRENE_UPDATE_SPTFLX (ITYP_OLD,IOLD,WGHTSP,1)
        IF (WGHTSC.GT.0..AND.ISSPTC.EQ.0)
     .    CALL EIRENE_UPDATE_SPTFLX (ITYP_OLD,IOLD,WGHTSC,1)
C
        IF (WGHTSP.GT.0..AND.ISSPTP.GT.0) THEN
C  PHYSICAL SPUTTERING, RESTORE PHYSICALLY SPUTTERED PARTICLE PARAMETERS
C  SCORE ALL RELEVANT TALLIES
C
          SPLFLG=SPLFLG+1.
          ISPZ=ISSPTP
          ITYP=ISPEZI(ISPZ,-1)
          IPHOT=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTP
          WEIGHT=WGHTSP
          VEL=VSPTP
          VELX=VXSPTP
          VELY=VYSPTP
          VELZ=VZSPTP
C
C  ITYP_OLD.NE.ITNEW=ITYP POSSIBLE
C
          CALL EIRENE_UPDATE_SPTFLX (ITYP_OLD,IOLD,WGHTSP,2)

          IND=2
          IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WGHTSP,IND)
          IF (NADSPC_S.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WGHTSP,IND,0)

          IF (IGASP_OLD.EQ.0) GOTO 4711 ! SCORE SPUTTERED PARTICLES
                                        ! ON SURFACE TALLIES ONLY
C                                         IF THEY ARE FOLLOWED.
C                                         OTHERWISE: ONLY ON SPUTTER TALLIES
          CALL EIRENE_UPDATE_SURFACE (ITYP_OLD,IOLD,WGHTSP,IND)
C
C.....................................................................
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
C  SPLITTING
C
          IF (NLEVEL.GE.MAXLEVEL) THEN
            WRITE (iunout,*) 'SPLITTING ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NLEVEL: ',NLEVEL
C   SPUTTERED PARTICLE HAS SCORED, BUT WILL NOT BE FOLLOWED
            GOTO 4711
          ENDIF
          NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
          DO 533 J=1,NPARTC
            RSPLST(J,NLEVEL)=RPST(J)
  533     CONTINUE
          DO 534 J=1,MPARTC
            ISPLST(J,NLEVEL)=IPST(J)
  534     CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
          NODES(NLEVEL)=2
C
        ENDIF
C
C  SPLITTING FOR PHYSICAL SPUTTERING DONE

C
 4711   IF (WGHTSC.GT.0..AND.ISSPTC.GT.0) THEN
C
C  CHEMICAL SPUTTERING, RESTORE CHEMICALLY SPUTTERED PARTICLE PARAMETERS
C  SCORE ALL RELEVANT TALLIES
C
          SPLFLG=SPLFLG+1.
          ISPZ=ISSPTC
          ITYP=ISPEZI(ISPZ,-1)
          IPHOT=ISPEZI(ISPZ,0)
          IATM=ISPEZI(ISPZ,1)
          IMOL=ISPEZI(ISPZ,2)
          IION=ISPEZI(ISPZ,3)
          IPLS=ISPEZI(ISPZ,4)
          E0=ESPTC
          WEIGHT=WGHTSC
          VEL=VSPTC
          VELX=VXSPTC
          VELY=VYSPTC
          VELZ=VZSPTC
C  ITYP_OLD.NE.ITNEW=ITYP POSSIBLE
C
          CALL EIRENE_UPDATE_SPTFLX (ITYP_OLD,IOLD,WGHTSC,2)

          IND=2
          IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WGHTSC,IND)
          IF (NADSPC_S.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WGHTSC,IND,0)

          IF (IGASC_OLD.EQ.0) GOTO 4712  ! SCORE SPUTTERED PARTICLES
                                         ! ON SURFACE TALLIES ONLY
C                                          IF THEY ARE FOLLOWED.
C                               OTHERWISE: ONLY ON SPUTTER TALLIES
          CALL EIRENE_UPDATE_SURFACE (ITYP_OLD,IOLD,WGHTSC,IND)
C
C.....................................................................
C  FOLLOW SPUTTERED PARTICLES LATER. PUT THEM INTO STATISTICAL CELLAR
C  SPLITTING
C

          IF (NLEVEL.GE.MAXLEVEL) THEN
            WRITE (iunout,*) 'SPLITTING ABANDONED FOR PART. NO. ',NPANU
            WRITE (iunout,*) 'CASCADE OVERFLOW: NLEVEL: ',NLEVEL
C   SPUTTERED PARTICLE HAS SCORED, BUT WILL NOT BE FOLLOWED
            GOTO 4712
          ENDIF

          NLEVEL=NLEVEL+1
C  SAVE LOCATION, WEIGHT AND OTHER PARAMETERS AT CURRENT LEVEL
          DO 535 J=1,NPARTC
            RSPLST(J,NLEVEL)=RPST(J)
  535     CONTINUE
          DO 536 J=1,MPARTC
            ISPLST(J,NLEVEL)=IPST(J)
  536     CONTINUE
C  NUMBER OF NODES AT THIS LEVEL
          NODES(NLEVEL)=2
C
        ENDIF
C
C  SPLITTING FOR CHEMICAL SPUTTERING DONE.
 4712   CONTINUE
C.................................................................
C
C  RESTORE INCIDENT PARTICLE, FOR SURFACE REFLECTION ROUTINE
C
        IF (SPLFLG.NE.0) THEN
          E0=E0S
          WEIGHT=WEIGHS
          VEL=VELS
          VELX=VELXS
          VELY=VELYS
          VELZ=VELZS
          ISPZ=ISPZS
          LGPART=.FALSE.
        ENDIF
C
cdr  surface model for "atomic" particles of type 1,3,4
cdr  tbd.: rename to : reflc1_atomic
        CALL EIRENE_REFLC1 (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .               ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))

        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

!  PARTICLE TYPE AND SPECIES MIGHT HAVE CHANGED
!  PREPARE POINTER FOR UNIFIED SUBROUTINES
        IF (LGPART) THEN
           CALL EIRENE_SWITCH_PARTINFO
        ELSE
          WEIGHT=0.
        ENDIF
C
C  ............................................
C  .                                          .
C  .  REFLECTION MODEL 700...799 FOR PHOTONS  .
C  ............................................
C
      ELSEIF (ITYP.EQ.0) THEN

cdr  surface model for particles type 0 ("photons")
        CALL EIRENE_REFLC1_PHOTON (WMINS,FMASS,FCHAR,NPRT(ISPZ),
     .               ISRF(ISPZ,MSURF),ISRT(ISPZ,MSURF))
        ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
        IF (LGPART) THEN
          CALL EIRENE_SWITCH_PARTINFO
        ELSE
          WEIGHT=0.
        ENDIF
      ENDIF
C
C  UPDATE REFLECTED PARTICLE AND ENERGY FLUX
C
      IRET = 0
      IF (.NOT.LGPART) RETURN
C
C  ITYP_OLD.NE.ITNEW=ITYP POSSIBLE
C
      IND=2
      CALL EIRENE_UPDATE_SURFACE (ITYP_OLD,IOLD,WEIGHT,IND)

      IF (NADSI.GE.1) CALL EIRENE_UPSUSR (WEIGHT,IND)
      IF (NADSPC_S.GE.1) CALL EIRENE_UPDATE_SPECTRUM (WEIGHT,IND,0)
C
C     NLTRJ = .FALSE.
C     TRAJ(ITRJ)%TRJ%NO_SURF = MSURF

      IF  (ITYP.EQ.ITYP_OLD) THEN
        IRET = 1
        RETURN
      END IF
      IF ((ITYP.EQ.1.AND.ITYP_OLD.EQ.2).OR.
     .    (ITYP.EQ.2.AND.ITYP_OLD.EQ.1)) THEN
        IRET = 1
        RETURN
      END IF
      IRET = 0
      RETURN
C
  995 CONTINUE
      WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE IN ESCAPE '
      WRITE (iunout,*) 'IMOL, MSURF ',IMOL,MSURF
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_ESCAPE
