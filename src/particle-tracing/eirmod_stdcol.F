      MODULE EIRMOD_STDCOL
c  eirene_stdcol   CALLED FROM FOLNEUT, FOLION, FOLSTAT_NEUT, FOLSTAT_ION, LININT: idimm=1,2,3
c  eirene_stdnor   CALLED FROM LOCATE: idimm=1,2,3,4

cdr  summer 2019: TK: remove ALTERNATE RETURN, and ENTRY.
cdr  dec. 2019        comments, cleanup, add missing features: stdcol_ass

cdr  msurf=ists in case nltet, nltri ? later used_ nlim+ists ?
CDR  IDIMM=5. also CALL STDCOL_XYMIX. Difference to IDIMM=4 ? removed.
cdr  Some comments are now in wrong places, misleading.
cdr  Periodicity in nltet option: call norusr. (make optional?)
cdr  Periodicity legveo=2 and nlcrc: use polyg. data? Hidden link
cdr  bug re NLSRFA in case idimm=3 ?

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CUPD
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CTETRA
      USE EIRMOD_COMPRT
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_CTRIG
      use EIRMOD_cfplk
      use EIRMOD_csdvi
      use EIRMOD_LEARC1, ONLY: EIRENE_LEARC1
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_STDCOL, EIRENE_STDNOR

c stdcol and stdcol_absorb
      REAL(DP) :: X0SA, Y0SA, Z0SA
      INTEGER :: NACLLS, ICOS, IWEI, MSURFS, IRET_OUT

      REAL(DP) :: COSROT, SINROT, VELX_OLD

c stdcol and stdcol_x_rad
      INTEGER :: NDUM

c stdcol, stdcol_x_rad, stdcol_y_pol, stdcol_z_tor
      INTEGER, EXTERNAL :: EIRENE_IDEZ

      INTEGER :: IDIMM

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(X0SA,Y0SA,Z0SA,
!$OMP&   NACLLS,ICOS,IWEI,MSURFS,IRET_OUT,
!$OMP&   COSROT,SINROT,VELX_OLD,
!$OMP&   NDUM,IDIMM)
#endif

      CONTAINS

!  18.08.06: after particle has been moved to surface correct particle
!            positions used in timep and timet
!pb  05.10.06: option for single-sided switching of non-default standard
!pb            surfaces introduced
!pb  07.12.06: cell number corrected for nltet option
!pb  22.03.07: LEVGEO=6 --> LEVGEO=10
!pb  18.04.08: typo corrected: NLSRFZ => NLSRFY in ELSEIF (IDIMM==3) block
!pb  25.07.07: periodicity in y-direction for LEVGEO=3 introduced
!pb  07.07.09: setting of NLSRFA added
cdr  29.07.17: added lgpart=false on absorbing surfaces (otherwise problems
cdr            with trace ions onto absorbing surfaces. corresponding fix in folion.
cdr  Nov. 17 : lmetspw arguments corrected
C
      SUBROUTINE EIRENE_STDCOL (ISTS,IDIMM_IN,SG,IRET)
C
C  1.) ADVANCE TRAJECTORY TO A NON-DEFAULT STANDARD SURFACE (DISTANCE: ZT)
C  2.) SWITCHING (CELL NUMBERS, FLAGS, ETC...)
C  3A.) IF TRANSPARENT AND NO SURFACE TALLIES: CONTINUE FLIGHT
C  3B.) ELSE: IF NEEDED, PREPARE REFLECTION (SURFACE NORMALS),
C             THEN CALL ESCAPE
C
C  SG = SIGN OF COSINE OF ANGLE OF INCIDENCE
C  ISTS = SURFACE INDEX IN NSTSI ARRAYS
C  IDIMM = INDEX (1,2,3) FOR: RADIAL, POLOIDAL OR TOROIDAL SURFACE
C  RETURN IRET=1: NO SURFACE TALLIES, FLIGHT CONTINUES
C  RETURN IRET=2: SURFACE TALLIES,
C             THEN ABSORPTION, REFLECTION MODEL OR CONTINUATION OF FLIGHT
C             (CALL SUBR. ESCAPE)
C

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: SG
      INTEGER, INTENT(IN) :: ISTS, IDIMM_IN
      INTEGER, INTENT(INOUT) :: IRET

      REAL(DP) :: WINK
      INTEGER :: IDUM, IAN, IEN, EIRENE_LEARCA, EIRENE_LEARC2
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_LEARCA, EIRENE_LEARC2

C
C   "COLLISION" WITH STANDARD SURFACE NO. MSURF=NLIM+ISTS
C   OF THE RADIAL   (OR X-) GRID: IDIMM=1
C   OR THE POLOIDAL (OR Y-) GRID: IDIMM=2
C   OR THE TOROIDAL (OR Z-) GRID: IDIMM=3
C
C  SAVE DATA OF OLD POINT FOR DIAGNOSTICS
      IRET = 0
      IRET_OUT = 0
      IDIMM = IDIMM_IN
      X0SA=X0
      Y0SA=Y0
      Z0SA=Z0
      MSURFS=MSURF
      NACLLS=NACELL
C  SET NEW POINT ON NON-DEFAULT STANDARD SURFACE ISTS. FLIGHT TIME: ZT
      X0=X0+VELX*ZT
      Y0=Y0+VELY*ZT
      Z0=Z0+VELZ*ZT
      TIME=TIME+ZT/VEL
      if (nlfem.or.nltet) then
cdr  no additional surfaces (NLIM) allowed in nlfem nor in nltet?
cpb  for NLFEM and NLTET MSURF already contains NLIM
        MSURF=ISTS
      else
        MSURF=NLIM+ISTS
      endif
      ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
      SCOS=SG
      ICOS=INT(SCOS)

      IF (IMETWL(MSURF) == 0) THEN
        NWLMT = NWLMT+1
        IWLMT(NWLMT) = MSURF
        IMETWL(MSURF) = NWLMT
      END IF

C  SET IPOLG, PHI AND OTHER ASSISTANT COORDINATES
      IPOLG=IPOLGN
      IF (NLTRA) THEN
        PHI=MOD(PHI-ATAN2(Z01,X01)+ATAN2(Z0,(RMTOR+X0)),PI2A)
      ENDIF
      CALL EIRENE_STDCOL_ASS
C
C  CARRY OUT DEFAULT CELL SWITCHING AT TRANSPARENT STANDARD SURFACES
      select case (IDIMM)

      case (1)
        NLSRFX=.TRUE.
        NLSRFY=.FALSE.
        NLSRFZ=.FALSE.
        NLSRFA=.FALSE.
        IF (ILIIN(MSURF).LE.0) NRCELL=NRCELL+ICOS
C
C  IN NLFEM, NLTET AND NLGEN OPTIONS: ALL SURFACES ARE IDIMM=1 SURFACES
        IF (NLFEM) NRCELL=MRSURF
        IF (NLTET) NRCELL=MRSURF
C       IF (NLGEN) NRCELL=MRSURF
C
      case (2)
        NLSRFX=.FALSE.
        NLSRFY=.TRUE.
        NLSRFZ=.FALSE.
        NLSRFA=.FALSE.
        IF (ILIIN(MSURF).LE.0) NPCELL=NPCELL+ICOS
C
        IF (NLPLG.AND.ILIIN(MSURF).GT.0.AND.SCOS.GT.0) IPOLG=IPOLG-1
C
      case (3)
        NLSRFX=.FALSE.
        NLSRFY=.FALSE.
        NLSRFZ=.TRUE.
CDR NEXT STATEMENT IS PROBABLY WRONG.
        NLSRFA=.TRUE.
        IF (ILIIN(MSURF).LE.0) NTCELL=NTCELL+ICOS

      case default
        WRITE (iunout,*) 'ERROR EXIT IN STDCOL, IDIMM ',IDIMM
        CALL EIRENE_EXIT_OWN(1)

      end select
C
      IWEI=ILSIDE(MSURF)*ICOS
!pb   IF (IWEI.LT.0) call stdcol_absorb
cdr
cdr: july 17: by removing this statement from here
cdr           the ILSIDE options for geometry debugging are partially disabled,
cdr           at least for absorbing surfaces, for which now code segment
cdr           stdcol_absorb... is bypassed.

      IF (ILIIN(MSURF).EQ.2) THEN
c  Regular absorbing surface (from both sides).
C  SCORE FLUXES LATER IN ESCAPE (DO NOT SET WEIGHT=0.D0 HERE),
C  AND STOP THEN.
C
        IF (NLTRC) CALL EIRENE_CHCTRC(X0,Y0,Z0,16,8)
        LGPART=.FALSE.
        IRET_OUT = 2
        IRET = IRET_OUT
        RETURN
      ENDIF
C
C  OPERATE A SWITCH
C
      IF (ILSWCH(MSURF).NE.0) THEN
C
C  TURN ON OR OFF THE STANDARD GRID CALCULATION
        IF (ISWICH(1,MSURF).NE.0) ITIME=ICOS*ISWICH(1,MSURF)

C  TURN ON OR OFF MFP AND REACTION RATES: PARTICLE ENTERS VACUUM
        IF (ISWICH(2,MSURF).NE.0) IFPATH=ICOS*ISWICH(2,MSURF)

C  TURN ON OR OFF VOLUME-AVERAGED TALLIES
        IF (ISWICH(3,MSURF).NE.0) IUPDTE=ICOS*ISWICH(3,MSURF)

C  NEW ADD. CELL INDEX NACELL
C  NOTE: STANDARD SURFACES CANNOT SWITCH NACELL, NBLOCK, NRCELL,....
        IF (ISWICH(4,MSURF).NE.0) THEN
          GOTO 999
C
C  ENTRANCE INTO STANDARD MESH, INTO BLOCK NBLOCK=ILBLCK
C  OR
C  EXIT FROM STANDARD MESH, INTO CELL NACELL=ILACLL
        ELSEIF (ISWICH(5,MSURF).NE.0) THEN
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO ILACLL
            IF (ILACLL(MSURF) > 0) THEN
              NACELL=ILACLL(MSURF)
              NBLOCK=NBMLTP
              NRCELL=0
              NPCELL=1
              NTCELL=1
              IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
                IWEI=-10
                CALL EIRENE_STDCOL_ABSORB
                IRET = IRET_OUT
                RETURN
              ENDIF
            ENDIF
          ELSEIF (NACELL.GT.0) THEN
C  ENTRANCE INTO STANDARD MESH, INTO NBLOCK=ILBLCK
            NBLOCK=ILBLCK(MSURF)
            NACELL=0
C  FIND  NRCELL,IPOLG IN STANDARD MESH, BLOCK NBLOCK
            select case (idimm)

            case (1)
              NRCELL=MIN0(NR1STM,MRSURF)
              IAN=MRSURF
              IEN=MRSURF
              NDUM=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,IAN,IEN,NLSRFX,NLSRFY,
     .                           NPANU,'STDCOL      ')
            case (2)
              IPOLG=MIN0(NP2NDM,MPSURF)
              IAN=MPSURF
              IEN=MPSURF
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IDUM,IAN,IEN,NLSRFX,NLSRFY,
     .                             NPANU,'STDCOL      ')
            case default
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,
     .                     .FALSE.,.FALSE.,NPANU,'STDCOL      ')
            end select

C  FIND NTCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=EIRENE_LEARCA(Z0,ZSURF,1,NT3RD,1,'STDCOL    ')
              ELSEIF (NLTRA) THEN
                NTCELL=EIRENE_LEARCA(PHI,ZSURF,1,NT3RD,1,'STDCOL    ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF

C  FIND NPCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLPOL) THEN
              select case (LEVGEO)
              case (1)
                NPCELL=EIRENE_LEARCA(Y0,PSURF,1,NP2ND,1,'STDCOL')
              case (2)
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=EIRENE_LEARCA(WINK,PSURF,1,NP2ND,1,'STDCOL')
                ELSE
                  NPCELL=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'STDCOL')
                ENDIF
              case (3)
                NPCELL=IPOLG
              case default
                WRITE (iunout,*) 'ERROR EXIT IN STDCOL, NLPOL ',LEVGEO
                CALL EIRENE_EXIT_OWN(1)
              end select
            ELSE
              NPCELL=1
            ENDIF
          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
C  ENTRANCE INTO STANDARD MESH, BLOCK ILBLCK
C  OR
C  EXIT FROM STANDARD MESH, INTO NACELL=NBLOCK+ILACLL
        ELSEIF (ISWICH(6,MSURF).NE.0) THEN
          IF (NACELL.EQ.0) THEN
C  SET CELL INDEX EQUAL TO ILACLL
            NACELL=NBLOCK+ICOS*ISWICH(6,MSURF)*ILACLL(MSURF)
            NBLOCK=NBMLTP
C
            NRCELL=0
            NPCELL=1
            NTCELL=1
            IF (.NOT.NLADD.OR.NACELL.GT.NRADD.OR.NACELL.LT.1) THEN
              IWEI=-10
              CALL EIRENE_STDCOL_ABSORB
              IRET = IRET_OUT
              RETURN
            ENDIF
          ELSEIF (NACELL.GT.0) THEN
C  ENTRANCE INTO STANDARD MESH, INTO NBLOCK=NACELL+ILBLCK
            NBLOCK=NACELL+ICOS*ISWICH(6,MSURF)*ILBLCK(MSURF)
            NACELL=0
C  FIND  NRCELL,IPOLG IN STANDARD MESH, BLOCK NBLOCK
            select case (idimm)
            case (1)
              NRCELL=MIN0(NR1STM,MRSURF)
              IAN=MRSURF
              IEN=MRSURF
              NDUM=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,IAN,IEN,NLSRFX,NLSRFY,
     .                           NPANU,'STDCOL      ')
            case (2)
              IPOLG=MIN0(NP2NDM,MPSURF)
              IAN=MPSURF
              IEN=MPSURF
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IDUM,IAN,IEN,NLSRFX,NLSRFY,
     .                             NPANU,'STDCOL      ')
            case default
              NRCELL=EIRENE_LEARC1(X0,Y0,Z0,IPOLG,1,NR1STM,
     .                     .FALSE.,.FALSE.,NPANU,'STDCOL      ')
            end select

C  FIND NTCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLTOR) THEN
              IF (NLTRZ) THEN
                NTCELL=EIRENE_LEARCA(Z0,ZSURF,1,NT3RD,1,'STDCOL    ')
              ELSEIF (NLTRA) THEN
                NTCELL=EIRENE_LEARCA(PHI,ZSURF,1,NT3RD,1,'STDCOL    ')
              ENDIF
            ELSE
              NTCELL=1
            ENDIF

C  FIND NPCELL IN STANDARD MESH, BLOCK NBLOCK
            IF (NLPOL) THEN
              select case (LEVGEO)
              case (1)
                NPCELL=EIRENE_LEARCA(Y0,PSURF,1,NP2ND,1,'STDCOL')
              case (2)
                IF (NLCRC) THEN
                  WINK=MOD(ATAN2(Y0,X0)+PI2A-PSURF(1),PI2A)+PSURF(1)
                  NPCELL=EIRENE_LEARCA(WINK,PSURF,1,NP2ND,1,'STDCOL')
                ELSE
                  NPCELL=EIRENE_LEARC2(X0,Y0,NRCELL,NPANU,'STDCOL')
                ENDIF
              case (3)
                NPCELL=IPOLG
              case default
                WRITE (iunout,*) 'ERROR EXIT IN STDCOL, NLPOL ',LEVGEO
                CALL EIRENE_EXIT_OWN(1)
              end select
            ELSE
              NPCELL=1
            ENDIF

          ENDIF
          NBLCKA=NSTRD*(NBLOCK-1)+NACELL
        ENDIF
      ENDIF
C
C  SWITCHING DONE
C
      IF (IWEI.LT.0) then
         CALL EIRENE_STDCOL_ABSORB
         IRET = IRET_OUT
         RETURN
      ENDIF
      IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
        CALL EIRENE_CHCTRC(X0,Y0,Z0,16,8)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
      ENDIF
C
C
      IF (ILIIN(MSURF).LT.0) THEN
        IF (ILIIN(MSURF).EQ.-1) then
          colflag = .true.
          IRET = 1
          RETURN
        endif
        IRET = 2
        RETURN
      ENDIF
C
C  ILIIN(MSURF) .GT. 0, AND  ILIIN(MSURF) .NE. 2
C  PREPARE REFLECTION, OR CARRY OUT PERIODICITY.
C  I.E. SET OUTER NORMAL
C
      select case (IDIMM)
      case (1)
         CALL EIRENE_STDCOL_X_RAD
      case (2)
         CALL EIRENE_STDCOL_Y_POL
      case (3)
         CALL EIRENE_STDCOL_Z_TOR
      end select
      IRET = IRET_OUT
      RETURN
C
  999 CONTINUE
      WRITE (iunout,*) 'ERROR IN STDCOL, WRONG CELL SWITCHING '
      WRITE (iunout,*) 'NPANU, MSURF = ',NPANU,MSURF
      CALL EIRENE_EXIT_OWN(1)
      RETURN
      END SUBROUTINE EIRENE_STDCOL

      SUBROUTINE EIRENE_STDNOR (X0E,Y0E,Z0E,IDIMM_IN,SCOSE,MSURFE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IDIMM_IN
      REAL(DP) :: X0E, Y0E, Z0E, SCOSE
      INTEGER :: MSURFE
      EXTERNAL :: EIRENE_EXIT_OWN
c
c  a) find (outer) surface normal vector at point of intersection X0E,Y0E,Z0E
c     on non-default standard surface MSURFE, (e.g. for surface reflection routines)
c  b) if iliin.gt.3: this is a periodicity surface.
c
c     In this case additionally change particle position, cell number,
c     and velocity of particle accordingly.
c     Note: this periodicity is ready only for some geometry options
c           and surfaces
c
c
c  idimm=1: x or radial surface
c  idimm=2: y or poloidal surface
c  idimm=3: z or toroidal surface
c  idimm=4: mixed x or y grid surface

      IDIMM = IDIMM_IN
      X0=X0E
      Y0=Y0E
      Z0=Z0E
      SCOS=SCOSE
      MSURF=MSURFE
      IRCELL=NRCELL
      IPCELL=NPCELL
      ITCELL=NTCELL
C     ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)

C  IN PARTICULAR BECAUSE OF PERIODICITY.
C     DISTINGUISH NEUTRAL AND IONISED PARTICLES

C
      select case (IDIMM)
      case (1)
         CALL EIRENE_STDCOL_X_RAD
      case (2)
         CALL EIRENE_STDCOL_Y_POL
      case (3)
         CALL EIRENE_STDCOL_Z_TOR
      case (4)
         CALL EIRENE_STDCOL_XYMIX
      case default
        WRITE (iunout,*) 'ERROR EXIT IN STDNOR, IDIMM ',IDIMM
        CALL EIRENE_EXIT_OWN(1)
      end select

      END SUBROUTINE EIRENE_STDNOR


      SUBROUTINE EIRENE_STDCOL_X_RAD
C
C  RADIAL (or X- ) SURFACE, ALSO: SURFACES IN TRIANGULAR,
C                           TETRAHEDRAL AND GENERAL (usr) GRIDS
      IMPLICIT NONE
      REAL(DP) :: TANPHI, PHINM, VCOS, VELS
      INTEGER :: IP, IST
      EXTERNAL :: EIRENE_NORUSR, EIRENE_NEWFIELD, EIRENE_EXIT_OWN,
     .            EIRENE_LEER, EIRENE_MASAGE
C
      select case (LEVGEO)
      case (1)
c  cartesian x coordinate
        CRTX=SCOS
        CRTY=0.
        CRTZ=0.
c
C  PERIODICITY SURFACE IN X DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          MRSURF=EIRENE_IDEZ(ILIIN(MSURF),2,2)
C  NEW X0, AND KEEP Y AND Z (NLTRZ) OR Y AND PHI (NLTRA) CONSTANT
          IF (NLTRA) THEN
            TANPHI=Z0/(X0+RMTOR)
            X0=RSURF(MRSURF)
            Z0=TANPHI*(X0+RMTOR)
          ELSEIF (NLTRZ) THEN
            X0=RSURF(MRSURF)
          ELSEIF (NLTRT) THEN
            WRITE (iunout,*) 'EXIT IN STDCOL_X_RAD'
            CALL EIRENE_EXIT_OWN(1)
          ENDIF
C  NEW CELL NUMBERS
          IF (SCOS.GT.0) NRCELL=MRSURF
          IF (SCOS.LT.0) NRCELL=MRSURF-1
          CALL EIRENE_STDCOL_ASS
          IRET_OUT = 2
          RETURN
        ENDIF

      case (2)
c  radial coordinate in polar-toroidal- system
        CRTX=(X0-EP1(MRSURF))*ELLQ(MRSURF)
        CRTY=Y0
        PHINM=SQRT(CRTX*CRTX+CRTY*CRTY)
        CRTX=CRTX/PHINM*SCOS
        CRTY=CRTY/PHINM*SCOS
        CRTZ=0.
c
C  PERIODICITY SURFACE IN RADIAL DIRECTION, LEVGEO=2
        IF (ILIIN(MSURF).GT.3) THEN
          GOTO 998  ! OPTION NOT AVAILABLE
        ENDIF

      case (3)
c  1st ("radial") coordinate in 2D polygonal grid
        IF (NLPOL)      IP=IPCELL
        IF (.NOT.NLPOL) IP=IPOLG
        CRTX=PLNX(MRSURF,IP)*SCOS
        CRTY=PLNY(MRSURF,IP)*SCOS
        CRTZ=0.
c
C  PERIODICITY SURFACE IN RADIAL DIRECTION, POLYGON GRID
        IF (ILIIN(MSURF).GT.3) THEN
          GOTO 998  ! OPTION NOT AVAILABLE
        ENDIF

      case (4)
c  triangular grid
        IP=IPOLG
        CRTX=PTRIX(IP,MRSURF)*SCOS
        CRTY=PTRIY(IP,MRSURF)*SCOS
        CRTZ=0.
c
C  PERIODICITY SURFACE IN TRIANGULAR GRID
        IF (ILIIN(MSURF).GT.3) THEN
          GOTO 998  ! OPTION NOT AVAILABLE
        ENDIF

      case (5)
c  tetrahedral grid
        IP=IPOLG
        CRTX=PTETX(IP,MRSURF)*SCOS
        CRTY=PTETY(IP,MRSURF)*SCOS
        CRTZ=PTETZ(IP,MRSURF)*SCOS
c
C  PERIODICITY SURFACE IN TETRAHEDRAL GRID
        IF (iliin(msurf).gt.3) THEN
          IST=MSURF-NLIM
          CALL EIRENE_NORUSR(IST,X0,Y0,Z0,CRTX,CRTY,CRTZ,SCOS,
     .                           VELX,VELY,VELZ,NRCELL)
        ENDIF

      case (10)
C
C  GENERAL GEOMETRY OPTION: PROVIDE OUTER SURFACE NORMAL UNIT VECTOR
C                           CRTX,CRTY,CRTZ
        IST=MSURF-NLIM
        CALL EIRENE_NORUSR(IST,X0,Y0,Z0,CRTX,CRTY,CRTZ,SCOS,
     .                     VELX,VELY,VELZ,NDUM)
C  IN CASE OF PERIODICITY: ALSO NEW POSITION, SPEED, SURFACE- AND CELL NUMBERS
        if (iliin(msurf).gt.3) then
          NRCELL=NDUM
C  IN CASE OF PERIODICITY: ALSO NEW POSITION, SPEED, SURFACE- AND CELL NUMBERS
C  NORUSR MAY SET NEW POINT, SURFACE NORMAL AND NEW VELOCITY

C  IN CASE OF TEST IONS: AT THIS POINT VELX,VELZ,VELX ARE THE REDUCED SPEED UNIT VECTOR COMPONENTS
C  E.G. THE PARALLEL TO B OR THE GUIDING CENTRE VELOCITY
          IF (ITYP.EQ.3) THEN
C  IN CASE OF TEST IONS: FIND NEW SIGPAR, LOCAL B FIELD, RETAIN REDUCED SPEEDS VELPER, VELPAR
C  FIND B FIELD IN CELL NCELL
!pb VELS is not used in NEWFIELD with option 0
!pb but for the sake of decent programming set VELS
            VELS = VEL
            CALL EIRENE_NEWFIELD(X0,Y0,Z0,VELS,0)
C  SIGPAR: SIGN OF PARALLEL VELOCITY WITH RESPECT TO B
            VCOS=VELX*BX(0)+VELY*BY(0)+VELZ*BZ(0)
            IF (ABS(VCOS).LT.EPS30) GOTO 992
            SIGPAR=SIGN(1._DP,VCOS)
C  VELPER, VELPAR: AS BEFORE, BUT SPEED UNIT VECTOR:
C  NOW REDUCED VELOCITY: GUIDING CENTRE APPROXIMATION

C  APPROXIMATION A)
C  USE B FIELD LINE AS TRAJECTORY
            VLXPAR=SIGPAR*BX(0)
            VLYPAR=SIGPAR*BY(0)
            VLZPAR=SIGPAR*BZ(0)
C  VL_PAR: PARALLEL UNIT SPEED VECTOR, VL_PAR=SIG*B
C     VL_PAR=(/VLXPAR,VLYPAR,VLZPAR/)
            VELX=VLXPAR
            VELY=VLYPAR
            VELZ=VLZPAR
            VEL =VELPAR
            LCART=.FALSE.
          ENDIF  ! TRACE ION CORRECTIONS AT PERIODICITY SURFACE: DONE
        ENDIF
      case default
       WRITE (iunout,*)
     .  'STDCOL IDIMM=1, BUT LEVGEO OUT OF RANGE. CALL EIRENE_EXIT '
        CALL EIRENE_EXIT_OWN(1)
      end select

      IRET_OUT = 2
      RETURN  ! CONTINUE IN CALLING PROGRAM (FOLNEUT, FOLION)
              ! AT "PARTICLE ESCAPE TO SURFACE"

  992 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('STDCOL IDIMM=1, PROJECTION TO V_PAR, V_PERP    ')
      CALL EIRENE_MASAGE
     .  ('PROBABLY ILL-DEFINED B FIELD WRT. PARTICLE SPEED')
      WRITE (iunout,*) 'BX(0),BY(0),BZ(0) ',BX(0),BY(0),BZ(0)
      IRET_OUT = 0
      RETURN
  998 CONTINUE
      CALL EIRENE_LEER(1)
      CALL EIRENE_MASAGE
     .  ('STDCOL IDIMM=1, UNWRITTEN PERIODICITY OPTION   ')
      WRITE (iunout,*) 'MSURF ',MSURF
      IRET_OUT = 0
      RETURN

      END SUBROUTINE EIRENE_STDCOL_X_RAD


      SUBROUTINE EIRENE_STDCOL_Y_POL
      IMPLICIT NONE
      REAL(DP) :: DST0, DSTT, FR, SINPHI, COSPHI
      INTEGER :: IR, MP_NEW
      EXTERNAL :: EIRENE_EXIT_OWN

C  POLOIDAL SURFACE  MPSURF
C
      select case (LEVGEO)
      case (1)
        CRTX=0.
        CRTY=SCOS
        CRTZ=0.
c
C  PERIODICITY SURFACE IN Y DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          MPSURF=EIRENE_IDEZ(ILIIN(MSURF),2,2)
          Y0=PSURF(MPSURF)
C  NEW CELL NUMBERS, NEW SURFACE NUMBER
          IF (SCOS.GT.0) NPCELL=MPSURF
          IF (SCOS.LT.0) NPCELL=MPSURF-1
          CALL EIRENE_STDCOL_ASS
          IRET_OUT = 2
          RETURN
        ENDIF

C  PERIODICITY SURFACE IN POLOIDAL DIRECTION:
C  USE POLYGON GRID EVEN IN CASE LEVGEO=2, NLCRC
      case (2:3)
        IF (NLRAD)      IR=IRCELL
C       IF (.NOT.NLRAD) IR=???
C
        IF (ILIIN(MSURF).GT.3) THEN
cdr Bug: in case levgeo=2 and nlcrc: xpol,... are undefined
          MP_NEW=EIRENE_IDEZ(ILIIN(MSURF),2,2)
          DST0 = SQRT((XPOL(IR,MPSURF)-X0)**2+(YPOL(IR,MPSURF)-Y0)**2)
          DSTT = BGLP(IR+1,MPSURF)-BGLP(IR,MPSURF)
          FR = DST0 / DSTT
          X0 = XPOL(IR,MP_NEW) + FR*VVTX(IR,MP_NEW)
          Y0 = YPOL(IR,MP_NEW) + FR*VVTY(IR,MP_NEW)
          CALL EIRENE_STDCOL_ASS
          IRET_OUT = 2

! NEW VELOCITY
! (PPLNX,PPLNY) ARE ALREADY NORMALIZED
          COSPHI=PPLNX(IR,MPSURF)*PPLNX(IR,MP_NEW) +
     .           PPLNY(IR,MPSURF)*PPLNY(IR,MP_NEW)
          SINPHI=SQRT(1._DP-COSPHI*COSPHI)
! COS(-PHI)=COS(PHI) ; SIN(-PHI)=-SIN(PHI)
          COSROT=COSPHI
          SINROT= -SINPHI
          VELX_OLD = VELX
          VELX = VELX_OLD*COSROT - VELY*SINROT
          VELY = VELX_OLD*SINROT + VELY*COSROT
C  NEW SURFACE AND CELL NUMBERS
          MPSURF=MP_NEW
          IF (SCOS.GT.0) THEN
            NPCELL = MPSURF
            IPOLG = MPSURF
            IPOLGN = MPSURF
          END IF
          IF (SCOS.LT.0) THEN
            NPCELL = MPSURF-1
            IPOLG = MPSURF-1
            IPOLGN = MPSURF-1
          END IF
        END IF
!PB USE OUTER NORMAL OF PERIODICITY COUNTERPART
        CRTX=PPLNX(IR,MPSURF)*SCOS
        CRTY=PPLNY(IR,MPSURF)*SCOS
        CRTZ=0.
      case default
        WRITE (iunout,*)
     .  'STDCOL IDIMM=2, BUT LEVGEO.GT.3 CALL EIRENE_EXIT '
        CALL EIRENE_EXIT_OWN(1)
      end select
      IRET_OUT = 2
      RETURN
      END SUBROUTINE EIRENE_STDCOL_Y_POL


      SUBROUTINE EIRENE_STDCOL_Z_TOR
      IMPLICIT NONE
      EXTERNAL :: EIRENE_EXIT_OWN
C
C  TOROIDAL SURFACE
C
      IF (NLTRZ) THEN
        CRTX=0.
        CRTY=0.
        CRTZ=SCOS

C  PERIODICITY SURFACE IN Z-DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          MTSURF=EIRENE_IDEZ(ILIIN(MSURF),2,2)
          Z0=ZSURF(MTSURF)
          IF (SCOS.GT.0) NTCELL=MTSURF
          IF (SCOS.LT.0) NTCELL=MTSURF-1
          CALL EIRENE_STDCOL_ASS
          IRET_OUT = 2
          RETURN
        ENDIF

      ELSEIF (NLTRA) THEN
        CRTX=-SIN(ALPHA)
        CRTY=0.
        CRTZ=SCOS*COS(ALPHA)

C  PERIODICITY SURFACE IN TOROIDAL-DIRECTION
        IF (ILIIN(MSURF).GT.3) THEN
          MTSURF=EIRENE_IDEZ(ILIIN(MSURF),2,2)
C  NEW PHI POSITION. KEEP X0,Y0 FIXED
          PHI=ZSURF(MTSURF)
C  ROTATE VELOCITY BECAUSE OF LOCAL COORDINATE SYSTEM ROTATION
          IF (ITYP.LE.2) THEN
C  ANGLE FOR ROTATION IS ONLY 2*ALPHA, BECAUSE THE REST OF THE TOROIDAL
C  PERIODICITY IS CONTAINED IN ROTATED LOCAL COORDINATE SYSTEMS.
            COSROT=COSAL
            SINROT=SINAL*FLOAT(NINCZ)
            VELX_OLD=VELX
            VELX= VELX_OLD*COSROT +VELZ*SINROT
            VELZ=-VELX_OLD*SINROT +VELZ*COSROT
          ENDIF
          Z01=-Z01
          Z0=Z01
C  NEW CELL NUMBERS
          IF (SCOS.GT.0) NTCELL=MTSURF
          IF (SCOS.LT.0) NTCELL=MTSURF-1
          IPERID=NTCELL
          CALL EIRENE_STDCOL_ASS
          IRET_OUT = 2
          RETURN
        ELSE
C         CRTX=?
          CRTY=0.
C         CRTZ=?
          WRITE (iunout,*)
     .      'STDCOL IDIMM=3, BUT NOT PERIODIC. CALL EIRENE_EXIT '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ELSE
        WRITE (iunout,*)
     .    'STDCOL IDIMM=3, BUT NEITHER NLTRZ NOR NLTRA: EXIT '
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
      IRET_OUT = 2
      RETURN
      END SUBROUTINE EIRENE_STDCOL_Z_TOR


      SUBROUTINE EIRENE_STDCOL_XYMIX
      IMPLICIT NONE
      EXTERNAL :: EIRENE_EXIT_OWN
C
C  MIXED X-Y GRID SURFACE,  idimm=4
C
      IF (NLSRFX) CALL EIRENE_STDCOL_X_RAD
      IF (NLSRFY) CALL EIRENE_STDCOL_Y_POL
      IF (NLSRFX .OR. NLSRFY) RETURN
      WRITE (iunout,*) ' ERROR IN STDCOL_XY-MIX,',
     .            ' IDIMM=4 AND NEITHER NLSRFX NOR NLSRFY '
      CALL EIRENE_EXIT_OWN(1)
C
C
      END SUBROUTINE EIRENE_STDCOL_XYMIX


      SUBROUTINE EIRENE_STDCOL_ABSORB
      IMPLICIT NONE
      INTEGER :: MSURFO
      EXTERNAL :: EIRENE_MASJ1, EIRENE_MASR3
C
C  IWEI.LT.0, I.E., ILIIN OPTION IS OVERRULED FROM THIS SIDE
C
      IF (IWEI.EQ.-1) THEN
C  PARTICLE HAS HIT A SURFACE FROM AN ABSORBING SIDE
C  SCORE FLUXES (DO NOT SET WEIGHT=0.D0, IRET_OUT=2)
C  AND ABSORB PARTICLE LATER
        IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,8)
          WRITE (iunout,*) 'ABSORB PARTICLE: NPANU ',NPANU
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        ENDIF
        IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WEIGHT
          ENDIF
        END IF
        IF (LSPUMP) LMETSPW(ISPZ) = .TRUE.
        LGPART=.FALSE.
        IRET_OUT = 2
        RETURN
      ELSEIF (IWEI.EQ.-2) THEN
C  KILL THIS PARTICLE BECAUSE IT COMES FROM WRONG SIDE
C  DO NOT SCORE FLUXES (SET WEIGHT=0.D0, IRET_OUT=2)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
        IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        ENDIF
        WRITE (iunout,*) 'ERROR DETECTED IN SUBR. STDCOL'
        WRITE (iunout,*) 'PARTICLE COMES FROM WRONG SIDE'
        CALL EIRENE_MASJ1 ('NPANU=  ',NPANU)
        MSURFO=MSURF
        IF (MSURFO.GT.NLIM) MSURFO=-(MSURFO-NLIM)
        CALL EIRENE_MASJ1 ('MSURF NW',MSURFO)
        IF (MSURFS.GT.NLIM) MSURFS=-(MSURFS-NLIM)
        CALL EIRENE_MASJ1 ('MSURF OD',MSURFS)
        CALL EIRENE_MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL EIRENE_MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL EIRENE_MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL EIRENE_MASR3 ('VEL,WEIGHT,E0           ',VEL,WEIGHT,E0)
        IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WEIGHT
          ENDIF
        END IF
        IF (LSPUMP) LMETSPW(ISPZ) = .TRUE.
        WEIGHT=0.
        LGPART=.FALSE.
        IRET_OUT = 2
        RETURN
      ELSEIF (IWEI.EQ.-3) THEN
C  SURFACE IS NOT SEEN BY THE PARTICLE BECAUSE OF ILSIDE OPTION
C  I.E., SURFACE IS TRANSPARENT FROM THIS SIDE
C  ACTS AS ILIIN=0 OPTION (NO SURFACE TALLIES, NO SWITCHES)
!pb for levgeo=4 or levgeo=5 increasing of cell numbers is wrong
        IF (LEVGEO <= 3) THEN
          IF (IDIMM.EQ.1.AND.ILIIN(MSURF).GT.0) NRCELL=NRCELL+ICOS
          IF (IDIMM.EQ.2.AND.ILIIN(MSURF).GT.0) NPCELL=NPCELL+ICOS
          IF (IDIMM.EQ.3.AND.ILIIN(MSURF).GT.0) NTCELL=NTCELL+ICOS
        END IF
        IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,8)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        ENDIF
        colflag = .true.
        IRET_OUT = 1
        RETURN
      ELSEIF (IWEI.EQ.-10) THEN
C  KILL THIS PARTICLE BECAUSE CELL NUMBER OUT OF RANGE DUE TO SWITCHING
C  DO NOT SCORE FLUXES (SET WEIGHT=0.D0)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
        PTRASH(ISTRA)=PTRASH(ISTRA)-WEIGHT
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
        ETRASH(ISTRA)=ETRASH(ISTRA)-WEIGHT*E0
        IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
          CALL EIRENE_CHCTRC(X0,Y0,Z0,16,18)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
        ENDIF
        WRITE (iunout,*) 'ERROR DETECTED IN SUBR. STDCOL'
        WRITE (iunout,*) 'NACELL OUT OF RANGE'
        CALL EIRENE_MASJ1 ('NPANU=  ',NPANU)
        WRITE (iunout,*) 'NLMLT,NLADD ',NLMLT,NLADD
        WRITE (iunout,*) 'NBMLT,NRADD ',NBMLT,NRADD
        MSURFO=MSURF
        IF (MSURFO.GT.NLIM) MSURFO=-(MSURFO-NLIM)
        CALL EIRENE_MASJ1 ('MSURF NW',MSURFO)
        IF (MSURFS.GT.NLIM) MSURFS=-(MSURFS-NLIM)
        CALL EIRENE_MASJ1 ('MSURF OD',MSURFS)
        CALL EIRENE_MASJ1 ('NACL NEW',NACELL)
        CALL EIRENE_MASJ1 ('NACL OLD',NACLLS)
        CALL EIRENE_MASR3 ('X0,Y0,Z0 (NEW)          ',X0,Y0,Z0)
        CALL EIRENE_MASR3 ('X0,Y0,Z0 (OLD)          ',X0SA,Y0SA,Z0SA)
        CALL EIRENE_MASR3 ('VELX,VELY,VELZ          ',VELX,VELY,VELZ)
        CALL EIRENE_MASR3 ('VEL,WEIGHT,E0           ',VEL,WEIGHT,E0)
        IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
          SPUMP(ISPZ,MSURF)=SPUMP(ISPZ,MSURF)+WEIGHT
        ENDIF
        IF (MSURFG.GT.0) THEN
          IF (LSPUMP) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC
#endif
            SPUMP(ISPZ,MSURFG)=SPUMP(ISPZ,MSURFG)+WEIGHT
          ENDIF
        END IF
        IF (LSPUMP) LMETSPW(ISPZ) = .TRUE.
        WEIGHT=0.
        LGPART=.FALSE.
        IRET_OUT = 2
        RETURN
      ENDIF
      END SUBROUTINE EIRENE_STDCOL_ABSORB


      SUBROUTINE EIRENE_STDCOL_ASS
c  set further assistant coordinates for curvilinear 2nd and or 3rd grid
c        x01,z01      (for tracing in toroidal grid)
c        x00,y00,z00  (for tracing in poloidal grid)

      IMPLICIT NONE

      X00=X0
      Y00=Y0
      Z00=Z0
      IF (NLTRA) X01=X0+RMTOR
      Z01=Z0
      END SUBROUTINE EIRENE_STDCOL_ASS

      END MODULE EIRMOD_STDCOL
