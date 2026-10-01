c  nov 03:  use relative distances to find neighbor segment,
c           otherwise sometimes problems with non-closing polygons encountered.
cdr june 17:  separate WRMESH (WRITING) and PLMESH (PLOTTING).
cdr Nov. 18:  fixes from ITER branch

      SUBROUTINE EIRENE_WRMESH
c  create closed polygonal contours, from the eirene standard and
c  additional surfaces
c  use ILPLG(isurf) flag, from input blocks 3A LEVGEO=3 OR LEVGEO=4,
C                         or certain additional surfaces, input block 3B,
c                         0<RLB<2.
c  This set of closed contours, together with their orientation, can be used
c  in 2D grid generators to produce multiply connected triangular grids.
c  The orientation indicates whether the inner or outer part of a closed contour
c  is a valid computational volume for triangular grid generation (not needed
c  for plotting)


c  EIRENE_WRMESH: Write closed contours onto output stream 78+ifoff.
c  EIRENE_PLMESH: plots these contours, using GR plot software.


      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CADGEO
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CINIT, ONLY: FORT_LC
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGEOM
      USE EIRMOD_CLGIN
      USE EIRMOD_CGRPTL
      USE EIRMOD_CTRIG
      USE EIRMOD_CGRID
      USE EIRMOD_CTRCEI
      IMPLICIT NONE


      REAL(DP), ALLOCATABLE :: partcont(:,:,:)
      INTEGER, ALLOCATABLE :: IDIAG(:), irip(:,:)
      INTEGER, SAVE :: MAXPOIN=2000
      REAL(DP) :: XPE, YPE, HELP, XT, YT, PHI1, X1, X2, Y1, Y2, PHI2
      REAL(DP) :: DISTQI, DISTQJ1, DISTQJ2, YMN, maxlen
      INTEGER  :: ICONT, IPOIN, IWST, IWEN, IWL, IWP,
     .            IWAN, IMN, I, NCONT, J, IUHR, ISTORE, IP, IH, IFOUND,
     .            ICO, IPO, IN, IS, IS1, ITRI, INBT, INBS
      LOGICAL  :: LFOUND, FORT78EXISTS
      LOGICAL, ALLOCATABLE :: FOUND(:,:)
      EXTERNAL :: EIRENE_LEER

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT BLOCK 3 EINGELESEN
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'SUBROUTINE WRMESH CALLED'
      CALL EIRENE_LEER(1)

      ALLOCATE (partcont(maxpoin,2,2))
      ALLOCATE (idiag(maxpoin))
      ALLOCATE (irip(maxpoin,2))

      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      IF (TRCSUR) THEN
        WRITE (iunout,*) 'NUMBER OF CONTOURS FOR FEM MESH: ',NCONT
        CALL EIRENE_LEER(1)
      END IF

      if (ncont == 0) then
        call EIRENE_leer(1)
        write (iunout,*) 'No contours specified in blocks 3a,3b'
        write (iunout,*)
     .    'No input file '//fort_lc//
     .    '78 for FEM mesh generator written'
        call EIRENE_leer(2)
        return
      endif
      fort78exists = .false.
      inquire(file='fort.78', exist=fort78exists)
      if (fort78exists) then 
        write(iunout, *) 'EIRENE_WRMESH: fort.78 file already ' // 
     .   'exists, exiting...'
        return 
      end if 

      IF (.NOT.ALLOCATED(NCONPOINT)) THEN
        ALLOCATE (NCONPOINT(NCONT))
        ALLOCATE (XCONTOUR(MAXPOIN,NCONT))
        ALLOCATE (YCONTOUR(MAXPOIN,NCONT))
      ENDIF
      NCONPOINT = 0
      ICO = 0

      DO ICONT = 1,NCONT
        IPOIN = 0
        MAXLEN = 0.
        irip=0
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN KONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT

        select case (LEVGEO)
        case (3)
c  ADDITIONAL SURFACES
          DO I=1,NLIMI
            IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
              IUHR=ILPLG(I)
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
              IF ((RLB(I) .GT. 0.) .AND. (RLB(I) .LT. 2.) .AND.
     >            (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >            .OR. P3(3,I) .EQ. 1.D55)) THEN
                IPOIN = IPOIN + 1
                IF (IPOIN.GT.MAXPOIN) THEN
                  WRITE(IUNOUT,*)
     .             'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .              ICONT
                  WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                  CALL EIRENE_EXTEND_ARRAYS
                  WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                ENDIF
                IF (A3LM(I) .EQ. 0.) THEN
C               X,Y-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(1,I)
                  PARTCONT(IPOIN,1,2) = P1(2,I)
                  PARTCONT(IPOIN,2,1) = P2(1,I)
                  PARTCONT(IPOIN,2,2) = P2(2,I)
                  idiag(ipoin)=i
                ELSEIF (A2LM(I) .EQ. 0.) THEN
C               X,Z-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(1,I)
                  PARTCONT(IPOIN,1,2) = P1(3,I)
                  PARTCONT(IPOIN,2,1) = P2(1,I)
                  PARTCONT(IPOIN,2,2) = P2(3,I)
                  idiag(ipoin)=i
                ELSEIF (A1LM(I) .EQ. 0.) THEN
C               Y,Z-KOORDINATEN
                  PARTCONT(IPOIN,1,1) = P1(2,I)
                  PARTCONT(IPOIN,1,2) = P1(3,I)
                  PARTCONT(IPOIN,2,1) = P2(2,I)
                  PARTCONT(IPOIN,2,2) = P2(3,I)
                  idiag(ipoin)=i
                ENDIF
                maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                IF (TRCGRD)
     >           WRITE(iunout,
     >            '(a,g14.7,a,g14.7,a,g14.7,a,g14.7,a,i4)')
     >            'Grabbed segment (',
     >            PARTCONT(IPOIN,1,1),',',PARTCONT(IPOIN,1,2),
     >            ') to (',
     >            PARTCONT(IPOIN,2,1),',',PARTCONT(IPOIN,2,2),
     >            ') from wall ',I
              ELSE
C  ERROR
                WRITE(iunout,'(a,f11.4,2i4)')
     >           'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),ILPLG(I),I
              ENDIF
            ENDIF
          ENDDO

          DO I=1,NSTSI
            IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
              IUHR=ILPLG(NLIM+I)
              IF (INUMP(I,2) .NE. 0) THEN
C  POLOIDAL SURFACES
                DO J=IRPTA(I,1),IRPTE(I,1)-1
                  IF ((XPOL(J,INUMP(I,2)).NE.XPOL(J+1,INUMP(I,2))) .OR.
     >                (YPOL(J,INUMP(I,2)).NE.YPOL(J+1,INUMP(I,2)))) THEN
                    IPOIN = IPOIN + 1
                    IF (IPOIN.GT.MAXPOIN) THEN
                      WRITE(IUNOUT,*)
     .                 'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                  ICONT
                      WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                      CALL EIRENE_EXTEND_ARRAYS
                      WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                    ENDIF
                    PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                    PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                    PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                    PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                    idiag(ipoin)=-i
                    irip(ipoin,1)=j
                    irip(ipoin,2)=INUMP(I,2)
                    maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                  ENDIF
                ENDDO
              ELSEIF (INUMP(I,1) .NE. 0) THEN
C  RADIAL SURFACES
                DO J=IRPTA(I,2),IRPTE(I,2)-1
                  IF ((XPOL(INUMP(I,1),J).NE.XPOL(INUMP(I,1),J+1)) .OR.
     >                (YPOL(INUMP(I,1),J).NE.YPOL(INUMP(I,1),J+1))) THEN
                    IPOIN = IPOIN + 1
                    IF (IPOIN.GT.MAXPOIN) THEN
                      WRITE(IUNOUT,*)
     .                 'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                  ICONT
                      WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                      CALL EIRENE_EXTEND_ARRAYS
                      WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                    ENDIF
                    PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                    PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                    PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                    PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                    idiag(ipoin)=-i
                    irip(ipoin,1)=INUMP(I,1)
                    irip(ipoin,2)=j
                    maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                  ENDIF
                ENDDO
              ELSE
C  ERROR
                WRITE(iunout,*) 'CASE NOT FORESEEN: INUMP: ',
     >                     (INUMP(I,J),J=1,3)
              ENDIF
            ENDIF
          ENDDO

        case (4)
C  TRIANGLE SIDES
          ALLOCATE (FOUND(1:3,1:NTRII))
          FOUND = .FALSE.
          DO ITRI = 1, NTRII
            DO IS = 1, 3
              IN=INMTI(IS,ITRI)
              LFOUND=FOUND(IS,ITRI)
              IF (IN /= 0 .AND. .NOT.LFOUND) THEN
                IF (ABS(ILPLG(IN)) == ICONT) THEN
                  IUHR=ILPLG(IN)
                  IS1 = IS+1
                  IF (IS1 > 3) IS1=1
                  IPOIN = IPOIN + 1
                  IF (IPOIN.GT.MAXPOIN) THEN
                    WRITE(IUNOUT,*)
     .               'INSUFFICIENT NUMBER OF POINTS FOR CONTOUR ',
     .                ICONT
                    WRITE(IUNOUT,*) 'INCREASE VALUE OF MAXPOIN'
                    CALL EIRENE_EXTEND_ARRAYS
                    WRITE(IUNOUT,*) 'MAXPOIN SET TO = ', MAXPOIN
                  ENDIF
                  PARTCONT(IPOIN,1,1) = XTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,1,2) = YTRIAN(NECKE(IS,ITRI))
                  PARTCONT(IPOIN,2,1) = XTRIAN(NECKE(IS1,ITRI))
                  PARTCONT(IPOIN,2,2) = YTRIAN(NECKE(IS1,ITRI))
                  idiag(ipoin)=IN
                  irip(ipoin,1)=itri
                  irip(ipoin,2)=is
                  maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                  FOUND(IS,ITRI)=.TRUE.
                  INBT=NCHBAR(IS,ITRI)
cwdk Make sure the corresponding side of the neighboring triangle is not found again
                  IF (INBT.GT.0) THEN
                    INBS=NSEITE(IS,ITRI)
                    FOUND(INBS,INBT)=.TRUE.
                  ENDIF
                END IF
              END IF
            END DO
          END DO
          DEALLOCATE (FOUND)
        end select

        IF (IPOIN.LE.0) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' NOT FOUND'
          CALL EIRENE_LEER(1)
          GOTO 1000
        ENDIF

C STUECKE DER AKTUELLEN KONTOUR WERDEN SORTIERT
        DO I=1,IPOIN-1
          XPE = PARTCONT(I,2,1)
          YPE = PARTCONT(I,2,2)
          DISTQI=(PARTCONT(I,2,1)-PARTCONT(I,1,1))**2+
     .           (PARTCONT(I,2,2)-PARTCONT(I,1,2))**2
          IFOUND=0
          DO J=I+1,IPOIN
            IF (IFOUND.EQ.1) CYCLE
            DISTQJ1=(XPE-PARTCONT(J,1,1))**2+
     .              (YPE-PARTCONT(J,1,2))**2
            DISTQJ2=(XPE-PARTCONT(J,2,1))**2+
     .              (YPE-PARTCONT(J,2,2))**2
            IF (DISTQJ1/DISTQI.LE.1.D-10) THEN
              IFOUND=1
              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ELSEIF (DISTQJ2/DISTQI.LE.1.D-10) THEN

              IFOUND=1
              HELP = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih

              ih=irip(i+1,1)
              irip(i+1,1)=irip(j,1)
              irip(j,1)=ih
              ih=irip(i+1,2)
              irip(i+1,2)=irip(j,2)
              irip(j,2)=ih
            ENDIF
          ENDDO
          IF (IFOUND.EQ.0) THEN
            WRITE(iunout,*) 'NO MATCHING POINT FOUND FOR CONTOUR ',
     >                       ICONT
            write(iunout,'(I4,3(1X,I4),4(1X,G14.7))') i,idiag(i),
     >                    irip(i,1),irip(i,2),
     >                    partcont(i,1,1),partcont(i,1,2),
     >                    partcont(i,2,1),partcont(i,2,2)
            WRITE(iunout,*) 'USE NEXT POINT'
            IP=I+1
            write(iunout,'(I4,3(1X,I4),4(1X,G14.7))') iP,idiag(iP),
     >                    irip(ip,1),irip(ip,2),
     >                    partcont(iP,1,1),partcont(iP,1,2),
     >                    partcont(iP,2,1),partcont(iP,2,2)
            CALL EIRENE_LEER(1)
          ENDIF
        ENDDO

        IF ((PARTCONT(1,1,1) .NE. PARTCONT(IPOIN,2,1)) .OR.
     >      (PARTCONT(1,1,2) .NE. PARTCONT(IPOIN,2,2))) THEN
          WRITE(iunout,*) 'CONTOUR ',ICONT,' IS NOT CLOSED'
        ELSE
          WRITE(iunout,*) 'CLOSED CONTOUR ',ICONT
        ENDIF

        IF (TRCSUR) THEN
          do i=1,ipoin
            write(iunout,'(1X,I4,3(1X,I4),4(1X,G14.7))') i,idiag(i),
     >                   irip(i,1),irip(i,2),
     >                   partcont(i,1,1),partcont(i,1,2),
     >                   partcont(i,2,1),partcont(i,2,2)
          enddo
          CALL EIRENE_LEER(1)
        END IF  !  contour no. icont done.

C  BERECHNUNG VON DELTA ALS MITTLERE LAENGE DER TEILSTUECKE
C  DELTA IST MASS FUER DIE GROESSE DER DREIECKE

        IF (ICONT .EQ. 1) THEN
          maxlen = maxlen / REAL(IPOIN,DP)
          WRITE(78+ifoff,*) maxlen
          WRITE(78+ifoff,*)
        endif

C BESTIMMUNG DES UHRZEIGERSINNS DER KONTOUR
C KONTOUR MUSS FUER DIE TRIANGULIERUNG FOLGENDERMASSEN AUSGEGEBEN
C WERDEN:
C  - IM UHRZEIGERSINN FUER INNERE BEGRENZUNGEN DES GEBIETES (POSITIV)
C  - GEGEN UHRZEIGERSINN FUER AEUSSERE BEGRENZUNGEN DES GEBIETES (NEGATIV)
        IMN=0
        YMN=PARTCONT(1,1,2)
        DO I=1,IPOIN
          IF (PARTCONT(I,2,2) .LT. YMN) THEN
            YMN = PARTCONT(I,2,2)
            IMN = I
          ENDIF
        ENDDO

        IF (IPOIN > 1) THEN
C  SONDERFALL IMN=IPOIN ENTFAELLT, DA ERSTER PUNKT GLEICH LETZTER
C  PUNKT GILT
          IF (IMN .EQ. 0) THEN
            XT = PARTCONT(1,1,1)
            YT = PARTCONT(1,1,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
            X1 = PARTCONT(IPOIN,1,1)
            Y1 = PARTCONT(IPOIN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
            X2 = PARTCONT(1,2,1)
            Y2 = PARTCONT(1,2,2)
          ELSE
            XT = PARTCONT(IMN,2,1)
            YT = PARTCONT(IMN,2,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
            X1 = PARTCONT(IMN,1,1)
            Y1 = PARTCONT(IMN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
            IF (IMN.LT.IPOIN) THEN
              X2 = PARTCONT(IMN+1,2,1)
              Y2 = PARTCONT(IMN+1,2,2)
            ELSE
              X2 = PARTCONT(1,2,1)
              Y2 = PARTCONT(1,2,2)
            ENDIF
          ENDIF

C  BESTIMME POLARWINKEL VON (X1,Y1) UND (X2,Y2) MIT (XT,YT) ALS URSPRUNG
          PHI1 = ATAN2 (Y1-YT,X1-XT)
          PHI2 = ATAN2 (Y2-YT,X2-XT)

          IF (PHI2 .GT. PHI1) THEN
C  ABSPEICHERUNG ERFOLGTE IM UHRZEIGERSINN
            ISTORE = 1
          ELSE
            ISTORE = -1
          ENDIF

        ELSE

          ISTORE = 1

        ENDIF

C  IUHR=ILPLG > 0 ==> IM UHRZEIGERSINN AUSGEBEN
C  IUHR=ILPLG < 0 ==> ENTGEGEN DEM UHRZEIGERSINN AUSGEBEN
        IWAN=1
        IWEN=IPOIN
        IWST=1
        IWP=1
        IWL=2
        IF (ISTORE*IUHR .LT. 0.) THEN
          IWAN=IPOIN
          IWEN=1
          IWST=-1
          IWP=2
          IWL=1
        ENDIF

        WRITE(78+ifoff,*) IPOIN+1
        IF (IUHR > 0) THEN
          ICO=ICO+1
          NCONPOINT(ICO)=IPOIN+1
          IPO=0
        END IF
        DO I=iwan,iwen,iwst
          WRITE(78+ifoff,'(1P,2(2X,E21.14))')
     >          PARTCONT(I,IWP,1),PARTCONT(I,IWP,2)
          IF (IUHR > 0) THEN
            IPO=IPO+1
            XCONTOUR(IPO,ICO) = PARTCONT(I,IWP,1)
            YCONTOUR(IPO,ICO) = PARTCONT(I,IWP,2)
          END IF
        ENDDO
        WRITE(78+ifoff,'(1P,2(2X,E21.14))') PARTCONT(IWEN,IWL,1),
     >                                      PARTCONT(IWEN,IWL,2)
        IF (IUHR > 0) THEN
          IPO=IPO+1
          XCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,1)
          YCONTOUR(IPO,ICO) = PARTCONT(IWEN,IWL,2)
        END IF

 1000   CONTINUE
      ENDDO    ! END OF DO ICONT.... LOOP
      NCONTOUR=ICO

      call EIRENE_leer(1)
      write (iunout,*)
     .  'input file ', FORT_LC, '78 for FEM mesh generator written'
      call EIRENE_leer(2)
      close(78+ifoff)

cdr
      if (allocated(nconpoint)) then
        DEALLOCATE (NCONPOINT)
        DEALLOCATE (XCONTOUR)
        DEALLOCATE (YCONTOUR)
      endif

      DEALLOCATE (partcont)
      DEALLOCATE (idiag, irip)

      RETURN

      CONTAINS

      SUBROUTINE EIRENE_EXTEND_ARRAYS

      IMPLICIT NONE
      REAL(DP), ALLOCATABLE :: pc(:,:,:)
      INTEGER, ALLOCATABLE :: ID(:), ir(:,:)
      INTEGER :: NEWPOIN

      ALLOCATE(PC(MAXPOIN,2,2))
      ALLOCATE(ID(MAXPOIN))
      ALLOCATE(IR(MAXPOIN,2))

      PC = PARTCONT
      ID = IDIAG
      IR = IRIP

      DEALLOCATE (PARTCONT)
      DEALLOCATE (IDIAG)
      DEALLOCATE (IRIP)

      NEWPOIN = MAXPOIN + 2000
      ALLOCATE (partcont(newpoin,2,2))
      ALLOCATE (idiag(newpoin))
      ALLOCATE (irip(newpoin,2))

      partcont(1:maxpoin,:,:) = pc(1:maxpoin,:,:)
      idiag(1:maxpoin) = id(1:maxpoin)
      irip(1:maxpoin,:) = ir(1:maxpoin,:)

      maxpoin = newpoin

      deallocate(pc)
      deallocate(id)
      deallocate(ir)

      return
      END SUBROUTINE EIRENE_EXTEND_ARRAYS

      END SUBROUTINE EIRENE_WRMESH
