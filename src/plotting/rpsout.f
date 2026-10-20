C  ifoff offset for output stream introduced
C
      SUBROUTINE EIRENE_RPSOUT

CDR  produce formatted output files for RAPS graphics system
C
C   options available so far: documentation to be written
C
C   nov. 15: something added for levgeo=1, y-z plot at given x position
CDR
C
C  ANZ: NUMBER OF CELLS
C  WRITE (17,...) LABELED COORDINATES OF VERTICES
C  WRITE (18,...) LABELING NUMBER OF VERTICES FOR EACH CELL
C  WRITE (19,...) VALUE OF TALLY AT VERTEX (BY LABEL)
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CRECH
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPLOT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CLGIN
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CGRPTL
      USE EIRMOD_CCONA
cdr   USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP) :: YWERT(2*NPTAL)
      INTEGER :: NCELL, NCO, IR, IP, IGR, IPPLG, NPUNKT, ISTS, IDIMP,
     .           NST, I, NSTAB, NRPS, IT, IF, IA, IE, ANZ
      integer :: ird1,ird2,ird3,ird4,ird5,ird6,ird7,ird8
      integer :: igroups, ipoints, ipl
      real(DP), ALLOCATABLE :: phelp(:,:),valcont(:), xcont(:), ycont(:)
      real(DP) :: xm, ym, length_m_p, xgeomin, xgeomax,
     .            ygeomin, ygeomax, epsrel, ahelp, vecax, vecay, vecbx,
     .            vecby,length_veca, length_vecb, del_eps,
     .            x1,x2,x3,x4,y1,y2,y3,y4,atri1,atri2
      integer :: icont, ipoint, iloop,ncont, ip_start, idel, ifc
      logical :: del_point

      TYPE(PPOINT), POINTER :: CUR
C
      allocate(valcont(iraps))
      valcont = huge(1._DP)
      igroups = 0
      I=0
      NSTAB=0
      CUR => FIRST_POINT
      IF (ASSOCIATED(CUR)) THEN
        DO WHILE (ASSOCIATED(CUR%NXTPNT))
          IF (CUR%NXTPNT%NPL2D.NE.0) THEN
C  START OF A NEW LINE
            NSTAB=NSTAB+1
          ENDIF
          CUR => CUR%NXTPNT
        END DO
      END IF
C
      IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        NST=NSTSI
        IF (NTIME.GT.0) NST=NST-1
        DO ISTS=1,NST
          IF (ILIIN(NLIM+ISTS) > 0) THEN
            IDIMP=TRANSFER(MAXLOC(INUMP(ISTS,1:3)),1)
            IF (IDIMP == 1) THEN
              IA=IRPTA(ISTS,2)
              IE=IRPTE(ISTS,2)
            ELSEIF (IDIMP == 2) THEN
              IA=IRPTA(ISTS,1)
              IE=IRPTE(ISTS,1)
            END IF
            IF ((IDIMP>0) .AND. (IDIMP<=2))  NSTAB=NSTAB+(IE-IA)
          END IF
        END DO
      END IF
C
      DO 5 IF=1,IRAPS
        NRPS=IUNRAPSVEC+IF
        OPEN (UNIT=NRPS,ACCESS='SEQUENTIAL',FORM='FORMATTED')
        REWIND NRPS
    5 CONTINUE
!pb
!pb find index of first RAPS plot, to be used for scaling with FCABS

      IFC = 1
      DO I = 1, NVOLPL
        IF (LRAPS3(I)) THEN
           IFC = I
           EXIT
        END IF
      END DO
C
C  3D
      if ((levgeo.eq.1.and.nlrad.and.nlpol.and.nltor.and.nltrz) .or.
     .    (levgeo.eq.5.and..not.lrpscut).or.lraps3d) then
         WRITE(17+ifoff,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)')
     .        '-1111','NPCO','1','3','1','1'
C  2D
      else
         WRITE(17+ifoff,'(1X,A5,8X,A4,11X,A1,11X,A1,11X,A1,11X,A1)')
     .        '-1111','NPCO','1','2','1','1'
      endif
C
      WRITE(19+ifoff,'(1X,A5,8X,A4,50(9X,I3))') '-1111',
     .        'NPST',1,IRAPS,1,(1,I=1,IRAPS)
C
      if (levgeo.eq.1.and.nlrad.and.nlpol.and.nltor.and.nltrz
     .    .and..not.lptorr) then
         ANZ=(nr1st-1)*(np2nd-1)*(nt3rd-1)
         i=0
         do ir=1,nr1st
            do ip=1,np2nd
               do it=1,nt3rd
                  DO IF=1,IRAPS
                     READ(IUNRAPSVEC+IF,*) YWERT(IF)
                  enddo
                  i = i+1
                  WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .                  I,(YWERT(IF),IF=1,IRAPS)
                  WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .                  I,rsurf(ir),psurf(ip),zsurf(it)
               enddo
            enddo
         enddo


         WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
         WRITE(18+ifoff,'(A8,I6,I9,I6)') 'HEXE8   ',1,ANZ,8

         do ir=1,nr1st-1
            do ip=1,np2nd-1
               do it=1,nt3rd-1
                  IRD1=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IRD2=IR+1+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IRD3=IR+1+((IP)+(IT-1)*NP2T3)*NR1P2
                  IRD4=IR+((IP)+(IT-1)*NP2T3)*NR1P2
                  IRD5=IR+((IP-1)+(IT)*NP2T3)*NR1P2
                  IRD6=IR+1+((IP-1)+(IT)*NP2T3)*NR1P2
                  IRD7=IR+1+((IP)+(IT)*NP2T3)*NR1P2
                  IRD8=IR+((IP)+(IT)*NP2T3)*NR1P2
                  WRITE(18+ifoff,'(1X,A1,8I10)')'0',ird1,ird2,ird3,ird4,
     .                                              ird5,ird6,ird7,ird8
               enddo
            enddo
         enddo

      ELSEIF (LEVGEO.EQ.1.AND.LPTORR) THEN
        ANZ = 0
        IT=IPTORR
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            IF (NSTGRD(NCELL).EQ.0) ANZ=ANZ+1
          ENDDO
        ENDDO
C
        WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
        WRITE(18+ifoff,'(A8,I6,I9,I6)') 'QUAM4   ',1,ANZ,4

C  DO 100, DO 110: RETAIN SAME SEQUENCE FOR READING FROM FORT(60+IF)
C                  AS IT WAS THE CASE FOR WRITING (RPSCOL,RPSVEC)

        I=0
        IT=IPTORR
cdr     write (iunout,*) ' in rpsout levego=1 and lptorr'
        DO 100 IR=1,NR1ST
           DO 110 IP=1,NP2ND

C  FORT 60+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
             DO 105 IF=1,IRAPS
               READ (IUNRAPSVEC+IF,*) YWERT(IF)
  105        CONTINUE

             IF (IP .NE. NP2ND) THEN
               I = I + 1
               WRITE(17+ifoff,'(I6,1P,2E12.4)') I,RSURF(IR)*FCABS1(IFC),
     .                                      PSURF(IP)*FCABS2(IFC)
C  EXCLUDE DEAD CELLS ON FORT.18
               NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
               IF (IR.LT.NR1ST.AND.NSTGRD(NCELL).EQ.0) THEN
                 WRITE(18+ifoff,'(1X,A1,4I10)') '0',
     .                     (IR-1)*NP2ND+IP,
     .                     (IR-1)*NP2ND+IP+1,
     .                     IR*NP2ND+IP+1,
     .                     IR*NP2ND+IP
               ELSE
cdr              write (iunout,*) ir,ip,it,ncell,nstgrd(ncell)
               ENDIF
             ENDIF
             IF (IP .NE. NP2ND) THEN
               WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .               I,  (YWERT(IF),IF=1,IRAPS)
             ELSE
               WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .               I+1,(YWERT(IF),IF=1,IRAPS)
             ENDIF
  110      CONTINUE
           I = I + 1
           WRITE(17+ifoff,'(I6,1P,2E12.4)') I,RSURF(IR)*FCABS1(IFC),
     .                                  PSURF(NP2ND)*FCABS2(IFC)
  100   CONTINUE
        NCO=I
C
      ELSEIF (LEVGEO.LE.2.AND.LPPOLR) THEN
        ANZ = 0
        IP=IPPOLR
        DO IR=1,NR1STM
          DO IT=1,NT3RDM
            NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
            IF (NSTGRD(NCELL).EQ.0) ANZ=ANZ+1
          ENDDO
        ENDDO
C
!pb        WRITE(18+ifoff,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
        WRITE(18+ifoff,'(A8,I6,I9,I6)') 'QUAM4   ',1,ANZ,4

        I=0
        IP=IPPOLR
        DO 1100 IR=1,NR1ST
          DO 2100 IT=1,NT3RD

C  FORT IUNRAPSVEC+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
            DO 2105 IF=1,IRAPS
              READ (IUNRAPSVEC+IF,*) YWERT(IF)
 2105       CONTINUE

            IF (IT .NE. NT3RD) THEN
              I = I + 1
              WRITE(17+ifoff,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(IT)
C  EXCLUDE DEAD CELLS ON FORT.18
              NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              IF (IR.LT.NR1ST.AND.NSTGRD(NCELL).EQ.0) THEN
                WRITE(18+ifoff,'(1X,A1,4I10)') '0',
     .                    (IR-1)*NT3RD+IT,
     .                    (IR-1)*NT3RD+IT+1,
     .                     IR*NT3RD+IT+1,
     .                     IR*NT3RD+IT
              ENDIF
            ENDIF
            IF (IT .NE. NT3RD) THEN
              WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .              I,  (YWERT(IF),IF=1,IRAPS)
            ELSE
              WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .              I+1,(YWERT(IF),IF=1,IRAPS)
            ENDIF
 2100     CONTINUE
          I = I + 1
          WRITE(17+ifoff,'(I6,1P,2E12.4)') I,RSURF(IR),ZSURF(NT3RD)
 1100   CONTINUE
        NCO=I
C
      ELSEIF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
        ANZ = 0
        IT=IPTORR
        DO  IR=1,NR1ST-1
          DO  IPPLG=1,NPPLG
            DO  IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)-1
              NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
              IF((NSTGRD(NCELL).EQ.0).or.(NSTGRD(NCELL).EQ.3)) ANZ=ANZ+1
            ENDDO
          ENDDO
        ENDDO

        if (lraps3d.and.lr3dcon) then
!pb          WRITE(18+ifoff,'(A3,2X,A14,45X,I8)')
!pb     .         'PSS','1BEISPIELDATEN',(ANZ+nstab)*(iplane-1)
          WRITE(18+ifoff,'(A3,I3,A60,I6)')
     .         'PSS',1,'BEISPIELDATEN',-3334
        else
!pb          WRITE(18+ifoff,'(A3,2X,A14,45X,I8)')'PSS','1BEISPIELDATEN',
!pb     .                                   (ANZ+NSTAB)*iplane
          WRITE(18+ifoff,'(A3,I3,A60,I6)')'PSS',1,'BEISPIELDATEN',-3334
        endif

        I=0
        ipoints=0
        do ipl=0,iplane-1
          DO IF=1,IRAPS
            open (IUNRAPSVEC+IF)
          enddo

          if (lraps3d.and.lr3dcon) then
!pb            if (ipl < iplane-1)
!pb     .        WRITE(18+ifoff,'(A8,2I6,5X,A1)') 'HEXE8        ',ipl+1,ANZ,'8'
            if (ipl < iplane-1)
     .        WRITE(18+ifoff,'(A8,I6,I9,I6)') 'HEXE8   ',ipl+1,ANZ,8
            igroups = iplane-1
          else
!pb            WRITE(18+ifoff,'(A8,2I6,5X,A1)') 'QUAM4        ',ipl+1,ANZ,'4'
            WRITE(18+ifoff,'(A8,I6,I9,I6)') 'QUAM4   ',ipl+1,ANZ,4
            igroups = iplane
          endif

          IT=IPTORR
          xgeomin = huge(1._dp)
          xgeomax = -xgeomin
          ygeomin = huge(1._dp)
          ygeomax = -ygeomin
          DO 10 IR=1,NR1ST
            NPUNKT = 0
            DO 11 IPPLG=1,NPPLG
              NPUNKT = NPUNKT+NPOINT(2,IPPLG)-NPOINT(1,IPPLG)+1
   11       CONTINUE
            DO 20 IPPLG=1,NPPLG
              xgeomin = min(xgeomin,minval(
     .                  xpol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              xgeomax = max(xgeomax,maxval(
     .                  xpol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              ygeomin = min(ygeomin,minval(
     .                  ypol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              ygeomax = max(ygeomax,maxval(
     .                  ypol(ir,NPOINT(1,IPPLG):NPOINT(2,IPPLG))))
              DO 30 IP=NPOINT(1,IPPLG),NPOINT(2,IPPLG)

C  FORT IUNRAPSVEC+IF WAS WRITTEN IN RPSCOL OR RPSVEC IN SAME DO LOOPS
                DO 25 IF=1,IRAPS
                  READ (IUNRAPSVEC+IF,*) YWERT(IF)
   25           CONTINUE

                IF (IP .NE. NPOINT(2,IPPLG)) THEN
                  I = I + 1
                  if (lraps3d.and.nltra) then
                    WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .                    I,XPOL(IR,IP)*cos(ipl*rapsdel)*FCABS1(IFC),
     .                      YPOL(IR,IP)*FCABS2(IFC),
     .                      xpol(ir,ip)*sin(ipl*rapsdel)*FCABS1(IFC)
                  elseif (lraps3d.and.nltrz) then
                    WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .                    I,XPOL(IR,IP)*FCABS1(IFC),
     .                      YPOL(IR,IP)*FCABS2(IFC),ipl*rapsdel
                  else
                    WRITE(17+ifoff,'(I6,1P,2E12.4)')
     .                    I,XPOL(IR,IP)*FCABS1(IFC),
     .                      YPOL(IR,IP)*FCABS2(IFC)
                  endif
C  EXCLUDE DEAD CELLS ON FORT.18
                  NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  IF (IR.LT.NR1ST.AND.((NSTGRD(NCELL).EQ.0) .OR.
     .               (NSTGRD(NCELL).EQ.3))) THEN
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1) then
!                        WRITE(18+ifoff,'(1X,A1,8I10)') '0',
!     .                       (IR-1)*NPUNKT+IP+ipl*ipoints,
!     .                       (IR-1)*NPUNKT+IP+1+ipl*ipoints,
!     .                        IR*NPUNKT+IP+1+ipl*ipoints,
!     .                        IR*NPUNKT+IP+ipl*ipoints,
!     .                       (IR-1)*NPUNKT+IP+(ipl+1)*ipoints,
!     .                       (IR-1)*NPUNKT+IP+1+(ipl+1)*ipoints,
!     .                        IR*NPUNKT+IP+1+(ipl+1)*ipoints,
!     .                        IR*NPUNKT+IP+(ipl+1)*ipoints
                        WRITE(18+ifoff,'(1X,A1,8I10)') '0',
     .                       NOPNT(INDPOINT(IR,IP))+ipl*ipoints,
     .                       NOPNT(INDPOINT(IR,IP+1))+ipl*ipoints,
     .                       NOPNT(INDPOINT(IR+1,IP+1))+ipl*ipoints,
     .                       NOPNT(INDPOINT(IR+1,IP))+ipl*ipoints,
     .                       NOPNT(INDPOINT(IR,IP))+(ipl+1)*ipoints,
     .                       NOPNT(INDPOINT(IR,IP+1))+(ipl+1)*ipoints,
     .                       NOPNT(INDPOINT(IR+1,IP+1))+(ipl+1)*ipoints,
     .                       NOPNT(INDPOINT(IR+1,IP))+(ipl+1)*ipoints
                      endif
                    else
!                      WRITE(18+ifoff,'(1X,A1,4I10)') '0',
!     .                     (IR-1)*NPUNKT+IP+ipl*ipoints,
!     .                     (IR-1)*NPUNKT+IP+1+ipl*ipoints,
!     .                      IR*NPUNKT+IP+1+ipl*ipoints,
!     .                      IR*NPUNKT+IP+ipl*ipoints
                      WRITE(18,'(1X,A1,4I10)') '0',
     .                     NOPNT(INDPOINT(IR,IP))+ipl*ipoints,
     .                     NOPNT(INDPOINT(IR,IP+1))+ipl*ipoints,
     .                     NOPNT(INDPOINT(IR+1,IP+1))+ipl*ipoints,
     .                     NOPNT(INDPOINT(IR+1,IP))+ipl*ipoints
                   endif
                  ENDIF
                ENDIF
                IF (IP .NE. NPOINT(2,IPPLG)) THEN
                  WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .                  I,(YWERT(IF),IF=1,IRAPS)
                ELSE
                  WRITE(19+ifoff,'(I6,1P,50E12.4)')
     .                  I+1,(YWERT(IF),IF=1,IRAPS)
                ENDIF
   30         CONTINUE
              I = I + 1
              if (lraps3d.and.nltra) then
                WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .          I,XPOL(IR,NPOINT(2,IPPLG))*cos(ipl*rapsdel)*FCABS1(IFC),
     .            YPOL(IR,NPOINT(2,IPPLG))*FCABS2(IFC),
     .            xpol(ir,NPOINT(2,IPPLG))*sin(ipl*rapsdel)*FCABS1(IFC)
              elseif (lraps3d.and.nltrz) then
                WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .            I,XPOL(IR,NPOINT(2,IPPLG))*FCABS1(IFC),
     .              YPOL(IR,NPOINT(2,IPPLG))*FCABS2(IFC),ipl*rapsdel
              else
                WRITE(17+ifoff,'(I6,1P,2E12.4)')
     .            I,XPOL(IR,NPOINT(2,IPPLG))*FCABS1(IFC),
     .              YPOL(IR,NPOINT(2,IPPLG))*FCABS2(IFC)
              endif
   20       CONTINUE
   10     CONTINUE
          if (ipl.eq.0) ipoints=i
          DO IF=1,IRAPS
            close (IUNRAPSVEC+IF)
          enddo
        enddo
        NCO=I
C
      ELSEIF ((LEVGEO.EQ.4.AND.LPTORR) .OR.
     .        (LEVGEO.EQ.5.AND.LRPSCUT)) THEN
C NSTGRD.NE.0 AUSBLENDEN, ANZ NEU BERECHENEN.
        IT=IPTORR
        ANZ=0
        DO I=1,NTRII
          IF (NSTGRD(I).EQ.0) ANZ=ANZ+1
        END DO

        xgeomin = minval(xtrian(1:nrknot))
        xgeomax = maxval(xtrian(1:nrknot))
        ygeomin = minval(ytrian(1:nrknot))
        ygeomax = maxval(ytrian(1:nrknot))
        if (lraps3d.and.lr3dcon) then
!pb          WRITE(18+ifoff,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',
!pb     .                                   (ANZ+NSTAB)*(iplane-1)
          WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
        else
!pb          WRITE(18+ifoff,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',
!pb     .                                   (ANZ+NSTAB)*iplane
          WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
        endif
        do ipl=0,iplane-1
          DO IF=1,IRAPS
            open(IUNRAPSVEC+IF)
          enddo
          DO 40 I=1,NRKNOT
            DO 50 IF=1,IRAPS
              READ(IUNRAPSVEC+IF,*) YWERT(IF)
              if (ywert(if) < valcont(if)) valcont(if) = ywert(if)
   50       CONTINUE
            WRITE(19+ifoff,'(I6,1P,50E12.4)') I+ipl*nrknot,
     .                                  (YWERT(IF),IF=1,IRAPS)
            if (lraps3d.and.nltra) then
              WRITE(17+ifoff,'(I6,1P,3E12.4)') I+ipl*nrknot,
     .              XTRIAN(I)*cos(ipl*rapsdel)*FCABS1(IFC),
     .              YTRIAN(I)*FCABS2(IFC),
     .              XTRIAN(I)*sin(ipl*rapsdel)*FCABS1(IFC)
            elseif (lraps3d.and.nltrz) then
              WRITE(17+ifoff,'(I6,1P,3E12.4)') I+ipl*nrknot,
     .           XTRIAN(I)*FCABS1(IFC),YTRIAN(I)*FCABS2(IFC),ipl*rapsdel
            else
              WRITE(17+ifoff,'(I6,1P,2E12.4)') I,XTRIAN(I)*FCABS1(IFC),
     .                                           YTRIAN(I)*FCABS2(IFC)
            endif
   40     CONTINUE
          DO IF=1,IRAPS
            close(IUNRAPSVEC+IF)
          enddo
          if (lraps3d.and.lr3dcon) then
            if (ipl < iplane-1) then
!pb              WRITE(18+ifoff,'(A8,2I6,5X,A1)') 'PENTA6  ',ipl+1,ANZ,'6'
              WRITE(18+ifoff,'(A8,I6,I9,I6)') 'PENTA6  ',ipl+1,ANZ,6

              DO I=1,NTRII
                IF (NSTGRD(I).EQ.0)
     .            WRITE(18+ifoff,'(1X,A1,6I10)')
     .                                  '0',NECKE(1,I)+ipl*nrknot,
     .                                      NECKE(2,I)+ipl*nrknot,
     .                                      NECKE(3,I)+ipl*nrknot,
     .                                      NECKE(1,I)+(ipl+1)*nrknot,
     .                                      NECKE(2,I)+(ipl+1)*nrknot,
     .                                      NECKE(3,I)+(ipl+1)*nrknot
              enddo
            endif
            igroups = iplane-1
          else
!pb            WRITE(18+ifoff,'(A8,2I6,5X,A1)') 'TRIM3   ',ipl+1,ANZ,'3'
            WRITE(18+ifoff,'(A8,I6,I9,I6)') 'TRIM3   ',ipl+1,ANZ,3

            DO 60 I=1,NTRII
              IF (NSTGRD(I).EQ.0)
     .          WRITE(18+ifoff,'(1X,A1,3I10)')
     .                                '0',NECKE(1,I)+ipl*nrknot,
     .                                    NECKE(2,I)+ipl*nrknot,
     .                                    NECKE(3,I)+ipl*nrknot
   60       CONTINUE
            igroups = iplane
          endif
        enddo
        NCO=NRKNOT*iplane
C
      ELSEIF (LEVGEO.EQ.5.AND..NOT.LRPSCUT) THEN
C TO BE DONE: NSTGRD.NE.0 AUSBLENDEN, ANZ NEU BERECHNEN.
        anz=ntet-ntet_collaps
        do i=1,ncoor
          do if=1,iraps
            read(IUNRAPSVEC+if,*) ywert(if)
          enddo
          WRITE(19+ifoff,'(I6,1P,50E12.4)') I,(YWERT(IF),IF=1,IRAPS)
          WRITE(17+ifoff,'(I6,1P,3E12.4)')
     .          I,XTETRA(I),YTETRA(I),ZTETRA(I)
        enddo
!pb        WRITE(18+ifoff,'(A3,2X,A14,45X,I8)') 'PSS','1BEISPIELDATEN',ANZ+NSTAB
        WRITE(18+ifoff,'(A3,I3,A60,I6)') 'PSS',1,'BEISPIELDATEN',-3334
        WRITE(18+ifoff,'(A8,I6,I9,I6)') 'TET4    ',1,ANZ,4
        do i=1,ntet
          if (ntbar(1,i) >= 0)
     .      WRITE(18+ifoff,'(1X,A1,4I10)') '0',NTECK(1,I),NTECK(2,I),
     .                                  NTECK(3,I),NTECK(4,I)
        enddo
        nco = anz
      ELSE
        nco = 0
      ENDIF
C
      if (nstab > 0) then
c ist das jemals getestet worden ???
        IGR=igroups+1
        do ipl=0,iplane-1
          if (lraps3d.and.lr3dcon) then
             if (ipl < iplane-1) then
                 WRITE(18+ifoff,'(A8,I6,I9,I6)')
     .                'QUAM4   ',IGR+ipl,NSTAB,4
                 igroups = igroups+1
              endif
          else
            WRITE(18+ifoff,'(A8,I6,I9,I6)') 'FLA2    ',IGR+ipl,NSTAB,2
            igroups = igroups+1
          endif
C
          CUR => FIRST_POINT
          IF (ASSOCIATED(CUR)) THEN
            DO WHILE (ASSOCIATED(CUR%NXTPNT))
! check if last point of line has been reached
              IF (CUR%NXTPNT%NPL2D.NE.0) THEN
                if (lraps3d.and.nltra) then
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                   CUR%XPL2D*cos(ipl*rapsdel)*FCABS1(IFC),
     .                   CUR%YPL2D*FCABS2(IFC),
     .                   CUR%XPL2D*sin(ipl*rapsdel)*FCABS1(IFC)
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D*cos(ipl*rapsdel)*FCABS1(IFC),
     .                   CUR%NXTPNT%YPL2D*FCABS2(IFC),
     .                   CUR%NXTPNT%XPL2D*sin(ipl*rapsdel)*FCABS1(IFC)
                elseif (lraps3d.and.nltrz) then
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                   CUR%XPL2D*FCABS1(IFC),
     .                   CUR%YPL2D*FCABS2(IFC),
     .                   ipl*rapsdel
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                   CUR%NXTPNT%XPL2D*FCABS1(IFC),
     .                   CUR%NXTPNT%YPL2D*FCABS2(IFC),
     .                   ipl*rapsdel
                else
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,2E12.4)') NCO,
     .                   CUR%XPL2D*FCABS1(IFC),
     .                   CUR%YPL2D*FCABS2(IFC)
                  NCO=NCO+1
                  WRITE (17+ifoff,'(I6,1P,2E12.4)')
     .                   NCO,CUR%NXTPNT%XPL2D*FCABS1(IFC),
     .                       CUR%NXTPNT%YPL2D*FCABS2(IFC)
                endif
                if (lraps3d.and.lr3dcon) then
                  if (ipl < iplane-1)
     .              WRITE (18+ifoff,'(1X,A1,4I10)') '0',
     .                     NCO-1,NCO,nco+ipl*2*nstab,
     .                     nco-1+ipl*2*nstab
                else
                  WRITE (18+ifoff,'(1X,A1,2I10)') '0',NCO-1,NCO
                endif
              ENDIF
              CUR => CUR%NXTPNT
            END DO
          END IF
C
          IF ((LEVGEO.EQ.2.OR.LEVGEO.EQ.3).AND.LPTORR) THEN
            NST=NSTSI
            IF (NTIME.GT.0) NST=NST-1
            DO ISTS=1,NST
              IF (ILIIN(NLIM+ISTS) > 0) THEN
                IDIMP=TRANSFER(MAXLOC(INUMP(ISTS,1:3)),1)

                IF (IDIMP == 1) THEN
                  IR=INUMP(ISTS,1)
                  IA=IRPTA(ISTS,2)
                  IE=IRPTE(ISTS,2)
                  DO IP=IA,IE-1
                    if (lraps3d.and.nltra) then
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*cos(ipl*rapsdel)*FCABS1(IFC),
     .                       YPOL(IR,IP)*FCABS2(IFC),
     .                       XPOL(IR,IP)*sin(ipl*rapsdel)*FCABS1(IFC)
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP+1)*cos(ipl*rapsdel)*FCABS1(IFC),
     .                       YPOL(IR,IP+1)*FCABS2(IFC),
     .                       XPOL(IR,IP+1)*sin(ipl*rapsdel)*FCABS1(IFC)
                    elseif (lraps3d.and.nltrz) then
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*FCABS1(IFC),
     .                       YPOL(IR,IP)*FCABS2(IFC),
     .                       ipl*rapsdel
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP+1)*FCABS1(IFC),
     .                       YPOL(IR,IP+1)*FCABS2(IFC),
     .                       ipl*rapsdel
                    else
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP)*FCABS1(IFC),
     .                           YPOL(IR,IP)*FCABS2(IFC)
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP+1)*FCABS1(IFC),
     .                           YPOL(IR,IP+1)*FCABS2(IFC)
                    endif
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1)
     .                  WRITE (18+ifoff,'(1X,A1,4I10)') '0',
     .                         NCO-1,NCO,nco+ipl*2*nstab,
     .                         nco-1+ipl*2*nstab
                    else
                      WRITE (18+ifoff,'(1X,A1,2I10)') '0',NCO-1,NCO
                    endif
                  END DO

                ELSEIF (IDIMP ==2) THEN
                  IP=INUMP(ISTS,IDIMP)
                  IA=IRPTA(ISTS,1)
                  IE=IRPTE(ISTS,1)
                  DO IR=IA,IE-1
                    if (lraps3d.and.nltra) then
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*cos(ipl*rapsdel)*FCABS1(IFC),
     .                       YPOL(IR,IP)*FCABS2(IFC),
     .                       XPOL(IR,IP)*sin(ipl*rapsdel)*FCABS1(IFC)
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR+1,IP)*cos(ipl*rapsdel)*FCABS1(IFC),
     .                       YPOL(IR+1,IP)*FCABS2(IFC),
     .                       XPOL(IR+1,IP)*sin(ipl*rapsdel)*FCABS1(IFC)
                    elseif (lraps3d.and.nltrz) then
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR,IP)*FCABS1(IFC),
     .                       YPOL(IR,IP)*FCABS2(IFC),
     .                       ipl*rapsdel
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,3E12.4)') NCO,
     .                       XPOL(IR+1,IP)*FCABS1(IFC),
     .                       YPOL(IR+1,IP)*FCABS2(IFC),
     .                       ipl*rapsdel
                    else
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR,IP)*FCABS1(IFC),
     .                           YPOL(IR,IP)*FCABS2(IFC)
                      NCO=NCO+1
                      WRITE (17+ifoff,'(I6,1P,2E12.4)')
     .                       NCO,XPOL(IR+1,IP)*FCABS1(IFC),
     .                           YPOL(IR+1,IP)*FCABS2(IFC)
                    endif
                    if (lraps3d.and.lr3dcon) then
                      if (ipl < iplane-1)
     .                  WRITE (18+ifoff,'(1X,A1,4I10)') '0',
     .                         NCO-1,NCO,nco+ipl*2*nstab,
     .                         nco-1+ipl*2*nstab
                    else
                      WRITE (18+ifoff,'(1X,A1,2I10)') '0',NCO-1,NCO
                    endif
                  END DO
                END IF
              END IF
            END DO
          END IF
        enddo
      endif
C
      if (lraps3d  .and. ((levgeo == 3) .or. (levgeo == 4)) ) then
        if (ncontour > 0) then
          ALLOCATE (phelp(maxval(nconpoint(1:ncontour))-1,2))
          ALLOCATE (xcont(maxval(nconpoint(1:ncontour))))
          ALLOCATE (ycont(maxval(nconpoint(1:ncontour))))
          igr=igroups+1
          epsrel = (xgeomax-xgeomin+ygeomax-ygeomin)/2.*eps5
          do icont=1,ncontour
            write(IUNRAPS+ifoff,*) 'xcontour ycontour icont: ',icont
            do ipoint = 1,nconpoint(icont)
               write(IUNRAPS+ifoff,'(2es12.4)') xcontour(ipoint,icont),
     .              ycontour(ipoint,icont)
            enddo
            write(IUNRAPS+ifoff,*)
c           bereinigte contour erstellen
            xcont(1:nconpoint(icont))=xcontour(1:nconpoint(icont),icont)
            ycont(1:nconpoint(icont))=ycontour(1:nconpoint(icont),icont)
            ncont = nconpoint(icont)
            del_point = .true.
            ip_start = 1
            do while (del_point)
               do ipoint=ip_start,ncont-1
                  vecax = xcont(ipoint)-xcont(ipoint+1)
                  vecay = ycont(ipoint)-ycont(ipoint+1)
                  length_veca = sqrt(vecax**2+vecay**2)
                  if (length_veca <
     .                (100./10**ceiling(log10(abs(xcont(ipoint)))) +
     .                 100./10**ceiling(log10(abs(xcont(ipoint)))))/2.)
     .            then
                     if (ipoint == 1) then
                        x1 = xcont(ncont-1)
                        y1 = ycont(ncont-1)
                     else
                        x1 = xcont(ipoint-1)
                        y1 = ycont(ipoint-1)
                     endif
                     x2 = xcont(ipoint)
                     y2 = ycont(ipoint)
                     x3 = xcont(ipoint+1)
                     y3 = ycont(ipoint+1)
                     if (ipoint == (ncont-1)) then
                        x4 = xcont(ipoint+2)
                        y4 = ycont(ipoint+2)
                     else
                        x4 = xcont(2)
                        y4 = ycont(2)
                     endif
                     atri1 = x1*(y3-y2)+x3*(y2-y1)+x2*(y1-y3)
                     atri2 = x4*(y3-y2)+x3*(y2-y4)+x2*(y4-y3)
                     if ((atri1 < -eps5) .and. (atri2 < -eps5)) then
                        ncont = ncont - 1
                        xcont(ipoint) = (xcont(ipoint)+
     .                                   xcont(ipoint+1))/2.
                        ycont(ipoint) = (ycont(ipoint)+
     .                                   ycont(ipoint+1))/2.
                        do idel=ipoint+1, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=max(1,ipoint-1)
                        exit
                     elseif (((atri1 > 0) .and. (atri2 < 0)) .or.
     .                       (abs(atri1) < eps5)) then
                        ncont = ncont - 1
                        do idel=ipoint, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=max(1,ipoint-1)
                        exit
                     elseif (((atri1 < 0) .and. (atri2 > 0)) .or.
     .                       (abs(atri2) < eps5)) then
                        ncont = ncont - 1
                        do idel=ipoint+1, ncont
                           xcont(idel) = xcont(idel+1)
                           ycont(idel) = ycont(idel+1)
                        enddo
                        xcont(ncont+1) = 0.
                        ycont(ncont+1) = 0.
                        ip_start=ipoint
                        exit
                     endif
                  endif
                  if (ipoint == (ncont-1)) then
                     del_point = .false.
                  endif
               enddo
             enddo
             do ipoint=1,NCONT-1
c     punkt berechnen
               if (ipoint == 1) then
                  vecax=xcont(ncont-1)-xcont(ipoint)
                  vecay=ycont(ncont-1)-ycont(ipoint)
               else
                  vecax=xcont(ipoint-1)-xcont(ipoint)
                  vecay=ycont(ipoint-1)-ycont(ipoint)
               endif
               vecbx = xcont(ipoint+1)-xcont(ipoint)
               vecby = ycont(ipoint+1)-ycont(ipoint)
               length_veca = sqrt(vecax**2+vecay**2)
               length_vecb = sqrt(vecbx**2+vecby**2)
               vecax = vecax/length_veca
               vecay = vecay/length_veca
               vecbx = vecbx/length_vecb
               vecby = vecby/length_vecb
               del_eps = min(epsrel,length_veca/10.,length_vecb/10.)
               xm = xcont(ipoint) + vecax + vecbx
               ym = ycont(ipoint) + vecay + vecby
               length_m_p = sqrt((xm-xcont(ipoint))**2+
     .              (ym-ycont(ipoint))**2)
               ahelp = xm*(ycont(ipoint+1)-ycont(ipoint)) +
     .              xcont(ipoint+1) * (ycont(ipoint)-ym) +
     .              xcont(ipoint) * (ym - ycont(ipoint+1))
c     punkt in richtung m verschieben
               if ((length_m_p < eps10) .or. (abs(ahelp) < eps10)) then
                  xm = -1.*ycont(ipoint)+ycont(ipoint+1)
                  ym = xcont(ipoint)-xcont(ipoint+1)
                  length_m_p = sqrt(xm**2+ym**2)
                  phelp(ipoint,1) = xcont(ipoint)+xm/length_m_p*del_eps
                  phelp(ipoint,2) = ycont(ipoint)+ym/length_m_p*del_eps
               iloop = 1
               do while (
     .              (abs(xcont(ipoint)-phelp(ipoint,1))<
     .              100./10**ceiling(log10(abs(xcont(ipoint)))))
     .              .and.
     .              (abs(ycont(ipoint)-phelp(ipoint,2))<
     .              100./10**ceiling(log10(abs(ycont(ipoint))))))
                  iloop = iloop + 1
                  phelp(ipoint,1) = xcont(ipoint)+xm/
     .                 length_m_p*del_eps*iloop
                  phelp(ipoint,2) = ycont(ipoint)+ym/
     .                 length_m_p*del_eps*iloop
               enddo
               else
               phelp(ipoint,1) = xcont(ipoint)+
     .              ahelp/abs(ahelp+eps30) *
     .              (xm-xcont(ipoint))/length_m_p*del_eps
               phelp(ipoint,2) = ycont(ipoint)+
     .              ahelp/abs(ahelp+eps30) *
     .              (ym-ycont(ipoint))/length_m_p*del_eps
               iloop = 1
               do while (
     .              (abs(xcont(ipoint)-phelp(ipoint,1))<
     .              100./10**ceiling(log10(abs(xcont(ipoint)))))
     .              .and.
     .              (abs(ycont(ipoint)-phelp(ipoint,2))<
     .              100./10**ceiling(log10(abs(ycont(ipoint))))))
                  iloop = iloop + 1
                  phelp(ipoint,1) = xcont(ipoint)+
     .                 ahelp/abs(ahelp+eps30) *
     .                 (xm-xcont(ipoint))/length_m_p*
     .                 del_eps*iloop
                  phelp(ipoint,2) = ycont(ipoint)+
     .                 ahelp/abs(ahelp+eps30) *
     .                 (ym-ycont(ipoint))/length_m_p*
     .                 del_eps*iloop
               enddo
             endif
            enddo
            do ipl=0,iplane-1
               write(IUNRAPS+ifoff,*) 'origx origy neux neuy'
               do ipoint=1,ncont-1
                  write(IUNRAPS+ifoff,'(4es12.4)')
     .                 xcont(ipoint), ycont(ipoint),
     .                 phelp(ipoint,1),phelp(ipoint,2)
c     punkte schreiben
                  if (nltra) then
                     WRITE(17+ifoff,'(I6,1P,3E12.4)') nco+ipoint,
     .                    phelp(Ipoint,1)*cos(ipl*rapsdel),
     .                    phelp(Ipoint,2),
     .                    phelp(Ipoint,1)*sin(ipl*rapsdel)
                     WRITE(19+ifoff,'(I6,1P,50E12.4)') nco+ipoint,
     .                                          (valcont(IF),IF=1,IRAPS)
                  elseif (nltrz) then
                     WRITE(17+ifoff,'(I6,1P,3E12.4)') nco+ipoint,
     .                    phelp(Ipoint,1),
     .                    phelp(Ipoint,2),
     .                    ipl*rapsdel
                     WRITE(19+ifoff,'(I6,1P,50E12.4)') nco+ipoint,
     .                                          (valcont(IF),IF=1,IRAPS)
                  endif
               enddo
               if (ipl < iplane-1) then
c     elementgruppe anfangen
                  WRITE(18+ifoff,'(A8,I6,I9,I6)')
     .                 'QUAM4   ',IGR+ipl,ncont-1,4
                  do ipoint=1,ncont-1
c     element schreiben
                     if (ipoint < ncont-1) then
                        WRITE(18+ifoff,'(1X,A1,4I10)') '0',
     .                       nco+ipoint,nco+ipoint+1,
     .                       nco+ipoint+1+ncont-1,
     .                       nco+ipoint+ncont-1
                     else
                        WRITE(18+ifoff,'(1X,A1,4I10)') '0',
     .                       nco+ipoint,nco+1,
     .                       nco+1+ncont-1,
     .                       nco+ipoint+ncont-1
                     endif
                  enddo
               endif
               nco = nco + ncont-1
            enddo
            igr=igr+iplane-1
          enddo
          deallocate(phelp)
          deallocate(xcont)
          deallocate(ycont)
        endif
      endif

      WRITE(17+ifoff,'(1X,A5,8X,A3,12X,A1,11X,A1,11X,A1,11X,A1)')
     .           '-9999','FIN','0','0','0','0'
      WRITE(19+ifoff,'(1X,A5,8X,A3,50(11X,I1))') '-9999',
     .           'FIN',(0,IF=1,IRAPS)
      deallocate(valcont)

      RETURN
      END SUBROUTINE EIRENE_RPSOUT
