C
c  written by P. Boerner, for FZJ proprietary IDL plotting tool.
c  not intended for 3rd party use.
c  last modified: jan 2017


cdr correction: 3 digits rather than 2 digits for I0 in outtal file name
cdr jan 22: add std dev. tallies SIGMA
cdr mar 22: try to document....
cdr aug 22: separate output files for individual strata
cdr         tbd: make pltstr(istra) flags also here.

      SUBROUTINE EIRENE_OUTIDLTAL
cdr write all volumetric output tallies, each one on a separate file.

cdr abuse of language: tallies are 2d arrays A(N1,N2).
cdr                    N1 "species index", (not necessarily a species, can also be a vector component, etc..)
cdr                    N2  cell number index.


cdr for each tally: write data for all "species" (1st dimension of tally array), ordered
cdr                 such that each "species" has its own column.
cdr for each tally skip cells in which all "species" entries for that tally are zero.
cdr consequently: the tally files written here may have different lengths,
cdr               depending on how many non-zero rows a tally has.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI

      IMPLICIT NONE
C
      REAL(DP), ALLOCATABLE :: VECTOR(:,:), TALAV(:), TALTOT(:)
      REAL(DP) :: OUTAUI, SMEAN
      INTEGER :: NFTI, NFTE, K, ITAL, I, ISTR, MXSPZ, IOUT, ICELL, KK
      INTEGER :: ISP, ISIG
      LOGICAL :: LFIRST
      CHARACTER(6) :: CISTR, CITAL, CISP, CNFTI, CNFTE
C
      CHARACTER(50) :: FNAME, FORMA, FORME, FORME2
      EXTERNAL :: EIRENE_RSTRT, EIRENE_SYMET, EIRENE_LEER
C
!      IF (NSBOX_TAL /= NSBOX) THEN
!         WRITE (IUNOUT,*) ' ERROR IN OUTIDLTAL '
!         WRITE (IUNOUT,*) ' NSBOX_TAL /= NSBOX '
!         WRITE (IUNOUT,*) ' THIS CASE IS NOT YET FORESEEN '
!         WRITE (IUNOUT,*) ' NO DATA WRITTEN '
!         RETURN
!      END IF

      MXSPZ = MAXVAL(NFSTVI(1:NTALV))
      MXSPZ = MAX(MXSPZ, NATM, NMOL, NION, NPHOT, NADV, NALV)
      ALLOCATE (VECTOR(NRAD,MXSPZ))
      ALLOCATE (TALTOT(MXSPZ))
      ALLOCATE (TALAV(MXSPZ))

      IOUT = 3 + IFOFF

      FORMA=REPEAT(' ',50)
      FORMA='(6X,   A25)'
      WRITE (FORMA(5:7),'(I3)') MXSPZ

      FORME=REPEAT(' ',50)
      FORME='(I10,   ES25.7)'
      WRITE (FORME(6:8),'(I3)') MXSPZ

      FORME2=REPEAT(' ',50)
      FORME2='(6X,   ES25.7)'
      WRITE (FORME2(5:7),'(I3)') MXSPZ

      LFIRST=.TRUE.
C
C ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA

      DO ISTR = 0, NSTRAI
        ISTRA=ISTR

        IF ((ISTRA == 0) .AND. (NSMSTRA /= 1)) CYCLE
C
        IF (XMCP(ISTRA).LT.1.) CYCLE
C
C
        IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
        ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
cdr  read ISTRA tallies from fort.10.
cdr  Totals (outau, fort.11) are still on storage.
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,
     .             NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .             TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
cdr  read only ISTRA=0 (sum over strata) tallies from fort.10.
cdr  Totals (outau, fort.11) are still on storage.
        ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
          IESTR=ISTRA
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,
     .             NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .             TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSE
          WRITE (iunout,*) 'ERROR IN OUTIDLTAL: STRATUM ISTRA= ',
     .                     ISTRA
          WRITE (iunout,*) 'DATA ARE NOT AVAILABLE. PRINTOUT ABANDONED'
          CALL EIRENE_LEER(1)
          CYCLE
        ENDIF
C
C
C  PRINT VOLUME-AVERAGED TALLIES, STRATUM: ISTRA
        write (cistr,'(I0)') istra

        DO 100 ITAL = 1, NTALV

          IF (.NOT.LIVTALV(ITAL)) THEN
            WRITE (iunout,*) ' TALLY NOT AVAILABLE (OUTIDLTAL)',
     .                       ' ITAL = ', ITAL
            CALL EIRENE_LEER(1)
            CYCLE
          ENDIF
C
          NFTI=1
          NFTE=NFSTVI(ITAL)

          IF (ITAL == NTALA) THEN
            NFTE = SUM(VERIFY(TXTSPC(1:NADV,NTALA),' '))
          END IF

          KK = 0
          DO 119 K=NFTI,NFTE
            IF (NEXTVI(ITAL) > 0) THEN
              KK = KK + 1
              IF ((K > NEXTVI(ITAL) .AND.
     .             (MOD(K,NEXTVI(ITAL)) == 1))) KK = KK + 1
            ELSE
              KK = K
            END IF
            CALL EIRENE_FETCH_OUTAUI (OUTAUI,ITAL,KK,ISTRA,IUNOUT)
C
            IF (NSBOX_TAL /= NSBOX) THEN
              DO I=1,NSBOX
                ICELL = NCLTAL(I)
                IF (ICELL > 0) THEN
                  VECTOR(I,K)=ESTIMV(NADDV(ITAL)+K,ICELL)
                ELSE
                  VECTOR(I,K)=0._DP
                END IF
              END DO
            ELSE
              DO 110 I=1,NSBOX_TAL
                VECTOR(I,K)=ESTIMV(NADDV(ITAL)+K,I)
  110         CONTINUE
          END IF

            TALTOT(K)=OUTAUI
            TALAV(K)=TALTOT(K)/VOLTOT
  119     CONTINUE
C
          write (cital,'(I0)') ital
          write (cnfti,'(I0)') nfti
          write (cnfte,'(I0)') nfte
          FNAME =
     .     'outtal_'//trim(cistr)//'_'//trim(cital)//
     .           '_'//trim(cnfti)//'-'//trim(cnfte)

          IF (LFIRST) THEN
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL')
          ELSE
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL',POSITION='APPEND')
          END IF

          WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'
          WRITE (IOUT,'(A,I6)') 'ISTRA = ',ISTRA
          WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'

cdr  tally name is the same for all "species"
          WRITE (IOUT,'(A)') TXTTAL(1,ITAL)

          WRITE (IOUT,'(A,I10)') 'NCELLS:   ',NSBOX
          WRITE (IOUT,'(A,I10)') 'NSPECIES: ',NFTE
          WRITE (IOUT,'(A)') 'SPECIES'
          WRITE (IOUT,FORMA) (TRIM(TXTSPC(K,ITAL)), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'UNITS'
          WRITE (IOUT,FORMA) (TRIM(TXTUNT(K,ITAL)), K=NFTI, NFTE)

          WRITE (IOUT,'(A)')
     .        'TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") '
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) (TALTOT(K), K=NFTI, NFTE)
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) (TALAV(K), K=NFTI, NFTE)

          WRITE (IOUT,'(A)')
     .      '================================================='
          DO I=1, NSBOX
cdr  only write cells with at least one non-zero contribution
            IF (ANY(ABS(VECTOR(I,NFTI:NFTE)) > EPS30))
     .        WRITE (IOUT,FORME) I,(VECTOR(I,K), K=NFTI, NFTE)
          END DO
          WRITE (IOUT,'(A)')
     .      '================================================='

          CLOSE (UNIT=IOUT)
C
  100   CONTINUE  ! ITAL

C  CHECK IF STANDARD DEVIATION TALLIES ARE AVAILABLE
        DO 101 ISIG=1,NSIGVI
          ITAL=IIH(ISIG)
          ISP =IGH(ISIG)
          IF (ISP.EQ.0) THEN
          ENDIF
          DO 112 I=1,NSBOX_TAL
            VECTOR(I,1)=SIGMA(ISIG,I)
  112     CONTINUE
          SMEAN=SGMS(ISIG)
C
          write (cital,'(I0)') ital
          write (cisp,'(I0)') isp
          FNAME = 'outtal_'//trim(cital)//'_'//trim(cisp)//'_STDV'

          IF (LFIRST) THEN
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL')
          ELSE
            OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .            ACCESS='SEQUENTIAL',POSITION='APPEND')
          END IF

          WRITE (IOUT,'(A)')
     .    '+++++++++++++++++++++++++++++++++++++++++++++++++'
          WRITE (IOUT,'(A,I6)') 'ISTRA = ',ISTRA
          WRITE (IOUT,'(A)')
     .    '+++++++++++++++++++++++++++++++++++++++++++++++++'

cdr  feb.22: special treatment for vectorial output tallies mapl_vec, mmpl_vec,....
cdr          "species index" may be just a particular vector component instead.
          if (ital.ge.97 .and. ital.le.100) then
            WRITE (IOUT,'(A)') TXTTAL(isp,ITAL)
          else
cdr  tally name is the same for all "species"
            WRITE (IOUT,'(A)') TXTTAL(1,ITAL)
          endif
          WRITE (IOUT,'(A,I10)') 'NCELLS:   ',NSBOX
          WRITE (IOUT,'(A,I10)') 'NSPECIES: ',1
          WRITE (IOUT,'(A)') 'SPECIES'
          WRITE (IOUT,FORMA) TRIM(TXTSPC(ISP,ITAL))
          WRITE (IOUT,'(A)') 'UNITS'
          WRITE (IOUT,FORMA) TRIM(TXTUNT(ISP,ITAL))

          WRITE (IOUT,'(A)') 'SDV OF TOTAL (%) '
          WRITE (IOUT,'(A)') 'TOTAL'
          WRITE (IOUT,FORME2) SMEAN
          WRITE (IOUT,'(A)') 'MEAN'
          WRITE (IOUT,FORME2) SMEAN

          WRITE (IOUT,'(A)')
     .      '================================================='

cdr  print only rows with at least one non-zero entry. Skip all the others.
cdr  This is the same as is done in OUTIDLPLA for input (field particle) tallies.
          DO I=1, NSBOX
            IF (ANY(ABS(VECTOR(I,1:1)) > EPS30))
     .        WRITE (IOUT,FORME) I,VECTOR(I,1)
          END DO
          WRITE (IOUT,'(A)')
     .      '================================================='

          CLOSE (UNIT=IOUT)

  101   CONTINUE  ! ISIG

        LFIRST = .FALSE.

      END DO  ! ISTR

      DEALLOCATE (VECTOR)
      DEALLOCATE (TALTOT)
      DEALLOCATE (TALAV)

      RETURN
      END SUBROUTINE EIRENE_OUTIDLTAL
