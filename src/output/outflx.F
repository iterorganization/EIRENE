C  printout surface fluxes (incl. sputter fluxes) for stratum 'istra'
C  loop over all surfaces selected in input block 11 for printout

cdr Feb. 2015:  total sputter tallies now included, resolved wrt. incident type
cdr these total tallies may include sputtering of unidentified wall material,
cdr  hence may be different from the totals obtained by sum over sputtered species-resolved fluxes
cdr  e.g. sptatot may be larger than summt, etc....

cdr SEPT.2014:  PRINTOUT OF SPUTTERED FLUXES REVISED
cdr  total sputter fluxes spttot(msurf) added.
CDR  TO BE DONE:
cdr  total sputter flux tallies resolved wrt. incident type to be added:
cdr  sptatot(msurf),sptmtot(msurf),sptitot(msurf),sptpltot(msurf),stpphtot(msurf)
cdr  many local variables have a redundant istra index, to enable MASYR1 printout format.

C
C
C 27.03.04; iliin=-3 option (only net fluxes on transp. surfaces) re-enforced
C           simultaneous changes in escape.f.
C           not active for bulk particle fluxes updated in subr. locate
C 05.05.04:  printout of spectrum: text improved
C
c 12.12.05:  some more lines shifted 2 columns to the right
c            (still not complete, done until line 467),
c            This is for surface do loop "do 10000"
c 16.12.05:  variance for adds tally: printout activated
c
c 16.01.06:  bug fix: suma1, suma2, etc... initialized (=0)
c            otherwise problems due to new options for deactivation of tallies
C 07.12.06:  some comments introduced to clarify status with iliin=-3 option
C 18.04.16:  reduced string length to match variable, J.Lore
cdr nov.21:  typo corrected, re. reflected photons
C
      SUBROUTINE EIRENE_OUTFLX(A,ISTRA)

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_CSDVI
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_CTRIG
      USE EIRMOD_CLOGAU, ONLY: NLSPCSPT_ATM, NLSPCSPT_MOL, NLSPCSPT_ION,
     .                         NLSPCSPT_PHOT, NLSPCSPT_PLS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      CHARACTER(32), INTENT(IN) :: A
      REAL(DP) :: SUMA1(0:NATM,0:NSTRA), VARA1(0:NATM,0:NSTRA),
     .          SUMM1(0:NMOL,0:NSTRA), VARM1(0:NMOL,0:NSTRA),
     .          SUMI1(0:NION,0:NSTRA), VARI1(0:NION,0:NSTRA),
     .          SUMPH1(0:NPHOT,0:NSTRA), VARPH1(0:NPHOT,0:NSTRA),
     .          SUMP1(0:NPLS,0:NSTRA), VARP1(0:NPLS,0:NSTRA),
     .          SUMA2(0:NATM,0:NSTRA), VARA2(0:NATM,0:NSTRA),
     .          SUMM2(0:NMOL,0:NSTRA), VARM2(0:NMOL,0:NSTRA),
     .          SUMI2(0:NION,0:NSTRA), VARI2(0:NION,0:NSTRA),
     .          SUMPH2(0:NPHOT,0:NSTRA), VARPH2(0:NPHOT,0:NSTRA),
     .          SUMP2(0:NPLS,0:NSTRA), VARP2(0:NPLS,0:NSTRA),
     .          SUMS(0:NADS,0:NSTRA),  VARS(0:NADS,0:NSTRA),
     .          SUML(0:NALS,0:NSTRA),  VARL(0:NALS,0:NSTRA)
      REAL(DP), ALLOCATABLE, SAVE :: HELP(:)
      REAL(DP), ALLOCATABLE :: HELP2(:,:), TOTAR(:)
      REAL(DP) :: HELPP(NLMPGS)
      REAL(DP) :: SUM1, DUMMY, SUMMT, SUMMEM, SUMMTI, SUMMEI,
     .          SUMMTA, SUMMTM, SUMMA, SUMMS, SUMMM, SUMMI, SUMMP,
     .          SUMA, SUMM, SUMI, SUMME, SUMMTP, SUMMEP, SUMP, TTTT,
     .          SUMMEA, SUMMEPH, SUMMTPH, SUMMPH, SUMPH, EN
      REAL(DP) :: TTSPTA, TTSPTM, TTSPTI, TTSPTPH, TTSPTP, TTSPT
      INTEGER :: NR, NP, NT, MSURFG, J, NCELL, NTOTAL, N1, N2, N3, I,
     .           ITRII, IPLGN, ISPR, ITEXT, ISTS, NFTI, NFTE,
     .           K, ITALS, IADS, IALS, IION, IPLS, IATM, N,
     .           IMOL, IPHOT, ISPC, IOUT, ISF, IE, IT, ILIM, ID, N1D,
     .           IRA, IRE, IPA, IPE, IS
      INTEGER :: IADTYP(0:4)
      LOGICAL :: LGVRA1(0:NATM,0:NSTRA),
     .           LGVRM1(0:NMOL,0:NSTRA),
     .           LGVRI1(0:NION,0:NSTRA),
     .           LGVRPH1(0:NPHOT,0:NSTRA),
     .           LGVRP1(0:NPLS,0:NSTRA)
      LOGICAL :: LGVRA2(0:NATM,0:NSTRA),
     .           LGVRM2(0:NMOL,0:NSTRA),
     .           LGVRI2(0:NION,0:NSTRA),
     .           LGVRPH2(0:NPHOT,0:NSTRA),
     .           LGVRP2(0:NPLS,0:NSTRA)
      LOGICAL :: LGVARS(0:NADS,0:NSTRA),
     .           LGVARL(0:NALS,0:NSTRA)
      LOGICAL :: LOGADS(0:NADS,0:NSTRA),
     .           LOGALS(0:NALS,0:NSTRA)
      LOGICAL :: PRINTED(NLIMPS)
      CHARACTER(8) :: TEXTA(NADS), TEXTL(NALS)
      CHARACTER(10) :: TEXTYP(0:4)
      CHARACTER(8) :: UNITINT(1:3), UNITOUT
      EXTERNAL :: EIRENE_INTVOL, EIRENE_PRTTLS, EIRENE_LEER,
     .            EIRENE_MASAGE, EIRENE_MASBOX, EIRENE_MASR1,
     .            EIRENE_MASYR1

      IF (.NOT.ALLOCATED(HELP)) ALLOCATE (HELP(NRAD))

      LGVRA1=.FALSE.
      LGVRM1=.FALSE.
      LGVRI1=.FALSE.
      LGVRPH1=.FALSE.
      LGVRP1=.FALSE.
      LGVRA2=.FALSE.
      LGVRM2=.FALSE.
      LGVRI2=.FALSE.
      LGVRPH2=.FALSE.
      LGVRP2=.FALSE.
      LGVARS=.FALSE.
      LGVARL=.FALSE.
      LOGADS=.FALSE.
      LOGALS=.FALSE.

      PRINTED = .FALSE.
      CALL EIRENE_LEER(1)
      WRITE (iunout,9999) A
C
C  SURFACE LOOP
C
      DO 10000 ISPR=1,NSURPR
C
C
        ITEXT=0
        I=NPRSRF(ISPR)
cdr remove disabled additional surfaces from printout
cdr All tallies are empty at these surfaces anyway
        if (i.le.nlim .and. iliin(i) == 0) cycle
        CALL EIRENE_LEER(2)
        IF (IGJUM0(I).NE.0) THEN
C  CHECK, IF SURFACE "I" IS STILL THERE AS ONE SIDE OF A TRIANGLE
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
            DO IPLGN=1,3
              ISTS=ABS(INMTI(IPLGN,ITRII))
              IF (ISTS.EQ.I) GOTO 5
            ENDDO
            ENDDO
          ENDIF
          IF (I.GT.NLIM) THEN
            ILIM = -(I-NLIM)
          ELSE
            ILIM = I
          ENDIF
          WRITE (iunout,*) ' SURFACE NO. ',ILIM,' : OUT '
          PRINTED(I) = .TRUE.
          GOTO 10000
        ENDIF
    5   CONTINUE
C  PRINT SURFACE AREA (NOT FOR "TIME SURFACE")
        CALL EIRENE_MASBOX(TXTSFL(I))
        IF (NTIME.GE.1.AND.I.EQ.NLIM+NSTSI) GOTO 1
        IF (SAREA(I).NE.666.) THEN
          WRITE (iunout,'(A22,1P,1E12.4)') ' SURFACE AREA (CM**2) ',
     .           SAREA(I)
        ELSE
          WRITE (iunout,*) 'SURFACE AREA (CM**2) ','?'
        ENDIF
    1   CONTINUE
C
C
        IF (I.GT.NLIM.AND.LEVGEO.LE.4.AND.NLMPGS.NE.NLIMPS) THEN
C
C  SPATIAL RESOLUTION ON NON-DEFAULT STANDARD SURFACE?
          HELP=0._DP
          ISTS=I-NLIM
C
          ITALS=NPRTLS(ISPR)
          IF (ITALS.LE.0.OR.ITALS.GT.NTALS) GOTO 11
          IF (.not.LIVTALS(ITALS)) GOTO 11
          NFTI=1
          NFTE=NFSTWI(ITALS)
          IF (NSPEZS(ISPR,1).GT.0) THEN
            NFTI=NSPEZS(ISPR,1)
            NFTE=MAX(NFTI,NSPEZS(ISPR,2))
          ENDIF
          IF ((LEVGEO == 4) .AND. (NFLAGS(ISPR) > 10)) THEN
            ALLOCATE(HELP2(SURF_TRIAN_ORDERED(I)%NUMTR+1,NFTI:NFTE))
            ALLOCATE(TOTAR(NFTI:NFTE))
            HELP2 = 0._DP
            TOTAR = 0._DP
          END IF
          DO 10 K=NFTI,NFTE
            IF (K.GT.NFSTWI(ITALS)) THEN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE IN OUTFLX'
              WRITE (iunout,*) 'ISTS,ITALS, K, ',ISTS,ITALS,K
              CALL EIRENE_LEER(1)
              GOTO 10
            ENDIF
C
            DO J=1,NLMPGS
              HELPP(J)=ESTIMS(NADDW(ITALS)+K,J)
            END DO
C
            sum1 = 0._dp
            n1 = 0
            n2 = 0
            n3 = 0
            select case (LEVGEO)
            case (:3)
              IF (INUMP(ISTS,2).NE.0) then
C  POLOIDAL SURFACE
                sum1=0
                np=1
                do nr=1,nr1st
                  do nt=1,nt3rd
                    MSURFG=NR+(NT-1)*NR1P2
                    MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                    NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                    HELP(ncell)=HELPP(MSURFG)
                    sum1=sum1+HELPP(msurfg)
                  ENDDO
                ENDDO
                N1=NR1ST
                N2=1
                N3=NT3RD
              ELSEIF (INUMP(ISTS,1).NE.0) then
C  RADIAL SURFACE
                sum1=0
                NR=1
                do np=1,np2nd
                  do nt=1,nt3rd
                    MSURFG=NP+(NT-1)*NP2T3
                    MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                    NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                    HELP(ncell)=HELPP(MSURFG)
                    sum1=sum1+HELPP(msurfg)
                  ENDDO
                ENDDO
                N1=1
                N2=NP2ND
                N3=NT3RD
              ELSEIF (INUMP(ISTS,3).NE.0) then
C  TOROIDAL SURFACE
                sum1=0
                nt=1
                do nr=1,nr1st
                  do np=1,np2nd
                    MSURFG=NR+(NP-1)*NR1P2
                    MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                    NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                    HELP(ncell)=HELPP(MSURFG)
                    sum1=sum1+HELPP(msurfg)
                  ENDDO
                ENDDO
                N1=NR1ST
                N2=NP2ND
                N3=1
              ENDIF
              IRA = IRPTA(ISTS,1)
              IRE = IRPTE(ISTS,1)
              IPA = IRPTA(ISTS,2)
              IPE = IRPTE(ISTS,2)
C  TRIANGULAR GRID
            case (4)
              sum1 = 0._dp
              nt = surf_trian_ordered(i)%numtr
              do j = 1, nt
                 it = surf_trian_ordered(i)%itrias(j)
                 is = surf_trian_ordered(i)%itrisi(j)
                 msurfg=nlim+nsts+inspat(is,it)
                 help(j) = helpp(msurfg)
                 sum1 = sum1 + helpp(msurfg)
              end do
              if (nflags(ispr) > 10) help2(1:nt,k) = help(1:nt)
              N1=NT+1
              N2=1
              N3=1
              IRA = 1
              IRE = N1
              IPA = 1
              IPE = N2
            end select

            NTOTAL=N1*N2*N3
            IF (NTOTAL > 0) THEN
              CALL EIRENE_INTVOL (HELP,1,1,NTOTAL,DUMMY,N1,N2,N3,1)
              IF (NFLAGS(ISPR) > 10) THEN
                 HELP2(NTOTAL,K) = DUMMY
                 TOTAR(K) = DUMMY
              END IF
              IF (ABS(DUMMY) > EPS60) THEN
                CALL EIRENE_PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,MOD(NFLAGS(ISPR),10),
     .                  NTLSFL(ISPR),
     .                  IRA,IRE,IPA,IPE,1,1)
              ELSE
                CALL EIRENE_PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,-1,
     .                  NTLSFL(ISPR),
     .                  IRA,IRE,IPA,IPE,1,1)
                CALL EIRENE_MASAGE
     .            ('IDENTICALLY ZERO, NOT PRINTED                ')
                CALL EIRENE_LEER(2)
              END IF
            END IF
   10     CONTINUE

          IF ((LEVGEO == 4) .AND. (NFLAGS(ISPR) > 10)) THEN
             CALL EIRENE_PRINT_SURF_TRIAN (IUNOUT)
             IF (NTLSFL(ISPR) > 0) THEN
                OPEN(UNIT=NTLSFL(ISPR),POSITION='APPEND')
                CALL EIRENE_PRINT_SURF_TRIAN (NTLSFL(ISPR))
                CLOSE(UNIT=NTLSFL(ISPR))
             END IF
             DEALLOCATE(HELP2)
             DEALLOCATE(TOTAR)
          END IF
   11     CONTINUE
        ENDIF  ! non-def. standard surf. with spatial resolution: done


C  SPECTRA
        IF (NTLSFL(ISPR) > 0) THEN
          ISF=I
          IF (I < 0) ISF=ABS(I)+NLIM
          IOUT = NTLSFL(ISPR)+ifoff

          TEXTYP(0) = 'PHOTONS   '
          TEXTYP(1) = 'ATOMS     '
          TEXTYP(2) = 'MOLECULES '
          TEXTYP(3) = 'TEST IONS '
          TEXTYP(4) = 'BULK IONS '

cdr units of energy integrated totals
          UNITINT(1)= '(AMP)   '
          UNITINT(2)= '(WATT)  '
          UNITINT(3)= '(DYN)   '
          UNITOUT   = '(?)     '

          IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

          DO ISPC=1,NADSPC
            IF ((ESTIML(ISPC)%ISRFCLL == 0) .AND.
     .          (ESTIML(ISPC)%ISPCSRF == ISF)) THEN
              WRITE (IOUT,*)
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' SPECTRUM CALCULATED FOR SURFACE ',I
              IT = ESTIML(ISPC)%ISPCTYP
              UNITOUT=UNITINT(IT)
              IF (IT == 1) THEN
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'PARTICLE FLUX IN AMP'
              ELSE
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'ENERGY FLUX IN WATT '
              END IF
c................................................................

              WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .              TEXTYP(ESTIML(ISPC)%IPRTYP)
              IF (ESTIML(ISPC)%IPRSP == 0) THEN
                WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
              ELSE
                WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .                 TEXTS(IADTYP(ESTIML(ISPC)%IPRTYP)+
     .                       ESTIML(ISPC)%IPRSP)
              END IF
              IF (ESTIML(ISPC)%LOG) THEN
                WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .               10._DP**ESTIML(ISPC)%SPCMIN
                WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .               10._DP**ESTIML(ISPC)%SPCMAX
                WRITE (IOUT,'(A)') ' LOGARITHMIC SPACING'
              ELSE
                WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .               ESTIML(ISPC)%SPCMIN
                WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .               ESTIML(ISPC)%SPCMAX
                WRITE (IOUT,'(A)') ' LINEAR SPACING'
              END IF
              WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ',
     .               ESTIML(ISPC)%NSPC
              WRITE (IOUT,*)
              IF (ESTIML(ISPC)%SPCS > EPS60) THEN
                IF (NSIGI_SPC == 0) THEN
C  STANDARD DEVIATION IS NOT AVAILABLE
                  DO IE=1, ESTIML(ISPC)%NSPC
C  DEAL WITH ENERGY BIN NO. IE
C  central energy bin value
                    EN = ESTIML(ISPC)%SPCMIN +
     .                   (IE-0.5_DP)*ESTIML(ISPC)%SPCDEL
                    IF (ESTIML(ISPC)%LOG) THEN
                      WRITE (IOUT,'(I6,2ES12.4)') IE,10._DP**EN,
     .                 ESTIML(ISPC)%SPC(IE)
                    ELSE
                      WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .                 ESTIML(ISPC)%SPC(IE)
                    END IF
                  END DO
                ELSE
c
C  STANDARD DEVIATION IS AVAILABLE
C
                  DO IE=1, ESTIML(ISPC)%NSPC
                    EN = ESTIML(ISPC)%SPCMIN +
     .                   (IE-0.5_DP)*ESTIML(ISPC)%SPCDEL
                    IF (ESTIML(ISPC)%LOG) THEN
                      WRITE (IOUT,'(I6,3ES12.4)') IE,10._DP**EN,
     .                   ESTIML(ISPC)%SPC(IE),
     .                   ESTIML(ISPC)%SDV(IE)
                    ELSE
                      WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
     .                   ESTIML(ISPC)%SPC(IE),
     .                   ESTIML(ISPC)%SDV(IE)
                    END IF
                  END DO
                END IF
chk Legendre polynomial expansion tallies
                IF (ESTIML(ISPC)%ISPCOPT==2) THEN
                  DO IE=1, ESTIML(ISPC)%NSPC
                    EN = ESTIML(ISPC)%SPCMIN +
     .                   (IE-0.5_DP)*ESTIML(ISPC)%SPCDEL
                    IF (ESTIML(ISPC)%LOG) THEN
                      WRITE (IOUT,'(I6,1ES12.4)') IE,10._DP**EN
                    ELSE
                      WRITE (IOUT,'(I6,1ES12.4)') IE,EN
                    END IF
                    WRITE (IOUT,'(99(1ES12.4))')
     .           (ESTIML(ISPC)%SPCAN(ID,IE),ID=1,ESTIML(ISPC)%ISPLDEG)
                  END DO
                END IF
              ELSE
                WRITE (IOUT,*) ' SPECTRUM IDENTICALLY 0 '
              END IF
C
C  PRINTOUT OF ENERGY INTEGRAL OVER SPECTRA
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' INTEGRAL OF SPECTRUM ',
     .               ESTIML(ISPC)%SPCS
              IF (NSIGI_SPC > 0)
     .          WRITE (IOUT,*) ' STANDARD DEVIATION  ',
     .               ESTIML(ISPC)%SGMS
            END IF
          END DO
        END IF

        IF (PRINTED(I)) CYCLE
C
C  *****************************************************
C   INCIDENT FLUXES, POSITIVE PARTIAL CURRENTS, ALSO: NET FLUXES
C  *****************************************************
C
C
C   SURFACE-AVERAGED TALLY NO.1 AND NO.26
C
        SUMMT=0.
        SUMME=0.
        SUMMTP=0.
        SUMMEP=0.
        SUMA=0.
        SUMM=0.
        SUMI=0.
        SUMP=0.
        SUMPH=0.
        SUMA1(:,ISTRA)=0._DP
        SUMA2(:,ISTRA)=0._DP
        DO 20 IATM=1,NATMI
          LGVRA1(IATM,ISTRA)=.FALSE.
          LGVRA2(IATM,ISTRA)=.FALSE.
          IF (LPOTAT) SUMA1(IATM,ISTRA)=POTAT(IATM,I)
          SUMMT=SUMMT+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
          SUMA=SUMA+SUMA1(IATM,ISTRA)
          IF (LEOTAT) SUMA2(IATM,ISTRA)=EOTAT(IATM,I)
          SUMME=SUMME+SUMA2(IATM,ISTRA)
   20   CONTINUE
C
        DO 21 N=1,NSIGSI
          IF (IIHW(N).EQ.1) THEN
            DO 22 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 22
              VARA1(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
   22       CONTINUE
          ELSEIF (IIHW(N).EQ.26) THEN
            DO 522 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 522
              VARA2(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
  522       CONTINUE
          ENDIF
   21   CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.7 AND NO.32
C
        SUMM1(:,ISTRA)=0._DP
        SUMM2(:,ISTRA)=0._DP
        DO 30 IMOL=1,NMOLI
          LGVRM1(IMOL,ISTRA)=.FALSE.
          LGVRM2(IMOL,ISTRA)=.FALSE.
          IF (LPOTML) SUMM1(IMOL,ISTRA)=POTML(IMOL,I)
          SUMMT=SUMMT+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
          SUMM=SUMM+SUMM1(IMOL,ISTRA)
          IF (LEOTML) SUMM2(IMOL,ISTRA)=EOTML(IMOL,I)
          SUMME=SUMME+SUMM2(IMOL,ISTRA)
   30   CONTINUE
C
        DO 23 N=1,NSIGSI
          IF (IIHW(N).EQ.7) THEN
            DO 24 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 24
              VARM1(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
   24       CONTINUE
          ELSEIF (IIHW(N).EQ.32) THEN
            DO 524 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 524
              VARM2(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
  524       CONTINUE
          ENDIF
   23   CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.13 AND NO.38
C
        SUMI1(:,ISTRA)=0._DP
        SUMI2(:,ISTRA)=0._DP
        DO 40 IION=1,NIONI
          LGVRI1(IION,ISTRA)=.FALSE.
          LGVRI2(IION,ISTRA)=.FALSE.
          IF (LPOTIO) SUMI1(IION,ISTRA)=POTIO(IION,I)
          SUMMT=SUMMT+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
          SUMI=SUMI+SUMI1(IION,ISTRA)
          IF (LEOTIO) SUMI2(IION,ISTRA)=EOTIO(IION,I)
          SUMME=SUMME+SUMI2(IION,ISTRA)
   40   CONTINUE
C
        DO 25 N=1,NSIGSI
          IF (IIHW(N).EQ.13) THEN
            DO 26 IION=1,NIONI
              IF (IGHW(N).NE.IION) GOTO 26
              VARI1(IION,ISTRA)=SIGMAW(N,I)
              LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
   26       CONTINUE
          ELSEIF (IIHW(N).EQ.38) THEN
            DO 526 IION=1,NIONI
              IF (IGHW(N).NE.IION) GOTO 526
              VARI2(IION,ISTRA)=SIGMAW(N,I)
              LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
  526       CONTINUE
          ENDIF
   25   CONTINUE
C
C   SURFACE-AVERAGED TALLY NO. 19 and NO. 44
C
        SUMPH1(:,ISTRA)=0._DP
        SUMPH2(:,ISTRA)=0._DP
        DO IPHOT=1,NPHOTI
          LGVRPH1(IPHOT,ISTRA)=.FALSE.
          LGVRPH2(IPHOT,ISTRA)=.FALSE.
          IF (LPOTPHT) SUMPH1(IPHOT,ISTRA)=POTPHT(IPHOT,I)
          SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
          SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
          IF (LEOTPHT) SUMPH2(IPHOT,ISTRA)=EOTPHT(IPHOT,I)
          SUMME=SUMME+SUMPH2(IPHOT,ISTRA)
        ENDDO

        DO N=1,NSIGSI
          IF (IIHW(N).EQ.19) THEN
            DO IPHOT=1,NPHOTI
              IF (IGHW(N).NE.IPHOT) CYCLE
              VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
              LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
            ENDDO
          ELSEIF (IIHW(N).EQ.44) THEN
            DO IPHOT=1,NPHOTI
              IF (IGHW(N).NE.IPHOT) CYCLE
              VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
              LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
            ENDDO
          ENDIF
        ENDDO
C
C   SURFACE-AVERAGED TALLY NO.25 AND NO.50
C
        SUMP1(:,ISTRA)=0._DP
        SUMP2(:,ISTRA)=0._DP
        DO 50 IPLS=1,NPLSI
         LGVRP1(IPLS,ISTRA)=.FALSE.
         LGVRP2(IPLS,ISTRA)=.FALSE.
         IF (LPOTPL) SUMP1(IPLS,ISTRA)=POTPL(IPLS,I)
         SUMMTP=SUMMTP+SUMP1(IPLS,ISTRA)*NPRT(NSPAMI+IPLS)
         SUMP=SUMP+SUMP1(IPLS,ISTRA)
         IF (LEOTPL) SUMP2(IPLS,ISTRA)=EOTPL(IPLS,I)
         SUMMEP=SUMMEP+SUMP2(IPLS,ISTRA)
   50   CONTINUE
C
        DO 27 N=1,NSIGSI
          IF (IIHW(N).EQ.25) THEN
            DO 28 IPLS=1,NPLSI
              IF (IGHW(N).NE.IPLS) GOTO 28
              VARP1(IPLS,ISTRA)=SIGMAW(N,I)
              LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
   28       CONTINUE
          ELSEIF (IIHW(N).EQ.50) THEN
            DO 528 IPLS=1,NPLSI
              IF (IGHW(N).NE.IPLS) GOTO 528
              VARP2(IPLS,ISTRA)=SIGMAW(N,I)
              LGVRP2(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
  528       CONTINUE
          ENDIF
   27   CONTINUE
C
C  printout starts here
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMP)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP) THEN
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASAGE
     .  ('NO FLUXES INCIDENT ON THIS SURFACE             ')
        CALL EIRENE_LEER(1)
      ELSE
        CALL EIRENE_LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (iunout,*) 'FLUX INCIDENT ON SURFACE:'
C  SURFACE-AVERAGED TALLY NO. 1
          IF (SUMA.NE.0._DP) THEN
            WRITE (iunout,*) 'INCIDENT: ATOMS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 26
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 7
          IF (SUMM.NE.0._DP) THEN
            WRITE (iunout,*) 'INCIDENT: MOLECULES'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 32
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 13
          IF (SUMI.NE.0._DP) THEN
            WRITE (iunout,*) 'INCIDENT: TEST IONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 38
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0._DP) THEN
            WRITE (iunout,*) 'INCIDENT: PHOTONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 44
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 25
          IF (SUMP.NE.0._DP) THEN
            WRITE (iunout,*) 'INCIDENT: BULK IONS (RECYCLING SOURCE)'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
C  SURFACE-AVERAGED TALLY NO. 50
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMP2,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP2,LGVRP2,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
          ENDIF
          CALL EIRENE_LEER (1)
          WRITE (iunout,*)
     .      'TOTAL INCIDENT "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0._DP)
     .    WRITE (iunout,*) 'EXCLUDING BULK IONS (RECYCLING SOURCE)'
          CALL EIRENE_MASR1 ('TOT.PFLX',SUMMT)
          CALL EIRENE_MASR1 ('TOT.EFLX',SUMME)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          IF (SUMA.NE.0._DP.OR.SUMM.NE.0._DP.OR.SUMI.NE.0._DP
     .                    .OR.SUMPH.NE.0._DP)
     .    WRITE (iunout,*)
     .      'PARTIAL PARTICLE AND ENERGY CURRENTS, POSITIVE'
C  SURFACE-AVERAGED TALLY NO. 1
          IF (SUMA.NE.0._DP) THEN
            WRITE (iunout,*) 'ATOMS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 26
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 7
          IF (SUMM.NE.0._DP) THEN
            WRITE (iunout,*) 'MOLECULES'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 32
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 13
          IF (SUMI.NE.0._DP) THEN
            WRITE (iunout,*) 'TEST IONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 38
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0._DP) THEN
            WRITE (iunout,*) 'PHOTONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 44
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 25
          IF (SUMP.NE.0._DP) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*)
     .        'PRIMARY FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (iunout,*) 'BULK IONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
C  SURFACE-AVERAGED TALLY NO. 50
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMP2,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP2,LGVRP2,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
          ENDIF
          CALL EIRENE_LEER (1)
          WRITE (iunout,*)
     .      'TOTAL POSITIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0._DP)
     .    WRITE (iunout,*) 'EXCLUDING BULK IONS (RECYCLING SOURCE)'
          CALL EIRENE_MASR1 ('POS.PFLX',SUMMT)
          CALL EIRENE_MASR1 ('POS.EFLX',SUMME)
C
C.....................................................................
C
        ELSEIF (ILIIN(I).EQ.-3) THEN
          IF (SUMA.NE.0._DP.OR.SUMM.NE.0._DP.OR.SUMI.NE.0._DP
     .                    .OR.SUMPH.NE.0._DP)
     .    WRITE (iunout,*) 'NET PARTICLE AND ENERGY CURRENTS'
C  SURFACE-AVERAGED TALLY NO. 1
          IF (SUMA.NE.0._DP) THEN
            WRITE (iunout,*) 'ATOMS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 26
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 7
          IF (SUMM.NE.0._DP) THEN
            WRITE (iunout,*) 'MOLECULES'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 32
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 13
          IF (SUMI.NE.0._DP) THEN
            WRITE (iunout,*) 'TEST IONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 38
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
          ENDIF
C  SURFACE-AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0._DP) THEN
            WRITE (iunout,*) 'PHOTONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 44
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
          ENDIF

C  SURFACE-AVERAGED TALLY NO. 25
          IF (SUMP.NE.0._DP) THEN
            CALL EIRENE_LEER(1)
            WRITE (iunout,*)
     .        'PRIMARY FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (iunout,*) 'BULK IONS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
C  SURFACE-AVERAGED TALLY NO. 50
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMP2,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP2,LGVRP2,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
          ENDIF
          CALL EIRENE_LEER (1)
          WRITE (iunout,*) 'TOTAL NET "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0._DP)
     .    WRITE (iunout,*) 'EXCLUDING BULK IONS (RECYCLING SOURCE)'
          CALL EIRENE_MASR1 ('NET PFLX',SUMMT)
          CALL EIRENE_MASR1 ('NET EFLX',SUMME)

        ENDIF
C
C  INDEPENDENT OF VALUE AND SIGN OF ILIIN:
C
        IF (SUMMTP.GT.0._DP) THEN
          CALL EIRENE_LEER (1)
          WRITE (iunout,*)
     .      'TOTAL INCIDENT RECYCLING SOURCE "ATOMIC" FLUXES'
          WRITE (iunout,*) 'BULK IONS'
          CALL EIRENE_MASR1 ('SRC.PFLX',SUMMTP)
          CALL EIRENE_MASR1 ('SRC.EFLX',SUMMEP)
        ENDIF
      ENDIF
C
C  ******************************************
C   REEMITTED FLUXES, NEGATIVE PARTIAL FLUXES
C  ******************************************
C
C   FIRST: FROM INCIDENT ATOMS
C
      SUMMTA=0._DP
      SUMMEA=0._DP
      SUMA=0._DP
      SUMM=0._DP
      SUMI=0._DP
      SUMPH=0._DP
C
C   SURFACE-AVERAGED TALLY NO.2 AND NO.27
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFAAT) SUMA1(IATM,ISTRA)=PRFAAT(IATM,I)
        SUMMTA=SUMMTA+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFAAT) SUMA2(IATM,ISTRA)=ERFAAT(IATM,I)
        SUMMEA=SUMMEA+SUMA2(IATM,ISTRA)
  102 CONTINUE
C
      DO 121 N=1,NSIGSI
        IF (IIHW(N).EQ.2) THEN
          DO 122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
  122     CONTINUE
        ELSEIF (IIHW(N).EQ.27) THEN
          DO 622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
  622     CONTINUE
        ENDIF
  121 CONTINUE
C
C
C   SURFACE-AVERAGED TALLY NO.8 AND NO.33
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFAML) SUMM1(IMOL,ISTRA)=PRFAML(IMOL,I)
        SUMMTA=SUMMTA+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFAML) SUMM2(IMOL,ISTRA)=ERFAML(IMOL,I)
        SUMMEA=SUMMEA+SUMM2(IMOL,ISTRA)
  103 CONTINUE
C
      DO 123 N=1,NSIGSI
        IF (IIHW(N).EQ.8) THEN
          DO 124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
  124     CONTINUE
        ELSEIF (IIHW(N).EQ.33) THEN
          DO 624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
  624     CONTINUE
        ENDIF
  123 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.14 AND NO.39
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFAIO) SUMI1(IION,ISTRA)=PRFAIO(IION,I)
        SUMMTA=SUMMTA+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFAIO) SUMI2(IION,ISTRA)=ERFAIO(IION,I)
        SUMMEA=SUMMEA+SUMI2(IION,ISTRA)
  104 CONTINUE
C
      DO 125 N=1,NSIGSI
        IF (IIHW(N).EQ.14) THEN
          DO 126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
  126     CONTINUE
        ELSEIF (IIHW(N).EQ.39) THEN
          DO 626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
  626     CONTINUE
        ENDIF
  125 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.20 AND NO.45
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFAPHT) SUMPH1(IPHOT,ISTRA)=PRFAPHT(IPHOT,I)
        SUMMTA=SUMMTA+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFAPHT) SUMPH2(IPHOT,ISTRA)=ERFAPHT(IPHOT,I)
        SUMMEA=SUMMEA+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.20) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.45) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP.AND.ILIIN(I).GT.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT ATOMS'
        CALL EIRENE_LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL EIRENE_LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT ATOMS:'
C  SURFACE-AVERAGED TALLY NO. 2
        IF (SUMA.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 27
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 8
        IF (SUMM.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 33
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 14
        IF (SUMI.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 39
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 20
        IF (SUMPH.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 45
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        ENDIF
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTA)
        CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEA)
C
        ELSEIF (ILIIN(I).LT.0) THEN  ! STILL ILIIN.NE.-3 AT THIS POINT
C  SURFACE-AVERAGED TALLY NO. 2
          IF (SUMA.NE.0._DP) THEN
            WRITE (iunout,*)
     .        'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (iunout,*) 'ATOMS'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 27
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
            CALL EIRENE_LEER (1)
            CALL EIRENE_MASR1 ('NEG.PFLX',SUMMTA)
            CALL EIRENE_MASR1 ('NEG.EFLX',SUMMEA)
          ENDIF
        ENDIF
C     ELSEIF(ILIIN = -3): NOTHING TO BE DONE HERE, NET FLUXES ARE ALREADY PRINTED.
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT MOLECULES
      SUMMTM=0.
      SUMMEM=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE-AVERAGED TALLY NO.3 AND NO.28
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 1102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFMAT) SUMA1(IATM,ISTRA)=PRFMAT(IATM,I)
        SUMMTM=SUMMTM+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFMAT) SUMA2(IATM,ISTRA)=ERFMAT(IATM,I)
        SUMMEM=SUMMEM+SUMA2(IATM,ISTRA)
 1102 CONTINUE
C
      DO 1121 N=1,NSIGSI
        IF (IIHW(N).EQ.3) THEN
          DO 1122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 1122     CONTINUE
        ELSEIF (IIHW(N).EQ.28) THEN
          DO 1622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 1622     CONTINUE
        ENDIF
 1121 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.9 AND NO.34
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 1103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFMML) SUMM1(IMOL,ISTRA)=PRFMML(IMOL,I)
        SUMMTM=SUMMTM+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFMML) SUMM2(IMOL,ISTRA)=ERFMML(IMOL,I)
        SUMMEM=SUMMEM+SUMM2(IMOL,ISTRA)
 1103 CONTINUE
C
      DO 1123 N=1,NSIGSI
        IF (IIHW(N).EQ.9) THEN
          DO 1124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 1124     CONTINUE
        ELSEIF (IIHW(N).EQ.34) THEN
          DO 1624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 1624     CONTINUE
        ENDIF
 1123 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.15 AND NO.40
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 1104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFMIO) SUMI1(IION,ISTRA)=PRFMIO(IION,I)
        SUMMTM=SUMMTM+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFMIO) SUMI2(IION,ISTRA)=ERFMIO(IION,I)
        SUMMEM=SUMMEM+SUMI2(IION,ISTRA)
 1104 CONTINUE
C
      DO 1125 N=1,NSIGSI
        IF (IIHW(N).EQ.15) THEN
          DO 1126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
 1126     CONTINUE
        ELSEIF (IIHW(N).EQ.40) THEN
          DO 1626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
 1626     CONTINUE
        ENDIF
 1125 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.21 AND NO.46
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFMPHT) SUMPH1(IPHOT,ISTRA)=PRFMPHT(IPHOT,I)
        SUMMTM=SUMMTM+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFMPHT) SUMPH2(IPHOT,ISTRA)=ERFMPHT(IPHOT,I)
        SUMMEM=SUMMEM+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.21) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.46) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP.AND.ILIIN(I).GT.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT MOLECULES'
        CALL EIRENE_LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL EIRENE_LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT MOLECULES:'
C  SURFACE-AVERAGED TALLY NO. 3
        IF (SUMA.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 28
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 9
        IF (SUMM.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 34
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 15
        IF (SUMI.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 40
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 21
        IF (SUMPH.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 46
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        ENDIF
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTM)
        CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEM)
C
        ELSEIF (ILIIN(I).LT.0) THEN  ! STILL ILIIN.NE.-3 AT THIS POINT
C  SURFACE-AVERAGED TALLY NO. 9
          IF (SUMM.NE.0._DP) THEN
            IF (ITEXT.EQ.0) WRITE (iunout,*)
     .         'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (iunout,*) 'MOLECULES'
            CALL EIRENE_MASYR1('P-FLUX:  ',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 34
            CALL EIRENE_MASYR1('E-FLUX:  ',
     .       SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
            CALL EIRENE_LEER (1)
            CALL EIRENE_MASR1 ('NEG.PFLX',SUMMTM)
            CALL EIRENE_MASR1 ('NEG.EFLX',SUMMEM)
          ENDIF
        ENDIF
C     ELSEIF(ILIIN = -3): NOTHING TO BE DONE HERE, NET FLUXES ARE ALREADY PRINTED.
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT TEST IONS
      SUMMTI=0.
      SUMMEI=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE-AVERAGED TALLY NO.4 AND NO.29
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 2102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFIAT) SUMA1(IATM,ISTRA)=PRFIAT(IATM,I)
        SUMMTI=SUMMTI+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFIAT) SUMA2(IATM,ISTRA)=ERFIAT(IATM,I)
        SUMMEI=SUMMEI+SUMA2(IATM,ISTRA)
 2102 CONTINUE
C
      DO 2121 N=1,NSIGSI
        IF (IIHW(N).EQ.4) THEN
          DO 2122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 2122     CONTINUE
        ELSEIF (IIHW(N).EQ.29) THEN
          DO 2622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 2622     CONTINUE
        ENDIF
 2121 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.10 AND NO.35
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 2103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFIML) SUMM1(IMOL,ISTRA)=PRFIML(IMOL,I)
        SUMMTI=SUMMTI+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFIML) SUMM2(IMOL,ISTRA)=ERFIML(IMOL,I)
        SUMMEI=SUMMEI+SUMM2(IMOL,ISTRA)
 2103 CONTINUE
C
      DO 2123 N=1,NSIGSI
        IF (IIHW(N).EQ.10) THEN
          DO 2124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 2124     CONTINUE
        ELSEIF (IIHW(N).EQ.35) THEN
          DO 2624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 2624     CONTINUE
        ENDIF
 2123 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.16 AND NO.41
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 2104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFIIO) SUMI1(IION,ISTRA)=PRFIIO(IION,I)
        SUMMTI=SUMMTI+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFIIO) SUMI2(IION,ISTRA)=ERFIIO(IION,I)
        SUMMEI=SUMMEI+SUMI2(IION,ISTRA)
 2104 CONTINUE
C
      DO 2125 N=1,NSIGSI
        IF (IIHW(N).EQ.16) THEN
          DO 2126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
 2126     CONTINUE
        ELSEIF (IIHW(N).EQ.41) THEN
          DO 2626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
 2626     CONTINUE
        ENDIF
 2125 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO 22 AND NO.47
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFIPHT) SUMPH1(IPHOT,ISTRA)=PRFIPHT(IPHOT,I)
        SUMMTI=SUMMTI+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFIPHT) SUMPH2(IPHOT,ISTRA)=ERFIPHT(IPHOT,I)
        SUMMEI=SUMMEI+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.22) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.47) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP.AND.ILIIN(I).GT.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT TEST IONS'
        CALL EIRENE_LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL EIRENE_LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT TEST IONS:'
C  SURFACE-AVERAGED TALLY NO. 4
        IF (SUMA.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 29
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 10
        IF (SUMM.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 35
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 16
        IF (SUMI.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 41
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 22
        IF (SUMPH.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 47
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        ENDIF
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTI)
        CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEI)
C
        ELSEIF (ILIIN(I).LT.0) THEN  ! STILL ILIIN.NE.-3 AT THIS POINT
C  SURFACE-AVERAGED TALLY NO. 16
        IF (SUMI.NE.0._DP) THEN
        IF (ITEXT.EQ.0) WRITE (iunout,*)
     .    'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
        ITEXT=1
        WRITE (iunout,*) 'TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 41
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('NEG.PFLX',SUMMTI)
        CALL EIRENE_MASR1 ('NEG.EFLX',SUMMEI)
        ENDIF
        ENDIF
C     ELSEIF(ILIIN = -3): NOTHING TO BE DONE HERE, NET FLUXES ARE ALREADY PRINTED.
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT PHOTONS
      SUMMTPH=0.
      SUMMEPH=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE-AVERAGED TALLY NO.5 AND NO.30
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 4102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFPHAT) SUMA1(IATM,ISTRA)=PRFPHAT(IATM,I)
        SUMMTPH=SUMMTPH+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFPHAT) SUMA2(IATM,ISTRA)=ERFPHAT(IATM,I)
        SUMMEPH=SUMMEPH+SUMA2(IATM,ISTRA)
 4102 CONTINUE
C
      DO 4121 N=1,NSIGSI
        IF (IIHW(N).EQ.5) THEN
          DO 4122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 4122     CONTINUE
        ELSEIF (IIHW(N).EQ.30) THEN
          DO 4622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 4622     CONTINUE
        ENDIF
 4121 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.11 AND NO.36
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 4103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFPHML) SUMM1(IMOL,ISTRA)=PRFPHML(IMOL,I)
        SUMMTPH=SUMMTPH+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFPHML) SUMM2(IMOL,ISTRA)=ERFPHML(IMOL,I)
        SUMMEPH=SUMMEPH+SUMM2(IMOL,ISTRA)
 4103 CONTINUE
C
      DO 4123 N=1,NSIGSI
        IF (IIHW(N).EQ.11) THEN
          DO 4124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 4124     CONTINUE
        ELSEIF (IIHW(N).EQ.36) THEN
          DO 4624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 4624     CONTINUE
        ENDIF
 4123 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.17 AND NO.42
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 4104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFPHIO) SUMI1(IION,ISTRA)=PRFPHIO(IION,I)
        SUMMTPH=SUMMTPH+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFPHIO) SUMI2(IION,ISTRA)=ERFPHIO(IION,I)
        SUMMEPH=SUMMEPH+SUMI2(IION,ISTRA)
 4104  CONTINUE
C
      DO 4125 N=1,NSIGSI
        IF (IIHW(N).EQ.17) THEN
          DO 4126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
 4126     CONTINUE
        ELSEIF (IIHW(N).EQ.42) THEN
          DO 4626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
 4626     CONTINUE
        ENDIF
 4125 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO 23 AND NO.48
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFPHPHT) SUMPH1(IPHOT,ISTRA)=PRFPHPHT(IPHOT,I)
        SUMMTPH=SUMMTPH+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFPHPHT) SUMPH2(IPHOT,ISTRA)=ERFPHPHT(IPHOT,I)
        SUMMEPH=SUMMEPH+SUMPH2(IPHOT,ISTRA)
      ENDDO

      DO N=1,NSIGSI
        IF (IIHW(N).EQ.23) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.48) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP.AND.ILIIN(I).GT.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT PHOTONS'
        CALL EIRENE_LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL EIRENE_LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT PHOTONS:'
C  SURFACE-AVERAGED TALLY NO. 5
        IF (SUMA.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 30
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 11
        IF (SUMM.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 36
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 17
        IF (SUMI.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 42
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 48
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        ENDIF
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTPH)
        CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEPH)
C
        ELSEIF (ILIIN(I).LT.0) THEN  ! STILL ILIIN.NE.-3 AT THIS POINT
C  SURFACE-AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0._DP) THEN
        IF (ITEXT.EQ.0) WRITE (iunout,*)
     .    'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
        ITEXT=1
        WRITE (iunout,*) 'PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 48
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_LEER (1)
        CALL EIRENE_MASR1 ('NEG.PFLX',SUMMTPH)
        CALL EIRENE_MASR1 ('NEG.EFLX',SUMMEPH)
        ENDIF
        ENDIF
C     ELSEIF(ILIIN = -3): NOTHING TO BE DONE HERE, NET FLUXES ARE ALREADY PRINTED.
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT BULK IONS
      SUMMTP=0.
      SUMMEP=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE-AVERAGED TALLY NO.6 AND NO.31
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 3102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFPAT) SUMA1(IATM,ISTRA)=PRFPAT(IATM,I)
        SUMMTP=SUMMTP+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFPAT) SUMA2(IATM,ISTRA)=ERFPAT(IATM,I)
        SUMMEP=SUMMEP+SUMA2(IATM,ISTRA)
 3102 CONTINUE
C
      DO 3121 N=1,NSIGSI
        IF (IIHW(N).EQ.6) THEN
          DO 3122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 3122     CONTINUE
        ELSEIF (IIHW(N).EQ.31) THEN
          DO 3622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
 3622     CONTINUE
        ENDIF
 3121 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.12 AND NO.37
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 3103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFPML) SUMM1(IMOL,ISTRA)=PRFPML(IMOL,I)
        SUMMTP=SUMMTP+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFPML) SUMM2(IMOL,ISTRA)=ERFPML(IMOL,I)
        SUMMEP=SUMMEP+SUMM2(IMOL,ISTRA)
 3103 CONTINUE
C
      DO 3123 N=1,NSIGSI
        IF (IIHW(N).EQ.12) THEN
          DO 3124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 3124     CONTINUE
        ELSEIF (IIHW(N).EQ.37) THEN
          DO 3624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
 3624     CONTINUE
        ENDIF
 3123 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.18 AND NO.43
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 3104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFPIO) SUMI1(IION,ISTRA)=PRFPIO(IION,I)
        SUMMTP=SUMMTP+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFPIO) SUMI2(IION,ISTRA)=ERFPIO(IION,I)
        SUMMEP=SUMMEP+SUMI2(IION,ISTRA)
 3104 CONTINUE
C
      DO 3125 N=1,NSIGSI
        IF (IIHW(N).EQ.18) THEN
          DO 3126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
 3126     CONTINUE
        ELSEIF (IIHW(N).EQ.43) THEN
          DO 3626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
 3626     CONTINUE
        ENDIF
 3125 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.24 AND NO.49
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFPPHT) SUMPH1(IPHOT,ISTRA)=PRFPPHT(IPHOT,I)
        SUMMTP=SUMMTP+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFPPHT) SUMPH2(IPHOT,ISTRA)=ERFPPHT(IPHOT,I)
        SUMMEP=SUMMEP+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.24) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.49) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0._DP) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT BULK IONS'
        CALL EIRENE_LEER(1)
      ELSE
        CALL EIRENE_LEER(1)
C
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT BULK IONS:'
        WRITE (iunout,*) '(RECYCLING SOURCE) '
C  SURFACE-AVERAGED TALLY NO. 6
        IF (SUMA.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
C  SURFACE-AVERAGED TALLY NO. 31
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMA2,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARA2,LGVRA2,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 12
        IF (SUMM.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
C  SURFACE-AVERAGED TALLY NO. 37
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMM2,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARM2,LGVRM2,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 18
        IF (SUMI.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
C  SURFACE-AVERAGED TALLY NO. 43
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMI2,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARI2,LGVRI2,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        ENDIF
C  SURFACE-AVERAGED TALLY NO. 24
        IF (SUMPH.NE.0._DP) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL EIRENE_MASYR1('P-FLUX:  ',
     .    SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
C  SURFACE-AVERAGED TALLY NO. 49
        CALL EIRENE_MASYR1('E-FLUX:  ',
     .    SUMPH2,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .    VARPH2,LGVRPH2,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        ENDIF
        CALL EIRENE_LEER (1)
        WRITE (iunout,*)
     .    'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
        WRITE (iunout,*)
     .    'CONTRIB. FROM INCIDENT BULK IONS (REC. SOURCE)'
        CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTP)
        CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEP)
C
      ENDIF
C
      IF (ABS(SUMMTA)+ABS(SUMMTM)+ABS(SUMMTI)+
     .    ABS(SUMMTPH).NE.0._DP) THEN
        CALL EIRENE_LEER (1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (iunout,*)
     .      'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0._DP) WRITE (iunout,*)
     .      '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
          CALL EIRENE_LEER (1)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          WRITE (iunout,*)
     .      'TOTAL NEGATIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0._DP) WRITE (iunout,*)
     .      '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL EIRENE_MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL EIRENE_MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
        ENDIF
C     ELSEIF(ILIIN = -3): NOTHING TO BE DONE HERE, NET FLUXES ARE ALREADY PRINTED.
      ENDIF
C
C
C   SPUTTERED FLUXES BY INCIDENT ATOMS
C
      TTSPTA  = 0._DP
      TTSPTM  = 0._DP
      TTSPTI  = 0._DP
      TTSPTPH = 0._DP
      TTSPTP  = 0._DP
      TTSPT   = 0._DP
C
C   SURFACE-AVERAGED TALLY NO.51
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO 202 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTAAT) SUMA1(IATM,ISTRA)=SPTAAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
  202 CONTINUE
C
      DO 221 N=1,NSIGSI
        IF (IIHW(N).EQ.51) THEN
          DO 222 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 222
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
  222     CONTINUE
        ENDIF
  221 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.56
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO 203 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTAML) SUMM1(IMOL,ISTRA)=SPTAML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
  203 CONTINUE
C
      DO 223 N=1,NSIGSI
        IF (IIHW(N).EQ.56) THEN
          DO 224 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 224
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
  224     CONTINUE
        ENDIF
  223 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.61
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO 204 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTAIO) SUMI1(IION,ISTRA)=SPTAIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
  204 CONTINUE
C
      DO 225 N=1,NSIGSI
        IF (IIHW(N).EQ.61) THEN
          DO 226 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 226
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
  226     CONTINUE
        ENDIF
  225 CONTINUE
C
C   SURFACE-AVERAGED TALLY NO.66
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTAPHT) SUMPH1(IPHOT,ISTRA)=SPTAPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.66) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE-AVERAGED TALLY NO.71
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO 205 IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTAPL) SUMP1(IPLS,ISTRA)=SPTAPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
  205 CONTINUE
C
      DO 227 N=1,NSIGSI
        IF (IIHW(N).EQ.71) THEN
          DO 228 IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) GOTO 228
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
  228     CONTINUE
        ENDIF
  227 CONTINUE
C
      TTSPTA  = TTSPTA + SUMMA
      TTSPTM  = TTSPTM + SUMMM
      TTSPTI  = TTSPTI + SUMMI
      TTSPTPH = TTSPTPH + SUMMPH
      TTSPTP  = TTSPTP + SUMMP

      TTTT = ABS(SUMMA)+ABS(SUMMM)+ABS(SUMMI)+ABS(SUMMPH)+ABS(SUMMP)
!      IF (TTTT.EQ.0._DP) THEN
!        CALL EIRENE_LEER(1)
!        CALL EIRENE_MASAGE
!     .  ('NO FLUXES SPUTTERED FROM THIS SURFACE BY INCIDENT ATOMS')
!        CALL EIRENE_LEER(1)
!        GOTO 206
!      ENDIF
      IF ( TTTT > 0._DP) THEN
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'FLUX SPUTTERED FROM SURFACE BY INCIDENT ATOMS:'
      IF (SUMMA.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 51
        CALL EIRENE_MASYR1('ATOMS    ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        IF (NLSPCSPT_ATM) THEN
          N1D=NATMP*NATM
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NATMP,N1D
            IF(SPTAAT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,51),SPTAAT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMM.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 56
        CALL EIRENE_MASYR1('MOLECULES',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        IF (NLSPCSPT_ATM) THEN
          N1D=NMOLP*NATM
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NMOLP,N1D
            IF(SPTAML(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,56),SPTAML(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMI.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 61
        CALL EIRENE_MASYR1('TEST IONS',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        IF (NLSPCSPT_ATM) THEN
          N1D=NIONP*NATM
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NIONP,N1D
            IF(SPTAIO(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,61),SPTAIO(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMPH.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 66
        CALL EIRENE_MASYR1('PHOTONS  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        IF (NLSPCSPT_ATM) THEN
          N1D=NPHOTP*NATM
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPHOTP,N1D
            IF(SPTAPHT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,66),SPTAPHT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMP.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 71
        CALL EIRENE_MASYR1('BULK IONS',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        IF (NLSPCSPT_ATM) THEN
          N1D=NPLSP*NATM
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPLSP,N1D
            IF(SPTAPL(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,71),SPTAPL(K,I)
          ENDDO
        ENDIF
      ENDIF
      CALL EIRENE_LEER (1)
      WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT ATOMS'
!pb      CALL EIRENE_MASR1 ('TOT. FLX',SUMMT)
!pb      CALL EIRENE_LEER(1)
      CALL EIRENE_MASR1 ('SPTATOT ',SPTATOT(I))
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.76) THEN
           CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
        ENDIF
      END DO
      CALL EIRENE_LEER(1)

      ELSEIF (SPTATOT(I) > 0._DP) THEN

        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT ATOMS'
        CALL EIRENE_MASR1 ('SPTATOT ',SPTATOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.76) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)

      END IF   ! TTTT
! 206 CONTINUE
C
C
C   SPUTTERED FLUXES BY INCIDENT MOLECULES
C
C   SURFACE-AVERAGED TALLY NO.52
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTMAT) SUMA1(IATM,ISTRA)=SPTMAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.52) THEN
          DO IATM=1,NATMI
            IF (IGHW(N).NE.IATM) CYCLE
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.57
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTMML) SUMM1(IMOL,ISTRA)=SPTMML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.57) THEN
          DO IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) CYCLE
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.62
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTMIO) SUMI1(IION,ISTRA)=SPTMIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.62) THEN
          DO IION=1,NIONI
            IF (IGHW(N).NE.IION) CYCLE
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.67
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTMPHT) SUMPH1(IPHOT,ISTRA)=SPTMPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.67) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE-AVERAGED TALLY NO.72
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTMPL) SUMP1(IPLS,ISTRA)=SPTMPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.72) THEN
          DO IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) CYCLE
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
          END DO
        ENDIF
      END DO
C
      TTSPTA  = TTSPTA + SUMMA
      TTSPTM  = TTSPTM + SUMMM
      TTSPTI  = TTSPTI + SUMMI
      TTSPTPH = TTSPTPH + SUMMPH
      TTSPTP  = TTSPTP + SUMMP

      TTTT = ABS(SUMMA)+ABS(SUMMM)+ABS(SUMMI)+ABS(SUMMPH)+ABS(SUMMP)

!     IF (TTTT.EQ.0._DP) THEN
!       CALL EIRENE_LEER(1)
!       CALL EIRENE_MASAGE
!    .  ('NO FLUXES SPUTTERED FROM THIS SURFACE BY INCIDENT MOLECULES')
!       CALL EIRENE_LEER(1)
!       GOTO 207
!     ENDIF

      IF (TTTT > 0) THEN
      CALL EIRENE_LEER(1)
      WRITE (iunout,*)
     .   'FLUX SPUTTERED FROM SURFACE BY INCIDENT MOLECULES:'
      IF (SUMMA.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 52
        CALL EIRENE_MASYR1('ATOMS    ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        IF (NLSPCSPT_MOL) THEN
          N1D=NATMP*NMOL
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NATMP,N1D
            IF(SPTMAT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,52),SPTMAT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMM.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 57
        CALL EIRENE_MASYR1('MOLECULES',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        IF (NLSPCSPT_MOL) THEN
          N1D=NMOLP*NMOL
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NMOLP,N1D
            IF(SPTMML(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,57),SPTMML(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMI.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 62
        CALL EIRENE_MASYR1('TEST IONS',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        IF (NLSPCSPT_MOL) THEN
          N1D=NIONP*NMOL
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NIONP,N1D
            IF(SPTMIO(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,62),SPTMIO(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMPH.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 67
        CALL EIRENE_MASYR1('PHOTONS  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        IF (NLSPCSPT_MOL) THEN
          N1D=NPHOTP*NMOL
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPHOTP,N1D
            IF(SPTMPHT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,67),SPTMPHT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMP.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 72
        CALL EIRENE_MASYR1('BULK IONS',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        IF (NLSPCSPT_MOL) THEN
          N1D=NPLSP*NMOL
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPLSP,N1D
            IF(SPTMPL(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,72),SPTMPL(K,I)
          ENDDO
        ENDIF
      ENDIF
      CALL EIRENE_LEER (1)
      WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT MOLECULES '
!pb      CALL EIRENE_MASR1 ('TOT. FLX',SUMMT)
!pb      CALL EIRENE_LEER(1)
      CALL EIRENE_MASR1 ('SPTMTOT ',SPTMTOT(I))
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.77) THEN
           CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
        ENDIF
      END DO
      CALL EIRENE_LEER(1)

      ELSEIF (SPTMTOT(I) > 0._DP) THEN

        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT MOLECULES'
        CALL EIRENE_MASR1 ('SPTMTOT ',SPTMTOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.77) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)

      END IF   !  TTTT
! 207 CONTINUE
C
C
C   SPUTTERED FLUXES BY INCIDENT TEST IONS
C
C   SURFACE-AVERAGED TALLY NO.53
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTIAT) SUMA1(IATM,ISTRA)=SPTIAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.53) THEN
          DO IATM=1,NATMI
            IF (IGHW(N).NE.IATM) CYCLE
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.58
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTIML) SUMM1(IMOL,ISTRA)=SPTIML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.58) THEN
          DO IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) CYCLE
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.63
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTIIO) SUMI1(IION,ISTRA)=SPTIIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.63) THEN
          DO IION=1,NIONI
            IF (IGHW(N).NE.IION) CYCLE
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.68
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTIPHT) SUMPH1(IPHOT,ISTRA)=SPTIPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.68) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE-AVERAGED TALLY NO.73
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTIPL) SUMP1(IPLS,ISTRA)=SPTIPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.73) THEN
          DO IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) CYCLE
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
          END DO
        ENDIF
      END DO
C
      TTSPTA  = TTSPTA + SUMMA
      TTSPTM  = TTSPTM + SUMMM
      TTSPTI  = TTSPTI + SUMMI
      TTSPTPH = TTSPTPH + SUMMPH
      TTSPTP  = TTSPTP + SUMMP

      TTTT = ABS(SUMMA)+ABS(SUMMM)+ABS(SUMMI)+ABS(SUMMPH)+ABS(SUMMP)

!      IF (TTTT.EQ.0._DP) THEN
!        CALL EIRENE_LEER(1)
!        CALL EIRENE_MASAGE
!     .  ('NO FLUXES SPUTTERED FROM THIS SURFACE BY INCIDENT TEST IONS')
!        CALL EIRENE_LEER(1)
!        GOTO 208
!      ENDIF

      IF (TTTT > 0._DP) THEN
      CALL EIRENE_LEER(1)
      WRITE (iunout,*)
     .   'FLUX SPUTTERED FROM SURFACE BY INCIDENT TEST IONS:'
      IF (SUMMA.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 53
        CALL EIRENE_MASYR1('ATOMS    ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        IF (NLSPCSPT_ION) THEN
          N1D=NATMP*NION
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NATMP,N1D
            IF(SPTIAT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,53),SPTIAT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMM.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 58
        CALL EIRENE_MASYR1('MOLECULES',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        IF (NLSPCSPT_ION) THEN
          N1D=NMOLP*NION
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NMOLP,N1D
            IF(SPTIML(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,58),SPTIML(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMI.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 63
        CALL EIRENE_MASYR1('TEST IONS',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        IF (NLSPCSPT_ION) THEN
          N1D=NIONP*NION
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NIONP,N1D
            IF(SPTIIO(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,63),SPTIIO(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMPH.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 68
        CALL EIRENE_MASYR1('PHOTONS  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        IF (NLSPCSPT_ION) THEN
          N1D=NPHOTP*NION
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPHOTP,N1D
            IF(SPTIPHT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,68),SPTIPHT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMP.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 73
        CALL EIRENE_MASYR1('BULK IONS',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        IF (NLSPCSPT_ION) THEN
          N1D=NPLSP*NION
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPLSP,N1D
            IF(SPTIPL(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,73),SPTIPL(K,I)
          ENDDO
        ENDIF
      ENDIF
      CALL EIRENE_LEER (1)
      WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT TEST IONS'
!pb      CALL EIRENE_MASR1 ('TOT. FLX',SUMMT)
!pb      CALL EIRENE_LEER(1)
      CALL EIRENE_MASR1 ('SPTITOT ',SPTITOT(I))
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.78) THEN
           CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
        ENDIF
      END DO
      CALL EIRENE_LEER(1)

      ELSEIF (SPTITOT(I) > 0._DP) THEN

        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT TEST IONS'
        CALL EIRENE_MASR1 ('SPTITOT ',SPTITOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.78) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)

      END IF  !   TTTT
! 208 CONTINUE
C
C
C   SPUTTERED FLUXES BY INCIDENT PHOTONS
C
C   SURFACE-AVERAGED TALLY NO.54
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTPHAT) SUMA1(IATM,ISTRA)=SPTPHAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.54) THEN
          DO IATM=1,NATMI
            IF (IGHW(N).NE.IATM) CYCLE
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.59
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTPHML) SUMM1(IMOL,ISTRA)=SPTPHML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.59) THEN
          DO IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) CYCLE
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.64
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTPHIO) SUMI1(IION,ISTRA)=SPTPHIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.64) THEN
          DO IION=1,NIONI
            IF (IGHW(N).NE.IION) CYCLE
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.69
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTPHPHT) SUMPH1(IPHOT,ISTRA)=SPTPHPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.69) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE-AVERAGED TALLY NO.74
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTPHPL) SUMP1(IPLS,ISTRA)=SPTPHPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.74) THEN
          DO IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) CYCLE
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
          END DO
        ENDIF
      END DO
C
      TTSPTA  = TTSPTA + SUMMA
      TTSPTM  = TTSPTM + SUMMM
      TTSPTI  = TTSPTI + SUMMI
      TTSPTPH = TTSPTPH + SUMMPH
      TTSPTP  = TTSPTP + SUMMP

      TTTT = ABS(SUMMA)+ABS(SUMMM)+ABS(SUMMI)+ABS(SUMMPH)+ABS(SUMMP)

!      IF (TTTT.EQ.0._DP) THEN
!        CALL EIRENE_LEER(1)
!        CALL EIRENE_MASAGE
!     .  ('NO FLUXES SPUTTERED FROM THIS SURFACE BY INCIDENT PHOTONS')
!        CALL EIRENE_LEER(1)
!        GOTO 209
!      ENDIF

      IF (TTTT > 0._DP) THEN
      CALL EIRENE_LEER(1)
      WRITE (iunout,*)
     .   'FLUX SPUTTERED FROM SURFACE BY INCIDENT PHOTONS:'
      IF (SUMMA.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 54
        CALL EIRENE_MASYR1('ATOMS    ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        IF (NLSPCSPT_PHOT) THEN
          N1D=NATMP*NPHOT
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NATMP,N1D
            IF(SPTPHAT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,54),SPTPHAT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMM.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 59
        CALL EIRENE_MASYR1('MOLECULES',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        IF (NLSPCSPT_PHOT) THEN
          N1D=NMOLP*NPHOT
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NMOLP,N1D
            IF(SPTPHML(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,59),SPTPHML(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMI.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 64
        CALL EIRENE_MASYR1('TEST IONS',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        IF (NLSPCSPT_PHOT) THEN
          N1D=NIONP*NPHOT
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NIONP,N1D
            IF(SPTPHIO(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,64),SPTPHIO(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMPH.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 69
        CALL EIRENE_MASYR1('PHOTONS  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        IF (NLSPCSPT_PHOT) THEN
          N1D=NPHOTP*NPHOT
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPHOTP,N1D
            IF(SPTPHPHT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,69),SPTPHPHT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMP.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 74
        CALL EIRENE_MASYR1('BULK IONS',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        IF (NLSPCSPT_PHOT) THEN
          N1D=NPLSP*NPHOT
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPLSP,N1D
            IF(SPTPHPL(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,74),SPTPHPL(K,I)
          ENDDO
        ENDIF
      ENDIF
      CALL EIRENE_LEER (1)
      WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT PHOTONS '
!pb      CALL EIRENE_MASR1 ('TOT. FLX',SUMMT)
!pb      CALL EIRENE_LEER(1)
      CALL EIRENE_MASR1 ('SPTPHTOT',SPTPHTOT(I))
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.79) THEN
           CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
        ENDIF
      END DO
      CALL EIRENE_LEER(1)

      ELSEIF (SPTPHTOT(I) > 0._DP) THEN

        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT PHOTONS'
        CALL EIRENE_MASR1 ('SPTPHTOT',SPTPHTOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.79) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)

      END IF  ! TTTT
! 209 CONTINUE
C
C
C   SPUTTERED FLUXES BY INCIDENT BULK IONS
C
C   SURFACE-AVERAGED TALLY NO.55
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTPAT) SUMA1(IATM,ISTRA)=SPTPAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.55) THEN
          DO IATM=1,NATMI
            IF (IGHW(N).NE.IATM) CYCLE
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.60
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTPML) SUMM1(IMOL,ISTRA)=SPTPML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.60) THEN
          DO IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) CYCLE
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.65
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTPIO) SUMI1(IION,ISTRA)=SPTPIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.65) THEN
          DO IION=1,NIONI
            IF (IGHW(N).NE.IION) CYCLE
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
          END DO
        ENDIF
      END DO
C
C   SURFACE-AVERAGED TALLY NO.70
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTPPHT) SUMPH1(IPHOT,ISTRA)=SPTPPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.70) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE-AVERAGED TALLY NO.75
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTPPL) SUMP1(IPLS,ISTRA)=SPTPPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
      END DO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.75) THEN
          DO IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) CYCLE
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
          END DO
        ENDIF
      END DO
C
      TTSPTA  = TTSPTA + SUMMA
      TTSPTM  = TTSPTM + SUMMM
      TTSPTI  = TTSPTI + SUMMI
      TTSPTPH = TTSPTPH + SUMMPH
      TTSPTP  = TTSPTP + SUMMP

      TTTT = ABS(SUMMA)+ABS(SUMMM)+ABS(SUMMI)+ABS(SUMMPH)+ABS(SUMMP)

!      IF (TTTT.EQ.0._DP) THEN
!        CALL EIRENE_LEER(1)
!        CALL EIRENE_MASAGE
!     .  ('NO FLUXES SPUTTERED FROM THIS SURFACE BY INCIDENT BULK IONS')
!        CALL EIRENE_LEER(1)
!        GOTO 210
!      ENDIF

      IF (TTTT > 0._DP) THEN
      CALL EIRENE_LEER(1)
      WRITE (iunout,*)
     .   'FLUX SPUTTERED FROM SURFACE BY INCIDENT BULK IONS'
      IF (SUMMA.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 55
        CALL EIRENE_MASYR1('ATOMS    ',
     .       SUMA1,LOGATM,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARA1,LGVRA1,ISTRA,0,NATM,NATM,0,NSTRA,TEXTS(NSPH+1))
        IF (NLSPCSPT_PLS) THEN
          N1D=NATMP*NPLS
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NATMP,N1D
            IF(SPTPAT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,55),SPTPAT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMM.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 60
        CALL EIRENE_MASYR1('MOLECULES',
     .       SUMM1,LOGMOL,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARM1,LGVRM1,ISTRA,0,NMOL,NMOL,0,NSTRA,TEXTS(NSPA+1))
        IF (NLSPCSPT_PLS) THEN
          N1D=NMOLP*NPLS
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NMOLP,N1D
            IF(SPTPML(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,60),SPTPML(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMI.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 65
        CALL EIRENE_MASYR1('TEST IONS',
     .       SUMI1,LOGION,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARI1,LGVRI1,ISTRA,0,NION,NION,0,NSTRA,TEXTS(NSPAM+1))
        IF (NLSPCSPT_PLS) THEN
          N1D=NIONP*NPLS
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NIONP,N1D
            IF(SPTPIO(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,65),SPTPIO(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMPH.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 70
        CALL EIRENE_MASYR1('PHOTONS  ',
     .       SUMPH1,LOGPHOT,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARPH1,LGVRPH1,ISTRA,0,NPHOT,NPHOT,0,NSTRA,TEXTS(0+1))
        IF (NLSPCSPT_PLS) THEN
          N1D=NPHOTP*NPLS
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPHOTP,N1D
            IF(SPTPPHT(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                        TXTSPW(K,70),SPTPPHT(K,I)
          ENDDO
        ENDIF
      ENDIF
      IF (SUMMP.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 75
        CALL EIRENE_MASYR1('BULK IONS',
     .       SUMP1,LOGPLS,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARP1,LGVRP1,ISTRA,0,NPLS,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        IF (NLSPCSPT_PLS) THEN
          N1D=NPLSP*NPLS
          WRITE(iunout,*) 'SPECIES RESOLVED CONTRIBUTIONS (not scaled):'
          DO K=NPLSP,N1D
            IF(SPTPPL(K,I)>0._DP) WRITE(iunout,'(1X,1A24,1X,1ES12.4)')
     .                                       TXTSPW(K,75),SPTPPL(K,I)
          ENDDO
        ENDIF
      ENDIF
      CALL EIRENE_LEER (1)
      WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT BULK IONS'
!pb      CALL EIRENE_MASR1 ('TOT. FLX',SUMMT)
!pb      CALL EIRENE_LEER(1)
      CALL EIRENE_MASR1 ('SPTPLTOT',SPTPLTOT(I))
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.80) THEN
           CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
        ENDIF
      END DO
      CALL EIRENE_LEER(1)

      ELSEIF (SPTPLTOT(I) > 0._DP) THEN

        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED BY INCIDENT BULK IONS'
        CALL EIRENE_MASR1 ('SPTPLTOT',SPTPLTOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.80) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)

      END IF  ! TTTT
! 210 CONTINUE

      TTSPT = TTSPTA + TTSPTM + TTSPTI + TTSPTPH + TTSPTP

      IF (ABS(TTSPT) < EPS10.AND.SPTTOT(I) < EPS10) THEN
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASAGE('NO FLUXES SPUTTERED FROM THIS SURFACE')
        CALL EIRENE_LEER(1)
      ELSEIF (ABS(TTSPT) > EPS10) THEN
        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'FLUX SPUTTERED AS'
        IF (ABS(TTSPTA) > EPS10) CALL EIRENE_MASR1 ('ATOMS   ',TTSPTA)
        IF (ABS(TTSPTM) > EPS10) CALL EIRENE_MASR1 ('MOLEC.  ',TTSPTM)
        IF (ABS(TTSPTI) > EPS10) CALL EIRENE_MASR1 ('TESTIONS',TTSPTI)
        IF (ABS(TTSPTPH) > EPS10)CALL EIRENE_MASR1 ('PHOTONS ',TTSPTPH)
        IF (ABS(TTSPTP) > EPS10) CALL EIRENE_MASR1 ('BULKIONS',TTSPTP)

        CALL EIRENE_LEER(1)
        WRITE (IUNOUT,*) 'TOTAL FLUX SPUTTERED FROM SURFACE'
        CALL EIRENE_MASR1 ('TOT. FLX',TTSPT)
      END IF

      IF (SPTTOT(I) > 0._DP.and.spttot(i).ne.ttspt) THEN
        CALL EIRENE_LEER (1)
        WRITE (IUNOUT,*) 'TOT. FLX SPUTTERED'
        write (iunout,*) '(not scaled by NLSCL option)'
        write (iunout,*) 'May include unidentified sputtered species'
        CALL EIRENE_MASR1 ('TOT. FLX',SPTTOT(I))
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.81) THEN
             CALL EIRENE_MASR1 ('ST.DEV.%',SIGMAW(N,I))
          ENDIF
        END DO
        CALL EIRENE_LEER(1)
      ENDIF

C
C   ADDITIONAL SURFACE-AVERAGED TALLIES
C
C   SURFACE-AVERAGED TALLY NO. 82 (=NTLSA)
C
      SUMMS=0.
      SUMS(:,ISTRA)=0._DP
      DO 302 IADS=1,NADSI
! TEXTA is only 8 characters long
        TEXTA(IADS)=TXTSPW(IADS,NTLSA)(1:8)
        LGVARS(IADS,ISTRA)=.FALSE.
        IF (LADDS) SUMS(IADS,ISTRA)=ADDS(IADS,I)
        LOGADS(IADS,ISTRA)=LADDS .AND. (ADDS(IADS,I).NE.0.)
        SUMMS=SUMMS+SUMS(IADS,ISTRA)
  302 CONTINUE
C
      DO 321 N=1,NSIGSI
        IF (IIHW(N).NE.82) GOTO 321
        IF (IGHW(N).EQ.0) GOTO 321
        DO 322 IADS=1,NADSI
          IF (IGHW(N).NE.IADS) GOTO 322
          VARS(IADS,ISTRA)=SIGMAW(N,I)
          LGVARS(IADS,ISTRA)=LOGADS(IADS,ISTRA)
  322   CONTINUE
  321 CONTINUE
C
      IF (SUMMS.EQ.0._DP) THEN
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASAGE
     .  ('NO ADDITIONAL SURFACE TALLIES AT THIS SURFACE')
        CALL EIRENE_LEER(1)
        GOTO 305
      ENDIF
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'ADDITIONAL SURFACE TALLIES:'
      IF (SUMMS.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 82
        CALL EIRENE_MASYR1('ADD.TALLY',
     .       SUMS,LOGADS,ISTRA,0,NADS,NADSI,0,NSTRA,TEXTA)
        CALL EIRENE_MASYR1('ST.DEV.% ',
     .       VARS,LGVARS,ISTRA,0,NADS,NADSI,0,NSTRA,TEXTA)
      ENDIF
      CALL EIRENE_LEER (1)
C
  305 CONTINUE
C
C
C   ALGEBRAIC SURFACE-AVERAGED TALLIES
C
C   SURFACE-AVERAGED TALLY NO. 83
C
      SUMMS=0.
      SUML(:,ISTRA)=0._DP
      DO 402 IALS=1,NALSI
! TEXTL is only 8 characters long
        TEXTL(IALS)=TXTSPW(IALS,NTLSR)(1:8)
        LGVARL(IALS,ISTRA)=.FALSE.
        IF (LALGS) SUML(IALS,ISTRA)=ALGS(IALS,I)
        LOGALS(IALS,ISTRA)=LALGS .AND. (ALGS(IALS,I).NE.0.)
        SUMMS=SUMMS+SUML(IALS,ISTRA)
  402 CONTINUE
C
      DO 421 N=1,NSIGSI
        IF (IIHW(N).NE.83) GOTO 421
        IF (IGHW(N).EQ.0) GOTO 421
        DO 422 IALS=1,NALSI
          IF (IGHW(N).NE.IALS) GOTO 422
          VARL(IALS,ISTRA)=SIGMAW(N,I)
          LGVARL(IALS,ISTRA)=LOGALS(IALS,ISTRA)
  422   CONTINUE
  421 CONTINUE
C
      IF (SUMMS.EQ.0._DP) THEN
        CALL EIRENE_LEER(1)
        CALL EIRENE_MASAGE
     .  ('NO ALGEBRAIC SURFACE TALLIES AT THIS SURFACE')
        CALL EIRENE_LEER(1)
        GOTO 405
      ENDIF
      CALL EIRENE_LEER(1)
      WRITE (iunout,*) 'ALGEBRAIC SURFACE TALLIES:'
      DO IALS=1,NALSI
        WRITE (iunout,*) 'NO ',IALS,'  ',TXTTLW(IALS,NTLSR)
      ENDDO
      IF (SUMMS.NE.0._DP) THEN
C  SURFACE-AVERAGED TALLY NO. 83
        CALL EIRENE_MASYR1
     .     ('         ',SUML,LOGALS,ISTRA,0,NALS,NALS,0,NSTRA,TEXTL)
      ENDIF
      CALL EIRENE_LEER (1)
C
  405 CONTINUE
C
C   SURFACE AVERAGED SPECTRA, INTEGRALS, SURFACE I, ISP
      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO ISPC=1,NADSPC
        IF ((ESTIML(ISPC)%ISRFCLL == 0) .AND.
     .      (ESTIML(ISPC)%ISPCSRF == I)) THEN
          CALL EIRENE_LEER (1)
          WRITE (iunout,'(A33)') ' SPECTRUM CALCULATED FOR SURFACE '
          IT = ESTIML(ISPC)%ISPCTYP
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
          END IF

          WRITE (iunout,'(A20,A8)') ' TYPE OF PARTICLE : ',
     .                       TEXTYP(ESTIML(ISPC)%IPRTYP)
          IF (ESTIML(ISPC)%IPRSP == 0) THEN
            WRITE (iunout,'(A10,10X,A16)') ' SPECIES :',
     .                                     'SUM OVER SPECIES'
          ELSE
            WRITE (iunout,'(A10,10X,A8)') ' SPECIES :',
     .             TEXTS(IADTYP(ESTIML(ISPC)%IPRTYP)+
     .                   ESTIML(ISPC)%IPRSP)
          END IF
          WRITE (iunout,'(A22,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .           ESTIML(ISPC)%SPCS
          IF (NSIGI_SPC > 0)
     .      WRITE (iunout,'(A22,ES12.4)') ' STANDARD DEVIATION   ',
     .           ESTIML(ISPC)%SGMS
        END IF
      END DO

        PRINTED(I) = .TRUE.
10000 CONTINUE  ! END OF LOOP OVER SURFACES I OR IST,
                ! FOR WHICH PRINTOUT WAS REQUESTED
      RETURN
 9999 FORMAT (1X,A32)

      CONTAINS

      SUBROUTINE EIRENE_PRINT_SURF_TRIAN(IUN)

      INTEGER, INTENT(IN) :: IUN
      INTEGER :: II, KK, IT, IS, IS1, I1, I2
      CHARACTER(1) :: TL(72)
      CHARACTER(14) :: FORM, FORMA

      DATA TL/72*'='/

      CALL EIRENE_LEER(3)
      WRITE (iun,*) TL
      WRITE (iun,*) 'TALLY:   ',TXTTLW(NFTI,ITALS)
      WRITE (iun,*) TL
      CALL EIRENE_LEER(1)

      FORM = '(I3,   ES12.4)'
      WRITE (FORM(5:7),'(I3)') NFTE-NFTI+1+5

      WRITE (IUN,'(63x,A)') 'SPECIES'
      FORMA = '(A3,   A12)'
      WRITE (FORMA(5:7),'(I3)') NFTE-NFTI+1+5
      WRITE (IUN,FORMA) 'I','XA','YA','XE','YE','ARC LENGTH',
     .                  (TRIM(TXTSPW(KK,ITALS)),KK=NFTI,NFTE)
      WRITE (IUN,FORMA) ' ','CM','CM','CM','CM','CM',
     .                  (TRIM(TXTUNW(KK,ITALS)),KK=NFTI,NFTE)

      DO II = 1, NTOTAL-1
         IT = SURF_TRIAN_ORDERED(I)%ITRIAS(II)
         IS = SURF_TRIAN_ORDERED(I)%ITRISI(II)
         IS1 = MOD(IS,3) + 1
         I1 = NECKE(IS,IT)
         I2 = NECKE(IS1,IT)
         WRITE (IUN,FORM)
     .      II, XTRIAN(I1), YTRIAN(I1), XTRIAN(i2), YTRIAN(I2),
     .      SURF_TRIAN_ORDERED(I)%BGLT(II+1),
     .      (HELP2(II,KK),KK=NFTI,NFTE)
      END DO

      WRITE (IUN,'(A)') 'TOTALS'
      FORM = '(63X,  ES12.4)'
      WRITE (FORM(6:7),'(I2)') NFTE-NFTI+1
      WRITE (IUN,FORM) (TOTAR(KK),KK=NFTI,NFTE)

      RETURN
      END SUBROUTINE EIRENE_PRINT_SURF_TRIAN

      END SUBROUTINE EIRENE_OUTFLX
