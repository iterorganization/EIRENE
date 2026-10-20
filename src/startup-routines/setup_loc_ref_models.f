
      SUBROUTINE EIRENE_SETUP_LOC_REF_MODELS
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CREF
      USE EIRMOD_CLGIN
      USE EIRMOD_COMPRT, ONLY : IUNOUT
      USE EIRMOD_PRESSURELOOP

      IMPLICIT NONE

      INTEGER :: J, JSPZ, NLJ, JJ, EIRENE_IDEZ
      REAL(DP) :: RSAVE
      TYPE(TSURFACE), POINTER :: SURFCUR, SURFCUR2
      TYPE(REFMODEL), POINTER :: REFCUR
      EXTERNAL :: EIRENE_EXIT_OWN, EIRENE_IDEZ
C
C  SOURCE PARAMETERS AND (REFLECTING) BOUNDARY CONDITIONS,
C  ON ADDITIONAL AND NON-DEFAULT STANDARD SURFACES
C
      DO J=0,NLIMPS
C
        RINTEG(J)=RINTEG(1)
        EINTEG(J)=EINTEG(1)
        AINTEG(J)=AINTEG(1)
        DO JSPZ=1,NSPZ
Crc changes to ilref, not sure if correct here
          ILREF(JSPZ,J)=ILREF(1,J)
          ISRS(JSPZ,J)=ISRS(1,J)
          ISRC(JSPZ,J)=ISRC(1,J)
          LCHSPNWL(JSPZ,J)=LCHSPNWL(1,J)
          TRANSP(JSPZ,1,J)=TRANSP(1,1,J)
          TRANSP(JSPZ,2,J)=TRANSP(1,2,J)
          RECYCF(JSPZ,J)=RECYCF(1,J)
          RECYCT(JSPZ,J)=RECYCT(1,J)
          RECPRM(JSPZ,J)=RECPRM(1,J)
          EXPPL(JSPZ,J)=EXPPL(1,J)
          EXPEL(JSPZ,J)=EXPEL(1,J)
          EXPIL(JSPZ,J)=EXPIL(1,J)
          RECYCS(JSPZ,J)=RECYCS(1,J)
          RECYCC(JSPZ,J)=RECYCC(1,J)
          SPTPRM(JSPZ,J)=SPTPRM(1,J)
          ESPUTS(JSPZ,J)=ESPUTS(1,J)
          ESPUTC(JSPZ,J)=ESPUTC(1,J)
        end do
      end do

!PB  REFLIST NEEDS TO BE KEPT FOR JSON OUTPUT
      REFCUR => REFLIST
      DO WHILE (ASSOCIATED(REFCUR))
        NULLIFY(SURFCUR2)
        SURFCUR => SURFLIST
        DO WHILE (ASSOCIATED(SURFCUR))
          IF (SURFCUR%MODNAME == REFCUR%REFNAME) THEN
            NLJ = SURFCUR%NOSURF
            ILREF(:,NLJ) = REFCUR%JLREF
            ILSPT(NLJ) = REFCUR%JLSPT
            ISRS(:,NLJ) = REFCUR%JSRS
            ISRC(:,NLJ) = REFCUR%JSRC
            LCHSPNWL(:,NLJ) = REFCUR%JLCHSPNWL
            ZNML(NLJ) = REFCUR%ZNMLR
            EWALL(NLJ) = REFCUR%EWALLR
            EWBIN(NLJ) = REFCUR%EWBINR
            TRANSP(:,1,NLJ) = REFCUR%TRANSPR(:,1)
            TRANSP(:,2,NLJ) = REFCUR%TRANSPR(:,2)
            FSHEAT(NLJ) = REFCUR%FSHEATR
            RECYCF(:,NLJ) = REFCUR%RCYCFR
            RECYCT(:,NLJ) = REFCUR%RCYCTR
            RECPRM(:,NLJ) = REFCUR%RCPRMR
            EXPPL(:,NLJ) = REFCUR%EXPPLR
            EXPEL(:,NLJ) = REFCUR%EXPELR
            EXPIL(:,NLJ) = REFCUR%EXPILR
            RECYCS(:,NLJ) = REFCUR%RCYCSR
            RECYCC(:,NLJ) = REFCUR%RCYCCR
            SPTPRM(:,NLJ) = REFCUR%STPRMR
            ESPUTS(:,NLJ) = REFCUR%ESPTSR
            ESPUTC(:,NLJ) = REFCUR%ESPTCR
crc changes to ilref, not sure if correct here, temporary solution
            !Initialize pressure feedback loop
            IF (ILREF(1,NLJ) == 4)THEN
              CALL initPressureFeedback(RPRESSFED(NLJ), NLJ,
     .                                  REFCUR%REFCELL,
     .                                  REFCUR%REFPRESS)
            END IF

            IF (.NOT.ASSOCIATED(SURFCUR2)) THEN
              SURFLIST => SURFCUR%NEXT
              DEALLOCATE(SURFCUR)
              SURFCUR => SURFLIST
            ELSE
              SURFCUR2%NEXT => SURFCUR%NEXT
              DEALLOCATE(SURFCUR)
              SURFCUR => SURFCUR2%NEXT
            END IF
          ELSE
            SURFCUR2 => SURFCUR
            SURFCUR => SURFCUR%NEXT
          END IF
        END DO
        REFCUR => REFCUR%NEXT
      ENDDO

      IF (ASSOCIATED(SURFLIST)) THEN
        WRITE (iunout,*)
     .    ' SURFACE DATA HAVE NOT BEEN DEFINED FOR MODEL:'
        DO WHILE (ASSOCIATED(SURFLIST))
          WRITE (iunout,*) SURFLIST%MODNAME
          SURFCUR => SURFLIST
          SURFLIST => SURFLIST%NEXT
          DEALLOCATE(SURFCUR)
        END DO
        WRITE (iunout,*) ' EXECUTION IS STOPPED'
        CALL EIRENE_EXIT_OWN(1)
      END IF

      DO 2000 J=0,NLIMPS
c  non-default standard surfaces have a negative surface index on printout
        JJ=J
        IF (JJ.GT.NLIM) JJ=-(J-NLIM)

        IF (ILCOL(J).LT.0) IGFIL(J)=1
        ILCOL(J)=MAX0(1,IABS(ILCOL(J)))

        IF ((ILIIN(J).LE.0.OR.ILIIN(J).GE.4).AND.(ILSPT(J).NE.0)) THEN
          WRITE (IUNOUT,'(1X,A,I4,A)') 'WARNING: SURFACE NO. ',JJ,
     .      ' IS TRANSPARENT OR PERIODIC BUT IS USED TO SPUTTER'
        END IF

        ISPUT(1,J)=EIRENE_IDEZ(ILSPT(J),1,2)
        ISPUT(2,J)=EIRENE_IDEZ(ILSPT(J),2,2)

        IF (ILIIN(J).LE.0.OR.ILIIN(J).GE.3) ILSPT(J)=0

        IF (ILIIN(J).LE.0) TRANSP(:,1,J)=0.D0
        IF (ILIIN(J).LE.0) TRANSP(:,2,J)=0.D0

        IF (ISPUT(1,J).NE.0.AND.ALL(RECYCS(:,J).EQ.0._DP)) THEN
          WRITE (IUNOUT,'(1X,A,I4,1x,A,A)') 'WARNING: SURFACE NO. ',JJ,
     .                     'HAS A PHYSICAL SPUTTER MODEL, BUT ',
     .                     'SPUTTER YIELD IS SCALED TO ZERO (RECYCS=0)'
        ENDIF
        IF (ISPUT(2,J).NE.0.AND.ALL(RECYCC(:,J).EQ.0._DP)) THEN
          WRITE (IUNOUT,'(1X,A,I4,1x,A,A)') 'WARNING: SURFACE NO. ',JJ,
     .                     'HAS A CHEMICAL SPUTTER MODEL, BUT ',
     .                     'SPUTTER YIELD IS SCALED TO ZERO (RECYCC=0)'
        ENDIF
        IF (ILIIN(J).EQ.2) RECYCF(:,J)=0.
        IF (ILIIN(J).EQ.2) RECYCT(:,J)=0.
C
        RSAVE=ZNML(J)
        ZNML(J)=DBLE(INT(RSAVE/100.D0))
        ZNCL(J)=RSAVE-100.*ZNML(J)
        XMLIM(J)=ZNML(J) !VK - for consistency with AK code
        XCLIM(J)=ZNCL(J) !VK - ditto
        DO 2001 JSPZ=1,NSPZ
          ISRF(JSPZ,J)=ISRF(JSPZ,1)
          ISRT(JSPZ,J)=ISRT(JSPZ,1)
 2001   CONTINUE
 2000 CONTINUE

      RETURN
      END SUBROUTINE EIRENE_SETUP_LOC_REF_MODELS
