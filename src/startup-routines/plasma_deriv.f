c  Set derived plasma data, such as electron density, vacuum flags, etc.
c  also:
c  "density models" to contruct background data from other given data :
c      Saha, Boltzmann, Planck,
c      corona, col-rad, file (fort.13 or fort.10)
c
c  presently:  Options "File" and "Boltzmann": may affect electron density.
c              hence: done prior to electron density, etc...
c              "Corona", "Colrad", "Saha", "Planck": need electron density as
c                            input, or, at least, do not affect n_e
c                            hence: done after electron density, etc...
C  may05
c  1) additional density model only for neutrals?  removed
c  2) Boltzmann factor only if Ti gt Tvac
c     to be checked: correct low T limit: everything in lower level?
c  3) new ti only if nlmlti
c  4) new vi only if nlmlv (nlmlv is new, in cinit.f, and set in input.f)
c  june05
c     bug fix: in colrad model: vnew=max(vvac,vold) nonsense, because
c          vold<0 possible. replaced by vnew=vold
c
c  march 06
c     new option: icall > 0, and call base_parameter
c       allows to use output (test particle) tallies and special "density model" to
c       construct new input (field particle) tallies (densities, temperatures, drift velocities)
c       e.g. for postprocessing (diagno), or for iterations (bgk).
c
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  06.03.07: new density models 'CONSTANT' and 'MULTIPLY' introduced
!pb            'CONSTANT' sets constant plasma profiles,
cdr             to a value specified in input block 5
!pb            'MULTIPLY' creates a new field particle density by multiplying an
!pb            existing plasma density
cdr            with a factor specified in input block 5

!pb  11.01.10: interpolation of plasma profiles to cell vertices added

cdr: may 2015
cdr: output tallies for new background: in case of multiple strata: how to get sum over strata?
cdr: do we need fort.10 ?
cdr: in case indpro(4)=8:  edrift, vdion:  only for ipls=1 available?
cdr: warnings in case of missing edrift removed: have been too many (one per cell)
cdr: jan 2016: automated resetting of nfilel to =3 or =9 removed.
cdr:           (had caused problems with t-dep mode)
cdr:           should be done more explicitly, by problem-specific routines,
cdr:           or mod_bgk, mod_timstep,.....
cdr:           postprocessing Balmer lines, etc:  to be confirmed that this now
cdr:           still works properly
cdr: jan 2017: generalized asymptotics options for A&M data structures included
cdr:           cleanup. logical FOUND seems to be redundant
cdr  oct 18  : NREAC+1 --> IRC  (A&M data for density models)
cdr            Calls to slreac from block 5 and 12 removed,
cdr            this is now already done in block 4 (plus call to SLREAC there).
cdr            After that,
cdr            A&M data are already on data structure REACDAT(IRC......)
cdr            i.e. A&M data for density models, such as corona or colrad
cdr            excitation rates (EI), or population levels (OT)....
cdr            These must now already be specified in block 4,
cdr            with the new, unified, parser for
cdr            reaction decks: READ_REACLINES.f.
cdr  may 2019: spectra here ? probaby wrong place. (corona, colrad ?).
cdr            unused? remove ?
cdr  Nov. 19:  add planckian (photon gas) density model (unfinished)
cdr  Oct. 20:  remove call to wrplam, i.e. decouple plasma background
cdr            preparation from data file handling, reduce complexity
cdr  Nov. 20:  remove call to slreac from here. All reaction parameters
cdr            for "density models" (Corona, Colrad) are already read
cdr            in block 4. Only use reaction number IRC here.
c
      SUBROUTINE EIRENE_PLASMA_DERIV (ICALL)

c  input:
c    nlmlti (via cinit.f): all bulk ions have own temperature Ti, on Ti(iplsti),
c                          set new Ti for ipls
c    nlmlv  (via cinit.f): all bulk ions have own velocity, Vx,Vy,Vz, on V*(iplsv)
c                          set new flow velocity for ipls

c    icall:

c    icall=0
c      pre-processing:
c      called PRIOR to Monte Carlo loop (from subr. input).
c      In this call all "density models" referring to output tallies
c      are ignored (e.g. 'fort.10').
c    icall=1
c      post-processing:
c      called AFTER Monte Carlo loop and sum over strata.
c        This allows to put output tallies (from fort.10) from a run onto the
c        background (input tallies) for a next iteration or for post-processing.
c        In this call all "density models" referring to input tallies are
c        ignored, because they are already done in the previous call with ICALL=0

c  for appropriate values of nfilel: = 1,3,4,6,8,9
c      write fort.13 (CALL WRPLAM) after all density models are done.

c   carry out specific "background models",
c   for bulk species IPLS
c      'fort.13': take background data from fort.13, species: IOLD
c      'fort.10': take test particle data from fort.10, species: IOLD

c  set derived plasma parameters:
c   DEIN             : electron density (from quasi-neutrality)
c   DEINL            : log electron density (with cutoffs)
c   TEINL            : log electron temperature (with cutoffs)
c   LGVAC(...,NPLS+1): electron vacuum flag
c   LGVAC(...,IPLS)  : bulk "ion" vacuum flag
c   LGVAC(...,0)     : background vacuum flag
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CTRCEI
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CGEOM
      USE EIRMOD_CSPEI
      USE EIRMOD_COMXS
      USE EIRMOD_CESTIM
      use EIRMOD_csdvi
      use EIRMOD_comsou
      use EIRMOD_cspei
      USE EIRMOD_SECOND_OWN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICALL
      REAL(DP) :: ZTII, ZTNI, FCT2, FCRG, FCT1, EIRENE_VDION, ZTEI,
     .            ZTNE,EMPLS, FCT0, TEPLS, DEPLS, DIPLS, AM1, TEF, DEF,
     .            BOLTZFAC, RCORONA, RCOLRAD, DELTAE,
     .            G_BOLTZ, G_PLANCK,
! rates and population coefficients for density models
     .            EIRENE_RATE_COEFF, EIRENE_OTHER_RATE_COEFF,
! asymptotics thereof
     .            RC1MIN, RC1MAX, RC2MIN, RC2MAX,
     .            BXP, BYP, BNORM, TE, DE
      REAL(DP) :: FP1(6), FP2(6)
cdr  Cartesian (contravariant) components of
cdr  parallel, grad-PSI and diamagn. B field, resp.
      REAL(DP) :: BPAR(3), BPER(3), BCROSS(3)
      REAL(DP) :: tpb1, tpb2
      REAL(DP), ALLOCATABLE :: DEINTF(:), SUMNI(:), SUMMNI(:),
     .                         BASE_DENSITY(:), BASE_TEMP(:),
     .                         BASE_VELOC(:,:),
     .                         TALLY(:)
      INTEGER :: IR, IN, IP, IPM, IPLS, IOLD, IRE,
     .           I, J, IAIN, ISPZ,
     .           ICO, KK,
     .           IPLSTI, IPLSV, IOLDTI, IOLDV, IBS, IFLG,
     .           JFEX1MN, JFEX1MX, JFEX2MN, JFEX2MX,
     .           ITAL, K, NFTI, NFTE, JPLS, IRC
cdr june 23: testing b_perp via input, rather than default
      real(dp) :: bxpp,bypp,bpp,dist1,dist2
      integer  :: ico1,ico2
      EXTERNAL :: EIRENE_RATE_COEFF, EIRENE_OTHER_RATE_COEFF,
     .            EIRENE_VDION

      TYPE(EIRENE_SPECTRUM) :: SPEC
      LOGICAL :: FOUND

      interface
        subroutine eirene_cell_to_corner (f,fcorner)
          use eirmod_precision
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fcorner(:)
        end subroutine eirene_cell_to_corner

        subroutine eirene_calc_grad (f, fdx, fdy, fdz, lfdx, lfdy, lfdz)
          use eirmod_precision
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fdx(:), fdy(:), fdz(:)
          logical, intent(in) :: lfdx, lfdy, lfdz
        end subroutine eirene_calc_grad

      end interface

      EXTERNAL :: EIRENE_RPLAM, EIRENE_EXIT_OWN

      FP1 = 0._DP
      FP2 = 0._DP
      RC1MIN = -HUGE(1._DP)
      RC1MAX =  HUGE(1._DP)
      RC2MIN = -HUGE(1._DP)
      RC2MAX =  HUGE(1._DP)
      JFEX1MN = 0
      JFEX1MX = 0
      JFEX2MN = 0
      JFEX2MX = 0

cdr
      write (iunout,*) 'plasma deriv called ', icall
      tpb1 = EIRENE_second_own()
      IBS = 0
cdr  hidden link: this seems to be assuming that nlshrt13=f.
cdr  But nlshrt13 is set true e.g. in couple_solps_iter
cdr  In that case only the virtual species IPLS=NPLS_FIX+1,NPLS are read.
cdr  IFLG=10:  set pointers DIINTF,... in RPLAM
      IFLG=10
      IF (ANY(CDENMODEL == FORT//'13'))
     .     CALL EIRENE_RPLAM(TRCFLE,IFLG,'PLASMA_DERIV')

      DO JPLS=1,NPLSI
        IPLS = JPLS
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
c       write (iunout,*) jpls, cdenmodel(IPLS)

        SELECT CASE (CDENMODEL(IPLS))

          CASE ('FORT.13','FTN13')

cdr fort.13 has been read, but do we have the long or short version?

            IOLD=TDMPAR(IPLS)%TDM%ISP(1)
c           ITOLD=TDMPAR(IPLS)%TDM%ITP(1) =4,  hard-wired
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)

cdr Feb 2020:
cdr hidden link: This next coding assumes that pointers diintf, tiinft,...
cdr              for IOLD are set. This is not necessarily the case, e.g
cdr              in stand alone cases, and if NLSHRT13=.false.
cdr              Or in case NLSHRT13 and IOLD < NPLS_FIX+1.
cdr Check if pointer DIINTF(IOLD,..) is set at all.
cdr Taken to be representative for species IOLD.

cdr to be checked first: allocated: plasma_bckgrnd, and: associated: diintf

            IF (any(diintf(iold,:).ne. 0._DP))  then

              IF (NLMLTI) THEN
                TIIN(IPLSTI,:)=MAX(TVAC,TIINTF(IOLDTI,:))
CDR           ELSE
CDR   in single TI cases: all TIIN(IPLS,:) are the same, and == TIIN(1,:)
              ENDIF
              DIIN(IPLS,:)=MAX(DVAC,DIINTF(IOLD,:))

              IF (NLMLV) THEN
                VXIN(IPLSV,:)=VXINTF(IOLDV,:)
                VYIN(IPLSV,:)=VYINTF(IOLDV,:)
                VZIN(IPLSV,:)=VZINTF(IOLDV,:)
CDR           ELSE
CDR  in single V.IN cases: all V.IN(IPLS,:) are the same, and == V.IN(1,:)
              END IF

            else
              write (iunout,*) 'Warning from PLASMA_DERIV'
              write (iunout,*) '"Density-model '//FORT//'13"'
              write (iunout,*) 'Species IPLS= ',iold, ' is either == 0,'
              write (iunout,*)
     .             'or could not be read from '//fort_lc//'13'
              write (iunout,*) 'Species IPLS left unchanged'
            ENDIF

          CASE ('FORT.10','FTN10')

c   itold = ?, type of particle on fort.10?
c   check: itold ge 0 and itold le 3
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)    ! species
cdr         itold=?                         ! type
            IOLDTI=MPLSTI(IOLD)             ! ?
            IOLDV=MPLSV(IOLD)

            ALLOCATE (BASE_DENSITY(NRAD))
            ALLOCATE (BASE_TEMP(NRAD))
            ALLOCATE (BASE_VELOC(3,NRAD))

            CALL EIRENE_GET_BASE_PARAMETERS(1)

            DIIN(IPLS,1:NSBOX)=MAX(DVAC,BASE_DENSITY(1:NSBOX))

            IF (NLMLTI) THEN
              TIIN(IPLSTI,:)=MAX(TVAC,BASE_TEMP(:))
cdr ???     ELSE
            ENDIF
            IF (NLMLV) THEN
              VXIN(IPLSV,:)=VXIN(IOLDV,:)
              VYIN(IPLSV,:)=VYIN(IOLDV,:)
              VZIN(IPLSV,:)=VZIN(IOLDV,:)
cdr ???     ELSE

            END IF

            DEALLOCATE (BASE_DENSITY)
            DEALLOCATE (BASE_TEMP)
            DEALLOCATE (BASE_VELOC)

            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              FOUND = .TRUE.
              DO IR=1,NSBOX
                IF (LSPCCLL(IR)) THEN
                  CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                  IF (FOUND) THEN
                    IBS = IBS + 1
                    SPEC%IPRTYP = 4
                    SPEC%IPRSP = IPLS
                    BACK_SPEC(IBS) = SPEC
                  END IF
                END IF
              ENDDO
            END IF

          CASE ('CONSTANT')

            IF (NLMLTI) THEN
              TIIN(IPLSTI,:)=MAX(TVAC,TDMPAR(IPLS)%TDM%TVAL)
            ELSE
CDR  ??
            ENDIF
            DIIN(IPLS,:)=MAX(DVAC,TDMPAR(IPLS)%TDM%DVAL)
            IF (NLMLV) THEN
              VXIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VXVAL
              VYIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VYVAL
              VZIN(IPLSV,:)=TDMPAR(IPLS)%TDM%VZVAL
            ELSE
CDR  ??
            END IF

          CASE ('MULTIPLY')
c           ITOLD=TDMPAR(IPLS)%TDM%ITP(1) =4,  hard-wired
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)

            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)

            ALLOCATE (BASE_DENSITY(NRAD))
            ALLOCATE (BASE_TEMP(NRAD))
            ALLOCATE (BASE_VELOC(3,NRAD))

            CALL EIRENE_GET_BASE_PARAMETERS(1)

            DO IR=1,NSBOX
              IF (NLMLTI) TIIN(IPLSTI,IR)=MAX(TVAC,BASE_TEMP(IR)*
     .                                        TDMPAR(IPLS)%TDM%TFACTOR)
              DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR)*
     .                               TDMPAR(IPLS)%TDM%DFACTOR)
              IF (NLMLV) THEN
                VXIN(IPLSV,IR)=VXIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
                VYIN(IPLSV,IR)=VYIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
                VZIN(IPLSV,IR)=VZIN(IOLDV,IR)*TDMPAR(IPLS)%TDM%VFACTOR
              END IF
            ENDDO

            DEALLOCATE (BASE_DENSITY)
            DEALLOCATE (BASE_TEMP)
            DEALLOCATE (BASE_VELOC)

          CASE ('BOLTZMANN')
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)

            ALLOCATE (BASE_DENSITY(NRAD))
            ALLOCATE (BASE_TEMP(NRAD))
            ALLOCATE (BASE_VELOC(3,NRAD))

            CALL EIRENE_GET_BASE_PARAMETERS(1)

            IF (NLMLTI) THEN
              TIIN(IPLSTI,:)=BASE_TEMP(:)
            ELSE
CDR  ??
            ENDIF
            IF (NLMLV) THEN
              VXIN(IPLSV,:)=VXIN(IOLDV,:)
              VYIN(IPLSV,:)=VYIN(IOLDV,:)
              VZIN(IPLSV,:)=VZIN(IOLDV,:)
            ELSE
CDR  ??
            END IF

            FOUND = .TRUE.
            G_BOLTZ=TDMPAR(IPLS)%TDM%G_BOLTZ
            DELTAE=TDMPAR(IPLS)%TDM%DELTAE
            DO IR=1,NSBOX
              BOLTZFAC=0._DP
              IF (TIIN(IOLDTI,IR).GT.TVAC)
     .          BOLTZFAC=G_BOLTZ*EXP(-DELTAE/BASE_TEMP(IR))
              DIIN(IPLS,IR)=MAX(DVAC,BASE_DENSITY(IR)*BOLTZFAC)

              IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                IF (LSPCCLL(IR)) THEN
                  CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                  IF (FOUND) THEN
                    IBS = IBS + 1
                    SPEC%SPC = SPEC%SPC * BOLTZFAC
                    SPEC%IPRTYP = 4
                    SPEC%IPRSP = IPLS
                    BACK_SPEC(IBS) = SPEC
                  END IF
                END IF
              END IF

            ENDDO
            DEALLOCATE (BASE_DENSITY)
            DEALLOCATE (BASE_TEMP)
            DEALLOCATE (BASE_VELOC)

          CASE ('PLANCK')
cdr  use planck function, grey body: const. emissivity factor
cdr  and set "grey" photon density here
            IOLD=TDMPAR(IPLS)%TDM%ISP(1)
            IOLDTI=MPLSTI(IOLD)

            ALLOCATE (BASE_DENSITY(NRAD))
            ALLOCATE (BASE_TEMP(NRAD))
            ALLOCATE (BASE_VELOC(3,NRAD))

            CALL EIRENE_GET_BASE_PARAMETERS(1)

            IF (NLMLTI) TIIN(IPLSTI,:)=BASE_TEMP(:)

            FOUND = .TRUE.
            G_PLANCK=TDMPAR(IPLS)%TDM%G_PLANCK
            DO IR=1,NSBOX
            ENDDO

            DEALLOCATE (BASE_DENSITY)
            DEALLOCATE (BASE_TEMP)
            DEALLOCATE (BASE_VELOC)
c         CASE DEFAULT
c
        END SELECT
      END DO  ! NPLSI

      IF (ALLOCATED(DEINTF)) DEALLOCATE(DEINTF)

c......................................................................

C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
C  SET ELECTRON DENSITY FROM QUASI-NEUTRALITY, FURTHER: TEINL, DEINL, LGVAC(..,0:NPLS+1)
      LGVAC=.TRUE.
      !DEIN=0._DP already initialised in plasma.f. Don't recompute if filled in prousr

      DO J=1,NSBOX
        IF(DEIN(J)==0._DP) THEN
        DO JPLS=1,NPLSI
cnh       28.10.2019
          IF(ZIIN(JPLS,J).NE.ZVAC) THEN
            DEIN(J)=DEIN(J)+ZIIN(JPLS,J)*DIIN(JPLS,J)
          ELSE
            DEIN(J)=DEIN(J)+DBLE(NCHRGP(JPLS))*DIIN(JPLS,J)
          ENDIF
        END DO
        ENDIF
C  SET 'LOG OF TEMPERATURE AND DENSITY' ARRAYS
        ZTEI=MAX(TVAC,MIN(TEIN(J),1.E10_DP))
        TEINL(J)=LOG(ZTEI)
        ZTNE=MAX(DVAC,MIN(DEIN(J),1.E20_DP))
        DEINL(J)=LOG(ZTNE)

cdr  set "vacuum flags": lgvac(j,ipls) turns off all reactions with
cdr  background species ipls, in cell j
cdr  ipls=npls+1:  electrons
cdr  ipls=0     :  all species, i.e.: all reactions are turned off in this cell.
cdr  ipls       :  set below, after special "density models" are done.
        TEPLS=TEIN(J)
        DEPLS=DEIN(J)
        LGVAC(J,NPLS+1)=TEPLS.LE.TVAC.OR.DEPLS.LE.DVAC
        LGVAC(J,0)     =LGVAC(J,0).AND.LGVAC(J,NPLS+1)
      END DO

      tpb2 = EIRENE_second_own()
      IF (TRCTIM) write (iunout,*) ' CPU time for log values ',tpb2-tpb1
      tpb1 = tpb2
c.....................................................................
C
C   FURTHER SPECIAL DENSITY MODELS, AFTER ELECTRON DENSITY DEIN IS SET:
C   SAHA  (NOT READY)
C   CORONA
C   COLRAD
C
      ALLOCATE (BASE_DENSITY(NRAD))
      ALLOCATE (BASE_TEMP(NRAD))
      ALLOCATE (BASE_VELOC(3,NRAD))

      DO JPLS=1,NPLSI
        IPLS=JPLS

        IF (LEN_TRIM(CDENMODEL(IPLS)) == 0) CYCLE

        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)

        SELECT CASE (CDENMODEL(IPLS))

        CASE ('SAHA      ')
!PB   TO BE WRITTEN
          WRITE (iunout,*) ' DENSITY PROFILE ACCORDING TO SAHA IS NOT ',
     .                'AVAILABLE '
          WRITE (iunout,*) ' CHOOSE DIFFERENT OPTION FOR ION ',
     .                'DENSITY ',IPLS
          CALL EIRENE_EXIT_OWN(1)
c...............................................................saha: done
        CASE ('CORONA    ')
cdr  density of a background "isotope" IPLS is derived from balance between
cdr  a single step electron impact excitation or ionisation
cdr  from a donor state IOLD=TDMPAR(IPLS)%TDM%ISP(1) ("gain"= RCORONA*DEIN)
cdr  followed by a radiative decay A_CORONA of that "isotope" IPLS ("loss" = A_CORONA).
cdr  The balance: n_iold * gain  = n_ipls * loss  provides density n_ipls
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)

          CALL EIRENE_GET_BASE_PARAMETERS(1)

          IF (NLMLTI) TIIN(IPLSTI,:)=BASE_TEMP(:)
cdr  why not: =TIIN(IOLDTI,:) ?
c
          IF (NLMLV) THEN
            VXIN(IPLSV,:)=VXIN(IOLDV,:)
            VYIN(IPLSV,:)=VYIN(IOLDV,:)
            VZIN(IPLSV,:)=VZIN(IOLDV,:)
          END IF

cdr inverse of A_CORONA
cdr tbd:  avoid division by zero here
          AM1=1._DP/TDMPAR(IPLS)%TDM%A_CORONA
          IRC = TDMPAR(IPLS)%TDM%IRC(1)

          DO IR=1,NSBOX
cdr nsbox rather than nsurf? Otherwise this will not work in additional cells
            RCORONA=0.0
            IF (.NOT.LGVAC(IR,NPLS+1)) THEN
cdr  Te cut off safety at 0.1 eV
              TEF=max(-2.3_DP,TEINL(IR))
              RCORONA = EIRENE_RATE_COEFF(IRC,IR,TEF,0._DP,.TRUE.,0)
            END IF
c  now RCORONA contains the excitation rate coefficient (cm**3/s),
c  and AMI is the inverse of the radiative decay rate (s).
c  Compute new density of species IPLS from equilibrium between
c  these two processes for the given "ground state" density BASE_DENSITY

            DIIN(IPLS,IR)=BASE_DENSITY(IR)*RCORONA*DEIN(IR)*AM1
            DIIN(IPLS,IR)=MAX(DVAC,DIIN(IPLS,IR))

cdr  what is this?  background spectrum ?
cdr  if so, why in corona part?
            IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
              IF (LSPCCLL(IR)) THEN
                CALL EIRENE_GET_SPECTRUM (IR,1,SPEC,FOUND)
                IF (FOUND) THEN
                  IBS = IBS + 1
                  SPEC%IPRTYP = 4
                  SPEC%IPRSP = IPLS
                  SPEC%SPC = SPEC%SPC * RCORONA*DEIN(IR)*AM1
                  BACK_SPEC(IBS) = SPEC
                END IF
              END IF
            END IF

          END DO
c...............................................................corona: done

          tpb2 = EIRENE_second_own()
          IF (TRCTIM)
     .     write (iunout,*) ' CPU time for corona ',ipls,tpb2-tpb1
          tpb1 = tpb2

        CASE ('COLRAD    ')
cdr  density of a background "isotope" IPLS is derived from collision radiative
cdr  models, as an CR equilibrium population of excited states.
cdr  From one or several donor states (components/contributions)
cdr  TEIN and DEIN (electron parameters) are given already.

cdr  IOLD=TDMPAR(IPLS)%TDM%ISP(IRE)
          IF (.NOT.ALLOCATED(SUMNI)) THEN
            ALLOCATE (SUMNI(NRAD))
            ALLOCATE (SUMMNI(NRAD))
          END IF
          SUMNI = 0._DP
          SUMMNI = 0._DP
          DIIN(IPLS,:)=0._DP
C  ARE THERE MULTIPLE ION TEMPERATURES?
          if (nlmlti) then
            TIIN(IPLSTI,:)=0._DP
cdr       else ??
cdr  one common temperature for all IPLS. Nothing to be done here.
          endif
C  ARE THERE MULTIPLE ION DRIFT VELOCITIES?
          if (nlmlv) then
            VXIN(IPLSV,:)=0._DP
            VYIN(IPLSV,:)=0._DP
            VZIN(IPLSV,:)=0._DP
cdr       else ??
cdr  one common flow velocity field for all IPLS. Nothing to be done here.
          endif

          FOUND = .TRUE.
          DO IRE=1,TDMPAR(IPLS)%TDM%NRE
c  sum up contributions coupled to one or more (NRE) base-species
            IOLD=TDMPAR(IPLS)%TDM%ISP(IRE)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
c           write (iunout,*) 'test colrad ',iold,ioldti,ioldv
            CALL EIRENE_GET_BASE_PARAMETERS(IRE)

c  temperature and density dependence in reduced population coefficient
cdr  Jan 2019: generalized from SELECT CASE (ISW): AMJUEL H.11, H.12
cdr            now to more general reaction card IRC
cdr            (e.g. population coefficients from internal CR codes
cdr                  or from AMJUEL H.11, or H.12 data fits)
            IRC = TDMPAR(IPLS)%TDM%IRC(IRE)
            DO IR=1,NSBOX
cdr nsbox rather than nsurf? Otherwise this will not work in additional cells
              IF (LGVAC(IR,NPLS+1)) CYCLE

              TE=TEIN(IR)
              DE=DEIN(IR)
              DEF=LOG(DE)
              TEF=max(-2.3_DP,LOG(TE)) ! cut off at 0.1 eV

              RCOLRAD=EIRENE_OTHER_RATE_COEFF(IRC,IR,TEF,DEF,.TRUE.,1)
              DIIN(IPLS,IR)=DIIN(IPLS,IR)+BASE_DENSITY(IR)*RCOLRAD
              if (nlmlti) then
                TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                            BASE_DENSITY(IR)*BASE_TEMP(IR)
              endif
              if (nlmlv) then
                VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VXIN(IOLDV,IR)
                VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VYIN(IOLDV,IR)
                VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                 NMASSP(IOLD)*BASE_DENSITY(IR)*VZIN(IOLDV,IR)
              endif
              SUMNI(IR) = SUMNI(IR) + BASE_DENSITY(IR)
              SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*BASE_DENSITY(IR)

cdr  oct 18: this must be wrong. What have spectra to do with CR models ?
cdr          Perhaps a left over from the obsolete corresponding option
cdr          in diagno, i.e. attempts to prepare line of sight scored
cdr          spectrally resolved densities, e.g. for photonic options ?

cdr identical code as above for H.11 ?
              IF ((ICALL > 0) .AND. (NBACK_SPEC > 0)) THEN
                IF (LSPCCLL(IR)) THEN
                  CALL EIRENE_GET_SPECTRUM (IR,IRE,SPEC,FOUND)
                  IF (FOUND) THEN
                    IF (IRE == 1) THEN
                      IBS = IBS + 1
                      SPEC%IPRTYP = 4
                      SPEC%IPRSP = IPLS
                      SPEC%SPC = SPEC%SPC * RCOLRAD
                      BACK_SPEC(IBS) = SPEC
                    ELSE
                      IF (  (BACK_SPEC(IBS)%NSPC == SPEC%NSPC)
     .                 .AND.(BACK_SPEC(IBS)%SPCMIN==SPEC%SPCMIN)
     .                 .AND.(BACK_SPEC(IBS)%SPCMAX==SPEC%SPCMAX))
     .                THEN
                        SPEC%SPC = SPEC%SPC * RCOLRAD
                        BACK_SPEC(IBS)%SPC =
     .                       BACK_SPEC(IBS)%SPC + SPEC%SPC
                      ELSE
                        WRITE (IUNOUT,*) ' ERROR IN PLASMA_DERIV,',
     .                       ' DENSITY MODEL COLRAD '
                        WRITE (IUNOUT,*) ' TWO SPECTRA CONTRIBUTING ',
     .                       ' TO THE SAME BACKGROUND SPECTRUM DO NOT',
     .                       ' MATCH '
                        WRITE (IUNOUT,*) ' IPLS, IRE ',IPLS, IRE
                        CALL EIRENE_EXIT_OWN(1)
                      END IF
                    END IF
                  END IF
                END IF
              END IF

            END DO  ! IR
          END DO  ! IRE

c  scale merged contributions from all contributing base-species
          DIIN(IPLS,:)=MAX(DVAC,DIIN(IPLS,:))
          if (nlmlti) then
            TIIN(IPLSTI,:)=MAX(TVAC,TIIN(IPLSTI,:)/(SUMNI(:)+eps60))
          endif
          if (nlmlv) then
            VXIN(IPLSV,:)=VXIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VYIN(IPLSV,:)=VYIN(IPLSV,:)/(SUMMNI(:)+eps60)
            VZIN(IPLSV,:)=VZIN(IPLSV,:)/(SUMMNI(:)+eps60)
          endif
c .................................................................colrad done

        CASE DEFAULT
!  NOTHING TO BE DONE HERE, ALREADY COMPLETED
        END SELECT ! density model

        tpb2 = EIRENE_second_own()
        IF (TRCTIM)
     .   write (iunout,*) ' CPU time for colrad ',ipls,tpb2-tpb1
        tpb1 = tpb2

      END DO   ! ipls

      IF (ALLOCATED(SUMNI)) THEN
        DEALLOCATE (SUMNI)
        DEALLOCATE (SUMMNI)
      END IF
      DEALLOCATE (BASE_DENSITY)
      DEALLOCATE (BASE_TEMP)
      DEALLOCATE (BASE_VELOC)

      NBACK_SPEC = IBS

      tpb2 = EIRENE_second_own()
      IF (TRCTIM)
     . write (iunout,*) ' CPU time for density models ',tpb2-tpb1
      tpb1 = tpb2
C
C  SPECIAL PLASMA BACKGROUND MODELS DONE
C
C  NEXT: SET SOME "DERIVED" FIELDS: EDRIFT, BPERP, BV_VEC, PMOM_VEC,
C                                   LGVAC, TIINL, DIINL, ZT1, ZRG

C  SET DRIFT ENERGY (EV)
      IF (LEDRIFT .AND. NLDRFT) THEN
        DO J=1,NSBOX
          DO IPLS=1,NPLSI
            IPLSV=MPLSV(IPLS)
            IF (INDPRO(4) == 8) THEN
              IF(IPLS.EQ.1) THEN
                EDRIFT(IPLS,J)=CVRSSP(IPLS)*EIRENE_VDION(J)**2
              ELSE
C               WRITE(iunout,*)'WARNING PLASMA_DERIV: IPLS>1 NO DRIFT!'
                EDRIFT(IPLS,J)=0._DP
              END IF
            ELSE
              EDRIFT(IPLS,J)=CVRSSP(IPLS)*
     .              (VXIN(IPLSV,J)**2+VYIN(IPLSV,J)**2+VZIN(IPLSV,J)**2)
            END IF
          END DO
        END DO
      ELSEIF (LEDRIFT .AND. .NOT. NLDRFT) THEN
        EDRIFT(1:NPLSI,1:NSBOX)=0.D0
      END IF
C
C  SET A DEFAULT B_PERP UNIT VECTOR in POL PLANE (X,Y) FOR 2D cases,
C      i.e. the 2D vector grad(PSI(X,Y)) in 2D cases.
C      B_PAR IS ALREADY GIVEN AS INPUT TALLY BXIN,BYIN,BZIN
C
c  IF BXIN AND BYIN ARE NOT AVAILABLE: ALSO: LBXPERP=LBYPERP=.FALSE.

cdr note: by default: lbxperp = lbyperp = false.
cdr       These "derived tallies" must be turned on explicitly,
cdr       at the end of input block 5, "optional" tallies.
cdr       Otherwise the perp and
cdr       diamagn. components of BV and PMOM are zero.
cdr probably these default settings need to be revised?
cdr June 23: Done: LBXPERP = LBYPERP = TRUE by default now.
cdr          SETPRM_INTAL modified accordingly for tallies -16 and -17

      IF (LBXPERP .AND. LBYPERP) THEN

       ico1=0
       ico2=0
       DO J=1,NSBOX
cdr  lets see if these are already set:
        bxpp=bxperp(j)
        bypp=byperp(j)
        bpp=sqrt(bxperp(j)**2+byperp(j)**2)
        if (bpp .gt. eps5) then
c         write (iunout,*) 'b_perp already available ,j ',bpp,j
          ico1=ico1+1
        endif

        IF (ABS(BXIN(J)) > EPS10) THEN
           BYP = 1._DP
           BXP = -BYIN(J)/BXIN(J)
        ELSEIF (ABS(BYIN(J)) > EPS10) THEN
           BXP = 1._DP
           BYP = -BXIN(J)/BYIN(J)
        ELSE
           BXP = 0._DP
           BYP = 0._DP
        END IF

C  CHECK ORIENTATION, SET B_PERP SUCH THAT B_PERP CROSS B_PAR > 0
C  (this coincides with the z-component of B_DIA being > 0)
        IF (BXIN(J)*BYP-BXP*BYIN(J) < 0._DP) THEN
           BXP = -BXP
           BYP = -BYP
        END IF
C  NORMALIZE
        BNORM=SQRT(BXP*BXP+BYP*BYP)
        if (bnorm.gt.eps60) then
          BXP=BXP/BNORM
          BYP=BYP/BNORM
        endif
cdr   and... bzperp(j)=0, by construction in 2D symmetric cases

        if (bpp .ne. 0.0) then
cdr compare with default
          dist1=(bxpp-bxp)**2 + (bypp-byp)**2
cdr  try opposite sign
          dist2=(bxpp+bxp)**2 + (bypp+byp)**2
          if (dist1 .gt. eps5 .and. dist2 .gt. eps5) then
            write (iunout,*) 'j available but different',j,dist1,dist2
          elseif (dist1 .gt. eps5) then
c           write (iunout,*) 'j available but different sign',j
            ico2=ico2+1
          endif
cdr  use the B_PERP transferred from outside, rather than the default
          bxperp(j)=bxpp
          byperp(j)=bypp
        else
cdr  use the default setting
          bxperp(j)=bxp
          byperp(j)=byp
        endif
       END DO  ! NSBOX
       write (iunout,*) 'B_PERP available in ICO1 cells ',ico1,nsbox
       write (iunout,*) 'B_PERP sign changed in ICO2 cells ',ico2
      END IF

      DO 5103 J=1,NSBOX
C  SET 'VACUUM REGION FLAGS'
C  LGVAC(...,0)  VACUUM, NO REACTION RATES AT ALL
C  LGVAC(...,IPLS)       NO REACTION RATES FOR BACKGROUND SPECIES IPLS
C  LGVAC(...,NPLS+1)     NO REACTION RATES FOR BACKGROUND ELECTRONS
C                        BUT PERHAPS FOR NEUTRAL BACKGROUND
        DO 5106 IPLS=1,NPLSI
          IPLSTI=MPLSTI(IPLS)
          EMPLS=1.5*TIIN(IPLSTI,J)
          IF (LEDRIFT) EMPLS=EMPLS+EDRIFT(IPLS,J)
          DIPLS=DIIN(IPLS,J)
          LGVAC(J,IPLS)=EMPLS.LE.TVAC.OR.DIPLS.LE.DVAC
          LGVAC(J,0)   =LGVAC(J,0).AND.LGVAC(J,IPLS)
 5106   CONTINUE
 5103 CONTINUE

      IF (LEVGEO.EQ.3) THEN
cdr set vacuum flags in polygonal grid cut cells (if any)
        DO 5161 I=1,NPPLG-1
          DO 5162 IP=NPOINT(2,I),NPOINT(1,I+1)-1
            IPM=IP-1
            DO 5163 IPLS=0,NPLS+1
              DO IR=1,NR1STM
                IN=IR+IPM*NR1ST
                LGVAC(IN,IPLS)=.TRUE.
              END DO
 5163       CONTINUE
 5162     CONTINUE
 5161   CONTINUE
      ENDIF
C
      DO 5205 IPLS=1,NPLSI
C  FACTOR FOR MOST PROBABLE SPEED, ALSO: THERMAL SPEED
        FCT0=1./RMASSP(IPLS)*2.*CVEL2A*CVEL2A
C  FACTOR FOR MEAN SPEED
        FCT1=1./RMASSP(IPLS)*8./PIA*CVEL2A*CVEL2A
C  FACTOR FOR ROOT MEAN SQUARE SPEED
        FCT2=1./RMASSP(IPLS)*3.*CVEL2A*CVEL2A

        FCRG=CVEL2A/SQRT(RMASSP(IPLS))
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)

        IF (LBV_VEC) THEN
cdr  tally -23, 1:3*nplsv
          DO ICO=0,2
            KK=ICO*NPLSV
            BV_VEC(KK+IPLSV,:)=0._DP
          ENDDO
        ENDIF
        IF (LPMOM_VEC) THEN
cdr   tally -24, 1:3*npls
          DO ICO=0,2
            KK=ICO*NPLS
            PMOM_VEC(KK+IPLS,:)=0._DP
          ENDDO
        ENDIF

        DO J=1,NSBOX
          ZTII=MAX(TVAC,MIN(TIIN(IPLSTI,J),1.E10_DP))
          TIINL(IPLSTI,J)=LOG(ZTII)
CDR  BPAR: PARALLEL COMPONENT UNIT VECTOR, DEFAULT: Z COORDINATE ONLY
          bpar(1)=0._dp
          bpar(2)=0._dp
          bpar(3)=1._dp
          if (lbxin) bpar(1)=bxin(j)
          if (lbyin) bpar(2)=byin(j)
          if (lbzin) bpar(3)=bzin(j)
          IF (LBV_VEC)
     .      BV_VEC(IPLSV,J)=Bpar(1)*VXIN(IPLSV,J)+
     .                      Bpar(2)*VYIN(IPLSV,J)+
     .                      Bpar(3)*VZIN(IPLSV,J)
          IF (LPMOM_VEC.AND.LBV_VEC)
     .      PMOM_VEC(IPLS,J)=BV_VEC(IPLSV,J)*
     .                       AMUA*RMASSP(IPLS)

CDR  BPER: PERP COMPONENT, GRAD PSI, DEFAULT: 0
cdr  PSI function input tally is not fully ready, I believe.
cdr  So I use input tallies bxperp,byperp.  Needs to be checked.
          bper(1)=0._dp
          bper(2)=0._dp
          bper(3)=0._dp
          if (lbxperp) bper(1)=bxperp(j)
          if (lbyperp) bper(2)=byperp(j)
          IF (LBV_VEC)
     .      BV_VEC(NPLSV+IPLSV,J)=Bper(1)*VXIN(IPLSV,J)+
     .                            Bper(2)*VYIN(IPLSV,J)+
     .                            Bper(3)*VZIN(IPLSV,J)
          IF (LPMOM_VEC.AND.LBV_VEC)
     .      PMOM_VEC(NPLS+IPLS,J)=BV_VEC(NPLSV+IPLSV,J)*
     .                            AMUA*RMASSP(IPLS)

CDR  DIAMAGN. COMPONENT (B X GRAD PSI), DEFAULT: 0
          bcross(1)=0._dp
          bcross(2)=0._dp
          bcross(3)=0._dp
          bcross(1) = bpar(2)*bper(3) - bpar(3)*bper(2)
          bcross(2) =-bpar(3)*bper(1) + bpar(1)*bper(3)
          bcross(3) = bpar(1)*bper(2) - bpar(2)*bper(1)
          IF (LBV_VEC)
     .      BV_VEC(2*NPLSV+IPLSV,J)=BCROSS(1)*VXIN(IPLSV,J)+
     .                              BCROSS(2)*VYIN(IPLSV,J)+
     .                              BCROSS(3)*VZIN(IPLSV,J)
          IF (LPMOM_VEC.AND.LBV_VEC)
     .      PMOM_VEC(2*NPLS+IPLS,J)=BV_VEC(2*NPLSV+IPLSV,J)*
     .                              AMUA*RMASSP(IPLS)

C
C  ZT1: FOR "EFFECTIVE" PLASMA PARTICLE VELOCITY IN CROSS-SECTIONS
C       FOR HEAVY PARTICLE INTERACTIONS
C       SQRT(ZT1) IS THE MEAN VELOCITY V_M AT TI=ZTII, TAKEN AS
C       ROOT MEAN SQUARE SPEED, M/2 V_M^2 = 3/2 KT
C
          ZT1(IPLS,J)=FCT2*ZTII
C
C  ZRGQ: VARIANCE FOR SAMPLING FROM MAXWELLIAN VELOCITY DISTRIBUTION
C  ZRG=SQRT(ZRGQ) = STANDARD DEVIATION
C
          ZRG(IPLS,J)=FCRG*SQRT(ZTII)
C
          ZTNI=MAX(DVAC,MIN(DIIN(IPLS,J),1.E20_DP))
          DIINL(IPLS,J)=LOG(ZTNI)
        END DO
 5205 CONTINUE
C
!     INTERPOLATE PLASMA PROFILES TO CELL VERTICES

      CALL EIRENE_ALLOC_CORNERS
      IF (LTESMO) THEN
        call eirene_cell_to_corner(TEIN,TEINCORNER)
      END IF

      IF (LTISMO) THEN
        do iplsti = 1, nplsti
          call eirene_cell_to_corner(TIIN(iplsti,:),
     .                               TIINCORNER(:,iplsti))
        end do
      END IF

      IF (LDESMO) THEN
        call eirene_cell_to_corner(DEIN,DEINCORNER)
      END IF

      IF (LDISMO) THEN
        do ipls = 1, npls
          call eirene_cell_to_corner(DIIN(ipls,:),DIINCORNER(:,ipls))
        end do
      ENDIF

      IF (LEDRIFTSMO) THEN
        do ipls = 1, npls
         call eirene_cell_to_corner(EDRIFT(ipls,:),EDRIFTCORNER(:,ipls))
        ENDDO
      ENDIF
      IF (LPMOM_VECSMO) THEN
        do ipls = 1, npls
         call eirene_cell_to_corner(PMOM_VEC(ipls,:),
     .          PMOM_VECCORNER(:,ipls))
         call eirene_cell_to_corner(PMOM_VEC(npls+ipls,:),
     .          PMOM_VECCORNER(:,npls+ipls))
         call eirene_cell_to_corner(PMOM_VEC(2*npls+ipls,:),
     .          PMOM_VECCORNER(:,2*npls+ipls))
        end do
      ENDIF

      IF (LVSMO) THEN
        do iplsv = 1, nplsv
          if (lvxsmo)
     .     call eirene_cell_to_corner(VXIN(iplsv,:),VXINCORNER(:,iplsv))
          if (lvysmo)
     .     call eirene_cell_to_corner(VYIN(iplsv,:),VYINCORNER(:,iplsv))
          if (lvzsmo)
     .     call eirene_cell_to_corner(VZIN(iplsv,:),VZINCORNER(:,iplsv))
          if (lbv_vecsmo) THEN
            call eirene_cell_to_corner(BV_VEC(iplsv,:),
     .          BV_VECCORNER(:,iplsv))
            call eirene_cell_to_corner(BV_VEC(nplsv+iplsv,:),
     .          BV_VECCORNER(:,nplsv+iplsv))
            call eirene_cell_to_corner(BV_VEC(2*nplsv+iplsv,:),
     .          BV_VECCORNER(:,2*nplsv+iplsv))
          ENDIF
        end do
      END IF

      IF (LBSMO) THEN
C  SMOOTH B FIELD
        if (lbxsmo) call eirene_cell_to_corner(BXIN,BXINCORNER)
        if (lbysmo) call eirene_cell_to_corner(BYIN,BYINCORNER)
        if (lbzsmo) call eirene_cell_to_corner(BZIN,BZINCORNER)
        if (lbfsmo) call eirene_cell_to_corner(BFIN,BFINCORNER)
      END IF

      IF (LESMO) THEN
C  SMOOTH E FIELD
        if (lexsmo)  call eirene_cell_to_corner(EXIN,EXCORNER)
        if (leysmo)  call eirene_cell_to_corner(EYIN,EYCORNER)
        if (lezsmo)  call eirene_cell_to_corner(EZIN,EZCORNER)
        if (lefsmo)  call eirene_cell_to_corner(EFIN,EFCORNER)
        if (lpotsmo) call eirene_cell_to_corner(POT,POTCORNER)
      END IF

      IF (LADSMO) THEN
        do iain = 1, nain
          call eirene_cell_to_corner(ADIN(iain,:),ADCORNER(:,iain))
        end do
      END IF

      IF (LVOLSMO) THEN
        call eirene_cell_to_corner(VOL,VOLCORNER)
      END IF

      IF (LWGHTSMO) THEN
        do ispz = 1, nspzmc
          call eirene_cell_to_corner(WGHT(ispz,:),WGHTCORNER(:,ispz))
        end do
      END IF

      IF (LBXPSMO) THEN
        call eirene_cell_to_corner(BXPERP,BXPERPCORNER)
      END IF

      IF (LBYPSMO) THEN
        call eirene_cell_to_corner(BYPERP,BYPERPCORNER)
      END IF

      IF (LPOTSMO) THEN
        call eirene_cell_to_corner(POT,POTCORNER)
      END IF

      IF (LPSISMO) THEN
        call eirene_cell_to_corner(PSI,PSICORNER)
      END IF

      IF (LZISMO) THEN
        do ipls = 1, npls
          call eirene_cell_to_corner(ZIIN(ipls,:),ZIINCORNER(:,ipls))
        ENDDO
      END IF

      IF (LFREE27SMO) THEN
        call eirene_cell_to_corner(FREE27,FREE27CORNER)
      END IF

      IF (LFREE28SMO) THEN
        call eirene_cell_to_corner(FREE28,FREE28CORNER)
      END IF

      IF (LFREE29SMO) THEN
        call eirene_cell_to_corner(FREE29,FREE29CORNER)
      END IF

      IF (LFREE30SMO) THEN
        call eirene_cell_to_corner(FREE30,FREE30CORNER)
      END IF
C
C  GRADIENTS
cdr       1:NTALG regular input tallies
cdr NTALG+1:NTALI gradients of these regular input tallies
      IF (ANY(LIVTALI(NTALG+1:NTALI))) THEN
        ALLOCATE (TALLY(NCORNER))
        DO ITAL = 1, NTALG
C  ITAL : NO. OF TALLY OF WHICH GRADIENT IS TO BE CALCULATED
C  KK   : INDEX OF GRADIENT TALLY COMPONENTS: KK+1, KK+2, KK+3
          KK = NTALG + (ITAL-1)*3
          IF (ANY(LIVTALI(KK+1:KK+3))) THEN
            NFTI = 1
            NFTE=NFSTPI(KK+1)  ! same as for KK+2, KK+3
            DO K=NFTI, NFTE
              TALLY(:)=CORNER_PROFILES(:,NADDCOR(ITAL)+K)
              CALL EIRENE_CALC_GRAD(TALLY,
     .                              PLSTLS(NADDP(KK+1)+K,:),
     .                              PLSTLS(NADDP(KK+2)+K,:),
     .                              PLSTLS(NADDP(KK+3)+K,:),
     .                              LIVTALI(KK+1),
     .                              LIVTALI(KK+2),
     .                              LIVTALI(KK+3))
            END DO
          END IF
        END DO
        DEALLOCATE (TALLY)
      END IF
C
C  SAVE PLASMA DATA AND ATOMIC DATA ON FORT.13
C
      IF ((NFILEL ==1) .OR. (NFILEL ==3) .OR. (NFILEL ==4)) THEN
cdr   calls either wrplam_long or wrplam_shrt, see flag: NLSHRT13
cdr     CALL EIRENE_WRPLAM(TRCFLE,'PLASMA_DERIV')
      END IF

      RETURN

      CONTAINS

cdr Set arrays Density, Temperature, from previous run (fort.10) or
cdr at the end of the present run, for background no. IPLS
cdr So far: temperature is defined as ratio-tally from energy density
cdr         and particle density, ignoring flow velocities.
cdr Also called with ICALL=0, then uses background field tallies,
cdr e.g., for COLRAD model with coupling to continuum.

      SUBROUTINE EIRENE_GET_BASE_PARAMETERS(IRE)
c  input: ire: number of density that contributes to the
c              evaluation of the expression for the selected species
c              ipls with special density/temperature option
c         ipls:
c         icall:
c  output: base_density, base_temp (missing: base_veloc)
c
      INTEGER, INTENT(IN) :: IRE
      INTEGER :: IG, IT, ISTRA, ITYP
      EXTERNAL :: EIRENE_RSTRT, EIRENE_SYMET, EIRENE_EXIT_OWN

      BASE_DENSITY = 0._DP
      BASE_TEMP = 0._DP
      BASE_VELOC = 0._DP

cdr fix ph1: moved here, because needed also for icall=0
      ISTRA = TDMPAR(IPLS)%TDM%ISTR(IRE)
      ITYP  = TDMPAR(IPLS)%TDM%ITP(IRE)
c     write (iunout,*) 'base_parameters',istra,ityp,ire

      IF (ICALL > 0) THEN

C FOR CALLS AFTER PARTICLE TRACING:
cdr: only needed if ityp .lt. 4 (test particles, not bulk particles)
cdr  fetch MC output tallies.
cdr  Particle densities, energy densities and momentum densities
cdr  would suffice. No need to read entire fort.10, fort.11
cdr  for this stratum.

        IF (ISTRA.EQ.IESTR.OR.ITYP.EQ.4) THEN
C  NOTHING TO BE DONE
        ELSEIF ((NFILEN.EQ.1.OR.NFILEN.EQ.2).OR.
     .         ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0)) THEN
          IESTR=ISTRA
          IF (TRCFLE) WRITE (IUNOUT,*) 'FROM PLASMA_DERIV:'
          CALL EIRENE_RSTRT(ISTRA,NSTRAI,
     .               NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .               TRCFLE)
          IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
            CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                 NLSYMP(ISTRA),NLSYMT(ISTRA))
          ENDIF
        ELSE
          WRITE (iunout,*)
     .      'ERROR IN PLASMA_DERIV: DATA FOR STRATUM ISTRA= ',ISTRA
          WRITE (iunout,*) 'ARE NOT AVAILABLE.'
          RETURN
        ENDIF
      END IF   ! ICALL > 0

      SELECT CASE (ITYP)
      CASE(0)
cdr photons
        IF (ASSOCIATED(PDENPH)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENPH(IOLD,IT)
              BASE_TEMP(IG) = EDENPH(IOLD,IT)/(PDENPH(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(1)
cdr atoms
        IF (ASSOCIATED(PDENA)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENA(IOLD,IT)
              BASE_TEMP(IG) = EDENA(IOLD,IT)/(PDENA(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(2)
cdr molecules
        IF (ASSOCIATED(PDENM)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENM(IOLD,IT)
              BASE_TEMP(IG) = EDENM(IOLD,IT)/(PDENM(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(3)
cdr test ions
        IF (ASSOCIATED(PDENI)) THEN
          DO IG=1,NRAD
            IT = NCLTAL(IG)
            IF (IT > 0) THEN
              BASE_DENSITY(IG) = PDENI(IOLD,IT)
              BASE_TEMP(IG) = EDENI(IOLD,IT)/(PDENI(IOLD,IT)+EPS60)
     .                        /1.5_DP
            END IF
          END DO
        END IF
      CASE(4)
cdr  bulk particles
        BASE_DENSITY(1:NRAD) = DIIN(IOLD,1:NRAD)
        BASE_TEMP(1:NRAD) = TIIN(IOLDTI,1:NRAD)
      CASE DEFAULT
        WRITE (IUNOUT,*) 'BULK SPECIES ', IPLS,
     .                   ' DENSITY MODEL: ',CDENMODEL(IPLS)
        WRITE (IUNOUT,*) ' WRONG PARTICLE TYPE SPECIFIED '
        WRITE (IUNOUT,*) ' CHECK INPUT AND RERUN CASE '
        CALL EIRENE_EXIT_OWN(1)
      END SELECT

      END SUBROUTINE EIRENE_GET_BASE_PARAMETERS



      SUBROUTINE EIRENE_GET_SPECTRUM (ICELL,IRE,SPEC,FOUND)
cdr return an (energy-resolved) spectrum SPEC in cell no. ICELL
cdr IRE: TDMPAR...(IRE)
cdr Probably unfinished code from somewhere?


      INTEGER, INTENT(IN) :: ICELL, IRE
      TYPE(EIRENE_SPECTRUM), INTENT(OUT) :: SPEC
      LOGICAL, INTENT(OUT) :: FOUND
      INTEGER :: ISPC

      FOUND = .FALSE.

      DO ISPC = 1, NADSPC
        IF ((ESTIML(ISPC)%ISRFCLL == 2) .AND.
     .      (ESTIML(ISPC)%ISPCSRF == ICELL) .AND.
     .      (ESTIML(ISPC)%IPRTYP == TDMPAR(IPLS)%TDM%ITP(IRE)).AND.
     .      (ESTIML(ISPC)%IPRSP == TDMPAR(IPLS)%TDM%ISP(IRE))) THEN
          IF (ASSOCIATED(SPEC%SPC)) THEN
            DEALLOCATE (SPEC%SPC)
            DEALLOCATE (SPEC%SDV)
            DEALLOCATE (SPEC%SGM)
          END IF
          ALLOCATE (SPEC%SPC(0:ESTIML(ISPC)%NSPC+1))
          ALLOCATE (SPEC%SDV(0:ESTIML(ISPC)%NSPC+1))
          ALLOCATE (SPEC%SGM(0:ESTIML(ISPC)%NSPC+1))
          IF (ESTIML(ISPC)%ISPCOPT==2) THEN
            ALLOCATE(SPEC%SPCAN(0:ESTIML(ISPC)%ISPLDEG,
     .                          0:ESTIML(ISPC)%NSPC+1))
          END IF
          SPEC = ESTIML(ISPC)
          FOUND = .TRUE.
        END IF
      END DO
      RETURN
      END SUBROUTINE EIRENE_GET_SPECTRUM

      END SUBROUTINE EIRENE_PLASMA_DERIV
