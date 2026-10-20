cdr  comments

cdr  may 18: some comments tried......NOT FINISHED


      subroutine eirene_emissivity(istr, lstart, lend, icall)

cdr  This routine is a generalization of old routines Ba_alpha,....,Ly-Beta.

cdr  Fill ADDV tallies with emissivities, stratum ISTR
cdr  for a set of lines ILINE = LSTART, LEND.
c
c    ADDV must have been allocated properly for (not checked here):
c    ADDV(iads),  iads = emis_lines(iline)%iadv_total,
c                        summed over lines iline=lstart,lend
c    ADDV(iadv),  iadv = emis_lines(iline)%compo(j)%iadv
c                        for iline=lstart,lend                     ! lines
c                           for j=1,emis_lines(iline)%num_compo    ! components per line
c
cdr: icall option: (perhaps not ready)
c    icall=0:  called from MCARLO, after each stratum, if NLEMIS
cdr                                and once again: for sum over strata.
cdr
c    icall=1:  called from SIGLINE in storage saving mode (MOD_ADDV=0)
c              and also if line of sight integral is defined in input block 12.
c
c
cdr  Line emissivity rates are defined in input block 4,
cdr  via ordinary "reaction decks".
c
cdr  Population coefficient for upper states (components)
cdr  must be defined via a reaction deck in block 4
cdr  either from databases (e.g.: AMJUEL, H.11 or H.12),
cdr  or via internal CR codes (H-colrad, He-colrad).
cdr  These population coefficients are transfered into here via call to OTHER_RATE_COEFF.f
cdr
cdr  NFILEN flag:
cdr  Write the newly defined tallies ADDV onto stream fort.10, fort.11, stratum ISTR



      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comsig
      use eirmod_ccona
      use eirmod_comusr, only : lgvac, tein, dein, diin, vol, nfilen
      use eirmod_comsou
      use eirmod_cgeom
      use eirmod_cgrid
      use eirmod_ctext
      use eirmod_coutau
      use eirmod_cspei
      use eirmod_ctrcei
      USE EIRMOD_CESTIM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      implicit none

      integer, intent(in) :: istr, lstart, lend, icall
      integer :: i, j, k, iline,
     .           iads, iadv, isp(3), itp(3), iratio, irc,
     .           irc_rat(2), icell, ncelc, ndens, idens
      real(dp) :: density(3), sigadd, add, powalf, powalfs,
     .            einstein, trans_en, DE, TE, TEF, DEF, popcf,
     .            EIRENE_OTHER_RATE_COEFF,
     .            ratio1, ratio2, fpop_esc
      REAL(DP) :: DUMMY(NRTAL)
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      logical :: lwrite
      CHARACTER(6) :: CISTRA
      character(len=80) :: ctest2
      EXTERNAL :: EIRENE_FTCRI, EIRENE_INTTAL, EIRENE_WRSTRT,
     .            EIRENE_LEER, EIRENE_MASBOX,
     .            EIRENE_OTHER_RATE_COEFF

      IF (TRCSIG .AND. ICALL.EQ.0) THEN
        CALL EIRENE_LEER(2)
        CALL EIRENE_FTCRI(ISTR,CISTRA)
        IF (ISTR.GT.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR STRATUM NO. '//CISTRA)
        IF (ISTR.EQ.0) CALL EIRENE_MASBOX
     .   ('SUBR. EMISSIVITY CALLED, FOR SUM OVER STRATA')
        CALL EIRENE_LEER(1)

        WRITE (iunout,*) 'AFTER INTEGRATION OVER COMPUTATIONAL DOMAIN'
      ENDIF

      do i = lstart, lend
        ILINE=I
        IF (TRCSIG .AND. ICALL.EQ.0) THEN
          ctest2 = emis_lines(iline)%line_name
          WRITE (iunout,'(1X,A,I2,3A)') 'LINE no. ',
     .                                   ILINE,', ',TRIM(CTEST2),':'
          write (iunout,'(1X,A,ES12.4)') 'EINSTEIN COEFFICIENT',
     .                                    emis_lines(i)%einstein
          write (iunout,'(1X,A,ES12.4/1x)') 'TRANSITION ENERGY   ',
     .                                    emis_lines(i)%trans_en

          WRITE (iunout,*) 'FLUX (AMP) AND POWER (WATT) BY '
        ENDIF

        einstein = emis_lines(i)%einstein
C  ENERGY FACTOR FOR POWER LOSS (W)
        trans_en = emis_lines(i)%trans_en * elcha

cdr initialize sum over components
        iads = emis_lines(i)%iadv_total
        addv(iads,:) = 0._dp
        powalfs = 0._dp

cdr run over components
        do j = 1, emis_lines(i)%num_compo
cdr  iadv: tally number on ADDV
          iadv = emis_lines(i)%compo(j)%iadv
cdr  irc:  emissivity line reaction label, as read from block 4
cdr        or set from default_emissivity.f
          irc  = emis_lines(i)%compo(j)%irc
c
          addv(iadv,:) = 0._dp
          sigadd = 0._dp
          powalf = 0._dp

cdr run over contributions: density models, isotopes, QSS states
          do k = 1, emis_lines(i)%compo(j)%num_contrib
            isp(1) = emis_lines(i)%compo(j)%contrib(k)%isp
            itp(1) = emis_lines(i)%compo(j)%contrib(k)%itp
            isp(2:3) = emis_lines(i)%compo(j)%contrib(k)%isp_rat
            itp(2:3) = emis_lines(i)%compo(j)%contrib(k)%itp_rat
            iratio = emis_lines(i)%compo(j)%contrib(k)%iratio

            irc_rat = emis_lines(i)%compo(j)%contrib(k)%irc_rat

            ndens = count(itp >= 0)
            lwrite = .true.

!  account for pop_esc
            fpop_esc = 1._dp
cdr  Option for internal CR code models only (ifit=5)
cdr  Each transition (line) can be
cdr  assigned a population escape factor.
cdr  This is then used for all calls to this CR code during the run,
cdr  e.g. for both effective rate coefficients and line emission densities
crc turned off for now
!            if (reacdat(irc)%oth%ifit == 5) then
!              fpop_esc = reacdat(irc)%oth%crm%pop_esc
!            end if

            DO ICELL=1,NSBOX
C
C  LOCAL BACKGROUND DATA ARE IN CELL ICELL
C  LOCAL TEST PARTICLE DATA ARE IN (PERHAPS COARSER) SCORING CELL NCELC
C  ACCUMULATE THE EMISSIVITIES ALSO ON THE COARSER "SCORING" GRID.
C
              NCELC=NCLTAL(ICELL)
C
              IF (NSTGRD(ICELL) > 0) CYCLE
              IF (LGVAC(ICELL,NPLS+1)) CYCLE

              TE=TEIN(ICELL)
              DE=DEIN(ICELL)

              DEF=LOG(DE)
              TEF=max(-2.30_DP,LOG(TE)) ! cut-off at 0.1 eV

              do idens = 1, ndens
                select case (itp(idens))
                  case (0)
                    density(idens) = pdenph(isp(idens),ncelc)
                  case (1)
                    density(idens) = pdena(isp(idens),ncelc)
                  case (2)
                    density(idens) = pdenm(isp(idens),ncelc)
                  case (3)
                    density(idens) = pdeni(isp(idens),ncelc)
                  case (4)
                    density(idens) = diin(isp(idens),icell)
                  case (5)
                    density(idens) = dein(icell)
                  case default
                    density(idens) = 0._dp
                    if (lwrite) then
                      write (iunout,*) ' ERROR IN EMISSIVITY'
                      write (iunout,*)
     .                  ' WRONG PARTICLE TYPE SPECIFIED FOR'
                      write (iunout,*) ' line ',i,
     .                   emis_lines(i)%line_name
                      write (iunout,*) ' component ',j,
     .                   emis_lines(i)%compo(j)%compo_name
                      write (iunout,*) ' contribution ',k
                      lwrite = .false.
                    end if
                end select
              end do

c  population coefficient, relative to density(1)
              popcf= EIRENE_OTHER_RATE_COEFF(IRC,ICELL,TEF,DEF,.TRUE.,1)
              add = popcf*density(1)

c  density ratio, if true parent density is not available (or in QSS mode)
c  then: ratio1 converts from density(1) to density
c  density is the "true" parent density for this component.
c  density(1) is taken as "intermediate" parent density. Fetch reduced population coefficent
c  and density ratio  ratio1="density"/"density(1)" will be applied,
c  to turn density(1) into "density"
c  e.g. density    = H2+
c       density(1) = H2
c       ratio1     = [H2+]/[H2]
c  this works when the second species involved in loss and gain rates
c  for species H2+ from H2 is the same, here: electron density, and hence cancels.
              if (iratio > 0) then

                ratio1 = EIRENE_OTHER_RATE_COEFF(IRC_RAT(1),ICELL,
     .                                          TEF,DEF,.TRUE.,1)
                add = add*ratio1
c
c  second conversion to yet another parent density
c  e.g: density    = H3+.   = [H2+] * [H2/ne] *ratio2 = [H2] * ratio1 * [H2/ne] *ratio2
c       density(1) = H2
c       ratio1     = H2+/H2(Te,ne) (CR equilibrium)
c       ratio2     = prod[H3+] from H2 impact/loss[H3+] from elec. impact (DR)
c  this works when the second species involved in loss and gain rate
c  for species H3+ from H2+ is not the same,
c  here: electron density, and H2 density, hence: does not cancel.
                if (iratio == 2) then
                  ratio2 = EIRENE_OTHER_RATE_COEFF(IRC_RAT(2),ICELL,
     .                                             TEF,DEF,.TRUE.,1)
                  add = add * density(2) / density(3) *ratio2
                end if
              end if

cdr so far: add is scored on the fine grid cell "icell".
cdr         add volume-weighted contribution to coarse cell "ncelc"
              sigadd = add * einstein * fpop_esc * vol(icell)

              addv(iadv,ncelc) = addv(iadv,ncelc) + sigadd
              addv(iads,ncelc) = addv(iads,ncelc) + sigadd

              powalf = powalf + sigadd
            end do               ! icell

          end do ! k contributions (summed) of component j of line iline

cdr ADDV was volume-weighted (extensive) sum. Now divide by coarse cell volume
cdr      to turn it into an intensive score: [...] per cm**3
          addv(iadv,1:nsbox_tal) = addv(iadv,1:nsbox_tal)
     .                             / voltal(1:nsbox_tal)

          powalf = powalf * trans_en
          powalfs = powalfs + powalf

          if (TRCSIG .AND. ICALL.EQ.0)
     .      WRITE (iunout,'(A50,2ES16.7)') ' COUPL. TO ' //
     .                     TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
     .                    ,POWALF/TRANS_EN*ELCHA,POWALF

          DUMMY(1:NSBOX_TAL) = ADDV(IADV,1:NSBOX_TAL)
          CALL EIRENE_INTTAL
     .         (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADV,ISTR),
     .          NR1TAL,NP2TAL,NT3TAL,NBMLT)
          ADDV(IADV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

          TXTTAL(IADV,NTALA) =REPEAT(' ',72)
          TXTTAL(IADV,NTALA) =TRIM(EMIS_LINES(I)%LINE_NAME) // ', ' //
     .                    ' SOURCE RATE '
          TXTSPC(IADV,NTALA) =TRIM(EMIS_LINES(I)%COMPO(J)%COMPO_NAME)
          TXTUNT(IADV,NTALA) ='PHOTONS/S/CM**3         '
          IF (TRCSIG) THEN
            WRITE (iunout,*) ' TALLY ADDV(IADV) prepared. IADV=',IADV
            CALL EIRENE_LEER(1)
          END IF

        end do ! j components of line ILINE are done

cdr  now sum over components: on tally ADDV(IADS)

        addv(iads,1:nsbox_tal) = addv(iads,1:nsbox_tal)
     .                           / voltal(1:nsbox_tal)
        IF (TRCSIG .AND. ICALL.EQ.0)
     .   WRITE (iunout,'(A50,2ES16.7)')
     .                  ' TOTAL FLUX (AMP) AND POWER (WATT) '
     .                  ,POWALFS/TRANS_EN*ELCHA,POWALFS

        DUMMY(1:NSBOX_TAL) = ADDV(IADS,1:NSBOX_TAL)
        CALL EIRENE_INTTAL
     .       (DUMMY,VOLTAL,1,1,NSBOX_TAL,ADDVI(IADS,ISTR),
     .        NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(IADS,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)

        TXTTAL(IADS,NTALA) =REPEAT(' ',72)
        TXTTAL(IADS,NTALA) ='SUM OVER COMPONENTS  '
        TXTSPC(IADS,NTALA) ='ALL COMPONENTS'
        TXTUNT(IADS,NTALA) ='PHOTONS/S/CM**3         '
        IF (TRCSIG) THEN
          WRITE (iunout,*) ' TALLY ADDV(IADV) prepared. IADV=',IADS
          CALL EIRENE_LEER(2)
        END IF

      end do ! line no. ILINE

C
C  WRITE ON STREAM 10 DATA FOR STRATUM NO. ISTR
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTR
        CALL EIRENE_WRSTRT(ISTR,NSTRAI,
     .              NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .              TRCFLE)
C
C  WRITE ON STREAM 11 (TOTALS, SUMMED OVER GRID) FOR STRATUM NO. ISTR
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 11  IRC= ',IRC

C  WRITE ON STREAM 10 ONLY DATA FOR SUM OVER STRATA
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_WRSTRT(ISTR,NSTRAI,
     .              NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .              TRCFLE)
C
C  WRITE ON STREAM 11 (TOTALS, SUMMED OVER GRID) ONLY DATA FOR SUM OVER STRATA
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ENDIF
C
      RETURN

      end subroutine eirene_emissivity
