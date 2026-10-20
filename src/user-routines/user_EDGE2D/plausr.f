c
c modified: s.wiesen@fz-juelich.de
c
C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| reads target plasma information from casename.zplasma         |
C| sets step-functions                                           |
C| reads also neutral particle fluxes from 'eirene.chemFluxDep'  |
C| and sets FLXOUT for flux dependency of chemical sputtering.   |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added reading of neutral fluxes from|
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C| 22/03/2011   D.Harting    Do not exit anymore if one of the   |
C|                           geometrical parameters in           |
C|                           eirene.chemFluxDep are not matching.|
C|                           Just ignore the neutral flux from   |
C|                           the previous run and continue.      |
C+---------------------------------------------------------------+
      SUBROUTINE EIRENE_PLAUSR
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CSTEP
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CINIT
      USE EIRMOD_COMSOU
      USE EIRMOD_CTRIG
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_LEARC1, ONLY: EIRENE_LEARC1
      USE EIRMOD_USRDATA

      IMPLICIT NONE

      REAL(DP) :: FACTOR, EIRENE_STEP
      INTEGER :: NLINES, ITRI, ISIDE, I, NBIN, ISTRA, ISRFS, ISOR, JJJ,
     .           ITEC1, ITEC2, ITEC3, ISTEP, INDSRF, IS1, IERROR, IPLS
      INTEGER :: EIRENE_IDEZ, MSURFG
      INTEGER, ALLOCATABLE :: KSTEP(:), INOSRC(:), IPLAN(:), IPLEN(:)
      REAL(DP) :: FLX, TE, TI, DE, MC, FE, FI, FSH, VP, FEL, ZI, DUM
      REAL(DP) :: DELR, FL, MCC, FFEL, CS, vx,vy,vz,di,usrval
      real(dp) :: xref, yref, bzref, facbz, x, y, rad, bx, by, bz, bf,
     .            errbx, errby, errbz, errbf, bxmax, bymax, bzmax,
     .            bfmax, dfdx, dfdy, dfdz, xref2, yref2, bzref2, facbz2,
     .            xref3, yref3, bzref3, facbz3
      integer :: nref, icell, nplcll, ipolg, nref2, nref3
      CHARACTER(256) :: line,sstr,filename
      character(2) :: cstr2
      integer, allocatable :: indextmp(:),itritmp(:),isidetmp(:)
      integer, allocatable :: isegtmp(:), offsets(:)
      real(dp), allocatable :: tetmp(:),fetmp(:),fshtmp(:)
      REAL(DP), ALLOCATABLE :: PSI_VALS(:), PSI_VALS_CORNER(:), copy(:)

      integer, parameter :: fp=31
      integer :: ll,ier,j,iseg,ind,fp2
      integer :: ftell

      integer :: ntr, NLIM_tmp, NSTS_tmp, NGITT_tmp, NGSTAL_tmp,
     &     NATM_tmp, NMOL_tmp, IPLS_tmp, NTR_tmp, NLIMPS_tmp, ISTEP_SPEZ
      REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: hydIonFLX, EIRENE_wall_area
      REAL(DP),ALLOCATABLE,DIMENSION(:)   :: hydNeutFLX
      INTEGER, ALLOCATABLE,DIMENSION(:,:) :: hydNeutFLX_info
      logical :: lex,dbg_out
      real(dp):: leng, tmp

      real(dp) :: neutralFluxFileVersion

      interface
        subroutine EIRENE_cell_to_corner (f, fcorner)
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: f(:)
          real(dp), intent(out) :: fcorner(:)
        end subroutine EIRENE_cell_to_corner

        subroutine EIRENE_df_dxyz (fcorner, icell, x, y, z,
     .                      dfdx, dfdy, dfdz)
          use eirmod_precision
          implicit none
          real(dp), intent(in) :: fcorner(:), x, y, z
          real(dp), intent(out) :: dfdx, dfdy, dfdz
          integer, intent(in) :: icell
        end subroutine EIRENE_df_dxyz
      end interface

      CALL EIRENE_ALLOC_CSTEP
      ALLOCATE (KSTEP(NSTEP))
      ALLOCATE (INOSRC(NSTEP))
      ALLOCATE (IPLAN(NSTEP))
      ALLOCATE (IPLEN(NSTEP))

      KSTEP = 0
      INOSRC = 0
      IPLAN = 0
      IPLEN = 0

      dbg_out=.false.
      NeutralFluxFileVersion = 1.2_DP

c begin dmh added 21.06.2010
      ll=len_trim(casename)
!
      filename=casename(1:ll) // '.zplasma'
      open (unit=fp+ifoff,file=filename,access='stream',
     .     form='formatted')

c first get number of triangles from misc plasma data:
        write(sstr,'(a20)')
     .          '*** MISC PLASMA DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
           write(*,*) 'PLAUSR: ',sstr,' not found'
           close(fp+ifoff)
           call EIRENE_exit_own(1)
        endif
        do j=1,2
           read(fp+ifoff,'(a)') line
        enddo

        read(fp+ifoff,'(i8)') ntr
        if (ntr /= nr1st-1) then
           write (*,*) 'PLAUSR:', sstr
           write (*,*) ' wrong number of triangles in plasma file'
           write (*,*) ' check for correct number in file ',filename
           call EIRENE_exit_own(1)
        endif
c allocate temporary array to store plasma flux (hydrogen isotope)
        allocate(hydIonFLX(3,ntr))
        hydIonFLX(:,:) = 0.0_DP
c allocate temporary array to store wall area of elements
        allocate(EIRENE_wall_area(3,ntr))
        EIRENE_wall_area(:,:) = 0.0_DP
c end dmh added 21.06.2010

c firstly, read misc target data
      write(sstr,'(a20)') '*** MISC TARGET DATA'
      CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
      if(ier /=0) then
         write(*,*) 'PLAUSR: ',sstr,' not found'
         close(fp+ifoff)
         CALL EIRENE_exit_own(1)
      endif
      do j=1,2
         read(fp+ifoff,'(a)') line
      enddo
      READ (fp+ifoff,*) NLINES
cswx 24sep07
      if(nlines == 0) then
        close(fp+ifoff)
        return
      endif
cswx
      if(nlines <= 0) then
         write(*,*) ' PLAUSR: error, nlines=',nlines
         close(fp+ifoff)
         call EIRENE_exit_own(1)
      endif

      allocate(indextmp(nlines))
      allocate(itritmp(nlines))
      allocate(isidetmp(nlines))
      allocate(tetmp(nlines))
      allocate(fetmp(nlines))
      allocate(fshtmp(nlines))
      allocate(isegtmp(nlines))
      allocate(offsets(npls_fix))
      do j=1,nlines
         read(fp+ifoff,'(3i7,3(1x,e14.7))')
     .        indextmp(j),
     .        itritmp(j),
     .        isidetmp(j),
     .        tetmp(j),
     .        fetmp(j),
     .        fshtmp(j)
         isegtmp(indextmp(j)) = j
      enddo
      ! scan to find offsets of ions strata data

      do ipls = 1, npls_fix
        write(cstr2,'(i2.2)') ipls
        write(sstr,'(a23)')
     .        '*** ION #'//cstr2//' TARGET DATA'
        CALL EIRENE_locstr_usr(fp+ifoff,sstr,ier)
        if(ier /=0) then
            write(*,*) 'PLAUSR: tag ',sstr,'not found'
            close(fp+ifoff)
            call EIRENE_exit_own(1)
        endif
        do j=1,2
            read(fp+ifoff,'(a)') line
        enddo
        inquire(unit=fp+ifoff, pos=offsets(ipls))
      enddo


c now read species-dependent target data
      DO ISTRA=1,NSTRAI
         IF (.NOT.NLSRF(ISTRA)) CYCLE
         DO ISRFS=1,NSRFSI(ISTRA)

c     get step function index ISTEP
            ISOR=SORLIM(ISRFS,ISTRA)
            ITEC1=EIRENE_IDEZ(ISOR,1,4)
            ITEC2=EIRENE_IDEZ(ISOR,2,4)
            ITEC3=EIRENE_IDEZ(ISOR,3,4)
            IF ((ITEC1 /= 4).AND.(ITEC2 /= 4).AND.(ITEC3 /= 4)) CYCLE
            ISTEP=SORIND(ISRFS,ISTRA)
            ISTEP=MOD(NINT(SORIND(ISRFS,ISTRA)),100)
            ISTEP_SPEZ=INT(SORIND(ISRFS,ISTRA) /100)
            IF(ISTEP_SPEZ==0) ISTEP_SPEZ = nspez(istra)
            IF (ISTEP.EQ.0) THEN
              WRITE (*,*) 'ERROR IN PRIMARY SOURCE DATA '
              WRITE (*,*) 'STEPFUNCTION REQUESTED FOR SOURCE SURFACE '
              WRITE (*,*) 'NO. ',INSOR(ISRFS,ISTRA),' BUT SORIND.EQ.0.'
              CALL EIRENE_EXIT_own(1)
            ELSEIF (ISTEP.GT.NSTEP) THEN
              CALL EIRENE_MASPRM('NSTEP',5,NSTEP,'ISTEP',5,ISTEP,IERROR)
              CALL EIRENE_EXIT_own(1)
            ENDIF
c     initialise stratum strength
            IF(ISRFS==1) FLUX(ISTRA) = 0.0_DP

c     get surface index INDSRF
            INDSRF=INSOR(ISRFS,ISTRA)
            IF (INDSRF < 0) INDSRF=NLIM+ABS(INDSRF)

c     get species index/indices IPLAN(ISTEP) --> IPLEN(ISTEP)
            IF (NSPEZ(ISTRA) <= 0) THEN
               IPLAN(ISTEP)=1
               IPLEN(ISTEP)=npls_fix
               ipls = 0
            ELSE
               IPLAN(ISTEP)=NSPEZ(ISTRA)
               IPLEN(ISTEP)=NSPEZ(ISTRA)
               IPLAN(ISTEP)=ISTEP_SPEZ
               IPLEN(ISTEP)=ISTEP_SPEZ
               ipls=ISTEP_SPEZ
            END IF
c     fudge species index for atomic impurity flux (get it from NEMODS index K)
            IF (NLATM(ISTRA)) THEN
              IPLS_tmp = EIRENE_IDEZ(NEMODS(ISTRA),4,4)
              IF (IPLS_tmp.gt.0) then
                ipls=IPLS_tmp
              ENDIF
            ENDIF
c     search target tag in .zplasma file
           IPLS_tmp = IPLS
           rewind(fp+ifoff)

c     read step functions
            KSTEP(ISTEP) = 0
            READ (fp+ifoff,*,pos=offsets(ipls)) NLINES
            DO I=1, NLINES
               read(fp+ifoff,'(i7,12(1x,e14.7))')
     .              ind,
     .              flx,ti,di,
     .              vx,vy,vz,
     .              fi,fel,
     .              vp,mc,zi,
     .              usrval

               iseg = isegtmp(ind)
               itri = itritmp(iseg)
               iside= isidetmp(iseg)
               te   = tetmp(iseg)
               fe   = fetmp(iseg)
               fsh  = fshtmp(iseg)

c begin added dmh 21.06.2010
c     store ion flux in [A] of main plasma
               if(ipls_tmp.eq.1) then
                  hydIonFLX(iside,itri) = flx
               endif
c end added dmh 21.06.2010

               IF (INMTI(ISIDE,ITRI) == INDSRF) THEN
                  IF (KSTEP(ISTEP) == 0) RRSTEP(ISTEP,1) = 0.0_DP
                  KSTEP(ISTEP) = KSTEP(ISTEP) + 1
                  INOSRC(ISTEP) = ISTRA
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1
                  IRSTEP(ISTEP,KSTEP(ISTEP))=ITRI
                  IPSTEP(ISTEP,KSTEP(ISTEP))=ISIDE
                  ITSTEP(ISTEP,KSTEP(ISTEP))=1
                  IASTEP(ISTEP,KSTEP(ISTEP))=0
                  IBSTEP(ISTEP,KSTEP(ISTEP))=1
                  DELR =  SQRT(
     .                 (XTRIAN(NECKE(ISIDE,ITRI))
     .                 -XTRIAN(NECKE(IS1,ITRI)))**2+
     .                 (YTRIAN(NECKE(ISIDE,ITRI))
     .                 -YTRIAN(NECKE(IS1,ITRI)))**2)
                  RRSTEP(ISTEP,KSTEP(ISTEP)+1)=
     .                 RRSTEP(ISTEP,KSTEP(ISTEP)) + DELR
                  TESTEP(ISTEP,KSTEP(ISTEP)) = TE
                  FESTEP(ISTEP,KSTEP(ISTEP)) = FE
C     IF NO SHEATH POTENTIAL SPECIFIED, DERIVE IT FROM ELECTRON ENERGY
C     FLUX BY SUBTRACTING THE KINETIC CONTRIBUTION 2.0*TE
                  IF (FSH.EQ.0..AND.FE.GE.2.0) FSH=FE-2.0
                  SHSTEP(ISTEP,KSTEP(ISTEP)) = FSH
                  DO IPLS=IPLAN(ISTEP), IPLEN(ISTEP)
                     ZISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ZI ! NCHARP(IPLS)
                     TISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = TI ! eV
                     DISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = DI ! 1/cm**3
                     FISTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FI !   1
                     VPSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(VP) ! cm/s
C     VP OVERRULES MC, IF VP IS GIVEN and MC=0
                     MCC=0.0_DP
                     IF (VP.NE.0.0_DP) THEN
                        CS=SQRT((TI+TE)/RMASSP(IPLS))*CVEL2A
                        IF (CS.GT.0.0_DP) MCC=VP/CS
                     ENDIF
                     jjj=kstep(istep)
                     IF (MC.EQ.0.) MC=MCC
!pb                     MCSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = abs(MC) ! 1
C     THIS NEXT VECTOR IS V-PARALLEL, IN CARTESIAN COORDINATES
                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VXIN(IPLS,ITRI)
                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VYIN(IPLS,ITRI)
                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZIN(IPLS,ITRI)
c                     VXSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VX
c                     VYSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VY
c                     VZSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = VZ
                     FLSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = ABS(FLX)/DELR
                     FLUX(ISTRA) = FLUX(ISTRA) + ABS(FLX)
C     IF NO ION KINETIC ENERGY FLUX IS SPECIFIED, DERIVE IT FROM  FI,MC,VP
                     ffel=0.0_DP
                     IF (FI.gt.0.0_DP.or.mc.gt.0.0_DP) then
                        ffel=(FI*TI+0.5_DP*MC*MC*(TE+TI))*abs(FLX)
                     endif
                     IF (FEL.EQ.0.0_DP) FEL=FFEL
                     IF (eirene_use_elstepdat_bug) THEN
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL
                     ELSE
                        ELSTEP(IPLS,ISTEP,KSTEP(ISTEP)) = FEL/DELR
                     ENDIF
csw
c                     usrstep(ipls,istep,kstep(istep)) = usrval
csw
c                 enddo ipls,iplan
                  ENDDO
C
C                 SET PLASMA PARAMETERS FOR ATOMIC STRATA
                  IF (NLATM(ISTRA)) THEN
                     TISTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) = TI ! eV
                     DISTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) = DI ! 1/cm**3
                     VPSTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) = abs(VP) ! cm/s
                     VXSTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) =
     &                    VXIN(IPLS_TMP,ITRI)
                     VYSTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) =
     &                    VYIN(IPLS_TMP,ITRI)
                     VZSTEP(IPLS_TMP,ISTEP,KSTEP(ISTEP)) =
     &                    VZIN(IPLS_TMP,ITRI)
                  ENDIF

c              endif inmti
               ENDIF

c           enddo nlines
            ENDDO
c        enddo isrfs
         ENDDO
c     enddo istra
      ENDDO



      DO ISTEP = 1, NSTEP
         FL=0.0_DP
         IF (KSTEP(ISTEP) > 0) THEN
            NBIN=KSTEP(ISTEP)+1
            IPLAN(ISTEP) = 1
            IPLEN(ISTEP) = npls_fix
            FL=EIRENE_STEP(IPLAN(ISTEP),IPLEN(ISTEP),NBIN,ISTEP,4)
         END IF
      END DO

c begin added dmh 21.06.2010 for flux dependency of chemical sputtering
c     read neutral flux [A] to target and walls from last EIRENE run

      allocate(hydNeutFLX(NGITT))
      hydNeutFLX(:) = 0.0_DP
      allocate(hydNeutFLX_info(NGITT,3))
      hydNeutFLX_info = 0
      NLIM_tmp   = NLIM
      NSTS_tmp   = NSTS
      NGITT_tmp  = NGITT
      NGSTAL_tmp = NGSTAL
      NATM_tmp   = NATM
      NMOL_tmp   = NMOL
      NLIMPS_tmp = NLIMPS
      NTR_tmp    = NTR
      fp2 = 4999
      inquire(file=trim(eirene_fstoreneutflux),exist=lex)
      if (eirene_wallFluxModel.ge.2.and.lex) then
          open(unit=fp2,file=trim(eirene_fstoreneutflux),
     &        access='sequential')
         read(fp2,'(a)') line
         read(line,'(a28)') sstr

c        check if netral flux file is compatible with actual code version
         if (index(sstr,"* Neutral flux file version:").ne.1) then
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
            write(IUNOUT,*) "Found obsolete neutral flux file: ",
     &           trim(eirene_fstoreneutflux)
            write(IUNOUT,*) "Ignoring neutral flux from previous run"
            WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
            lex =.false.
         else
            read(line,'(a28,f14.6)') sstr,tmp
            if (tmp.ne.NeutralFluxFileVersion) then
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
               write(IUNOUT,'(a,a,a,f14.6)')
     &              "Found obsolete neutral flux file: ",
     &              trim(eirene_fstoreneutflux),"; version:",tmp
               write(IUNOUT,'(a,f14.6)')
     &              "but actual neutral flux file should be version:",
     &              NeutralFluxFileVersion
               write(IUNOUT,*) "Ignoring neutral flux from previous run"
               WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
               lex =.false.
            else
c              neutral flux file is compatible with actual code version, so read it
               read(fp2,'(a)') line
               read(fp2,'(8i8)') NLIM_tmp, NSTS_tmp, NGITT_tmp,
     &              NGSTAL_tmp, NATM_tmp, NMOL_tmp, NLIMPS_tmp, NTR_tmp
               if ((NLIMPS_tmp-NLIM_tmp-NSTS_tmp.ne.NLIMPS-NLIM-NSTS)
     &              .or.(NSTS_tmp.ne.NSTS).or.(NTR_tmp.ne.NTR).or.
     &              (NGITT_tmp.ne.NGITT).or.(NGSTAL_tmp.ne.NGSTAL)) then
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning"
                  WRITE(IUNOUT,*) "Neutral flux file ",
     &                 trim(eirene_fstoreneutflux),
     &                 " does not fit to simulation!"
                  WRITE(IUNOUT,*)"NLIM = ",NLIM,"; NLIM_tmp = ",NLIM_tmp
                  WRITE(IUNOUT,*)"NSTS = ",NSTS,"; NSTS_tmp = ",NSTS_tmp
                  WRITE(IUNOUT,*)"NGITT = ",NGITT,
     &                 "; NGITT_tmp = ",NGITT_tmp
                  WRITE(IUNOUT,*)"NGSTAL = ",NGSTAL,
     &                 "; NGSTAL_tmp = ",NGSTAL_tmp
                  WRITE(IUNOUT,*)"NATM = ",NATM,"; NATM_tmp = ",NATM_tmp
                  WRITE(IUNOUT,*)"NMOL = ",NMOL,"; NMOL_tmp = ",NMOL_tmp
                  WRITE(IUNOUT,*)"NTRII= ",NTR,"; NTRII_tmp = ",NTR_tmp
                  WRITE(IUNOUT,*)"NLIMPS_tmp-NLIM_tmp-NSTS_tmp = ",
     &                 NLIMPS_tmp-NLIM_tmp-NSTS_tmp,
     &                 "NLIMPS-NLIM-NSTS = ",NLIMPS-NLIM-NSTS
                  write(IUNOUT,*)
     &                 "Ignoring neutral flux from previous run"
                  WRITE(IUNOUT,*) "* EIRENE_PLAUSR: Warning end"
                  lex =.false.
               endif

               if (lex) then
                  read(fp2,'(a)') line
                  read(fp2,'(a)') line
                  do ISRFS=1,NLIMPS
                    do i= 1, surf_trian(ISRFS)%numtr
                      itri = surf_trian(ISRFS)%itrias(i)
                      iside = surf_trian(ISRFS)%itrisi(i)

                      read(fp2,'(i6,1x,2(e14.6,1x),i6,1x,i6,1x,i6)')
     &                    j, hydNeutFLX(j), tmp,
     &                    hydNeutFLX_info(j,1), hydNeutFLX_info(j,2),
     &                    hydNeutFLX_info(j,3)
                     if (INSPAT(iside,itri).ne.j) then
                        WRITE(IUNOUT,*) "* EIRENE_PLAUSR:"
                        write(IUNOUT,*) "Error reading from file: ",
     &                       trim(eirene_fstoreneutflux)
                        write(IUNOUT,*) "INSPAT=",INSPAT(iside,itri),
     &                                  ";  in file j = ",j
                        CALL EIRENE_EXIT_own(1)
                     endif
                   enddo        ! isrfs
                 enddo          ! i
               endif            !(lex)
            endif               !(tmp.ne.NeutralFluxFileVersion)
         endif                  !(index(sstr,"* Neutral flux file version:").ne.1)
         close(fp2)
      endif

c     calculate target flux for chemical sputtering of Roth-formula
      FLXOUT(:) = 0.0_DP

c     use hydrogen ion flux [A] to wall
c     FLXOUT is needed in #/(cm^2 s) (convert from A to #/(cm^2 s))
      if (eirene_wallFluxModel.ge.1) then
        do ISRFS=1,NLIMPS
          do i= 1, surf_trian(ISRFS)%numtr
            itri = surf_trian(ISRFS)%itrias(i)
            iside = surf_trian(ISRFS)%itrisi(i)
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1
            leng  = SQRT(
     .             (XTRIAN(NECKE(ISIDE,ITRI))
     .             -XTRIAN(NECKE(IS1,ITRI)))**2+
     .             (YTRIAN(NECKE(ISIDE,ITRI))
     .             -YTRIAN(NECKE(IS1,ITRI)))**2)
            EIRENE_wall_area(iside,itri) = PIA*leng
     &             *( XTRIAN(NECKE(ISIDE,ITRI))
     &             +  XTRIAN(NECKE(IS1,ITRI)) )
            if (hydIonFLX(iside,itri).eq.0) cycle
            MSURFG=NLIM+NSTS+INSPAT(iside,itri)
            FLXOUT(MSURFG) = FLXOUT(MSURFG)+DABS(hydIonFLX(iside,itri))
     &              /(1.6022D-19*EIRENE_wall_area(iside,itri))
          enddo               ! iside
        enddo                  ! itri
      endif

      if ((eirene_wallFluxModel.ge.2).and.lex) then
c     use hydrogen neutral flux [A] to wall
         do i = 1, NGITT
            if (hydNeutFLX(i).eq.0) cycle
            if (hydNeutFLX_info(i,1).eq.0) cycle
            itri  = hydNeutFLX_info(i,1)
            iside = hydNeutFLX_info(i,2)
            MSURFG=NLIM+NSTS+INSPAT(iside,itri)
            FLXOUT(MSURFG) = FLXOUT(MSURFG) + DABS(hydNeutFLX(i))
     &              /(1.6022D-19*EIRENE_wall_area(iside,itri))
         enddo                  ! i
      endif

c end added dmh 21.06.2010 for flux dependency of chemical sputtering

c     Debug output of boundary
      if (dbg_out) then
      open(unit=fp2, file="eirene.chemSput_dbgOut",access='sequential')
      write(fp2,*) NLIMPS_tmp
      write(fp2,'(a,a,a)')
     &     "* index, nsurf,           R1,           R2,",
     &     "           Z1,           Z2,       FLXOUT,",
     &     "      NeutFLX,       IonFLX,         area"
      do i=1,NGITT
         lex=.false.
         do itri=1,ntr
            do iside=1,3
               if (INSPAT(iside,itri).eq.i) then
                  if (lex) write(fp2,*)"*Edge twice found"
                  lex=.true.
                  IS1 = ISIDE + 1
                  IF (IS1.GT.3) IS1=1
                  write(fp2,'(I8,I7,8(x,e13.6))')
     &                 i,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI)),
     &                 FLXOUT(i)/
     &                 (1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydNeutFLX(i)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 hydIonFLX(iside,itri)
     &                 /(1.6022D-19*EIRENE_wall_area(iside,itri)),
     &                 EIRENE_wall_area(iside,itri)
               endif
            enddo
         enddo
         if (.not.lex) then
            write(fp2,'(I8,I7,8(x,e13.6))') i,
     &     0,0.0_DP,0.0_DP,0.0_DP,0.0_DP,FLXOUT(i),hydNeutFLX(i),0.0_DP
         endif
      enddo
      close(fp2)


      open(unit=fp2, file="eirene.chemSput_dbgOut2",access='sequential')
      write(fp2,*) ntr
      write(fp2,'(a,a)')
     &     "*   itri, iside, nsurf,           R1,           R2,",
     &     "           Z1,           Z2"
      do itri=1,ntr
         do iside=1,3
            IS1 = ISIDE + 1
            IF (IS1.GT.3) IS1=1
            write(fp2,'(I8,x,I6,x,I6,4(x,e13.6))')
     &                 itri,iside,INMTI(iside,itri),
     &                 XTRIAN(NECKE(ISIDE,ITRI)),
     &                 XTRIAN(NECKE(IS1,ITRI)),
     &                 YTRIAN(NECKE(ISIDE,ITRI)),
     &                 YTRIAN(NECKE(IS1,ITRI))
         enddo
      enddo
      close(fp2)

c dmh begin added output of step function data
      open(unit=fp2, file="eirene.stepdat_dbg",access='sequential')
      write(fp2,"(A,A,A)") "ipls istep i  IRSTEP(ITRI)  IPSTEP(ISIDE)",
     &     "  X1   Y1  DR  FLSTEP  ELSTEP  SHSTEP  TISTEP",
     &     "  TESTEP"
      do ipls=1,npls_fix
         do istep=1,nstep
            do i=1,KSTEP(ISTEP)
               ITRI = IRSTEP(ISTEP,I)
               ISIDE= IPSTEP(ISTEP,I)
               IS1 = ISIDE+1
               if(IS1.GT.3) IS1=1
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE,
     &              XTRIAN(NECKE(ISIDE,ITRI)),YTRIAN(NECKE(ISIDE,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i),
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i),
     &              TESTEP(ISTEP,i)
               write(fp2,"(2(I3,x),I4,x,I6,x,I1,8(x,e13.6))")
     &              ipls, istep, i, ITRI, ISIDE,
     &              XTRIAN(NECKE(IS1,ITRI)), YTRIAN(NECKE(IS1,ITRI)),
     &              (RRSTEP(istep,i+1)-RRSTEP(istep,i)),
     &              FLSTEP(IPLS,ISTEP,i), ELSTEP(IPLS,ISTEP,i),
     &              SHSTEP(ISTEP,i), TISTEP(IPLS,ISTEP,i),
     &              TESTEP(ISTEP,i)
               write(fp2,*)
            enddo
         enddo
      enddo
      close(fp2)
      endif
c dmh begin added output of step function data

      deallocate(hydNeutFLX)
      deallocate(hydNeutFLX_info)

c     cleanup
      DEALLOCATE (KSTEP)
      DEALLOCATE (INOSRC)
      DEALLOCATE (IPLAN)
      DEALLOCATE (IPLEN)
      deallocate(indextmp)
      deallocate(itritmp)
      deallocate(isidetmp)
      deallocate(tetmp)
      deallocate(fetmp)
      deallocate(fshtmp)
      deallocate(isegtmp)
      deallocate(hydIonFLX)
      deallocate(EIRENE_wall_area)
      deallocate(offsets)

      close(fp+ifoff)

      return

      END SUBROUTINE EIRENE_PLAUSR
